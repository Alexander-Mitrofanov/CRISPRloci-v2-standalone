"""Immutable package resources, identities and model compatibility checks."""
import hashlib
import json
from importlib.resources import files

import xgboost as xgb

from .errors import ModelError
from .features import FEATURE_NAMES

PROVENANCE_SHA256 = "dd7d3765e5eb9e32fbab712dee21a26a5c6be63aad100d5515b7f628fd13bc43"


def load_model() -> tuple[xgb.Booster, dict]:
    try:
        data = files("repeattyper").joinpath("data")
        provenance_bytes = data.joinpath("provenance.json").read_bytes()
        if hashlib.sha256(provenance_bytes).hexdigest() != PROVENANCE_SHA256:
            raise ModelError("Bundled provenance checksum mismatch")
        provenance = json.loads(provenance_bytes)
        resources = {}
        for filename, expected in provenance["files"].items():
            payload = data.joinpath(filename).read_bytes()
            if hashlib.sha256(payload).hexdigest() != expected:
                raise ModelError(f"Bundled model resource checksum mismatch: {filename}")
            resources[filename] = payload
        booster = xgb.Booster({"nthread": 1, "device": "cpu"})
        booster.load_model(bytearray(resources["xgb_repeats.json"]))
        rows = [row.split(":") for row in resources["type_dict.tab"].decode("ascii").splitlines()]
        labels = [row[0] for row in rows]
        indices = [int(row[1]) for row in rows]
        configuration = json.loads(booster.save_config())["learner"]
        if (indices != list(range(37)) or labels != provenance["labels"]
                or booster.num_features() != len(FEATURE_NAMES) or booster.attr("best_iteration") != "41"
                or configuration["objective"]["name"] != "multi:softprob"
                or configuration["learner_model_param"]["num_class"] != "37"):
            raise ModelError("Model feature, label or iteration contract mismatch")
    except ModelError:
        raise
    except (xgb.core.XGBoostError, OSError, ValueError, KeyError, TypeError, IndexError) as exc:
        raise ModelError("Bundled classifier is incompatible") from exc
    return booster, provenance
