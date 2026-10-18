"""Pure-NumPy inference for CRISPRidentify's compact classifier."""

from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np

from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.feature_schema import validate_compact_feature_schema


FEATURE_SCHEMA = COMPACT_FEATURE_SCHEMA
MODEL_ROOT = Path(__file__).resolve().parents[1] / "trained_models"
MANIFEST_PATH = MODEL_ROOT / "compact_model.json"


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class CompactHistGradientModel:
    def __init__(self, archive):
        self.baseline = float(archive["baseline"][0])
        self.learning_rate = float(archive["learning_rate"][0])
        self.n_features = int(archive["n_features"][0])
        self.tree_offsets = np.asarray(
            archive["tree_offsets"], dtype=np.int64
        )
        self.values = np.asarray(archive["values"], dtype=float)
        self.features = np.asarray(archive["features"], dtype=np.int64)
        self.thresholds = np.asarray(archive["thresholds"], dtype=float)
        self.left = np.asarray(archive["left"], dtype=np.int64)
        self.right = np.asarray(archive["right"], dtype=np.int64)
        self.is_leaf = np.asarray(archive["is_leaf"], dtype=bool)
        self.missing_left = np.asarray(
            archive["missing_left"], dtype=bool
        )

    def predict(self, feature_matrix):
        matrix = np.asarray(feature_matrix, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != self.n_features:
            raise ValueError("portable booster input has incompatible shape")
        prediction = np.full(matrix.shape[0], self.baseline, dtype=float)
        row_indexes = np.arange(matrix.shape[0])
        for offset, stop in zip(
                self.tree_offsets[:-1], self.tree_offsets[1:]):
            nodes = np.full(matrix.shape[0], offset, dtype=np.int64)
            active = np.ones(matrix.shape[0], dtype=bool)
            while np.any(active):
                active_rows = row_indexes[active]
                active_nodes = nodes[active]
                leaf = self.is_leaf[active_nodes]
                if np.any(leaf):
                    leaf_rows = active_rows[leaf]
                    prediction[leaf_rows] += (
                        self.values[active_nodes[leaf]]
                    )
                    active[leaf_rows] = False
                branch_rows = active_rows[~leaf]
                if not len(branch_rows):
                    continue
                branch_nodes = active_nodes[~leaf]
                feature_values = matrix[
                    branch_rows, self.features[branch_nodes]
                ]
                go_left = feature_values <= self.thresholds[branch_nodes]
                missing = np.isnan(feature_values)
                if np.any(missing):
                    go_left[missing] = self.missing_left[
                        branch_nodes[missing]
                    ]
                next_local = np.where(
                    go_left,
                    self.left[branch_nodes],
                    self.right[branch_nodes],
                )
                nodes[branch_rows] = offset + next_local
                if np.any(nodes[branch_rows] >= stop):
                    raise ValueError("portable booster contains invalid links")
        return prediction


@lru_cache(maxsize=1)
def load_compact_model():
    manifest = json.loads(MANIFEST_PATH.read_text())
    if (
            manifest.get("format")
            != "CRISPRidentify-portable-hist-gradient-booster"
            or manifest.get("version") != 1):
        raise ValueError("compact model manifest is incompatible")
    if (
            validate_compact_feature_schema(
                manifest.get("feature_schema", ())
            ) != FEATURE_SCHEMA):
        raise ValueError("compact model policy is incompatible")
    archive_path = MODEL_ROOT / manifest["archive"]
    if _sha256_file(archive_path) != manifest.get("archive_sha256"):
        raise ValueError("compact model archive checksum mismatch")
    archive = np.load(str(archive_path), allow_pickle=False)
    model = CompactHistGradientModel(archive)
    if model.n_features != len(FEATURE_SCHEMA):
        raise ValueError("compact model width is incompatible")
    return model


def predict_compact_scores(feature_matrix):
    scores = np.asarray(
        load_compact_model().predict(feature_matrix), dtype=float
    )
    if not np.isfinite(scores).all():
        raise ValueError("compact model returned invalid scores")
    return np.clip(scores, 0.0, 1.0)

