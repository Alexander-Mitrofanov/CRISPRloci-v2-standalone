"""Deduplicated repeat folds and local support, without dense probability arrays."""

import hashlib
import json

from .common import fields, register, sources


def repeat_data(read):
    doc = read.json("repeat/repeat-analysis.json")
    folds = {}

    def fold_id(sequence, isolated):
        fold = {
            "rna_sequence": sequence,
            **fields(
                isolated,
                "status reason method structure mfe_kcal_mol ensemble_kcal_mol",
            ),
        }
        identity = (
            "fold_"
            + hashlib.sha256(
                json.dumps(fold, sort_keys=True, allow_nan=False).encode()
            ).hexdigest()[:24]
        )
        register(folds, identity, {"id": identity, **fold})
        return identity

    instances = []
    for row in doc["instances"]:
        read.step()
        instances.append(
            {
                **fields(
                    row,
                    "id source_id array_id repeat_index start end terminal orientation context_id context_start",
                ),
                "fold_id": fold_id(row["rna_sequence"], row["isolated"]),
                "comparison": fields(
                    row["comparison"],
                    "status reason motif_origin isolated_pair_support context_pair_support support_change",
                ),
            }
        )
    data = {
        "coordinate_system": "zero_based_half_open",
        "support_semantics": "descriptive-thermodynamic-model-evidence",
        "molecule": doc["molecule"],
        "orientation_policy": doc["orientation_policy"],
        "model": fields(
            doc["model"],
            "engine parameters parameter_sha256 temperature dangles salt_molar minimum_loop_size allow_GU pseudoknots no_lonely_pairs window span",
        ),
        "sources": sources(doc["sources"]),
        "classification": {
            key: fields(value, "status reason")
            for key, value in doc["classification"].items()
            if key in {"strand", "subtype", "family"}
        },
        "instances": instances,
    }
    if doc["input_kind"] == "arrays":
        data["context_semantics"] = "analysis-window-not-biological-boundary"
        data["arrays"] = [
            fields(row, "id native_id source_id start end category detector_strand")
            for row in doc["arrays"]
        ]
        data["contexts"] = [
            fields(row, "id array_id orientation status reason window span")
            for row in doc["contexts"]
        ]
        data["consensus_folds"] = [
            {
                **fields(row, "array_id orientation"),
                "fold_id": fold_id(row["rna_sequence"], row["isolated"]),
            }
            for row in doc["consensus_folds"]
        ]
    data["folds"] = list(folds.values())
    return data
