"""Biologically routed cassette subtype classification."""

from __future__ import annotations

from casandra.architecture_runtime import classify_architecture_local


def classify_hybrid_cassette(
    cassette,
    architecture_model_dir,
    genes_by_id,
    predictions_by_id,
):
    """Use architecture learning for Type III; direct profiles elsewhere."""
    direct = cassette["classification"]
    if direct["type"] != "III":
        return {
            "class": direct["class"],
            "type": direct["type"],
            "subtype": direct["subtype"],
            "method": "direct_profile_aggregation",
            "confidence": direct["confidence"]["subtype_vote_fraction"],
            "direct_profile_result": direct,
        }
    subtype, confidence = classify_architecture_local(
        architecture_model_dir, cassette, genes_by_id, predictions_by_id
    )
    return {
        "class": direct["class"],
        "type": direct["type"],
        "subtype": subtype,
        "method": "type_iii_ordered_profile_architecture_extratrees",
        "confidence": confidence,
        "direct_profile_result": direct,
    }
