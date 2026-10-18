"""Biologically routed cassette subtype classification."""

from __future__ import annotations

from casandra.architecture_runtime import classify_architecture_local

ARCHITECTURE_METHODS = {
    "II": "type_ii_ordered_profile_architecture_extratrees",
    "III": "type_iii_ordered_profile_architecture_extratrees",
}


def classify_hybrid_cassette(cassette, architecture_model, genes_by_id, predictions_by_id):
    """Use ordered architecture evidence for Types II/III and direct profiles elsewhere."""
    direct = cassette["classification"]
    method = ARCHITECTURE_METHODS.get(direct["type"])
    if method is None:
        return {
            "class": direct["class"],
            "type": direct["type"],
            "subtype": direct["subtype"],
            "method": "direct_profile_aggregation",
            "confidence": direct["confidence"]["subtype_vote_fraction"],
            "direct_profile_result": direct,
        }
    subtype, confidence = classify_architecture_local(
        architecture_model, cassette, genes_by_id, predictions_by_id
    )
    return {
        "class": direct["class"],
        "type": direct["type"],
        "subtype": subtype,
        "method": method,
        "confidence": confidence,
        "direct_profile_result": direct,
    }
