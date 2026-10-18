"""Self-contained production routing for cassette subtype calls."""

from __future__ import annotations

from casandra.architecture_runtime import classify_architecture_local
from casandra.classification_v4 import normalize_classification


def _normalize_final_route(raw_route):
    """Attach the raw route claim and emit its hierarchy-consistent view."""

    raw_classification = {
        key: raw_route[key] for key in ("class", "type", "subtype")
    }
    classification, hierarchy_conflict_masked = normalize_classification(
        raw_classification
    )
    return {
        **raw_route,
        **classification,
        "evidence": {
            "raw_classification": raw_classification,
            "hierarchy_conflict_masked": hierarchy_conflict_masked,
        },
    }


def classify_hybrid_cassette_v2(cassette, architecture_model_dir, genes_by_id, predictions_by_id):
    direct = cassette["classification"]
    if direct["type"] != "III":
        has_subtype = direct["subtype"] is not None
        confidence_key = (
            "subtype_vote_fraction" if has_subtype else "type_vote_fraction"
        )
        return _normalize_final_route({
            "class": direct["class"],
            "type": direct["type"],
            "subtype": direct["subtype"],
            "method": (
                "direct_profile_aggregation"
                if has_subtype else "direct_profile_type_aggregation"
            ),
            "confidence": direct["confidence"][confidence_key],
            "direct_profile_result": direct,
        })
    subtype, confidence = classify_architecture_local(
        architecture_model_dir, cassette, genes_by_id, predictions_by_id
    )
    return _normalize_final_route({
        "class": direct["class"],
        "type": direct["type"],
        "subtype": subtype,
        "method": "type_iii_ordered_profile_architecture_extratrees",
        "confidence": confidence,
        "direct_profile_result": direct,
    })
