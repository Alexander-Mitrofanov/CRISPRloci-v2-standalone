"""Frozen hierarchy normalization for independent CasAndra v4 heads."""

from __future__ import annotations


CAS_TYPE_TO_CLASS = {
    "I": "1",
    "II": "2",
    "III": "1",
    "IV": "1",
    "V": "2",
    "VI": "2",
    "VII": "1",
}


def normalize_classification(raw_classification):
    """Return a hierarchy-consistent view of independent-head predictions.

    The v4 protein model scores class, type, and subtype independently. Keep
    the type winner, derive its class from the frozen Cas hierarchy, and keep
    the subtype winner only when it belongs to that type. No lower-ranked
    label is substituted for a conflicting winner.
    """

    raw_class = raw_classification["class"]
    raw_type = raw_classification["type"]
    raw_subtype = raw_classification["subtype"]
    normalized = {
        "class": (
            CAS_TYPE_TO_CLASS[raw_type]
            if raw_type is not None else raw_class
        ),
        "type": raw_type,
        "subtype": (
            raw_subtype
            if (
                raw_type is not None
                and raw_subtype is not None
                and raw_subtype.startswith(raw_type + "-")
            )
            else None
        ),
    }
    return normalized, raw_classification != normalized
