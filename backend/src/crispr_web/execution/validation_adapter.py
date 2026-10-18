"""Adapter membership and orientation-result validation."""

from __future__ import annotations

import math
from pathlib import PurePosixPath
from typing import Any, Mapping

from .validation_common import _finite_number, _nonempty_string
from .validation_newick import _NEWICK_LABEL


def _adapter_group_leaves(
    manifest: Mapping[str, Any], group_count: int
) -> dict[str, frozenset[str]] | None:
    groups = manifest.get("groups")
    if not isinstance(groups, list) or len(groups) != group_count:
        return None
    group_leaves: dict[str, frozenset[str]] = {}
    for entry in groups:
        if not isinstance(entry, Mapping) or not _nonempty_string(entry.get("name")):
            return None
        name = str(entry["name"]).strip()
        arrays = entry.get("arrays")
        array_count = entry.get("array_count")
        if (
            name in group_leaves
            or not isinstance(arrays, list)
            or isinstance(array_count, bool)
            or not isinstance(array_count, int)
            or array_count < 2
            or array_count != len(arrays)
        ):
            return None
        leaves: set[str] = set()
        for array in arrays:
            if not isinstance(array, Mapping) or not _nonempty_string(
                array.get("ccdb_file")
            ):
                return None
            relative = PurePosixPath(str(array["ccdb_file"]))
            if (
                relative.is_absolute()
                or ".." in relative.parts
                or len(relative.parts) != 4
                or relative.parts[0] != "ccdb"
                or relative.parts[1] != name
                or relative.parts[2] not in {"pos_strand", "neg_strand"}
                or relative.suffix.lower() != ".fasta"
                or _NEWICK_LABEL.fullmatch(relative.stem) is None
                or relative.stem in leaves
            ):
                return None
            leaves.add(relative.stem)
        if len(leaves) != array_count:
            return None
        group_leaves[name] = frozenset(leaves)
    return group_leaves if len(group_leaves) == group_count else None


def _adapter_group_names(
    manifest: Mapping[str, Any], group_count: int
) -> set[str] | None:
    group_leaves = _adapter_group_leaves(manifest, group_count)
    return None if group_leaves is None else set(group_leaves)


def _valid_orientation_result(entry: object) -> tuple[str, dict[str, object]] | None:
    if not isinstance(entry, Mapping) or not _nonempty_string(entry.get("group")):
        return None
    prediction = entry.get("prediction")
    if prediction not in {"Forward", "Reverse", "ND"}:
        return None
    expected_reverse = {"Forward": False, "Reverse": True, "ND": None}[prediction]
    if entry.get("recommended_reverse") is not expected_reverse:
        return None
    decisive = entry.get("decisive")
    if not isinstance(decisive, bool) or decisive != (prediction != "ND"):
        return None
    threshold = entry.get("confidence_threshold")
    forward = entry.get("forward_ln_likelihood_bdm")
    reverse = entry.get("reverse_ln_likelihood_bdm")
    difference = entry.get("forward_minus_reverse_ln_likelihood_bdm")
    if (
        not _finite_number(threshold)
        or float(threshold) < 0
        or not _finite_number(forward)
        or not _finite_number(reverse)
        or not _finite_number(difference)
    ):
        return None
    threshold_value = float(threshold)
    forward_value = float(forward)
    reverse_value = float(reverse)
    difference_value = float(difference)
    if not math.isclose(
        difference_value,
        forward_value - reverse_value,
        rel_tol=1e-9,
        abs_tol=1e-9,
    ):
        return None
    expected_prediction = (
        "Reverse"
        if difference_value < -threshold_value
        else "Forward"
        if difference_value > threshold_value
        else "ND"
    )
    if prediction != expected_prediction:
        return None
    return str(entry["group"]).strip(), {
        "prediction": prediction,
        "recommended_reverse": expected_reverse,
        "confidence_threshold": threshold_value,
        "forward_ln_likelihood_bdm": forward_value,
        "reverse_ln_likelihood_bdm": reverse_value,
        "forward_minus_reverse_ln_likelihood_bdm": difference_value,
    }

