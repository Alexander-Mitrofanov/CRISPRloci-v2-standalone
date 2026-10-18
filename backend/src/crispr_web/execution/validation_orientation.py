"""CRISPR-evOr metadata, likelihood, and tree-policy validation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import math
from pathlib import Path, PurePosixPath
from typing import Any

from .files import _read_json, _safe_nonempty_output
from .validation_adapter import _valid_orientation_result
from .validation_common import (
    _csv_finite_float,
    _csv_nonnegative_float,
    _csv_nonnegative_integer,
    _nonempty_string,
    _read_result_csv,
)
from .validation_newick import _tree_has_exact_leaves, _tree_matches_result
from .validation_evolution import valid_metrics, valid_policy, correction, csv_correction, stage_corrections


def _validate_orientation_outputs(
    output_root: Path,
    downstream_stage: Mapping[str, object],
    expected_group_leaves: Mapping[str, frozenset[str]],
    checkpoint: Callable[[], None] | None = None,
    *, expected_orientation_boundary: float | None = None,
    expected_bias_corrections=True, expected_rho_failure_policy="fallback_uncorrected",
) -> bool:
    expected_groups = set(expected_group_leaves)
    group_count = len(expected_groups)
    evor_root = output_root / "crispr_evor"
    result_rows = _read_result_csv(
        output_root,
        evor_root / "0_results.csv",
        required_fields={
            "name",
            "predicted orientation",
            "recommend reversing array",
            "confidence threshold",
            "ln_lh_bdm",
            "reversed_ln_lh_bdm",
            "ln_lh_bdm - reversed_ln_lh_bdm",
            "Deletion model preferred by LRT",
            "nb of reconstructed insertions",
            "tree length",
        },
        max_rows=group_count,
        checkpoint=checkpoint,
    )
    if result_rows is None or set(result_rows) != expected_groups:
        return False
    metadata_path = evor_root / "crispr_evor_results.json"
    if not _safe_nonempty_output(output_root, metadata_path):
        return False
    metadata = _read_json(output_root, metadata_path, checkpoint=checkpoint)
    results = metadata.get("results")
    semantics = metadata.get("orientation_semantics")
    if (
        metadata.get("tool") != "CRISPR-evOr"
        or not _nonempty_string(metadata.get("version"))
        or metadata.get("status") != "completed"
        or not isinstance(results, list)
        or len(results) != group_count
        or metadata.get("tree_policy") != "estimated_separately"
        or downstream_stage.get("tree_policy") != "estimated_separately"
        or not isinstance(semantics, Mapping)
        or any(
            not _nonempty_string(semantics.get(label))
            for label in ("Forward", "Reverse", "ND")
        )
    ):
        return False

    policy = metadata.get("bias_correction_policy")
    if not valid_policy(policy, expected_bias_corrections, expected_rho_failure_policy, "per_group_per_orientation"):
        return False
    warnings = []
    result_entries: dict[str, dict[str, object]] = {}
    for raw_entry in results:
        validated = _valid_orientation_result(raw_entry)
        if validated is None:
            return False
        group, entry = validated
        if expected_orientation_boundary is not None and entry.get("confidence_threshold") != expected_orientation_boundary:
            return False
        if group in result_entries:
            return False
        corrections = raw_entry.get("bias_corrections")
        if not isinstance(corrections, Mapping):
            return False
        for direction in ("forward", "reverse"):
            effective = corrections.get(direction)
            if not correction(effective, expected_bias_corrections, expected_rho_failure_policy, allow_not_applicable=True):
                return False
            if effective["rho"] == "fallback_uncorrected":
                warnings.append({"group": group, "orientation": direction, "code": effective["warning_code"], "message": effective["warning"]})
        selected = corrections["reverse" if entry["prediction"] == "Reverse" else "forward"]
        row = result_rows.get(group, {})
        if csv_correction(row) != selected or row.get("rho_bias_failure_policy") != expected_rho_failure_policy:
            return False
        result_entries[group] = entry
    if not stage_corrections(downstream_stage, policy, warnings) or set(result_entries) != expected_groups:
        return False
    if not _validate_orientation_rows(result_rows, result_entries):
        return False
    return _validate_orientation_trees(
        evor_root,
        metadata,
        result_rows,
        result_entries,
        expected_group_leaves,
        checkpoint,
    )


def _validate_orientation_rows(
    result_rows: Mapping[str, Mapping[str, str]],
    result_entries: Mapping[str, Mapping[str, object]],
) -> bool:
    for group, entry in result_entries.items():
        row = result_rows[group]
        if (not valid_metrics(row) or not _csv_nonnegative_integer(row, "nb of reconstructed deletions")
                or not _csv_nonnegative_float(row, "deletion_rate_bdm")
                or not _csv_nonnegative_float(row, "insertion_rate_bdm")):
            return False
        recommendation = (row.get("recommend reversing array") or "").strip().lower()
        expected_recommendation = {
            False: "false",
            True: "true",
            None: "nd",
        }[entry["recommended_reverse"]]
        if (row.get("predicted orientation") or "").strip() != entry[
            "prediction"
        ] or recommendation != expected_recommendation:
            return False
        if (
            (row.get("Deletion model preferred by LRT") or "").strip()
            not in {"BDM", "IDM"}
            or not _csv_nonnegative_integer(row, "nb of reconstructed insertions")
            or not _csv_nonnegative_float(row, "tree length")
        ):
            return False
        selected_likelihood = (
            "reverse_ln_likelihood_bdm"
            if entry["prediction"] == "Reverse"
            else "forward_ln_likelihood_bdm"
        )
        alternate_likelihood = (
            "forward_ln_likelihood_bdm"
            if entry["prediction"] == "Reverse"
            else "reverse_ln_likelihood_bdm"
        )
        comparisons = (
            ("confidence threshold", "confidence_threshold"),
            ("ln_lh_bdm", selected_likelihood),
            ("reversed_ln_lh_bdm", alternate_likelihood),
            (
                "ln_lh_bdm - reversed_ln_lh_bdm",
                "forward_minus_reverse_ln_likelihood_bdm",
            ),
        )
        for csv_field, metadata_field in comparisons:
            csv_value = _csv_finite_float(row, csv_field)
            if csv_value is None or not math.isclose(
                csv_value,
                float(entry[metadata_field]),
                rel_tol=1e-9,
                abs_tol=1e-9,
            ):
                return False
    return True


def _validate_orientation_trees(
    evor_root: Path,
    metadata: Mapping[str, Any],
    result_rows: Mapping[str, Mapping[str, str]],
    result_entries: Mapping[str, Mapping[str, object]],
    expected_group_leaves: Mapping[str, frozenset[str]],
    checkpoint: Callable[[], None] | None,
) -> bool:
    expected_groups = set(expected_group_leaves)
    reference = metadata.get("orientation_trees")
    if not _nonempty_string(reference):
        return False
    relative = PurePosixPath(str(reference))
    if relative.is_absolute() or ".." in relative.parts:
        return False
    tree_path = evor_root.joinpath(*relative.parts)
    if not _safe_nonempty_output(evor_root, tree_path):
        return False
    tree_document = _read_json(evor_root, tree_path, checkpoint=checkpoint)
    tree_maps: dict[str, Mapping[str, Any]] = {}
    for policy in ("forward", "reverse", "selected"):
        tree_map = tree_document.get(policy)
        if (
            not isinstance(tree_map, Mapping)
            or set(tree_map) != expected_groups
            or any(
                not _tree_has_exact_leaves(
                    tree_map[group],
                    expected_group_leaves[group],
                )
                for group in expected_groups
            )
        ):
            return False
        tree_maps[policy] = tree_map
    for group, entry in result_entries.items():
        expected_policy = "reverse" if entry["prediction"] == "Reverse" else "forward"
        selected_tree = tree_maps["selected"][group]
        tree_length = _csv_finite_float(result_rows[group], "tree length")
        if (
            selected_tree != tree_maps[expected_policy][group]
            or tree_length is None
            or not _tree_matches_result(
                selected_tree,
                expected_group_leaves[group],
                tree_length,
            )
        ):
            return False
    return True
