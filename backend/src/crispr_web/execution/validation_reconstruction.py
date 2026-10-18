"""SpacerPlacer result and generated-tree validation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path

from .files import _read_json, _safe_nonempty_output
from .validation_common import (
    _csv_finite_float,
    _csv_nonnegative_float,
    _csv_nonnegative_integer,
    _read_result_csv,
)
from .validation_newick import _tree_matches_result
from .validation_evolution import valid_metrics, valid_policy, correction, csv_correction, stage_corrections


def _validate_reconstruction_outputs(
    output_root: Path,
    expected_group_leaves: Mapping[str, frozenset[str]],
    checkpoint: Callable[[], None] | None = None,
    *, expected_bias_corrections=True, expected_rho_failure_policy="fallback_uncorrected", downstream_stage=None,
) -> bool:
    expected_groups = set(expected_group_leaves)
    result_rows = _read_result_csv(
        output_root,
        output_root / "spacerplacer" / "0_results.csv",
        required_fields={
            "name",
            "deletion_rate_bdm",
            "insertion_rate_bdm",
            "nb of spacers in alignment",
            "ln_lh_bdm",
            "Deletion model preferred by LRT",
            "nb of reconstructed insertions",
            "nb of reconstructed deletions",
            "tree length",
            "run_time",
        },
        max_rows=len(expected_groups),
        checkpoint=checkpoint,
    )
    tree_path = output_root / "spacerplacer" / "additional_data" / "dict_nwk_trees.json"
    if result_rows is None or set(result_rows) != expected_groups:
        return False
    warnings = []
    for group, row in result_rows.items():
        effective = csv_correction(row)
        if (not valid_metrics(row) or row.get("rho_bias_failure_policy") != expected_rho_failure_policy
                or not correction(effective, expected_bias_corrections, expected_rho_failure_policy)):
            return False
        if effective["rho"] == "fallback_uncorrected":
            warnings.append({"group": group, "orientation": "input_order", "code": effective["warning_code"], "message": effective["warning"]})
        if (
            any(
                not _csv_nonnegative_float(row, field)
                for field in (
                    "deletion_rate_bdm",
                    "insertion_rate_bdm",
                    "tree length",
                    "run_time",
                )
            )
            or any(
                not _csv_nonnegative_integer(row, field)
                for field in (
                    "nb of spacers in alignment",
                    "nb of reconstructed insertions",
                    "nb of reconstructed deletions",
                )
            )
            or _csv_finite_float(row, "ln_lh_bdm") is None
            or (row.get("Deletion model preferred by LRT") or "").strip()
            not in {"BDM", "IDM"}
        ):
            return False
    if downstream_stage is not None:
        policy = downstream_stage.get("bias_correction_policy")
        if (not valid_policy(policy, expected_bias_corrections, expected_rho_failure_policy, "per_group")
                or not stage_corrections(downstream_stage, policy, warnings)):
            return False
    if not _safe_nonempty_output(output_root, tree_path):
        return False
    tree_map = _read_json(output_root, tree_path, checkpoint=checkpoint)
    return set(tree_map) == expected_groups and all(
        _tree_matches_result(
            tree_map[group],
            expected_group_leaves[group],
            float(result_rows[group]["tree length"]),
        )
        for group in expected_groups
    )
