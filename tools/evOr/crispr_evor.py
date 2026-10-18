"""Dedicated, web-friendly CLI for the CRISPR-evOr orientation workflow."""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import sys
import tempfile

import spacerplacer


DESCRIPTION = (
    "CRISPR-evOr: predict CRISPR array acquisition orientation by comparing "
    "forward and reverse evolutionary reconstruction likelihoods."
)
RHO_BIAS_FALLBACK_WARNING_CODE = "rho_bias_fit_unavailable_used_uncorrected"
RHO_BIAS_FAILURE_POLICIES = frozenset({"fail", "fallback_uncorrected"})


def prepare_arguments(arguments):
    """Force orientation comparison while preserving every supported core option."""
    prepared = list(arguments)
    if "--determine_orientation" not in prepared:
        insertion_index = prepared.index("--") if "--" in prepared else len(prepared)
        prepared.insert(insertion_index, "--determine_orientation")
    return prepared


def _prediction(recommendation):
    normalized = str(recommendation).strip().lower()
    if normalized == "true":
        return "Reverse", True
    if normalized == "false":
        return "Forward", False
    if normalized == "nd":
        return "ND", None
    raise ValueError(f"Unknown orientation recommendation: {recommendation!r}")


def _finite_float(row, field):
    group = row.get("name")
    raw_value = row.get(field)
    if raw_value is None:
        raise ValueError(f"Missing {field} for group {group!r}")
    try:
        value = float(raw_value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid {field} for group {group!r}: {raw_value!r}") from error
    if not math.isfinite(value):
        raise ValueError(f"Non-finite {field} for group {group!r}")
    return value


def _validated_threshold(value):
    if isinstance(value, bool):
        raise ValueError("confidence_threshold must be finite and non-negative")
    try:
        threshold = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError("confidence_threshold must be finite and non-negative") from error
    if not math.isfinite(threshold) or threshold < 0:
        raise ValueError("confidence_threshold must be finite and non-negative")
    return threshold


def _bias_correction_result(
    row,
    *,
    prefix,
    alpha_requested,
    rho_requested,
    rho_failure_policy,
):
    alpha_status = (row.get(f"{prefix}alpha_bias_correction_status") or "").strip()
    rho_status = (row.get(f"{prefix}rho_bias_correction_status") or "").strip()
    row_policy = (row.get(f"{prefix}rho_bias_failure_policy") or "").strip()
    warning_code = (row.get(f"{prefix}rho_bias_correction_warning_code") or "").strip()
    warning = (row.get(f"{prefix}rho_bias_correction_warning") or "").strip()

    # Older detailed CSVs did not carry these columns. A successful strict run
    # implies that requested corrections were applied; a disabled correction
    # remains unambiguous from the invocation policy.
    if not alpha_status:
        alpha_status = "applied" if alpha_requested else "disabled"
    if not rho_status:
        if rho_requested and rho_failure_policy != 'fail':
            raise ValueError('Rho bias fallback policy requires explicit effective correction provenance')
        rho_status = "applied" if rho_requested else "disabled"
    if row_policy and row_policy != rho_failure_policy:
        raise ValueError(
            f"Detailed orientation result reports rho policy {row_policy!r}; "
            f"invocation requested {rho_failure_policy!r}"
        )

    valid_alpha = {"applied", "not_applicable"} if alpha_requested else {"disabled", "not_applicable"}
    valid_rho = (
        {"applied", "fallback_uncorrected", "not_applicable"}
        if rho_requested
        else {"disabled", "not_applicable"}
    )
    if alpha_status not in valid_alpha:
        raise ValueError(f"Invalid effective alpha bias-correction status: {alpha_status!r}")
    if rho_status not in valid_rho:
        raise ValueError(f"Invalid effective rho bias-correction status: {rho_status!r}")
    if rho_status == "fallback_uncorrected":
        if rho_failure_policy != "fallback_uncorrected":
            raise ValueError("Rho bias fallback occurred under a strict failure policy")
        if warning_code != RHO_BIAS_FALLBACK_WARNING_CODE or not warning:
            raise ValueError("Rho bias fallback is missing its machine-readable warning provenance")
    elif warning_code == RHO_BIAS_FALLBACK_WARNING_CODE:
        raise ValueError("Rho bias fallback warning disagrees with the effective correction status")

    return {
        "alpha": alpha_status,
        "rho": rho_status,
        "warning_code": warning_code or None,
        "warning": warning or None,
    }


def write_orientation_json(
    output_path,
    confidence_threshold,
    tree_path=None,
    *,
    alpha_bias_correction=True,
    rho_bias_correction=True,
    rho_bias_failure_policy="fail",
):
    """Create a compact, validated JSON result from the detailed comparison CSV."""
    threshold = _validated_threshold(confidence_threshold)
    if rho_bias_failure_policy not in RHO_BIAS_FAILURE_POLICIES:
        raise ValueError(
            f"rho_bias_failure_policy must be one of {sorted(RHO_BIAS_FAILURE_POLICIES)}"
        )
    output_directory = Path(output_path).expanduser().resolve()
    comparison_path = output_directory / "detailed_results" / "0_detailed_results_w_orientation.csv"
    if not comparison_path.is_file():
        raise FileNotFoundError(f"Orientation comparison output was not produced: {comparison_path}")

    required_fields = {
        "name",
        "ln_lh_bdm",
        "reversed_ln_lh_bdm",
        "ln_lh_bdm - reversed_ln_lh_bdm",
        "recommend reversing array",
    }
    results = []
    seen_groups = set()
    with comparison_path.open("r", encoding="utf-8", newline="") as comparison_file:
        reader = csv.DictReader(comparison_file)
        fieldnames = reader.fieldnames or []
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("Orientation comparison CSV contains duplicate columns")
        missing_fields = required_fields.difference(fieldnames)
        if missing_fields:
            raise ValueError(f"Orientation comparison CSV is missing fields: {sorted(missing_fields)}")

        for line_number, row in enumerate(reader, start=2):
            group = (row.get("name") or "").strip()
            if not group:
                raise ValueError(f"Orientation comparison CSV has an empty group on line {line_number}")
            if group in seen_groups:
                raise ValueError(f"Duplicate orientation result for group {group!r}")
            seen_groups.add(group)
            row["name"] = group

            prediction, recommended_reverse = _prediction(row["recommend reversing array"])
            forward_log_likelihood = _finite_float(row, "ln_lh_bdm")
            reverse_log_likelihood = _finite_float(row, "reversed_ln_lh_bdm")
            difference = _finite_float(row, "ln_lh_bdm - reversed_ln_lh_bdm")
            if not math.isclose(
                difference,
                forward_log_likelihood - reverse_log_likelihood,
                rel_tol=1e-9,
                abs_tol=1e-9,
            ):
                raise ValueError(f"Inconsistent likelihood difference for group {group!r}")
            expected_prediction = (
                "Reverse" if difference < -threshold else "Forward" if difference > threshold else "ND"
            )
            if prediction != expected_prediction:
                raise ValueError(
                    f"Orientation recommendation for group {group!r} disagrees with its likelihood difference"
                )
            bias_corrections = {
                "forward": _bias_correction_result(
                    row,
                    prefix="",
                    alpha_requested=alpha_bias_correction,
                    rho_requested=rho_bias_correction,
                    rho_failure_policy=rho_bias_failure_policy,
                ),
                "reverse": _bias_correction_result(
                    row,
                    prefix="reversed_",
                    alpha_requested=alpha_bias_correction,
                    rho_requested=rho_bias_correction,
                    rho_failure_policy=rho_bias_failure_policy,
                ),
            }
            results.append(
                {
                    "group": group,
                    "prediction": prediction,
                    "recommended_reverse": recommended_reverse,
                    "decisive": prediction != "ND",
                    "confidence_threshold": threshold,
                    "forward_ln_likelihood_bdm": forward_log_likelihood,
                    "reverse_ln_likelihood_bdm": reverse_log_likelihood,
                    "forward_minus_reverse_ln_likelihood_bdm": difference,
                    "bias_corrections": bias_corrections,
                }
            )

    orientation_trees_path = output_directory / "additional_data" / "orientation_trees.json"
    if not orientation_trees_path.is_file():
        raise FileNotFoundError(f"Orientation tree output was not produced: {orientation_trees_path}")
    try:
        with orientation_trees_path.open("r", encoding="utf-8") as tree_file:
            orientation_trees = json.load(tree_file)
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Invalid orientation tree output: {orientation_trees_path}") from error
    if not isinstance(orientation_trees, dict):
        raise ValueError("Orientation tree output must be a JSON object")
    for policy in ("forward", "reverse", "selected"):
        trees = orientation_trees.get(policy)
        if not isinstance(trees, dict):
            raise ValueError(f"Orientation tree output is missing the {policy!r} tree map")
        missing_groups = seen_groups.difference(trees)
        if missing_groups:
            raise ValueError(f"{policy.capitalize()} tree map is missing groups: {sorted(missing_groups)}")
        extra_groups = set(trees).difference(seen_groups)
        if extra_groups:
            raise ValueError(f"{policy.capitalize()} tree map has unexpected groups: {sorted(extra_groups)}")
        invalid_groups = [
            group for group in seen_groups if not isinstance(trees[group], str) or not trees[group].strip()
        ]
        if invalid_groups:
            raise ValueError(f"{policy.capitalize()} tree map has invalid trees: {sorted(invalid_groups)}")

    for result in results:
        group = result['group']
        selected_policy = 'reverse' if result['prediction'] == 'Reverse' else 'forward'
        if orientation_trees['selected'][group] != orientation_trees[selected_policy][group]:
            raise ValueError(f"Selected tree disagrees with the orientation prediction for group {group!r}")

    document = {
        "tool": "CRISPR-evOr",
        "version": spacerplacer.__version__,
        "status": "completed" if results else "completed_no_results",
        "orientation_semantics": {
            "Forward": "the input spacer order is the predicted acquisition orientation",
            "Reverse": "the reverse of the input spacer order is predicted",
            "ND": "the absolute log-likelihood difference does not exceed the confidence threshold",
        },
        "bias_correction_policy": {
            "alpha_requested": bool(alpha_bias_correction),
            "rho_requested": bool(rho_bias_correction),
            "rho_failure": rho_bias_failure_policy,
            "fallback_scope": "per_group_per_orientation",
            "orientation_likelihood_effect": (
                "none; rho fallback changes rate estimates, not optimizer likelihoods or the IDM/BDM LRT; "
                "orientation uses the BDM likelihood"
            ),
        },
        "tree_policy": "provided_shared" if tree_path is not None else "estimated_separately",
        "orientation_trees": "additional_data/orientation_trees.json",
        "results": results,
    }
    result_path = output_directory / "crispr_evor_results.json"
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=output_directory,
            prefix=".crispr_evor_results.",
            suffix=".tmp",
            delete=False,
        ) as result_file:
            temporary_path = Path(result_file.name)
            json.dump(document, result_file, indent=2, allow_nan=False)
            result_file.write(chr(10))
        os.replace(temporary_path, result_path)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return result_path


def main(argv=None):
    arguments = prepare_arguments(sys.argv[1:] if argv is None else argv)
    parsed = spacerplacer.main(arguments, prog="crispr-evor", description=DESCRIPTION)
    return write_orientation_json(
        parsed.output_path,
        parsed.orientation_decision_boundary,
        tree_path=parsed.tree_path,
        alpha_bias_correction=not parsed.no_alpha_bias_correction,
        rho_bias_correction=not parsed.no_rho_bias_correction,
        rho_bias_failure_policy=parsed.rho_bias_failure_policy,
    )


def cli():
    """Console entry point: successful runs must return process status zero."""
    main()


if __name__ == "__main__":
    cli()
