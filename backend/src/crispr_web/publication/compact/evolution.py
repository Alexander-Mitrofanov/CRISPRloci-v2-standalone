"""Comparable groups, likelihoods, rates and trees with explicit order semantics."""

import csv
import io
from pathlib import PurePosixPath

from ...summaries.result_tables import _coerce_csv
from .common import fields


METRICS = {
    "deletion_rate_bdm": "deletion_rate_bdm",
    "insertion_rate_bdm": "insertion_rate_bdm",
    "per_spacer_deletion_rate_bdm": "per_spacer_deletion_rate_bdm",
    "alpha_bdm": "alpha_bdm",
    "deletion_rate_idm": "deletion_rate_idm",
    "insertion_rate_idm": "insertion_rate_idm",
    "ln_lh_bdm": "ln_likelihood_bdm",
    "ln_lh_idm": "ln_likelihood_idm",
    "Deletion model preferred by LRT": "preferred_deletion_model",
    "test_statistic (-2*ln_lh_ratio)": "lrt_statistic",
    "significance level": "lrt_significance_level",
    "chi2_quantile": "lrt_chi2_quantile",
    "test result": "lrt_result",
    "nb of spacers in alignment": "aligned_spacers",
    "nb of unique spacers": "unique_spacers",
    "nb of reconstructed insertions": "reconstructed_insertions",
    "nb of reconstructed deletions": "reconstructed_deletions",
    "nb of reconstructed duplications": "reconstructed_duplications",
    "nb of reconstructed rearrangements": "reconstructed_rearrangements",
    "nb of reconstructed reacquisitions": "reconstructed_reacquisitions",
    "nb of reconstructed independent gains": "reconstructed_independent_gains",
    "nb of reconstructed other dup. events": "reconstructed_other_duplications",
    "nb of reconstructed dup/rearr candidates": "duplication_rearrangement_candidates",
    "nb of reconstructed ind. acquisitions, not duplicates": "independent_acquisitions",
    "avg array length": "mean_array_length",
    "tree length": "tree_length",
    "tree height": "tree_height",
    "alpha_bias_correction_status": "alpha_bias_status",
    "rho_bias_correction_status": "rho_bias_status",
    "rho_bias_correction_warning_code": "rho_warning_code",
}


def reconstruction_metrics(read, directory):
    payload = read.bytes(f"{directory}/0_results.csv").decode("utf-8")
    metrics = {}
    try:
        table = csv.DictReader(io.StringIO(payload, newline=""), strict=True)
        names = table.fieldnames or []
        if "name" not in names or len(names) != len(set(names)):
            raise ValueError("missing or duplicate reconstruction columns")
        for row in table:
            read.step()
            if None in row or len(metrics) >= 100_000:
                raise ValueError("invalid or oversized reconstruction table")
            identity = row["name"]
            if not identity or identity in metrics:
                raise ValueError("missing or duplicate reconstruction group")
            selected = {}
            for old, new in METRICS.items():
                if old not in row:
                    continue
                value = row[old]
                if value is not None and len(value) > 512:
                    raise ValueError("oversized reconstruction value")
                selected[new] = (
                    value.lower() == "true"
                    if value and value.lower() in {"true", "false"}
                    else _coerce_csv(value)
                )
            metrics[identity] = selected
    except csv.Error as error:
        raise ValueError("invalid reconstruction CSV") from error
    return metrics


def evolution_data(read, mode, status):
    manifest = read.json("adapter/manifest.json")
    groups = []
    for raw in manifest["groups"]:
        read.step()
        members = []
        for member in raw["arrays"]:
            leaf = PurePosixPath(member["ccdb_file"]).stem
            members.append(
                {
                    "leaf": leaf,
                    **fields(
                        member, "source_id array_id strand input_sequence_orientation"
                    ),
                }
            )
        groups.append({"id": raw["name"], "members": members})
    data = {
        "order_semantics": "array-acquisition-order",
        "history_semantics": "model-estimated",
        "groups": groups,
        "excluded_arrays": [
            fields(row, "source_id array_id reason")
            for row in manifest["skipped_arrays"]
        ],
    }
    if status == "completed_no_eligible_groups":
        if groups:
            raise ValueError("no-group result has groups")
        return data

    directory = "crispr_evor" if mode == "orientation" else "spacerplacer"
    metrics = reconstruction_metrics(read, directory)
    if len(groups) != len(metrics) or set(metrics) != {group["id"] for group in groups}:
        raise ValueError("reconstruction group mismatch")
    if mode == "reconstruction":
        trees = read.json("spacerplacer/additional_data/dict_nwk_trees.json")
        if set(trees) != set(metrics):
            raise ValueError("reconstruction tree group mismatch")
        for group in groups:
            group.update(
                hypothesis="input",
                tree=trees[group["id"]],
                reconstruction=metrics[group["id"]],
            )
        return data

    doc = read.json("crispr_evor/crispr_evor_results.json")
    comparisons = {row["group"]: row for row in doc["results"]}
    trees = read.json("crispr_evor/additional_data/orientation_trees.json")
    if (
        set(comparisons) != set(metrics)
        or len(comparisons) != len(doc["results"])
        or any(
            set(trees[policy]) != set(metrics)
            for policy in ("forward", "reverse", "selected")
        )
    ):
        raise ValueError("orientation group mismatch")
    data["tree_policy"] = doc["tree_policy"]
    data["bias_policy"] = fields(
        doc["bias_correction_policy"],
        "alpha_requested rho_requested rho_failure fallback_scope",
    )
    for group in groups:
        row = comparisons[group["id"]]
        prediction = row["prediction"]
        if prediction not in {"Forward", "Reverse", "ND"}:
            raise ValueError("unsupported orientation prediction")
        group["orientation"] = {
            "prediction": prediction,
            "threshold": row["confidence_threshold"],
            "input_ln_likelihood_bdm": row["forward_ln_likelihood_bdm"],
            "reverse_ln_likelihood_bdm": row["reverse_ln_likelihood_bdm"],
            "input_minus_reverse_ln_likelihood_bdm": row[
                "forward_minus_reverse_ln_likelihood_bdm"
            ],
        }
        group["selected_hypothesis"] = "reverse" if prediction == "Reverse" else "input"
        group["selection_basis"] = (
            "unresolved_default" if prediction == "ND" else "supported"
        )
        group["trees"] = {
            "input": trees["forward"][group["id"]],
            "reverse": trees["reverse"][group["id"]],
        }
        group["reconstruction"] = metrics[group["id"]]
        if "bias_corrections" in row:
            group["bias_corrections"] = {
                direction: fields(
                    row["bias_corrections"][direction], "alpha rho warning_code"
                )
                for direction in ("forward", "reverse")
            }
    return data
