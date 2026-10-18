"""Deterministic construction and hierarchical classification of Cas cassettes."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path


DEFAULT_MAX_NON_CAS_GENES = 3
DEFAULT_MAX_CAS_GAP_NT = 20_000
DEFAULT_MIN_SUBTYPE_VOTE_FRACTION = 0.60


def _weighted_label(predictions, head):
    votes: dict[str, float] = defaultdict(float)
    supporters: dict[str, list[str]] = defaultdict(list)
    for prediction in predictions:
        label = prediction["classification"].get(head)
        if label is None:
            continue
        weight = max(1.0, float(prediction["positive_profile_score"]))
        votes[str(label)] += weight
        supporters[str(label)].append(str(prediction["sequence_id"]))
    if not votes:
        return None, 0.0, {}, {}
    ordered = sorted(votes, key=lambda label: (-votes[label], label))
    total = sum(votes.values())
    winner = ordered[0]
    return winner, votes[winner] / total, dict(sorted(votes.items())), dict(supporters)


def classify_cassette(predictions, min_subtype_vote_fraction=DEFAULT_MIN_SUBTYPE_VOTE_FRACTION):
    cas_predictions = [row for row in predictions if row["is_cas"]]
    cas_class, class_fraction, class_votes, _ = _weighted_label(cas_predictions, "class")
    cas_type, type_fraction, type_votes, _ = _weighted_label(cas_predictions, "type")
    subtype, subtype_fraction, subtype_votes, supporters = _weighted_label(
        cas_predictions, "subtype"
    )
    accepted_subtype = subtype if subtype_fraction >= min_subtype_vote_fraction else None
    if accepted_subtype is not None and cas_type is not None:
        if not accepted_subtype.startswith(f"{cas_type}-"):
            accepted_subtype = None
    return {
        "class": cas_class,
        "type": cas_type,
        "subtype": accepted_subtype,
        "raw_subtype_winner": subtype,
        "confidence": {
            "class_vote_fraction": class_fraction,
            "type_vote_fraction": type_fraction,
            "subtype_vote_fraction": subtype_fraction,
        },
        "votes": {"class": class_votes, "type": type_votes, "subtype": subtype_votes},
        "subtype_supporting_proteins": supporters.get(subtype, []) if subtype else [],
        "abstained_subtype": subtype is not None and accepted_subtype is None,
    }


def construct_cassettes(
    genes: list[dict[str, object]],
    predictions: list[dict[str, object]],
    *,
    max_non_cas_genes: int = DEFAULT_MAX_NON_CAS_GENES,
    max_cas_gap_nt: int = DEFAULT_MAX_CAS_GAP_NT,
    min_subtype_vote_fraction: float = DEFAULT_MIN_SUBTYPE_VOTE_FRACTION,
) -> list[dict[str, object]]:
    """Join nearby Cas-positive genes without using CRISPR-array evidence."""
    if max_non_cas_genes < 0 or max_cas_gap_nt < 0:
        raise ValueError("Cassette gap limits must be non-negative")
    prediction_by_id = {str(row["sequence_id"]): row for row in predictions}
    missing = sorted(str(row["protein_id"]) for row in genes if str(row["protein_id"]) not in prediction_by_id)
    if missing:
        raise ValueError(f"Missing predictions for {len(missing)} genes")
    by_contig: dict[str, list[dict[str, object]]] = defaultdict(list)
    for gene in genes:
        by_contig[str(gene["contig_id"])].append(gene)
    output = []
    for contig_id, contig_genes in sorted(by_contig.items()):
        ordered = sorted(
            contig_genes,
            key=lambda row: (int(row["start_1based"]), int(row["end_1based_inclusive"]), str(row["protein_id"])),
        )
        index_by_id = {str(row["protein_id"]): index for index, row in enumerate(ordered)}
        positives = [row for row in ordered if prediction_by_id[str(row["protein_id"])]["is_cas"]]
        groups: list[list[dict[str, object]]] = []
        for gene in positives:
            if not groups:
                groups.append([gene])
                continue
            previous = groups[-1][-1]
            intervening = index_by_id[str(gene["protein_id"])] - index_by_id[str(previous["protein_id"])] - 1
            gap_nt = int(gene["start_1based"]) - int(previous["end_1based_inclusive"]) - 1
            if intervening <= max_non_cas_genes and gap_nt <= max_cas_gap_nt:
                groups[-1].append(gene)
            else:
                groups.append([gene])
        for ordinal, group in enumerate(groups, start=1):
            first_index = index_by_id[str(group[0]["protein_id"])]
            last_index = index_by_id[str(group[-1]["protein_id"])]
            context = ordered[first_index : last_index + 1]
            group_predictions = [prediction_by_id[str(row["protein_id"])] for row in group]
            classification = classify_cassette(
                group_predictions, min_subtype_vote_fraction=min_subtype_vote_fraction
            )
            start = min(int(row["start_1based"]) for row in group)
            end = max(int(row["end_1based_inclusive"]) for row in group)
            output.append(
                {
                    "cassette_id": f"casandra|contig={contig_id}|cassette={ordinal:04d}|loc={start}-{end}",
                    "contig_id": contig_id,
                    "start_1based": start,
                    "end_1based_inclusive": end,
                    "cas_gene_count": len(group),
                    "context_gene_count": len(context),
                    "cas_protein_ids": [str(row["protein_id"]) for row in group],
                    "context_protein_ids": [str(row["protein_id"]) for row in context],
                    "classification": classification,
                    "construction": {
                        "max_non_cas_genes": max_non_cas_genes,
                        "max_cas_gap_nt": max_cas_gap_nt,
                        "crispr_array_evidence_used": False,
                    },
                }
            )
    return output


def write_cassettes(path: Path, cassettes: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for cassette in cassettes:
            json.dump(cassette, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
