"""Evidence aggregation and transparent, uncalibrated ranking."""

from __future__ import annotations

import math
import re
from typing import Any

from .adapters import nearest_cas
from .config import Settings
from .models import Array, Candidate, Hit, Interval, Sequence
from .sequence import subsequence
from .short_rescue_policy import SHORT_RESCUE_POLICY


def best_hit(hits: list[Hit], forward_only: bool = False) -> Hit | None:
    eligible = [h for h in hits if not forward_only or h.strand == "+"]
    return (
        min(eligible, key=lambda h: (-h.score, h.evalue, h.model, h.interval)) if eligible else None
    )


def type_v_motifs(sequence: str) -> dict[str, Any]:
    """V2 masks with corrected IUPAC classes, overlap search, and end coordinates."""
    motif = re.search("CC[CT]CC[ACGT]{6}GG[AG]GG", sequence)
    paired = None
    for first in re.finditer("(?=TTT)", sequence):
        first_end = first.start() + 3
        # Search only the allowed 12–28-base gap, including the whole second motif.
        second_start = sequence.find("CTTTC", first_end + 12, first_end + 28 + 5)
        if second_start >= 0:
            paired = {
                "TTT": Interval(first.start(), first_end).public(),
                "CTTTC": Interval(second_start, second_start + 5).public(),
                "intervening_bases": second_start - first_end,
            }
            break
    return {
        "paired_motif": paired,
        "mask_motif": None
        if motif is None
        else {"interval": Interval(*motif.span()).public(), "sequence": motif.group()},
        "coordinate_reference": "oriented_prediction_sequence",
    }


def candidate_evidence(
    candidate: Candidate,
    source: Sequence,
    models: list[Hit],
    terms: list[Hit],
    interaction: dict[str, Any] | None,
    cas: dict[str, Any],
    settings: Settings,
) -> dict[str, Any]:
    # Candidate windows are already in transcript orientation. Reverse hits are
    # retained as raw evidence but cannot support this transcript hypothesis.
    model = best_hit(models, forward_only=True)
    forward_terms = [t for t in terms if t.strand == "+"]
    term = (
        min(forward_terms, key=lambda h: (h.evalue, -h.score, h.interval))
        if forward_terms
        else None
    )
    target = (
        None
        if interaction is None
        else Interval.from_public(
            interaction["target_interval"]["start"], interaction["target_interval"]["end"]
        )
    )
    anti_tail = float(
        target is not None
        and model is not None
        and target.overlap(model.interval) < 10
        and target.start <= model.interval.start
    )
    tail_term = float(
        model is not None
        and term is not None
        and model.interval.overlap(term.interval) < 10
        and model.interval.start <= term.interval.start
    )
    expected_orientation = (
        candidate.anti.array.strand is not None
        and candidate.anti.repeat_strand == candidate.anti.array.strand
    )
    features = {
        "array": candidate.anti.array.score,
        "identity": candidate.anti.identity,
        "coverage": candidate.anti.coverage,
        "identity_coverage": candidate.anti.identity * candidate.anti.coverage,
        "interaction": max(0.0, -interaction["energy_kcal_mol"] / 10) if interaction else 0.0,
        "model": max(0.0, model.score / 10) if model else 0.0,
        "terminator": max(0.0, -math.log10(max(term.evalue, 1e-300))) if term else 0.0,
        "orientation": float(
            expected_orientation and candidate.anti.complementary and term is not None
        ),
        "anti_tail": anti_tail,
        "tail_terminator": tail_term,
    }
    contributions = {key: value * settings.weights[key] for key, value in features.items()}

    def hit_record(hit: Hit | None) -> dict[str, Any] | None:
        if hit is None:
            return None
        genomic_strand = (
            candidate.anti.strand
            if hit.strand == "+"
            else ("+" if candidate.anti.strand == "-" else "-")
        )
        return {
            **hit.public(),
            "genomic_interval": candidate.to_genome(hit.interval).public(),
            "genomic_strand": genomic_strand,
        }

    return {
        "id": candidate.id,
        "candidate_origin": "anti_repeat",
        "admission_branch": "standard",
        "evidence_class": evidence_class(
            interaction is not None, model is not None, term is not None
        ),
        "boundary_kind": "anti_repeat_search_window",
        "interval": candidate.interval.public(),
        "strand": candidate.anti.strand,
        "sequence": candidate.sequence,
        "array": {
            "id": candidate.anti.array.id,
            "interval": candidate.anti.array.interval.public(),
            "category": candidate.anti.array.category,
            "score": candidate.anti.array.score,
            "predicted_strand": candidate.anti.array.strand,
            "repeat": candidate.anti.repeat,
            "tested_repeat_strand": candidate.anti.repeat_strand,
        },
        "anti_repeat": {
            "interval": candidate.anti.interval.public(),
            "local_interval": candidate.anti_local.public(),
            "sequence": subsequence(source.bases, candidate.anti.interval, candidate.anti.strand),
            "identity": candidate.anti.identity,
            "coverage": candidate.anti.coverage,
            "complementary_hypothesis": candidate.anti.complementary,
            "array_distance_bases": candidate.anti.interval.distance(candidate.anti.array.interval),
        },
        "interaction": interaction,
        "model_hit": hit_record(model),
        "terminator": hit_record(term),
        "model_hits": [hit_record(h) for h in models],
        "terminator_hits": [hit_record(h) for h in terms],
        "poly_t_signals": [
            Interval(*m.span()).public() for m in re.finditer("T{4,}", candidate.sequence)
        ],
        "cas_context": nearest_cas(candidate.interval, cas, settings.model_type),
        "ranking": {
            "score": sum(contributions.values()),
            "features": features,
            "contributions": contributions,
            "calibrated_probability": False,
        },
        "type_v_motifs": type_v_motifs(candidate.sequence) if settings.type_v_motifs else None,
    }


def merge_model_hits(
    source: Sequence, hits: list[Hit], gap: int, motifs: bool = False
) -> list[dict[str, Any]]:
    """Merge by source/strand; preserve every contributing model and raw hit."""
    groups: list[list[Hit]] = []
    for strand in ("+", "-"):
        active: list[Hit] = []
        end = -1
        for hit in sorted((h for h in hits if h.strand == strand), key=lambda h: h.interval):
            if active and hit.interval.start - end > gap:
                groups.append(active)
                active = []
            active.append(hit)
            end = max(end, hit.interval.end) if len(active) > 1 else hit.interval.end
        if active:
            groups.append(active)
    rows = []
    for group in groups:
        interval = Interval(
            min(h.interval.start for h in group), max(h.interval.end for h in group)
        )
        sequence = subsequence(source.bases, interval, group[0].strand)
        rows.append(
            {
                "interval": interval.public(),
                "strand": group[0].strand,
                "sequence": sequence,
                "best_evalue": min(h.evalue for h in group),
                "best_score": max(h.score for h in group),
                "hits": [h.public() for h in group],
                "type_v_motifs": type_v_motifs(sequence) if motifs else None,
            }
        )
    return sorted(rows, key=lambda r: (r["interval"]["start"], r["strand"]))


def model_rescue_predictions(
    source: Sequence,
    hits: list[Hit],
    arrays: list[Array],
    existing: list[dict[str, Any]],
    cas: dict[str, Any],
    settings: Settings,
) -> list[dict[str, Any]]:
    """Retain strong genomic CM evidence even when no anti-repeat was proposed.

    These intervals are model-supported regions, not inferred full transcripts.
    Do not relax sequence thresholds or invent an array/anti-repeat association.
    """
    eligible = []
    for hit in hits:
        if hit.evalue > settings.model_rescue_evalue:
            continue
        if any(hit.interval.overlap(array.interval) for array in arrays):
            continue
        covered = any(
            p["strand"] == hit.strand
            and p["interval"]["start"] - 1 <= hit.interval.start
            and p["interval"]["end"] >= hit.interval.end
            for p in existing
        )
        if not covered:
            eligible.append(hit)
    rows = merge_model_hits(source, eligible, 0, settings.type_v_motifs)
    for index, row in enumerate(rows):
        interval = Interval.from_public(row["interval"]["start"], row["interval"]["end"])
        records = []
        for hit in row["hits"]:
            genomic = Interval.from_public(hit["interval"]["start"], hit["interval"]["end"])
            local = (
                Interval(genomic.start - interval.start, genomic.end - interval.start)
                if row["strand"] == "+"
                else Interval(interval.end - genomic.end, interval.end - genomic.start)
            )
            records.append(
                {
                    **hit,
                    "interval": local.public(),
                    "strand": "+",
                    "genomic_interval": genomic.public(),
                    "genomic_strand": row["strand"],
                }
            )
        model = min(records, key=lambda h: (-h["score"], h["evalue"]))
        features = dict.fromkeys(settings.weights, 0.0)
        features["model"] = max(0.0, model["score"] / 10)
        contributions = {k: v * settings.weights[k] for k, v in features.items()}
        row.update(
            {
                "id": f"model_rescue{index:06d}",
                "candidate_origin": "covariance_model",
                "admission_branch": "genomic_cm_rescue",
                "evidence_class": "cm_only",
                "boundary_kind": "model_supported_region",
                "array": None,
                "anti_repeat": None,
                "interaction": None,
                "terminator": None,
                "terminator_hits": [],
                "model_hit": model,
                "model_hits": records,
                "cas_context": nearest_cas(interval, cas, settings.model_type),
                "poly_t_signals": [
                    Interval(*m.span()).public() for m in re.finditer("T{4,}", row["sequence"])
                ],
                "ranking": {
                    "score": sum(contributions.values()),
                    "features": features,
                    "contributions": contributions,
                    "calibrated_probability": False,
                },
            }
        )
    return rows


def short_rescue_admission(
    candidate: Candidate,
    models: list[Hit],
    terms: list[Hit],
    interaction: dict[str, Any] | None,
    settings: Settings,
) -> dict[str, Any]:
    """Frozen experimental Type II gate, using actual antiparallel paired residues.

    Candidate-window CM E-values have a different search context from genomic
    CM rescue. These predeclared thresholds are hypotheses, not probabilities.
    """
    reasons = []
    anchored: list[int] = []
    pairs: list[dict[str, Any]] = []
    supported_models: list[Hit] = []
    supported_terms: list[Hit] = []
    if interaction is None:
        reasons.append("no_anchored_interaction")
    else:
        target = Interval.from_public(
            **{key: interaction["target_interval"][key] for key in ("start", "end")}
        )
        target_positions = interaction["target_paired_positions"]
        repeat_positions = interaction["repeat_paired_positions"]
        if interaction["energy_kcal_mol"] >= SHORT_RESCUE_POLICY.energy_strictly_below:
            reasons.append("nonnegative_energy")
        for t, q in zip(target_positions, reversed(repeat_positions), strict=True):
            bases = (candidate.sequence[t - 1] + candidate.anti.repeat[q - 1]).replace("T", "U")
            pairs.append({"target": t, "repeat": q, "bases": bases})
            if candidate.anti_local.start <= t - 1 < candidate.anti_local.end:
                anchored.append(t)
        if any(p["bases"] not in {"AU", "UA", "GC", "CG", "GU", "UG"} for p in pairs):
            reasons.append("noncanonical_or_ambiguous_pair")
        if len(anchored) < SHORT_RESCUE_POLICY.min_anchored_paired_bases:
            reasons.append("fewer_than_10_anchored_paired_bases")
        if (
            len(anchored) / candidate.anti_local.length
            < SHORT_RESCUE_POLICY.min_anchored_paired_fraction
        ):
            reasons.append("less_than_half_anchor_paired")
        supported_models = [
            h
            for h in models
            if h.strand == "+"
            and h.evalue <= SHORT_RESCUE_POLICY.cm_evalue_max
            and h.interval.start >= target.start
            and h.interval.overlap(target) < SHORT_RESCUE_POLICY.cm_duplex_overlap_strictly_below
        ]
        supported_terms = [
            h
            for h in terms
            if h.strand == "+"
            and h.evalue <= SHORT_RESCUE_POLICY.terminator_evalue_max
            and h.interval.start >= target.end
            and h.interval.end <= len(candidate.sequence)
        ]
    if not supported_models and not supported_terms:
        reasons.append("no_independent_regional_support")
    return {
        "accepted": not reasons,
        "rejection_reasons": reasons,
        "anchored_paired_positions": anchored,
        "anchored_paired_count": len(anchored),
        "anchored_paired_fraction": len(anchored) / candidate.anti_local.length,
        "antiparallel_pairs": pairs,
        "supporting_model_hits": [h.public() for h in supported_models],
        "supporting_terminator_hits": [h.public() for h in supported_terms],
        "support_search_context": "separate short-proposal candidate windows; oriented transcript",
        "model_family_responsibility": "user-supplied custom models"
        if settings.models
        else "bundled II",
        "seed_target_range": candidate.anti_local.public(),
        "thresholds": SHORT_RESCUE_POLICY.public(),
    }


def evidence_class(
    duplex: bool,
    cm: bool,
    terminator: bool,
    *,
    sequence: bool = True,
) -> str:
    if not sequence:
        return "cm_only"
    labels = ["sequence"]
    if duplex:
        labels.append("duplex")
    if cm:
        labels.append("cm")
    if terminator:
        labels.append("terminator")
    return "_".join(labels) if len(labels) > 1 else "sequence_only"
