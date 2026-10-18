"""Pure source-backed planning for periodic CRISPR lattice repair.

This module is intentionally not wired into the production pipeline.  It
answers one narrow question: can an already accepted repeat core be expanded
to a complete, source-supported periodic lattice without using its classifier
score, legacy consensus, or comparator coordinates?

The planner is fail-closed:

* every input and output interval is replayed against the original source;
* period, phase, internal holes, and terminal occurrences are resolved in one
  atomic search rather than by sequential weak-repeat insertions;
* all near-optimal phase alternatives are retained;
* occurrence count and repeat-start locations must be identical across those
  alternatives, while exact repeat ends may remain explicitly ambiguous;
* coding-context and coding-tandem assessments are mandatory external
  authorization inputs;
* configured search caps cause abstention and never return a partial result.

Coordinates are zero-based and end-exclusive.
"""

from dataclasses import dataclass
import hashlib
import math
from typing import Optional, Tuple

from components.components_evaluated_arrays_enhancement import (
    _longest_exact_source_anchor,
    _source_phase_boundary_evidence,
    _source_edit_identity,
    _source_only_alignment_seed,
    _source_phase_geometry_score,
    _source_spacer_statistics,
    build_source_backed_candidate_from_boundaries,
    validate_candidate_source_reconstruction,
)


COORDINATE_SYSTEM = "zero_based_half_open"
RULE_NAME = "source-periodic-lattice-phase-then-terminal"
RULE_VERSION = "0.1.0-isolated"


def source_sequence_sha256(sequence):
    """Return the canonical uppercase source-sequence digest."""

    if not isinstance(sequence, str):
        raise ValueError("source DNA must be a string")
    try:
        encoded = sequence.upper().encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("source DNA must contain ASCII characters") from error
    return hashlib.sha256(encoded).hexdigest()


def _slice_sha256(source, interval):
    return source_sequence_sha256(source[interval[0]:interval[1]])


@dataclass(frozen=True)
class SourcePeriodicLatticeRule:
    """Frozen conservative limits for the isolated proposal search."""

    minimum_seed_repeats: int = 8
    minimum_supported_lattice_edges: int = 7
    minimum_repeat_length: int = 20
    maximum_repeat_length: int = 55
    minimum_spacer_length: int = 18
    maximum_spacer_length: int = 80
    maximum_start_lattice_residual: int = 2
    maximum_internal_missing_slots: int = 8
    maximum_boundary_shift: int = 12
    maximum_phase_hypotheses: int = 625
    maximum_terminal_extensions_per_side: int = 4
    maximum_total_repeat_count: int = 64
    minimum_phase_score_improvement: float = 30.0
    plausible_phase_score_delta: float = 3.0
    minimum_repeat_identity: float = 0.70
    strong_repeat_identity: float = 0.85
    minimum_mean_repeat_identity: float = 0.89
    minimum_representative_mean_repeat_identity: float = 0.90
    minimum_strong_repeat_count: int = 8
    maximum_weak_repeat_fraction: float = 1.0 / 3.0
    minimum_exact_anchor_length: int = 6
    maximum_spacer_length_cv: float = 0.05
    maximum_distinct_spacer_identity: float = 0.80
    terminal_kmer_length: int = 5
    maximum_terminal_kmer_dominance: float = 0.80
    minimum_parent_interval_coverage: float = 1.0
    minimum_column_coverage_fraction: float = 0.80
    spacer_edge_width: int = 3
    background_flank: int = 5000


SOURCE_PERIODIC_LATTICE_RULE = SourcePeriodicLatticeRule()


@dataclass(frozen=True)
class SourcePeriodicLatticeSafetyEvidence:
    """External fail-closed authorization for source-only geometry repair."""

    source_sha256: str
    coding_context_assessed: bool
    tandem_evidence_assessed: bool
    coding_tandem_supported: bool
    tandem_decision_reason: str


@dataclass(frozen=True)
class LatticePeriodFit:
    period: int
    observed_slots: Tuple[int, ...]
    signed_residuals: Tuple[int, ...]
    maximum_absolute_residual: int
    total_absolute_residual: int
    internal_missing_slot_count: int


@dataclass(frozen=True)
class SourcePeriodicLatticeGeometryPrefilterResult:
    """Cheap source-geometry screen that carries no safety authorization."""

    eligible: bool
    reason: str
    source_sha256: str
    parent_repeat_intervals: Tuple[Tuple[int, int], ...]
    period_fits: Tuple[LatticePeriodFit, ...]
    selected_period: Optional[int]
    observed_lattice_slots: Tuple[int, ...]
    phase_hypothesis_capacity: int
    cap_exhausted: bool


@dataclass(frozen=True)
class SourcePeriodicLatticePrefilterResult:
    eligible: bool
    reason: str
    source_sha256: str
    parent_repeat_intervals: Tuple[Tuple[int, int], ...]
    period_fits: Tuple[LatticePeriodFit, ...]
    selected_period: Optional[int]
    observed_lattice_slots: Tuple[int, ...]
    phase_hypothesis_capacity: int
    cap_exhausted: bool


@dataclass(frozen=True)
class ParentOccurrenceMapping:
    parent_index: int
    parent_interval: Tuple[int, int]
    lattice_slot: int
    proposed_interval: Tuple[int, int]
    parent_interval_coverage: float
    parent_source_sha256: str
    proposed_source_sha256: str


@dataclass(frozen=True)
class FlankingLatticeStopEvidence:
    side: str
    lattice_slot: int
    interval: Optional[Tuple[int, int]]
    status: str
    repeat_identity: Optional[float]
    exact_anchor_length: Optional[int]
    source_sha256: Optional[str]
    gate_reasons: Tuple[str, ...]


@dataclass(frozen=True)
class PhaseBoundarySideEvidence:
    varies: bool
    passes: bool
    included_dominance: float
    immediately_excluded_dominance: float
    dominance_improvement: float
    included_discordant_canonical_count: int
    evidence_branch: Optional[str]


@dataclass(frozen=True)
class PhaseBoundaryEvidence:
    left: PhaseBoundarySideEvidence
    right: PhaseBoundarySideEvidence


@dataclass(frozen=True)
class SourcePeriodicLatticeAlternative:
    """One complete and independently source-replayable phase alternative."""

    left_shift: int
    right_shift: int
    period: int
    repeat_length: int
    phase_source_score: float
    phase_score_improvement: float
    phase_score_delta_from_best: float
    lattice_slots: Tuple[int, ...]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    repeat_source_sha256s: Tuple[str, ...]
    spacer_intervals: Tuple[Tuple[int, int], ...]
    spacer_source_sha256s: Tuple[str, ...]
    envelope: Tuple[int, int]
    envelope_source_sha256: str
    source_medoid: str
    source_medoid_sha256: str
    repeat_identities: Tuple[float, ...]
    exact_anchor_lengths: Tuple[int, ...]
    mean_repeat_identity: float
    strong_repeat_count: int
    weak_repeat_slots: Tuple[int, ...]
    spacer_lengths: Tuple[int, ...]
    spacer_length_cv: float
    distinct_spacer_count: int
    maximum_distinct_spacer_identity: float
    terminal_kmer_dominance: float
    internal_missing_slots: Tuple[int, ...]
    terminal_added_slots: Tuple[int, ...]
    parent_mappings: Tuple[ParentOccurrenceMapping, ...]
    flanking_stops: Tuple[FlankingLatticeStopEvidence, ...]


@dataclass(frozen=True)
class SourcePeriodicLatticeProposalResult:
    """Immutable output of one exhaustive bounded proposal search."""

    status: str
    reason: str
    coordinate_system: str
    rule_name: str
    rule_version: str
    rule: SourcePeriodicLatticeRule
    source_sha256: str
    safety_evidence: SourcePeriodicLatticeSafetyEvidence
    parent_repeat_intervals: Tuple[Tuple[int, int], ...]
    parent_repeat_source_sha256s: Tuple[str, ...]
    period_fits: Tuple[LatticePeriodFit, ...]
    selected_period: Optional[int]
    observed_lattice_slots: Tuple[int, ...]
    phase_hypothesis_capacity: int
    evaluated_phase_hypothesis_count: int
    geometrically_valid_phase_hypothesis_count: int
    structurally_evaluated_phase_hypothesis_count: int
    structurally_valid_phase_hypothesis_count: int
    structurally_pruned_phase_hypothesis_count: int
    plausible_phase_alternative_count: int
    background_interval: Optional[Tuple[int, int]]
    background_source_sha256: Optional[str]
    background_probabilities: Tuple[Tuple[str, float], ...]
    baseline_phase_source_score: Optional[float]
    best_phase_source_score: Optional[float]
    phase_boundary_evidence: Optional[PhaseBoundaryEvidence]
    alternatives: Tuple[SourcePeriodicLatticeAlternative, ...]
    representative_index: Optional[int]
    repeat_count_invariant: bool
    repeat_start_locations_invariant: bool
    invariant_repeat_starts: Tuple[int, ...]
    exact_boundary_status: str
    classification_ceiling: Optional[str]
    cap_exhausted: bool
    gate_reason_counts: Tuple[Tuple[str, int], ...]

    @property
    def representative(self):
        if self.representative_index is None:
            return None
        return self.alternatives[self.representative_index]


@dataclass(frozen=True)
class _GeometryMetrics:
    medoid: str
    repeat_identities: Tuple[float, ...]
    exact_anchor_lengths: Tuple[int, ...]
    mean_repeat_identity: float
    strong_repeat_count: int
    weak_indexes: Tuple[int, ...]
    spacers: Tuple[str, ...]
    spacer_lengths: Tuple[int, ...]
    spacer_length_cv: float
    distinct_spacer_count: int
    maximum_distinct_spacer_identity: float
    terminal_kmer_dominance: float
    gate_reasons: Tuple[str, ...]


@dataclass(frozen=True)
class _PhaseDraft:
    left_shift: int
    right_shift: int
    lattice_slots: Tuple[int, ...]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    parent_mappings: Tuple[ParentOccurrenceMapping, ...]
    metrics: Optional[_GeometryMetrics]
    source_score: Optional[float] = None


@dataclass(frozen=True)
class _ExtendedDraft:
    phase: _PhaseDraft
    lattice_slots: Tuple[int, ...]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    metrics: _GeometryMetrics
    terminal_added_slots: Tuple[int, ...]
    flanking_stops: Tuple[FlankingLatticeStopEvidence, ...]
    cap_exhausted: bool


def _validate_rule(rule):
    if not isinstance(rule, SourcePeriodicLatticeRule):
        raise TypeError("rule must be SourcePeriodicLatticeRule")
    integer_fields = (
        "minimum_seed_repeats",
        "minimum_supported_lattice_edges",
        "minimum_repeat_length",
        "maximum_repeat_length",
        "minimum_spacer_length",
        "maximum_spacer_length",
        "maximum_start_lattice_residual",
        "maximum_internal_missing_slots",
        "maximum_boundary_shift",
        "maximum_phase_hypotheses",
        "maximum_terminal_extensions_per_side",
        "maximum_total_repeat_count",
        "minimum_strong_repeat_count",
        "minimum_exact_anchor_length",
        "terminal_kmer_length",
        "spacer_edge_width",
        "background_flank",
    )
    for field in integer_fields:
        value = getattr(rule, field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise ValueError("{} must be a non-negative integer".format(field))
    if rule.minimum_seed_repeats < 3:
        raise ValueError("minimum_seed_repeats must be at least three")
    if rule.minimum_supported_lattice_edges < 2:
        raise ValueError(
            "minimum_supported_lattice_edges must be at least two"
        )
    if rule.minimum_repeat_length < 1:
        raise ValueError("minimum_repeat_length must be positive")
    if rule.maximum_repeat_length < rule.minimum_repeat_length:
        raise ValueError("repeat-length bounds are inconsistent")
    if rule.maximum_spacer_length < rule.minimum_spacer_length:
        raise ValueError("spacer-length bounds are inconsistent")
    if rule.maximum_phase_hypotheses < 1:
        raise ValueError("maximum_phase_hypotheses must be positive")
    if rule.maximum_total_repeat_count < rule.minimum_seed_repeats:
        raise ValueError("maximum_total_repeat_count is too small")
    if rule.minimum_strong_repeat_count > rule.maximum_total_repeat_count:
        raise ValueError("minimum_strong_repeat_count is too large")
    if rule.terminal_kmer_length < 1:
        raise ValueError("terminal_kmer_length must be positive")
    if rule.spacer_edge_width < 1:
        raise ValueError("spacer_edge_width must be positive")

    probability_fields = (
        "minimum_repeat_identity",
        "strong_repeat_identity",
        "minimum_mean_repeat_identity",
        "minimum_representative_mean_repeat_identity",
        "maximum_weak_repeat_fraction",
        "maximum_distinct_spacer_identity",
        "maximum_terminal_kmer_dominance",
        "minimum_parent_interval_coverage",
        "minimum_column_coverage_fraction",
    )
    for field in probability_fields:
        value = getattr(rule, field)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not 0.0 <= value <= 1.0
        ):
            raise ValueError("{} must be finite and in [0, 1]".format(field))
    if rule.strong_repeat_identity < rule.minimum_repeat_identity:
        raise ValueError("strong repeat identity is below the weak minimum")
    if not 0.0 < rule.minimum_column_coverage_fraction <= 1.0:
        raise ValueError(
            "minimum_column_coverage_fraction must be in (0, 1]"
        )
    numeric_fields = (
        "minimum_phase_score_improvement",
        "plausible_phase_score_delta",
        "maximum_spacer_length_cv",
    )
    for field in numeric_fields:
        value = getattr(rule, field)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0.0
        ):
            raise ValueError("{} must be finite and non-negative".format(field))


def _validate_safety_evidence(source, evidence):
    if not isinstance(evidence, SourcePeriodicLatticeSafetyEvidence):
        raise TypeError(
            "safety_evidence must be SourcePeriodicLatticeSafetyEvidence"
        )
    observed_hash = source_sequence_sha256(source)
    if evidence.source_sha256 != observed_hash:
        raise ValueError("safety evidence source hash does not match source DNA")
    if not isinstance(evidence.tandem_decision_reason, str):
        raise ValueError("tandem_decision_reason must be a string")
    for field in (
        "coding_context_assessed",
        "tandem_evidence_assessed",
        "coding_tandem_supported",
    ):
        if not isinstance(getattr(evidence, field), bool):
            raise ValueError("{} must be boolean".format(field))


def _parent_intervals(candidate):
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts,
            candidate.list_repeats,
        )
    )


def _period_bounds(rule):
    return (
        rule.minimum_repeat_length + rule.minimum_spacer_length,
        rule.maximum_repeat_length + rule.maximum_spacer_length,
    )


def _infer_periods(starts, rule):
    minimum_period, maximum_period = _period_bounds(rule)
    fits = []
    anchor = starts[0]
    for period in range(minimum_period, maximum_period + 1):
        slots = []
        residuals = []
        previous_slot = -1
        for start in starts:
            delta = start - anchor
            slot = int(math.floor(delta / float(period) + 0.5))
            residual = delta - slot * period
            if (
                slot <= previous_slot
                or abs(residual) > rule.maximum_start_lattice_residual
            ):
                break
            slots.append(slot)
            residuals.append(residual)
            previous_slot = slot
        else:
            missing = slots[-1] + 1 - len(slots)
            if (
                len(starts) - 1
                >= rule.minimum_supported_lattice_edges
                and missing <= rule.maximum_internal_missing_slots
                and slots[-1] + 1 <= rule.maximum_total_repeat_count
            ):
                fits.append(LatticePeriodFit(
                    period=period,
                    observed_slots=tuple(slots),
                    signed_residuals=tuple(residuals),
                    maximum_absolute_residual=max(map(abs, residuals)),
                    total_absolute_residual=sum(map(abs, residuals)),
                    internal_missing_slot_count=missing,
                ))
    return tuple(sorted(
        fits,
        key=lambda fit: (
            fit.maximum_absolute_residual,
            fit.total_absolute_residual,
            fit.internal_missing_slot_count,
            fit.period,
        ),
    ))


def source_periodic_lattice_geometry_prefilter(
    full_dna,
    crispr_candidate,
    rule=SOURCE_PERIODIC_LATTICE_RULE,
):
    """Reject geometry that cannot enter the bounded lattice search.

    This screen accepts no safety evidence and cannot authorize a proposal.
    It lets callers skip coding-context work for geometrically ineligible
    accepted rows without fabricating an assessed tandem decision.
    """

    _validate_rule(rule)
    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    validate_candidate_source_reconstruction(
        full_dna,
        crispr_candidate.list_repeats,
        crispr_candidate.list_repeat_starts,
        crispr_candidate.list_spacers,
    )
    parent = _parent_intervals(crispr_candidate)
    source_hash = source_sequence_sha256(full_dna)
    capacity = (2 * rule.maximum_boundary_shift + 1) ** 2

    def result(eligible, reason, period_fits=(), cap_exhausted=False):
        selected = period_fits[0] if len(period_fits) == 1 else None
        return SourcePeriodicLatticeGeometryPrefilterResult(
            eligible=eligible,
            reason=reason,
            source_sha256=source_hash,
            parent_repeat_intervals=parent,
            period_fits=tuple(period_fits),
            selected_period=(selected.period if selected else None),
            observed_lattice_slots=(
                selected.observed_slots if selected else ()
            ),
            phase_hypothesis_capacity=capacity,
            cap_exhausted=cap_exhausted,
        )

    if len(parent) < rule.minimum_seed_repeats:
        return result(False, "insufficient_seed_repeat_count")

    period_fits = _infer_periods(
        tuple(start for start, _end in parent), rule
    )
    if not period_fits:
        return result(False, "period_not_supported")
    if len(period_fits) != 1:
        return result(False, "period_not_unique", period_fits)
    if capacity > rule.maximum_phase_hypotheses:
        return result(
            False,
            "phase_hypothesis_cap_exceeded",
            period_fits,
            cap_exhausted=True,
        )
    return result(True, "eligible", period_fits)


def source_periodic_lattice_prefilter(
    full_dna,
    crispr_candidate,
    safety_evidence,
    rule=SOURCE_PERIODIC_LATTICE_RULE,
):
    """Cheaply reject calls that cannot enter the authorized lattice search."""

    _validate_safety_evidence(full_dna, safety_evidence)
    geometry = source_periodic_lattice_geometry_prefilter(
        full_dna,
        crispr_candidate,
        rule=rule,
    )

    def result(eligible, reason):
        return SourcePeriodicLatticePrefilterResult(
            eligible=eligible,
            reason=reason,
            source_sha256=geometry.source_sha256,
            parent_repeat_intervals=geometry.parent_repeat_intervals,
            period_fits=geometry.period_fits,
            selected_period=geometry.selected_period,
            observed_lattice_slots=geometry.observed_lattice_slots,
            phase_hypothesis_capacity=geometry.phase_hypothesis_capacity,
            cap_exhausted=geometry.cap_exhausted,
        )

    if not safety_evidence.coding_context_assessed:
        return result(False, "coding_context_unassessed")
    if not safety_evidence.tandem_evidence_assessed:
        return result(False, "tandem_evidence_unassessed")
    if safety_evidence.coding_tandem_supported:
        return result(False, "coding_tandem_supported")
    return result(geometry.eligible, geometry.reason)


def _interval_overlap_fraction(inner, outer):
    overlap = max(
        0,
        min(inner[1], outer[1]) - max(inner[0], outer[0]),
    )
    return overlap / float(inner[1] - inner[0])


def _score_rule(rule):
    return {
        "maximum_boundary_shift": rule.maximum_boundary_shift,
        "minimum_repeats": rule.minimum_seed_repeats,
        "minimum_repeat_length": rule.minimum_repeat_length,
        "maximum_repeat_length": rule.maximum_repeat_length,
        "minimum_mean_repeat_length": rule.minimum_repeat_length,
        "maximum_mean_repeat_length": rule.maximum_repeat_length,
        "minimum_spacer_length": rule.minimum_spacer_length,
        "maximum_spacer_length": rule.maximum_spacer_length,
        "minimum_column_coverage_fraction": (
            rule.minimum_column_coverage_fraction
        ),
        "spacer_edge_width": rule.spacer_edge_width,
        "background_flank": rule.background_flank,
        "minimum_score_improvement": rule.minimum_phase_score_improvement,
        "minimum_winner_margin": rule.plausible_phase_score_delta,
    }


def _geometry_metrics(source, intervals, rule):
    repeats = tuple(source[start:end] for start, end in intervals)
    medoid = _source_only_alignment_seed(repeats)
    identities = tuple(
        _source_edit_identity(repeat, medoid) for repeat in repeats
    )
    anchors = tuple(
        _longest_exact_source_anchor(repeat, medoid) for repeat in repeats
    )
    mean_identity = sum(identities) / float(len(identities))
    weak_indexes = tuple(
        index for index, identity in enumerate(identities)
        if identity < rule.strong_repeat_identity
    )
    strong_count = len(identities) - len(weak_indexes)
    spacers = tuple(
        source[left_end:right_start]
        for (_left_start, left_end), (right_start, _right_end)
        in zip(intervals, intervals[1:])
    )
    spacer_stats = _source_spacer_statistics(
        spacers, rule.terminal_kmer_length
    )
    reasons = []
    repeat_lengths = tuple(end - start for start, end in intervals)
    if not all(
        rule.minimum_repeat_length <= length <= rule.maximum_repeat_length
        for length in repeat_lengths
    ):
        reasons.append("repeat_length_out_of_bounds")
    if min(identities) < rule.minimum_repeat_identity:
        reasons.append("repeat_identity_below_minimum")
    if min(anchors) < rule.minimum_exact_anchor_length:
        reasons.append("repeat_anchor_below_minimum")
    if mean_identity < rule.minimum_mean_repeat_identity:
        reasons.append("mean_repeat_identity_below_minimum")
    if strong_count < rule.minimum_strong_repeat_count:
        reasons.append("strong_repeat_count_below_minimum")
    if len(weak_indexes) / float(len(identities)) > (
        rule.maximum_weak_repeat_fraction
    ):
        reasons.append("weak_repeat_fraction_above_maximum")
    if not spacer_stats["lengths"]:
        reasons.append("spacers_missing")
    elif not all(
        rule.minimum_spacer_length <= length <= rule.maximum_spacer_length
        for length in spacer_stats["lengths"]
    ):
        reasons.append("spacer_length_out_of_bounds")
    if spacer_stats["cv"] > rule.maximum_spacer_length_cv:
        reasons.append("spacer_length_cv_above_maximum")
    if spacer_stats["distinct_count"] != len(spacers):
        reasons.append("duplicate_spacers")
    if (
        spacer_stats["maximum_identity"]
        > rule.maximum_distinct_spacer_identity
    ):
        reasons.append("spacer_identity_above_maximum")
    if (
        spacer_stats["terminal_kmer_dominance"]
        >= rule.maximum_terminal_kmer_dominance
    ):
        reasons.append("terminal_kmer_dominance_above_maximum")
    return _GeometryMetrics(
        medoid=medoid,
        repeat_identities=identities,
        exact_anchor_lengths=anchors,
        mean_repeat_identity=mean_identity,
        strong_repeat_count=strong_count,
        weak_indexes=weak_indexes,
        spacers=spacers,
        spacer_lengths=tuple(spacer_stats["lengths"]),
        spacer_length_cv=spacer_stats["cv"],
        distinct_spacer_count=spacer_stats["distinct_count"],
        maximum_distinct_spacer_identity=spacer_stats["maximum_identity"],
        terminal_kmer_dominance=spacer_stats[
            "terminal_kmer_dominance"
        ],
        gate_reasons=tuple(reasons),
    )


def _background(source, geometries, flank):
    starts = tuple(
        start for geometry in geometries for start, _end in geometry
    )
    ends = tuple(
        end for geometry in geometries for _start, end in geometry
    )
    interval = (
        max(0, min(starts) - flank),
        min(len(source), max(ends) + flank),
    )
    counts = {
        base: source[interval[0]:interval[1]].upper().count(base)
        for base in "ACGT"
    }
    size = sum(counts.values())
    probabilities = {
        base: (counts[base] + 1.0) / (size + 4.0)
        for base in "ACGT"
    }
    return interval, probabilities


def _phase_drafts(source, parent, period_fit, rule):
    phase_capacity = (2 * rule.maximum_boundary_shift + 1) ** 2
    if phase_capacity > rule.maximum_phase_hypotheses:
        return (), phase_capacity, True

    observed_slots = period_fit.observed_slots
    first_parent_start, first_parent_end = parent[0]
    all_slots = tuple(range(observed_slots[-1] + 1))
    drafts = []
    for left_shift in range(
        -rule.maximum_boundary_shift, rule.maximum_boundary_shift + 1
    ):
        lattice_anchor = first_parent_start + left_shift
        for right_shift in range(
            -rule.maximum_boundary_shift, rule.maximum_boundary_shift + 1
        ):
            repeat_end = first_parent_end + right_shift
            repeat_length = repeat_end - lattice_anchor
            if not (
                rule.minimum_repeat_length
                <= repeat_length
                <= rule.maximum_repeat_length
            ):
                continue
            intervals = tuple(
                (
                    lattice_anchor + slot * period_fit.period,
                    lattice_anchor + slot * period_fit.period + repeat_length,
                )
                for slot in all_slots
            )
            if (
                intervals[0][0] < 0
                or intervals[-1][1] > len(source)
                or any(
                    left[1] > right[0]
                    for left, right in zip(intervals, intervals[1:])
                )
            ):
                continue

            mappings = []
            for parent_index, (parent_interval, lattice_slot) in enumerate(
                zip(parent, observed_slots)
            ):
                proposed_interval = intervals[lattice_slot]
                coverage = _interval_overlap_fraction(
                    parent_interval, proposed_interval
                )
                if coverage < rule.minimum_parent_interval_coverage:
                    break
                mappings.append(ParentOccurrenceMapping(
                    parent_index=parent_index,
                    parent_interval=parent_interval,
                    lattice_slot=lattice_slot,
                    proposed_interval=proposed_interval,
                    parent_interval_coverage=coverage,
                    parent_source_sha256=_slice_sha256(
                        source, parent_interval
                    ),
                    proposed_source_sha256=_slice_sha256(
                        source, proposed_interval
                    ),
                ))
            else:
                drafts.append(_PhaseDraft(
                    left_shift=left_shift,
                    right_shift=right_shift,
                    lattice_slots=all_slots,
                    repeat_intervals=intervals,
                    parent_mappings=tuple(mappings),
                    metrics=None,
                ))
    return tuple(drafts), phase_capacity, False


def _with_scores(source, drafts, background_probabilities, rule):
    score_rule = _score_rule(rule)
    scored = []
    for draft in drafts:
        intervals = draft.repeat_intervals
        score = _source_phase_geometry_score(
            source,
            tuple(start for start, _end in intervals),
            tuple(end for _start, end in intervals),
            score_rule,
            background_probabilities,
        )
        scored.append(_PhaseDraft(
            left_shift=draft.left_shift,
            right_shift=draft.right_shift,
            lattice_slots=draft.lattice_slots,
            repeat_intervals=draft.repeat_intervals,
            parent_mappings=draft.parent_mappings,
            metrics=draft.metrics,
            source_score=score,
        ))
    return tuple(scored)


def _terminal_gate(source, interval, medoid, rule):
    if interval[0] < 0 or interval[1] > len(source):
        return None, None, ("source_boundary",)
    sequence = source[interval[0]:interval[1]]
    identity = _source_edit_identity(sequence, medoid)
    anchor = _longest_exact_source_anchor(sequence, medoid)
    reasons = []
    if identity < rule.minimum_repeat_identity:
        reasons.append("repeat_identity_below_minimum")
    if anchor < rule.minimum_exact_anchor_length:
        reasons.append("repeat_anchor_below_minimum")
    return identity, anchor, tuple(reasons)


def _extend_terminals(source, draft, period, rule):
    slot_intervals = dict(zip(draft.lattice_slots, draft.repeat_intervals))
    initial_slots = set(slot_intervals)
    flanking_stops = []
    cap_exhausted = False
    repeat_length = (
        draft.repeat_intervals[0][1] - draft.repeat_intervals[0][0]
    )
    anchor = (
        draft.repeat_intervals[0][0] - draft.lattice_slots[0] * period
    )

    for side in ("left", "right"):
        added = 0
        while True:
            next_slot = (
                min(slot_intervals) - 1
                if side == "left"
                else max(slot_intervals) + 1
            )
            start = anchor + next_slot * period
            interval = (start, start + repeat_length)
            identity, exact_anchor, reasons = _terminal_gate(
                source, interval, draft.metrics.medoid, rule
            )
            if reasons:
                status = (
                    "contig_boundary"
                    if reasons == ("source_boundary",)
                    else "repeat_gate_failed"
                )
                flanking_stops.append(FlankingLatticeStopEvidence(
                    side=side,
                    lattice_slot=next_slot,
                    interval=(
                        interval
                        if 0 <= interval[0] < interval[1] <= len(source)
                        else None
                    ),
                    status=status,
                    repeat_identity=identity,
                    exact_anchor_length=exact_anchor,
                    source_sha256=(
                        _slice_sha256(source, interval)
                        if 0 <= interval[0] < interval[1] <= len(source)
                        else None
                    ),
                    gate_reasons=reasons,
                ))
                break
            if (
                added >= rule.maximum_terminal_extensions_per_side
                or len(slot_intervals) >= rule.maximum_total_repeat_count
            ):
                cap_exhausted = True
                flanking_stops.append(FlankingLatticeStopEvidence(
                    side=side,
                    lattice_slot=next_slot,
                    interval=interval,
                    status="cap_exhausted_with_passing_repeat",
                    repeat_identity=identity,
                    exact_anchor_length=exact_anchor,
                    source_sha256=_slice_sha256(source, interval),
                    gate_reasons=(),
                ))
                break
            slot_intervals[next_slot] = interval
            added += 1

    slots = tuple(sorted(slot_intervals))
    intervals = tuple(slot_intervals[slot] for slot in slots)
    return _ExtendedDraft(
        phase=draft,
        lattice_slots=slots,
        repeat_intervals=intervals,
        metrics=_geometry_metrics(source, intervals, rule),
        terminal_added_slots=tuple(
            slot for slot in slots if slot not in initial_slots
        ),
        flanking_stops=tuple(sorted(
            flanking_stops, key=lambda item: item.side
        )),
        cap_exhausted=cap_exhausted,
    )


def _alternative(
    source,
    extended,
    period_fit,
    baseline_score,
    best_score,
):
    intervals = extended.repeat_intervals
    slots = extended.lattice_slots
    spacers = tuple(
        (left_end, right_start)
        for (_left_start, left_end), (right_start, _right_end)
        in zip(intervals, intervals[1:])
    )
    metrics = extended.metrics
    observed_slots = set(period_fit.observed_slots)
    internal_slots = set(extended.phase.lattice_slots)
    return SourcePeriodicLatticeAlternative(
        left_shift=extended.phase.left_shift,
        right_shift=extended.phase.right_shift,
        period=period_fit.period,
        repeat_length=intervals[0][1] - intervals[0][0],
        phase_source_score=extended.phase.source_score,
        phase_score_improvement=(
            extended.phase.source_score - baseline_score
        ),
        phase_score_delta_from_best=(
            best_score - extended.phase.source_score
        ),
        lattice_slots=slots,
        repeat_intervals=intervals,
        repeat_source_sha256s=tuple(
            _slice_sha256(source, interval) for interval in intervals
        ),
        spacer_intervals=spacers,
        spacer_source_sha256s=tuple(
            _slice_sha256(source, interval) for interval in spacers
        ),
        envelope=(intervals[0][0], intervals[-1][1]),
        envelope_source_sha256=_slice_sha256(
            source, (intervals[0][0], intervals[-1][1])
        ),
        source_medoid=metrics.medoid,
        source_medoid_sha256=source_sequence_sha256(metrics.medoid),
        repeat_identities=metrics.repeat_identities,
        exact_anchor_lengths=metrics.exact_anchor_lengths,
        mean_repeat_identity=metrics.mean_repeat_identity,
        strong_repeat_count=metrics.strong_repeat_count,
        weak_repeat_slots=tuple(
            slots[index] for index in metrics.weak_indexes
        ),
        spacer_lengths=metrics.spacer_lengths,
        spacer_length_cv=metrics.spacer_length_cv,
        distinct_spacer_count=metrics.distinct_spacer_count,
        maximum_distinct_spacer_identity=(
            metrics.maximum_distinct_spacer_identity
        ),
        terminal_kmer_dominance=metrics.terminal_kmer_dominance,
        internal_missing_slots=tuple(sorted(
            internal_slots - observed_slots
        )),
        terminal_added_slots=extended.terminal_added_slots,
        parent_mappings=extended.phase.parent_mappings,
        flanking_stops=extended.flanking_stops,
    )


def _reason_counts(drafts):
    counts = {}
    for draft in drafts:
        for reason in draft.metrics.gate_reasons:
            counts[reason] = counts.get(reason, 0) + 1
    return tuple(sorted(counts.items()))


def _freeze_phase_boundary_evidence(raw_evidence):
    def freeze_side(side):
        evidence = raw_evidence[side]
        return PhaseBoundarySideEvidence(
            varies=evidence["varies"],
            passes=evidence["passes"],
            included_dominance=evidence["included"]["dominance"],
            immediately_excluded_dominance=(
                evidence["immediately_excluded"]["dominance"]
            ),
            dominance_improvement=evidence["dominance_improvement"],
            included_discordant_canonical_count=(
                evidence["included_discordant_canonical_count"]
            ),
            evidence_branch=evidence["evidence_branch"],
        )

    return PhaseBoundaryEvidence(
        left=freeze_side("left"),
        right=freeze_side("right"),
    )


def _result(
    *,
    status,
    reason,
    source_hash,
    safety_evidence,
    rule,
    parent,
    period_fits=(),
    selected_period=None,
    observed_slots=(),
    capacity=0,
    evaluated=0,
    geometric_count=0,
    structural_evaluated_count=0,
    structural_count=0,
    structural_pruned_count=0,
    background_interval=None,
    background_probabilities=(),
    baseline_score=None,
    best_score=None,
    phase_boundary_evidence=None,
    alternatives=(),
    representative_index=None,
    count_invariant=False,
    locations_invariant=False,
    invariant_starts=(),
    boundary_status="unassessed",
    classification_ceiling=None,
    cap_exhausted=False,
    gate_reason_counts=(),
):
    return SourcePeriodicLatticeProposalResult(
        status=status,
        reason=reason,
        coordinate_system=COORDINATE_SYSTEM,
        rule_name=RULE_NAME,
        rule_version=RULE_VERSION,
        rule=rule,
        source_sha256=source_hash,
        safety_evidence=safety_evidence,
        parent_repeat_intervals=tuple(parent),
        parent_repeat_source_sha256s=(),
        period_fits=tuple(period_fits),
        selected_period=selected_period,
        observed_lattice_slots=tuple(observed_slots),
        phase_hypothesis_capacity=capacity,
        evaluated_phase_hypothesis_count=evaluated,
        geometrically_valid_phase_hypothesis_count=geometric_count,
        structurally_evaluated_phase_hypothesis_count=(
            structural_evaluated_count
        ),
        structurally_valid_phase_hypothesis_count=structural_count,
        structurally_pruned_phase_hypothesis_count=structural_pruned_count,
        plausible_phase_alternative_count=len(alternatives),
        background_interval=background_interval,
        background_source_sha256=None,
        background_probabilities=tuple(background_probabilities),
        baseline_phase_source_score=baseline_score,
        best_phase_source_score=best_score,
        phase_boundary_evidence=phase_boundary_evidence,
        alternatives=tuple(alternatives),
        representative_index=representative_index,
        repeat_count_invariant=count_invariant,
        repeat_start_locations_invariant=locations_invariant,
        invariant_repeat_starts=tuple(invariant_starts),
        exact_boundary_status=boundary_status,
        classification_ceiling=classification_ceiling,
        cap_exhausted=cap_exhausted,
        gate_reason_counts=tuple(gate_reason_counts),
    )


def _finalize_result_source_hashes(source, result):
    """Populate hashes that depend only on the validated source string."""

    return SourcePeriodicLatticeProposalResult(
        status=result.status,
        reason=result.reason,
        coordinate_system=result.coordinate_system,
        rule_name=result.rule_name,
        rule_version=result.rule_version,
        rule=result.rule,
        source_sha256=result.source_sha256,
        safety_evidence=result.safety_evidence,
        parent_repeat_intervals=result.parent_repeat_intervals,
        parent_repeat_source_sha256s=tuple(
            _slice_sha256(source, interval)
            for interval in result.parent_repeat_intervals
        ),
        period_fits=result.period_fits,
        selected_period=result.selected_period,
        observed_lattice_slots=result.observed_lattice_slots,
        phase_hypothesis_capacity=result.phase_hypothesis_capacity,
        evaluated_phase_hypothesis_count=(
            result.evaluated_phase_hypothesis_count
        ),
        geometrically_valid_phase_hypothesis_count=(
            result.geometrically_valid_phase_hypothesis_count
        ),
        structurally_evaluated_phase_hypothesis_count=(
            result.structurally_evaluated_phase_hypothesis_count
        ),
        structurally_valid_phase_hypothesis_count=(
            result.structurally_valid_phase_hypothesis_count
        ),
        structurally_pruned_phase_hypothesis_count=(
            result.structurally_pruned_phase_hypothesis_count
        ),
        plausible_phase_alternative_count=(
            result.plausible_phase_alternative_count
        ),
        background_interval=result.background_interval,
        background_source_sha256=(
            _slice_sha256(source, result.background_interval)
            if result.background_interval is not None else None
        ),
        background_probabilities=result.background_probabilities,
        baseline_phase_source_score=result.baseline_phase_source_score,
        best_phase_source_score=result.best_phase_source_score,
        phase_boundary_evidence=result.phase_boundary_evidence,
        alternatives=result.alternatives,
        representative_index=result.representative_index,
        repeat_count_invariant=result.repeat_count_invariant,
        repeat_start_locations_invariant=(
            result.repeat_start_locations_invariant
        ),
        invariant_repeat_starts=result.invariant_repeat_starts,
        exact_boundary_status=result.exact_boundary_status,
        classification_ceiling=result.classification_ceiling,
        cap_exhausted=result.cap_exhausted,
        gate_reason_counts=result.gate_reason_counts,
    )


def propose_source_periodic_lattice_repair(
    full_dna,
    crispr_candidate,
    safety_evidence,
    rule=SOURCE_PERIODIC_LATTICE_RULE,
):
    """Return one immutable atomic lattice proposal or an explicit abstention.

    The result never mutates or classifies ``crispr_candidate``.  Invalid
    source/candidate contracts raise ``ValueError``; insufficient biological
    evidence returns ``status == "abstained"``.
    """

    prefilter = source_periodic_lattice_prefilter(
        full_dna, crispr_candidate, safety_evidence, rule=rule
    )
    parent = prefilter.parent_repeat_intervals
    source_hash = prefilter.source_sha256

    def abstain(reason, **kwargs):
        return _finalize_result_source_hashes(
            full_dna,
            _result(
                status="abstained",
                reason=reason,
                source_hash=source_hash,
                safety_evidence=safety_evidence,
                rule=rule,
                parent=parent,
                **kwargs
            ),
        )

    if not prefilter.eligible:
        return abstain(
            prefilter.reason,
            period_fits=prefilter.period_fits,
            selected_period=prefilter.selected_period,
            observed_slots=prefilter.observed_lattice_slots,
            capacity=prefilter.phase_hypothesis_capacity,
            cap_exhausted=prefilter.cap_exhausted,
        )

    period_fits = prefilter.period_fits
    period_fit = period_fits[0]

    drafts, capacity, cap_exhausted = _phase_drafts(
        full_dna, parent, period_fit, rule
    )
    if cap_exhausted:
        return abstain(
            "phase_hypothesis_cap_exceeded",
            period_fits=period_fits,
            selected_period=period_fit.period,
            observed_slots=period_fit.observed_slots,
            capacity=capacity,
            cap_exhausted=True,
        )
    if not drafts:
        return abstain(
            "no_geometrically_valid_phase_hypotheses",
            period_fits=period_fits,
            selected_period=period_fit.period,
            observed_slots=period_fit.observed_slots,
            capacity=capacity,
            evaluated=capacity,
        )

    background_interval, background_probabilities = _background(
        full_dna,
        tuple(draft.repeat_intervals for draft in drafts),
        rule.background_flank,
    )
    scored = _with_scores(
        full_dna, drafts, background_probabilities, rule
    )
    baseline = next(
        (
            draft for draft in scored
            if draft.left_shift == 0 and draft.right_shift == 0
        ),
        None,
    )
    common = {
        "period_fits": period_fits,
        "selected_period": period_fit.period,
        "observed_slots": period_fit.observed_slots,
        "capacity": capacity,
        "evaluated": capacity,
        "geometric_count": len(drafts),
        "background_interval": background_interval,
        "background_probabilities": tuple(sorted(
            background_probabilities.items()
        )),
        "gate_reason_counts": (),
    }
    if baseline is None or baseline.source_score is None:
        return abstain("baseline_phase_not_scorable", **common)

    score_ordered = tuple(sorted(
        (
            draft for draft in scored
            if draft.source_score is not None
        ),
        key=lambda draft: (
            -draft.source_score,
            abs(draft.left_shift) + abs(draft.right_shift),
            abs(draft.left_shift),
            abs(draft.right_shift),
            -draft.left_shift,
            -draft.right_shift,
        ),
    ))
    structurally_evaluated = []
    structurally_valid = []
    best_score = None
    for draft in score_ordered:
        if (
            best_score is not None
            and best_score - draft.source_score
            > rule.plausible_phase_score_delta
        ):
            break
        metrics = _geometry_metrics(
            full_dna, draft.repeat_intervals, rule
        )
        evaluated_draft = _PhaseDraft(
            left_shift=draft.left_shift,
            right_shift=draft.right_shift,
            lattice_slots=draft.lattice_slots,
            repeat_intervals=draft.repeat_intervals,
            parent_mappings=draft.parent_mappings,
            metrics=metrics,
            source_score=draft.source_score,
        )
        structurally_evaluated.append(evaluated_draft)
        if not metrics.gate_reasons:
            structurally_valid.append(evaluated_draft)
            if best_score is None:
                best_score = draft.source_score

    common.update({
        "structural_evaluated_count": len(structurally_evaluated),
        "structural_count": len(structurally_valid),
        "structural_pruned_count": (
            len(score_ordered) - len(structurally_evaluated)
        ),
        "baseline_score": baseline.source_score,
        "gate_reason_counts": _reason_counts(structurally_evaluated),
    })
    if not structurally_valid:
        return abstain("no_structurally_valid_phase_hypotheses", **common)

    ordered = tuple(structurally_valid)
    common["best_score"] = best_score
    if (
        best_score - baseline.source_score
        < rule.minimum_phase_score_improvement
    ):
        return abstain("phase_score_improvement_below_minimum", **common)
    plausible = ordered

    def geometry(draft):
        return (
            tuple(start for start, _end in draft.repeat_intervals),
            tuple(end for _start, end in draft.repeat_intervals),
        )

    winner_geometry = geometry(ordered[0])
    raw_boundary_evidence = _source_phase_boundary_evidence(
        full_dna,
        winner_geometry,
        tuple(geometry(draft) for draft in plausible),
    )
    phase_boundary_evidence = _freeze_phase_boundary_evidence(
        raw_boundary_evidence
    )
    plausible = tuple(
        draft for draft in plausible
        if (
            not phase_boundary_evidence.left.passes
            or geometry(draft)[0] == winner_geometry[0]
        )
        and (
            not phase_boundary_evidence.right.passes
            or geometry(draft)[1] == winner_geometry[1]
        )
    )
    common["phase_boundary_evidence"] = phase_boundary_evidence

    extended = tuple(
        _extend_terminals(
            full_dna, draft, period_fit.period, rule
        )
        for draft in plausible
    )
    if any(item.cap_exhausted for item in extended):
        return abstain(
            "terminal_extension_cap_exhausted",
            cap_exhausted=True,
            **common
        )
    if any(item.metrics.gate_reasons for item in extended):
        return abstain(
            "plausible_phase_terminal_structure_failed",
            **common
        )

    alternatives = tuple(
        _alternative(
            full_dna,
            item,
            period_fit,
            baseline.source_score,
            best_score,
        )
        for item in extended
    )
    if (
        alternatives[0].mean_repeat_identity
        < rule.minimum_representative_mean_repeat_identity
    ):
        return abstain(
            "representative_mean_repeat_identity_below_minimum",
            alternatives=alternatives, **common
        )
    counts = {len(alternative.repeat_intervals) for alternative in alternatives}
    starts_by_alternative = {
        tuple(start for start, _end in alternative.repeat_intervals)
        for alternative in alternatives
    }
    count_invariant = len(counts) == 1
    locations_invariant = len(starts_by_alternative) == 1
    invariant_starts = (
        next(iter(starts_by_alternative))
        if locations_invariant else ()
    )
    if not count_invariant or not locations_invariant:
        return abstain(
            "phase_alternatives_disagree_on_occurrence_lattice",
            alternatives=alternatives,
            count_invariant=count_invariant,
            locations_invariant=locations_invariant,
            invariant_starts=invariant_starts,
            **common
        )

    exact_geometry_count = len({
        alternative.repeat_intervals for alternative in alternatives
    })
    boundary_status = (
        "exact" if exact_geometry_count == 1 else "ambiguous"
    )
    parent_geometry = tuple(parent)
    if all(
        alternative.repeat_intervals == parent_geometry
        for alternative in alternatives
    ):
        return abstain(
            "no_geometry_change",
            alternatives=alternatives,
            count_invariant=True,
            locations_invariant=True,
            invariant_starts=invariant_starts,
            boundary_status=boundary_status,
            **common
        )

    result = _result(
        status="proposed",
        reason=(
            "source_lattice_supported_boundary_exact"
            if boundary_status == "exact"
            else "source_lattice_supported_boundary_ambiguous"
        ),
        source_hash=source_hash,
        safety_evidence=safety_evidence,
        rule=rule,
        parent=parent,
        alternatives=alternatives,
        representative_index=0,
        count_invariant=True,
        locations_invariant=True,
        invariant_starts=invariant_starts,
        boundary_status=boundary_status,
        classification_ceiling="Possible",
        **common
    )
    return _finalize_result_source_hashes(full_dna, result)


def replay_source_periodic_lattice_proposal(full_dna, result):
    """Reauthorize and materialize every complete proposal alternative.

    Reauthorization rebuilds the parent from recorded source intervals and
    reruns the bounded search.  Consequently a source mutation, forged hash,
    altered interval, incomplete alternative, or changed gate result fails
    before any candidate can be materialized.
    """

    if not isinstance(result, SourcePeriodicLatticeProposalResult):
        raise TypeError("result must be SourcePeriodicLatticeProposalResult")
    if result.status != "proposed":
        raise ValueError("only a proposed lattice result can be replayed")
    if source_sequence_sha256(full_dna) != result.source_sha256:
        raise ValueError("proposal source hash does not match source DNA")
    for interval, expected_hash in zip(
        result.parent_repeat_intervals,
        result.parent_repeat_source_sha256s,
    ):
        if _slice_sha256(full_dna, interval) != expected_hash:
            raise ValueError("parent repeat source hash mismatch")

    parent = build_source_backed_candidate_from_boundaries(
        full_dna,
        None,
        tuple(start for start, _end in result.parent_repeat_intervals),
        tuple(end for _start, end in result.parent_repeat_intervals),
    )
    regenerated = propose_source_periodic_lattice_repair(
        full_dna,
        parent,
        result.safety_evidence,
        rule=result.rule,
    )
    if regenerated != result:
        raise ValueError(
            "proposal does not match a fresh exhaustive source replay"
        )

    candidates = []
    for alternative in result.alternatives:
        for interval, expected_hash in zip(
            alternative.repeat_intervals,
            alternative.repeat_source_sha256s,
        ):
            if _slice_sha256(full_dna, interval) != expected_hash:
                raise ValueError("alternative repeat source hash mismatch")
        for interval, expected_hash in zip(
            alternative.spacer_intervals,
            alternative.spacer_source_sha256s,
        ):
            if _slice_sha256(full_dna, interval) != expected_hash:
                raise ValueError("alternative spacer source hash mismatch")
        if (
            _slice_sha256(full_dna, alternative.envelope)
            != alternative.envelope_source_sha256
        ):
            raise ValueError("alternative envelope source hash mismatch")
        candidate = build_source_backed_candidate_from_boundaries(
            full_dna,
            None,
            tuple(start for start, _end in alternative.repeat_intervals),
            tuple(end for _start, end in alternative.repeat_intervals),
        )
        validate_candidate_source_reconstruction(
            full_dna,
            candidate.list_repeats,
            candidate.list_repeat_starts,
            candidate.list_spacers,
        )
        candidates.append(candidate)
    return tuple(candidates)


