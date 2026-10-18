"""Bounded accepted-family rescue for two narrow fast-run failure modes.

The compact classifier occasionally rejects a source-valid satellite array
even though the same repeat family is already represented by an accepted call
on the same source.  This module deliberately does *not* lower the global
model threshold.  It exposes two hard-capped, comparator-blind construction
arms:

``family_terminal``
    A five-to-eight-repeat Low candidate is a high-identity prefix/suffix
    phase of a nearby Bona-fide repeat.  Extending that phase and adding one
    source terminal must change the unchanged strict source-intact decision
    from the single ``spacer_length_cv_above_maximum`` rejection to eligible.

``exact_template_satellite``
    A two-repeat candidate rejected only for repeat count contains the same
    exact template in both repeats as at least three nearby accepted arrays.
    The template is located independently in each parent repeat and one
    cadence-compatible terminal is added from source DNA.

Both arms are generation rules, not external-comparator rules.  They consume
only the current source, current scored rows, and current accepted calls; all
coordinates are replayed from source and every ambiguity/cap fails closed.
"""

from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import statistics

from components.components_evaluated_arrays_enhancement import (
    _longest_exact_source_anchor,
    _source_only_alignment_seed,
    _source_repeat_edit_distance,
    build_source_backed_candidate_from_boundaries,
    validate_candidate_source_reconstruction,
)
from components.fast_run_rescue import candidate_source_interval
from components.fast_two_repeat_terminal_repair import (
    MAXIMUM_ADMISSION_REPEAT_DERIVED_SPACER_COVERAGE,
    MINIMUM_REPEAT_FLANK_CONTINUATION_IDENTITY,
    MINIMUM_SPACER_FLANK_CONTINUATION_IDENTITY,
    flank_continuation_identities,
    repeat_derived_spacer_coverage,
    repeat_dust_triplet_score,
)
from components.source_intact_existence_rescue import (
    ELIGIBLE,
    SourceIntactEvaluationContext,
    evaluate_source_intact_existence_in_context,
)


CLASSIFICATION_METHOD = "source_backed_accepted_family_rescue"
FAMILY_TERMINAL = "family_terminal"
EXACT_TEMPLATE_SATELLITE = "exact_template_satellite"
LOW_SCORE_PARENT = "low_score_row"
RAW_TWO_REPEAT_SEED_PARENT = "raw_two_repeat_seed"
COORDINATE_SYSTEM = "zero_based_half_open"
SEARCH_IMPLEMENTATION_VERSION = "2.0.0"


@dataclass(frozen=True)
class AcceptedFamilyRescuePolicy:
    """Frozen policy for the opt-in family-supported admission stage."""

    name: str = "fast-accepted-family-rescue"
    version: str = "1.1.0"
    minimum_completed_score: float = 0.30
    maximum_completed_score_exclusive: float = 0.50

    terminal_parent_minimum_score: float = 0.20
    terminal_parent_minimum_repeats: int = 5
    terminal_parent_maximum_repeats: int = 8
    terminal_family_maximum_gap: int = 2000
    terminal_family_minimum_local_identity: float = 0.90
    terminal_family_minimum_length_delta: int = 1
    terminal_family_maximum_length_delta: int = 20
    terminal_minimum_identity: float = 0.80

    template_parent_minimum_score: float = 0.35
    template_family_maximum_gap: int = 25000
    template_minimum_distinct_accepted_arrays: int = 3
    template_minimum_exact_repeat_support_per_array: int = 8
    template_minimum_length: int = 20
    template_maximum_length: int = 55
    template_spacer_radius: int = 10
    template_terminal_minimum_identity: float = 0.80
    template_terminal_minimum_exact_anchor: int = 17
    template_minimum_repeat_identity: float = 0.80
    template_minimum_mean_repeat_identity: float = 0.90
    template_maximum_spacer_cv: float = 0.10
    template_maximum_spacer_identity: float = 0.60

    minimum_spacer_length: int = 18
    maximum_spacer_length: int = 80
    maximum_repeat_derived_spacer_coverage: float = (
        MAXIMUM_ADMISSION_REPEAT_DERIVED_SPACER_COVERAGE
    )
    minimum_repeat_flank_continuation_identity: float = (
        MINIMUM_REPEAT_FLANK_CONTINUATION_IDENTITY
    )
    minimum_spacer_flank_continuation_identity: float = (
        MINIMUM_SPACER_FLANK_CONTINUATION_IDENTITY
    )
    maximum_parent_rows_per_source: int = 4096
    maximum_raw_seeds_per_source: int = 1024
    maximum_raw_seed_parent_matches_per_source: int = 8
    maximum_raw_seed_terminal_windows_per_parent: int = 64
    maximum_family_anchors_per_parent: int = 64
    maximum_terminal_windows_per_parent: int = 512
    maximum_proposals_per_source: int = 64

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("family rescue policy name must be non-empty")
        if not isinstance(self.version, str) or not self.version:
            raise ValueError("family rescue policy version must be non-empty")
        for field_name in (
            "terminal_parent_minimum_repeats",
            "terminal_parent_maximum_repeats",
            "terminal_family_maximum_gap",
            "terminal_family_minimum_length_delta",
            "terminal_family_maximum_length_delta",
            "template_family_maximum_gap",
            "template_minimum_distinct_accepted_arrays",
            "template_minimum_exact_repeat_support_per_array",
            "template_minimum_length",
            "template_maximum_length",
            "template_spacer_radius",
            "template_terminal_minimum_exact_anchor",
            "minimum_spacer_length",
            "maximum_spacer_length",
            "maximum_parent_rows_per_source",
            "maximum_raw_seeds_per_source",
            "maximum_raw_seed_parent_matches_per_source",
            "maximum_raw_seed_terminal_windows_per_parent",
            "maximum_family_anchors_per_parent",
            "maximum_terminal_windows_per_parent",
            "maximum_proposals_per_source",
        ):
            value = getattr(self, field_name)
            if not isinstance(value, int) or isinstance(value, bool) or value < 1:
                raise ValueError("{} must be a positive integer".format(field_name))
        if self.terminal_parent_minimum_repeats > self.terminal_parent_maximum_repeats:
            raise ValueError("terminal parent repeat bounds are inconsistent")
        if self.terminal_family_minimum_length_delta > self.terminal_family_maximum_length_delta:
            raise ValueError("terminal family length deltas are inconsistent")
        if self.template_minimum_length > self.template_maximum_length:
            raise ValueError("template length bounds are inconsistent")
        if self.minimum_spacer_length > self.maximum_spacer_length:
            raise ValueError("family rescue spacer bounds are inconsistent")
        for field_name in (
            "minimum_completed_score",
            "maximum_completed_score_exclusive",
            "terminal_parent_minimum_score",
            "terminal_family_minimum_local_identity",
            "terminal_minimum_identity",
            "template_parent_minimum_score",
            "template_terminal_minimum_identity",
            "template_minimum_repeat_identity",
            "template_minimum_mean_repeat_identity",
            "template_maximum_spacer_cv",
            "template_maximum_spacer_identity",
            "maximum_repeat_derived_spacer_coverage",
            "minimum_repeat_flank_continuation_identity",
            "minimum_spacer_flank_continuation_identity",
        ):
            value = getattr(self, field_name)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
                or not 0.0 <= float(value) <= 1.0
            ):
                raise ValueError("{} must be finite in [0, 1]".format(field_name))
        if self.minimum_completed_score >= self.maximum_completed_score_exclusive:
            raise ValueError("completed-score bounds are inconsistent")

    def as_dict(self):
        return asdict(self)


ACCEPTED_FAMILY_RESCUE_POLICY = AcceptedFamilyRescuePolicy()


@dataclass(frozen=True)
class AcceptedFamilyAnchor:
    category_index: int
    score: float
    repeat_intervals: tuple
    consensus: str

    def as_dict(self):
        return {
            "category_index": self.category_index,
            "score": self.score,
            "repeat_intervals": [list(interval) for interval in self.repeat_intervals],
            "consensus": self.consensus,
        }


@dataclass(frozen=True)
class AcceptedFamilyRawSeedParent:
    """One source-rebuilt raw seed that passed the exact-family prefilter."""

    candidate: object
    seed_repeat_intervals: tuple
    cluster_indices: tuple

    def __post_init__(self):
        if _geometry(self.candidate) != self.seed_repeat_intervals:
            raise ValueError("raw-seed parent geometry differs from its seed")
        if (
            not self.cluster_indices
            or tuple(sorted(set(self.cluster_indices))) != self.cluster_indices
            or any(
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 0
                for value in self.cluster_indices
            )
        ):
            raise ValueError("raw-seed parent cluster indexes are invalid")


@dataclass(frozen=True)
class AcceptedFamilyRescueProposal:
    arm: str
    candidate: object
    parent_row: object
    parent_score: float
    parent_repeat_intervals: tuple
    family_template: str
    anchors: tuple
    phase_offset: object
    terminal_side: str
    terminal_start: int
    source_metrics: dict
    base_source_intact_decision: object = None
    completed_source_intact_decision: object = None
    parent_origin: str = LOW_SCORE_PARENT
    raw_seed_cluster_indices: tuple = ()

    @property
    def repeat_intervals(self):
        return _geometry(self.candidate)

    def source_rank(self):
        metrics = self.source_metrics
        return (
            round(float(metrics["terminal_identity"]), 12),
            int(metrics["terminal_exact_anchor"]),
            round(float(metrics["minimum_repeat_identity"]), 12),
            round(float(metrics["mean_repeat_identity"]), 12),
            -round(float(metrics["spacer_cv"]), 12),
            -round(float(metrics["maximum_spacer_identity"]), 12),
            -int(metrics["accepted_family_minimum_gap"]),
        )

    def as_evidence(self, source_sha256, fresh_score, policy):
        evidence = {
            "arm": self.arm,
            "coordinate_system": COORDINATE_SYSTEM,
            "source_sha256": source_sha256,
            "parent_score": float(self.parent_score),
            "parent_repeat_intervals": [
                list(interval) for interval in self.parent_repeat_intervals
            ],
            "parent_origin": self.parent_origin,
            "raw_seed_cluster_indices": list(self.raw_seed_cluster_indices),
            "family_template": self.family_template,
            "accepted_family_anchors": [
                anchor.as_dict() for anchor in self.anchors
            ],
            "phase_offset": self.phase_offset,
            "terminal_side": self.terminal_side,
            "terminal_start": self.terminal_start,
            "proposal_repeat_intervals": [
                list(interval) for interval in self.repeat_intervals
            ],
            "source_metrics": _json_metrics(self.source_metrics),
            "fresh_score": float(fresh_score),
            "policy": policy.as_dict(),
        }
        if self.base_source_intact_decision is not None:
            evidence["base_source_intact_decision"] = (
                self.base_source_intact_decision.as_dict()
            )
        if self.completed_source_intact_decision is not None:
            evidence["completed_source_intact_decision"] = (
                self.completed_source_intact_decision.as_dict()
            )
        return evidence


def _json_metrics(metrics):
    result = {}
    for key, value in sorted(metrics.items()):
        if isinstance(value, tuple):
            result[key] = list(value)
        elif isinstance(value, (int, float, str, bool)) or value is None:
            result[key] = value
        else:
            raise TypeError("family rescue metric is not JSON-safe: {}".format(key))
    return result


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _geometry_digest(geometries):
    payload = json.dumps(
        [
            [[int(start), int(end)] for start, end in geometry]
            for geometry in geometries
        ],
        separators=(",", ":"),
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _clustered_geometry_digest(geometries, grouped_clusters):
    payload = json.dumps(
        [
            {
                "geometry": [
                    [int(start), int(end)] for start, end in geometry
                ],
                "cluster_indices": sorted(grouped_clusters[geometry]),
            }
            for geometry in geometries
        ],
        separators=(",", ":"),
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN"))[::-1]


def _canonical_repeat(sequence):
    reverse = _reverse_complement(sequence)
    return min(sequence, reverse)


def _edit_identity(first, second):
    denominator = max(len(first), len(second))
    if denominator == 0:
        return 1.0
    return 1.0 - _source_repeat_edit_distance(first, second) / float(denominator)


def _geometry(candidate):
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(candidate.list_repeat_starts, candidate.list_repeats)
    )


def _source_gap(first, second):
    first_start, first_end = first
    second_start, second_end = second
    return max(0, second_start - first_end, first_start - second_end)


def _strict_overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _spacer_statistics(spacers):
    lengths = tuple(map(len, spacers))
    mean = sum(lengths) / float(len(lengths)) if lengths else 0.0
    cv = statistics.pstdev(lengths) / mean if mean else 0.0
    maximum_identity = max((
        _edit_identity(first, second)
        for index, first in enumerate(spacers[:-1])
        for second in spacers[index + 1:]
    ), default=0.0)
    return cv, maximum_identity


def _repeat_statistics(repeats):
    medoid = _source_only_alignment_seed(tuple(repeats))
    identities = tuple(_edit_identity(repeat, medoid) for repeat in repeats)
    return min(identities), sum(identities) / float(len(identities))


def _failed_structural_predicates(structural_filter, candidate):
    evidence = structural_filter.evidence(candidate)
    failed = tuple(
        name for name, status in evidence.predicate_results if status == "failed"
    )
    return evidence, failed


def _family_guard_metrics(source, candidate, terminal_template):
    repeats = tuple(candidate.list_repeats)
    spacers = tuple(candidate.list_spacers)
    minimum_identity, mean_identity = _repeat_statistics(repeats)
    spacer_cv, spacer_identity = _spacer_statistics(spacers)
    coverage = repeat_derived_spacer_coverage(candidate)
    dust = repeat_dust_triplet_score(candidate)
    flank = flank_continuation_identities(source, candidate)
    terminal = repeats[0] if candidate.list_repeat_starts[0] < min(
        candidate.list_repeat_starts[1:]
    ) else repeats[-1]
    # Callers overwrite terminal_identity/anchor with the actual terminal when
    # the terminal is not inferable solely from ordering.
    return {
        "minimum_repeat_identity": minimum_identity,
        "mean_repeat_identity": mean_identity,
        "spacer_cv": spacer_cv,
        "maximum_spacer_identity": spacer_identity,
        "maximum_repeat_derived_spacer_coverage": coverage,
        "repeat_dust_triplet_score": dust,
        "maximum_repeat_flank_identity": flank[
            "maximum_repeat_flank_identity"
        ],
        "maximum_spacer_flank_identity": flank[
            "maximum_spacer_flank_identity"
        ],
        "terminal_identity": _edit_identity(terminal, terminal_template),
        "terminal_exact_anchor": _longest_exact_source_anchor(
            terminal, terminal_template
        ),
    }


def _common_guards_pass(metrics, policy):
    if metrics["maximum_repeat_derived_spacer_coverage"] >= (
        policy.maximum_repeat_derived_spacer_coverage
    ):
        return False
    source_continues = (
        metrics["maximum_repeat_flank_identity"]
        >= policy.minimum_repeat_flank_continuation_identity
        and metrics["maximum_spacer_flank_identity"]
        >= policy.minimum_spacer_flank_continuation_identity
    )
    return not source_continues


def _best_unique_local_alignment(query, reference, minimum_identity):
    if not query or len(query) > len(reference):
        return None
    scored = []
    for offset in range(len(reference) - len(query) + 1):
        window = reference[offset:offset + len(query)]
        identity = sum(first == second for first, second in zip(query, window)) / float(len(query))
        scored.append((identity, offset))
    best = max((record[0] for record in scored), default=-1.0)
    offsets = tuple(offset for identity, offset in scored if identity == best)
    if best < minimum_identity or len(offsets) != 1:
        return None
    return best, offsets[0]


def _candidate_from_intervals(source, intervals):
    return build_source_backed_candidate_from_boundaries(
        source,
        None,
        tuple(start for start, _end in intervals),
        tuple(end for _start, end in intervals),
    )


def _accepted_records(categories):
    records = []
    for category_index in (0, 2):
        for rows in categories[category_index].values():
            for row in rows:
                candidate = row[1]
                records.append({
                    "category_index": category_index,
                    "row": row,
                    "score": float(row[0]),
                    "candidate": candidate,
                    "interval": candidate_source_interval(candidate),
                    "geometry": _geometry(candidate),
                })
    return tuple(records)


def _anchor(record):
    candidate = record["candidate"]
    return AcceptedFamilyAnchor(
        category_index=record["category_index"],
        score=record["score"],
        repeat_intervals=record["geometry"],
        consensus=candidate.consensus,
    )


def _exact_consensus_support(record):
    candidate = record["candidate"]
    return sum(
        repeat == candidate.consensus for repeat in candidate.list_repeats
    )


def _nonoverlapping_supported_family_records(records, policy):
    """Return one deterministic, strongly supported anchor per source locus."""

    supported = tuple(
        record
        for record in records
        if _exact_consensus_support(record)
        >= policy.template_minimum_exact_repeat_support_per_array
    )
    if not supported:
        return ()
    ordered = sorted(supported, key=lambda record: record["interval"])
    components = []
    current = []
    current_end = None
    for record in ordered:
        start, end = record["interval"]
        if not current or start < current_end:
            current.append(record)
            current_end = end if current_end is None else max(current_end, end)
        else:
            components.append(current)
            current = [record]
            current_end = end
    if current:
        components.append(current)
    return tuple(
        sorted(component, key=lambda record: (
            -_exact_consensus_support(record),
            record["category_index"],
            -record["score"],
            -len(record["candidate"].list_repeats),
            record["interval"],
            record["geometry"],
        ))[0]
        for component in components
    )


def _supported_template_family_records(accepted, policy):
    """Index the monotonic global prerequisite for the template arm.

    Every admissible local template family must draw at least three records
    from one exact-supported canonical family on the source.  A local distance
    filter can only remove records, so canonical families with fewer than the
    required global record count can never admit a parent.  We intentionally
    use record count rather than overlap-component count: adding a spanning
    record can merge components, so component count is not monotonic.
    """

    families = {}
    for record in accepted:
        consensus = record["candidate"].consensus
        if not (
            policy.template_minimum_length <= len(consensus)
            <= policy.template_maximum_length
            and _exact_consensus_support(record)
            >= policy.template_minimum_exact_repeat_support_per_array
        ):
            continue
        families.setdefault(_canonical_repeat(consensus), []).append(record)
    return {
        canonical: tuple(records)
        for canonical, records in families.items()
        if len(records) >= policy.template_minimum_distinct_accepted_arrays
    }


def _terminal_family_proposals(
    source,
    context_cache,
    parent_row,
    accepted,
    structural_filter,
    policy,
    counters,
):
    parent = parent_row[1]
    score = float(parent_row[0])
    repeats = tuple(parent.list_repeats)
    if not (
        policy.terminal_parent_minimum_score <= score
        < policy.maximum_completed_score_exclusive
        and policy.terminal_parent_minimum_repeats <= len(repeats)
        <= policy.terminal_parent_maximum_repeats
    ):
        return ()
    evidence, failed = _failed_structural_predicates(structural_filter, parent)
    if not evidence.structural_valid or failed:
        return ()
    parent_interval = candidate_source_interval(parent)
    if any(_strict_overlap(parent_interval, record["interval"]) for record in accepted):
        return ()
    query = parent.consensus
    if not isinstance(query, str) or not query:
        return ()

    anchors = []
    for record in accepted:
        if record["category_index"] != 0:
            continue
        gap = _source_gap(parent_interval, record["interval"])
        if gap > policy.terminal_family_maximum_gap:
            continue
        accepted_consensus = record["candidate"].consensus
        delta = len(accepted_consensus) - len(query)
        if not (
            policy.terminal_family_minimum_length_delta <= delta
            <= policy.terminal_family_maximum_length_delta
        ):
            continue
        for reference in tuple(dict.fromkeys((
            accepted_consensus,
            _reverse_complement(accepted_consensus),
        ))):
            alignment = _best_unique_local_alignment(
                query, reference, policy.terminal_family_minimum_local_identity
            )
            if alignment is None:
                continue
            anchors.append((record, reference, alignment[0], alignment[1], gap))
    counters["terminal_family_anchor_matches"] += len(anchors)
    if not anchors or len(anchors) > policy.maximum_family_anchors_per_parent:
        if len(anchors) > policy.maximum_family_anchors_per_parent:
            counters["family_anchor_cap_exhausted"] += 1
        return ()

    proposals = {}
    windows_tested = 0
    for record, reference, local_identity, phase_offset, gap in anchors:
        length = len(reference)
        starts = tuple(start - phase_offset for start in parent.list_repeat_starts)
        intervals = tuple((start, start + length) for start in starts)
        if (
            any(start < 0 or end > len(source) for start, end in intervals)
            or any(left[1] > right[0] for left, right in zip(intervals, intervals[1:]))
        ):
            continue
        try:
            base = _candidate_from_intervals(source, intervals)
        except ValueError:
            continue
        context = context_cache.get("source_intact")
        if context is None:
            context = SourceIntactEvaluationContext(source)
            context_cache["source_intact"] = context
            counters["source_intact_context_constructions"] += 1
        base_decision = evaluate_source_intact_existence_in_context(context, base)
        base_codes = tuple(reason.code for reason in base_decision.reasons)
        if base_decision.decision == ELIGIBLE or base_codes != (
            "spacer_length_cv_above_maximum",
        ):
            continue

        for side in ("left", "right"):
            if side == "right":
                lower = intervals[-1][1] + policy.minimum_spacer_length
                upper = intervals[-1][1] + policy.maximum_spacer_length
            else:
                lower = intervals[0][0] - policy.maximum_spacer_length - length
                upper = intervals[0][0] - policy.minimum_spacer_length - length
            for terminal_start in range(max(0, lower), min(len(source) - length, upper) + 1):
                windows_tested += 1
                counters["terminal_windows_tested"] += 1
                if windows_tested > policy.maximum_terminal_windows_per_parent:
                    counters["terminal_window_cap_exhausted"] += 1
                    return ()
                terminal = source[terminal_start:terminal_start + length]
                terminal_identity = _edit_identity(terminal, reference)
                if terminal_identity < policy.terminal_minimum_identity:
                    continue
                completed_intervals = (
                    intervals + ((terminal_start, terminal_start + length),)
                    if side == "right"
                    else ((terminal_start, terminal_start + length),) + intervals
                )
                try:
                    candidate = _candidate_from_intervals(source, completed_intervals)
                except ValueError:
                    continue
                if any(
                    _strict_overlap(candidate_source_interval(candidate), item["interval"])
                    for item in accepted
                ):
                    continue
                completed = evaluate_source_intact_existence_in_context(context, candidate)
                if completed.decision != ELIGIBLE or not structural_filter(candidate):
                    continue
                metrics = _family_guard_metrics(source, candidate, reference)
                terminal_index = 0 if side == "left" else len(candidate.list_repeats) - 1
                terminal_repeat = candidate.list_repeats[terminal_index]
                metrics["terminal_identity"] = _edit_identity(terminal_repeat, reference)
                metrics["terminal_exact_anchor"] = _longest_exact_source_anchor(
                    terminal_repeat, reference
                )
                metrics["accepted_family_minimum_gap"] = gap
                metrics["parent_family_local_identity"] = local_identity
                if not _common_guards_pass(metrics, policy):
                    continue
                proposal = AcceptedFamilyRescueProposal(
                    arm=FAMILY_TERMINAL,
                    candidate=candidate,
                    parent_row=parent_row,
                    parent_score=score,
                    parent_repeat_intervals=_geometry(parent),
                    family_template=reference,
                    anchors=(_anchor(record),),
                    phase_offset=phase_offset,
                    terminal_side=side,
                    terminal_start=terminal_start,
                    source_metrics=metrics,
                    base_source_intact_decision=base_decision,
                    completed_source_intact_decision=completed,
                )
                geometry = proposal.repeat_intervals
                prior = proposals.get(geometry)
                if prior is None or proposal.source_rank() > prior.source_rank():
                    proposals[geometry] = proposal
    return tuple(proposals.values())


def _template_satellite_proposals(
    source,
    parent_row,
    accepted,
    structural_filter,
    policy,
    counters,
    *,
    parent_origin=LOW_SCORE_PARENT,
    raw_seed_cluster_indices=(),
    maximum_terminal_windows=None,
    supported_template_families=None,
):
    if parent_origin not in (LOW_SCORE_PARENT, RAW_TWO_REPEAT_SEED_PARENT):
        raise ValueError("accepted-family parent origin is unsupported")
    raw_seed_cluster_indices = tuple(raw_seed_cluster_indices)
    if parent_origin == LOW_SCORE_PARENT and raw_seed_cluster_indices:
        raise ValueError("Low-score parent cannot carry raw-seed clusters")
    if parent_origin == RAW_TWO_REPEAT_SEED_PARENT and not raw_seed_cluster_indices:
        raise ValueError("raw-seed parent lost its source clusters")
    if maximum_terminal_windows is None:
        maximum_terminal_windows = policy.maximum_terminal_windows_per_parent
    if (
        not isinstance(maximum_terminal_windows, int)
        or isinstance(maximum_terminal_windows, bool)
        or maximum_terminal_windows < 1
    ):
        raise ValueError("template terminal-window cap must be positive")
    parent = parent_row[1]
    score = float(parent_row[0])
    if not (
        policy.template_parent_minimum_score <= score
        < policy.maximum_completed_score_exclusive
        and len(parent.list_repeats) == 2
        and len(parent.list_spacers) == 1
    ):
        return ()
    evidence, failed = _failed_structural_predicates(structural_filter, parent)
    # ``structural_valid`` means that all non-cardinality morphology checks
    # are usable; the canonical filter still records the sole failed
    # predicate separately. This is the same entrance contract used by the
    # existing two-repeat terminal repair.
    if not evidence.structural_valid or failed != ("minimum_repeat_count",):
        return ()
    parent_interval = candidate_source_interval(parent)
    if any(_strict_overlap(parent_interval, record["interval"]) for record in accepted):
        return ()

    if supported_template_families is None:
        supported_template_families = _supported_template_family_records(
            accepted, policy
        )
    families = {
        canonical: [
            record for record in records
            if _source_gap(parent_interval, record["interval"])
            <= policy.template_family_maximum_gap
        ]
        for canonical, records in supported_template_families.items()
    }

    proposals = {}
    windows_tested = 0
    parent_geometry = _geometry(parent)
    for canonical, records in families.items():
        records = _nonoverlapping_supported_family_records(records, policy)
        if len(records) < policy.template_minimum_distinct_accepted_arrays:
            continue
        if len(records) > policy.maximum_family_anchors_per_parent:
            counters["family_anchor_cap_exhausted"] += 1
            continue
        for template in tuple(dict.fromkeys((canonical, _reverse_complement(canonical)))):
            offsets = []
            for repeat in parent.list_repeats:
                first = repeat.find(template)
                if first < 0 or repeat.find(template, first + 1) >= 0:
                    offsets = []
                    break
                offsets.append(first)
            if len(offsets) != 2:
                continue
            core_starts = tuple(
                start + offset
                for start, offset in zip(parent.list_repeat_starts, offsets)
            )
            length = len(template)
            core = tuple((start, start + length) for start in core_starts)
            first_spacer = core[1][0] - core[0][1]
            spacer_lower = max(
                policy.minimum_spacer_length,
                first_spacer - policy.template_spacer_radius,
            )
            spacer_upper = min(
                policy.maximum_spacer_length,
                first_spacer + policy.template_spacer_radius,
            )
            if spacer_lower > spacer_upper:
                continue
            for side in ("left", "right"):
                if side == "right":
                    lower = core[-1][1] + spacer_lower
                    upper = core[-1][1] + spacer_upper
                else:
                    lower = core[0][0] - spacer_upper - length
                    upper = core[0][0] - spacer_lower - length
                for terminal_start in range(max(0, lower), min(len(source) - length, upper) + 1):
                    windows_tested += 1
                    counters["template_terminal_windows_tested"] += 1
                    if windows_tested > maximum_terminal_windows:
                        counters["terminal_window_cap_exhausted"] += 1
                        return ()
                    terminal = source[terminal_start:terminal_start + length]
                    terminal_identity = _edit_identity(terminal, template)
                    terminal_anchor = _longest_exact_source_anchor(terminal, template)
                    if (
                        terminal_identity < policy.template_terminal_minimum_identity
                        or terminal_anchor < policy.template_terminal_minimum_exact_anchor
                    ):
                        continue
                    intervals = (
                        core + ((terminal_start, terminal_start + length),)
                        if side == "right"
                        else ((terminal_start, terminal_start + length),) + core
                    )
                    try:
                        candidate = _candidate_from_intervals(source, intervals)
                    except ValueError:
                        continue
                    if any(
                        _strict_overlap(candidate_source_interval(candidate), item["interval"])
                        for item in accepted
                    ):
                        continue
                    if not structural_filter(candidate):
                        continue
                    metrics = _family_guard_metrics(source, candidate, template)
                    terminal_index = 0 if side == "left" else 2
                    terminal_repeat = candidate.list_repeats[terminal_index]
                    metrics["terminal_identity"] = _edit_identity(terminal_repeat, template)
                    metrics["terminal_exact_anchor"] = _longest_exact_source_anchor(
                        terminal_repeat, template
                    )
                    metrics["accepted_family_minimum_gap"] = min(
                        _source_gap(parent_interval, record["interval"])
                        for record in records
                    )
                    metrics["accepted_family_array_count"] = len(records)
                    metrics["template_offsets"] = tuple(offsets)
                    if (
                        metrics["minimum_repeat_identity"]
                        < policy.template_minimum_repeat_identity
                        or metrics["mean_repeat_identity"]
                        < policy.template_minimum_mean_repeat_identity
                        or metrics["spacer_cv"] > policy.template_maximum_spacer_cv
                        or metrics["maximum_spacer_identity"]
                        > policy.template_maximum_spacer_identity
                        or not _common_guards_pass(metrics, policy)
                    ):
                        continue
                    proposal = AcceptedFamilyRescueProposal(
                        arm=EXACT_TEMPLATE_SATELLITE,
                        candidate=candidate,
                        parent_row=parent_row,
                        parent_score=score,
                        parent_repeat_intervals=parent_geometry,
                        family_template=template,
                        anchors=tuple(_anchor(record) for record in records),
                        phase_offset=tuple(offsets),
                        terminal_side=side,
                        terminal_start=terminal_start,
                        source_metrics=metrics,
                        parent_origin=parent_origin,
                        raw_seed_cluster_indices=raw_seed_cluster_indices,
                    )
                    geometry = proposal.repeat_intervals
                    prior = proposals.get(geometry)
                    if prior is None or proposal.source_rank() > prior.source_rank():
                        proposals[geometry] = proposal
    return tuple(proposals.values())


def _select_unique_source_winner(proposals, counters):
    proposals = tuple(proposals)
    if not proposals:
        return None
    ordered = sorted(proposals, key=lambda item: item.source_rank(), reverse=True)
    if len(ordered) > 1 and ordered[0].source_rank() == ordered[1].source_rank():
        counters["ambiguous_source_rank"] += 1
        return None
    return ordered[0]


def select_exact_template_raw_seed_parents(
    source,
    categories,
    raw_seeds,
    structural_filter,
    *,
    accepted_categories=None,
    policy=ACCEPTED_FAMILY_RESCUE_POLICY,
):
    """Preselect raw two-repeat seeds without broadly scoring rejected rows.

    A seed reaches compact-model scoring only if the source-rebuilt parent can
    already produce a unique exact-template-satellite proposal at the frozen
    lower parent-score boundary.  Thus scoring work is proportional to rare
    accepted-family matches rather than to every rejected two-repeat path.
    """

    if not isinstance(source, str) or not source:
        raise ValueError("raw-seed family rescue source must be non-empty")
    source = source.upper()
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("raw-seed family rescue requires five categories")
    if accepted_categories is None:
        accepted_categories = categories
    if (
        not isinstance(accepted_categories, (list, tuple))
        or len(accepted_categories) != 5
    ):
        raise ValueError("raw-seed family anchor snapshot is malformed")
    if not callable(structural_filter) or not hasattr(structural_filter, "evidence"):
        raise TypeError("raw-seed family rescue requires the structural filter")
    try:
        raw_seeds = tuple(raw_seeds)
    except TypeError as error:
        raise ValueError("raw two-repeat seeds must be iterable") from error

    counters = Counter(
        raw_seeds_input=len(raw_seeds),
        raw_seed_cap=policy.maximum_raw_seeds_per_source,
        raw_seed_parent_match_cap=(
            policy.maximum_raw_seed_parent_matches_per_source
        ),
    )
    counters["search_implementation_version"] = (
        SEARCH_IMPLEMENTATION_VERSION
    )
    if len(raw_seeds) > policy.maximum_raw_seeds_per_source:
        counters["raw_seed_cap_exhausted"] = 1
        return (), dict(counters)

    accepted = _accepted_records(accepted_categories)
    supported_template_families = _supported_template_family_records(
        accepted, policy
    )
    categorized_geometries = set()
    for category in categories:
        for rows in category.values():
            for row in rows:
                if not isinstance(row, (tuple, list)) or len(row) < 2:
                    raise ValueError("raw-seed prefilter saw a malformed row")
                categorized_geometries.add(_geometry(row[1]))

    grouped_clusters = {}
    for seed in raw_seeds:
        try:
            starts = tuple(seed.repeat_starts)
            ends = tuple(seed.repeat_ends)
            cluster_index = seed.cluster_index
        except AttributeError as error:
            raise ValueError("raw two-repeat seed is malformed") from error
        if (
            len(starts) != 2
            or len(ends) != 2
            or any(
                not isinstance(value, int) or isinstance(value, bool)
                for value in starts + ends
            )
            or not (0 <= starts[0] < ends[0] <= starts[1] < ends[1] <= len(source))
            or not isinstance(cluster_index, int)
            or isinstance(cluster_index, bool)
            or cluster_index < 0
        ):
            raise ValueError("raw two-repeat seed geometry is invalid")
        geometry = tuple(zip(starts, ends))
        if geometry in categorized_geometries:
            counters["raw_seed_already_categorized"] += 1
            continue
        grouped_clusters.setdefault(geometry, set()).add(cluster_index)

    counters["raw_seed_unique_uncategorized_geometries"] = len(grouped_clusters)
    counters["raw_seed_geometry_duplicates"] = (
        len(raw_seeds)
        - counters["raw_seed_already_categorized"]
        - len(grouped_clusters)
    )
    counters["raw_seed_supported_template_families"] = len(
        supported_template_families
    )
    ordered_geometries = tuple(sorted(grouped_clusters))
    if not supported_template_families:
        counters["raw_seed_global_exact_family_impossible"] = 1
        counters["raw_seed_global_family_skips"] = len(
            grouped_clusters
        )
        counters["raw_seed_skipped_geometry_count"] = len(
            ordered_geometries
        )
        counters["raw_seed_skipped_geometry_sha256"] = _geometry_digest(
            ordered_geometries
        )
        counters["raw_seed_skipped_clustered_geometry_sha256"] = (
            _clustered_geometry_digest(
                ordered_geometries, grouped_clusters
            )
        )
        counters["raw_seed_prefilter_terminal_windows_tested"] = 0
        counters["raw_seed_parent_matches"] = 0
        counters["raw_seed_parent_match_count_status"] = "complete"
        return (), dict(counters)

    selected = []
    scratch_terminal_windows = 0
    for geometry_index, geometry in enumerate(ordered_geometries):
        parent = _candidate_from_intervals(source, geometry)
        validate_candidate_source_reconstruction(
            source,
            parent.list_repeats,
            parent.list_repeat_starts,
            parent.list_spacers,
        )
        evidence, failed = _failed_structural_predicates(
            structural_filter, parent
        )
        if not evidence.structural_valid or failed != ("minimum_repeat_count",):
            counters["raw_seed_structural_rejections"] += 1
            continue
        parent_interval = candidate_source_interval(parent)
        if any(
            _strict_overlap(parent_interval, record["interval"])
            for record in accepted
        ):
            counters["raw_seed_accepted_overlap_rejections"] += 1
            continue
        cluster_indices = tuple(sorted(grouped_clusters[geometry]))
        scratch = Counter()
        provisional = _template_satellite_proposals(
            source,
            [policy.template_parent_minimum_score, parent, {}],
            accepted,
            structural_filter,
            policy,
            scratch,
            parent_origin=RAW_TWO_REPEAT_SEED_PARENT,
            raw_seed_cluster_indices=cluster_indices,
            maximum_terminal_windows=(
                policy.maximum_raw_seed_terminal_windows_per_parent
            ),
            supported_template_families=supported_template_families,
        )
        scratch_terminal_windows += scratch["template_terminal_windows_tested"]
        if _select_unique_source_winner(provisional, scratch) is None:
            counters["raw_seed_exact_family_prefilter_rejections"] += 1
            continue
        selected.append(AcceptedFamilyRawSeedParent(
            candidate=parent,
            seed_repeat_intervals=geometry,
            cluster_indices=cluster_indices,
        ))
        if len(selected) > policy.maximum_raw_seed_parent_matches_per_source:
            counters["raw_seed_prefilter_terminal_windows_tested"] = (
                scratch_terminal_windows
            )
            counters["raw_seed_parent_matches"] = len(selected)
            counters["raw_seed_parent_matches_at_cap"] = len(selected)
            counters["raw_seed_parent_matches_is_lower_bound"] = 1
            counters["raw_seed_parent_match_count_status"] = (
                "lower_bound_cap_exhausted"
            )
            unexamined = ordered_geometries[geometry_index + 1:]
            counters["raw_seed_unexamined_geometry_count"] = len(
                unexamined
            )
            counters["raw_seed_unexamined_geometry_sha256"] = (
                _geometry_digest(unexamined)
            )
            counters["raw_seed_unexamined_clustered_geometry_sha256"] = (
                _clustered_geometry_digest(unexamined, grouped_clusters)
            )
            counters["raw_seed_parent_match_cap_exhausted"] = 1
            return (), dict(counters)

    counters["raw_seed_prefilter_terminal_windows_tested"] = (
        scratch_terminal_windows
    )
    counters["raw_seed_parent_matches"] = len(selected)
    counters["raw_seed_parent_match_count_status"] = "complete"
    if len(selected) > policy.maximum_raw_seed_parent_matches_per_source:
        counters["raw_seed_parent_match_cap_exhausted"] = 1
        return (), dict(counters)
    selected.sort(key=lambda item: item.seed_repeat_intervals)
    return tuple(selected), dict(counters)


def generate_accepted_family_rescue_proposals(
    source,
    categories,
    structural_filter,
    *,
    raw_seed_parent_rows=(),
    accepted_categories=None,
    policy=ACCEPTED_FAMILY_RESCUE_POLICY,
):
    """Return source-only family proposals plus deterministic telemetry."""

    if not isinstance(source, str) or not source:
        raise ValueError("family rescue source must be a non-empty string")
    source = source.upper()
    if not isinstance(policy, AcceptedFamilyRescuePolicy):
        raise TypeError("family rescue policy has the wrong type")
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("family rescue requires five category mappings")
    if accepted_categories is None:
        accepted_categories = categories
    if (
        not isinstance(accepted_categories, (list, tuple))
        or len(accepted_categories) != 5
    ):
        raise ValueError("family rescue accepted-anchor snapshot is malformed")
    if not callable(structural_filter) or not hasattr(structural_filter, "evidence"):
        raise TypeError("family rescue requires the canonical structural filter")

    accepted = _accepted_records(accepted_categories)
    supported_template_families = _supported_template_family_records(
        accepted, policy
    )
    has_terminal_family_anchor = any(
        record["category_index"] == 0 for record in accepted
    )
    low_rows = tuple(
        row
        for rows in categories[4].values()
        for row in rows
    )
    try:
        raw_seed_parent_rows = tuple(raw_seed_parent_rows)
    except TypeError as error:
        raise ValueError("raw-seed parent rows must be iterable") from error
    counters = Counter(
        low_rows=len(low_rows),
        raw_seed_parent_rows=len(raw_seed_parent_rows),
        accepted_rows=len(accepted),
        parent_row_cap=policy.maximum_parent_rows_per_source,
        proposal_cap=policy.maximum_proposals_per_source,
        supported_template_families=len(supported_template_families),
    )
    counters["search_implementation_version"] = (
        SEARCH_IMPLEMENTATION_VERSION
    )
    if (
        len(low_rows) + len(raw_seed_parent_rows)
        > policy.maximum_parent_rows_per_source
    ):
        counters["parent_row_cap_exhausted"] = 1
        return (), dict(counters)

    context_cache = {}
    winners = []
    for row in low_rows:
        if not isinstance(row, (tuple, list)) or len(row) < 2:
            raise ValueError("family rescue received a malformed Low row")
        counters["terminal_parents_considered"] += 1
        terminal = (
            _terminal_family_proposals(
                source,
                context_cache,
                row,
                accepted,
                structural_filter,
                policy,
                counters,
            )
            if has_terminal_family_anchor else ()
        )
        template = (
            _template_satellite_proposals(
                source,
                row,
                accepted,
                structural_filter,
                policy,
                counters,
                supported_template_families=(
                    supported_template_families
                ),
            )
            if supported_template_families else ()
        )
        counters["family_terminal_raw_proposals"] += len(terminal)
        counters["exact_template_raw_proposals"] += len(template)
        winner = _select_unique_source_winner(terminal + template, counters)
        if winner is not None:
            winners.append(winner)

    for raw_record in raw_seed_parent_rows:
        if (
            not isinstance(raw_record, (tuple, list))
            or len(raw_record) != 2
            or not isinstance(raw_record[1], AcceptedFamilyRawSeedParent)
        ):
            raise ValueError("raw-seed scored parent record is malformed")
        row, parent_record = raw_record
        if (
            not isinstance(row, (tuple, list))
            or len(row) < 2
            or row[1] is not parent_record.candidate
        ):
            raise ValueError("raw-seed parent scoring changed identity")
        counters["raw_seed_parents_considered"] += 1
        template = _template_satellite_proposals(
            source,
            row,
            accepted,
            structural_filter,
            policy,
            counters,
            parent_origin=RAW_TWO_REPEAT_SEED_PARENT,
            raw_seed_cluster_indices=parent_record.cluster_indices,
            maximum_terminal_windows=(
                policy.maximum_raw_seed_terminal_windows_per_parent
            ),
            supported_template_families=supported_template_families,
        )
        counters["exact_template_raw_proposals"] += len(template)
        winner = _select_unique_source_winner(template, counters)
        if winner is not None:
            winners.append(winner)

    # One admission per strictly overlapping unresolved source component.
    winners.sort(key=lambda item: (candidate_source_interval(item.candidate), item.arm))
    components = []
    current = []
    current_end = None
    for proposal in winners:
        start, end = candidate_source_interval(proposal.candidate)
        if not current or start < current_end:
            current.append(proposal)
            current_end = end if current_end is None else max(current_end, end)
        else:
            components.append(current)
            current = [proposal]
            current_end = end
    if current:
        components.append(current)
    selected = []
    for component in components:
        winner = _select_unique_source_winner(component, counters)
        if winner is not None:
            selected.append(winner)
    if len(selected) > policy.maximum_proposals_per_source:
        counters["proposal_cap_exhausted"] = 1
        return (), dict(counters)
    selected.sort(key=lambda item: candidate_source_interval(item.candidate))
    counters["source_components"] = len(components)
    counters["proposals"] = len(selected)
    counters["family_terminal_proposals"] = sum(
        item.arm == FAMILY_TERMINAL for item in selected
    )
    counters["exact_template_proposals"] = sum(
        item.arm == EXACT_TEMPLATE_SATELLITE for item in selected
    )
    counters["source_sha256"] = _sha256(source)
    return tuple(selected), dict(counters)


def completed_score_is_admissible(score, policy=ACCEPTED_FAMILY_RESCUE_POLICY):
    if (
        isinstance(score, bool)
        or not isinstance(score, (int, float))
        or not math.isfinite(float(score))
    ):
        return False
    return (
        policy.minimum_completed_score <= float(score)
        < policy.maximum_completed_score_exclusive
    )


def classification_basis(source, proposal, fresh_score, *, policy=ACCEPTED_FAMILY_RESCUE_POLICY):
    """Build the runtime-object-free basis attached to one Possible row."""

    source = source.upper()
    if not completed_score_is_admissible(fresh_score, policy):
        raise ValueError("family rescue fresh score is outside the policy band")
    if not isinstance(proposal, AcceptedFamilyRescueProposal):
        raise TypeError("family rescue proposal has the wrong type")
    return {
        "method": CLASSIFICATION_METHOD,
        "classification_ceiling": "Possible",
        "boundary_status": "AMBIGUOUS",
        "exact_boundary_certified": False,
        "accepted_family_rescue": proposal.as_evidence(
            _sha256(source), fresh_score, policy
        ),
    }


def _canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _record_intervals(value, source_length, label):
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError("{} must be a non-empty interval sequence".format(label))
    intervals = []
    for index, raw in enumerate(value):
        if not isinstance(raw, (list, tuple)) or len(raw) != 2:
            raise ValueError("{}[{}] is malformed".format(label, index))
        start, end = raw
        if (
            not isinstance(start, int)
            or isinstance(start, bool)
            or not isinstance(end, int)
            or isinstance(end, bool)
            or not 0 <= start < end <= source_length
        ):
            raise ValueError("{}[{}] is outside source".format(label, index))
        intervals.append((start, end))
    if any(left[1] > right[0] for left, right in zip(intervals, intervals[1:])):
        raise ValueError("{} is not monotonic".format(label))
    return tuple(intervals)


def _same_candidate(first, second):
    return (
        _geometry(first) == _geometry(second)
        and tuple(first.list_repeats) == tuple(second.list_repeats)
        and tuple(first.list_spacers) == tuple(second.list_spacers)
        and first.consensus == second.consensus
    )


def replay_classification_evidence(source, candidate, evidence, row_score):
    """Recompute one serialized family-rescue record from source DNA.

    This validator intentionally does not trust a provenance boolean.  It
    rebuilds parent, anchors, and proposal, then replays the family relation,
    terminal placement, strict-existence decision, and all source guards.
    The accepted category indexes/scores are transaction attestations: report
    validation can source-replay those anchor geometries but does not rerun the
    compact model for the historical anchor rows.
    """

    if not isinstance(source, str) or not source:
        raise ValueError("family rescue replay source is invalid")
    source = source.upper()
    if not isinstance(evidence, dict):
        raise ValueError("accepted_family_rescue must be a mapping")
    common_keys = {
        "arm",
        "coordinate_system",
        "source_sha256",
        "parent_score",
        "parent_repeat_intervals",
        "parent_origin",
        "raw_seed_cluster_indices",
        "family_template",
        "accepted_family_anchors",
        "phase_offset",
        "terminal_side",
        "terminal_start",
        "proposal_repeat_intervals",
        "source_metrics",
        "fresh_score",
        "policy",
    }
    arm = evidence.get("arm")
    required = set(common_keys)
    if arm == FAMILY_TERMINAL:
        required.update((
            "base_source_intact_decision",
            "completed_source_intact_decision",
        ))
    elif arm != EXACT_TEMPLATE_SATELLITE:
        raise ValueError("family rescue arm is unsupported")
    if set(evidence) != required:
        raise ValueError("family rescue evidence keys are stale")
    if evidence["coordinate_system"] != COORDINATE_SYSTEM:
        raise ValueError("family rescue coordinate system is stale")
    source_hash = _sha256(source)
    if evidence["source_sha256"] != source_hash:
        raise ValueError("family rescue source hash is stale")
    policy_mapping = evidence["policy"]
    if not isinstance(policy_mapping, dict):
        raise ValueError("family rescue policy is malformed")
    try:
        policy = AcceptedFamilyRescuePolicy(**policy_mapping)
    except (TypeError, ValueError) as error:
        raise ValueError("family rescue policy is invalid") from error
    if policy.as_dict() != ACCEPTED_FAMILY_RESCUE_POLICY.as_dict():
        raise ValueError("family rescue policy is not the supported version")
    fresh_score = evidence["fresh_score"]
    if (
        isinstance(row_score, bool)
        or not isinstance(row_score, (int, float))
        or float(fresh_score) != float(row_score)
        or not completed_score_is_admissible(float(row_score), policy)
    ):
        raise ValueError("family rescue fresh score is stale")

    proposal_intervals = _record_intervals(
        evidence["proposal_repeat_intervals"],
        len(source),
        "proposal_repeat_intervals",
    )
    rebuilt = _candidate_from_intervals(source, proposal_intervals)
    if not _same_candidate(candidate, rebuilt):
        raise ValueError("family rescue proposal differs from reported candidate")
    parent_intervals = _record_intervals(
        evidence["parent_repeat_intervals"],
        len(source),
        "parent_repeat_intervals",
    )
    parent = _candidate_from_intervals(source, parent_intervals)
    parent_score = evidence["parent_score"]
    if (
        isinstance(parent_score, bool)
        or not isinstance(parent_score, (int, float))
        or not math.isfinite(float(parent_score))
    ):
        raise ValueError("family rescue parent score is invalid")
    parent_origin = evidence["parent_origin"]
    raw_seed_cluster_indices = evidence["raw_seed_cluster_indices"]
    if parent_origin not in (LOW_SCORE_PARENT, RAW_TWO_REPEAT_SEED_PARENT):
        raise ValueError("family rescue parent origin is invalid")
    if not isinstance(raw_seed_cluster_indices, list) or any(
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
        for value in raw_seed_cluster_indices
    ):
        raise ValueError("family rescue raw-seed clusters are malformed")
    if raw_seed_cluster_indices != sorted(set(raw_seed_cluster_indices)):
        raise ValueError("family rescue raw-seed clusters are not canonical")
    if (
        (parent_origin == LOW_SCORE_PARENT and raw_seed_cluster_indices)
        or (
            parent_origin == RAW_TWO_REPEAT_SEED_PARENT
            and not raw_seed_cluster_indices
        )
    ):
        raise ValueError("family rescue parent origin evidence is inconsistent")
    template = evidence["family_template"]
    if not isinstance(template, str) or not template:
        raise ValueError("family rescue template is invalid")
    anchors_raw = evidence["accepted_family_anchors"]
    if not isinstance(anchors_raw, list) or not anchors_raw:
        raise ValueError("family rescue anchors are missing")
    anchors = []
    seen_anchor_intervals = set()
    for index, raw in enumerate(anchors_raw):
        if not isinstance(raw, dict) or set(raw) != {
            "category_index", "score", "repeat_intervals", "consensus"
        }:
            raise ValueError("family rescue anchor is malformed")
        if raw["category_index"] not in (0, 2):
            raise ValueError("family rescue anchor category is not accepted")
        if (
            isinstance(raw["score"], bool)
            or not isinstance(raw["score"], (int, float))
            or not math.isfinite(float(raw["score"]))
        ):
            raise ValueError("family rescue anchor score is invalid")
        intervals = _record_intervals(
            raw["repeat_intervals"],
            len(source),
            "accepted_family_anchors[{}].repeat_intervals".format(index),
        )
        anchor_candidate = _candidate_from_intervals(source, intervals)
        if anchor_candidate.consensus != raw["consensus"]:
            raise ValueError("family rescue anchor consensus is stale")
        interval = candidate_source_interval(anchor_candidate)
        if interval in seen_anchor_intervals:
            raise ValueError("family rescue anchor interval is duplicated")
        seen_anchor_intervals.add(interval)
        anchors.append((raw, anchor_candidate, interval))

    side = evidence["terminal_side"]
    terminal_start = evidence["terminal_start"]
    if side not in ("left", "right"):
        raise ValueError("family rescue terminal side is invalid")
    if not isinstance(terminal_start, int) or isinstance(terminal_start, bool):
        raise ValueError("family rescue terminal start is invalid")
    parent_interval = candidate_source_interval(parent)
    terminal_index = 0 if side == "left" else len(rebuilt.list_repeats) - 1
    terminal_repeat = rebuilt.list_repeats[terminal_index]
    if rebuilt.list_repeat_starts[terminal_index] != terminal_start:
        raise ValueError("family rescue terminal start differs from proposal")

    expected_metrics = _family_guard_metrics(source, rebuilt, template)
    expected_metrics["terminal_identity"] = _edit_identity(
        terminal_repeat, template
    )
    expected_metrics["terminal_exact_anchor"] = _longest_exact_source_anchor(
        terminal_repeat, template
    )
    minimum_gap = min(
        _source_gap(parent_interval, anchor_interval)
        for _raw, _anchor_candidate, anchor_interval in anchors
    )
    expected_metrics["accepted_family_minimum_gap"] = minimum_gap

    if arm == FAMILY_TERMINAL:
        if not (
            parent_origin == LOW_SCORE_PARENT
            and not raw_seed_cluster_indices
            and
            policy.terminal_parent_minimum_score <= float(parent_score)
            < policy.maximum_completed_score_exclusive
            and policy.terminal_parent_minimum_repeats
            <= len(parent.list_repeats)
            <= policy.terminal_parent_maximum_repeats
            and len(anchors) == 1
            and anchors[0][0]["category_index"] == 0
            and minimum_gap <= policy.terminal_family_maximum_gap
        ):
            raise ValueError("family-terminal parent/anchor contract is stale")
        anchor_consensus = anchors[0][1].consensus
        if template not in (anchor_consensus, _reverse_complement(anchor_consensus)):
            raise ValueError("family-terminal template is not its anchor family")
        delta = len(template) - len(parent.consensus)
        if not (
            policy.terminal_family_minimum_length_delta <= delta
            <= policy.terminal_family_maximum_length_delta
        ):
            raise ValueError("family-terminal phase length delta is invalid")
        alignment = _best_unique_local_alignment(
            parent.consensus,
            template,
            policy.terminal_family_minimum_local_identity,
        )
        if alignment is None or evidence["phase_offset"] != alignment[1]:
            raise ValueError("family-terminal phase alignment is stale")
        expected_metrics["parent_family_local_identity"] = alignment[0]
        length = len(template)
        base_intervals = tuple(
            (start - alignment[1], start - alignment[1] + length)
            for start in parent.list_repeat_starts
        )
        base = _candidate_from_intervals(source, base_intervals)
        expected_proposal = (
            base_intervals + ((terminal_start, terminal_start + length),)
            if side == "right"
            else ((terminal_start, terminal_start + length),) + base_intervals
        )
        if proposal_intervals != expected_proposal:
            raise ValueError("family-terminal proposal is not a one-terminal phase completion")
        context = SourceIntactEvaluationContext(source)
        base_decision = evaluate_source_intact_existence_in_context(context, base)
        completed_decision = evaluate_source_intact_existence_in_context(
            context, rebuilt
        )
        if (
            tuple(reason.code for reason in base_decision.reasons)
            != ("spacer_length_cv_above_maximum",)
            or completed_decision.decision != ELIGIBLE
            or _canonical_json(evidence["base_source_intact_decision"])
            != _canonical_json(base_decision.as_dict())
            or _canonical_json(evidence["completed_source_intact_decision"])
            != _canonical_json(completed_decision.as_dict())
        ):
            raise ValueError("family-terminal strict source-intact replay failed")
        if expected_metrics["terminal_identity"] < policy.terminal_minimum_identity:
            raise ValueError("family-terminal terminal identity is too weak")
    else:
        anchor_intervals = tuple(item[2] for item in anchors)
        anchors_are_nonoverlapping = not any(
            _strict_overlap(first, second)
            for index, first in enumerate(anchor_intervals[:-1])
            for second in anchor_intervals[index + 1:]
        )
        anchors_have_exact_support = all(
            sum(
                repeat == anchor.consensus
                for repeat in anchor.list_repeats
            )
            >= policy.template_minimum_exact_repeat_support_per_array
            for _raw, anchor, _interval in anchors
        )
        if not (
            policy.template_parent_minimum_score <= float(parent_score)
            < policy.maximum_completed_score_exclusive
            and len(parent.list_repeats) == 2
            and len(anchors) >= policy.template_minimum_distinct_accepted_arrays
            and anchors_are_nonoverlapping
            and anchors_have_exact_support
            and minimum_gap <= policy.template_family_maximum_gap
            and policy.template_minimum_length <= len(template)
            <= policy.template_maximum_length
        ):
            raise ValueError("exact-template parent/anchor contract is stale")
        canonical = _canonical_repeat(template)
        if any(_canonical_repeat(anchor.consensus) != canonical for _r, anchor, _i in anchors):
            raise ValueError("exact-template anchors do not share one family")
        offsets = []
        for repeat in parent.list_repeats:
            offset = repeat.find(template)
            if offset < 0 or repeat.find(template, offset + 1) >= 0:
                raise ValueError("exact-template parent containment is ambiguous")
            offsets.append(offset)
        try:
            recorded_offsets = tuple(evidence["phase_offset"])
        except TypeError as error:
            raise ValueError("exact-template offsets are malformed") from error
        if recorded_offsets != tuple(offsets):
            raise ValueError("exact-template independent offsets are stale")
        length = len(template)
        core = tuple(
            (start + offset, start + offset + length)
            for start, offset in zip(parent.list_repeat_starts, offsets)
        )
        expected_proposal = (
            core + ((terminal_start, terminal_start + length),)
            if side == "right"
            else ((terminal_start, terminal_start + length),) + core
        )
        if proposal_intervals != expected_proposal:
            raise ValueError("exact-template proposal geometry is stale")
        first_spacer = core[1][0] - core[0][1]
        terminal_spacer = (
            terminal_start - core[-1][1]
            if side == "right"
            else core[0][0] - (terminal_start + length)
        )
        if not (
            policy.minimum_spacer_length <= terminal_spacer
            <= policy.maximum_spacer_length
            and abs(terminal_spacer - first_spacer)
            <= policy.template_spacer_radius
            and expected_metrics["terminal_identity"]
            >= policy.template_terminal_minimum_identity
            and expected_metrics["terminal_exact_anchor"]
            >= policy.template_terminal_minimum_exact_anchor
        ):
            raise ValueError("exact-template terminal cadence/evidence is stale")
        expected_metrics["accepted_family_array_count"] = len(anchors)
        expected_metrics["template_offsets"] = tuple(offsets)
        if (
            expected_metrics["minimum_repeat_identity"]
            < policy.template_minimum_repeat_identity
            or expected_metrics["mean_repeat_identity"]
            < policy.template_minimum_mean_repeat_identity
            or expected_metrics["spacer_cv"] > policy.template_maximum_spacer_cv
            or expected_metrics["maximum_spacer_identity"]
            > policy.template_maximum_spacer_identity
        ):
            raise ValueError("exact-template source morphology is stale")

    if not _common_guards_pass(expected_metrics, policy):
        raise ValueError("family rescue common source guards reject")
    if _canonical_json(evidence["source_metrics"]) != _canonical_json(
        _json_metrics(expected_metrics)
    ):
        raise ValueError("family rescue source metrics are stale")
    return evidence


__all__ = (
    "ACCEPTED_FAMILY_RESCUE_POLICY",
    "AcceptedFamilyRawSeedParent",
    "AcceptedFamilyRescuePolicy",
    "AcceptedFamilyRescueProposal",
    "CLASSIFICATION_METHOD",
    "EXACT_TEMPLATE_SATELLITE",
    "FAMILY_TERMINAL",
    "LOW_SCORE_PARENT",
    "RAW_TWO_REPEAT_SEED_PARENT",
    "classification_basis",
    "completed_score_is_admissible",
    "generate_accepted_family_rescue_proposals",
    "replay_classification_evidence",
    "select_exact_template_raw_seed_parents",
)
