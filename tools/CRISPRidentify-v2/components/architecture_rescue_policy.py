"""Conservative, source-backed architecture-rescue eligibility policy.

This module evaluates an already enhanced CRISPR candidate.  It deliberately
does not enhance, phase-collapse, classify, categorize, promote, remove, or
mutate candidates.  The caller supplies the classifier score and named
features; this module independently reconstructs the candidate from the
original sequence and applies the frozen-300 audit gates to that final
geometry.

Optional coding-context inputs are paired.  Supplying neither member means
that evidence was not assessed and does not block eligibility.  Supplying
only one member forces review.  A complete coding pair is not a standalone
veto: strong CDS coverage must agree with both explicit in-frame repeat
periodicity and tandem-HMM support.  This avoids rejecting genuine arrays
merely because an ORF caller spans them or their cadence happens to be
divisible by three.
"""

from collections import Counter
from dataclasses import asdict, dataclass, field
import math
from numbers import Real
from typing import Any, Dict, Optional, Tuple

from components.reconstruction_validator import (
    ReconstructionResult,
    validate_candidate_reconstruction,
)


ELIGIBLE = "eligible"
REVIEW = "review"
INELIGIBLE = "ineligible"

REJECT_REASON = "reject"
REVIEW_REASON = "review"


@dataclass(frozen=True)
class ArchitectureRescuePolicy:
    """All thresholds used by the audited post-enhancement policy."""

    name: str = "frozen300-post-enhancement-architecture-rescue"
    version: str = "1.1.0"
    min_repeat_count: int = 5
    min_classifier_score: float = 0.40
    min_consensus_length: int = 20
    max_consensus_length: int = 55
    min_spacer_length: int = 18
    max_spacer_length: int = 80
    max_deleted_spacer_count: int = 0
    max_spacer_length_cv: float = 0.10
    min_repeat_similarity: float = 0.70
    max_spacer_similarity: float = 0.50
    max_distinct_spacer_edit_identity: float = 0.80
    duplicate_pair_min_spacer_count: int = 7
    duplicate_pair_min_other_distinct_spacers: int = 5
    terminal_kmer_length: int = 5
    terminal_kmer_review_dominance: float = 0.80
    coding_overlap_veto_fraction: float = 0.80
    context_cds_score_veto_minimum: float = 20.0

    def __post_init__(self):
        def require_nonempty_string(attribute):
            value = getattr(self, attribute)
            if not isinstance(value, str) or not value:
                raise ValueError(
                    "{} must be a non-empty string".format(attribute)
                )

        def require_integer(attribute, minimum):
            value = getattr(self, attribute)
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < minimum
            ):
                raise ValueError(
                    "{} must be an integer >= {}".format(attribute, minimum)
                )

        def require_probability(attribute):
            value = getattr(self, attribute)
            if (
                not isinstance(value, Real)
                or isinstance(value, bool)
                or not math.isfinite(float(value))
                or not 0.0 <= float(value) <= 1.0
            ):
                raise ValueError(
                    "{} must be a finite number in [0, 1]".format(attribute)
                )

        require_nonempty_string("name")
        require_nonempty_string("version")
        require_integer("min_repeat_count", 2)
        require_integer("min_consensus_length", 1)
        require_integer("max_consensus_length", 1)
        require_integer("min_spacer_length", 1)
        require_integer("max_spacer_length", 1)
        require_integer("max_deleted_spacer_count", 0)
        require_integer("duplicate_pair_min_spacer_count", 2)
        require_integer("duplicate_pair_min_other_distinct_spacers", 1)
        require_integer("terminal_kmer_length", 1)

        for attribute in (
            "min_classifier_score",
            "max_spacer_length_cv",
            "min_repeat_similarity",
            "max_spacer_similarity",
            "max_distinct_spacer_edit_identity",
            "terminal_kmer_review_dominance",
            "coding_overlap_veto_fraction",
        ):
            require_probability(attribute)

        if (
            not isinstance(self.context_cds_score_veto_minimum, Real)
            or isinstance(self.context_cds_score_veto_minimum, bool)
            or not math.isfinite(float(self.context_cds_score_veto_minimum))
            or self.context_cds_score_veto_minimum < 0.0
        ):
            raise ValueError(
                "context_cds_score_veto_minimum must be a finite "
                "non-negative number"
            )

        if self.min_consensus_length > self.max_consensus_length:
            raise ValueError(
                "min_consensus_length must not exceed max_consensus_length"
            )
        if self.min_spacer_length > self.max_spacer_length:
            raise ValueError(
                "min_spacer_length must not exceed max_spacer_length"
            )
        if self.terminal_kmer_length > self.min_spacer_length:
            raise ValueError(
                "terminal_kmer_length must not exceed min_spacer_length"
            )
        minimum_duplicate_geometry = (
            self.duplicate_pair_min_other_distinct_spacers + 2
        )
        if (
            self.duplicate_pair_min_spacer_count
            < minimum_duplicate_geometry
        ):
            raise ValueError(
                "duplicate_pair_min_spacer_count must allow the duplicate "
                "pair plus every required other distinct spacer"
            )

    def as_dict(self):
        return asdict(self)


AUDITED_ARCHITECTURE_RESCUE_POLICY = ArchitectureRescuePolicy()


@dataclass(frozen=True)
class ArchitectureRescueReason:
    """One stable, machine-readable rejection or review reason."""

    code: str
    severity: str
    message: str
    details: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self):
        value = {
            "code": self.code,
            "severity": self.severity,
            "message": self.message,
        }
        if self.details:
            value["details"] = dict(self.details)
        return value


@dataclass(frozen=True)
class ArchitectureRescueStats:
    """Geometry and evidence values measured for one final candidate."""

    reconstruction_valid: bool
    reconstruction_issue_codes: Tuple[str, ...]
    array_start: Optional[int]
    array_end: Optional[int]
    repeat_count: int
    spacer_count: int
    nondeleted_spacer_count: int
    deleted_spacer_count: int
    consensus_length: Optional[int]
    nondeleted_spacer_lengths: Tuple[int, ...]
    mean_spacer_length: Optional[float]
    spacer_length_cv: Optional[float]
    unique_spacer_count: int
    exact_duplicate_pair_count: int
    max_distinct_spacer_edit_identity: Optional[float]
    max_distinct_spacer_pair: Optional[Tuple[str, str]]
    dominant_prefix_5mer: Optional[str]
    prefix_5mer_dominance: Optional[float]
    dominant_suffix_5mer: Optional[str]
    suffix_5mer_dominance: Optional[float]
    classifier_score: Optional[float]
    repeat_similarity: Optional[float]
    spacer_similarity: Optional[float]
    afsf_pass: Optional[bool]
    coding_veto_assessed: bool
    coding_veto_triggered: bool
    tandem_veto_assessed: bool
    tandem_veto_triggered: bool

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class ArchitectureRescueDecision:
    """Structured result; it carries no pipeline category or mutation."""

    decision: str
    reasons: Tuple[ArchitectureRescueReason, ...]
    stats: ArchitectureRescueStats
    reconstruction: Optional[ReconstructionResult]
    policy: ArchitectureRescuePolicy

    @property
    def eligible(self):
        return self.decision == ELIGIBLE

    @property
    def requires_review(self):
        return self.decision == REVIEW

    @property
    def rejection_codes(self):
        return tuple(
            reason.code
            for reason in self.reasons
            if reason.severity == REJECT_REASON
        )

    @property
    def review_codes(self):
        return tuple(
            reason.code
            for reason in self.reasons
            if reason.severity == REVIEW_REASON
        )

    def as_dict(self):
        return {
            "decision": self.decision,
            "eligible": self.eligible,
            "requires_review": self.requires_review,
            "rejection_codes": list(self.rejection_codes),
            "review_codes": list(self.review_codes),
            "reasons": [reason.as_dict() for reason in self.reasons],
            "stats": self.stats.as_dict(),
            "reconstruction": (
                self.reconstruction.as_dict()
                if self.reconstruction is not None
                else None
            ),
            "policy": self.policy.as_dict(),
        }


def _candidate_tuple(candidate, attribute):
    value = getattr(candidate, attribute, ())
    if value is None or isinstance(value, (str, bytes)):
        return ()
    try:
        return tuple(value)
    except TypeError:
        return ()


def _finite_real(value):
    return (
        isinstance(value, Real)
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _unit_interval(value):
    return _finite_real(value) and 0.0 <= float(value) <= 1.0


def _population_cv(lengths):
    if not lengths:
        return None, None
    mean = sum(lengths) / float(len(lengths))
    if mean == 0:
        return mean, None
    variance = sum((length - mean) ** 2 for length in lengths) / len(lengths)
    return mean, math.sqrt(variance) / mean


def _global_edit_distance(first, second):
    """Return Levenshtein distance with full-sequence end-gap penalties."""

    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for first_index, first_base in enumerate(first, 1):
        current = [first_index]
        for second_index, second_base in enumerate(second, 1):
            current.append(min(
                current[-1] + 1,
                previous[second_index] + 1,
                previous[second_index - 1] + (first_base != second_base),
            ))
        previous = current
    return previous[-1]


def _global_edit_identity(first, second):
    denominator = max(len(first), len(second))
    if denominator == 0:
        return 1.0
    return 1.0 - _global_edit_distance(first, second) / float(denominator)


def _maximum_distinct_spacer_identity(spacers):
    distinct = sorted(set(spacers))
    if len(distinct) < 2:
        return None, None

    best_identity = -1.0
    best_pair = None
    for left_index, left in enumerate(distinct[:-1]):
        for right in distinct[left_index + 1:]:
            identity = _global_edit_identity(left, right)
            if identity > best_identity:
                best_identity = identity
                best_pair = (left, right)
    return best_identity, best_pair


def _dominant_kmer(spacers, kmer_length, from_end=False):
    if not spacers:
        return None, None
    if from_end:
        kmers = [spacer[-kmer_length:] for spacer in spacers]
    else:
        kmers = [spacer[:kmer_length] for spacer in spacers]
    counts = Counter(kmers)
    dominant, count = sorted(
        counts.items(),
        key=lambda item: (-item[1], item[0]),
    )[0]
    return dominant, count / float(len(spacers))


def _reason(reasons, code, severity, message, **details):
    reasons.append(ArchitectureRescueReason(
        code=code,
        severity=severity,
        message=message,
        details=details,
    ))


def evaluate_architecture_rescue_eligibility(
    genome_sequence,
    candidate,
    *,
    classifier_score,
    repeat_similarity,
    spacer_similarity,
    afsf_pass,
    coding_overlap_fraction=None,
    context_cds_score=None,
    in_frame_repeat_periodicity=None,
    tandem_hmm_support=None,
    coordinate_offset=0,
    policy=AUDITED_ARCHITECTURE_RESCUE_POLICY,
):
    """Evaluate automatic rescue eligibility for one enhanced candidate.

    Similarities and the classifier score are supplied explicitly because the
    policy must judge the final, post-enhancement feature vector.  Empty spacer
    strings are treated as deleted spacers: they remain part of exact source
    reconstruction, are checked against ``max_deleted_spacer_count``, and are
    excluded from length, CV, diversity, and terminal k-mer statistics.

    Global spacer identity is ``1 - Levenshtein/max(lengths)`` and therefore
    penalizes all end gaps.  Exact duplicates are handled separately: at most
    one sequence may occur twice, and that exception requires at least seven
    nondeleted spacers plus five *other* distinct spacers.
    """

    if not isinstance(policy, ArchitectureRescuePolicy):
        raise TypeError("policy must be an ArchitectureRescuePolicy")

    reasons = []
    repeats = _candidate_tuple(candidate, "list_repeats")
    spacers = _candidate_tuple(candidate, "list_spacers")
    consensus = getattr(candidate, "consensus", None)

    reconstruction = None
    reconstruction_error = None
    try:
        reconstruction = validate_candidate_reconstruction(
            genome_sequence,
            candidate,
            coordinate_offset=coordinate_offset,
        )
    except (AttributeError, TypeError, ValueError) as error:
        reconstruction_error = str(error)

    if reconstruction is None:
        reconstruction_valid = False
        reconstruction_issue_codes = ("validation_error",)
        array_start = None
        array_end = None
        _reason(
            reasons,
            "source_reconstruction_failed",
            REJECT_REASON,
            "candidate could not be reconstructed from the original sequence",
            error=reconstruction_error,
        )
    else:
        reconstruction_valid = reconstruction.valid
        reconstruction_issue_codes = tuple(
            issue.code for issue in reconstruction.issues
        )
        array_start = reconstruction.array_start
        array_end = reconstruction.array_end
        if not reconstruction_valid:
            _reason(
                reasons,
                "source_reconstruction_failed",
                REJECT_REASON,
                "candidate sequences or coordinates do not reconstruct the original sequence",
                issue_codes=reconstruction_issue_codes,
            )

    normalized_spacers = tuple(
        spacer.upper()
        for spacer in spacers
        if isinstance(spacer, str) and spacer
    )
    deleted_spacer_count = sum(
        isinstance(spacer, str) and not spacer for spacer in spacers
    )
    spacer_lengths = tuple(len(spacer) for spacer in normalized_spacers)
    mean_spacer_length, spacer_length_cv = _population_cv(spacer_lengths)
    counts = Counter(normalized_spacers)
    duplicate_pair_count = sum(
        count * (count - 1) // 2 for count in counts.values()
    )
    max_edit_identity, max_edit_pair = _maximum_distinct_spacer_identity(
        normalized_spacers
    )
    prefix_kmer, prefix_dominance = _dominant_kmer(
        normalized_spacers,
        policy.terminal_kmer_length,
    )
    suffix_kmer, suffix_dominance = _dominant_kmer(
        normalized_spacers,
        policy.terminal_kmer_length,
        from_end=True,
    )

    score_value = float(classifier_score) if _finite_real(classifier_score) else None
    repeat_similarity_value = (
        float(repeat_similarity) if _finite_real(repeat_similarity) else None
    )
    spacer_similarity_value = (
        float(spacer_similarity) if _finite_real(spacer_similarity) else None
    )
    afsf_value = afsf_pass if isinstance(afsf_pass, bool) else None

    if len(repeats) < policy.min_repeat_count:
        _reason(
            reasons,
            "repeat_count_below_minimum",
            REJECT_REASON,
            "candidate has fewer repeats than the rescue minimum",
            observed=len(repeats),
            minimum=policy.min_repeat_count,
        )

    if not _unit_interval(classifier_score):
        _reason(
            reasons,
            "invalid_classifier_score",
            REJECT_REASON,
            "classifier_score must be a finite number in [0, 1]",
            observed=repr(classifier_score),
        )
    elif score_value < policy.min_classifier_score:
        _reason(
            reasons,
            "classifier_score_below_minimum",
            REJECT_REASON,
            "classifier score is below the rescue minimum",
            observed=score_value,
            minimum=policy.min_classifier_score,
        )

    consensus_length = len(consensus) if isinstance(consensus, str) else None
    if consensus_length is None:
        _reason(
            reasons,
            "invalid_repeat_consensus",
            REJECT_REASON,
            "candidate consensus must be a string",
            observed=repr(consensus),
        )
    elif not (
        policy.min_consensus_length
        <= consensus_length
        <= policy.max_consensus_length
    ):
        _reason(
            reasons,
            "consensus_length_out_of_range",
            REJECT_REASON,
            "repeat consensus length is outside the rescue interval",
            observed=consensus_length,
            minimum=policy.min_consensus_length,
            maximum=policy.max_consensus_length,
        )

    if not normalized_spacers:
        _reason(
            reasons,
            "no_nondeleted_spacers",
            REJECT_REASON,
            "candidate has no nondeleted spacers to assess",
        )
    else:
        invalid_spacer_indexes = tuple(
            index
            for index, spacer in enumerate(spacers)
            if (
                isinstance(spacer, str)
                and spacer
                and not (
                    policy.min_spacer_length
                    <= len(spacer)
                    <= policy.max_spacer_length
                )
            )
        )
        if invalid_spacer_indexes:
            _reason(
                reasons,
                "spacer_length_out_of_range",
                REJECT_REASON,
                "one or more nondeleted spacers are outside the rescue interval",
                spacer_indexes=invalid_spacer_indexes,
                observed_lengths=tuple(
                    len(spacers[index]) for index in invalid_spacer_indexes
                ),
                minimum=policy.min_spacer_length,
                maximum=policy.max_spacer_length,
            )
        if (
            spacer_length_cv is not None
            and spacer_length_cv > policy.max_spacer_length_cv
        ):
            _reason(
                reasons,
                "spacer_length_cv_above_maximum",
                REJECT_REASON,
                "population spacer-length CV exceeds the rescue maximum",
                observed=spacer_length_cv,
                maximum=policy.max_spacer_length_cv,
            )

    if deleted_spacer_count > policy.max_deleted_spacer_count:
        _reason(
            reasons,
            "spacer_deletion_count_above_maximum",
            REJECT_REASON,
            "candidate has more deleted spacers than the rescue policy allows",
            observed=deleted_spacer_count,
            maximum=policy.max_deleted_spacer_count,
        )

    if not _unit_interval(repeat_similarity):
        _reason(
            reasons,
            "invalid_repeat_similarity",
            REJECT_REASON,
            "repeat_similarity must be a finite number in [0, 1]",
            observed=repr(repeat_similarity),
        )
    elif repeat_similarity_value < policy.min_repeat_similarity:
        _reason(
            reasons,
            "repeat_similarity_below_minimum",
            REJECT_REASON,
            "repeat similarity is below the rescue minimum",
            observed=repeat_similarity_value,
            minimum=policy.min_repeat_similarity,
        )

    if not _unit_interval(spacer_similarity):
        _reason(
            reasons,
            "invalid_spacer_similarity",
            REJECT_REASON,
            "spacer_similarity must be a finite number in [0, 1]",
            observed=repr(spacer_similarity),
        )
    elif spacer_similarity_value > policy.max_spacer_similarity:
        _reason(
            reasons,
            "spacer_similarity_above_maximum",
            REJECT_REASON,
            "spacer similarity exceeds the rescue maximum",
            observed=spacer_similarity_value,
            maximum=policy.max_spacer_similarity,
        )

    if not isinstance(afsf_pass, bool):
        _reason(
            reasons,
            "invalid_afsf_pass",
            REJECT_REASON,
            "afsf_pass must be an explicit boolean",
            observed=repr(afsf_pass),
        )
    elif not afsf_pass:
        _reason(
            reasons,
            "afsf_failed",
            REJECT_REASON,
            "candidate did not pass the supplied AFSF result",
        )

    duplicate_multiplicities = sorted(
        count for count in counts.values() if count > 1
    )
    duplicate_rule_passes = not duplicate_multiplicities
    if duplicate_multiplicities == [2]:
        duplicate_rule_passes = (
            len(normalized_spacers)
            >= policy.duplicate_pair_min_spacer_count
            and len(counts) - 1
            >= policy.duplicate_pair_min_other_distinct_spacers
        )
    if not duplicate_rule_passes:
        _reason(
            reasons,
            "exact_spacer_duplicate_rule_failed",
            REJECT_REASON,
            "exact spacer duplication exceeds the audited one-pair exception",
            duplicate_multiplicities=tuple(duplicate_multiplicities),
            nondeleted_spacer_count=len(normalized_spacers),
            other_distinct_spacer_count=(
                len(counts) - 1 if duplicate_multiplicities == [2] else None
            ),
            required_spacer_count=policy.duplicate_pair_min_spacer_count,
            required_other_distinct_spacers=(
                policy.duplicate_pair_min_other_distinct_spacers
            ),
        )

    if (
        max_edit_identity is not None
        and max_edit_identity > policy.max_distinct_spacer_edit_identity
    ):
        _reason(
            reasons,
            "distinct_spacer_identity_above_maximum",
            REJECT_REASON,
            "maximum global edit identity among distinct spacers is too high",
            observed=max_edit_identity,
            maximum=policy.max_distinct_spacer_edit_identity,
            spacer_pair=max_edit_pair,
        )

    terminal_dominance = max(
        value
        for value in (prefix_dominance, suffix_dominance)
        if value is not None
    ) if normalized_spacers else None
    if (
        terminal_dominance is not None
        and terminal_dominance >= policy.terminal_kmer_review_dominance
    ):
        _reason(
            reasons,
            "terminal_5mer_phase_review",
            REVIEW_REASON,
            "a shared terminal 5-mer requires source-phase review",
            dominant_prefix=prefix_kmer,
            prefix_dominance=prefix_dominance,
            dominant_suffix=suffix_kmer,
            suffix_dominance=suffix_dominance,
            review_threshold=policy.terminal_kmer_review_dominance,
        )

    coding_veto_assessed = False
    coding_veto_triggered = False
    coding_pair_supplied = (
        coding_overlap_fraction is not None or context_cds_score is not None
    )
    if coding_pair_supplied and (
        coding_overlap_fraction is None or context_cds_score is None
    ):
        _reason(
            reasons,
            "incomplete_coding_veto_inputs",
            REVIEW_REASON,
            "coding overlap fraction and context CDS score must be supplied together",
        )
    elif coding_pair_supplied:
        coding_overlap_valid = _unit_interval(coding_overlap_fraction)
        context_score_valid = (
            _finite_real(context_cds_score) and float(context_cds_score) >= 0.0
        )
        if not coding_overlap_valid:
            _reason(
                reasons,
                "invalid_coding_overlap_fraction",
                REJECT_REASON,
                "coding_overlap_fraction must be a finite number in [0, 1]",
                observed=repr(coding_overlap_fraction),
            )
        if not context_score_valid:
            _reason(
                reasons,
                "invalid_context_cds_score",
                REJECT_REASON,
                "context_cds_score must be a finite non-negative number",
                observed=repr(context_cds_score),
            )
        if coding_overlap_valid and context_score_valid:
            coding_veto_assessed = True
            coding_context_signal = (
                float(coding_overlap_fraction)
                >= policy.coding_overlap_veto_fraction
                and float(context_cds_score)
                >= policy.context_cds_score_veto_minimum
            )
            coding_veto_triggered = (
                coding_context_signal
                and in_frame_repeat_periodicity is True
                and tandem_hmm_support is True
            )
            if coding_veto_triggered:
                _reason(
                    reasons,
                    "coding_context_veto",
                    REJECT_REASON,
                    "context CDS, in-frame periodicity, and tandem-HMM evidence agree",
                    coding_overlap_fraction=float(coding_overlap_fraction),
                    context_cds_score=float(context_cds_score),
                    in_frame_repeat_periodicity=True,
                    tandem_hmm_support=True,
                    overlap_threshold=policy.coding_overlap_veto_fraction,
                    score_threshold=policy.context_cds_score_veto_minimum,
                )

    tandem_veto_assessed = False
    tandem_veto_triggered = False
    tandem_pair_supplied = (
        in_frame_repeat_periodicity is not None
        or tandem_hmm_support is not None
    )
    if tandem_pair_supplied and (
        in_frame_repeat_periodicity is None or tandem_hmm_support is None
    ):
        _reason(
            reasons,
            "incomplete_tandem_veto_inputs",
            REVIEW_REASON,
            "in-frame periodicity and tandem-HMM support must be supplied together",
        )
    elif tandem_pair_supplied:
        periodicity_valid = isinstance(in_frame_repeat_periodicity, bool)
        hmm_valid = isinstance(tandem_hmm_support, bool)
        if not periodicity_valid:
            _reason(
                reasons,
                "invalid_in_frame_repeat_periodicity",
                REJECT_REASON,
                "in_frame_repeat_periodicity must be an explicit boolean",
                observed=repr(in_frame_repeat_periodicity),
            )
        if not hmm_valid:
            _reason(
                reasons,
                "invalid_tandem_hmm_support",
                REJECT_REASON,
                "tandem_hmm_support must be an explicit boolean",
                observed=repr(tandem_hmm_support),
            )
        if periodicity_valid and hmm_valid:
            tandem_veto_assessed = True
            tandem_veto_triggered = (
                in_frame_repeat_periodicity and tandem_hmm_support
            )
            if tandem_veto_triggered:
                _reason(
                    reasons,
                    "coding_tandem_veto",
                    REJECT_REASON,
                    "in-frame periodicity and tandem-HMM support agree",
                )

    stats = ArchitectureRescueStats(
        reconstruction_valid=reconstruction_valid,
        reconstruction_issue_codes=reconstruction_issue_codes,
        array_start=array_start,
        array_end=array_end,
        repeat_count=len(repeats),
        spacer_count=len(spacers),
        nondeleted_spacer_count=len(normalized_spacers),
        deleted_spacer_count=deleted_spacer_count,
        consensus_length=consensus_length,
        nondeleted_spacer_lengths=spacer_lengths,
        mean_spacer_length=mean_spacer_length,
        spacer_length_cv=spacer_length_cv,
        unique_spacer_count=len(counts),
        exact_duplicate_pair_count=duplicate_pair_count,
        max_distinct_spacer_edit_identity=max_edit_identity,
        max_distinct_spacer_pair=max_edit_pair,
        dominant_prefix_5mer=prefix_kmer,
        prefix_5mer_dominance=prefix_dominance,
        dominant_suffix_5mer=suffix_kmer,
        suffix_5mer_dominance=suffix_dominance,
        classifier_score=score_value,
        repeat_similarity=repeat_similarity_value,
        spacer_similarity=spacer_similarity_value,
        afsf_pass=afsf_value,
        coding_veto_assessed=coding_veto_assessed,
        coding_veto_triggered=coding_veto_triggered,
        tandem_veto_assessed=tandem_veto_assessed,
        tandem_veto_triggered=tandem_veto_triggered,
    )

    if any(reason.severity == REJECT_REASON for reason in reasons):
        decision = INELIGIBLE
    elif any(reason.severity == REVIEW_REASON for reason in reasons):
        decision = REVIEW
    else:
        decision = ELIGIBLE

    return ArchitectureRescueDecision(
        decision=decision,
        reasons=tuple(reasons),
        stats=stats,
        reconstruction=reconstruction,
        policy=policy,
    )
