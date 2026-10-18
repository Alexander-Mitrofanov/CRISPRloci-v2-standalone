"""Pure source-intact evidence for existence of a CRISPR array.

This module deliberately answers only whether one already-detected candidate
has a conservative, classifier-independent CRISPR architecture.  It does not
choose a repeat phase, certify exact boundaries, classify a candidate, mutate
categories, or promote output.

Candidate coordinates and sequences are treated as claims.  Reconstruction is
checked against the complete original source, and every sequence used for a
metric is sliced back from that source.  An eligible result therefore carries
``AMBIGUOUS`` boundary status: existence evidence is not boundary evidence.
"""

from collections import Counter
from dataclasses import asdict, dataclass, field
import hashlib
import math
from numbers import Real
from typing import Optional, Tuple

from components.reconstruction_validator import (
    IUPAC_DNA_ALPHABET,
    validate_reconstruction,
)


ELIGIBLE = "eligible"
INELIGIBLE = "ineligible"

BOUNDARY_UNASSESSED = "UNASSESSED"
BOUNDARY_AMBIGUOUS = "AMBIGUOUS"


@dataclass(frozen=True)
class SourceIntactExistencePolicy:
    """Versioned thresholds for the conservative source-intact gate."""

    name: str = "source-intact-existence-rescue"
    version: str = "1.0.0"
    minimum_repeat_count: int = 5
    minimum_consensus_length: int = 20
    maximum_consensus_length: int = 55
    minimum_repeat_length: int = 20
    maximum_repeat_length: int = 55
    maximum_repeat_length_cv: float = 0.08
    minimum_spacer_length: int = 18
    maximum_spacer_length: int = 80
    maximum_spacer_length_cv: float = 0.08
    minimum_mean_repeat_identity: float = 0.90
    minimum_repeat_identity: float = 0.80
    duplicate_pair_minimum_spacer_count: int = 7
    duplicate_pair_minimum_other_distinct_spacers: int = 5
    maximum_distinct_spacer_identity: float = 0.65
    terminal_kmer_length: int = 5
    maximum_prefix_kmer_dominance: float = 0.50
    maximum_suffix_kmer_dominance: float = 0.50

    def __post_init__(self):
        for field_name in ("name", "version"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value:
                raise ValueError(
                    "{} must be a non-empty string".format(field_name)
                )

        integer_minima = {
            "minimum_repeat_count": 2,
            "minimum_consensus_length": 1,
            "maximum_consensus_length": 1,
            "minimum_repeat_length": 1,
            "maximum_repeat_length": 1,
            "minimum_spacer_length": 1,
            "maximum_spacer_length": 1,
            "duplicate_pair_minimum_spacer_count": 2,
            "duplicate_pair_minimum_other_distinct_spacers": 1,
            "terminal_kmer_length": 1,
        }
        for field_name, minimum in integer_minima.items():
            value = getattr(self, field_name)
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < minimum
            ):
                raise ValueError(
                    "{} must be an integer >= {}".format(
                        field_name, minimum
                    )
                )

        for field_name in (
            "maximum_repeat_length_cv",
            "maximum_spacer_length_cv",
            "minimum_mean_repeat_identity",
            "minimum_repeat_identity",
            "maximum_distinct_spacer_identity",
            "maximum_prefix_kmer_dominance",
            "maximum_suffix_kmer_dominance",
        ):
            value = getattr(self, field_name)
            if (
                not isinstance(value, Real)
                or isinstance(value, bool)
                or not math.isfinite(float(value))
                or not 0.0 <= float(value) <= 1.0
            ):
                raise ValueError(
                    "{} must be a finite number in [0, 1]".format(
                        field_name
                    )
                )

        if self.minimum_consensus_length > self.maximum_consensus_length:
            raise ValueError("consensus-length bounds are inconsistent")
        if self.minimum_repeat_length > self.maximum_repeat_length:
            raise ValueError("repeat-length bounds are inconsistent")
        if self.minimum_spacer_length > self.maximum_spacer_length:
            raise ValueError("spacer-length bounds are inconsistent")
        if self.terminal_kmer_length > self.minimum_spacer_length:
            raise ValueError(
                "terminal_kmer_length cannot exceed the minimum spacer length"
            )
        if (
            self.duplicate_pair_minimum_spacer_count
            < self.duplicate_pair_minimum_other_distinct_spacers + 2
        ):
            raise ValueError(
                "the duplicate-pair spacer minimum must allow the pair and "
                "all required other distinct spacers"
            )

    def as_dict(self):
        return asdict(self)


SOURCE_INTACT_EXISTENCE_POLICY = SourceIntactExistencePolicy()


@dataclass(frozen=True)
class SourceIntactEvaluationContext:
    """One normalized, source-bound context reusable across candidates.

    Construction validates and uppercases the complete IUPAC source, then
    computes its SHA-256 exactly once.  The digest is not caller-supplied, so
    a context cannot accidentally bind candidates to a stale source hash.
    """

    source_sequence: str
    source_sha256: str = field(init=False)

    def __post_init__(self):
        source = self.source_sequence
        if (
            not isinstance(source, str)
            or not source
            or set(source).difference(IUPAC_DNA_ALPHABET)
        ):
            raise ValueError(
                "source_sequence must be one non-empty IUPAC DNA string"
            )
        normalized = source.upper()
        object.__setattr__(self, "source_sequence", normalized)
        object.__setattr__(self, "source_sha256", _sha256(normalized))


@dataclass(frozen=True)
class SourceIntactExistenceReason:
    """One stable rejection reason."""

    code: str
    message: str

    def as_dict(self):
        return {"code": self.code, "message": self.message}


@dataclass(frozen=True)
class SourceIntactExistenceMetrics:
    """Immutable source-derived measurements, including unavailable values."""

    source_length: Optional[int]
    source_sha256: Optional[str]
    reconstruction_valid: bool
    reconstruction_issue_codes: Tuple[str, ...]
    array_start: Optional[int]
    array_end: Optional[int]
    array_source_slice_sha256: Optional[str]
    repeat_count: Optional[int]
    supplied_spacer_count: Optional[int]
    nonempty_supplied_spacer_count: Optional[int]
    consensus_length: Optional[int]
    source_repeat_sequences: Tuple[str, ...]
    source_spacer_sequences: Tuple[str, ...]
    source_repeat_lengths: Tuple[int, ...]
    source_spacer_lengths: Tuple[int, ...]
    mean_repeat_length: Optional[float]
    repeat_length_cv: Optional[float]
    mean_spacer_length: Optional[float]
    spacer_length_cv: Optional[float]
    source_repeat_medoid: Optional[str]
    source_repeat_medoid_total_edit_distance: Optional[int]
    mean_repeat_identity_to_medoid: Optional[float]
    minimum_repeat_identity_to_medoid: Optional[float]
    distinct_spacer_count: Optional[int]
    duplicate_spacer_multiplicities: Tuple[int, ...]
    exact_duplicate_spacer_pair_count: Optional[int]
    maximum_distinct_spacer_identity: Optional[float]
    maximum_distinct_spacer_pair: Optional[Tuple[str, str]]
    dominant_prefix_kmer: Optional[str]
    prefix_kmer_dominance: Optional[float]
    dominant_suffix_kmer: Optional[str]
    suffix_kmer_dominance: Optional[float]

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourceIntactExistenceDecision:
    """A non-mutating existence decision with no exact-boundary state."""

    decision: str
    reasons: Tuple[SourceIntactExistenceReason, ...]
    metrics: SourceIntactExistenceMetrics
    boundary_status: str
    policy: SourceIntactExistencePolicy

    @property
    def eligible(self):
        return self.decision == ELIGIBLE

    @property
    def rejection_codes(self):
        return tuple(reason.code for reason in self.reasons)

    def as_dict(self):
        return {
            "decision": self.decision,
            "eligible": self.eligible,
            "rejection_codes": list(self.rejection_codes),
            "reasons": [reason.as_dict() for reason in self.reasons],
            "metrics": self.metrics.as_dict(),
            "boundary_status": self.boundary_status,
            "policy": self.policy.as_dict(),
        }


_REASON_MESSAGES = {
    "invalid_source_sequence": (
        "source must be one non-empty string containing only IUPAC DNA symbols"
    ),
    "malformed_candidate": (
        "candidate must expose iterable repeat starts, repeats, and spacers"
    ),
    "source_reconstruction_failed": (
        "candidate coordinates or sequences do not reconstruct the source"
    ),
    "source_slicing_unavailable": (
        "safe repeat and spacer intervals could not be sliced from the source"
    ),
    "repeat_count_below_minimum": (
        "candidate has fewer source repeats than the existence minimum"
    ),
    "spacer_count_or_nonempty_rule_failed": (
        "candidate must supply exactly one non-empty spacer between repeats"
    ),
    "invalid_consensus": (
        "candidate consensus must be one non-empty ungapped IUPAC DNA string"
    ),
    "consensus_length_out_of_range": (
        "candidate consensus length is outside the calibrated interval"
    ),
    "repeat_length_out_of_range": (
        "one or more source-sliced repeats are outside the calibrated interval"
    ),
    "repeat_length_cv_above_maximum": (
        "source repeat-length population CV exceeds the calibrated maximum"
    ),
    "spacer_length_out_of_range": (
        "one or more source-sliced spacers are outside the calibrated interval"
    ),
    "spacer_length_cv_above_maximum": (
        "source spacer-length population CV exceeds the calibrated maximum"
    ),
    "mean_repeat_identity_below_minimum": (
        "mean source-repeat identity to the deterministic medoid is too low"
    ),
    "minimum_repeat_identity_below_minimum": (
        "minimum source-repeat identity to the deterministic medoid is too low"
    ),
    "exact_spacer_duplicate_rule_failed": (
        "exact source-spacer duplication exceeds the one-pair exception"
    ),
    "distinct_spacer_identity_above_maximum": (
        "two distinct source spacers are too similar"
    ),
    "prefix_kmer_dominance_above_maximum": (
        "source-spacer prefix k-mer dominance exceeds the calibrated maximum"
    ),
    "suffix_kmer_dominance_above_maximum": (
        "source-spacer suffix k-mer dominance exceeds the calibrated maximum"
    ),
}


def _reason(code):
    return SourceIntactExistenceReason(
        code=code,
        message=_REASON_MESSAGES[code],
    )


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _coefficient_of_variation(values):
    if not values:
        return None, None
    mean = sum(values) / float(len(values))
    if mean == 0.0:
        return mean, 0.0
    variance = sum((value - mean) ** 2 for value in values) / len(values)
    return mean, math.sqrt(variance) / mean


def _edit_distance(first, second):
    """Return deterministic Levenshtein distance without optional packages."""

    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for first_index, first_character in enumerate(first, 1):
        current = [first_index]
        for second_index, second_character in enumerate(second, 1):
            current.append(min(
                current[-1] + 1,
                previous[second_index] + 1,
                previous[second_index - 1]
                + (first_character != second_character),
            ))
        previous = current
    return previous[-1]


def _edit_identity(first, second):
    denominator = max(len(first), len(second))
    if denominator == 0:
        return 1.0
    return 1.0 - (_edit_distance(first, second) / float(denominator))


def _deterministic_medoid(sequences):
    unique = tuple(sorted(set(sequences)))
    ranked = tuple(
        (
            sum(_edit_distance(candidate, sequence) for sequence in sequences),
            candidate,
        )
        for candidate in unique
    )
    total_distance, medoid = min(ranked)
    return medoid, total_distance


def _maximum_distinct_identity(sequences):
    distinct = tuple(sorted(set(sequences)))
    records = tuple(
        (_edit_identity(left, right), (left, right))
        for left_index, left in enumerate(distinct)
        for right in distinct[left_index + 1 :]
    )
    if not records:
        return None, None
    maximum = max(record[0] for record in records)
    pair = min(record[1] for record in records if record[0] == maximum)
    return maximum, pair


def _dominant_terminal_kmer(sequences, length, from_end=False):
    if not sequences:
        return None, None
    kmers = tuple(
        sequence[-length:] if from_end else sequence[:length]
        for sequence in sequences
    )
    counts = Counter(kmers)
    maximum = max(counts.values())
    dominant = min(kmer for kmer, count in counts.items() if count == maximum)
    return dominant, maximum / float(len(kmers))


def _as_tuple_field(candidate, name):
    value = getattr(candidate, name)
    if isinstance(value, (str, bytes)) or value is None:
        raise TypeError("{} must be an iterable, not a string".format(name))
    return tuple(value)


def _safe_source_slices(source, starts, repeats):
    """Slice repeat/spacer sequences using only safe source intervals."""

    if (
        not repeats
        or len(starts) != len(repeats)
        or any(
            not isinstance(start, int) or isinstance(start, bool)
            for start in starts
        )
        or any(not isinstance(repeat, str) or not repeat for repeat in repeats)
    ):
        return None
    ends = tuple(
        start + len(repeat) for start, repeat in zip(starts, repeats)
    )
    if any(
        start < 0 or end <= start or end > len(source)
        for start, end in zip(starts, ends)
    ):
        return None
    if any(
        right_start <= left_start or right_start < left_end
        for left_start, left_end, right_start in zip(
            starts, ends, starts[1:]
        )
    ):
        return None
    source_repeats = tuple(
        source[start:end] for start, end in zip(starts, ends)
    )
    source_spacers = tuple(
        source[end:next_start] for end, next_start in zip(ends, starts[1:])
    )
    return (
        source_repeats,
        source_spacers,
        starts[0],
        ends[-1],
    )


def _empty_metrics(
    source=None,
    *,
    source_sha256=None,
    repeat_count=None,
    supplied_spacer_count=None,
    nonempty_supplied_spacer_count=None,
    consensus_length=None,
    reconstruction_issue_codes=(),
):
    source_length = len(source) if isinstance(source, str) else None
    if source_sha256 is None and isinstance(source, str):
        source_sha256 = _sha256(source)
    return SourceIntactExistenceMetrics(
        source_length=source_length,
        source_sha256=source_sha256,
        reconstruction_valid=False,
        reconstruction_issue_codes=tuple(reconstruction_issue_codes),
        array_start=None,
        array_end=None,
        array_source_slice_sha256=None,
        repeat_count=repeat_count,
        supplied_spacer_count=supplied_spacer_count,
        nonempty_supplied_spacer_count=nonempty_supplied_spacer_count,
        consensus_length=consensus_length,
        source_repeat_sequences=(),
        source_spacer_sequences=(),
        source_repeat_lengths=(),
        source_spacer_lengths=(),
        mean_repeat_length=None,
        repeat_length_cv=None,
        mean_spacer_length=None,
        spacer_length_cv=None,
        source_repeat_medoid=None,
        source_repeat_medoid_total_edit_distance=None,
        mean_repeat_identity_to_medoid=None,
        minimum_repeat_identity_to_medoid=None,
        distinct_spacer_count=None,
        duplicate_spacer_multiplicities=(),
        exact_duplicate_spacer_pair_count=None,
        maximum_distinct_spacer_identity=None,
        maximum_distinct_spacer_pair=None,
        dominant_prefix_kmer=None,
        prefix_kmer_dominance=None,
        dominant_suffix_kmer=None,
        suffix_kmer_dominance=None,
    )


def _closed_decision(reasons, metrics, policy):
    reason_tuple = tuple(reasons)
    eligible = not reason_tuple
    return SourceIntactExistenceDecision(
        decision=ELIGIBLE if eligible else INELIGIBLE,
        reasons=reason_tuple,
        metrics=metrics,
        boundary_status=(
            BOUNDARY_AMBIGUOUS if eligible else BOUNDARY_UNASSESSED
        ),
        policy=policy,
    )


def evaluate_source_intact_existence_in_context(
    context,
    candidate,
    *,
    policy=SOURCE_INTACT_EXISTENCE_POLICY,
):
    """Evaluate one candidate using a prevalidated, prehashed source."""

    if not isinstance(policy, SourceIntactExistencePolicy):
        raise TypeError(
            "policy must be a SourceIntactExistencePolicy instance"
        )

    if not isinstance(context, SourceIntactEvaluationContext):
        raise TypeError("context must be a SourceIntactEvaluationContext")

    reasons = []
    source = context.source_sequence

    try:
        starts = _as_tuple_field(candidate, "list_repeat_starts")
        supplied_repeats = _as_tuple_field(candidate, "list_repeats")
        supplied_spacers = _as_tuple_field(candidate, "list_spacers")
    except (AttributeError, TypeError):
        reasons.append(_reason("malformed_candidate"))
        return _closed_decision(
            reasons,
            _empty_metrics(
                source,
                source_sha256=context.source_sha256,
            ),
            policy,
        )

    repeat_count = len(supplied_repeats)
    supplied_spacer_count = len(supplied_spacers)
    nonempty_supplied_spacer_count = sum(
        isinstance(spacer, str) and bool(spacer)
        for spacer in supplied_spacers
    )
    consensus = getattr(candidate, "consensus", None)
    consensus_valid = (
        isinstance(consensus, str)
        and bool(consensus)
        and not set(consensus).difference(IUPAC_DNA_ALPHABET)
    )
    consensus_length = len(consensus) if consensus_valid else None

    reconstruction = validate_reconstruction(
        source,
        starts,
        supplied_repeats,
        supplied_spacers,
    )
    if not reconstruction.valid:
        reasons.append(_reason("source_reconstruction_failed"))

    if repeat_count < policy.minimum_repeat_count:
        reasons.append(_reason("repeat_count_below_minimum"))

    if (
        supplied_spacer_count != max(0, repeat_count - 1)
        or nonempty_supplied_spacer_count != supplied_spacer_count
    ):
        reasons.append(_reason("spacer_count_or_nonempty_rule_failed"))

    if not consensus_valid:
        reasons.append(_reason("invalid_consensus"))
    elif not (
        policy.minimum_consensus_length
        <= consensus_length
        <= policy.maximum_consensus_length
    ):
        reasons.append(_reason("consensus_length_out_of_range"))

    sliced = _safe_source_slices(source, starts, supplied_repeats)
    if sliced is None:
        reasons.append(_reason("source_slicing_unavailable"))
        metrics = _empty_metrics(
            source,
            source_sha256=context.source_sha256,
            repeat_count=repeat_count,
            supplied_spacer_count=supplied_spacer_count,
            nonempty_supplied_spacer_count=(
                nonempty_supplied_spacer_count
            ),
            consensus_length=consensus_length,
            reconstruction_issue_codes=tuple(
                issue.code for issue in reconstruction.issues
            ),
        )
        return _closed_decision(reasons, metrics, policy)

    source_repeats, source_spacers, array_start, array_end = sliced
    repeat_lengths = tuple(len(sequence) for sequence in source_repeats)
    spacer_lengths = tuple(len(sequence) for sequence in source_spacers)
    mean_repeat_length, repeat_length_cv = _coefficient_of_variation(
        repeat_lengths
    )
    mean_spacer_length, spacer_length_cv = _coefficient_of_variation(
        spacer_lengths
    )

    medoid, medoid_total_distance = _deterministic_medoid(source_repeats)
    repeat_identities = tuple(
        _edit_identity(sequence, medoid) for sequence in source_repeats
    )
    mean_repeat_identity = sum(repeat_identities) / len(repeat_identities)
    minimum_repeat_identity = min(repeat_identities)

    spacer_counts = Counter(source_spacers)
    duplicate_multiplicities = tuple(sorted(
        count for count in spacer_counts.values() if count > 1
    ))
    exact_duplicate_pair_count = sum(
        count * (count - 1) // 2 for count in spacer_counts.values()
    )
    distinct_spacer_count = len(spacer_counts)
    maximum_distinct_identity, maximum_distinct_pair = (
        _maximum_distinct_identity(source_spacers)
    )
    prefix_kmer, prefix_dominance = _dominant_terminal_kmer(
        source_spacers,
        policy.terminal_kmer_length,
    )
    suffix_kmer, suffix_dominance = _dominant_terminal_kmer(
        source_spacers,
        policy.terminal_kmer_length,
        from_end=True,
    )

    if any(
        length < policy.minimum_repeat_length
        or length > policy.maximum_repeat_length
        for length in repeat_lengths
    ):
        reasons.append(_reason("repeat_length_out_of_range"))
    if repeat_length_cv > policy.maximum_repeat_length_cv:
        reasons.append(_reason("repeat_length_cv_above_maximum"))

    if (
        not spacer_lengths
        or any(
            length < policy.minimum_spacer_length
            or length > policy.maximum_spacer_length
            for length in spacer_lengths
        )
    ):
        reasons.append(_reason("spacer_length_out_of_range"))
    if (
        spacer_length_cv is None
        or spacer_length_cv > policy.maximum_spacer_length_cv
    ):
        reasons.append(_reason("spacer_length_cv_above_maximum"))

    if (
        mean_repeat_identity
        < policy.minimum_mean_repeat_identity
    ):
        reasons.append(_reason("mean_repeat_identity_below_minimum"))
    if minimum_repeat_identity < policy.minimum_repeat_identity:
        reasons.append(_reason("minimum_repeat_identity_below_minimum"))

    duplicate_rule_passes = not duplicate_multiplicities
    if duplicate_multiplicities == (2,):
        duplicate_rule_passes = (
            len(source_spacers)
            >= policy.duplicate_pair_minimum_spacer_count
            and distinct_spacer_count - 1
            >= policy.duplicate_pair_minimum_other_distinct_spacers
        )
    if not duplicate_rule_passes:
        reasons.append(_reason("exact_spacer_duplicate_rule_failed"))

    if (
        maximum_distinct_identity is not None
        and maximum_distinct_identity
        > policy.maximum_distinct_spacer_identity
    ):
        reasons.append(_reason("distinct_spacer_identity_above_maximum"))
    if (
        prefix_dominance is None
        or prefix_dominance
        > policy.maximum_prefix_kmer_dominance
    ):
        reasons.append(_reason("prefix_kmer_dominance_above_maximum"))
    if (
        suffix_dominance is None
        or suffix_dominance
        > policy.maximum_suffix_kmer_dominance
    ):
        reasons.append(_reason("suffix_kmer_dominance_above_maximum"))

    metrics = SourceIntactExistenceMetrics(
        source_length=len(source),
        source_sha256=context.source_sha256,
        reconstruction_valid=reconstruction.valid,
        reconstruction_issue_codes=tuple(
            issue.code for issue in reconstruction.issues
        ),
        array_start=array_start,
        array_end=array_end,
        array_source_slice_sha256=_sha256(source[array_start:array_end]),
        repeat_count=repeat_count,
        supplied_spacer_count=supplied_spacer_count,
        nonempty_supplied_spacer_count=nonempty_supplied_spacer_count,
        consensus_length=consensus_length,
        source_repeat_sequences=source_repeats,
        source_spacer_sequences=source_spacers,
        source_repeat_lengths=repeat_lengths,
        source_spacer_lengths=spacer_lengths,
        mean_repeat_length=mean_repeat_length,
        repeat_length_cv=repeat_length_cv,
        mean_spacer_length=mean_spacer_length,
        spacer_length_cv=spacer_length_cv,
        source_repeat_medoid=medoid,
        source_repeat_medoid_total_edit_distance=medoid_total_distance,
        mean_repeat_identity_to_medoid=mean_repeat_identity,
        minimum_repeat_identity_to_medoid=minimum_repeat_identity,
        distinct_spacer_count=distinct_spacer_count,
        duplicate_spacer_multiplicities=duplicate_multiplicities,
        exact_duplicate_spacer_pair_count=exact_duplicate_pair_count,
        maximum_distinct_spacer_identity=maximum_distinct_identity,
        maximum_distinct_spacer_pair=maximum_distinct_pair,
        dominant_prefix_kmer=prefix_kmer,
        prefix_kmer_dominance=prefix_dominance,
        dominant_suffix_kmer=suffix_kmer,
        suffix_kmer_dominance=suffix_dominance,
    )
    return _closed_decision(reasons, metrics, policy)


def evaluate_source_intact_existence(
    source_sequence,
    candidate,
    *,
    policy=SOURCE_INTACT_EXISTENCE_POLICY,
):
    """Return a pure, exhaustive source-intact existence decision.

    This compatibility wrapper preserves the historical fail-closed behavior
    for malformed biological source input.  Callers evaluating many
    candidates from one genome should construct one
    :class:`SourceIntactEvaluationContext` and call
    :func:`evaluate_source_intact_existence_in_context` instead.
    """

    if not isinstance(policy, SourceIntactExistencePolicy):
        raise TypeError(
            "policy must be a SourceIntactExistencePolicy instance"
        )
    if (
        not isinstance(source_sequence, str)
        or not source_sequence
        or set(source_sequence).difference(IUPAC_DNA_ALPHABET)
    ):
        return _closed_decision(
            (_reason("invalid_source_sequence"),),
            _empty_metrics(),
            policy,
        )
    context = SourceIntactEvaluationContext(source_sequence)
    return evaluate_source_intact_existence_in_context(
        context,
        candidate,
        policy=policy,
    )
