"""Conservative recovery of regular CRISPR subpaths from irregular parents.

This module is intentionally isolated from classification and category
mutation.  It treats every parent coordinate and sequence as an untrusted
claim, replays it against the supplied original source, and can only emit
source-derived candidate geometries.

The recovery scope is deliberately narrow: a reconstructed diagnostic parent
must fail the source-intact architecture gate *only* because of spacer length
and spacer-length variation, and the path must contain an extreme long-spacer
break.  The maximal source path on each side of every such break is evaluated
without searching smaller windows.  This avoids turning arbitrary repetitive
sequence into a subarray search problem.

Returned proposals are not classifications.  They carry ambiguous biological
boundary status and a source-bound payload for the existing coding-tandem
validator.  ``admission_ready`` is always false until a separate transaction
has run that veto and a fresh classifier evaluation.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
import statistics
from typing import Optional, Tuple

from components.reconstruction_validator import (
    IUPAC_DNA_ALPHABET,
    validate_reconstruction,
)
from components.source_intact_existence_rescue import (
    SOURCE_INTACT_EXISTENCE_POLICY,
    SourceIntactEvaluationContext,
    SourceIntactExistenceDecision,
    SourceIntactExistencePolicy,
    evaluate_source_intact_existence_in_context,
)


COORDINATE_SYSTEM = "zero_based_half_open"

DIAGNOSTIC_PARENT = "diagnostic"
ACCEPTED_PARENT = "accepted"

ABSTAINED = "abstained"
PROPOSAL_PENDING_TANDEM_VETO = "proposal_pending_tandem_veto"

BOUNDARY_AMBIGUOUS_PARENT_PHASE = "ambiguous_parent_phase"

_SPACER_ONLY_PARENT_REJECTIONS = frozenset({
    "spacer_length_out_of_range",
    "spacer_length_cv_above_maximum",
})


def _recovery_existence_policy():
    values = SOURCE_INTACT_EXISTENCE_POLICY.as_dict()
    values.update({
        "name": "source-regular-subpath-existence-gate",
        "version": "1.0.0",
        "minimum_repeat_count": 8,
    })
    return SourceIntactExistencePolicy(**values)


@dataclass(frozen=True)
class SourceRegularSubpathRecoveryPolicy:
    """Frozen thresholds for the conservative split proposal."""

    name: str = "source-regular-subpath-recovery"
    version: str = "1.0.0"
    existence_policy: SourceIntactExistencePolicy = field(
        default_factory=_recovery_existence_policy
    )
    minimum_break_length_over_spacer_maximum: int = 10
    minimum_break_to_reference_ratio: float = 2.5
    minimum_break_excess_over_reference: int = 40

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("policy name must be a non-empty string")
        if not isinstance(self.version, str) or not self.version:
            raise ValueError("policy version must be a non-empty string")
        if not isinstance(
            self.existence_policy, SourceIntactExistencePolicy
        ):
            raise ValueError(
                "existence_policy must be a SourceIntactExistencePolicy"
            )
        if self.existence_policy.minimum_repeat_count < 8:
            raise ValueError(
                "regular-subpath recovery requires at least eight repeats"
            )
        for field_name in (
            "minimum_break_length_over_spacer_maximum",
            "minimum_break_excess_over_reference",
        ):
            value = getattr(self, field_name)
            if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 1
            ):
                raise ValueError(
                    "{} must be a positive integer".format(field_name)
                )
        ratio = self.minimum_break_to_reference_ratio
        if (
            isinstance(ratio, bool)
            or not isinstance(ratio, (int, float))
            or not math.isfinite(float(ratio))
            or ratio <= 1.0
        ):
            raise ValueError(
                "minimum_break_to_reference_ratio must be finite and > 1"
            )

    def as_dict(self):
        return asdict(self)


SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY = (
    SourceRegularSubpathRecoveryPolicy()
)


@dataclass(frozen=True)
class SourceRegularSubpathContext:
    """Reusable source window bound to one hash and existence context."""

    source_sequence: str
    source_coordinate_offset: int = 0
    source_sha256: str = field(init=False)
    existence_context: SourceIntactEvaluationContext = field(
        init=False, repr=False
    )

    def __post_init__(self):
        if (
            not isinstance(self.source_coordinate_offset, int)
            or isinstance(self.source_coordinate_offset, bool)
            or self.source_coordinate_offset < 0
        ):
            raise ValueError(
                "source_coordinate_offset must be a non-negative integer"
            )
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
        existence_context = SourceIntactEvaluationContext(normalized)
        object.__setattr__(self, "source_sequence", normalized)
        object.__setattr__(
            self, "source_sha256", existence_context.source_sha256
        )
        object.__setattr__(
            self, "existence_context", existence_context
        )

    @property
    def source_sequence_interval(self):
        return (
            self.source_coordinate_offset,
            self.source_coordinate_offset + len(self.source_sequence),
        )


@dataclass(frozen=True)
class SourceRegularSubpathCandidate:
    """Immutable, source-derived candidate geometry in source coordinates."""

    list_repeat_starts: Tuple[int, ...]
    list_repeats: Tuple[str, ...]
    list_spacers: Tuple[str, ...]
    consensus: str

    @property
    def source_interval(self):
        return (
            self.list_repeat_starts[0],
            self.list_repeat_starts[-1] + len(self.list_repeats[-1]),
        )

    def as_dict(self):
        return {
            "repeat_starts": list(self.list_repeat_starts),
            "repeats": list(self.list_repeats),
            "spacers": list(self.list_spacers),
            "consensus": self.consensus,
            "source_interval": list(self.source_interval),
        }


@dataclass(frozen=True)
class SourceRegularSubpathBreak:
    """One extreme source spacer that separates maximal path segments."""

    spacer_index: int
    source_interval: Tuple[int, int]
    length: int
    reference_spacer_median: float
    length_to_reference_ratio: float
    excess_over_reference: float
    source_slice_sha256: str

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourceRegularSubpathTandemVetoPayload:
    """Source-bound geometry sufficient to run the coding-tandem veto."""

    required: bool
    source_coordinate_offset: int
    source_sequence_interval: Tuple[int, int]
    source_sha256: str
    candidate_interval: Tuple[int, int]
    candidate_source_slice_sha256: str
    candidate_geometry_sha256: str
    repeat_intervals: Tuple[Tuple[int, int], ...]
    spacer_intervals: Tuple[Tuple[int, int], ...]
    repeat_start_deltas: Tuple[int, ...]
    repeat_source_slice_sha256s: Tuple[str, ...]
    spacer_source_slice_sha256s: Tuple[str, ...]
    validator_contract: str

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourceRegularSubpathProposal:
    """One maximal regular segment, still pending external veto and scoring."""

    parent_identifier: str
    parent_interval: Tuple[int, int]
    parent_source_slice_sha256: str
    parent_index_span: Tuple[int, int]
    candidate: SourceRegularSubpathCandidate
    candidate_source_slice_sha256: str
    left_break: Optional[SourceRegularSubpathBreak]
    right_break: Optional[SourceRegularSubpathBreak]
    existence_decision: SourceIntactExistenceDecision
    boundary_status: str
    tandem_veto: SourceRegularSubpathTandemVetoPayload

    @property
    def source_interval(self):
        return self.candidate.source_interval

    @property
    def admission_ready(self):
        return False

    def as_dict(self):
        existence_record = self.existence_decision.as_dict()
        existence_metrics = dict(existence_record["metrics"])
        source_offset = self.tandem_veto.source_coordinate_offset
        for field_name in ("array_start", "array_end"):
            value = existence_metrics[field_name]
            if value is not None:
                existence_metrics[field_name] = value + source_offset
        existence_metrics["coordinate_system"] = COORDINATE_SYSTEM
        existence_metrics["source_sequence_interval"] = list(
            self.tandem_veto.source_sequence_interval
        )
        existence_record["metrics"] = existence_metrics
        return {
            "parent_identifier": self.parent_identifier,
            "parent_interval": list(self.parent_interval),
            "parent_source_slice_sha256": self.parent_source_slice_sha256,
            "parent_index_span": list(self.parent_index_span),
            "candidate": self.candidate.as_dict(),
            "candidate_source_slice_sha256": (
                self.candidate_source_slice_sha256
            ),
            "left_break": (
                None if self.left_break is None else self.left_break.as_dict()
            ),
            "right_break": (
                None
                if self.right_break is None
                else self.right_break.as_dict()
            ),
            "existence_decision": existence_record,
            "boundary_status": self.boundary_status,
            "tandem_veto": self.tandem_veto.as_dict(),
            "admission_ready": self.admission_ready,
            "classification_assigned": False,
        }


@dataclass(frozen=True)
class SourceRegularSubpathSegmentEvaluation:
    """Audit record for one maximal segment between strong breaks."""

    parent_index_span: Tuple[int, int]
    source_interval: Tuple[int, int]
    repeat_count: int
    eligible: bool
    rejection_codes: Tuple[str, ...]

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourceRegularSubpathSuppression:
    """A structurally valid proposal suppressed by accepted output."""

    proposal_interval: Tuple[int, int]
    overlapping_accepted_intervals: Tuple[Tuple[int, int], ...]

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourceRegularSubpathRecoveryResult:
    """Complete non-mutating decision for one parent candidate."""

    status: str
    reason_codes: Tuple[str, ...]
    parent_identifier: str
    parent_disposition: str
    source_coordinate_offset: int
    source_sequence_interval: Tuple[int, int]
    source_sha256: Optional[str]
    parent_interval: Optional[Tuple[int, int]]
    parent_source_slice_sha256: Optional[str]
    parent_reconstruction_issue_codes: Tuple[str, ...]
    parent_gate_rejection_codes: Tuple[str, ...]
    strong_breaks: Tuple[SourceRegularSubpathBreak, ...]
    segment_evaluations: Tuple[
        SourceRegularSubpathSegmentEvaluation, ...
    ]
    proposals: Tuple[SourceRegularSubpathProposal, ...]
    suppressions: Tuple[SourceRegularSubpathSuppression, ...]
    policy: SourceRegularSubpathRecoveryPolicy

    @property
    def admission_ready(self):
        return False

    def as_dict(self):
        return {
            "status": self.status,
            "reason_codes": list(self.reason_codes),
            "parent_identifier": self.parent_identifier,
            "parent_disposition": self.parent_disposition,
            "source_coordinate_offset": self.source_coordinate_offset,
            "source_sequence_interval": list(self.source_sequence_interval),
            "source_sha256": self.source_sha256,
            "parent_interval": (
                None
                if self.parent_interval is None
                else list(self.parent_interval)
            ),
            "parent_source_slice_sha256": self.parent_source_slice_sha256,
            "parent_reconstruction_issue_codes": list(
                self.parent_reconstruction_issue_codes
            ),
            "parent_gate_rejection_codes": list(
                self.parent_gate_rejection_codes
            ),
            "strong_breaks": [
                item.as_dict() for item in self.strong_breaks
            ],
            "segment_evaluations": [
                item.as_dict() for item in self.segment_evaluations
            ],
            "proposals": [item.as_dict() for item in self.proposals],
            "suppressions": [
                item.as_dict() for item in self.suppressions
            ],
            "policy": self.policy.as_dict(),
            "admission_ready": self.admission_ready,
            "classification_assigned": False,
        }


@dataclass(frozen=True)
class _SourcePath:
    starts: Tuple[int, ...]
    local_starts: Tuple[int, ...]
    repeats: Tuple[str, ...]
    spacers: Tuple[str, ...]
    source_interval: Tuple[int, int]
    source_slice_sha256: str


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _canonical_sha256(value):
    encoded = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _candidate_fields(candidate):
    try:
        starts = tuple(candidate.list_repeat_starts)
        repeats = tuple(candidate.list_repeats)
        spacers = tuple(candidate.list_spacers)
    except (AttributeError, TypeError):
        return None
    return starts, repeats, spacers


def _source_path(source, source_offset, candidate):
    fields = _candidate_fields(candidate)
    if fields is None:
        return None, ("malformed_candidate",)
    starts, claimed_repeats, claimed_spacers = fields
    local_starts = tuple(
        start - source_offset
        if isinstance(start, int) and not isinstance(start, bool)
        else start
        for start in starts
    )
    reconstruction = validate_reconstruction(
        source,
        local_starts,
        claimed_repeats,
        claimed_spacers,
    )
    if not reconstruction.valid:
        return None, tuple(issue.code for issue in reconstruction.issues)

    source_repeats = tuple(
        source[start:start + len(repeat)]
        for start, repeat in zip(local_starts, claimed_repeats)
    )
    source_spacers = tuple(
        source[start + len(repeat):next_start]
        for start, repeat, next_start in zip(
            local_starts,
            source_repeats,
            local_starts[1:],
        )
    )
    local_interval = (
        local_starts[0],
        local_starts[-1] + len(source_repeats[-1]),
    )
    source_interval = (
        local_interval[0] + source_offset,
        local_interval[1] + source_offset,
    )
    return _SourcePath(
        starts=tuple(starts),
        local_starts=local_starts,
        repeats=source_repeats,
        spacers=source_spacers,
        source_interval=source_interval,
        source_slice_sha256=_sha256(
            source[local_interval[0]:local_interval[1]]
        ),
    ), ()


def _candidate_from_path(path, start_index, end_index, consensus=None):
    repeats = path.repeats[start_index:end_index]
    if consensus is None:
        consensus = repeats[0]
    return SourceRegularSubpathCandidate(
        list_repeat_starts=path.starts[start_index:end_index],
        list_repeats=repeats,
        list_spacers=path.spacers[start_index:end_index - 1],
        consensus=consensus,
    )


def _local_candidate(candidate, source_offset):
    return SourceRegularSubpathCandidate(
        list_repeat_starts=tuple(
            start - source_offset
            for start in candidate.list_repeat_starts
        ),
        list_repeats=candidate.list_repeats,
        list_spacers=candidate.list_spacers,
        consensus=candidate.consensus,
    )


def _evaluate_candidate(context, candidate, source_offset, policy):
    return evaluate_source_intact_existence_in_context(
        context,
        _local_candidate(candidate, source_offset),
        policy=policy.existence_policy,
    )


def _strong_breaks(source, source_offset, path, policy):
    spacer_lengths = tuple(len(spacer) for spacer in path.spacers)
    lower = policy.existence_policy.minimum_spacer_length
    upper = policy.existence_policy.maximum_spacer_length
    reference_lengths = tuple(
        length for length in spacer_lengths if lower <= length <= upper
    )
    if not reference_lengths:
        return ()
    reference = float(statistics.median(reference_lengths))
    minimum_length = (
        upper + policy.minimum_break_length_over_spacer_maximum
    )
    records = []
    for index, length in enumerate(spacer_lengths):
        ratio = length / reference
        excess = length - reference
        if (
            length < minimum_length
            or ratio < policy.minimum_break_to_reference_ratio
            or excess < policy.minimum_break_excess_over_reference
        ):
            continue
        local_start = (
            path.local_starts[index] + len(path.repeats[index])
        )
        local_end = path.local_starts[index + 1]
        records.append(SourceRegularSubpathBreak(
            spacer_index=index,
            source_interval=(
                local_start + source_offset,
                local_end + source_offset,
            ),
            length=length,
            reference_spacer_median=reference,
            length_to_reference_ratio=ratio,
            excess_over_reference=excess,
            source_slice_sha256=_sha256(source[local_start:local_end]),
        ))
    return tuple(records)


def _segment_spans(repeat_count, breaks):
    boundaries = (0,) + tuple(
        item.spacer_index + 1 for item in breaks
    ) + (repeat_count,)
    return tuple(zip(boundaries, boundaries[1:]))


def _tandem_payload(
    source,
    source_sha256,
    source_offset,
    path,
    start_index,
    end_index,
    candidate,
):
    local_starts = path.local_starts[start_index:end_index]
    repeats = path.repeats[start_index:end_index]
    spacers = path.spacers[start_index:end_index - 1]
    repeat_intervals = tuple(
        (start + source_offset, start + source_offset + len(repeat))
        for start, repeat in zip(local_starts, repeats)
    )
    spacer_intervals = tuple(
        (
            local_starts[index] + len(repeats[index]) + source_offset,
            local_starts[index + 1] + source_offset,
        )
        for index in range(len(repeats) - 1)
    )
    candidate_interval = candidate.source_interval
    local_interval = (
        candidate_interval[0] - source_offset,
        candidate_interval[1] - source_offset,
    )
    geometry_record = {
        "coordinate_system": COORDINATE_SYSTEM,
        "repeat_starts": candidate.list_repeat_starts,
        "repeats": candidate.list_repeats,
        "spacers": candidate.list_spacers,
        "consensus": candidate.consensus,
    }
    return SourceRegularSubpathTandemVetoPayload(
        required=True,
        source_coordinate_offset=source_offset,
        source_sequence_interval=(
            source_offset,
            source_offset + len(source),
        ),
        source_sha256=source_sha256,
        candidate_interval=candidate_interval,
        candidate_source_slice_sha256=_sha256(
            source[local_interval[0]:local_interval[1]]
        ),
        candidate_geometry_sha256=_canonical_sha256(geometry_record),
        repeat_intervals=repeat_intervals,
        spacer_intervals=spacer_intervals,
        repeat_start_deltas=tuple(
            right - left for left, right in zip(
                candidate.list_repeat_starts,
                candidate.list_repeat_starts[1:],
            )
        ),
        repeat_source_slice_sha256s=tuple(
            _sha256(sequence) for sequence in repeats
        ),
        spacer_source_slice_sha256s=tuple(
            _sha256(sequence) for sequence in spacers
        ),
        validator_contract=(
            "accepted_family_tandem_validator."
            "validate_accepted_family_tandem"
        ),
    )


def _overlaps(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _result(
    *,
    status,
    reason_codes,
    parent_identifier,
    parent_disposition,
    source_offset,
    source_length,
    source_sha256,
    policy,
    parent_path=None,
    parent_issue_codes=(),
    parent_gate_codes=(),
    breaks=(),
    segment_evaluations=(),
    proposals=(),
    suppressions=(),
):
    return SourceRegularSubpathRecoveryResult(
        status=status,
        reason_codes=tuple(reason_codes),
        parent_identifier=parent_identifier,
        parent_disposition=parent_disposition,
        source_coordinate_offset=source_offset,
        source_sequence_interval=(
            source_offset,
            source_offset + source_length,
        ),
        source_sha256=source_sha256,
        parent_interval=(
            None if parent_path is None else parent_path.source_interval
        ),
        parent_source_slice_sha256=(
            None
            if parent_path is None
            else parent_path.source_slice_sha256
        ),
        parent_reconstruction_issue_codes=tuple(parent_issue_codes),
        parent_gate_rejection_codes=tuple(parent_gate_codes),
        strong_breaks=tuple(breaks),
        segment_evaluations=tuple(segment_evaluations),
        proposals=tuple(proposals),
        suppressions=tuple(suppressions),
        policy=policy,
    )


def propose_source_regular_subpaths(
    source_sequence,
    parent_candidate,
    *,
    parent_identifier,
    parent_disposition,
    accepted_candidates,
    source_coordinate_offset=0,
    policy=SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY,
):
    """Return source-replayed maximal regular subpaths, or abstain.

    ``source_sequence`` may be a full source record (the production case) or a
    pinned source window.  Candidate starts always remain in source-record
    coordinates; ``source_coordinate_offset`` gives the zero-based start of a
    supplied window.

    ``accepted_candidates`` is required rather than defaulted so an integration
    cannot accidentally omit the accepted-output protection.  Every accepted
    geometry is itself replayed against the same source before overlap is used.
    """

    if not isinstance(policy, SourceRegularSubpathRecoveryPolicy):
        raise TypeError(
            "policy must be a SourceRegularSubpathRecoveryPolicy"
        )
    if not isinstance(parent_identifier, str) or not parent_identifier:
        raise ValueError("parent_identifier must be a non-empty string")
    if parent_disposition not in (DIAGNOSTIC_PARENT, ACCEPTED_PARENT):
        raise ValueError(
            "parent_disposition must be diagnostic or accepted"
        )
    if not isinstance(accepted_candidates, tuple):
        raise TypeError("accepted_candidates must be a tuple")
    if (
        not isinstance(source_coordinate_offset, int)
        or isinstance(source_coordinate_offset, bool)
        or source_coordinate_offset < 0
    ):
        raise ValueError(
            "source_coordinate_offset must be a non-negative integer"
        )

    supplied_context = (
        source_sequence
        if isinstance(source_sequence, SourceRegularSubpathContext)
        else None
    )
    raw_source = (
        supplied_context.source_sequence
        if supplied_context is not None
        else source_sequence
    )
    if supplied_context is None and (
        not isinstance(raw_source, str)
        or not raw_source
        or set(raw_source).difference(IUPAC_DNA_ALPHABET)
    ):
        return _result(
            status=ABSTAINED,
            reason_codes=("invalid_source_sequence",),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=(
                len(raw_source)
                if isinstance(raw_source, str)
                else 0
            ),
            source_sha256=None,
            policy=policy,
        )

    if supplied_context is None:
        recovery_context = SourceRegularSubpathContext(
            raw_source,
            source_coordinate_offset,
        )
    else:
        recovery_context = supplied_context
        if (
            source_coordinate_offset
            != recovery_context.source_coordinate_offset
        ):
            raise ValueError(
                "context offset differs from source_coordinate_offset"
            )
    source = recovery_context.source_sequence
    source_sha = recovery_context.source_sha256
    context = recovery_context.existence_context
    parent_path, parent_issues = _source_path(
        source,
        source_coordinate_offset,
        parent_candidate,
    )
    if parent_path is None:
        return _result(
            status=ABSTAINED,
            reason_codes=("parent_source_reconstruction_failed",),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=len(source),
            source_sha256=source_sha,
            parent_issue_codes=parent_issues,
            policy=policy,
        )

    if parent_disposition == ACCEPTED_PARENT:
        return _result(
            status=ABSTAINED,
            reason_codes=("accepted_parent_protected",),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=len(source),
            source_sha256=source_sha,
            parent_path=parent_path,
            policy=policy,
        )

    accepted_paths = []
    for accepted_index, accepted_candidate in enumerate(accepted_candidates):
        accepted_path, accepted_issues = _source_path(
            source,
            source_coordinate_offset,
            accepted_candidate,
        )
        if accepted_path is None:
            return _result(
                status=ABSTAINED,
                reason_codes=(
                    "accepted_context_source_reconstruction_failed",
                    "accepted_context_index_{}".format(accepted_index),
                ) + tuple(
                    "accepted_context_issue_{}".format(code)
                    for code in accepted_issues
                ),
                parent_identifier=parent_identifier,
                parent_disposition=parent_disposition,
                source_offset=source_coordinate_offset,
                source_length=len(source),
                source_sha256=source_sha,
                parent_path=parent_path,
                policy=policy,
            )
        accepted_paths.append(accepted_path)

    parent_source_candidate = _candidate_from_path(
        parent_path,
        0,
        len(parent_path.repeats),
    )
    parent_decision = _evaluate_candidate(
        context,
        parent_source_candidate,
        source_coordinate_offset,
        policy,
    )
    parent_gate_codes = parent_decision.rejection_codes
    if parent_decision.eligible:
        return _result(
            status=ABSTAINED,
            reason_codes=("parent_is_already_regular",),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=len(source),
            source_sha256=source_sha,
            parent_path=parent_path,
            parent_gate_codes=parent_gate_codes,
            policy=policy,
        )
    if (
        not parent_gate_codes
        or not set(parent_gate_codes).issubset(
            _SPACER_ONLY_PARENT_REJECTIONS
        )
        or set(parent_gate_codes) != _SPACER_ONLY_PARENT_REJECTIONS
    ):
        return _result(
            status=ABSTAINED,
            reason_codes=("parent_has_non_spacer_architecture_failures",),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=len(source),
            source_sha256=source_sha,
            parent_path=parent_path,
            parent_gate_codes=parent_gate_codes,
            policy=policy,
        )

    breaks = _strong_breaks(
        source,
        source_coordinate_offset,
        parent_path,
        policy,
    )
    if not breaks:
        return _result(
            status=ABSTAINED,
            reason_codes=("extreme_source_spacer_break_not_found",),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=len(source),
            source_sha256=source_sha,
            parent_path=parent_path,
            parent_gate_codes=parent_gate_codes,
            policy=policy,
        )

    breaks_by_index = {item.spacer_index: item for item in breaks}
    evaluations = []
    structural_proposals = []
    for start_index, end_index in _segment_spans(
        len(parent_path.repeats), breaks
    ):
        repeat_count = end_index - start_index
        source_interval = (
            parent_path.starts[start_index],
            parent_path.starts[end_index - 1]
            + len(parent_path.repeats[end_index - 1]),
        )
        if repeat_count < policy.existence_policy.minimum_repeat_count:
            evaluations.append(SourceRegularSubpathSegmentEvaluation(
                parent_index_span=(start_index, end_index),
                source_interval=source_interval,
                repeat_count=repeat_count,
                eligible=False,
                rejection_codes=(
                    "repeat_count_below_recovery_minimum",
                ),
            ))
            continue

        preliminary = _candidate_from_path(
            parent_path, start_index, end_index
        )
        preliminary_decision = _evaluate_candidate(
            context,
            preliminary,
            source_coordinate_offset,
            policy,
        )
        if not preliminary_decision.eligible:
            evaluations.append(SourceRegularSubpathSegmentEvaluation(
                parent_index_span=(start_index, end_index),
                source_interval=source_interval,
                repeat_count=repeat_count,
                eligible=False,
                rejection_codes=preliminary_decision.rejection_codes,
            ))
            continue

        consensus = (
            preliminary_decision.metrics.source_repeat_medoid
        )
        candidate = _candidate_from_path(
            parent_path,
            start_index,
            end_index,
            consensus=consensus,
        )
        decision = _evaluate_candidate(
            context,
            candidate,
            source_coordinate_offset,
            policy,
        )
        if not decision.eligible:
            raise RuntimeError(
                "source-derived medoid changed an eligible subpath decision"
            )

        evaluations.append(SourceRegularSubpathSegmentEvaluation(
            parent_index_span=(start_index, end_index),
            source_interval=source_interval,
            repeat_count=repeat_count,
            eligible=True,
            rejection_codes=(),
        ))
        local_interval = (
            source_interval[0] - source_coordinate_offset,
            source_interval[1] - source_coordinate_offset,
        )
        structural_proposals.append(SourceRegularSubpathProposal(
            parent_identifier=parent_identifier,
            parent_interval=parent_path.source_interval,
            parent_source_slice_sha256=parent_path.source_slice_sha256,
            parent_index_span=(start_index, end_index),
            candidate=candidate,
            candidate_source_slice_sha256=_sha256(
                source[local_interval[0]:local_interval[1]]
            ),
            left_break=breaks_by_index.get(start_index - 1),
            right_break=breaks_by_index.get(end_index - 1),
            existence_decision=decision,
            boundary_status=BOUNDARY_AMBIGUOUS_PARENT_PHASE,
            tandem_veto=_tandem_payload(
                source,
                source_sha,
                source_coordinate_offset,
                parent_path,
                start_index,
                end_index,
                candidate,
            ),
        ))

    proposals = []
    suppressions = []
    accepted_intervals = tuple(
        item.source_interval for item in accepted_paths
    )
    for proposal in structural_proposals:
        overlaps = tuple(
            interval
            for interval in accepted_intervals
            if _overlaps(proposal.source_interval, interval)
        )
        if overlaps:
            suppressions.append(SourceRegularSubpathSuppression(
                proposal_interval=proposal.source_interval,
                overlapping_accepted_intervals=overlaps,
            ))
        else:
            proposals.append(proposal)

    if not proposals:
        reason = (
            "all_regular_subpaths_overlap_accepted_arrays"
            if suppressions
            else "no_eligible_maximal_regular_subpath"
        )
        return _result(
            status=ABSTAINED,
            reason_codes=(reason,),
            parent_identifier=parent_identifier,
            parent_disposition=parent_disposition,
            source_offset=source_coordinate_offset,
            source_length=len(source),
            source_sha256=source_sha,
            parent_path=parent_path,
            parent_gate_codes=parent_gate_codes,
            breaks=breaks,
            segment_evaluations=evaluations,
            suppressions=suppressions,
            policy=policy,
        )

    return _result(
        status=PROPOSAL_PENDING_TANDEM_VETO,
        reason_codes=(
            "source_regular_subpath_requires_tandem_veto_and_rescoring",
        ),
        parent_identifier=parent_identifier,
        parent_disposition=parent_disposition,
        source_offset=source_coordinate_offset,
        source_length=len(source),
        source_sha256=source_sha,
        parent_path=parent_path,
        parent_gate_codes=parent_gate_codes,
        breaks=breaks,
        segment_evaluations=evaluations,
        proposals=proposals,
        suppressions=suppressions,
        policy=policy,
    )


def propose_source_regular_subpaths_in_context(
    context,
    parent_candidate,
    *,
    parent_identifier,
    parent_disposition,
    accepted_candidates,
    policy=SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY,
):
    """Evaluate one parent while hashing and validating its source once."""

    if not isinstance(context, SourceRegularSubpathContext):
        raise TypeError("context must be a SourceRegularSubpathContext")
    return propose_source_regular_subpaths(
        context,
        parent_candidate,
        parent_identifier=parent_identifier,
        parent_disposition=parent_disposition,
        accepted_candidates=accepted_candidates,
        source_coordinate_offset=context.source_coordinate_offset,
        policy=policy,
    )


__all__ = (
    "ABSTAINED",
    "ACCEPTED_PARENT",
    "BOUNDARY_AMBIGUOUS_PARENT_PHASE",
    "COORDINATE_SYSTEM",
    "DIAGNOSTIC_PARENT",
    "PROPOSAL_PENDING_TANDEM_VETO",
    "SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY",
    "SourceRegularSubpathBreak",
    "SourceRegularSubpathCandidate",
    "SourceRegularSubpathContext",
    "SourceRegularSubpathProposal",
    "SourceRegularSubpathRecoveryPolicy",
    "SourceRegularSubpathRecoveryResult",
    "SourceRegularSubpathSegmentEvaluation",
    "SourceRegularSubpathSuppression",
    "SourceRegularSubpathTandemVetoPayload",
    "propose_source_regular_subpaths",
    "propose_source_regular_subpaths_in_context",
)
