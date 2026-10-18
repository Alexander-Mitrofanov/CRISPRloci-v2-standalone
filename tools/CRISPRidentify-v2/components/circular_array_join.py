"""Conservative, source-replayed joins for arrays crossing a circular origin.

The ordinary :class:`CrisprCandidate` representation is linear: repeat starts
must increase and one array has one contiguous interval.  Consequently, this
module does not manufacture a ``start > end`` candidate.  It constructs an
immutable, orientation-aware linearized *view* for scoring while retaining the
true source geometry as ordered, zero-based half-open segments.

This module is intentionally not wired into the pipeline.  A join is accepted
only when independently supplied topology and compatibility evidence is
complete and source-bound, both children replay exactly against the complete
source, the modulo-origin spacer is structurally plausible, explicit
structural and tandem gates pass, and a scorer evaluates the newly constructed
view at or above the caller's Possible threshold.  Scores attached to either
child are neither accepted nor inspected.
"""

from dataclasses import dataclass
from datetime import date
import hashlib
import json
import math
from numbers import Real
import re
from typing import Callable, Optional, Tuple

from components.reconstruction_validator import (
    IUPAC_DNA_ALPHABET,
    validate_reconstruction,
)
from components.sequence_utils import reverse_complement


RULE_NAME = "source-replayed-circular-boundary-join"
RULE_VERSION = "1.0.0"
COORDINATE_SYSTEM = "zero_based_half_open"

TOPOLOGY_CIRCULAR = "circular"
TOPOLOGY_LINEAR = "linear"
TOPOLOGY_UNKNOWN = "unknown"
TOPOLOGIES = frozenset((
    TOPOLOGY_CIRCULAR,
    TOPOLOGY_LINEAR,
    TOPOLOGY_UNKNOWN,
))

ORIENTATION_FORWARD = "forward"
ORIENTATION_REVERSE = "reverse"
ORIENTATIONS = frozenset((ORIENTATION_FORWARD, ORIENTATION_REVERSE))

STATUS_JOINED = "joined"
STATUS_ABSTAINED = "abstained"
STATUS_REJECTED = "rejected"

GATE_STRUCTURAL = "structural"
GATE_TANDEM = "tandem"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _canonical_sha256(value):
    encoded = json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    return _sha256(encoded)


def _require_text(value, field_name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("{} must be a non-empty string".format(field_name))


def _require_sha256(value, field_name):
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "{} must be a lowercase hexadecimal SHA-256".format(field_name)
        )


def _require_bool(value, field_name):
    if not isinstance(value, bool):
        raise ValueError("{} must be a boolean".format(field_name))


def _normalize_source(source):
    if not isinstance(source, str) or not source:
        raise ValueError("source must be one non-empty DNA string")
    try:
        source.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("source must contain ASCII DNA symbols") from error
    if set(source).difference(IUPAC_DNA_ALPHABET):
        raise ValueError("source contains non-IUPAC DNA symbols")
    return source.upper()


def source_sha256(source):
    """Return the hash used to bind evidence to a normalized source."""

    return _sha256(_normalize_source(source))


@dataclass(frozen=True)
class CircularTopologyEvidence:
    """Immutable external evidence for the source molecule's topology."""

    topology: str
    accession: str
    source_sha256: str
    evidence_method: str
    evidence_date: str
    evidence_sha256: str
    evidence_complete: bool
    evidence_ambiguous: bool

    def __post_init__(self):
        if self.topology not in TOPOLOGIES:
            raise ValueError("topology must be circular, linear, or unknown")
        _require_text(self.accession, "accession")
        _require_sha256(self.source_sha256, "source_sha256")
        _require_text(self.evidence_method, "evidence_method")
        _require_text(self.evidence_date, "evidence_date")
        try:
            date.fromisoformat(self.evidence_date)
        except ValueError as error:
            raise ValueError("evidence_date must be an ISO calendar date") from error
        _require_sha256(self.evidence_sha256, "evidence_sha256")
        _require_bool(self.evidence_complete, "evidence_complete")
        _require_bool(self.evidence_ambiguous, "evidence_ambiguous")

    def as_dict(self):
        return {
            "topology": self.topology,
            "accession": self.accession,
            "source_sha256": self.source_sha256,
            "evidence_method": self.evidence_method,
            "evidence_date": self.evidence_date,
            "evidence_sha256": self.evidence_sha256,
            "evidence_complete": self.evidence_complete,
            "evidence_ambiguous": self.evidence_ambiguous,
        }


@dataclass(frozen=True)
class CircularJoinCompatibilityEvidence:
    """Source-bound assessment that the two edge children are one family."""

    child_geometry_sha256s: Tuple[str, str]
    repeat_family_compatible: bool
    repeat_phase_compatible: bool
    orientation_compatible: bool
    orientation: str
    assessment_complete: bool
    assessment_ambiguous: bool
    evidence_method: str
    evidence_sha256: str

    def __post_init__(self):
        if not isinstance(self.child_geometry_sha256s, tuple):
            raise ValueError("child_geometry_sha256s must be an immutable tuple")
        if len(self.child_geometry_sha256s) != 2:
            raise ValueError("exactly two child geometry hashes are required")
        for index, value in enumerate(self.child_geometry_sha256s):
            _require_sha256(
                value,
                "child_geometry_sha256s[{}]".format(index),
            )
        if self.child_geometry_sha256s != tuple(
                sorted(self.child_geometry_sha256s)):
            raise ValueError("child geometry hashes must be sorted")
        for field_name in (
            "repeat_family_compatible",
            "repeat_phase_compatible",
            "orientation_compatible",
            "assessment_complete",
            "assessment_ambiguous",
        ):
            _require_bool(getattr(self, field_name), field_name)
        if self.orientation not in ORIENTATIONS:
            raise ValueError("orientation must be forward or reverse")
        _require_text(self.evidence_method, "evidence_method")
        _require_sha256(self.evidence_sha256, "evidence_sha256")

    def as_dict(self):
        return {
            "child_geometry_sha256s": list(self.child_geometry_sha256s),
            "repeat_family_compatible": self.repeat_family_compatible,
            "repeat_phase_compatible": self.repeat_phase_compatible,
            "orientation_compatible": self.orientation_compatible,
            "orientation": self.orientation,
            "assessment_complete": self.assessment_complete,
            "assessment_ambiguous": self.assessment_ambiguous,
            "evidence_method": self.evidence_method,
            "evidence_sha256": self.evidence_sha256,
        }


@dataclass(frozen=True)
class CircularArrayJoinPolicy:
    """Explicit structural spacer bounds and Possible score threshold."""

    minimum_junction_spacer_length: int
    maximum_junction_spacer_length: int
    possible_score_threshold: float

    def __post_init__(self):
        minimum = self.minimum_junction_spacer_length
        maximum = self.maximum_junction_spacer_length
        if (
            not isinstance(minimum, int)
            or isinstance(minimum, bool)
            or minimum < 0
        ):
            raise ValueError(
                "minimum_junction_spacer_length must be a non-negative integer"
            )
        if (
            not isinstance(maximum, int)
            or isinstance(maximum, bool)
            or maximum < minimum
        ):
            raise ValueError(
                "maximum_junction_spacer_length must be an integer >= minimum"
            )
        threshold = self.possible_score_threshold
        if (
            not isinstance(threshold, Real)
            or isinstance(threshold, bool)
            or not math.isfinite(float(threshold))
            or not 0.0 <= float(threshold) <= 1.0
        ):
            raise ValueError("possible_score_threshold must be finite in [0, 1]")

    def as_dict(self):
        return {
            "minimum_junction_spacer_length": (
                self.minimum_junction_spacer_length
            ),
            "maximum_junction_spacer_length": (
                self.maximum_junction_spacer_length
            ),
            "possible_score_threshold": float(self.possible_score_threshold),
        }


@dataclass(frozen=True)
class CircularSourceSegment:
    """One legal source interval, with source and oriented sequence hashes."""

    role: str
    interval: Tuple[int, int]
    traversal: str
    source_slice_sha256: str
    oriented_slice_sha256: str

    @property
    def length(self):
        return self.interval[1] - self.interval[0]

    def as_dict(self):
        return {
            "role": self.role,
            "interval": list(self.interval),
            "traversal": self.traversal,
            "length": self.length,
            "source_slice_sha256": self.source_slice_sha256,
            "oriented_slice_sha256": self.oriented_slice_sha256,
        }


@dataclass(frozen=True)
class CircularChildGeometry:
    """Runtime-object-free geometry independently replayed from one child."""

    geometry_sha256: str
    source_sha256: str
    source_interval: Tuple[int, int]
    source_slice_sha256: str
    repeat_intervals: Tuple[Tuple[int, int], ...]
    repeat_source_slice_sha256s: Tuple[str, ...]
    spacer_intervals: Tuple[Tuple[int, int], ...]
    spacer_source_slice_sha256s: Tuple[str, ...]
    repeat_count: int
    spacer_count: int

    def as_dict(self):
        return {
            "geometry_sha256": self.geometry_sha256,
            "source_sha256": self.source_sha256,
            "coordinate_system": COORDINATE_SYSTEM,
            "source_interval": list(self.source_interval),
            "source_slice_sha256": self.source_slice_sha256,
            "repeat_intervals": [
                list(interval) for interval in self.repeat_intervals
            ],
            "repeat_source_slice_sha256s": list(
                self.repeat_source_slice_sha256s
            ),
            "spacer_intervals": [
                list(interval) for interval in self.spacer_intervals
            ],
            "spacer_source_slice_sha256s": list(
                self.spacer_source_slice_sha256s
            ),
            "repeat_count": self.repeat_count,
            "spacer_count": self.spacer_count,
        }


@dataclass(frozen=True)
class SegmentedJunctionSpacer:
    """The origin-crossing spacer as two non-wrapping source segments."""

    orientation: str
    sequence: str
    sequence_sha256: str
    physical_source_sequence_sha256: str
    source_segments: Tuple[CircularSourceSegment, CircularSourceSegment]

    @property
    def length(self):
        return len(self.sequence)

    def as_dict(self):
        return {
            "orientation": self.orientation,
            "sequence": self.sequence,
            "length": self.length,
            "sequence_sha256": self.sequence_sha256,
            "physical_source_sequence_sha256": (
                self.physical_source_sequence_sha256
            ),
            "source_segments": [
                segment.as_dict() for segment in self.source_segments
            ],
        }


@dataclass(frozen=True)
class CircularLinearizedView:
    """Fresh candidate-like sequence view in the supplied array orientation."""

    orientation: str
    repeats: Tuple[str, ...]
    spacers: Tuple[str, ...]
    linear_repeat_starts: Tuple[int, ...]
    repeat_source_intervals: Tuple[Tuple[int, int], ...]
    repeat_source_slice_sha256s: Tuple[str, ...]
    repeat_oriented_sequence_sha256s: Tuple[str, ...]
    spacer_source_segments: Tuple[Tuple[Tuple[int, int], ...], ...]
    spacer_source_segment_sha256s: Tuple[Tuple[str, ...], ...]
    spacer_oriented_sequence_sha256s: Tuple[str, ...]
    child_geometry_sha256s_in_order: Tuple[str, str]
    source_sha256: str
    content_sha256: str
    view_sha256: str

    @property
    def list_repeats(self):
        """Candidate-like alias for callbacks that consume repeat sequences."""

        return self.repeats

    @property
    def list_spacers(self):
        """Candidate-like alias for callbacks that consume spacer sequences."""

        return self.spacers

    @property
    def list_repeat_starts(self):
        """Candidate-like alias in the artificial linear coordinate space."""

        return self.linear_repeat_starts

    def as_dict(self):
        return {
            "orientation": self.orientation,
            "repeats": list(self.repeats),
            "spacers": list(self.spacers),
            "linear_repeat_starts": list(self.linear_repeat_starts),
            "repeat_source_intervals": [
                list(interval) for interval in self.repeat_source_intervals
            ],
            "repeat_source_slice_sha256s": list(
                self.repeat_source_slice_sha256s
            ),
            "repeat_oriented_sequence_sha256s": list(
                self.repeat_oriented_sequence_sha256s
            ),
            "spacer_source_segments": [
                [list(interval) for interval in segments]
                for segments in self.spacer_source_segments
            ],
            "spacer_source_segment_sha256s": [
                list(hashes) for hashes in self.spacer_source_segment_sha256s
            ],
            "spacer_oriented_sequence_sha256s": list(
                self.spacer_oriented_sequence_sha256s
            ),
            "child_geometry_sha256s_in_order": list(
                self.child_geometry_sha256s_in_order
            ),
            "source_sha256": self.source_sha256,
            "content_sha256": self.content_sha256,
            "view_sha256": self.view_sha256,
        }


@dataclass(frozen=True)
class CircularJoinGateDecision:
    """One explicit source-view-bound structural or tandem gate decision."""

    gate: str
    passed: bool
    reason_code: str
    subject_view_sha256: str
    evidence_method: str
    evidence_sha256: str

    def __post_init__(self):
        if self.gate not in (GATE_STRUCTURAL, GATE_TANDEM):
            raise ValueError("gate must be structural or tandem")
        _require_bool(self.passed, "passed")
        _require_text(self.reason_code, "reason_code")
        _require_sha256(self.subject_view_sha256, "subject_view_sha256")
        _require_text(self.evidence_method, "evidence_method")
        _require_sha256(self.evidence_sha256, "evidence_sha256")

    def as_dict(self):
        return {
            "gate": self.gate,
            "passed": self.passed,
            "reason_code": self.reason_code,
            "subject_view_sha256": self.subject_view_sha256,
            "evidence_method": self.evidence_method,
            "evidence_sha256": self.evidence_sha256,
        }


@dataclass(frozen=True)
class CircularJoinScoreEvidence:
    """Score computed from the new view, never inherited from a child."""

    score: float
    subject_view_sha256: str
    evidence_method: str
    evidence_sha256: str

    def __post_init__(self):
        if (
            not isinstance(self.score, Real)
            or isinstance(self.score, bool)
            or not math.isfinite(float(self.score))
            or not 0.0 <= float(self.score) <= 1.0
        ):
            raise ValueError("score must be finite in [0, 1]")
        _require_sha256(self.subject_view_sha256, "subject_view_sha256")
        _require_text(self.evidence_method, "evidence_method")
        _require_sha256(self.evidence_sha256, "evidence_sha256")

    def as_dict(self):
        return {
            "score": float(self.score),
            "subject_view_sha256": self.subject_view_sha256,
            "evidence_method": self.evidence_method,
            "evidence_sha256": self.evidence_sha256,
        }


@dataclass(frozen=True)
class CircularJoinCheck:
    """One deterministic check evaluated while considering a join."""

    check: str
    passed: bool
    reason_code: str

    def as_dict(self):
        return {
            "check": self.check,
            "passed": self.passed,
            "reason_code": self.reason_code,
        }


@dataclass(frozen=True)
class CircularJoinProposal:
    """Segmented source geometry and fresh linearized view for one pair."""

    orientation: str
    ordered_source_segments: Tuple[
        CircularSourceSegment,
        CircularSourceSegment,
    ]
    locus_source_segments: Tuple[
        CircularSourceSegment,
        CircularSourceSegment,
    ]
    junction_spacer: SegmentedJunctionSpacer
    linearized_view: CircularLinearizedView
    proposal_sha256: str

    def as_dict(self):
        return {
            "orientation": self.orientation,
            "ordered_source_segments": [
                segment.as_dict() for segment in self.ordered_source_segments
            ],
            "locus_source_segments": [
                segment.as_dict() for segment in self.locus_source_segments
            ],
            "junction_spacer": self.junction_spacer.as_dict(),
            "linearized_view": self.linearized_view.as_dict(),
            "proposal_sha256": self.proposal_sha256,
        }


@dataclass(frozen=True)
class CircularArrayJoinDecision:
    """Complete, JSON-safe outcome of one conservative join attempt."""

    status: str
    reason_codes: Tuple[str, ...]
    source_sha256: str
    topology_evidence: CircularTopologyEvidence
    compatibility_evidence: CircularJoinCompatibilityEvidence
    policy: CircularArrayJoinPolicy
    checks: Tuple[CircularJoinCheck, ...]
    child_geometries: Tuple[CircularChildGeometry, ...]
    proposal: Optional[CircularJoinProposal]
    structural_decision: Optional[CircularJoinGateDecision]
    tandem_decision: Optional[CircularJoinGateDecision]
    score_evidence: Optional[CircularJoinScoreEvidence]

    @property
    def joined(self):
        return self.status == STATUS_JOINED

    def as_dict(self):
        record = {
            "rule": {
                "name": RULE_NAME,
                "version": RULE_VERSION,
                "coordinate_system": COORDINATE_SYSTEM,
                "automatic_join_requires_explicit_circular_topology": True,
                "child_scores_inherited": False,
                "junction_representation": "segmented_non_wrapping_intervals",
            },
            "status": self.status,
            "joined": self.joined,
            "reason_codes": list(self.reason_codes),
            "source_sha256": self.source_sha256,
            "topology_evidence": self.topology_evidence.as_dict(),
            "compatibility_evidence": self.compatibility_evidence.as_dict(),
            "policy": self.policy.as_dict(),
            "checks": [check.as_dict() for check in self.checks],
            "child_geometries": [
                geometry.as_dict() for geometry in self.child_geometries
            ],
            "proposal": None if self.proposal is None else self.proposal.as_dict(),
            "structural_decision": (
                None
                if self.structural_decision is None
                else self.structural_decision.as_dict()
            ),
            "tandem_decision": (
                None
                if self.tandem_decision is None
                else self.tandem_decision.as_dict()
            ),
            "score_evidence": (
                None
                if self.score_evidence is None
                else self.score_evidence.as_dict()
            ),
        }
        # This is an invariant check as well as a stable deep copy.  It makes
        # accidental runtime-object or NaN leakage fail at the module boundary.
        return json.loads(json.dumps(
            record,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ))


class CircularJoinValidationError(ValueError):
    """A source/candidate claim could not be replayed exactly."""

    def __init__(self, reason_codes):
        self.reason_codes = tuple(sorted(set(reason_codes)))
        super().__init__(", ".join(self.reason_codes))


@dataclass(frozen=True)
class _ChildData:
    geometry: CircularChildGeometry
    repeats: Tuple[str, ...]
    spacers: Tuple[str, ...]


def _source_interval(source, interval, claimed_sequence=None):
    start, end = interval
    if (
        not isinstance(start, int)
        or isinstance(start, bool)
        or not isinstance(end, int)
        or isinstance(end, bool)
        or not 0 <= start <= end <= len(source)
    ):
        raise CircularJoinValidationError(("source_interval_out_of_bounds",))
    observed = source[start:end]
    if claimed_sequence is not None and observed != claimed_sequence.upper():
        raise CircularJoinValidationError(("source_interval_sequence_mismatch",))
    return observed


def _snapshot_child_data(source, candidate, source_hash=None):
    try:
        starts = tuple(candidate.list_repeat_starts)
        repeats = tuple(candidate.list_repeats)
        spacers = tuple(candidate.list_spacers)
    except (AttributeError, TypeError) as error:
        raise CircularJoinValidationError(("invalid_child_candidate_fields",)) from error

    try:
        reconstruction = validate_reconstruction(
            source,
            starts,
            repeats,
            spacers,
        )
    except (TypeError, ValueError) as error:
        raise CircularJoinValidationError(("invalid_child_reconstruction",)) from error
    if not reconstruction.valid:
        raise CircularJoinValidationError(tuple(
            "child_{}".format(issue.code) for issue in reconstruction.issues
        ))

    normalized_repeats = tuple(repeat.upper() for repeat in repeats)
    normalized_spacers = tuple(spacer.upper() for spacer in spacers)
    genomic_starts = tuple(reconstruction.genomic_repeat_starts)
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(genomic_starts, normalized_repeats)
    )
    spacer_intervals = tuple(
        (repeat_intervals[index][1], repeat_intervals[index + 1][0])
        for index in range(len(repeat_intervals) - 1)
    )

    try:
        repeat_slices = tuple(
            _source_interval(source, interval, repeat)
            for interval, repeat in zip(repeat_intervals, normalized_repeats)
        )
        spacer_slices = tuple(
            _source_interval(source, interval, spacer)
            for interval, spacer in zip(spacer_intervals, normalized_spacers)
        )
        source_interval = (repeat_intervals[0][0], repeat_intervals[-1][1])
        source_slice = _source_interval(source, source_interval)
    except IndexError as error:
        raise CircularJoinValidationError(("child_has_no_repeat_intervals",)) from error

    if source_hash is None:
        source_hash = _sha256(source)
    payload = {
        "coordinate_system": COORDINATE_SYSTEM,
        "source_sha256": source_hash,
        "source_interval": source_interval,
        "source_slice_sha256": _sha256(source_slice),
        "repeat_intervals": repeat_intervals,
        "repeat_source_slice_sha256s": tuple(
            _sha256(sequence) for sequence in repeat_slices
        ),
        "spacer_intervals": spacer_intervals,
        "spacer_source_slice_sha256s": tuple(
            _sha256(sequence) for sequence in spacer_slices
        ),
    }
    geometry = CircularChildGeometry(
        geometry_sha256=_canonical_sha256(payload),
        source_sha256=source_hash,
        source_interval=source_interval,
        source_slice_sha256=_sha256(source_slice),
        repeat_intervals=repeat_intervals,
        repeat_source_slice_sha256s=tuple(
            _sha256(sequence) for sequence in repeat_slices
        ),
        spacer_intervals=spacer_intervals,
        spacer_source_slice_sha256s=tuple(
            _sha256(sequence) for sequence in spacer_slices
        ),
        repeat_count=len(repeat_intervals),
        spacer_count=len(spacer_intervals),
    )
    return _ChildData(
        geometry=geometry,
        repeats=repeat_slices,
        spacers=spacer_slices,
    )


def snapshot_circular_child(source, candidate):
    """Replay a child and return its source-bound immutable geometry."""

    normalized_source = _normalize_source(source)
    source_hash = _sha256(normalized_source)
    return _snapshot_child_data(
        normalized_source, candidate, source_hash=source_hash
    ).geometry


def snapshot_circular_children(source, candidates):
    """Replay children while hashing their immutable source exactly once."""

    normalized_source = _normalize_source(source)
    source_hash = _sha256(normalized_source)
    return tuple(
        _snapshot_child_data(
            normalized_source, candidate, source_hash=source_hash
        ).geometry
        for candidate in candidates
    )


def _segment(source, role, interval, orientation):
    source_slice = _source_interval(source, interval)
    if orientation == ORIENTATION_FORWARD:
        traversal = "increasing"
        oriented = source_slice
    else:
        traversal = "decreasing"
        oriented = reverse_complement(source_slice)
    return CircularSourceSegment(
        role=role,
        interval=interval,
        traversal=traversal,
        source_slice_sha256=_sha256(source_slice),
        oriented_slice_sha256=_sha256(oriented),
    )


def _orient_child(child, orientation):
    if orientation == ORIENTATION_FORWARD:
        repeats = child.repeats
        spacers = child.spacers
        repeat_intervals = child.geometry.repeat_intervals
        spacer_segments = tuple(
            (interval,) for interval in child.geometry.spacer_intervals
        )
    else:
        repeats = tuple(
            reverse_complement(sequence)
            for sequence in reversed(child.repeats)
        )
        spacers = tuple(
            reverse_complement(sequence)
            for sequence in reversed(child.spacers)
        )
        repeat_intervals = tuple(reversed(child.geometry.repeat_intervals))
        spacer_segments = tuple(
            (interval,)
            for interval in reversed(child.geometry.spacer_intervals)
        )
    return repeats, spacers, repeat_intervals, spacer_segments


def _linear_repeat_starts(repeats, spacers):
    cursor = 0
    starts = []
    for index, repeat in enumerate(repeats):
        starts.append(cursor)
        cursor += len(repeat)
        if index < len(spacers):
            cursor += len(spacers[index])
    return tuple(starts)


def _build_proposal(source, head, tail, orientation):
    tail_end = tail.geometry.source_interval[1]
    head_start = head.geometry.source_interval[0]
    tail_junction_interval = (tail_end, len(source))
    head_junction_interval = (0, head_start)
    tail_junction = _source_interval(source, tail_junction_interval)
    head_junction = _source_interval(source, head_junction_interval)
    physical_junction = tail_junction + head_junction

    if orientation == ORIENTATION_FORWARD:
        first_child = tail
        second_child = head
        junction_sequence = physical_junction
        junction_segments = (
            _segment(
                source,
                "junction_tail_segment",
                tail_junction_interval,
                orientation,
            ),
            _segment(
                source,
                "junction_head_segment",
                head_junction_interval,
                orientation,
            ),
        )
        ordered_source_segments = (
            _segment(
                source,
                "tail_child",
                tail.geometry.source_interval,
                orientation,
            ),
            _segment(
                source,
                "head_child",
                head.geometry.source_interval,
                orientation,
            ),
        )
        locus_source_segments = (
            _segment(
                source,
                "tail_locus_through_origin",
                (tail.geometry.source_interval[0], len(source)),
                orientation,
            ),
            _segment(
                source,
                "head_locus_after_origin",
                (0, head.geometry.source_interval[1]),
                orientation,
            ),
        )
    else:
        first_child = head
        second_child = tail
        junction_sequence = reverse_complement(physical_junction)
        junction_segments = (
            _segment(
                source,
                "junction_head_segment",
                head_junction_interval,
                orientation,
            ),
            _segment(
                source,
                "junction_tail_segment",
                tail_junction_interval,
                orientation,
            ),
        )
        ordered_source_segments = (
            _segment(
                source,
                "head_child",
                head.geometry.source_interval,
                orientation,
            ),
            _segment(
                source,
                "tail_child",
                tail.geometry.source_interval,
                orientation,
            ),
        )
        locus_source_segments = (
            _segment(
                source,
                "head_locus_before_origin",
                (0, head.geometry.source_interval[1]),
                orientation,
            ),
            _segment(
                source,
                "tail_locus_through_origin",
                (tail.geometry.source_interval[0], len(source)),
                orientation,
            ),
        )

    (
        first_repeats,
        first_spacers,
        first_repeat_intervals,
        first_spacer_segments,
    ) = _orient_child(first_child, orientation)
    (
        second_repeats,
        second_spacers,
        second_repeat_intervals,
        second_spacer_segments,
    ) = _orient_child(second_child, orientation)

    repeats = first_repeats + second_repeats
    spacers = first_spacers + (junction_sequence,) + second_spacers
    repeat_intervals = first_repeat_intervals + second_repeat_intervals
    spacer_segments = (
        first_spacer_segments
        + (tuple(segment.interval for segment in junction_segments),)
        + second_spacer_segments
    )
    starts = _linear_repeat_starts(repeats, spacers)
    content_payload = {
        "orientation": orientation,
        "repeats": repeats,
        "spacers": spacers,
    }
    content_hash = _canonical_sha256(content_payload)
    view_payload = {
        "content_sha256": content_hash,
        "source_sha256": _sha256(source),
        "repeat_source_intervals": repeat_intervals,
        "spacer_source_segments": spacer_segments,
        "child_geometry_sha256s_in_order": (
            first_child.geometry.geometry_sha256,
            second_child.geometry.geometry_sha256,
        ),
    }
    view = CircularLinearizedView(
        orientation=orientation,
        repeats=repeats,
        spacers=spacers,
        linear_repeat_starts=starts,
        repeat_source_intervals=repeat_intervals,
        repeat_source_slice_sha256s=tuple(
            _sha256(_source_interval(source, interval))
            for interval in repeat_intervals
        ),
        repeat_oriented_sequence_sha256s=tuple(
            _sha256(item) for item in repeats
        ),
        spacer_source_segments=spacer_segments,
        spacer_source_segment_sha256s=tuple(
            tuple(
                _sha256(_source_interval(source, interval))
                for interval in segments
            )
            for segments in spacer_segments
        ),
        spacer_oriented_sequence_sha256s=tuple(
            _sha256(item) for item in spacers
        ),
        child_geometry_sha256s_in_order=(
            first_child.geometry.geometry_sha256,
            second_child.geometry.geometry_sha256,
        ),
        source_sha256=_sha256(source),
        content_sha256=content_hash,
        view_sha256=_canonical_sha256(view_payload),
    )
    junction = SegmentedJunctionSpacer(
        orientation=orientation,
        sequence=junction_sequence,
        sequence_sha256=_sha256(junction_sequence),
        physical_source_sequence_sha256=_sha256(physical_junction),
        source_segments=junction_segments,
    )
    proposal_payload = {
        "orientation": orientation,
        "ordered_source_segments": tuple(
            (
                segment.role,
                segment.interval,
                segment.traversal,
                segment.source_slice_sha256,
                segment.oriented_slice_sha256,
            )
            for segment in ordered_source_segments
        ),
        "locus_source_segments": tuple(
            (
                segment.role,
                segment.interval,
                segment.traversal,
                segment.source_slice_sha256,
                segment.oriented_slice_sha256,
            )
            for segment in locus_source_segments
        ),
        "junction_spacer": junction.as_dict(),
        "view_sha256": view.view_sha256,
    }
    return CircularJoinProposal(
        orientation=orientation,
        ordered_source_segments=ordered_source_segments,
        locus_source_segments=locus_source_segments,
        junction_spacer=junction,
        linearized_view=view,
        proposal_sha256=_canonical_sha256(proposal_payload),
    )


def _check(check, passed, reason_code):
    return CircularJoinCheck(
        check=check,
        passed=bool(passed),
        reason_code=reason_code,
    )


def _reason_codes(checks):
    return tuple(sorted({
        check.reason_code for check in checks if not check.passed
    }))


def _result(
    status,
    source_hash,
    topology_evidence,
    compatibility_evidence,
    policy,
    checks,
    child_geometries=(),
    proposal=None,
    structural_decision=None,
    tandem_decision=None,
    score_evidence=None,
    success_reason=None,
):
    reasons = _reason_codes(checks)
    if not reasons and success_reason is not None:
        reasons = (success_reason,)
    return CircularArrayJoinDecision(
        status=status,
        reason_codes=reasons,
        source_sha256=source_hash,
        topology_evidence=topology_evidence,
        compatibility_evidence=compatibility_evidence,
        policy=policy,
        checks=tuple(checks),
        child_geometries=tuple(child_geometries),
        proposal=proposal,
        structural_decision=structural_decision,
        tandem_decision=tandem_decision,
        score_evidence=score_evidence,
    )


def _resolve_gate(expected_gate, view, callback, explicit_decision):
    if (callback is None) == (explicit_decision is None):
        return None, "{}_gate_evidence_configuration_invalid".format(expected_gate)
    if callback is not None:
        try:
            decision = callback(view)
        except Exception:  # A gate failure must fail closed, regardless of origin.
            return None, "{}_gate_callback_failed".format(expected_gate)
    else:
        decision = explicit_decision
    if not isinstance(decision, CircularJoinGateDecision):
        return None, "{}_gate_decision_invalid".format(expected_gate)
    if decision.gate != expected_gate:
        return None, "{}_gate_kind_mismatch".format(expected_gate)
    if decision.subject_view_sha256 != view.view_sha256:
        return None, "{}_gate_stale_view".format(expected_gate)
    return decision, None


def evaluate_circular_array_join(
    source,
    first_candidate,
    second_candidate,
    topology_evidence,
    compatibility_evidence,
    policy,
    *,
    structural_callback: Optional[
        Callable[[CircularLinearizedView], CircularJoinGateDecision]
    ] = None,
    tandem_callback: Optional[
        Callable[[CircularLinearizedView], CircularJoinGateDecision]
    ] = None,
    scoring_callback: Optional[
        Callable[[CircularLinearizedView], CircularJoinScoreEvidence]
    ] = None,
    structural_decision: Optional[CircularJoinGateDecision] = None,
    tandem_decision: Optional[CircularJoinGateDecision] = None,
):
    """Evaluate a two-child circular-origin join without mutating either child.

    ``first_candidate`` and ``second_candidate`` are deliberately unordered;
    their source-replayed intervals determine head and tail roles.  Structural
    and tandem decisions can be supplied directly when already bound to the
    exact view hash, or computed by callbacks.  A scoring callback is always
    required so a child score cannot be inherited accidentally.
    """

    if not isinstance(topology_evidence, CircularTopologyEvidence):
        raise TypeError("topology_evidence has the wrong type")
    if not isinstance(
            compatibility_evidence, CircularJoinCompatibilityEvidence):
        raise TypeError("compatibility_evidence has the wrong type")
    if not isinstance(policy, CircularArrayJoinPolicy):
        raise TypeError("policy has the wrong type")

    normalized_source = _normalize_source(source)
    source_hash = _sha256(normalized_source)
    checks = []

    topology_source_matches = topology_evidence.source_sha256 == source_hash
    checks.append(_check(
        "topology_source_binding",
        topology_source_matches,
        "topology_source_hash_mismatch",
    ))

    child_data = []
    child_failure_codes = []
    for candidate in (first_candidate, second_candidate):
        try:
            child_data.append(_snapshot_child_data(
                normalized_source, candidate, source_hash=source_hash
            ))
        except CircularJoinValidationError as error:
            child_failure_codes.extend(error.reason_codes)
    children_valid = len(child_data) == 2 and not child_failure_codes
    if child_failure_codes:
        for reason_code in sorted(set(child_failure_codes)):
            checks.append(_check(
                "child_source_replay",
                False,
                reason_code,
            ))
    else:
        checks.append(_check(
            "child_source_replay",
            children_valid,
            "child_source_replay_failed",
        ))
    child_geometries = tuple(sorted(
        (item.geometry for item in child_data),
        key=lambda item: (item.source_interval, item.geometry_sha256),
    ))
    if not topology_source_matches or not children_valid:
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
        )

    child_data.sort(key=lambda item: (
        item.geometry.source_interval,
        item.geometry.geometry_sha256,
    ))
    head, tail = child_data
    observed_hashes = tuple(sorted((
        head.geometry.geometry_sha256,
        tail.geometry.geometry_sha256,
    )))
    distinct = observed_hashes[0] != observed_hashes[1]
    checks.append(_check(
        "distinct_child_geometry",
        distinct,
        "duplicate_child_geometry",
    ))
    disjoint = head.geometry.source_interval[1] <= tail.geometry.source_interval[0]
    checks.append(_check(
        "disjoint_linear_children",
        disjoint,
        "child_source_intervals_overlap",
    ))

    topology_complete = topology_evidence.evidence_complete
    checks.append(_check(
        "topology_evidence_complete",
        topology_complete,
        "topology_evidence_incomplete",
    ))
    topology_unambiguous = not topology_evidence.evidence_ambiguous
    checks.append(_check(
        "topology_evidence_unambiguous",
        topology_unambiguous,
        "topology_evidence_ambiguous",
    ))
    is_circular = topology_evidence.topology == TOPOLOGY_CIRCULAR
    topology_reason = {
        TOPOLOGY_CIRCULAR: "explicit_circular_topology",
        TOPOLOGY_LINEAR: "explicit_linear_topology",
        TOPOLOGY_UNKNOWN: "unknown_topology",
    }[topology_evidence.topology]
    checks.append(_check("circular_topology", is_circular, topology_reason))

    compatibility_bound = (
        compatibility_evidence.child_geometry_sha256s == observed_hashes
    )
    checks.append(_check(
        "compatibility_source_binding",
        compatibility_bound,
        "compatibility_child_hash_mismatch",
    ))
    checks.append(_check(
        "compatibility_assessment_complete",
        compatibility_evidence.assessment_complete,
        "compatibility_assessment_incomplete",
    ))
    checks.append(_check(
        "compatibility_assessment_unambiguous",
        not compatibility_evidence.assessment_ambiguous,
        "compatibility_assessment_ambiguous",
    ))
    checks.append(_check(
        "repeat_family_compatibility",
        compatibility_evidence.repeat_family_compatible,
        "repeat_family_mismatch",
    ))
    checks.append(_check(
        "repeat_phase_compatibility",
        compatibility_evidence.repeat_phase_compatible,
        "repeat_phase_mismatch",
    ))
    checks.append(_check(
        "orientation_compatibility",
        compatibility_evidence.orientation_compatible,
        "orientation_mismatch",
    ))

    junction_length = (
        len(normalized_source)
        - tail.geometry.source_interval[1]
        + head.geometry.source_interval[0]
    )
    above_minimum = junction_length >= policy.minimum_junction_spacer_length
    below_maximum = junction_length <= policy.maximum_junction_spacer_length
    checks.append(_check(
        "junction_spacer_minimum",
        above_minimum,
        "junction_spacer_below_minimum",
    ))
    checks.append(_check(
        "junction_spacer_maximum",
        below_maximum,
        "junction_spacer_above_maximum",
    ))

    if not distinct or not disjoint:
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
        )

    evidence_checks = checks[2:]
    if any(not check.passed for check in evidence_checks):
        return _result(
            STATUS_ABSTAINED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
        )

    proposal = _build_proposal(
        normalized_source,
        head,
        tail,
        compatibility_evidence.orientation,
    )
    view = proposal.linearized_view

    structural, structural_error = _resolve_gate(
        GATE_STRUCTURAL,
        view,
        structural_callback,
        structural_decision,
    )
    if structural_error is not None:
        checks.append(_check("structural_gate", False, structural_error))
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
        )
    checks.append(_check(
        "structural_gate",
        structural.passed,
        structural.reason_code,
    ))
    if not structural.passed:
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
            structural_decision=structural,
        )

    tandem, tandem_error = _resolve_gate(
        GATE_TANDEM,
        view,
        tandem_callback,
        tandem_decision,
    )
    if tandem_error is not None:
        checks.append(_check("tandem_gate", False, tandem_error))
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
            structural_decision=structural,
        )
    checks.append(_check("tandem_gate", tandem.passed, tandem.reason_code))
    if not tandem.passed:
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
            structural_decision=structural,
            tandem_decision=tandem,
        )

    if scoring_callback is None:
        checks.append(_check(
            "new_join_score",
            False,
            "new_join_scoring_callback_missing",
        ))
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
            structural_decision=structural,
            tandem_decision=tandem,
        )
    try:
        score = scoring_callback(view)
    except Exception:
        score = None
        score_error = "new_join_scoring_callback_failed"
    else:
        score_error = None
        if not isinstance(score, CircularJoinScoreEvidence):
            score_error = "new_join_score_evidence_invalid"
        elif score.subject_view_sha256 != view.view_sha256:
            score_error = "new_join_score_stale_view"
    if score_error is not None:
        checks.append(_check("new_join_score", False, score_error))
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
            structural_decision=structural,
            tandem_decision=tandem,
        )

    score_passes = float(score.score) >= float(policy.possible_score_threshold)
    checks.append(_check(
        "new_join_score",
        score_passes,
        "new_join_score_below_possible_threshold",
    ))
    if not score_passes:
        return _result(
            STATUS_REJECTED,
            source_hash,
            topology_evidence,
            compatibility_evidence,
            policy,
            checks,
            child_geometries,
            proposal,
            structural_decision=structural,
            tandem_decision=tandem,
            score_evidence=score,
        )

    return _result(
        STATUS_JOINED,
        source_hash,
        topology_evidence,
        compatibility_evidence,
        policy,
        checks,
        child_geometries,
        proposal,
        structural_decision=structural,
        tandem_decision=tandem,
        score_evidence=score,
        success_reason="source_replayed_circular_join_accepted_as_possible",
    )
