"""Fail-closed category transaction for regular diagnostic subpaths.

The structural planner in :mod:`source_regular_subpath_recovery` deliberately
does not classify its proposals.  This module is the production boundary that
adds the missing safety gates without weakening that separation:

* every accepted row, reviewed diagnostic parent, family member, and proposal
  is replayed against the original source;
* supplied coding-tandem decisions are compared with fresh recomputation;
* overlapping proposals are rejected as complete collision components rather
  than ranked by score or iteration order;
* a scorer is called only after geometry planning, and its score cannot promote
  a recovered path beyond ``Possible``; and
* mutations are copy-on-write and occur only after every callback result and
  stale-state check has passed.

Coordinates are zero-based and end-exclusive.  Row and candidate membership is
identity-exact.  Public audit records contain no live row or candidate objects.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from numbers import Integral, Real
from typing import Optional, Tuple

from components.accepted_family_tandem_validator import (
    AcceptedFamilyTandemDecision,
    validate_accepted_family_tandem,
)
from components.architecture_rescue_coding_context import SourceCDS
from components.boundary_assessment import (
    AlternativeRelationshipInput,
    BoundaryAssessmentNormalizer,
    BoundaryDimensionInput,
    DIMENSION_OCCURRENCE_PATH,
    STATUS_POLICY_SELECTED,
    geometry_input_from_candidate,
)
from components.reconstruction_validator import (
    ReconstructionValidationError,
    validate_candidate_reconstruction,
    validate_reconstruction,
)
from components.source_regular_subpath_recovery import (
    BOUNDARY_AMBIGUOUS_PARENT_PHASE,
    DIAGNOSTIC_PARENT,
    SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY,
    SourceRegularSubpathContext,
    SourceRegularSubpathProposal,
    propose_source_regular_subpaths_in_context,
)


RULE_NAME = "source-regular-subpath-recovery-transaction"
RULE_VERSION = "1.2.0"
CLASSIFICATION_METHOD = "source_backed_regular_subpath_recovery"
COORDINATE_SYSTEM = "zero_based_half_open"
REPORT_REPLAY_SNAPSHOT_SCHEMA = (
    "crispridentify.source_regular_subpath_report_replay_snapshot"
)
REPORT_REPLAY_SNAPSHOT_VERSION = "1.0.0"
ROW_VALUE_NORMALIZER_VERSION = "1.0.0"

ACTION_ABSTAIN = "abstain"
ACTION_REJECT = "reject"
ACTION_SCORING_REQUIRED = "scoring_required"
ACTION_ADMIT_POSSIBLE = "admit_possible"

POSSIBLE_CATEGORY_INDEX = 2
POSSIBLE_CATEGORY_NAME = "Possible"

_CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
_ACCEPTED_CATEGORY_INDEXES = (0, 2)
_DIAGNOSTIC_CATEGORY_INDEXES = (1, 3, 4)


def _json_copy(value, label):
    try:
        encoded = json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        return json.loads(encoded)
    except (TypeError, ValueError) as error:
        raise ValueError("{} must be JSON-safe".format(label)) from error


def _json_key(value, label):
    try:
        return json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise ValueError("{} must be JSON-safe".format(label)) from error


def _canonical_hash(value):
    return hashlib.sha256(
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("ascii")
    ).hexdigest()


def _slice_hash(source, interval):
    return hashlib.sha256(
        source[interval[0]:interval[1]].encode("ascii")
    ).hexdigest()


def _source_hash(source):
    return hashlib.sha256(source.encode("ascii")).hexdigest()


def _finite_score(value, unit_interval=False):
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(float(value))
    ):
        return False
    return not unit_interval or 0.0 <= float(value) <= 1.0


def canonicalize_regular_subpath_row_value(value, label="row value"):
    """Return a deterministic JSON-safe copy of row features or metadata.

    NumPy is intentionally not imported.  Array-like values must expose a
    deterministic ``tolist`` conversion; mappings require string keys so JSON
    key coercion cannot create collisions.  No object identity or ``repr`` is
    admitted into the binding.
    """

    active = set()

    def normalize(item, path):
        if item is None or isinstance(item, (str, bool)):
            return item
        if isinstance(item, Integral):
            return int(item)
        if isinstance(item, Real):
            if not math.isfinite(float(item)):
                raise ValueError("{} must contain only finite numbers".format(path))
            return float(item)
        if isinstance(item, Mapping):
            identity = id(item)
            if identity in active:
                raise ValueError("{} contains a reference cycle".format(path))
            active.add(identity)
            try:
                normalized = {}
                keys = tuple(item)
                if any(not isinstance(key, str) for key in keys):
                    raise ValueError(
                        "{} mapping keys must be strings".format(path)
                    )
                for key in sorted(keys):
                    normalized[key] = normalize(
                        item[key],
                        "{}.{}".format(path, key),
                    )
                return normalized
            finally:
                active.remove(identity)
        if isinstance(item, (list, tuple)):
            identity = id(item)
            if identity in active:
                raise ValueError("{} contains a reference cycle".format(path))
            active.add(identity)
            try:
                return [
                    normalize(member, "{}[{}]".format(path, index))
                    for index, member in enumerate(item)
                ]
            finally:
                active.remove(identity)
        tolist = getattr(item, "tolist", None)
        if callable(tolist):
            identity = id(item)
            if identity in active:
                raise ValueError("{} contains a reference cycle".format(path))
            active.add(identity)
            try:
                converted = tolist()
                if converted is item:
                    raise ValueError(
                        "{} tolist conversion returned itself".format(path)
                    )
                return normalize(converted, path + ".tolist()")
            finally:
                active.remove(identity)
        raise ValueError(
            "{} contains unsupported value type {}".format(
                path,
                type(item).__name__,
            )
        )

    normalized = normalize(value, label)
    _json_key(normalized, label)
    return normalized


def regular_subpath_row_value_sha256(value, label="row value"):
    """Return the shared canonical SHA-256 for a row value binding."""

    return _canonical_hash(
        canonicalize_regular_subpath_row_value(value, label)
    )


def _complete_family_record(geometries):
    return {
        "caller_declared_complete": True,
        "candidate_count": len(geometries),
        "candidate_geometries": [
            geometry.as_dict() for geometry in geometries
        ],
    }


def _structural_planner_record(
    status,
    reason_codes,
    proposal=None,
    proposal_json=None,
):
    if proposal is not None and proposal_json is not None:
        raise ValueError("structural proposal has two representations")
    if proposal_json is not None:
        proposal_record = json.loads(proposal_json)
    else:
        proposal_record = (
            None
            if proposal is None
            else _json_copy(proposal.as_dict(), "structural proposal")
        )
    return {
        "status": status,
        "reason_codes": list(reason_codes),
        "policy": SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY.as_dict(),
        "proposal": proposal_record,
    }


@dataclass(frozen=True)
class RegularSubpathGeometry:
    """Source-replayed geometry with a deterministic runtime-free digest."""

    source_interval: Tuple[int, int]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    repeat_starts: Tuple[int, ...]
    repeat_sequences: Tuple[str, ...]
    spacer_sequences: Tuple[str, ...]
    consensus: str
    source_slice_sha256: str
    geometry_sha256: str

    @property
    def repeat_count(self):
        return len(self.repeat_intervals)

    def as_dict(self):
        return {
            "coordinate_system": COORDINATE_SYSTEM,
            "source_interval": list(self.source_interval),
            "repeat_count": self.repeat_count,
            "repeat_intervals": [
                list(interval) for interval in self.repeat_intervals
            ],
            "repeat_starts": list(self.repeat_starts),
            "repeat_sequences": list(self.repeat_sequences),
            "spacer_sequences": list(self.spacer_sequences),
            "consensus": self.consensus,
            "source_slice_sha256": self.source_slice_sha256,
            "geometry_sha256": self.geometry_sha256,
        }


def _candidate_state(candidate):
    try:
        starts = tuple(candidate.list_repeat_starts)
        repeats = tuple(candidate.list_repeats)
        spacers = tuple(candidate.list_spacers)
        consensus = candidate.consensus
    except (AttributeError, TypeError) as error:
        raise TypeError(
            "candidate does not expose standard CRISPRidentify fields"
        ) from error
    if not isinstance(consensus, str) or not consensus:
        raise ValueError("candidate consensus must be a non-empty string")
    return starts, repeats, spacers, consensus


def _source_geometry(source, candidate, label):
    state = _candidate_state(candidate)
    starts, claimed_repeats, claimed_spacers, consensus = state
    reconstruction = None
    if (
        starts
        and len(starts) == len(claimed_repeats)
        and all(
            isinstance(start, int) and not isinstance(start, bool)
            for start in starts
        )
        and all(
            isinstance(repeat, str) and repeat
            for repeat in claimed_repeats
        )
    ):
        window_start = starts[0]
        window_end = starts[-1] + len(claimed_repeats[-1])
        if 0 <= window_start < window_end <= len(source):
            reconstruction = validate_reconstruction(
                source[window_start:window_end],
                tuple(start - window_start for start in starts),
                claimed_repeats,
                claimed_spacers,
            )
            if not reconstruction.valid:
                # Preserve the validator's original-source coordinates and
                # issue details for malformed claims.  Valid common-case rows
                # avoid copying and uppercasing the whole chromosome.
                reconstruction = None
    try:
        if reconstruction is None:
            reconstruction = validate_candidate_reconstruction(source, candidate)
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError("{} is not source-reconstructable".format(label)) from error
    if not reconstruction.valid:
        raise ReconstructionValidationError(reconstruction)
    repeats = tuple(
        source[start:start + len(repeat)]
        for start, repeat in zip(starts, claimed_repeats)
    )
    spacers = tuple(
        source[start + len(repeat):next_start]
        for start, repeat, next_start in zip(starts, repeats, starts[1:])
    )
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    source_interval = (
        repeat_intervals[0][0], repeat_intervals[-1][1]
    )
    # Keep this canonical payload identical to the stable planner's tandem
    # veto digest so independently supplied evidence binds to one geometry.
    record = {
        "coordinate_system": COORDINATE_SYSTEM,
        "repeat_starts": starts,
        "repeats": repeats,
        "spacers": spacers,
        "consensus": consensus,
    }
    return RegularSubpathGeometry(
        source_interval=source_interval,
        repeat_intervals=repeat_intervals,
        repeat_starts=starts,
        repeat_sequences=repeats,
        spacer_sequences=spacers,
        consensus=consensus,
        source_slice_sha256=_slice_hash(source, source_interval),
        geometry_sha256=_canonical_hash(record),
    )


@dataclass(frozen=True)
class RegularSubpathProposalTandemEvidence:
    """Caller-supplied tandem decision bound to one proposal digest."""

    candidate_geometry_sha256: str
    decision: object

    def __post_init__(self):
        if (
            not isinstance(self.candidate_geometry_sha256, str)
            or not self.candidate_geometry_sha256
        ):
            raise ValueError(
                "candidate_geometry_sha256 must be a non-empty string"
            )


@dataclass(frozen=True)
class SourceRegularSubpathDiagnosticRowContext:
    """Complete family and tandem contract for one diagnostic row."""

    context_key: object
    diagnostic_row: object
    target_locus_key: object
    complete_family_candidates: tuple
    supplied_parent_tandem_decision: object
    supplied_proposal_tandem_evidence: tuple
    family_complete: bool = True

    def __post_init__(self):
        if not isinstance(self.complete_family_candidates, tuple):
            raise TypeError("complete_family_candidates must be a tuple")
        if not isinstance(self.supplied_proposal_tandem_evidence, tuple):
            raise TypeError(
                "supplied_proposal_tandem_evidence must be a tuple"
            )
        for index, evidence in enumerate(
            self.supplied_proposal_tandem_evidence
        ):
            if not isinstance(evidence, RegularSubpathProposalTandemEvidence):
                raise TypeError(
                    "supplied_proposal_tandem_evidence[{}] has the wrong "
                    "type".format(index)
                )
        if not isinstance(self.family_complete, bool):
            raise TypeError("family_complete must be a boolean")


@dataclass(frozen=True)
class RegularSubpathOccurrencePathRelationship:
    """Explicitly path-local boundary claim, never a biological phase claim."""

    relationship: str
    parent_interval: Tuple[int, int]
    proposal_interval: Tuple[int, int]
    parent_repeat_index_span: Tuple[int, int]
    shared_repeat_intervals: Tuple[Tuple[int, int], ...]
    exact_parent_occurrence_path_subpath: bool
    biological_boundary_status: str
    exact_biological_repeat_phase: bool
    exact_array_boundaries: bool
    selection_basis: str
    planner_policy_name: str
    planner_policy_version: str

    def as_dict(self):
        value = asdict(self)
        for field_name in (
            "parent_interval",
            "proposal_interval",
            "parent_repeat_index_span",
        ):
            value[field_name] = list(value[field_name])
        value["shared_repeat_intervals"] = [
            list(interval) for interval in self.shared_repeat_intervals
        ]
        return value


@dataclass(frozen=True)
class RegularSubpathTandemAudit:
    """Supplied-versus-fresh tandem evidence for one exact geometry."""

    subject: str
    candidate_geometry_sha256: str
    supplied: bool
    recomputed: bool
    matched: bool
    supported: Optional[bool]
    decision: Optional[str]
    reason: Optional[str]
    supplied_decision: Optional[AcceptedFamilyTandemDecision]
    recomputed_decision: Optional[AcceptedFamilyTandemDecision]

    def as_dict(self):
        return {
            "subject": self.subject,
            "candidate_geometry_sha256": self.candidate_geometry_sha256,
            "supplied": self.supplied,
            "recomputed": self.recomputed,
            "matched": self.matched,
            "supported": self.supported,
            "decision": self.decision,
            "reason": self.reason,
            "supplied_decision": (
                None
                if self.supplied_decision is None
                else self.supplied_decision.as_dict()
            ),
            "recomputed_decision": (
                None
                if self.recomputed_decision is None
                else self.recomputed_decision.as_dict()
            ),
        }


@dataclass(frozen=True)
class RegularSubpathRecoveryScoreRequest:
    """One exact candidate requiring a fresh classifier evaluation."""

    request_id: str
    candidate: object
    source_sha256: str
    parent_locator_json: str
    target_locus_key_json: str
    parent_geometry: RegularSubpathGeometry
    proposal_geometry: RegularSubpathGeometry
    occurrence_relationship: RegularSubpathOccurrencePathRelationship

    @property
    def parent_locator(self):
        return json.loads(self.parent_locator_json)

    @property
    def target_locus_key(self):
        return json.loads(self.target_locus_key_json)

    def as_dict(self):
        return {
            "request_id": self.request_id,
            "source_sha256": self.source_sha256,
            "parent_locator": self.parent_locator,
            "target_locus_key": self.target_locus_key,
            "parent_geometry": self.parent_geometry.as_dict(),
            "proposal_geometry": self.proposal_geometry.as_dict(),
            "occurrence_relationship": self.occurrence_relationship.as_dict(),
        }


@dataclass(frozen=True)
class RegularSubpathFreshScoreResult:
    """Fresh standard scored row for one exact request candidate."""

    request_id: str
    scored_row: object
    evidence_method: str

    def __post_init__(self):
        if not isinstance(self.request_id, str) or not self.request_id:
            raise ValueError("request_id must be a non-empty string")
        if not isinstance(self.evidence_method, str) or not self.evidence_method:
            raise ValueError("evidence_method must be a non-empty string")


@dataclass(frozen=True)
class _RowLocation:
    category_index: int
    locus_key: object
    locus_key_json: str
    row_ordinal: int
    row: object


@dataclass(frozen=True)
class _SnapshotRow:
    row: object
    row_length: int
    score: float
    candidate: object
    candidate_state: tuple
    features: object
    metadata: object
    metadata_json: Optional[str]


@dataclass(frozen=True)
class RegularSubpathRecoveryDraft:
    """One planned proposal or one parent-level abstention."""

    action: str
    reason_codes: Tuple[str, ...]
    context_key_json: str
    parent_locator_json: str
    target_locus_key_json: str
    parent_score: float
    parent_geometry: RegularSubpathGeometry
    proposal_id: Optional[str]
    proposal_geometry: Optional[RegularSubpathGeometry]
    occurrence_relationship: Optional[RegularSubpathOccurrencePathRelationship]
    parent_tandem_audit: RegularSubpathTandemAudit
    proposal_tandem_audit: Optional[RegularSubpathTandemAudit]
    collision_peer_proposal_ids: Tuple[str, ...]
    score_request: Optional[RegularSubpathRecoveryScoreRequest]
    planner_status: str
    planner_reason_codes: Tuple[str, ...]
    _runtime_parent_row: object
    _runtime_target_locus_key: object
    _runtime_family_candidates: tuple
    _runtime_family_geometries: tuple
    _runtime_proposal: Optional[SourceRegularSubpathProposal]
    _runtime_supplied_parent_tandem_decision: object
    _runtime_supplied_proposal_tandem_decision: object

    @property
    def context_key(self):
        return json.loads(self.context_key_json)

    @property
    def parent_locator(self):
        return json.loads(self.parent_locator_json)

    @property
    def target_locus_key(self):
        return json.loads(self.target_locus_key_json)

    def as_dict(self):
        return {
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "context_key": self.context_key,
            "parent_locator": self.parent_locator,
            "target_locus_key": self.target_locus_key,
            "parent_score": self.parent_score,
            "parent_geometry": self.parent_geometry.as_dict(),
            "proposal_id": self.proposal_id,
            "proposal_geometry": (
                None
                if self.proposal_geometry is None
                else self.proposal_geometry.as_dict()
            ),
            "occurrence_relationship": (
                None
                if self.occurrence_relationship is None
                else self.occurrence_relationship.as_dict()
            ),
            "parent_tandem_audit": self.parent_tandem_audit.as_dict(),
            "proposal_tandem_audit": (
                None
                if self.proposal_tandem_audit is None
                else self.proposal_tandem_audit.as_dict()
            ),
            "collision_peer_proposal_ids": list(
                self.collision_peer_proposal_ids
            ),
            "score_request_id": (
                None if self.score_request is None else self.score_request.request_id
            ),
            "planner_status": self.planner_status,
            "planner_reason_codes": list(self.planner_reason_codes),
            "complete_family": _complete_family_record(
                self._runtime_family_geometries
            ),
            "structural_planner": _structural_planner_record(
                self.planner_status,
                self.planner_reason_codes,
                self._runtime_proposal,
            ),
        }


@dataclass(frozen=True)
class SourceRegularSubpathRecoveryBatch:
    """Frozen bridge between structural/tandem planning and fresh scoring."""

    source_sha256: str
    source_cds: Tuple[SourceCDS, ...]
    drafts: Tuple[RegularSubpathRecoveryDraft, ...]
    score_requests: Tuple[RegularSubpathRecoveryScoreRequest, ...]
    accepted_geometry_sha256s: Tuple[str, ...]
    parent_context_count: int
    report_replay_snapshot_json: str
    _category_snapshot: tuple
    _source_cds_snapshot: tuple
    _runtime_accepted_candidates: tuple
    _runtime_accepted_geometries: tuple

    @property
    def proposed_candidates(self):
        return tuple(request.candidate for request in self.score_requests)

    @property
    def report_replay_snapshot(self):
        return json.loads(self.report_replay_snapshot_json)

    def as_dict(self):
        return {
            "source_sha256": self.source_sha256,
            "source_cds": [asdict(cds) for cds in self.source_cds],
            "parent_context_count": self.parent_context_count,
            "report_replay_snapshot": self.report_replay_snapshot,
            "accepted_geometry_sha256s": list(
                self.accepted_geometry_sha256s
            ),
            "accepted_geometries": [
                geometry.as_dict()
                for geometry in self._runtime_accepted_geometries
            ],
            "drafts": [draft.as_dict() for draft in self.drafts],
            "score_requests": [
                request.as_dict() for request in self.score_requests
            ],
        }


@dataclass(frozen=True)
class RegularSubpathRecoveryRecord:
    """Runtime-free final decision for one planned proposal or parent."""

    action: str
    reason_codes: Tuple[str, ...]
    source_sha256: str
    context_key_json: str
    parent_locator_json: str
    target_locus_key_json: str
    proposal_id: Optional[str]
    request_id: Optional[str]
    parent_score: float
    fresh_score: Optional[float]
    fresh_evidence_method: Optional[str]
    parent_geometry: RegularSubpathGeometry
    proposal_geometry: Optional[RegularSubpathGeometry]
    occurrence_relationship: Optional[RegularSubpathOccurrencePathRelationship]
    parent_tandem_audit: RegularSubpathTandemAudit
    proposal_tandem_audit: Optional[RegularSubpathTandemAudit]
    complete_family_geometries: Tuple[RegularSubpathGeometry, ...]
    accepted_geometry_sha256s: Tuple[str, ...]
    planner_status: str
    planner_reason_codes: Tuple[str, ...]
    structural_proposal_json: Optional[str]
    report_replay_snapshot_sha256: str
    boundary_assessment_id: Optional[str]
    collision_peer_proposal_ids: Tuple[str, ...]

    @property
    def context_key(self):
        return json.loads(self.context_key_json)

    @property
    def parent_locator(self):
        return json.loads(self.parent_locator_json)

    @property
    def target_locus_key(self):
        return json.loads(self.target_locus_key_json)

    def as_dict(self):
        return {
            "rule": {"name": RULE_NAME, "version": RULE_VERSION},
            "source_sha256": self.source_sha256,
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "context_key": self.context_key,
            "parent_locator": self.parent_locator,
            "target_locus_key": self.target_locus_key,
            "proposal_id": self.proposal_id,
            "request_id": self.request_id,
            "parent_score": self.parent_score,
            "fresh_score": self.fresh_score,
            "fresh_evidence_method": self.fresh_evidence_method,
            "parent_geometry": self.parent_geometry.as_dict(),
            "proposal_geometry": (
                None
                if self.proposal_geometry is None
                else self.proposal_geometry.as_dict()
            ),
            "occurrence_relationship": (
                None
                if self.occurrence_relationship is None
                else self.occurrence_relationship.as_dict()
            ),
            "parent_tandem_audit": self.parent_tandem_audit.as_dict(),
            "proposal_tandem_audit": (
                None
                if self.proposal_tandem_audit is None
                else self.proposal_tandem_audit.as_dict()
            ),
            "complete_family": _complete_family_record(
                self.complete_family_geometries
            ),
            "accepted_overlap_protection": {
                "accepted_rows_identity_complete": True,
                "accepted_geometry_sha256s": list(
                    self.accepted_geometry_sha256s
                ),
            },
            "structural_planner": _structural_planner_record(
                self.planner_status,
                self.planner_reason_codes,
                proposal_json=self.structural_proposal_json,
            ),
            "report_replay_snapshot_sha256": (
                self.report_replay_snapshot_sha256
            ),
            "boundary_assessment_id": self.boundary_assessment_id,
            "collision_peer_proposal_ids": list(
                self.collision_peer_proposal_ids
            ),
            "classification_cap": POSSIBLE_CATEGORY_NAME,
            "bona_fide_admission_allowed": False,
        }


@dataclass(frozen=True)
class SourceRegularSubpathRecoveryStats:
    total_rows_before: int
    total_rows_after: int
    accepted_rows_before: int
    accepted_rows_after: int
    possible_rows_before: int
    possible_rows_after: int
    diagnostic_rows_before: int
    diagnostic_rows_after: int
    parents_reviewed: int
    structural_proposals: int
    collision_rejections: int
    tandem_missing_rejections: int
    tandem_stale_rejections: int
    tandem_supported_rejections: int
    scoring_requests: int
    scorer_calls: int
    possible_admissions: int
    bona_fide_admissions: int

    def as_dict(self):
        return {
            field_name: getattr(self, field_name)
            for field_name in self.__dataclass_fields__
        }


def _validate_row(row, accepted=False, label="category row"):
    if not isinstance(row, (list, tuple)) or len(row) not in (3, 4):
        raise ValueError(
            "{} must contain score, candidate, features, and optional "
            "metadata".format(label)
        )
    if not _finite_score(row[0], unit_interval=accepted):
        suffix = " in [0, 1]" if accepted else ""
        raise ValueError("{} score must be finite{}".format(label, suffix))
    if row[1] is None:
        raise ValueError("{} must contain a candidate".format(label))
    if isinstance(row[2], (str, bytes)) or not isinstance(row[2], (list, tuple)):
        raise ValueError("{} features must be a sequence".format(label))
    if len(row) == 4:
        if not isinstance(row[3], Mapping):
            raise ValueError("{} metadata must be a mapping".format(label))
        _json_copy(row[3], "{} metadata".format(label))


def _snapshot_row(row):
    metadata = row[3] if len(row) == 4 else None
    return _SnapshotRow(
        row=row,
        row_length=len(row),
        score=float(row[0]),
        candidate=row[1],
        candidate_state=_candidate_state(row[1]),
        features=row[2],
        metadata=metadata,
        metadata_json=(
            None if metadata is None else _json_key(metadata, "row metadata")
        ),
    )


def _index_categories(categories):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("regular-subpath recovery requires five categories")
    locations = []
    accepted = []
    diagnostics = []
    snapshot = []
    row_ids = set()
    candidate_ids = set()
    for category_index, category in enumerate(categories):
        if not isinstance(category, Mapping):
            raise ValueError("each category must be a mapping")
        entries = []
        for locus_key, rows in category.items():
            locus_json = _json_key(locus_key, "category locus key")
            if not isinstance(rows, (list, tuple)) or not rows:
                raise ValueError(
                    "category locus values must be non-empty row sequences"
                )
            row_snapshots = []
            for ordinal, row in enumerate(rows):
                is_accepted = category_index in _ACCEPTED_CATEGORY_INDEXES
                _validate_row(row, accepted=is_accepted)
                if id(row) in row_ids:
                    raise ValueError("one row occurs more than once in categories")
                if id(row[1]) in candidate_ids:
                    raise ValueError(
                        "one candidate object has multiple category rows"
                    )
                row_ids.add(id(row))
                candidate_ids.add(id(row[1]))
                location = _RowLocation(
                    category_index=category_index,
                    locus_key=locus_key,
                    locus_key_json=locus_json,
                    row_ordinal=ordinal,
                    row=row,
                )
                locations.append(location)
                if is_accepted:
                    accepted.append(location)
                else:
                    diagnostics.append(location)
                row_snapshots.append(_snapshot_row(row))
            entries.append((locus_key, tuple(row_snapshots)))
        snapshot.append(tuple(entries))
    return tuple(locations), tuple(accepted), tuple(diagnostics), tuple(snapshot)


def _validate_snapshot(categories, snapshot):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("categories changed after regular-subpath planning")
    for category_index, expected_entries in enumerate(snapshot):
        category = categories[category_index]
        if not isinstance(category, Mapping) or len(category) != len(expected_entries):
            raise ValueError("categories changed after regular-subpath planning")
        for locus_key, expected_rows in expected_entries:
            if locus_key not in category:
                raise ValueError("a category locus changed after planning")
            observed_rows = category[locus_key]
            if (
                not isinstance(observed_rows, (list, tuple))
                or len(observed_rows) != len(expected_rows)
            ):
                raise ValueError("category row cardinality changed after planning")
            for observed, expected in zip(observed_rows, expected_rows):
                metadata = observed[3] if len(observed) == 4 else None
                if (
                    observed is not expected.row
                    or len(observed) != expected.row_length
                    or float(observed[0]) != expected.score
                    or observed[1] is not expected.candidate
                    or _candidate_state(observed[1]) != expected.candidate_state
                    or observed[2] is not expected.features
                    or metadata is not expected.metadata
                    or (
                        None
                        if metadata is None
                        else _json_key(metadata, "row metadata")
                    ) != expected.metadata_json
                ):
                    raise ValueError("a category row changed after planning")


def _validate_source_context(context):
    if not isinstance(context, SourceRegularSubpathContext):
        raise TypeError("context must be SourceRegularSubpathContext")
    if context.source_coordinate_offset != 0:
        raise ValueError(
            "production tandem validation requires a full-source context "
            "with offset zero"
        )
    observed_sha256 = _source_hash(context.source_sequence)
    if observed_sha256 != context.source_sha256:
        raise ValueError("source context sequence and hash disagree")
    existence_context = context.existence_context
    if (
        existence_context.source_sequence != context.source_sequence
        or existence_context.source_sha256 != observed_sha256
    ):
        raise ValueError("source and existence contexts disagree")


def _validate_source_cds(source, source_cds):
    if not isinstance(source_cds, tuple):
        raise TypeError("source_cds must be an explicit tuple")
    seen = set()
    for index, cds in enumerate(source_cds):
        if not isinstance(cds, SourceCDS):
            raise TypeError("source_cds[{}] is not SourceCDS".format(index))
        if cds.end > len(source):
            raise ValueError("source_cds[{}] lies outside source".format(index))
        if cds in seen:
            raise ValueError("source_cds must not contain duplicates")
        seen.add(cds)
    return tuple(sorted(
        source_cds,
        key=lambda cds: (
            cds.start,
            cds.end,
            cds.is_complement,
            cds.confidence,
        ),
    ))


def _index_accepted_rows(accepted_rows, accepted_locations):
    if not isinstance(accepted_rows, tuple):
        raise TypeError("accepted_rows must be an explicit tuple")
    if len({id(row) for row in accepted_rows}) != len(accepted_rows):
        raise ValueError("accepted_rows contains duplicate row identities")
    expected = {id(location.row): location.row for location in accepted_locations}
    observed = {id(row): row for row in accepted_rows}
    if set(observed) != set(expected) or any(
        observed[row_id] is not expected[row_id] for row_id in expected
    ):
        raise ValueError(
            "accepted_rows must identity-match every Bona-fide and Possible row"
        )
    return tuple(location.row[1] for location in accepted_locations)


def _index_diagnostic_contexts(contexts, diagnostic_locations):
    if not isinstance(contexts, tuple):
        raise TypeError("diagnostic_contexts must be an explicit tuple")
    by_row_id = {id(location.row): location for location in diagnostic_locations}
    indexed = []
    seen_rows = set()
    seen_keys = set()
    for index, context in enumerate(contexts):
        if not isinstance(context, SourceRegularSubpathDiagnosticRowContext):
            raise TypeError(
                "diagnostic_contexts[{}] has the wrong type".format(index)
            )
        if not context.family_complete:
            raise ValueError("every diagnostic family must be explicitly complete")
        location = by_row_id.get(id(context.diagnostic_row))
        if location is None or location.row is not context.diagnostic_row:
            raise ValueError(
                "diagnostic context refers to an accepted or foreign row"
            )
        if location.category_index not in _DIAGNOSTIC_CATEGORY_INDEXES:
            raise RuntimeError("diagnostic category indexing changed")
        if id(context.diagnostic_row) in seen_rows:
            raise ValueError("one diagnostic row has duplicate contexts")
        context_key_json = _json_key(context.context_key, "context_key")
        target_json = _json_key(context.target_locus_key, "target_locus_key")
        try:
            hash(context.target_locus_key)
        except TypeError as error:
            raise ValueError("target_locus_key must be hashable") from error
        if context_key_json in seen_keys:
            raise ValueError("diagnostic context keys must be unique")
        seen_rows.add(id(context.diagnostic_row))
        seen_keys.add(context_key_json)
        indexed.append((context, location, context_key_json, target_json))
    return tuple(indexed)


def _validate_family(source, parent_candidate, family_candidates):
    if not family_candidates:
        raise ValueError("complete_family_candidates must be non-empty")
    if len({id(item) for item in family_candidates}) != len(family_candidates):
        raise ValueError("complete_family_candidates contains duplicate objects")
    parent_occurrences = sum(item is parent_candidate for item in family_candidates)
    if parent_occurrences != 1:
        raise ValueError(
            "complete_family_candidates must contain the parent exactly once"
        )
    geometries = tuple(
        _source_geometry(source, candidate, "family candidate {}".format(index))
        for index, candidate in enumerate(family_candidates)
    )
    hashes = tuple(geometry.geometry_sha256 for geometry in geometries)
    if len(set(hashes)) != len(hashes):
        raise ValueError("complete_family_candidates contains duplicate geometry")
    return geometries


def _tandem_audit(subject, geometry, supplied, recomputed):
    supplied_valid = isinstance(supplied, AcceptedFamilyTandemDecision)
    recomputed_valid = isinstance(recomputed, AcceptedFamilyTandemDecision)
    matched = supplied_valid and recomputed_valid and supplied == recomputed
    selected = recomputed if recomputed_valid else None
    return RegularSubpathTandemAudit(
        subject=subject,
        candidate_geometry_sha256=geometry.geometry_sha256,
        supplied=supplied_valid,
        recomputed=recomputed_valid,
        matched=matched,
        supported=(None if selected is None else selected.supported),
        decision=(None if selected is None else selected.decision),
        reason=(None if selected is None else selected.reason),
        supplied_decision=(supplied if supplied_valid else None),
        recomputed_decision=(recomputed if recomputed_valid else None),
    )


def _tandem_reason(audit, prefix):
    if not audit.supplied:
        return ("{}_tandem_evidence_missing".format(prefix),)
    if not audit.recomputed or not audit.matched:
        return ("{}_tandem_evidence_stale_or_mismatched".format(prefix),)
    if audit.supported:
        return ("{}_coding_tandem_supported".format(prefix),)
    return ()


def _parent_locator(location, context_key_json):
    return _json_key({
        "category_index": location.category_index,
        "category_name": _CATEGORY_NAMES[location.category_index],
        "locus_key": json.loads(location.locus_key_json),
        "row_ordinal": location.row_ordinal,
        "context_key": json.loads(context_key_json),
    }, "parent locator")


def _proposal_id(source_sha256, parent_locator_json, geometry_sha256):
    digest = _canonical_hash({
        "source_sha256": source_sha256,
        "parent_locator": json.loads(parent_locator_json),
        "proposal_geometry_sha256": geometry_sha256,
    })
    return "source-regular-subpath-" + digest[:24]


def _relationship(parent_geometry, proposal, proposal_geometry):
    start_index, end_index = proposal.parent_index_span
    shared = parent_geometry.repeat_intervals[start_index:end_index]
    if shared != proposal_geometry.repeat_intervals:
        raise RuntimeError("planner proposal is not an exact parent occurrence subpath")
    if proposal.boundary_status != BOUNDARY_AMBIGUOUS_PARENT_PHASE:
        raise RuntimeError("planner made an unsupported boundary-phase claim")
    return RegularSubpathOccurrencePathRelationship(
        relationship="policy_selected_occurrence_path_subpath",
        parent_interval=parent_geometry.source_interval,
        proposal_interval=proposal_geometry.source_interval,
        parent_repeat_index_span=proposal.parent_index_span,
        shared_repeat_intervals=shared,
        exact_parent_occurrence_path_subpath=True,
        biological_boundary_status=BOUNDARY_AMBIGUOUS_PARENT_PHASE,
        exact_biological_repeat_phase=False,
        exact_array_boundaries=False,
        selection_basis=(
            "maximal_source_path_segment_between_extreme_spacer_breaks"
        ),
        planner_policy_name=SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY.name,
        planner_policy_version=SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY.version,
    )


def _collision_peers(drafts):
    # Interval graphs have a direct sweep decomposition into connected
    # components.  Avoid the former all-pairs scan for the overwhelmingly
    # common case of many disjoint diagnostic loci.  Materializing every peer
    # remains necessarily quadratic only inside a real collision component,
    # because the public record explicitly lists all of those peers.
    proposal_drafts = sorted(
        (draft for draft in drafts if draft.proposal_id),
        key=lambda draft: (
            draft.proposal_geometry.source_interval[0],
            draft.proposal_geometry.source_interval[1],
            draft.proposal_id,
        ),
    )
    peers = {}

    def record_component(component):
        if len(component) < 2:
            return
        proposal_ids = tuple(sorted(
            draft.proposal_id for draft in component
        ))
        for proposal_id in proposal_ids:
            peers[proposal_id] = tuple(
                peer for peer in proposal_ids if peer != proposal_id
            )

    component = []
    component_max_end = None
    for draft in proposal_drafts:
        start, end = draft.proposal_geometry.source_interval
        if component and start >= component_max_end:
            record_component(component)
            component = []
            component_max_end = None
        component.append(draft)
        component_max_end = (
            end if component_max_end is None else max(component_max_end, end)
        )
    record_component(component)
    return peers


def _canonical_row_manifest_entry(source, location, geometry=None):
    if geometry is None:
        geometry = _source_geometry(
            source,
            location.row[1],
            "pretransaction category row",
        )
    features = canonicalize_regular_subpath_row_value(
        location.row[2],
        "pretransaction row features",
    )
    metadata = (
        None
        if len(location.row) == 3
        else canonicalize_regular_subpath_row_value(
            location.row[3],
            "pretransaction row metadata",
        )
    )
    features_sha256 = regular_subpath_row_value_sha256(
        features,
        "canonical row features",
    )
    metadata_sha256 = regular_subpath_row_value_sha256(
        metadata,
        "canonical row metadata",
    )
    record = {
        "category_index": location.category_index,
        "category_name": _CATEGORY_NAMES[location.category_index],
        "locus_key": json.loads(location.locus_key_json),
        "row_ordinal": location.row_ordinal,
        "row_field_count": len(location.row),
        "score": float(location.row[0]),
        "geometry": geometry.as_dict(),
        "features": features,
        "features_sha256": features_sha256,
        "metadata": metadata,
        "metadata_sha256": metadata_sha256,
    }
    record["row_binding_sha256"] = _canonical_hash(record)
    return record


def _reviewed_parent_replay_records(drafts):
    indexed = {}
    for draft in drafts:
        record = {
            "context_key": draft.context_key,
            "parent_locator": draft.parent_locator,
            "target_locus_key": draft.target_locus_key,
            "parent_geometry": draft.parent_geometry.as_dict(),
            "complete_family": _complete_family_record(
                draft._runtime_family_geometries
            ),
        }
        previous = indexed.get(draft.parent_locator_json)
        if previous is not None and previous != record:
            raise RuntimeError(
                "one reviewed parent has inconsistent replay state"
            )
        indexed[draft.parent_locator_json] = record
    return [indexed[key] for key in sorted(indexed)]


def _build_report_replay_snapshot(
    source,
    source_sha256,
    source_cds,
    locations,
    accepted_geometries,
    drafts,
    score_requests,
    row_geometry_by_id=None,
):
    row_geometry_by_id = (
        {}
        if row_geometry_by_id is None
        else row_geometry_by_id
    )
    ordered_locations = sorted(
        locations,
        key=lambda location: (
            location.category_index,
            location.locus_key_json,
            location.row_ordinal,
        ),
    )
    rows = [
        _canonical_row_manifest_entry(
            source,
            location,
            row_geometry_by_id.get(id(location.row)),
        )
        for location in ordered_locations
    ]
    accepted = sorted(
        (geometry.as_dict() for geometry in accepted_geometries),
        key=lambda geometry: (
            geometry["source_interval"],
            geometry["geometry_sha256"],
        ),
    )
    cds_records = [
        canonicalize_regular_subpath_row_value(
            asdict(cds),
            "source CDS",
        )
        for cds in source_cds
    ]
    reviewed_parents = _reviewed_parent_replay_records(drafts)
    draft_records = sorted(
        (_json_copy(draft.as_dict(), "recovery draft") for draft in drafts),
        key=lambda draft: (
            _json_key(draft["parent_locator"], "draft parent locator"),
            "" if draft["proposal_id"] is None else draft["proposal_id"],
        ),
    )
    request_records = sorted(
        (
            _json_copy(request.as_dict(), "score request")
            for request in score_requests
        ),
        key=lambda request: request["request_id"],
    )
    payload = {
        "schema": {
            "name": REPORT_REPLAY_SNAPSHOT_SCHEMA,
            "version": REPORT_REPLAY_SNAPSHOT_VERSION,
        },
        "transaction_rule": {
            "name": RULE_NAME,
            "version": RULE_VERSION,
        },
        "coordinate_system": COORDINATE_SYSTEM,
        "source_sha256": source_sha256,
        "source_length": len(source),
        "source_cds": cds_records,
        "accepted_geometries": accepted,
        "reviewed_diagnostic_parents": reviewed_parents,
        "pretransaction_row_manifest": rows,
        "transaction_plan": {
            "drafts": draft_records,
            "score_requests": request_records,
        },
        "bindings": {
            "row_value_normalizer_version": ROW_VALUE_NORMALIZER_VERSION,
            "source_cds_sha256": _canonical_hash(cds_records),
            "accepted_geometries_sha256": _canonical_hash(accepted),
            "reviewed_diagnostic_parents_sha256": _canonical_hash(
                reviewed_parents
            ),
            "pretransaction_row_manifest_sha256": _canonical_hash(rows),
            "transaction_plan_sha256": _canonical_hash({
                "drafts": draft_records,
                "score_requests": request_records,
            }),
        },
        "invariants": {
            "accepted_rows_are_complete_category_0_and_2_snapshot": True,
            "all_pretransaction_rows_source_replayed": True,
            "all_reviewed_parent_locators_canonical": True,
            "features_and_metadata_embedded_and_sha256_bound": True,
            "runtime_object_identity_serialized": False,
        },
    }
    payload["snapshot_sha256"] = _canonical_hash(payload)
    return _json_copy(payload, "report replay snapshot")


def _validate_report_replay_snapshot_json(value, expected_source_sha256):
    if not isinstance(value, str) or not value:
        raise ValueError("report replay snapshot must be canonical JSON")
    try:
        snapshot = json.loads(value)
    except (TypeError, ValueError) as error:
        raise ValueError("report replay snapshot is invalid JSON") from error
    if _json_key(snapshot, "report replay snapshot") != value:
        raise ValueError("report replay snapshot JSON is not canonical")
    if snapshot.get("schema") != {
        "name": REPORT_REPLAY_SNAPSHOT_SCHEMA,
        "version": REPORT_REPLAY_SNAPSHOT_VERSION,
    }:
        raise ValueError("report replay snapshot schema changed")
    if snapshot.get("source_sha256") != expected_source_sha256:
        raise ValueError("report replay snapshot source changed")
    if snapshot.get("transaction_rule") != {
        "name": RULE_NAME,
        "version": RULE_VERSION,
    }:
        raise ValueError("report replay snapshot transaction rule changed")
    supplied_hash = snapshot.get("snapshot_sha256")
    unsigned = dict(snapshot)
    unsigned.pop("snapshot_sha256", None)
    if supplied_hash != _canonical_hash(unsigned):
        raise ValueError("report replay snapshot digest changed")
    return snapshot


def _regular_subpath_boundary_assessment(normalizer, draft):
    source_hash = normalizer.source_sha256
    parent = draft._runtime_parent_row[1]
    proposal = draft._runtime_proposal.candidate
    return normalizer.normalize(
        geometry_input_from_candidate(proposal, source_hash),
        alternative_relationships=(
            AlternativeRelationshipInput(
                geometry=geometry_input_from_candidate(parent, source_hash),
                label="diagnostic-parent-occurrence-path",
                dimensions=(DIMENSION_OCCURRENCE_PATH,),
                reason_codes=(
                    "parent_source_replayed",
                    "proposal_is_policy_selected_parent_subpath",
                ),
            ),
        ),
        occurrence_path=BoundaryDimensionInput(
            status=STATUS_POLICY_SELECTED,
            selection_policy_id=(
                "source-regular-subpath-recovery:"
                + SOURCE_REGULAR_SUBPATH_RECOVERY_POLICY.version
            ),
            reason_codes=(
                "maximal_regular_segment_selected_between_extreme_breaks",
                "parent_occurrence_path_retained_as_alternative",
            ),
        ),
        repeat_phase=BoundaryDimensionInput(),
    )


def plan_source_regular_subpath_recovery_batch(
    context,
    categories,
    diagnostic_contexts,
    accepted_rows,
    source_cds,
):
    """Plan source-backed proposals and expose only fully gated score requests."""

    _validate_source_context(context)
    source = context.source_sequence
    source_cds = _validate_source_cds(source, source_cds)
    locations, accepted_locations, diagnostics, snapshot = _index_categories(
        categories
    )
    row_geometry_by_id = {
        id(location.row): _source_geometry(
            source,
            location.row[1],
            "pretransaction category row {}:{}:{}".format(
                location.category_index,
                location.locus_key_json,
                location.row_ordinal,
            ),
        )
        for location in locations
    }
    accepted_candidates = _index_accepted_rows(
        accepted_rows, accepted_locations
    )
    accepted_geometries = tuple(
        row_geometry_by_id[id(location.row)]
        for location in accepted_locations
    )
    indexed_contexts = _index_diagnostic_contexts(
        diagnostic_contexts, diagnostics
    )

    drafts = []
    for diagnostic, location, context_key_json, target_json in indexed_contexts:
        parent = location.row[1]
        parent_geometry = row_geometry_by_id[id(location.row)]
        family_geometries = _validate_family(
            source, parent, diagnostic.complete_family_candidates
        )
        parent_locator_json = _parent_locator(location, context_key_json)
        supplied_parent = diagnostic.supplied_parent_tandem_decision
        recomputed_parent = validate_accepted_family_tandem(
            source,
            parent,
            diagnostic.complete_family_candidates,
            source_cds,
        )
        parent_audit = _tandem_audit(
            "parent", parent_geometry, supplied_parent, recomputed_parent
        )
        parent_reasons = _tandem_reason(parent_audit, "parent")

        planner_result = propose_source_regular_subpaths_in_context(
            context,
            parent,
            parent_identifier=_canonical_hash(
                json.loads(parent_locator_json)
            )[:24],
            parent_disposition=DIAGNOSTIC_PARENT,
            accepted_candidates=accepted_candidates,
        )
        supplied_by_hash = {}
        duplicate_evidence = False
        for evidence in diagnostic.supplied_proposal_tandem_evidence:
            if evidence.candidate_geometry_sha256 in supplied_by_hash:
                duplicate_evidence = True
            supplied_by_hash[evidence.candidate_geometry_sha256] = evidence.decision

        proposal_hashes = {
            proposal.tandem_veto.candidate_geometry_sha256
            for proposal in planner_result.proposals
        }
        evidence_hashes = set(supplied_by_hash)
        foreign_evidence = duplicate_evidence or bool(
            evidence_hashes.difference(proposal_hashes)
        )

        if not planner_result.proposals:
            evidence_reasons = (
                ("proposal_tandem_evidence_set_mismatched",)
                if foreign_evidence
                else ()
            )
            reasons = (
                tuple(parent_reasons)
                + evidence_reasons
                + tuple(planner_result.reason_codes)
            )
            drafts.append(RegularSubpathRecoveryDraft(
                action=(
                    ACTION_REJECT
                    if parent_reasons or evidence_reasons
                    else ACTION_ABSTAIN
                ),
                reason_codes=reasons,
                context_key_json=context_key_json,
                parent_locator_json=parent_locator_json,
                target_locus_key_json=target_json,
                parent_score=float(location.row[0]),
                parent_geometry=parent_geometry,
                proposal_id=None,
                proposal_geometry=None,
                occurrence_relationship=None,
                parent_tandem_audit=parent_audit,
                proposal_tandem_audit=None,
                collision_peer_proposal_ids=(),
                score_request=None,
                planner_status=planner_result.status,
                planner_reason_codes=tuple(planner_result.reason_codes),
                _runtime_parent_row=location.row,
                _runtime_target_locus_key=diagnostic.target_locus_key,
                _runtime_family_candidates=diagnostic.complete_family_candidates,
                _runtime_family_geometries=family_geometries,
                _runtime_proposal=None,
                _runtime_supplied_parent_tandem_decision=supplied_parent,
                _runtime_supplied_proposal_tandem_decision=None,
            ))
            continue

        for proposal in planner_result.proposals:
            proposal_geometry = _source_geometry(
                source, proposal.candidate, "planner proposal"
            )
            if (
                proposal_geometry.geometry_sha256
                != proposal.tandem_veto.candidate_geometry_sha256
            ):
                raise RuntimeError("planner proposal digest disagrees with source replay")
            relationship = _relationship(
                parent_geometry, proposal, proposal_geometry
            )
            supplied_proposal = supplied_by_hash.get(
                proposal.tandem_veto.candidate_geometry_sha256
            )
            recomputed_proposal = validate_accepted_family_tandem(
                source,
                proposal.candidate,
                diagnostic.complete_family_candidates + (proposal.candidate,),
                source_cds,
            )
            proposal_audit = _tandem_audit(
                "proposal",
                proposal_geometry,
                supplied_proposal,
                recomputed_proposal,
            )
            reasons = list(parent_reasons)
            reasons.extend(_tandem_reason(proposal_audit, "proposal"))
            if foreign_evidence:
                reasons.append("proposal_tandem_evidence_set_mismatched")
            proposal_id = _proposal_id(
                context.source_sha256,
                parent_locator_json,
                proposal_geometry.geometry_sha256,
            )
            request = None
            if reasons:
                action = ACTION_REJECT
            else:
                action = ACTION_SCORING_REQUIRED
                request = RegularSubpathRecoveryScoreRequest(
                    request_id=proposal_id,
                    candidate=proposal.candidate,
                    source_sha256=context.source_sha256,
                    parent_locator_json=parent_locator_json,
                    target_locus_key_json=target_json,
                    parent_geometry=parent_geometry,
                    proposal_geometry=proposal_geometry,
                    occurrence_relationship=relationship,
                )
                reasons.append("fresh_scoring_required")
            drafts.append(RegularSubpathRecoveryDraft(
                action=action,
                reason_codes=tuple(reasons),
                context_key_json=context_key_json,
                parent_locator_json=parent_locator_json,
                target_locus_key_json=target_json,
                parent_score=float(location.row[0]),
                parent_geometry=parent_geometry,
                proposal_id=proposal_id,
                proposal_geometry=proposal_geometry,
                occurrence_relationship=relationship,
                parent_tandem_audit=parent_audit,
                proposal_tandem_audit=proposal_audit,
                collision_peer_proposal_ids=(),
                score_request=request,
                planner_status=planner_result.status,
                planner_reason_codes=tuple(planner_result.reason_codes),
                _runtime_parent_row=location.row,
                _runtime_target_locus_key=diagnostic.target_locus_key,
                _runtime_family_candidates=diagnostic.complete_family_candidates,
                _runtime_family_geometries=family_geometries,
                _runtime_proposal=proposal,
                _runtime_supplied_parent_tandem_decision=supplied_parent,
                _runtime_supplied_proposal_tandem_decision=supplied_proposal,
            ))

    if len({draft.proposal_id for draft in drafts if draft.proposal_id}) != sum(
        draft.proposal_id is not None for draft in drafts
    ):
        raise RuntimeError("regular-subpath proposal IDs are not unique")

    collision_peers = _collision_peers(drafts)
    collision_rewritten = []
    for draft in drafts:
        peers = collision_peers.get(draft.proposal_id, ())
        if peers:
            reasons = tuple(
                code
                for code in draft.reason_codes
                if code != "fresh_scoring_required"
            ) + ("proposal_interval_collision",)
            draft = replace(
                draft,
                action=ACTION_REJECT,
                reason_codes=reasons,
                collision_peer_proposal_ids=peers,
                score_request=None,
            )
        collision_rewritten.append(draft)
    drafts = collision_rewritten
    requests = tuple(
        draft.score_request
        for draft in drafts
        if draft.action == ACTION_SCORING_REQUIRED
    )
    if len({request.request_id for request in requests}) != len(requests):
        raise RuntimeError("regular-subpath score request IDs are not unique")

    report_replay_snapshot = _build_report_replay_snapshot(
        source,
        context.source_sha256,
        source_cds,
        locations,
        accepted_geometries,
        drafts,
        requests,
        row_geometry_by_id=row_geometry_by_id,
    )
    report_replay_snapshot_json = _json_key(
        report_replay_snapshot,
        "report replay snapshot",
    )

    batch = SourceRegularSubpathRecoveryBatch(
        source_sha256=context.source_sha256,
        source_cds=source_cds,
        drafts=tuple(drafts),
        score_requests=requests,
        accepted_geometry_sha256s=tuple(
            geometry.geometry_sha256 for geometry in accepted_geometries
        ),
        parent_context_count=len(indexed_contexts),
        report_replay_snapshot_json=report_replay_snapshot_json,
        _category_snapshot=snapshot,
        _source_cds_snapshot=tuple(
            (cds.start, cds.end, cds.is_complement, cds.confidence)
            for cds in source_cds
        ),
        _runtime_accepted_candidates=accepted_candidates,
        _runtime_accepted_geometries=accepted_geometries,
    )
    _json_copy(batch.as_dict(), "regular-subpath recovery batch")
    return batch


def _validate_batch_structure(batch):
    if not isinstance(batch, SourceRegularSubpathRecoveryBatch):
        raise TypeError("batch must be SourceRegularSubpathRecoveryBatch")
    _validate_report_replay_snapshot_json(
        batch.report_replay_snapshot_json,
        batch.source_sha256,
    )
    request_by_id = {}
    for request in batch.score_requests:
        if not isinstance(request, RegularSubpathRecoveryScoreRequest):
            raise TypeError("batch contains a malformed score request")
        if request.request_id in request_by_id:
            raise ValueError("batch contains duplicate score request IDs")
        request_by_id[request.request_id] = request
        if request.source_sha256 != batch.source_sha256:
            raise ValueError("score request source hash disagrees with batch")
    seen_proposals = set()
    referenced_requests = set()
    for draft in batch.drafts:
        if not isinstance(draft, RegularSubpathRecoveryDraft):
            raise TypeError("batch contains a malformed recovery draft")
        if draft.proposal_id is not None:
            if draft.proposal_id in seen_proposals:
                raise ValueError("batch contains duplicate proposal IDs")
            seen_proposals.add(draft.proposal_id)
            if draft.proposal_geometry is None:
                raise ValueError("proposal draft is missing source geometry")
            expected_proposal_id = _proposal_id(
                batch.source_sha256,
                draft.parent_locator_json,
                draft.proposal_geometry.geometry_sha256,
            )
            if draft.proposal_id != expected_proposal_id:
                raise ValueError("proposal ID changed after planning")
        elif (
            draft.proposal_geometry is not None
            or draft.occurrence_relationship is not None
            or draft._runtime_proposal is not None
        ):
            raise ValueError("parent-only draft retained proposal state")
        try:
            hash(draft._runtime_target_locus_key)
        except TypeError as error:
            raise ValueError(
                "runtime target locus key changed after planning"
            ) from error
        if _json_key(
            draft._runtime_target_locus_key,
            "runtime target_locus_key",
        ) != draft.target_locus_key_json:
            raise ValueError("runtime target locus key changed after planning")
        if draft.action == ACTION_SCORING_REQUIRED:
            if draft.score_request is None:
                raise ValueError("scoring draft is missing its request")
            request = request_by_id.get(draft.score_request.request_id)
            if request is not draft.score_request:
                raise ValueError("draft score request disagrees with batch")
            if (
                request.request_id != draft.proposal_id
                or request.parent_locator_json != draft.parent_locator_json
                or request.target_locus_key_json
                != draft.target_locus_key_json
                or request.parent_geometry != draft.parent_geometry
                or request.proposal_geometry != draft.proposal_geometry
                or request.occurrence_relationship
                != draft.occurrence_relationship
                or draft._runtime_proposal is None
                or request.candidate is not draft._runtime_proposal.candidate
            ):
                raise ValueError("score request changed after planning")
            referenced_requests.add(request.request_id)
        elif draft.score_request is not None:
            raise ValueError("non-scoring draft retained a score request")
        if draft.collision_peer_proposal_ids and draft.action != ACTION_REJECT:
            raise ValueError("collision draft was not rejected")
    if referenced_requests != set(request_by_id):
        raise ValueError("batch has unreferenced score requests")


def _validate_parent_row_binding(categories, draft):
    try:
        locator = json.loads(draft.parent_locator_json)
        category_index = locator["category_index"]
        category_name = locator["category_name"]
        locus_key = locator["locus_key"]
        row_ordinal = locator["row_ordinal"]
        locator_context_key = locator["context_key"]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("parent locator changed after planning") from error
    if (
        not isinstance(category_index, int)
        or isinstance(category_index, bool)
        or not 0 <= category_index < len(_CATEGORY_NAMES)
        or category_name != _CATEGORY_NAMES[category_index]
        or locator_context_key != json.loads(draft.context_key_json)
        or not isinstance(row_ordinal, int)
        or isinstance(row_ordinal, bool)
        or row_ordinal < 0
    ):
        raise ValueError("parent locator changed after planning")

    locator_locus_json = _json_key(locus_key, "parent locator locus key")
    matching_rows = [
        rows
        for observed_locus_key, rows in categories[category_index].items()
        if _json_key(
            observed_locus_key,
            "current parent locus key",
        ) == locator_locus_json
    ]
    if len(matching_rows) != 1 or row_ordinal >= len(matching_rows[0]):
        raise ValueError("parent locator changed after planning")
    if matching_rows[0][row_ordinal] is not draft._runtime_parent_row:
        raise ValueError("runtime parent row changed after planning")


def _revalidate_runtime_state(context, categories, batch):
    _validate_source_context(context)
    if context.source_sha256 != batch.source_sha256:
        raise ValueError("source changed after regular-subpath planning")
    if context.source_coordinate_offset != 0:
        raise ValueError("source context offset changed after planning")
    _validate_snapshot(categories, batch._category_snapshot)
    _validate_source_cds(context.source_sequence, batch.source_cds)
    if tuple(
        (cds.start, cds.end, cds.is_complement, cds.confidence)
        for cds in batch.source_cds
    ) != batch._source_cds_snapshot:
        raise ValueError("source CDS changed after planning")
    observed_accepted = tuple(
        _source_geometry(
            context.source_sequence,
            candidate,
            "accepted candidate after planning",
        )
        for candidate in batch._runtime_accepted_candidates
    )
    if observed_accepted != batch._runtime_accepted_geometries:
        raise ValueError("an accepted candidate changed after planning")
    if tuple(
        geometry.geometry_sha256 for geometry in observed_accepted
    ) != batch.accepted_geometry_sha256s:
        raise ValueError("accepted geometry digests changed after planning")

    for draft in batch.drafts:
        _validate_parent_row_binding(categories, draft)
        parent = draft._runtime_parent_row[1]
        observed_parent = _source_geometry(
            context.source_sequence, parent, "diagnostic parent after planning"
        )
        if observed_parent != draft.parent_geometry:
            raise ValueError("a diagnostic parent changed after planning")
        observed_family = tuple(
            _source_geometry(
                context.source_sequence,
                candidate,
                "family candidate after planning",
            )
            for candidate in draft._runtime_family_candidates
        )
        if observed_family != draft._runtime_family_geometries:
            raise ValueError("a complete family changed after planning")
        parent_decision = validate_accepted_family_tandem(
            context.source_sequence,
            parent,
            draft._runtime_family_candidates,
            batch.source_cds,
        )
        parent_audit = _tandem_audit(
            "parent",
            observed_parent,
            draft._runtime_supplied_parent_tandem_decision,
            parent_decision,
        )
        if parent_audit != draft.parent_tandem_audit:
            raise ValueError("parent tandem evidence changed after planning")
        if draft._runtime_proposal is not None:
            proposal = draft._runtime_proposal
            observed_proposal = _source_geometry(
                context.source_sequence,
                proposal.candidate,
                "proposal after planning",
            )
            if observed_proposal != draft.proposal_geometry:
                raise ValueError("a proposal changed after planning")
            proposal_decision = validate_accepted_family_tandem(
                context.source_sequence,
                proposal.candidate,
                draft._runtime_family_candidates + (proposal.candidate,),
                batch.source_cds,
            )
            proposal_audit = _tandem_audit(
                "proposal",
                observed_proposal,
                draft._runtime_supplied_proposal_tandem_decision,
                proposal_decision,
            )
            if proposal_audit != draft.proposal_tandem_audit:
                raise ValueError("proposal tandem evidence changed after planning")

    locations, _accepted, _diagnostics, _snapshot = _index_categories(
        categories
    )
    row_geometry_by_id = {
        id(location.row): _source_geometry(
            context.source_sequence,
            location.row[1],
            "pretransaction category row after planning",
        )
        for location in locations
    }
    regenerated_snapshot = _build_report_replay_snapshot(
        context.source_sequence,
        context.source_sha256,
        batch.source_cds,
        locations,
        observed_accepted,
        batch.drafts,
        batch.score_requests,
        row_geometry_by_id=row_geometry_by_id,
    )
    if _json_key(
        regenerated_snapshot,
        "regenerated report replay snapshot",
    ) != batch.report_replay_snapshot_json:
        raise ValueError("report replay snapshot changed after planning")


def _validate_fresh_result(context, request, result):
    if not isinstance(result, RegularSubpathFreshScoreResult):
        raise TypeError("fresh_scorer must return RegularSubpathFreshScoreResult")
    if result.request_id != request.request_id:
        raise ValueError("fresh score result refers to the wrong request")
    if not isinstance(result.evidence_method, str) or not result.evidence_method:
        raise ValueError("fresh score evidence_method changed after creation")
    row = result.scored_row
    _validate_row(row, accepted=True, label="fresh scored row")
    if row[1] is not request.candidate:
        raise ValueError("fresh scored row must retain the exact request candidate")
    geometry = _source_geometry(
        context.source_sequence, row[1], "fresh scored proposal"
    )
    if geometry != request.proposal_geometry:
        raise ValueError("fresh scored proposal geometry changed")
    return result


def _prior_basis(row):
    if len(row) != 4:
        return None
    value = row[3].get("classification_basis")
    return None if value is None else _json_copy(value, "prior basis")


def _admitted_row(
    draft,
    result,
    accepted_geometry_sha256s,
    report_replay_snapshot,
    boundary_assessment,
):
    scored = result.scored_row
    metadata = {} if len(scored) == 3 else _json_copy(
        scored[3], "fresh scorer metadata"
    )
    prior_fresh_basis = metadata.get("classification_basis")
    if prior_fresh_basis is not None:
        prior_fresh_basis = _json_copy(
            prior_fresh_basis,
            "fresh scorer prior classification_basis",
        )
    prior_boundary_assessment = metadata.get("boundary_assessment")
    if prior_boundary_assessment is not None:
        prior_boundary_assessment = _json_copy(
            prior_boundary_assessment,
            "fresh scorer prior boundary_assessment",
        )
    metadata["classification_basis"] = {
        "method": CLASSIFICATION_METHOD,
        "rule": {"name": RULE_NAME, "version": RULE_VERSION},
        "source_sha256": draft.score_request.source_sha256,
        "classification_cap": POSSIBLE_CATEGORY_NAME,
        "bona_fide_admission_allowed": False,
        "exact_boundary_certified": False,
        "boundary_assessment_id": boundary_assessment.assessment_id,
        "report_replay_snapshot": _json_copy(
            report_replay_snapshot,
            "classification report replay snapshot",
        ),
        "parent": {
            "locator": json.loads(draft.parent_locator_json),
            "score": draft.parent_score,
            "prior_classification_basis": _prior_basis(
                draft._runtime_parent_row
            ),
            "preserved_in_diagnostic_category": True,
            "geometry": draft.parent_geometry.as_dict(),
        },
        "proposal": {
            "proposal_id": draft.proposal_id,
            "geometry": draft.proposal_geometry.as_dict(),
            "occurrence_relationship": (
                draft.occurrence_relationship.as_dict()
            ),
            "boundary_status": (
                "policy_selected_occurrence_path_repeat_phase_unassessed"
            ),
            "planner_boundary_status": BOUNDARY_AMBIGUOUS_PARENT_PHASE,
            "planner_boundary_status_interpretation": (
                "conservative_planner_limitation_not_normalized_phase_evidence"
            ),
            "exact_biological_repeat_phase": False,
            "exact_array_boundaries": False,
            "normalized_boundary_dimensions": {
                "occurrence_path": "policy_selected",
                "repeat_phase": "unassessed",
            },
            "structural_planner": _structural_planner_record(
                draft.planner_status,
                draft.planner_reason_codes,
                draft._runtime_proposal,
            ),
        },
        "tandem_validation": {
            "parent": draft.parent_tandem_audit.as_dict(),
            "proposal": draft.proposal_tandem_audit.as_dict(),
            "complete_family": _complete_family_record(
                draft._runtime_family_geometries
            ),
        },
        "accepted_overlap_protection": {
            "accepted_rows_identity_complete": True,
            "accepted_geometry_sha256s": list(
                accepted_geometry_sha256s
            ),
        },
        "fresh_classifier": {
            "request_id": result.request_id,
            "score": float(scored[0]),
            "evidence_method": result.evidence_method,
            "prior_classification_basis": prior_fresh_basis,
            "prior_boundary_assessment": prior_boundary_assessment,
        },
    }
    metadata["boundary_assessment"] = boundary_assessment.as_dict()
    _json_copy(metadata, "admitted row metadata")
    return [float(scored[0]), scored[1], scored[2], metadata]


def _count_rows(categories, indexes=None):
    selected = range(5) if indexes is None else indexes
    return sum(
        len(rows)
        for index in selected
        for rows in categories[index].values()
    )


def apply_source_regular_subpath_recovery_batch(
    context,
    categories,
    batch,
    fresh_scorer,
):
    """Fresh-score eligible proposals and append them to Possible atomically."""

    _validate_source_context(context)
    _validate_batch_structure(batch)
    if not callable(fresh_scorer):
        raise TypeError("fresh_scorer must be callable")
    _revalidate_runtime_state(context, categories, batch)

    results = {}
    scorer_calls = 0
    for request in batch.score_requests:
        result = fresh_scorer(request)
        scorer_calls += 1
        result = _validate_fresh_result(context, request, result)
        if request.request_id in results:
            raise RuntimeError("fresh scorer produced a duplicate request result")
        results[request.request_id] = result

    # A scorer can retain and mutate any object passed to it.  Replay all state
    # after the final callback and before constructing even the copied output.
    _validate_batch_structure(batch)
    _revalidate_runtime_state(context, categories, batch)
    for request in batch.score_requests:
        _validate_fresh_result(context, request, results[request.request_id])

    report_replay_snapshot = _validate_report_replay_snapshot_json(
        batch.report_replay_snapshot_json,
        batch.source_sha256,
    )
    boundary_normalizer = BoundaryAssessmentNormalizer(
        context.source_sequence,
        batch.source_sha256,
    )

    additions = []
    records = []
    for draft in batch.drafts:
        if draft.action == ACTION_SCORING_REQUIRED:
            result = results[draft.score_request.request_id]
            boundary_assessment = _regular_subpath_boundary_assessment(
                boundary_normalizer,
                draft,
            )
            admitted = _admitted_row(
                draft,
                result,
                batch.accepted_geometry_sha256s,
                report_replay_snapshot,
                boundary_assessment,
            )
            additions.append((draft, admitted))
            action = ACTION_ADMIT_POSSIBLE
            reasons = ("source_regular_subpath_admitted_possible",)
            fresh_score = float(result.scored_row[0])
            evidence_method = result.evidence_method
            request_id = result.request_id
            boundary_assessment_id = boundary_assessment.assessment_id
        else:
            action = draft.action
            reasons = draft.reason_codes
            fresh_score = None
            evidence_method = None
            request_id = None
            boundary_assessment_id = None
        records.append(RegularSubpathRecoveryRecord(
            action=action,
            reason_codes=tuple(reasons),
            source_sha256=batch.source_sha256,
            context_key_json=draft.context_key_json,
            parent_locator_json=draft.parent_locator_json,
            target_locus_key_json=draft.target_locus_key_json,
            proposal_id=draft.proposal_id,
            request_id=request_id,
            parent_score=draft.parent_score,
            fresh_score=fresh_score,
            fresh_evidence_method=evidence_method,
            parent_geometry=draft.parent_geometry,
            proposal_geometry=draft.proposal_geometry,
            occurrence_relationship=draft.occurrence_relationship,
            parent_tandem_audit=draft.parent_tandem_audit,
            proposal_tandem_audit=draft.proposal_tandem_audit,
            complete_family_geometries=draft._runtime_family_geometries,
            accepted_geometry_sha256s=batch.accepted_geometry_sha256s,
            planner_status=draft.planner_status,
            planner_reason_codes=draft.planner_reason_codes,
            structural_proposal_json=(
                None
                if draft._runtime_proposal is None
                else _json_key(
                    draft._runtime_proposal.as_dict(),
                    "structural proposal",
                )
            ),
            report_replay_snapshot_sha256=(
                report_replay_snapshot["snapshot_sha256"]
            ),
            boundary_assessment_id=boundary_assessment_id,
            collision_peer_proposal_ids=draft.collision_peer_proposal_ids,
        ))

    for record in records:
        _json_copy(record.as_dict(), "regular-subpath recovery record")

    # Every replay, callback, score, and metadata check above succeeds before
    # the copy exists.  The input category graph remains untouched.
    new_categories = [
        {locus_key: list(rows) for locus_key, rows in category.items()}
        for category in categories
    ]
    for draft, admitted in additions:
        rows = new_categories[POSSIBLE_CATEGORY_INDEX].setdefault(
            draft._runtime_target_locus_key, []
        )
        rows.append(admitted)

    total_before = _count_rows(categories)
    total_after = _count_rows(new_categories)
    accepted_before = _count_rows(categories, _ACCEPTED_CATEGORY_INDEXES)
    accepted_after = _count_rows(new_categories, _ACCEPTED_CATEGORY_INDEXES)
    possible_before = _count_rows(categories, (POSSIBLE_CATEGORY_INDEX,))
    possible_after = _count_rows(new_categories, (POSSIBLE_CATEGORY_INDEX,))
    diagnostic_before = _count_rows(categories, _DIAGNOSTIC_CATEGORY_INDEXES)
    diagnostic_after = _count_rows(new_categories, _DIAGNOSTIC_CATEGORY_INDEXES)
    admission_count = len(additions)
    if total_after != total_before + admission_count:
        raise RuntimeError("regular-subpath transaction lost or duplicated rows")
    if accepted_after != accepted_before + admission_count:
        raise RuntimeError("accepted row conservation failed")
    if possible_after != possible_before + admission_count:
        raise RuntimeError("Possible row conservation failed")
    if diagnostic_after != diagnostic_before:
        raise RuntimeError("diagnostic parents were not preserved")
    if _count_rows(new_categories, (0,)) != _count_rows(categories, (0,)):
        raise RuntimeError("regular-subpath recovery changed Bona-fide")
    _validate_snapshot(categories, batch._category_snapshot)

    all_original_row_ids = {
        id(row)
        for category in categories
        for rows in category.values()
        for row in rows
    }
    retained_row_ids = {
        id(row)
        for category in new_categories
        for rows in category.values()
        for row in rows
    }
    if not all_original_row_ids.issubset(retained_row_ids):
        raise RuntimeError("an original category row was not preserved")

    reason_sets = [set(draft.reason_codes) for draft in batch.drafts]
    stats = SourceRegularSubpathRecoveryStats(
        total_rows_before=total_before,
        total_rows_after=total_after,
        accepted_rows_before=accepted_before,
        accepted_rows_after=accepted_after,
        possible_rows_before=possible_before,
        possible_rows_after=possible_after,
        diagnostic_rows_before=diagnostic_before,
        diagnostic_rows_after=diagnostic_after,
        parents_reviewed=batch.parent_context_count,
        structural_proposals=sum(
            draft.proposal_id is not None for draft in batch.drafts
        ),
        collision_rejections=sum(
            "proposal_interval_collision" in reasons for reasons in reason_sets
        ),
        tandem_missing_rejections=sum(
            any(code.endswith("tandem_evidence_missing") for code in reasons)
            for reasons in reason_sets
        ),
        tandem_stale_rejections=sum(
            any("stale_or_mismatched" in code for code in reasons)
            or "proposal_tandem_evidence_set_mismatched" in reasons
            for reasons in reason_sets
        ),
        tandem_supported_rejections=sum(
            any(code.endswith("coding_tandem_supported") for code in reasons)
            for reasons in reason_sets
        ),
        scoring_requests=len(batch.score_requests),
        scorer_calls=scorer_calls,
        possible_admissions=admission_count,
        bona_fide_admissions=0,
    )
    _json_copy(stats.as_dict(), "regular-subpath recovery stats")
    return new_categories, tuple(records), stats


__all__ = (
    "ACTION_ABSTAIN",
    "ACTION_ADMIT_POSSIBLE",
    "ACTION_REJECT",
    "ACTION_SCORING_REQUIRED",
    "CLASSIFICATION_METHOD",
    "REPORT_REPLAY_SNAPSHOT_SCHEMA",
    "REPORT_REPLAY_SNAPSHOT_VERSION",
    "ROW_VALUE_NORMALIZER_VERSION",
    "RegularSubpathFreshScoreResult",
    "RegularSubpathGeometry",
    "RegularSubpathOccurrencePathRelationship",
    "RegularSubpathProposalTandemEvidence",
    "RegularSubpathRecoveryDraft",
    "RegularSubpathRecoveryRecord",
    "RegularSubpathRecoveryScoreRequest",
    "RegularSubpathTandemAudit",
    "SourceRegularSubpathDiagnosticRowContext",
    "SourceRegularSubpathRecoveryBatch",
    "SourceRegularSubpathRecoveryStats",
    "apply_source_regular_subpath_recovery_batch",
    "canonicalize_regular_subpath_row_value",
    "plan_source_regular_subpath_recovery_batch",
    "regular_subpath_row_value_sha256",
)
