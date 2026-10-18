"""Atomic accepted-row integration for source-periodic lattice repairs.

Only Bona-fide and Possible rows are reviewed.  A geometry-only screen runs
before coding work, source CDS records are supplied by the caller, every
proposal alternative is rebuilt through the public replay API, and every
known geometry is checked against the complete deduplicated family.  Safe
representatives are exposed together for one caller-controlled fresh scoring
pass.

The commit is copy-on-write and fail-closed.  Collision components abstain
without score arbitration or iteration-order choices.  A successful repair
always remains under the ``Possible`` ceiling and is stored under the actual
representative source envelope.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
from numbers import Real
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
    DIMENSION_REPEAT_PHASE,
    STATUS_AMBIGUOUS,
    STATUS_POLICY_SELECTED,
    geometry_input_from_candidate,
)
from components.reconstruction_validator import (
    ReconstructionValidationError,
    validate_candidate_reconstruction,
)
from components.source_periodic_lattice_repair import (
    SOURCE_PERIODIC_LATTICE_RULE,
    SourcePeriodicLatticeGeometryPrefilterResult,
    SourcePeriodicLatticeProposalResult,
    SourcePeriodicLatticeRule,
    SourcePeriodicLatticeSafetyEvidence,
    propose_source_periodic_lattice_repair,
    replay_source_periodic_lattice_proposal,
    source_periodic_lattice_geometry_prefilter,
    source_sequence_sha256,
)


RULE_NAME = "source-periodic-lattice-category-transaction"
RULE_VERSION = "1.0.0"
CLASSIFICATION_METHOD = "source_periodic_lattice_repair"
COORDINATE_SYSTEM = "zero_based_half_open"
CLASSIFICATION_CEILING = "Possible"

CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
ACCEPTED_CATEGORY_INDEXES = (0, 2)
POSSIBLE_CATEGORY_INDEX = 2

ACTION_PREFILTER_REJECTED = "geometry_prefilter_rejected"
ACTION_PARENT_TANDEM_REJECTED = "parent_tandem_rejected"
ACTION_PROPOSAL_ABSTAINED = "proposal_abstained"
ACTION_FAMILY_TANDEM_REJECTED = "complete_family_tandem_rejected"
ACTION_COLLISION_ABSTAINED = "collision_abstained"
ACTION_SCORING_REQUIRED = "scoring_required"
ACTION_REPAIRED = "repaired_possible"


@dataclass(frozen=True)
class SourcePeriodicLatticeFamilyTandemReview:
    """One direct decision bound to the complete replayed family."""

    role: str
    alternative_index: Optional[int]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    decision: AcceptedFamilyTandemDecision
    _runtime_candidate: object

    def as_dict(self):
        return {
            "role": self.role,
            "alternative_index": self.alternative_index,
            "repeat_intervals": [
                list(interval) for interval in self.repeat_intervals
            ],
            "decision": _json_copy(
                self.decision.as_dict(),
                "complete-family tandem decision",
            ),
        }


@dataclass(frozen=True)
class SourcePeriodicLatticeScoreRequest:
    """One exact representative that requires a fresh complete score."""

    request_id: str
    candidate: object
    source_category_index: int
    source_category_name: str
    source_locus_key_json: str
    source_row_ordinal: int
    target_locus_key_json: str
    source_sha256: str
    proposal: SourcePeriodicLatticeProposalResult

    @property
    def source_locus_key(self):
        return json.loads(self.source_locus_key_json)

    @property
    def target_locus_key(self):
        return json.loads(self.target_locus_key_json)

    def as_dict(self):
        return {
            "request_id": self.request_id,
            "source_category_index": self.source_category_index,
            "source_category_name": self.source_category_name,
            "source_locus_key": self.source_locus_key,
            "source_row_ordinal": self.source_row_ordinal,
            "target_locus_key": self.target_locus_key,
            "source_sha256": self.source_sha256,
            "representative_repeat_intervals": [
                list(interval)
                for interval in self.proposal.representative.repeat_intervals
            ],
        }


@dataclass(frozen=True)
class SourcePeriodicLatticeScoreResult:
    """Fresh normal-ensemble row for one exact request."""

    request_id: str
    scored_row: object


@dataclass(frozen=True)
class SourcePeriodicLatticeRepairDraft:
    """One immutable accepted-row planning decision."""

    action: str
    reason_codes: Tuple[str, ...]
    source_category_index: int
    source_category_name: str
    source_locus_key_json: str
    source_row_ordinal: int
    parent_row: object
    geometry_prefilter: SourcePeriodicLatticeGeometryPrefilterResult
    parent_direct_tandem: Optional[AcceptedFamilyTandemDecision]
    proposal: Optional[SourcePeriodicLatticeProposalResult]
    replayed_candidates: Tuple[object, ...]
    complete_family_tandem_reviews: Tuple[
        SourcePeriodicLatticeFamilyTandemReview, ...
    ]
    score_request: Optional[SourcePeriodicLatticeScoreRequest]
    collision_with: Tuple[str, ...]
    _runtime_source_locus_key: object
    _runtime_target_locus_key: object

    @property
    def source_locus_key(self):
        return json.loads(self.source_locus_key_json)

    def as_dict(self):
        return {
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "source_category_index": self.source_category_index,
            "source_category_name": self.source_category_name,
            "source_locus_key": self.source_locus_key,
            "source_row_ordinal": self.source_row_ordinal,
            "geometry_prefilter": _json_copy(
                asdict(self.geometry_prefilter),
                "geometry prefilter",
            ),
            "parent_direct_tandem": (
                None
                if self.parent_direct_tandem is None
                else _json_copy(
                    self.parent_direct_tandem.as_dict(),
                    "parent direct tandem",
                )
            ),
            "proposal": (
                None
                if self.proposal is None
                else _proposal_as_dict(self.proposal)
            ),
            "complete_family_tandem_reviews": [
                review.as_dict()
                for review in self.complete_family_tandem_reviews
            ],
            "score_request": (
                None
                if self.score_request is None
                else self.score_request.as_dict()
            ),
            "collision_with": list(self.collision_with),
        }


@dataclass(frozen=True)
class SourcePeriodicLatticeCostAudit:
    accepted_rows_reviewed: int
    geometry_prefilter_calls: int
    geometry_prefilter_eligible: int
    parent_direct_tandem_evaluations: int
    proposal_searches: int
    proposal_replays: int
    complete_family_tandem_evaluations: int
    scoring_requests: int
    source_cds_prediction_calls: int
    source_cds_reused_from_pipeline: bool
    estimated_planning_full_source_hashes: int
    cost_note: str

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourcePeriodicLatticeRepairBatch:
    """Frozen bridge between planning and one fresh scoring pass."""

    source_sha256: str
    drafts: Tuple[SourcePeriodicLatticeRepairDraft, ...]
    score_requests: Tuple[SourcePeriodicLatticeScoreRequest, ...]
    cost_audit: SourcePeriodicLatticeCostAudit
    _category_snapshot: tuple

    @property
    def proposed_candidates(self):
        return tuple(request.candidate for request in self.score_requests)

    def as_dict(self):
        return {
            "source_sha256": self.source_sha256,
            "drafts": [draft.as_dict() for draft in self.drafts],
            "score_requests": [
                request.as_dict() for request in self.score_requests
            ],
            "cost_audit": self.cost_audit.as_dict(),
        }


@dataclass(frozen=True)
class SourcePeriodicLatticeRepairRecord:
    """One detached final decision with complete proposal provenance."""

    action: str
    reason_codes: Tuple[str, ...]
    source_category_index: int
    source_category_name: str
    source_locus_key_json: str
    source_row_ordinal: int
    target_category_index: Optional[int]
    target_category_name: Optional[str]
    target_locus_key_json: Optional[str]
    parent_score: float
    fresh_score: Optional[float]
    parent_repeat_intervals: Tuple[Tuple[int, int], ...]
    representative_repeat_intervals: Tuple[Tuple[int, int], ...]
    geometry_prefilter: dict
    parent_direct_tandem: Optional[dict]
    proposal: Optional[dict]
    complete_family_tandem_reviews: Tuple[dict, ...]
    fresh_scoring_metadata: object
    boundary_assessment_id: Optional[str]
    collision_with: Tuple[str, ...]

    @property
    def source_locus_key(self):
        return json.loads(self.source_locus_key_json)

    @property
    def target_locus_key(self):
        if self.target_locus_key_json is None:
            return None
        return json.loads(self.target_locus_key_json)

    def as_dict(self):
        return {
            "rule": {"name": RULE_NAME, "version": RULE_VERSION},
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "coordinate_system": COORDINATE_SYSTEM,
            "classification_ceiling": CLASSIFICATION_CEILING,
            "source_category_index": self.source_category_index,
            "source_category_name": self.source_category_name,
            "source_locus_key": self.source_locus_key,
            "source_row_ordinal": self.source_row_ordinal,
            "target_category_index": self.target_category_index,
            "target_category_name": self.target_category_name,
            "target_locus_key": self.target_locus_key,
            "parent_score": self.parent_score,
            "fresh_score": self.fresh_score,
            "parent_repeat_intervals": [
                list(interval) for interval in self.parent_repeat_intervals
            ],
            "representative_repeat_intervals": [
                list(interval)
                for interval in self.representative_repeat_intervals
            ],
            "geometry_prefilter": _json_copy(
                self.geometry_prefilter,
                "record geometry prefilter",
            ),
            "parent_direct_tandem": _json_copy(
                self.parent_direct_tandem,
                "record parent direct tandem",
            ),
            "proposal": _json_copy(self.proposal, "record proposal"),
            "complete_family_tandem_reviews": _json_copy(
                list(self.complete_family_tandem_reviews),
                "record complete-family reviews",
            ),
            "fresh_scoring_metadata": _json_copy(
                self.fresh_scoring_metadata,
                "record fresh scoring metadata",
            ),
            "boundary_assessment_id": self.boundary_assessment_id,
            "boundary_status": "AMBIGUOUS",
            "exact_boundary_certified": False,
            "collision_with": list(self.collision_with),
        }


@dataclass(frozen=True)
class SourcePeriodicLatticeRepairStats:
    total_rows_before: int
    total_rows_after: int
    accepted_rows_before: int
    accepted_rows_after: int
    accepted_rows_reviewed: int
    geometry_prefilter_rejected: int
    parent_tandem_rejected: int
    proposal_abstained: int
    complete_family_tandem_rejected: int
    collision_abstained: int
    scoring_requests: int
    repaired_rows: int
    bona_fide_to_possible: int
    possible_to_possible: int
    apply_proposal_replays: int
    apply_direct_tandem_revalidations: int

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class _AcceptedLocation:
    category_index: int
    locus_key: object
    locus_key_json: str
    row_ordinal: int
    row: object
    repeat_intervals: tuple
    envelope: tuple


def _json_copy(value, label):
    try:
        return json.loads(json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ))
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


def _proposal_as_dict(proposal):
    return _json_copy(asdict(proposal), "periodic-lattice proposal")


def _candidate_geometry(source, candidate):
    result = validate_candidate_reconstruction(source, candidate)
    if not result.valid:
        raise ReconstructionValidationError(result)
    intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts,
            candidate.list_repeats,
        )
    )
    return intervals, (intervals[0][0], intervals[-1][1])


def _validate_row(source, row):
    if not isinstance(row, (tuple, list)) or len(row) not in (3, 4):
        raise ValueError(
            "category rows require score, candidate, features, and "
            "optional metadata"
        )
    if (
        isinstance(row[0], bool)
        or not isinstance(row[0], Real)
        or not math.isfinite(float(row[0]))
    ):
        raise ValueError("category row score must be finite")
    if not isinstance(row[2], (tuple, list)):
        raise ValueError("category row features must be a sequence")
    if len(row) == 4 and not isinstance(row[3], Mapping):
        raise ValueError("category row metadata must be a mapping")
    return _candidate_geometry(source, row[1])


def _validate_source_cds(source, source_cds):
    if not isinstance(source_cds, tuple):
        raise TypeError("source_cds must be a tuple of SourceCDS records")
    for index, cds in enumerate(source_cds):
        if not isinstance(cds, SourceCDS):
            raise TypeError(
                "source_cds[{}] is not a SourceCDS".format(index)
            )
        if cds.end > len(source):
            raise ValueError(
                "source_cds[{}] lies outside the source".format(index)
            )


def _index_categories(source, categories):
    if not isinstance(categories, (tuple, list)) or len(categories) != 5:
        raise ValueError("periodic-lattice repair requires five categories")
    row_ids = set()
    candidate_ids = set()
    accepted = []
    snapshot = []
    for category_index, category in enumerate(categories):
        if not isinstance(category, Mapping):
            raise ValueError("every category must be a mapping")
        category_snapshot = []
        for locus_key, rows in category.items():
            if not isinstance(rows, (tuple, list)) or not rows:
                raise ValueError(
                    "category locus values must be non-empty row sequences"
                )
            locus_json = _json_key(locus_key, "category locus key")
            for row_ordinal, row in enumerate(rows):
                intervals, envelope = _validate_row(source, row)
                if id(row) in row_ids:
                    raise ValueError("one row occurs more than once")
                if id(row[1]) in candidate_ids:
                    raise ValueError(
                        "one candidate has more than one category row"
                    )
                row_ids.add(id(row))
                candidate_ids.add(id(row[1]))
                category_snapshot.append((
                    locus_json,
                    row_ordinal,
                    row,
                    row[1],
                    intervals,
                    float(row[0]),
                ))
                if category_index in ACCEPTED_CATEGORY_INDEXES:
                    accepted.append(_AcceptedLocation(
                        category_index=category_index,
                        locus_key=locus_key,
                        locus_key_json=locus_json,
                        row_ordinal=row_ordinal,
                        row=row,
                        repeat_intervals=intervals,
                        envelope=envelope,
                    ))
        snapshot.append(tuple(category_snapshot))
    return tuple(accepted), tuple(snapshot)


def _snapshot_matches(source, categories, expected):
    _accepted, observed = _index_categories(source, categories)
    if len(observed) != len(expected):
        return False
    for observed_category, expected_category in zip(observed, expected):
        if len(observed_category) != len(expected_category):
            return False
        for current, prior in zip(observed_category, expected_category):
            if (
                current[:2] != prior[:2]
                or current[2] is not prior[2]
                or current[3] is not prior[3]
                or current[4:] != prior[4:]
            ):
                return False
    return True


def _intervals_overlap(first, second):
    return max(first[0], second[0]) < min(first[1], second[1])


def _deduplicated_family(parent, alternatives, source):
    family = []
    seen = set()
    for candidate in (parent,) + tuple(alternatives):
        geometry, _envelope = _candidate_geometry(source, candidate)
        if geometry in seen:
            continue
        seen.add(geometry)
        family.append(candidate)
    return tuple(family)


def _complete_family_reviews(source, parent, alternatives, source_cds):
    complete_family = _deduplicated_family(
        parent,
        alternatives,
        source,
    )
    entries = (("parent", None, parent),) + tuple(
        ("alternative", index, candidate)
        for index, candidate in enumerate(alternatives)
    )
    reviews = []
    for role, alternative_index, candidate in entries:
        intervals, _envelope = _candidate_geometry(source, candidate)
        decision = validate_accepted_family_tandem(
            source,
            candidate,
            complete_family,
            source_cds,
        )
        reviews.append(SourcePeriodicLatticeFamilyTandemReview(
            role=role,
            alternative_index=alternative_index,
            repeat_intervals=intervals,
            decision=decision,
            _runtime_candidate=candidate,
        ))
    return tuple(reviews)


def _request_id(source_hash, location, proposal):
    payload = {
        "rule": [RULE_NAME, RULE_VERSION],
        "source_sha256": source_hash,
        "source_category_index": location.category_index,
        "source_locus_key": json.loads(location.locus_key_json),
        "source_row_ordinal": location.row_ordinal,
        "representative_repeat_intervals": [
            list(interval)
            for interval in proposal.representative.repeat_intervals
        ],
    }
    return "source-periodic-lattice-request-sha256:" + hashlib.sha256(
        json.dumps(
            payload,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("ascii")
    ).hexdigest()


def _collision_label(location):
    return "accepted:{}:{}:{}".format(
        location.category_index,
        location.locus_key_json,
        location.row_ordinal,
    )


def _initial_draft(location, screen, **values):
    defaults = {
        "action": ACTION_PREFILTER_REJECTED,
        "reason_codes": (screen.reason,),
        "parent_direct_tandem": None,
        "proposal": None,
        "replayed_candidates": (),
        "complete_family_tandem_reviews": (),
        "score_request": None,
        "collision_with": (),
        "_runtime_target_locus_key": None,
    }
    defaults.update(values)
    return SourcePeriodicLatticeRepairDraft(
        source_category_index=location.category_index,
        source_category_name=CATEGORY_NAMES[location.category_index],
        source_locus_key_json=location.locus_key_json,
        source_row_ordinal=location.row_ordinal,
        parent_row=location.row,
        geometry_prefilter=screen,
        _runtime_source_locus_key=location.locus_key,
        **defaults
    )


def _proposal_envelopes(draft):
    if draft.proposal is None:
        return ()
    return tuple(
        alternative.envelope for alternative in draft.proposal.alternatives
    )


def _apply_collision_abstentions(drafts, accepted_locations):
    collision_reasons = {index: set() for index in range(len(drafts))}
    collision_with = {index: set() for index in range(len(drafts))}
    planned = tuple(
        index for index, draft in enumerate(drafts)
        if draft.action == ACTION_SCORING_REQUIRED
    )

    for index in planned:
        draft = drafts[index]
        for other in accepted_locations:
            if other.row is draft.parent_row:
                continue
            if any(
                _intervals_overlap(envelope, other.envelope)
                for envelope in _proposal_envelopes(draft)
            ):
                collision_reasons[index].add(
                    "proposal_overlaps_other_accepted_row"
                )
                collision_with[index].add(_collision_label(other))

    for left_position, left_index in enumerate(planned):
        for right_index in planned[left_position + 1:]:
            if any(
                _intervals_overlap(left, right)
                for left in _proposal_envelopes(drafts[left_index])
                for right in _proposal_envelopes(drafts[right_index])
            ):
                left_id = drafts[left_index].score_request.request_id
                right_id = drafts[right_index].score_request.request_id
                collision_reasons[left_index].add(
                    "proposal_overlaps_other_proposal"
                )
                collision_reasons[right_index].add(
                    "proposal_overlaps_other_proposal"
                )
                collision_with[left_index].add(right_id)
                collision_with[right_index].add(left_id)

    resolved = []
    for index, draft in enumerate(drafts):
        if not collision_reasons[index]:
            resolved.append(draft)
            continue
        resolved.append(replace(
            draft,
            action=ACTION_COLLISION_ABSTAINED,
            reason_codes=tuple(sorted(collision_reasons[index])),
            score_request=None,
            collision_with=tuple(sorted(collision_with[index])),
            _runtime_target_locus_key=None,
        ))
    return tuple(resolved)


def plan_source_periodic_lattice_repairs(
    source,
    categories,
    source_cds,
    rule=SOURCE_PERIODIC_LATTICE_RULE,
):
    """Review every accepted row and expose only fully authorized requests."""

    if not isinstance(source, str) or not source:
        raise ValueError("source must be a non-empty DNA string")
    if not isinstance(rule, SourcePeriodicLatticeRule):
        raise TypeError("rule must be SourcePeriodicLatticeRule")
    _validate_source_cds(source, source_cds)
    accepted_locations, snapshot = _index_categories(source, categories)
    source_hash = source_sequence_sha256(source)
    drafts = []
    cost = {
        "screen": 0,
        "eligible": 0,
        "parent_tandem": 0,
        "search": 0,
        "replay": 0,
        "family_tandem": 0,
    }

    for location in accepted_locations:
        parent = location.row[1]
        screen = source_periodic_lattice_geometry_prefilter(
            source,
            parent,
            rule=rule,
        )
        cost["screen"] += 1
        if not screen.eligible:
            drafts.append(_initial_draft(location, screen))
            continue
        cost["eligible"] += 1

        parent_direct = validate_accepted_family_tandem(
            source,
            parent,
            (parent,),
            source_cds,
        )
        cost["parent_tandem"] += 1
        if parent_direct.supported:
            drafts.append(_initial_draft(
                location,
                screen,
                action=ACTION_PARENT_TANDEM_REJECTED,
                reason_codes=("parent_direct_tandem_supported",),
                parent_direct_tandem=parent_direct,
            ))
            continue

        evidence = SourcePeriodicLatticeSafetyEvidence(
            source_sha256=source_hash,
            coding_context_assessed=True,
            tandem_evidence_assessed=True,
            coding_tandem_supported=False,
            tandem_decision_reason=parent_direct.reason,
        )
        proposal = propose_source_periodic_lattice_repair(
            source,
            parent,
            evidence,
            rule=rule,
        )
        cost["search"] += 1
        if proposal.status != "proposed":
            drafts.append(_initial_draft(
                location,
                screen,
                action=ACTION_PROPOSAL_ABSTAINED,
                reason_codes=(proposal.reason,),
                parent_direct_tandem=parent_direct,
                proposal=proposal,
            ))
            continue
        if (
            proposal.classification_ceiling != CLASSIFICATION_CEILING
            or proposal.representative is None
            or not proposal.repeat_count_invariant
            or not proposal.repeat_start_locations_invariant
        ):
            raise RuntimeError(
                "proposal violated the transaction acceptance contract"
            )

        replayed = replay_source_periodic_lattice_proposal(source, proposal)
        cost["replay"] += 1
        if len(replayed) != len(proposal.alternatives):
            raise RuntimeError("proposal replay changed alternative cardinality")
        reviews = _complete_family_reviews(
            source,
            parent,
            replayed,
            source_cds,
        )
        cost["family_tandem"] += len(reviews)
        supported = tuple(
            review for review in reviews if review.decision.supported
        )
        if supported:
            drafts.append(_initial_draft(
                location,
                screen,
                action=ACTION_FAMILY_TANDEM_REJECTED,
                reason_codes=tuple(sorted({
                    "complete_family_tandem_supported_{}".format(
                        review.role
                    )
                    for review in supported
                })),
                parent_direct_tandem=parent_direct,
                proposal=proposal,
                replayed_candidates=replayed,
                complete_family_tandem_reviews=reviews,
            ))
            continue

        representative = replayed[proposal.representative_index]
        representative_geometry, target_key = _candidate_geometry(
            source,
            representative,
        )
        if representative_geometry != proposal.representative.repeat_intervals:
            raise RuntimeError("replayed representative geometry is stale")
        target_key_json = _json_key(target_key, "target geometry key")
        request = SourcePeriodicLatticeScoreRequest(
            request_id=_request_id(source_hash, location, proposal),
            candidate=representative,
            source_category_index=location.category_index,
            source_category_name=CATEGORY_NAMES[location.category_index],
            source_locus_key_json=location.locus_key_json,
            source_row_ordinal=location.row_ordinal,
            target_locus_key_json=target_key_json,
            source_sha256=source_hash,
            proposal=proposal,
        )
        drafts.append(_initial_draft(
            location,
            screen,
            action=ACTION_SCORING_REQUIRED,
            reason_codes=("source_periodic_lattice_repair_scoring_required",),
            parent_direct_tandem=parent_direct,
            proposal=proposal,
            replayed_candidates=replayed,
            complete_family_tandem_reviews=reviews,
            score_request=request,
            _runtime_target_locus_key=target_key,
        ))

    drafts = _apply_collision_abstentions(
        tuple(drafts),
        accepted_locations,
    )
    requests = tuple(
        draft.score_request for draft in drafts
        if draft.action == ACTION_SCORING_REQUIRED
    )
    if len({request.request_id for request in requests}) != len(requests):
        raise RuntimeError("two lattice repairs have one request id")
    cost_audit = SourcePeriodicLatticeCostAudit(
        accepted_rows_reviewed=len(accepted_locations),
        geometry_prefilter_calls=cost["screen"],
        geometry_prefilter_eligible=cost["eligible"],
        parent_direct_tandem_evaluations=cost["parent_tandem"],
        proposal_searches=cost["search"],
        proposal_replays=cost["replay"],
        complete_family_tandem_evaluations=cost["family_tandem"],
        scoring_requests=len(requests),
        source_cds_prediction_calls=0,
        source_cds_reused_from_pipeline=True,
        estimated_planning_full_source_hashes=(
            1
            + cost["screen"]
            + cost["parent_tandem"]
            + 2 * cost["search"]
            + 3 * cost["replay"]
            + cost["family_tandem"]
        ),
        cost_note=(
            "source CDS prediction is caller-owned and never repeated here; "
            "the estimate counts explicit full-source hashes in geometry "
            "screens, tandem decisions, proposal authorization, and replay"
        ),
    )
    batch = SourcePeriodicLatticeRepairBatch(
        source_sha256=source_hash,
        drafts=drafts,
        score_requests=requests,
        cost_audit=cost_audit,
        _category_snapshot=snapshot,
    )
    _json_copy(batch.as_dict(), "periodic-lattice repair batch")
    return batch


def _index_score_results(batch, score_results):
    if isinstance(score_results, (str, bytes)) or score_results is None:
        raise TypeError("score_results must be an iterable")
    try:
        score_results = tuple(score_results)
    except TypeError as error:
        raise TypeError("score_results must be an iterable") from error
    expected = {request.request_id: request for request in batch.score_requests}
    indexed = {}
    for result in score_results:
        if not isinstance(result, SourcePeriodicLatticeScoreResult):
            raise TypeError(
                "score results must be SourcePeriodicLatticeScoreResult"
            )
        if result.request_id in indexed:
            raise ValueError("duplicate score result request id")
        if result.request_id not in expected:
            raise ValueError("score result does not belong to the batch")
        indexed[result.request_id] = result
    if set(indexed) != set(expected):
        raise ValueError("score result mapping is incomplete")
    return indexed


def _validate_scored_row(source, request, row):
    _validate_row(source, row)
    if row[1] is not request.candidate:
        raise ValueError(
            "fresh scorer changed representative identity or order"
        )
    intervals, envelope = _candidate_geometry(source, row[1])
    if intervals != request.proposal.representative.repeat_intervals:
        raise ValueError("fresh scorer changed representative geometry")
    if list(envelope) != request.target_locus_key:
        raise ValueError("score request target key is not actual geometry")
    if len(row) == 4:
        _json_copy(row[3], "fresh scoring metadata")


def _same_reviews(first, second):
    return tuple(review.as_dict() for review in first) == tuple(
        review.as_dict() for review in second
    )


def _boundary_assessment(
    normalizer,
    parent,
    representative,
    proposal,
    replayed,
):
    source_hash = normalizer.source_sha256
    relationships = []
    representative_alternative = proposal.representative
    parent_dimensions = []
    if (
        len(proposal.parent_repeat_intervals)
        != len(representative_alternative.repeat_intervals)
        or representative_alternative.terminal_added_slots
        or representative_alternative.internal_missing_slots
    ):
        parent_dimensions.append(DIMENSION_OCCURRENCE_PATH)
    if any(
        mapping.parent_interval != mapping.proposed_interval
        for mapping in representative_alternative.parent_mappings
    ):
        parent_dimensions.append(DIMENSION_REPEAT_PHASE)
    if parent_dimensions:
        relationships.append(AlternativeRelationshipInput(
            geometry=geometry_input_from_candidate(parent, source_hash),
            label="accepted-parent-path",
            dimensions=tuple(parent_dimensions),
            reason_codes=(
                "accepted_parent_source_replayed",
                "periodic_lattice_policy_selected_representative",
            ),
        ))

    for index, candidate in enumerate(replayed):
        if index == proposal.representative_index:
            continue
        relationships.append(AlternativeRelationshipInput(
            geometry=geometry_input_from_candidate(candidate, source_hash),
            label="periodic-lattice-phase-alternative-{}".format(index),
            dimensions=(DIMENSION_REPEAT_PHASE,),
            reason_codes=(
                "source_valid_plausible_repeat_phase",
                "repeat_starts_and_count_invariant",
            ),
        ))

    occurrence = (
        BoundaryDimensionInput(
            status=STATUS_POLICY_SELECTED,
            selection_policy_id=(
                "source-periodic-lattice-phase-then-terminal:"
                + proposal.rule_version
            ),
            reason_codes=(
                "occurrence_count_source_lattice_supported",
                "reported_path_selected_by_periodic_lattice_policy",
            ),
        )
        if any(
            DIMENSION_OCCURRENCE_PATH in relationship.dimensions
            for relationship in relationships
        )
        else BoundaryDimensionInput()
    )
    repeat_phase = (
        BoundaryDimensionInput(
            status=STATUS_AMBIGUOUS,
            reason_codes=(
                "plausible_source_valid_repeat_phases_retained",
                "reported_geometry_is_representative_not_exact",
            ),
        )
        if any(
            DIMENSION_REPEAT_PHASE in relationship.dimensions
            for relationship in relationships
        )
        else BoundaryDimensionInput()
    )
    return normalizer.normalize(
        geometry_input_from_candidate(representative, source_hash),
        alternative_relationships=tuple(relationships),
        occurrence_path=occurrence,
        repeat_phase=repeat_phase,
    )


def _record(
    draft,
    *,
    action=None,
    reason_codes=None,
    fresh_score=None,
    fresh_scoring_metadata=None,
    boundary_assessment_id=None,
):
    proposal = draft.proposal
    representative = (
        ()
        if proposal is None or proposal.representative is None
        else proposal.representative.repeat_intervals
    )
    request = draft.score_request
    return SourcePeriodicLatticeRepairRecord(
        action=draft.action if action is None else action,
        reason_codes=(
            draft.reason_codes if reason_codes is None else reason_codes
        ),
        source_category_index=draft.source_category_index,
        source_category_name=draft.source_category_name,
        source_locus_key_json=draft.source_locus_key_json,
        source_row_ordinal=draft.source_row_ordinal,
        target_category_index=(
            POSSIBLE_CATEGORY_INDEX if action == ACTION_REPAIRED else None
        ),
        target_category_name=(
            CATEGORY_NAMES[POSSIBLE_CATEGORY_INDEX]
            if action == ACTION_REPAIRED else None
        ),
        target_locus_key_json=(
            request.target_locus_key_json
            if action == ACTION_REPAIRED else None
        ),
        parent_score=float(draft.parent_row[0]),
        fresh_score=fresh_score,
        parent_repeat_intervals=draft.geometry_prefilter.parent_repeat_intervals,
        representative_repeat_intervals=representative,
        geometry_prefilter=_json_copy(
            asdict(draft.geometry_prefilter),
            "record geometry prefilter",
        ),
        parent_direct_tandem=(
            None
            if draft.parent_direct_tandem is None
            else _json_copy(
                draft.parent_direct_tandem.as_dict(),
                "record parent direct tandem",
            )
        ),
        proposal=(
            None if proposal is None else _proposal_as_dict(proposal)
        ),
        complete_family_tandem_reviews=tuple(
            review.as_dict()
            for review in draft.complete_family_tandem_reviews
        ),
        fresh_scoring_metadata=fresh_scoring_metadata,
        boundary_assessment_id=boundary_assessment_id,
        collision_with=draft.collision_with,
    )


def _replacement_row(parent_row, scored_row, record, assessment):
    original_metadata = (
        dict(parent_row[3]) if len(parent_row) == 4 else {}
    )
    prior_basis = original_metadata.get("classification_basis")
    if prior_basis is not None:
        prior_basis = _json_copy(
            prior_basis,
            "prior classification_basis",
        )
    prior_boundary = original_metadata.get("boundary_assessment")
    if prior_boundary is not None:
        prior_boundary = _json_copy(
            prior_boundary,
            "prior boundary_assessment",
        )
    basis = {
        "method": CLASSIFICATION_METHOD,
        "prior_classification_basis": prior_basis,
        "prior_boundary_assessment": prior_boundary,
        "source_periodic_lattice_repair": _json_copy(
            record.as_dict(),
            "periodic-lattice classification record",
        ),
        "transition": {
            "source_category_index": record.source_category_index,
            "source_category_name": record.source_category_name,
            "target_category_index": POSSIBLE_CATEGORY_INDEX,
            "target_category_name": CATEGORY_NAMES[POSSIBLE_CATEGORY_INDEX],
        },
        "classification_ceiling": CLASSIFICATION_CEILING,
        "boundary_status": "AMBIGUOUS",
        "exact_boundary_certified": False,
    }
    original_metadata["classification_basis"] = basis
    original_metadata["boundary_assessment"] = assessment.as_dict()
    fields = (
        scored_row[0],
        scored_row[1],
        scored_row[2],
        original_metadata,
    )
    return tuple(fields) if isinstance(parent_row, tuple) else list(fields)


def _count_rows(categories, indexes=None):
    selected = range(len(categories)) if indexes is None else indexes
    return sum(
        len(rows)
        for category_index in selected
        for rows in categories[category_index].values()
    )


def apply_source_periodic_lattice_repairs(
    source,
    categories,
    source_cds,
    batch,
    score_results,
):
    """Reauthorize a complete batch, then replace rows in one copied state."""

    if not isinstance(batch, SourcePeriodicLatticeRepairBatch):
        raise TypeError("batch must be SourcePeriodicLatticeRepairBatch")
    _validate_source_cds(source, source_cds)
    if source_sequence_sha256(source) != batch.source_sha256:
        raise ValueError("source changed after periodic-lattice planning")
    if not _snapshot_matches(source, categories, batch._category_snapshot):
        raise ValueError("categories changed after periodic-lattice planning")
    results = _index_score_results(batch, score_results)
    for request in batch.score_requests:
        _validate_scored_row(
            source,
            request,
            results[request.request_id].scored_row,
        )

    normalizer = BoundaryAssessmentNormalizer(
        source,
        batch.source_sha256,
    )
    replacements = {}
    records = []
    apply_replays = 0
    apply_tandem = 0
    for draft in batch.drafts:
        if draft.action != ACTION_SCORING_REQUIRED:
            records.append(_record(draft))
            continue
        request = draft.score_request
        proposal = draft.proposal
        parent = draft.parent_row[1]

        parent_direct = validate_accepted_family_tandem(
            source,
            parent,
            (parent,),
            source_cds,
        )
        apply_tandem += 1
        if (
            parent_direct != draft.parent_direct_tandem
            or parent_direct.supported
        ):
            raise RuntimeError(
                "parent tandem authorization changed before commit"
            )
        replayed = replay_source_periodic_lattice_proposal(source, proposal)
        apply_replays += 1
        fresh_reviews = _complete_family_reviews(
            source,
            parent,
            replayed,
            source_cds,
        )
        apply_tandem += len(fresh_reviews)
        if (
            not _same_reviews(
                fresh_reviews,
                draft.complete_family_tandem_reviews,
            )
            or any(review.decision.supported for review in fresh_reviews)
        ):
            raise RuntimeError(
                "complete-family tandem authorization changed before commit"
            )
        representative = request.candidate
        assessment = _boundary_assessment(
            normalizer,
            parent,
            representative,
            proposal,
            replayed,
        )
        scored_row = results[request.request_id].scored_row
        fresh_metadata = (
            None if len(scored_row) == 3 else _json_copy(
                scored_row[3],
                "fresh scoring metadata",
            )
        )
        record = _record(
            draft,
            action=ACTION_REPAIRED,
            reason_codes=("source_periodic_lattice_repaired_possible",),
            fresh_score=float(scored_row[0]),
            fresh_scoring_metadata=fresh_metadata,
            boundary_assessment_id=assessment.assessment_id,
        )
        replacement_row = _replacement_row(
            draft.parent_row,
            scored_row,
            record,
            assessment,
        )
        replacements[id(draft.parent_row)] = (
            draft,
            replacement_row,
        )
        records.append(record)

    new_categories = [
        {locus_key: list(rows) for locus_key, rows in category.items()}
        for category in categories
    ]
    for draft, replacement_row in replacements.values():
        rows = new_categories[draft.source_category_index].get(
            draft._runtime_source_locus_key
        )
        if rows is None:
            raise RuntimeError("parent locus disappeared during commit")
        matches = [
            index for index, row in enumerate(rows)
            if row is draft.parent_row
        ]
        if len(matches) != 1:
            raise RuntimeError("parent row identity changed during commit")
        rows.pop(matches[0])
        if not rows:
            del new_categories[draft.source_category_index][
                draft._runtime_source_locus_key
            ]
        new_categories[POSSIBLE_CATEGORY_INDEX].setdefault(
            draft._runtime_target_locus_key,
            [],
        ).append(replacement_row)

    total_before = _count_rows(categories)
    total_after = _count_rows(new_categories)
    accepted_before = _count_rows(categories, ACCEPTED_CATEGORY_INDEXES)
    accepted_after = _count_rows(new_categories, ACCEPTED_CATEGORY_INDEXES)
    if total_before != total_after:
        raise RuntimeError("lattice repair changed total row cardinality")
    if accepted_before != accepted_after:
        raise RuntimeError("lattice repair changed accepted cardinality")
    parent_ids = set(replacements)
    observed_parent_ids = {
        id(row)
        for category in new_categories
        for rows in category.values()
        for row in rows
        if id(row) in parent_ids
    }
    if observed_parent_ids:
        raise RuntimeError("a repaired parent remained beside its child")
    child_candidate_ids = [
        id(replacement[1][1]) for replacement in replacements.values()
    ]
    if len(child_candidate_ids) != len(set(child_candidate_ids)):
        raise RuntimeError("one repaired child was inserted more than once")
    if not _snapshot_matches(source, categories, batch._category_snapshot):
        raise RuntimeError("input categories mutated during transaction")

    repaired = len(replacements)
    stats = SourcePeriodicLatticeRepairStats(
        total_rows_before=total_before,
        total_rows_after=total_after,
        accepted_rows_before=accepted_before,
        accepted_rows_after=accepted_after,
        accepted_rows_reviewed=len(batch.drafts),
        geometry_prefilter_rejected=sum(
            draft.action == ACTION_PREFILTER_REJECTED
            for draft in batch.drafts
        ),
        parent_tandem_rejected=sum(
            draft.action == ACTION_PARENT_TANDEM_REJECTED
            for draft in batch.drafts
        ),
        proposal_abstained=sum(
            draft.action == ACTION_PROPOSAL_ABSTAINED
            for draft in batch.drafts
        ),
        complete_family_tandem_rejected=sum(
            draft.action == ACTION_FAMILY_TANDEM_REJECTED
            for draft in batch.drafts
        ),
        collision_abstained=sum(
            draft.action == ACTION_COLLISION_ABSTAINED
            for draft in batch.drafts
        ),
        scoring_requests=len(batch.score_requests),
        repaired_rows=repaired,
        bona_fide_to_possible=sum(
            draft.source_category_index == 0
            for draft, _replacement in replacements.values()
        ),
        possible_to_possible=sum(
            draft.source_category_index == 2
            for draft, _replacement in replacements.values()
        ),
        apply_proposal_replays=apply_replays,
        apply_direct_tandem_revalidations=apply_tandem,
    )
    _json_copy(
        [record.as_dict() for record in records],
        "periodic-lattice repair records",
    )
    _json_copy(stats.as_dict(), "periodic-lattice repair stats")
    return new_categories, tuple(records), stats


__all__ = (
    "ACTION_COLLISION_ABSTAINED",
    "ACTION_FAMILY_TANDEM_REJECTED",
    "ACTION_PARENT_TANDEM_REJECTED",
    "ACTION_PREFILTER_REJECTED",
    "ACTION_PROPOSAL_ABSTAINED",
    "ACTION_REPAIRED",
    "ACTION_SCORING_REQUIRED",
    "CLASSIFICATION_METHOD",
    "SourcePeriodicLatticeCostAudit",
    "SourcePeriodicLatticeFamilyTandemReview",
    "SourcePeriodicLatticeRepairBatch",
    "SourcePeriodicLatticeRepairDraft",
    "SourcePeriodicLatticeRepairRecord",
    "SourcePeriodicLatticeRepairStats",
    "SourcePeriodicLatticeScoreRequest",
    "SourcePeriodicLatticeScoreResult",
    "apply_source_periodic_lattice_repairs",
    "plan_source_periodic_lattice_repairs",
)
