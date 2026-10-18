"""Two-phase category transaction for accepted-core boundary completion.

The geometry planner in :mod:`accepted_core_boundary_completion` deliberately
separates array existence from terminal-boundary evidence.  This module keeps
that separation intact at the category boundary:

1. every currently accepted row receives one explicit strict-core and complete
   witness context;
2. source-backed geometry is planned without consulting a fresh model score;
3. all changed geometries are exposed together for one caller-controlled batch
   scoring pass; and
4. an exact, complete score-result mapping is validated before any replacement
   is made in a copy of the categories.

The fresh score and feature vectors are reported for the completed geometry,
but they never choose its category or decide whether the already accepted
array exists.  A successful completion stays in the exact Bona-fide/Possible
category and locus occupied by its authorizing core.

Coordinates are zero-based and end-exclusive.  Rows and candidates are
selected by object identity, never by value equality.
"""

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Real
from typing import Optional, Tuple

from components.accepted_core_boundary_completion import (
    BoundaryGeometryWitness,
    CoreAuthorizationEvidence,
    FinalBoundaryEvidence,
    SourceGeometry,
    WitnessContextEvidence,
    plan_accepted_core_boundary_completion,
)
from components.components_evaluated_arrays_enhancement import (
    source_backed_terminal_repeat_proposals,
    source_backed_weak_terminal_repeat_proposals,
)


RULE_NAME = "accepted-core-boundary-completion"
RULE_VERSION = "1.0.0"
CLASSIFICATION_METHOD = "source_backed_accepted_core_boundary_completion"
COORDINATE_SYSTEM = "zero_based_half_open"

_CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
_ACCEPTED_CATEGORY_INDEXES = (0, 2)
_GEOMETRY_READY_REASON = ("final_evidence_missing",)


@dataclass(frozen=True)
class AcceptedCoreBoundaryContext:
    """Complete caller contract for one identity-exact accepted row.

    ``witnesses`` must contain explicit :class:`BoundaryGeometryWitness`
    objects so that the accepted status of every supplied path is auditable.
    ``witness_context`` is an explicit completeness attestation; both of its
    fields must be true.  The transaction additionally cross-checks every
    supplied ``accepted=True`` witness against the accepted category rows.
    """

    accepted_row: object
    core_decision: object
    core_authorization: CoreAuthorizationEvidence
    witnesses: Tuple[BoundaryGeometryWitness, ...]
    witness_context: WitnessContextEvidence

    def __post_init__(self):
        if not isinstance(self.core_authorization, CoreAuthorizationEvidence):
            raise TypeError(
                "core_authorization must be CoreAuthorizationEvidence"
            )
        if not isinstance(self.witnesses, tuple):
            raise TypeError("witnesses must be an explicit tuple")
        for index, witness in enumerate(self.witnesses):
            if not isinstance(witness, BoundaryGeometryWitness):
                raise TypeError(
                    "witnesses[{}] is not a BoundaryGeometryWitness".format(index)
                )
        if not isinstance(self.witness_context, WitnessContextEvidence):
            raise TypeError(
                "witness_context must be WitnessContextEvidence"
            )


@dataclass(frozen=True)
class BoundaryCompletionScoreRequest:
    """One exact proposed candidate that must participate in batch scoring."""

    request_id: str
    candidate: object
    category_index: int
    category_name: str
    locus_key_json: str
    row_ordinal: int
    source_sha256: str
    core_geometry: SourceGeometry
    proposed_geometry: SourceGeometry

    @property
    def locus_key(self):
        return json.loads(self.locus_key_json)

    def as_dict(self):
        return {
            "request_id": self.request_id,
            "category_index": self.category_index,
            "category_name": self.category_name,
            "locus_key": self.locus_key,
            "row_ordinal": self.row_ordinal,
            "source_sha256": self.source_sha256,
            "core_geometry": _geometry_as_dict(self.core_geometry),
            "proposed_geometry": _geometry_as_dict(self.proposed_geometry),
        }


@dataclass(frozen=True)
class BoundaryCompletionScoreResult:
    """Fresh row and final evidence for one exact scoring request."""

    request_id: str
    scored_row: object
    final_evidence: FinalBoundaryEvidence

    def __post_init__(self):
        if not isinstance(self.request_id, str) or not self.request_id:
            raise ValueError("a score result needs a non-empty request_id")
        if not isinstance(self.final_evidence, FinalBoundaryEvidence):
            raise TypeError("final_evidence must be FinalBoundaryEvidence")


@dataclass(frozen=True)
class AcceptedCoreBoundaryDraft:
    """One immutable planning decision for an accepted category row."""

    action: str
    reason_codes: Tuple[str, ...]
    category_index: int
    category_name: str
    locus_key_json: str
    row_ordinal: int
    core_row: object
    planner_plan: object
    score_request: Optional[BoundaryCompletionScoreRequest]
    _runtime_locus_key: object
    _witness_geometry_snapshot: tuple

    @property
    def locus_key(self):
        return json.loads(self.locus_key_json)

    def as_dict(self):
        value = {
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "category_index": self.category_index,
            "category_name": self.category_name,
            "locus_key": self.locus_key,
            "row_ordinal": self.row_ordinal,
            "core_score": float(self.planner_plan.core_score),
            "core_geometry": _geometry_as_dict(
                self.planner_plan.core_geometry
            ),
            "proposed_geometry": _geometry_as_dict(
                self.planner_plan.proposed_geometry
            ),
        }
        if self.score_request is not None:
            value["score_request_id"] = self.score_request.request_id
        return value


@dataclass(frozen=True)
class AcceptedCoreBoundaryCompletionBatch:
    """Frozen bridge between source planning and one batch scoring pass."""

    source_sha256: str
    drafts: Tuple[AcceptedCoreBoundaryDraft, ...]
    score_requests: Tuple[BoundaryCompletionScoreRequest, ...]
    _category_snapshot: tuple

    @property
    def proposed_candidates(self):
        """Return exact candidates in the required batch-scoring order."""

        return tuple(request.candidate for request in self.score_requests)

    def as_dict(self):
        return {
            "source_sha256": self.source_sha256,
            "accepted_rows_planned": len(self.drafts),
            "score_requests": [
                request.as_dict() for request in self.score_requests
            ],
            "drafts": [draft.as_dict() for draft in self.drafts],
        }


@dataclass(frozen=True)
class AcceptedCoreBoundaryCompletionRecord:
    """One runtime-object-free final decision for an accepted core."""

    action: str
    reason_codes: Tuple[str, ...]
    category_index: int
    category_name: str
    locus_key_json: str
    row_ordinal: int
    request_id: Optional[str]
    core_score: float
    fresh_score: Optional[float]
    core_geometry: SourceGeometry
    proposed_geometry: SourceGeometry
    final_evidence: Optional[FinalBoundaryEvidence]

    @property
    def locus_key(self):
        return json.loads(self.locus_key_json)

    def as_dict(self):
        return {
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "category_index": self.category_index,
            "category_name": self.category_name,
            "locus_key": self.locus_key,
            "row_ordinal": self.row_ordinal,
            "request_id": self.request_id,
            "core_score": self.core_score,
            "fresh_score": self.fresh_score,
            "core_geometry": _geometry_as_dict(self.core_geometry),
            "proposed_geometry": _geometry_as_dict(self.proposed_geometry),
            "final_evidence": (
                None
                if self.final_evidence is None
                else _final_evidence_as_dict(self.final_evidence)
            ),
        }


@dataclass(frozen=True)
class AcceptedCoreBoundaryCompletionStats:
    """Cardinality and decision totals for the copy-on-write transaction."""

    total_rows_before: int
    total_rows_after: int
    accepted_rows_before: int
    accepted_rows_after: int
    accepted_rows_planned: int
    scoring_requests: int
    completed_rows: int
    unchanged_rows: int
    planning_abstentions: int
    final_evidence_abstentions: int
    bona_fide_completed: int
    possible_completed: int

    def as_dict(self):
        return {
            field_name: getattr(self, field_name)
            for field_name in self.__dataclass_fields__
        }


@dataclass(frozen=True)
class _RowLocation:
    category_index: int
    locus_key: object
    locus_key_json: Optional[str]
    row_ordinal: int
    row: object


@dataclass(frozen=True)
class _SnapshotRow:
    row: object
    row_length: int
    score: float
    candidate: object
    features: object
    metadata: object


def _source_hash(source):
    return hashlib.sha256(source.upper().encode("ascii")).hexdigest()


def _validate_source(source):
    if not isinstance(source, str):
        raise TypeError("source must be a string")
    if not source:
        raise ValueError("source must be non-empty")
    try:
        source.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("source must contain ASCII DNA symbols") from error


def _finite_score(value, unit_interval=False):
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(float(value))
    ):
        return False
    return not unit_interval or 0.0 <= float(value) <= 1.0


def _json_copy(value, label):
    try:
        serialized = json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        return json.loads(serialized)
    except (TypeError, ValueError) as error:
        raise ValueError("{} must be JSON-safe".format(label)) from error


def _locus_key_json(locus_key):
    try:
        return json.dumps(
            locus_key,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise ValueError("accepted-row locus keys must be JSON-safe") from error


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
        prior = row[3].get("classification_basis")
        if prior is not None:
            if not isinstance(prior, Mapping):
                raise ValueError(
                    "{} prior classification_basis must be a mapping".format(
                        label
                    )
                )
            _json_copy(prior, "{} prior classification_basis".format(label))


def _snapshot_row(row):
    return _SnapshotRow(
        row=row,
        row_length=len(row),
        score=float(row[0]),
        candidate=row[1],
        features=row[2],
        metadata=(row[3] if len(row) == 4 else None),
    )


def _index_categories(categories):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("boundary completion requires exactly five categories")

    all_locations = []
    accepted_locations = []
    row_ids = set()
    candidate_ids = set()
    snapshot = []
    for category_index, category in enumerate(categories):
        if not isinstance(category, Mapping):
            raise ValueError("each category must be a mapping")
        snapshot_entries = []
        for locus_key, rows in category.items():
            if not isinstance(rows, (list, tuple)) or not rows:
                raise ValueError(
                    "category locus values must be non-empty row sequences"
                )
            locus_json = (
                _locus_key_json(locus_key)
                if category_index in _ACCEPTED_CATEGORY_INDEXES
                else None
            )
            snapshot_rows = []
            for row_ordinal, row in enumerate(rows):
                accepted = category_index in _ACCEPTED_CATEGORY_INDEXES
                _validate_row(row, accepted=accepted)
                if id(row) in row_ids:
                    raise ValueError(
                        "one scored row occurs more than once in categories"
                    )
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
                    row_ordinal=row_ordinal,
                    row=row,
                )
                all_locations.append(location)
                if accepted:
                    accepted_locations.append(location)
                snapshot_rows.append(_snapshot_row(row))
            snapshot_entries.append((locus_key, tuple(snapshot_rows)))
        snapshot.append(tuple(snapshot_entries))
    return (
        tuple(all_locations),
        tuple(accepted_locations),
        tuple(snapshot),
    )


def _validate_snapshot(categories, snapshot):
    if not isinstance(categories, (list, tuple)) or len(categories) != len(snapshot):
        raise ValueError("categories changed after boundary planning")
    for category_index, expected_entries in enumerate(snapshot):
        category = categories[category_index]
        if not isinstance(category, Mapping) or len(category) != len(expected_entries):
            raise ValueError("categories changed after boundary planning")
        for locus_key, expected_rows in expected_entries:
            if locus_key not in category:
                raise ValueError("a category locus changed after boundary planning")
            observed_rows = category[locus_key]
            if (
                not isinstance(observed_rows, (list, tuple))
                or len(observed_rows) != len(expected_rows)
            ):
                raise ValueError("category row cardinality changed after planning")
            for observed, expected in zip(observed_rows, expected_rows):
                if (
                    observed is not expected.row
                    or len(observed) != expected.row_length
                    or float(observed[0]) != expected.score
                    or observed[1] is not expected.candidate
                    or observed[2] is not expected.features
                    or (
                        observed[3] if len(observed) == 4 else None
                    ) is not expected.metadata
                ):
                    raise ValueError("a category row changed after boundary planning")


def _index_contexts(contexts, accepted_locations):
    if contexts is None or isinstance(contexts, (str, bytes)):
        raise TypeError("accepted_core_contexts must be an iterable")
    try:
        contexts = tuple(contexts)
    except TypeError as error:
        raise TypeError("accepted_core_contexts must be an iterable") from error

    accepted_by_row_id = {
        id(location.row): location for location in accepted_locations
    }
    accepted_by_candidate_id = {
        id(location.row[1]): location for location in accepted_locations
    }
    indexed = {}
    for index, context in enumerate(contexts):
        if not isinstance(context, AcceptedCoreBoundaryContext):
            raise TypeError(
                "accepted_core_contexts[{}] is not an "
                "AcceptedCoreBoundaryContext".format(index)
            )
        row_id = id(context.accepted_row)
        if row_id in indexed:
            raise ValueError("an accepted row has duplicate boundary contexts")
        location = accepted_by_row_id.get(row_id)
        if location is None or location.row is not context.accepted_row:
            raise ValueError(
                "boundary context refers to a non-accepted or foreign row"
            )
        if not context.witness_context.complete:
            raise ValueError(
                "every accepted row requires a complete witness context"
            )

        witness_candidate_ids = tuple(
            id(witness.candidate) for witness in context.witnesses
        )
        witness_path_ids = tuple(
            witness.path_id for witness in context.witnesses
        )
        if len(set(witness_candidate_ids)) != len(witness_candidate_ids):
            raise ValueError("boundary witnesses must not contain duplicates")
        if len(set(witness_path_ids)) != len(witness_path_ids):
            raise ValueError("boundary witness path IDs must be unique")
        if id(context.accepted_row[1]) in witness_candidate_ids:
            raise ValueError("an accepted core must not witness itself")
        for witness in context.witnesses:
            actual = accepted_by_candidate_id.get(id(witness.candidate))
            actually_accepted = (
                actual is not None and actual.row[1] is witness.candidate
            )
            if witness.accepted != actually_accepted:
                raise ValueError(
                    "witness accepted status disagrees with categories"
                )
        indexed[row_id] = context

    missing_count = sum(
        id(location.row) not in indexed for location in accepted_locations
    )
    if missing_count:
        raise ValueError(
            "{} accepted row(s) are missing boundary context".format(
                missing_count
            )
        )
    if len(indexed) != len(accepted_locations):
        raise RuntimeError("accepted boundary context cardinality changed")
    return indexed, accepted_by_candidate_id


def _intervals_overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _request_id(source_sha256, location, plan):
    payload = {
        "source_sha256": source_sha256,
        "category_index": location.category_index,
        "locus_key": json.loads(location.locus_key_json),
        "row_ordinal": location.row_ordinal,
        "core_path": [list(value) for value in plan.core_geometry.repeat_intervals],
        "proposed_path": [
            list(value) for value in plan.proposed_geometry.repeat_intervals
        ],
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return "accepted-core-boundary-" + hashlib.sha256(encoded).hexdigest()[:24]


def _geometry_as_dict(geometry):
    return {
        "coordinate_system": COORDINATE_SYSTEM,
        "source_interval": list(geometry.source_interval),
        "repeat_count": geometry.repeat_count,
        "repeat_intervals": [
            list(interval) for interval in geometry.repeat_intervals
        ],
        "repeat_sequences": list(geometry.repeat_sequences),
        "spacer_sequences": list(geometry.spacer_sequences),
    }


def _step_as_dict(step):
    return {
        "tier": step.tier,
        "side": step.side,
        "repeat_interval": list(step.repeat_interval),
        "spacer_interval": list(step.spacer_interval),
        "witness_path_ids": list(step.witness_path_ids),
        "evidence": [
            {"name": name, "value": value} for name, value in step.evidence
        ],
    }


def _final_evidence_as_dict(evidence):
    return {
        "fresh_score": float(evidence.fresh_score),
        "source_reconstruction_valid": evidence.source_reconstruction_valid,
        "structural_valid": evidence.structural_valid,
        "phase_unambiguous": evidence.phase_unambiguous,
        "coding_context_assessed": evidence.coding_context_assessed,
        "tandem_evidence_assessed": evidence.tandem_evidence_assessed,
        "coding_tandem_clear": evidence.coding_tandem_clear,
        "reason_codes": list(evidence.reason_codes),
    }


def _final_failure_codes(evidence):
    failed = []
    for field_name in (
        "source_reconstruction_valid",
        "structural_valid",
        "phase_unambiguous",
        "coding_context_assessed",
        "tandem_evidence_assessed",
        "coding_tandem_clear",
    ):
        if not getattr(evidence, field_name):
            failed.append("final_{}_failed".format(field_name))
    failed.extend(evidence.reason_codes)
    return tuple(failed)


def plan_accepted_core_boundary_completion_batch(
    source,
    categories,
    accepted_core_contexts,
    strong_proposal_provider=source_backed_terminal_repeat_proposals,
    weak_proposal_provider=source_backed_weak_terminal_repeat_proposals,
):
    """Plan every accepted core and return exact candidates for batch scoring.

    This function never calls a scorer.  The underlying planner is invoked
    without a final-evidence provider, so its model-score field cannot
    influence geometry planning.  A changed plan whose only outstanding gate
    is ``final_evidence_missing`` becomes one score request.
    """

    _validate_source(source)
    if not callable(strong_proposal_provider):
        raise TypeError("strong_proposal_provider must be callable")
    if not callable(weak_proposal_provider):
        raise TypeError("weak_proposal_provider must be callable")

    _all_locations, accepted_locations, snapshot = _index_categories(categories)
    contexts, accepted_by_candidate_id = _index_contexts(
        accepted_core_contexts, accepted_locations
    )
    source_sha256 = _source_hash(source)
    accepted_geometries = {
        candidate_id: SourceGeometry.from_candidate(source, location.row[1])
        for candidate_id, location in accepted_by_candidate_id.items()
    }

    drafts = []
    requests = []
    proposed_paths = {}
    proposed_candidate_ids = set()
    for location in accepted_locations:
        context = contexts[id(location.row)]
        plan = plan_accepted_core_boundary_completion(
            source,
            location.row,
            _CATEGORY_NAMES[location.category_index],
            context.core_decision,
            context.core_authorization,
            witnesses=context.witnesses,
            witness_context=context.witness_context,
            final_evidence_provider=None,
            strong_proposal_provider=strong_proposal_provider,
            weak_proposal_provider=weak_proposal_provider,
        )
        witness_geometry_snapshot = tuple(
            (
                witness.candidate,
                SourceGeometry.from_candidate(source, witness.candidate),
            )
            for witness in context.witnesses
        )

        score_request = None
        if (
            plan.action == "abstain"
            and plan.changed
            and plan.reason_codes == _GEOMETRY_READY_REASON
        ):
            witness_accepted_ids = {
                id(witness.candidate)
                for witness in context.witnesses
                if witness.accepted
            }
            for candidate_id, other_geometry in accepted_geometries.items():
                if candidate_id == id(location.row[1]):
                    continue
                if (
                    _intervals_overlap(
                        plan.proposed_geometry.source_interval,
                        other_geometry.source_interval,
                    )
                    and candidate_id not in witness_accepted_ids
                ):
                    raise ValueError(
                        "an overlapping accepted path is missing from witnesses"
                    )

            proposed_path = plan.proposed_geometry.repeat_intervals
            previous = proposed_paths.get(proposed_path)
            if previous is not None:
                raise ValueError(
                    "two accepted cores proposed the same completed geometry"
                )
            candidate = plan.materialize_proposed_candidate(source)
            if id(candidate) in proposed_candidate_ids:
                raise RuntimeError("proposed candidate identity was reused")
            if id(candidate) in accepted_by_candidate_id:
                raise RuntimeError("proposed candidate is already accepted")
            request = BoundaryCompletionScoreRequest(
                request_id=_request_id(source_sha256, location, plan),
                candidate=candidate,
                category_index=location.category_index,
                category_name=_CATEGORY_NAMES[location.category_index],
                locus_key_json=location.locus_key_json,
                row_ordinal=location.row_ordinal,
                source_sha256=source_sha256,
                core_geometry=plan.core_geometry,
                proposed_geometry=plan.proposed_geometry,
            )
            proposed_paths[proposed_path] = request.request_id
            proposed_candidate_ids.add(id(candidate))
            requests.append(request)
            action = "scoring_required"
            reasons = ("source_geometry_completion_planned",)
            score_request = request
        elif plan.action == "unchanged":
            action = "unchanged"
            reasons = plan.reason_codes
        elif plan.action == "abstain":
            action = "abstained"
            reasons = plan.reason_codes
        else:
            raise RuntimeError(
                "geometry-only planner returned unexpected action {!r}".format(
                    plan.action
                )
            )

        drafts.append(
            AcceptedCoreBoundaryDraft(
                action=action,
                reason_codes=tuple(reasons),
                category_index=location.category_index,
                category_name=_CATEGORY_NAMES[location.category_index],
                locus_key_json=location.locus_key_json,
                row_ordinal=location.row_ordinal,
                core_row=location.row,
                planner_plan=plan,
                score_request=score_request,
                _runtime_locus_key=location.locus_key,
                _witness_geometry_snapshot=witness_geometry_snapshot,
            )
        )

    if len({request.request_id for request in requests}) != len(requests):
        raise RuntimeError("boundary score request IDs are not unique")
    return AcceptedCoreBoundaryCompletionBatch(
        source_sha256=source_sha256,
        drafts=tuple(drafts),
        score_requests=tuple(requests),
        _category_snapshot=snapshot,
    )


def _index_results(results, batch):
    if results is None or isinstance(results, (str, bytes)):
        raise TypeError("score_results must be an iterable")
    try:
        results = tuple(results)
    except TypeError as error:
        raise TypeError("score_results must be an iterable") from error

    requests = {request.request_id: request for request in batch.score_requests}
    indexed = {}
    for index, result in enumerate(results):
        if not isinstance(result, BoundaryCompletionScoreResult):
            raise TypeError(
                "score_results[{}] is not a "
                "BoundaryCompletionScoreResult".format(index)
            )
        if result.request_id in indexed:
            raise ValueError("a boundary score request has duplicate results")
        if result.request_id not in requests:
            raise ValueError("score result refers to a foreign request")
        indexed[result.request_id] = result

    missing = tuple(
        request_id for request_id in requests if request_id not in indexed
    )
    if missing:
        raise ValueError(
            "{} boundary score request(s) are missing results".format(len(missing))
        )
    if len(indexed) != len(requests):
        raise RuntimeError("boundary score-result cardinality changed")
    return indexed


def _validate_batch_structure(batch):
    if not isinstance(batch, AcceptedCoreBoundaryCompletionBatch):
        raise TypeError(
            "batch must be AcceptedCoreBoundaryCompletionBatch"
        )
    requests = tuple(batch.score_requests)
    if len({request.request_id for request in requests}) != len(requests):
        raise ValueError("batch contains duplicate score request IDs")
    request_by_id = {request.request_id: request for request in requests}
    seen = set()
    for draft in batch.drafts:
        if not isinstance(draft, AcceptedCoreBoundaryDraft):
            raise TypeError("batch contains a malformed boundary draft")
        if draft.action == "scoring_required":
            request = draft.score_request
            if (
                request is None
                or request_by_id.get(request.request_id) is not request
                or request.request_id in seen
            ):
                raise ValueError("batch draft/request mapping is inconsistent")
            seen.add(request.request_id)
        elif draft.score_request is not None:
            raise ValueError("a non-scoring draft contains a score request")
        if not isinstance(draft._witness_geometry_snapshot, tuple):
            raise ValueError("batch contains a malformed witness snapshot")
        for witness_snapshot in draft._witness_geometry_snapshot:
            if (
                not isinstance(witness_snapshot, tuple)
                or len(witness_snapshot) != 2
                or not isinstance(witness_snapshot[1], SourceGeometry)
            ):
                raise ValueError("batch contains a malformed witness snapshot")
    if seen != set(request_by_id):
        raise ValueError("batch contains unclaimed score requests")


def _validate_scored_result(source, request, result):
    row = result.scored_row
    _validate_row(row, accepted=True, label="fresh scored row")
    if row[1] is not request.candidate:
        raise ValueError(
            "fresh scorer must return the exact proposed candidate object"
        )
    observed_geometry = SourceGeometry.from_candidate(source, row[1])
    if observed_geometry != request.proposed_geometry:
        raise ValueError("fresh scored candidate geometry changed after planning")
    rematerialized = request.proposed_geometry.materialize(source)
    if SourceGeometry.from_candidate(source, rematerialized) != observed_geometry:
        raise RuntimeError("planned source geometry cannot be rematerialized")
    if float(row[0]) != float(result.final_evidence.fresh_score):
        raise ValueError(
            "fresh scored row and final evidence scores disagree"
        )
    scorer_metadata = row[3] if len(row) == 4 else None
    if scorer_metadata is not None:
        _json_copy(scorer_metadata, "fresh scoring metadata")


def _completion_basis(draft, result):
    plan = draft.planner_plan
    core_row = draft.core_row
    core_metadata = core_row[3] if len(core_row) == 4 else {}
    prior_basis = core_metadata.get("classification_basis")
    if prior_basis is not None:
        prior_basis = _json_copy(
            prior_basis, "prior classification_basis"
        )
    fresh_metadata = (
        result.scored_row[3] if len(result.scored_row) == 4 else None
    )
    if fresh_metadata is not None:
        fresh_metadata = _json_copy(
            fresh_metadata, "fresh scoring metadata"
        )
    return {
        "method": CLASSIFICATION_METHOD,
        "geometry_derivation": "accepted_core_source_terminal_fixed_point",
        "coordinate_system": COORDINATE_SYSTEM,
        "rule": {
            "name": RULE_NAME,
            "version": RULE_VERSION,
        },
        "category_inheritance": {
            "category_index": draft.category_index,
            "category_name": draft.category_name,
            "locus_key": draft.locus_key,
            "fresh_score_controls_existence": False,
            "fresh_score_controls_category": False,
        },
        "authorized_core_score": float(core_row[0]),
        "fresh_completed_score": float(result.scored_row[0]),
        "prior_classification_basis": prior_basis,
        "fresh_scoring_metadata": fresh_metadata,
        "source_sha256": plan.source_sha256,
        "core_policy": {
            "name": plan.core_policy_name,
            "version": plan.core_policy_version,
        },
        "core_geometry": _geometry_as_dict(plan.core_geometry),
        "completed_geometry": _geometry_as_dict(plan.proposed_geometry),
        "completion_steps": [_step_as_dict(step) for step in plan.steps],
        "final_boundary_evidence": _final_evidence_as_dict(
            result.final_evidence
        ),
    }


def _completed_row(draft, result):
    core_row = draft.core_row
    metadata = dict(core_row[3]) if len(core_row) == 4 else {}
    metadata["classification_basis"] = _completion_basis(draft, result)
    fields = (
        result.scored_row[0],
        result.scored_row[1],
        result.scored_row[2],
        metadata,
    )
    return tuple(fields) if isinstance(core_row, tuple) else list(fields)


def _count_rows(categories, indexes=None):
    category_indexes = range(len(categories)) if indexes is None else indexes
    return sum(
        len(rows)
        for category_index in category_indexes
        for rows in categories[category_index].values()
    )


def apply_accepted_core_boundary_completion_batch(
    source,
    categories,
    batch,
    score_results,
):
    """Validate one complete batch and copy-on-write accepted replacements.

    All batch, source, category-snapshot, score cardinality, candidate
    identity, source reconstruction, metadata, and evidence checks complete
    before the copied categories are edited.  Failed final biological evidence
    retains that core row by identity.  A low fresh score alone never blocks a
    replacement and never changes its accepted category.
    """

    _validate_source(source)
    _validate_batch_structure(batch)
    if _source_hash(source) != batch.source_sha256:
        raise ValueError("source changed after boundary planning")
    _validate_snapshot(categories, batch._category_snapshot)
    for draft in batch.drafts:
        observed_core = SourceGeometry.from_candidate(
            source, draft.core_row[1]
        )
        if observed_core != draft.planner_plan.core_geometry:
            raise ValueError("an accepted core changed after boundary planning")
        for witness, expected_geometry in draft._witness_geometry_snapshot:
            observed_witness = SourceGeometry.from_candidate(source, witness)
            if observed_witness != expected_geometry:
                raise ValueError("a boundary witness changed after planning")
    results = _index_results(score_results, batch)

    for request in batch.score_requests:
        if request.source_sha256 != batch.source_sha256:
            raise ValueError("score request source hash disagrees with batch")
        _validate_scored_result(
            source, request, results[request.request_id]
        )

    replacements = {}
    records = []
    final_abstentions = 0
    completed_by_category = {0: 0, 2: 0}
    for draft in batch.drafts:
        plan = draft.planner_plan
        if draft.action == "scoring_required":
            request = draft.score_request
            result = results[request.request_id]
            if result.final_evidence.passes:
                replacement = _completed_row(draft, result)
                replacements[id(draft.core_row)] = replacement
                action = "completed"
                reasons = ("accepted_core_boundary_completed",)
                completed_by_category[draft.category_index] += 1
            else:
                action = "abstained"
                reasons = _final_failure_codes(result.final_evidence)
                final_abstentions += 1
            fresh_score = float(result.scored_row[0])
            final_evidence = result.final_evidence
            request_id = request.request_id
        else:
            action = draft.action
            reasons = draft.reason_codes
            fresh_score = None
            final_evidence = None
            request_id = None
        records.append(
            AcceptedCoreBoundaryCompletionRecord(
                action=action,
                reason_codes=tuple(reasons),
                category_index=draft.category_index,
                category_name=draft.category_name,
                locus_key_json=draft.locus_key_json,
                row_ordinal=draft.row_ordinal,
                request_id=request_id,
                core_score=float(plan.core_score),
                fresh_score=fresh_score,
                core_geometry=plan.core_geometry,
                proposed_geometry=plan.proposed_geometry,
                final_evidence=final_evidence,
            )
        )

    # Metadata construction above is intentionally completed before this copy.
    new_categories = [
        {locus_key: list(rows) for locus_key, rows in category.items()}
        for category in categories
    ]
    replaced_ids = set()
    for draft in batch.drafts:
        replacement = replacements.get(id(draft.core_row))
        if replacement is None:
            continue
        rows = new_categories[draft.category_index][draft._runtime_locus_key]
        if (
            draft.row_ordinal >= len(rows)
            or rows[draft.row_ordinal] is not draft.core_row
        ):
            raise RuntimeError(
                "accepted row identity changed during boundary transaction"
            )
        rows[draft.row_ordinal] = replacement
        replaced_ids.add(id(draft.core_row))

    if replaced_ids != set(replacements):
        raise RuntimeError("not every completed core row was replaced exactly once")

    total_before = _count_rows(categories)
    total_after = _count_rows(new_categories)
    accepted_before = _count_rows(categories, _ACCEPTED_CATEGORY_INDEXES)
    accepted_after = _count_rows(new_categories, _ACCEPTED_CATEGORY_INDEXES)
    if total_after != total_before:
        raise RuntimeError("boundary completion changed total row cardinality")
    if accepted_after != accepted_before:
        raise RuntimeError("boundary completion changed accepted cardinality")
    for category_index, category in enumerate(categories):
        if set(new_categories[category_index]) != set(category):
            raise RuntimeError("boundary completion changed category locus keys")

    # The caller's input remains an exact snapshot after the transaction.
    _validate_snapshot(categories, batch._category_snapshot)
    stats = AcceptedCoreBoundaryCompletionStats(
        total_rows_before=total_before,
        total_rows_after=total_after,
        accepted_rows_before=accepted_before,
        accepted_rows_after=accepted_after,
        accepted_rows_planned=len(batch.drafts),
        scoring_requests=len(batch.score_requests),
        completed_rows=len(replacements),
        unchanged_rows=sum(
            draft.action == "unchanged" for draft in batch.drafts
        ),
        planning_abstentions=sum(
            draft.action == "abstained" for draft in batch.drafts
        ),
        final_evidence_abstentions=final_abstentions,
        bona_fide_completed=completed_by_category[0],
        possible_completed=completed_by_category[2],
    )
    return new_categories, tuple(records), stats


__all__ = (
    "AcceptedCoreBoundaryCompletionBatch",
    "AcceptedCoreBoundaryCompletionRecord",
    "AcceptedCoreBoundaryCompletionStats",
    "AcceptedCoreBoundaryContext",
    "AcceptedCoreBoundaryDraft",
    "BoundaryCompletionScoreRequest",
    "BoundaryCompletionScoreResult",
    "CLASSIFICATION_METHOD",
    "apply_accepted_core_boundary_completion_batch",
    "plan_accepted_core_boundary_completion_batch",
)
