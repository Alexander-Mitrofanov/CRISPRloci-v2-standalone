"""Conservative source-intact recovery planning and category transaction.

This module is the narrow bridge between classifier-discarded candidates and
the accepted ``Possible`` category.  It deliberately does not discover a
candidate, choose a repeat phase, compare model scores, or certify exact
boundaries.  Every coordinate used here is replayed against the complete
source sequence before interval arithmetic is performed.

Planning and mutation are separate.  The planner requires explicit complete
family, lineage, augmentation, accepted-overlap, and direct tandem evidence.
It returns a frozen plan.  The apply function rejects a stale plan and edits a
copy of the five category mappings while preserving unrelated row objects.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from numbers import Real
from typing import Optional, Tuple

from components.accepted_family_tandem_validator import (
    AcceptedFamilyTandemDecision,
    validate_accepted_family_tandem,
)
from components.accepted_source_intact_conflict import (
    ACCEPTED_RETAIN,
    ACCEPTED_RETIRE,
    SOURCE_ADMIT_POSSIBLE,
    SOURCE_ATTACH_ALTERNATIVE,
    AcceptedSourceConflictDecision,
    arbitrate_accepted_source_intact_conflict,
)
from components.architecture_rescue_coding_context import SourceCDS
from components.components_detection_refinement import (
    AdvancedFuzzySearchFilter,
    AdvancedFuzzySearchFilterEvidence,
)
from components.reconstruction_validator import (
    IUPAC_DNA_ALPHABET,
    validate_candidate_reconstruction,
)
from components.source_intact_existence_rescue import (
    BOUNDARY_AMBIGUOUS,
    SOURCE_INTACT_EXISTENCE_POLICY,
    SourceIntactEvaluationContext,
    SourceIntactExistenceDecision,
    SourceIntactExistencePolicy,
    evaluate_source_intact_existence_in_context,
)
from components.source_intact_terminal_authorization import (
    SourceIntactTerminalAuthorization,
    authorize_source_intact_terminal_evidence,
)


RULE_NAME = "source-intact-recovery-transaction"
RULE_VERSION = "1.1.0"
CLASSIFICATION_METHOD = "source_backed_intact_existence_rescue"
COORDINATE_SYSTEM = "zero_based_half_open"

CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
ACCEPTED_CATEGORY_INDEXES = (0, 2)
RECOVERY_SOURCE_CATEGORY_INDEXES = (3, 4)
POSSIBLE_CATEGORY_INDEX = 2
DISCARDED_CATEGORY_INDEX = 3

ACTION_ADMIT_POSSIBLE = "admit_possible"
ACTION_REPLACE_POSSIBLE = "replace_conflicting_possible"
ACTION_ATTACH_ALTERNATIVE = "attach_ambiguous_alternative"
ACTION_ABSTAIN = "abstain"
ACTION_REJECT = "reject"


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


def _source_hash(source):
    return hashlib.sha256(source.encode("ascii")).hexdigest()


def _slice_hash(source, interval):
    return _source_hash(source[interval[0] : interval[1]])


def _validate_source(source):
    if not isinstance(source, str) or not source:
        raise ValueError("source must be one non-empty DNA string")
    if set(source).difference(IUPAC_DNA_ALPHABET):
        raise ValueError("source contains non-IUPAC DNA symbols")
    return source.upper()


def _validate_source_cds(source, source_cds):
    if not isinstance(source_cds, tuple):
        raise TypeError("source_cds must be a tuple of SourceCDS objects")
    for index, cds in enumerate(source_cds):
        if not isinstance(cds, SourceCDS):
            raise TypeError(
                "source_cds[{}] is not a SourceCDS".format(index)
            )
        if cds.end > len(source):
            raise ValueError(
                "source_cds[{}] lies outside the source".format(index)
            )


@dataclass(frozen=True)
class SourceIntactRecoveryGeometry:
    """Source-replayed candidate geometry with no runtime identity."""

    geometry_id: str
    interval: Tuple[int, int]
    repeat_intervals: Tuple[Tuple[int, int], ...]
    spacer_intervals: Tuple[Tuple[int, int], ...]
    repeat_count: int
    spacer_count: int
    source_slice_sha256: str
    repeat_source_slice_sha256s: Tuple[str, ...]
    spacer_source_slice_sha256s: Tuple[str, ...]

    def as_dict(self):
        return asdict(self)


def _candidate_geometry(source, candidate):
    try:
        reconstruction = validate_candidate_reconstruction(source, candidate)
    except (AttributeError, TypeError, ValueError):
        return None
    if not reconstruction.valid:
        return None
    repeats = tuple(candidate.list_repeats)
    starts = tuple(reconstruction.genomic_repeat_starts)
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    if not repeat_intervals:
        return None
    spacer_intervals = tuple(
        (left_end, right_start)
        for (_left_start, left_end), (right_start, _right_end)
        in zip(repeat_intervals, repeat_intervals[1:])
    )
    interval = (repeat_intervals[0][0], repeat_intervals[-1][1])
    geometry_payload = {
        "coordinate_system": COORDINATE_SYSTEM,
        "repeat_intervals": repeat_intervals,
        "spacer_intervals": spacer_intervals,
        "repeat_hashes": tuple(
            _slice_hash(source, item) for item in repeat_intervals
        ),
        "spacer_hashes": tuple(
            _slice_hash(source, item) for item in spacer_intervals
        ),
    }
    geometry_id = _source_hash(
        json.dumps(
            geometry_payload,
            separators=(",", ":"),
            sort_keys=True,
        )
    )
    return SourceIntactRecoveryGeometry(
        geometry_id=geometry_id,
        interval=interval,
        repeat_intervals=repeat_intervals,
        spacer_intervals=spacer_intervals,
        repeat_count=len(repeat_intervals),
        spacer_count=len(spacer_intervals),
        source_slice_sha256=_slice_hash(source, interval),
        repeat_source_slice_sha256s=tuple(
            _slice_hash(source, item) for item in repeat_intervals
        ),
        spacer_source_slice_sha256s=tuple(
            _slice_hash(source, item) for item in spacer_intervals
        ),
    )


def _overlap(first, second):
    start = max(first[0], second[0])
    end = min(first[1], second[1])
    return None if end <= start else (start, end)


def _relationship_kind(first, second):
    if first.repeat_intervals == second.repeat_intervals:
        return "identical_geometry"
    first_repeats = set(first.repeat_intervals)
    second_repeats = set(second.repeat_intervals)
    if first_repeats.issubset(second_repeats):
        return "first_exact_repeat_subpath"
    if second_repeats.issubset(first_repeats):
        return "second_exact_repeat_subpath"
    return "overlapping_boundary_or_phase_alternatives"


@dataclass(frozen=True)
class SourceIntactAlternativeRelationship:
    """Explicit source-validated relationship; it selects no representative."""

    source_sha256: str
    relationship: str
    overlap_interval: Tuple[int, int]
    first_geometry: SourceIntactRecoveryGeometry
    second_geometry: SourceIntactRecoveryGeometry
    boundary_status: str
    representative_selected: bool
    _first_row: object
    _second_row: object

    def as_dict(self):
        return {
            "source_sha256": self.source_sha256,
            "relationship": self.relationship,
            "overlap_interval": list(self.overlap_interval),
            "first_geometry": self.first_geometry.as_dict(),
            "second_geometry": self.second_geometry.as_dict(),
            "boundary_status": self.boundary_status,
            "representative_selected": self.representative_selected,
            "classifier_score_compared": False,
            "candidate_length_preference_used": False,
        }


def validate_source_intact_alternative_relationship(
    source_sequence, first_row, second_row
):
    """Build an immutable relation after replaying both rows from source."""

    source = _validate_source(source_sequence)
    _validate_scored_row(first_row)
    _validate_scored_row(second_row)
    if first_row is second_row:
        raise ValueError("an alternative relationship needs two distinct rows")
    first = _candidate_geometry(source, first_row[1])
    second = _candidate_geometry(source, second_row[1])
    if first is None or second is None:
        raise ValueError("alternative candidates must reconstruct from source")
    overlap = _overlap(first.interval, second.interval)
    if overlap is None:
        raise ValueError("disjoint arrays are not biological alternatives")
    ordered = sorted(
        ((first.geometry_id, first_row, first),
         (second.geometry_id, second_row, second)),
        key=lambda item: item[0],
    )
    (_first_id, ordered_first_row, ordered_first), (
        _second_id,
        ordered_second_row,
        ordered_second,
    ) = ordered
    return SourceIntactAlternativeRelationship(
        source_sha256=_source_hash(source),
        relationship=_relationship_kind(ordered_first, ordered_second),
        overlap_interval=_overlap(
            ordered_first.interval, ordered_second.interval
        ),
        first_geometry=ordered_first,
        second_geometry=ordered_second,
        boundary_status=BOUNDARY_AMBIGUOUS,
        representative_selected=False,
        _first_row=ordered_first_row,
        _second_row=ordered_second_row,
    )


@dataclass(frozen=True)
class SourceBackedRowFamilyEvidence:
    """Direct tandem evidence for one row and its exhaustive source family."""

    scored_row: object
    complete_family_candidates: Tuple[object, ...]
    tandem_decision: AcceptedFamilyTandemDecision


@dataclass(frozen=True)
class SourceIntactRecoveryCandidateContext:
    """One classifier-discarded row proposed for conservative recovery."""

    scored_row: object
    target_locus_key: object
    complete_family_candidates: Tuple[object, ...]
    tandem_decision: AcceptedFamilyTandemDecision


@dataclass(frozen=True)
class SourceIntactRecoveryLocusContext:
    """Explicit completeness contract for one recovery review group."""

    context_key: object
    recovery_candidates: Tuple[SourceIntactRecoveryCandidateContext, ...]
    overlapping_accepted: Tuple[SourceBackedRowFamilyEvidence, ...]
    alternative_relationships: Tuple[
        SourceIntactAlternativeRelationship, ...
    ] = ()
    family_complete: bool = True
    lineage_complete: bool = True
    augmentation_complete: bool = True
    accepted_overlap_complete: bool = True

    def __post_init__(self):
        if not isinstance(self.recovery_candidates, tuple):
            raise TypeError("recovery_candidates must be an explicit tuple")
        if not isinstance(self.overlapping_accepted, tuple):
            raise TypeError("overlapping_accepted must be an explicit tuple")
        if not isinstance(self.alternative_relationships, tuple):
            raise TypeError(
                "alternative_relationships must be an explicit tuple"
            )
        for field_name in (
            "family_complete",
            "lineage_complete",
            "augmentation_complete",
            "accepted_overlap_complete",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise TypeError("{} must be boolean".format(field_name))


@dataclass(frozen=True)
class SourceIntactRecoveryEvidence:
    """Exhaustive source-only authorization evidence for one recovery row."""

    source_sha256: str
    geometry: Optional[SourceIntactRecoveryGeometry]
    family_geometries: Tuple[SourceIntactRecoveryGeometry, ...]
    existence: SourceIntactExistenceDecision
    structural_filter: Optional[AdvancedFuzzySearchFilterEvidence]
    terminal_authorization: Optional[SourceIntactTerminalAuthorization]
    supplied_tandem: object
    recomputed_tandem: Optional[AcceptedFamilyTandemDecision]
    invalid_evidence_codes: Tuple[str, ...]
    authorized: bool

    def as_dict(self):
        return {
            "source_sha256": self.source_sha256,
            "geometry": (
                None if self.geometry is None else self.geometry.as_dict()
            ),
            "family_geometries": [
                geometry.as_dict() for geometry in self.family_geometries
            ],
            "existence": self.existence.as_dict(),
            "structural_filter": (
                None
                if self.structural_filter is None
                else self.structural_filter.as_dict()
            ),
            "terminal_authorization": (
                None
                if self.terminal_authorization is None
                else self.terminal_authorization.as_dict()
            ),
            "supplied_tandem": (
                self.supplied_tandem.as_dict()
                if isinstance(
                    self.supplied_tandem, AcceptedFamilyTandemDecision
                )
                else None
            ),
            "recomputed_tandem": (
                None
                if self.recomputed_tandem is None
                else self.recomputed_tandem.as_dict()
            ),
            "invalid_evidence_codes": list(self.invalid_evidence_codes),
            "authorized": self.authorized,
            "terminal_status_required": "explicit_policy_authorization",
            "boundary_status": BOUNDARY_AMBIGUOUS,
            "exact_boundary_certified": False,
            "classifier_score_used": False,
            "candidate_length_preference_used": False,
        }


@dataclass(frozen=True)
class SourceIntactRecoveryAcceptedTarget:
    category_index: int
    category_name: str
    locus_key_json: str
    row_ordinal: int
    geometry: SourceIntactRecoveryGeometry

    @property
    def locus_key(self):
        return json.loads(self.locus_key_json)

    def as_dict(self):
        return {
            "category_index": self.category_index,
            "category_name": self.category_name,
            "locus_key": self.locus_key,
            "row_ordinal": self.row_ordinal,
            "geometry": self.geometry.as_dict(),
        }


@dataclass(frozen=True)
class SourceIntactRecoveryRecord:
    """One runtime-object-free recovery decision and transition record."""

    action: str
    reason_codes: Tuple[str, ...]
    context_key_json: str
    source_category_index: int
    source_category_name: str
    source_locus_key_json: str
    source_row_ordinal: int
    target_locus_key_json: str
    score: float
    evidence: SourceIntactRecoveryEvidence
    accepted_targets: Tuple[SourceIntactRecoveryAcceptedTarget, ...]
    conflict_decisions: Tuple[AcceptedSourceConflictDecision, ...]
    alternative_relationships: Tuple[
        SourceIntactAlternativeRelationship, ...
    ]

    @property
    def context_key(self):
        return json.loads(self.context_key_json)

    @property
    def source_locus_key(self):
        return json.loads(self.source_locus_key_json)

    @property
    def target_locus_key(self):
        return json.loads(self.target_locus_key_json)

    def as_dict(self):
        return {
            "rule": {"name": RULE_NAME, "version": RULE_VERSION},
            "action": self.action,
            "reason_codes": list(self.reason_codes),
            "context_key": self.context_key,
            "source_category_index": self.source_category_index,
            "source_category_name": self.source_category_name,
            "source_locus_key": self.source_locus_key,
            "source_row_ordinal": self.source_row_ordinal,
            "target_category_index": (
                POSSIBLE_CATEGORY_INDEX
                if self.action in (
                    ACTION_ADMIT_POSSIBLE,
                    ACTION_REPLACE_POSSIBLE,
                )
                else None
            ),
            "target_category_name": (
                CATEGORY_NAMES[POSSIBLE_CATEGORY_INDEX]
                if self.action in (
                    ACTION_ADMIT_POSSIBLE,
                    ACTION_REPLACE_POSSIBLE,
                )
                else None
            ),
            "target_locus_key": self.target_locus_key,
            "score": self.score,
            "evidence": self.evidence.as_dict(),
            "accepted_targets": [
                target.as_dict() for target in self.accepted_targets
            ],
            "conflict_decisions": [
                decision.as_dict() for decision in self.conflict_decisions
            ],
            "alternative_relationships": [
                relationship.as_dict()
                for relationship in self.alternative_relationships
            ],
            "boundary_status": BOUNDARY_AMBIGUOUS,
            "exact_boundary_certified": False,
            "classifier_score_used": False,
            "candidate_length_preference_used": False,
        }


@dataclass(frozen=True)
class _RowLocation:
    category_index: int
    locus_key: object
    locus_key_json: str
    row_ordinal: int
    row: object


@dataclass(frozen=True)
class _PlanDraft:
    record: SourceIntactRecoveryRecord
    recovery_row: object
    accepted_rows: Tuple[object, ...]
    _target_locus_key: object


@dataclass(frozen=True)
class SourceIntactRecoveryPlan:
    source_sha256: str
    drafts: Tuple[_PlanDraft, ...]
    cost_audit: object
    _category_snapshot: tuple

    @property
    def records(self):
        return tuple(draft.record for draft in self.drafts)

    def as_dict(self):
        return {
            "rule": {"name": RULE_NAME, "version": RULE_VERSION},
            "source_sha256": self.source_sha256,
            "cost_audit": self.cost_audit.as_dict(),
            "records": [record.as_dict() for record in self.records],
        }


@dataclass(frozen=True)
class SourceIntactRecoveryCostAudit:
    """Explicit full-source work remaining after context reuse."""

    source_context_full_hashes: int
    context_bound_strict_gate_evaluations: int
    direct_tandem_recomputations: int
    accepted_conflict_arbitrations: int
    estimated_residual_full_source_hashes_during_planning: int
    residual_cost_note: str

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class SourceIntactRecoveryStats:
    total_rows_before: int
    total_rows_after: int
    reviewed: int
    admitted_possible: int
    replaced_possible: int
    alternatives_attached: int
    abstained: int
    rejected: int
    accepted_rows_before: int
    accepted_rows_after: int

    def as_dict(self):
        return asdict(self)


def _validate_scored_row(row):
    if not isinstance(row, (list, tuple)) or len(row) not in (3, 4):
        raise ValueError(
            "category rows must contain score, candidate, features, "
            "and optional metadata"
        )
    if (
        isinstance(row[0], bool)
        or not isinstance(row[0], Real)
        or not math.isfinite(float(row[0]))
    ):
        raise ValueError("category row scores must be finite numbers")
    if row[1] is None:
        raise ValueError("category rows must contain a candidate")
    if not isinstance(row[2], (list, tuple)):
        raise ValueError("category row features must be a sequence")
    if len(row) == 4:
        if not isinstance(row[3], Mapping):
            raise ValueError("category row metadata must be a mapping")
        prior = row[3].get("classification_basis")
        if prior is not None:
            if not isinstance(prior, Mapping):
                raise ValueError("classification_basis must be a mapping")
            _json_copy(prior, "classification_basis")


def _candidate_state(candidate):
    fields = []
    for name in (
        "list_repeat_starts",
        "list_repeats",
        "list_spacers",
        "list_gaped_repeats",
        "consensus",
    ):
        value = getattr(candidate, name, None)
        if isinstance(value, list):
            value = tuple(value)
        fields.append((name, value))
    return tuple(fields)


def _index_categories(categories):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("source-intact recovery requires five categories")
    locations = []
    row_ids = set()
    candidate_ids = set()
    snapshot = []
    for category_index, category in enumerate(categories):
        if not isinstance(category, Mapping):
            raise ValueError("each category must be a mapping")
        category_snapshot = []
        for locus_key, rows in category.items():
            locus_json = _json_key(locus_key, "category locus key")
            if not isinstance(rows, (list, tuple)) or not rows:
                raise ValueError(
                    "category locus values must be non-empty row sequences"
                )
            row_snapshot = []
            for ordinal, row in enumerate(rows):
                _validate_scored_row(row)
                if id(row) in row_ids:
                    raise ValueError("one scored row occurs more than once")
                if id(row[1]) in candidate_ids:
                    raise ValueError(
                        "one candidate object has multiple category rows"
                    )
                row_ids.add(id(row))
                candidate_ids.add(id(row[1]))
                location = _RowLocation(
                    category_index,
                    locus_key,
                    locus_json,
                    ordinal,
                    row,
                )
                locations.append(location)
                row_snapshot.append((
                    row,
                    len(row),
                    float(row[0]),
                    row[1],
                    _candidate_state(row[1]),
                    _json_copy(dict(row[3]), "row metadata")
                    if len(row) == 4 else None,
                ))
            category_snapshot.append((locus_json, tuple(row_snapshot)))
        snapshot.append(tuple(category_snapshot))
    return tuple(locations), tuple(snapshot)


def _snapshot_matches(categories, expected):
    _locations, observed = _index_categories(categories)
    if len(observed) != len(expected):
        return False
    for expected_category, observed_category in zip(expected, observed):
        if len(expected_category) != len(observed_category):
            return False
        for expected_locus, observed_locus in zip(
            expected_category, observed_category
        ):
            if expected_locus[0] != observed_locus[0]:
                return False
            if len(expected_locus[1]) != len(observed_locus[1]):
                return False
            for expected_row, observed_row in zip(
                expected_locus[1], observed_locus[1]
            ):
                if (
                    expected_row[0] is not observed_row[0]
                    or expected_row[1:] != observed_row[1:]
                ):
                    return False
    return True


def _family_evidence(
    source,
    evaluation_context,
    row,
    complete_family_candidates,
    supplied_tandem,
    structural_filter,
    source_cds,
    existence_policy,
    cost_counts,
    *,
    recovery_authorization,
):
    existence = evaluate_source_intact_existence_in_context(
        evaluation_context,
        row[1],
        policy=existence_policy,
    )
    cost_counts["strict_gate"] += 1
    geometry = _candidate_geometry(source, row[1])
    invalid_codes = []
    if geometry is None:
        invalid_codes.append("candidate_source_reconstruction_failed")

    family_geometries = []
    family_valid = True
    if (
        not isinstance(complete_family_candidates, tuple)
        or not complete_family_candidates
    ):
        invalid_codes.append("complete_family_missing")
        family_valid = False
        family = ()
    else:
        family = complete_family_candidates
        family_ids = tuple(id(candidate) for candidate in family)
        if len(set(family_ids)) != len(family_ids):
            invalid_codes.append("complete_family_contains_duplicate_object")
            family_valid = False
        if sum(candidate is row[1] for candidate in family) != 1:
            invalid_codes.append(
                "complete_family_does_not_include_row_candidate_once"
            )
            family_valid = False
        for candidate in family:
            candidate_geometry = _candidate_geometry(source, candidate)
            if candidate_geometry is None:
                invalid_codes.append(
                    "complete_family_candidate_source_reconstruction_failed"
                )
                family_valid = False
            else:
                family_geometries.append(candidate_geometry)

    try:
        structural = structural_filter.evidence(row[1])
    except (AttributeError, IndexError, TypeError, ValueError, ZeroDivisionError):
        structural = None
        invalid_codes.append("structural_filter_evidence_failed")

    terminal_authorization = (
        None
        if structural is None
        else authorize_source_intact_terminal_evidence(
            existence,
            structural,
        )
    )
    should_recompute_tandem = (
        geometry is not None
        and family_valid
        and (
            not recovery_authorization
            or (
                terminal_authorization is not None
                and terminal_authorization.authorized
            )
        )
    )
    recomputed_tandem = None
    if should_recompute_tandem:
        if not isinstance(supplied_tandem, AcceptedFamilyTandemDecision):
            invalid_codes.append("direct_tandem_evidence_missing_or_wrong_type")
        try:
            recomputed_tandem = validate_accepted_family_tandem(
                source,
                row[1],
                family,
                source_cds,
            )
        except (AttributeError, TypeError, ValueError):
            invalid_codes.append("direct_tandem_recomputation_failed")
        else:
            cost_counts["tandem"] += 1
            if supplied_tandem != recomputed_tandem:
                invalid_codes.append("direct_tandem_evidence_stale_or_mismatched")
    authorized = (
        not invalid_codes
        and geometry is not None
        and existence.eligible
        and structural is not None
        and structural.passes
        and terminal_authorization is not None
        and terminal_authorization.authorized
        and recomputed_tandem is not None
        and not recomputed_tandem.supported
    )
    return SourceIntactRecoveryEvidence(
        source_sha256=evaluation_context.source_sha256,
        geometry=geometry,
        family_geometries=tuple(sorted(
            family_geometries,
            key=lambda item: (item.interval, item.geometry_id),
        )),
        existence=existence,
        structural_filter=structural,
        terminal_authorization=terminal_authorization,
        supplied_tandem=supplied_tandem,
        recomputed_tandem=recomputed_tandem,
        invalid_evidence_codes=tuple(sorted(set(invalid_codes))),
        authorized=authorized,
    )


def _authorization_reason_codes(evidence):
    reasons = list(evidence.invalid_evidence_codes)
    if not evidence.existence.eligible:
        reasons.append("strict_source_intact_existence_rejected")
    if evidence.structural_filter is None:
        reasons.append("structural_filter_evidence_unavailable")
    elif not evidence.structural_filter.passes:
        reasons.append("structural_filter_rejected")
    if evidence.terminal_authorization is None:
        reasons.append("terminal_evidence_authorization_unavailable")
    elif not evidence.terminal_authorization.authorized:
        reasons.append("terminal_evidence_not_authorized")
        reasons.extend(
            "terminal_authorization_{}".format(reason)
            for reason in evidence.terminal_authorization.reason_codes
        )
    if (
        evidence.recomputed_tandem is not None
        and evidence.recomputed_tandem.supported
    ):
        reasons.append("direct_coding_tandem_supported")
    return tuple(sorted(set(reasons)))


def _accepted_target(location, geometry):
    return SourceIntactRecoveryAcceptedTarget(
        category_index=location.category_index,
        category_name=CATEGORY_NAMES[location.category_index],
        locus_key_json=location.locus_key_json,
        row_ordinal=location.row_ordinal,
        geometry=geometry,
    )


def _relationship_pair_key(first_row, second_row):
    return frozenset((id(first_row), id(second_row)))


def _validated_context_relationships(
    source, source_sha256, context, recovery_rows
):
    invalid_codes = []
    relationships = []
    pair_keys = set()
    recovery_ids = {id(row): row for row in recovery_rows}
    for relationship in context.alternative_relationships:
        if not isinstance(
            relationship, SourceIntactAlternativeRelationship
        ):
            invalid_codes.append("alternative_relationship_wrong_type")
            continue
        first = recovery_ids.get(id(relationship._first_row))
        second = recovery_ids.get(id(relationship._second_row))
        if (
            first is not relationship._first_row
            or second is not relationship._second_row
        ):
            invalid_codes.append(
                "alternative_relationship_refers_to_foreign_row"
            )
            continue
        pair_key = _relationship_pair_key(first, second)
        if pair_key in pair_keys:
            invalid_codes.append("duplicate_alternative_relationship")
            continue
        pair_keys.add(pair_key)
        first_geometry = _candidate_geometry(source, first[1])
        second_geometry = _candidate_geometry(source, second[1])
        if first_geometry is None or second_geometry is None:
            invalid_codes.append("alternative_relationship_revalidation_failed")
            continue
        ordered = sorted(
            ((first_geometry.geometry_id, first_geometry),
             (second_geometry.geometry_id, second_geometry)),
            key=lambda item: item[0],
        )
        rebuilt_payload = {
            "source_sha256": source_sha256,
            "relationship": _relationship_kind(ordered[0][1], ordered[1][1]),
            "overlap_interval": list(_overlap(
                ordered[0][1].interval, ordered[1][1].interval
            ) or ()),
            "first_geometry": ordered[0][1].as_dict(),
            "second_geometry": ordered[1][1].as_dict(),
            "boundary_status": BOUNDARY_AMBIGUOUS,
            "representative_selected": False,
            "classifier_score_compared": False,
            "candidate_length_preference_used": False,
        }
        if relationship.as_dict() != rebuilt_payload:
            invalid_codes.append("alternative_relationship_stale_or_mismatched")
            continue
        relationships.append(relationship)
    return (
        tuple(sorted(
            relationships,
            key=lambda item: (
                item.first_geometry.geometry_id,
                item.second_geometry.geometry_id,
            ),
        )),
        frozenset(pair_keys),
        tuple(sorted(set(invalid_codes))),
    )


def _context_completeness_reasons(context):
    reasons = []
    for field_name, reason in (
        ("family_complete", "family_context_incomplete"),
        ("lineage_complete", "lineage_context_incomplete"),
        ("augmentation_complete", "augmentation_context_incomplete"),
        (
            "accepted_overlap_complete",
            "accepted_overlap_context_not_attested_complete",
        ),
    ):
        if not getattr(context, field_name):
            reasons.append(reason)
    return tuple(reasons)


def _record(
    *,
    action,
    reasons,
    context_key_json,
    location,
    target_locus_key_json,
    evidence,
    accepted_targets=(),
    conflicts=(),
    relationships=(),
):
    return SourceIntactRecoveryRecord(
        action=action,
        reason_codes=tuple(sorted(set(reasons))),
        context_key_json=context_key_json,
        source_category_index=location.category_index,
        source_category_name=CATEGORY_NAMES[location.category_index],
        source_locus_key_json=location.locus_key_json,
        source_row_ordinal=location.row_ordinal,
        target_locus_key_json=target_locus_key_json,
        score=float(location.row[0]),
        evidence=evidence,
        accepted_targets=tuple(accepted_targets),
        conflict_decisions=tuple(conflicts),
        alternative_relationships=tuple(relationships),
    )


def _accepted_category_source_state(source, locations):
    geometries = {}
    invalid = []
    for location in locations:
        if location.category_index not in ACCEPTED_CATEGORY_INDEXES:
            continue
        geometry = _candidate_geometry(source, location.row[1])
        if geometry is None:
            invalid.append(location)
        else:
            geometries[id(location.row)] = geometry
    return geometries, tuple(invalid)


def _expected_overlapping_accepted(
    recovery_geometries, accepted_locations, accepted_geometries
):
    expected = []
    for accepted in accepted_locations:
        accepted_geometry = accepted_geometries.get(id(accepted.row))
        if accepted_geometry is None:
            continue
        if any(
            geometry is not None
            and _overlap(geometry.interval, accepted_geometry.interval)
            is not None
            for geometry in recovery_geometries
        ):
            expected.append(accepted)
    return tuple(expected)


def plan_source_intact_recoveries(
    source_sequence,
    categories,
    locus_contexts,
    structural_filter,
    source_cds,
    *,
    existence_policy=SOURCE_INTACT_EXISTENCE_POLICY,
):
    """Return a frozen, non-mutating recovery plan.

    A recovery row must currently occur in ``Possible Discarded`` or ``Low
    score``.  Contexts need not cover unrelated rows, but each context must
    enumerate every accepted row overlapping any of its recovery candidates.
    """

    source = _validate_source(source_sequence)
    evaluation_context = SourceIntactEvaluationContext(source)
    _validate_source_cds(source, source_cds)
    if not isinstance(structural_filter, AdvancedFuzzySearchFilter):
        raise TypeError("structural_filter must be AdvancedFuzzySearchFilter")
    if not isinstance(existence_policy, SourceIntactExistencePolicy):
        raise TypeError("existence_policy must be SourceIntactExistencePolicy")
    locations, category_snapshot = _index_categories(categories)
    location_by_id = {id(location.row): location for location in locations}
    accepted_locations = tuple(
        sorted(
            (
                location
                for location in locations
                if location.category_index in ACCEPTED_CATEGORY_INDEXES
            ),
            key=lambda location: (
                location.category_index,
                location.locus_key_json,
                location.row_ordinal,
            ),
        )
    )
    accepted_geometries, invalid_accepted = (
        _accepted_category_source_state(source, locations)
    )

    if isinstance(locus_contexts, (str, bytes)) or locus_contexts is None:
        raise TypeError("locus_contexts must be an iterable of contexts")
    try:
        locus_contexts = tuple(locus_contexts)
    except TypeError as error:
        raise TypeError(
            "locus_contexts must be an iterable of contexts"
        ) from error

    seen_recovery_rows = set()
    drafts = []
    cost_counts = {"strict_gate": 0, "tandem": 0, "conflict": 0}
    ordered_contexts = []
    for index, context in enumerate(locus_contexts):
        if not isinstance(context, SourceIntactRecoveryLocusContext):
            raise TypeError(
                "locus_contexts[{}] has the wrong type".format(index)
            )
        if not context.recovery_candidates:
            raise ValueError("a recovery context needs at least one candidate")
        context_key_json = _json_key(
            context.context_key, "recovery context key"
        )
        ordered_contexts.append((context_key_json, context))

    for context_key_json, context in sorted(
        ordered_contexts, key=lambda item: item[0]
    ):
        context_reasons = list(_context_completeness_reasons(context))
        recovery_items = []
        recovery_rows = []
        for candidate_context in context.recovery_candidates:
            if not isinstance(
                candidate_context, SourceIntactRecoveryCandidateContext
            ):
                raise TypeError("recovery candidate context has wrong type")
            row = candidate_context.scored_row
            if id(row) in seen_recovery_rows:
                raise ValueError(
                    "a recovery row occurs in more than one context"
                )
            seen_recovery_rows.add(id(row))
            location = location_by_id.get(id(row))
            if location is None or location.row is not row:
                raise ValueError("recovery context refers to a foreign row")
            if location.category_index not in RECOVERY_SOURCE_CATEGORY_INDEXES:
                raise ValueError(
                    "recovery rows must be classifier-discarded rows"
                )
            target_locus_key_json = _json_key(
                candidate_context.target_locus_key,
                "recovery target locus key",
            )
            try:
                hash(candidate_context.target_locus_key)
            except TypeError as error:
                raise ValueError(
                    "recovery target locus key must be hashable"
                ) from error
            evidence = _family_evidence(
                source,
                evaluation_context,
                row,
                candidate_context.complete_family_candidates,
                candidate_context.tandem_decision,
                structural_filter,
                source_cds,
                existence_policy,
                cost_counts,
                recovery_authorization=True,
            )
            recovery_items.append((
                candidate_context,
                location,
                target_locus_key_json,
                evidence,
            ))
            recovery_rows.append(row)

        relationships, relation_pairs, relationship_errors = (
            _validated_context_relationships(
                source,
                evaluation_context.source_sha256,
                context,
                tuple(recovery_rows),
            )
        )
        context_reasons.extend(relationship_errors)

        accepted_evidence_by_id = {}
        for accepted_context in context.overlapping_accepted:
            if not isinstance(
                accepted_context, SourceBackedRowFamilyEvidence
            ):
                raise TypeError("accepted row family evidence has wrong type")
            row = accepted_context.scored_row
            location = location_by_id.get(id(row))
            if (
                location is None
                or location.row is not row
                or location.category_index not in ACCEPTED_CATEGORY_INDEXES
            ):
                raise ValueError(
                    "accepted overlap context refers to a non-accepted row"
                )
            if id(row) in accepted_evidence_by_id:
                context_reasons.append("duplicate_accepted_overlap_evidence")
                continue
            evidence = _family_evidence(
                source,
                evaluation_context,
                row,
                accepted_context.complete_family_candidates,
                accepted_context.tandem_decision,
                structural_filter,
                source_cds,
                existence_policy,
                cost_counts,
                recovery_authorization=False,
            )
            accepted_evidence_by_id[id(row)] = (
                accepted_context,
                location,
                evidence,
            )

        recovery_geometries = tuple(
            evidence.geometry for _ctx, _loc, _key, evidence in recovery_items
        )
        expected_accepted = _expected_overlapping_accepted(
            recovery_geometries,
            accepted_locations,
            accepted_geometries,
        )
        expected_ids = {id(location.row) for location in expected_accepted}
        supplied_ids = set(accepted_evidence_by_id)
        if expected_ids != supplied_ids:
            context_reasons.append("accepted_overlap_context_incomplete")
        if invalid_accepted:
            context_reasons.append(
                "accepted_category_contains_unreconstructable_candidate"
            )

        peer_overlaps = {id(row): [] for row in recovery_rows}
        for left_index, left_item in enumerate(recovery_items):
            left_row = left_item[1].row
            left_evidence = left_item[3]
            if not left_evidence.authorized:
                continue
            for right_item in recovery_items[left_index + 1 :]:
                right_row = right_item[1].row
                right_evidence = right_item[3]
                if not right_evidence.authorized:
                    continue
                overlap = _overlap(
                    left_evidence.geometry.interval,
                    right_evidence.geometry.interval,
                )
                if overlap is None:
                    continue
                explicit = (
                    _relationship_pair_key(left_row, right_row)
                    in relation_pairs
                )
                peer_overlaps[id(left_row)].append((right_row, explicit))
                peer_overlaps[id(right_row)].append((left_row, explicit))

        for (
            candidate_context,
            location,
            target_locus_key_json,
            evidence,
        ) in recovery_items:
            reasons = list(context_reasons)
            if not evidence.authorized:
                reasons.extend(_authorization_reason_codes(evidence))
            if reasons:
                action = ACTION_REJECT
                record = _record(
                    action=action,
                    reasons=reasons,
                    context_key_json=context_key_json,
                    location=location,
                    target_locus_key_json=target_locus_key_json,
                    evidence=evidence,
                    relationships=relationships,
                )
                drafts.append(_PlanDraft(
                    record,
                    location.row,
                    (),
                    candidate_context.target_locus_key,
                ))
                continue

            peers = peer_overlaps[id(location.row)]
            if peers:
                all_explicit = all(explicit for _row, explicit in peers)
                reason = (
                    "explicit_biological_alternatives_abstained_without_"
                    "authorized_representative"
                    if all_explicit
                    else "overlapping_recovery_geometries_unresolved"
                )
                record = _record(
                    action=ACTION_ABSTAIN,
                    reasons=(reason,),
                    context_key_json=context_key_json,
                    location=location,
                    target_locus_key_json=target_locus_key_json,
                    evidence=evidence,
                    relationships=relationships,
                )
                drafts.append(_PlanDraft(
                    record,
                    location.row,
                    (),
                    candidate_context.target_locus_key,
                ))
                continue

            overlapping = []
            for accepted_location in expected_accepted:
                accepted_geometry = accepted_geometries[id(
                    accepted_location.row
                )]
                if _overlap(
                    evidence.geometry.interval, accepted_geometry.interval
                ) is not None:
                    overlapping.append(accepted_location)

            accepted_targets = tuple(
                _accepted_target(
                    accepted_location,
                    accepted_geometries[id(accepted_location.row)],
                )
                for accepted_location in overlapping
            )
            invalid_accepted_evidence = []
            for accepted_location in overlapping:
                accepted_evidence = accepted_evidence_by_id[
                    id(accepted_location.row)
                ][2]
                if accepted_evidence.invalid_evidence_codes:
                    invalid_accepted_evidence.extend(
                        "accepted_{}".format(code)
                        for code in accepted_evidence.invalid_evidence_codes
                    )
                if (
                    accepted_evidence.recomputed_tandem is not None
                    and accepted_evidence.recomputed_tandem.supported
                ):
                    invalid_accepted_evidence.append(
                        "accepted_direct_tandem_demotion_must_run_first"
                    )
            if invalid_accepted_evidence:
                record = _record(
                    action=ACTION_REJECT,
                    reasons=invalid_accepted_evidence,
                    context_key_json=context_key_json,
                    location=location,
                    target_locus_key_json=target_locus_key_json,
                    evidence=evidence,
                    accepted_targets=accepted_targets,
                    relationships=relationships,
                )
                drafts.append(_PlanDraft(
                    record,
                    location.row,
                    (),
                    candidate_context.target_locus_key,
                ))
                continue

            if not overlapping:
                record = _record(
                    action=ACTION_ADMIT_POSSIBLE,
                    reasons=("independent_source_intact_array",),
                    context_key_json=context_key_json,
                    location=location,
                    target_locus_key_json=target_locus_key_json,
                    evidence=evidence,
                    relationships=relationships,
                )
                drafts.append(_PlanDraft(
                    record,
                    location.row,
                    (),
                    candidate_context.target_locus_key,
                ))
                continue

            conflicts = []
            for accepted_location in overlapping:
                accepted_context, _location, accepted_evidence = (
                    accepted_evidence_by_id[id(accepted_location.row)]
                )
                conflicts.append(arbitrate_accepted_source_intact_conflict(
                    source,
                    accepted_location.row[1],
                    location.row[1],
                    accepted_category=CATEGORY_NAMES[
                        accepted_location.category_index
                    ],
                    structural_filter=structural_filter,
                    accepted_tandem_decision=(
                        accepted_evidence.recomputed_tandem
                    ),
                    source_intact_tandem_decision=(
                        evidence.recomputed_tandem
                    ),
                    existence_policy=existence_policy,
                ))
                cost_counts["conflict"] += 1

            if len(overlapping) > 1:
                record = _record(
                    action=ACTION_ATTACH_ALTERNATIVE,
                    reasons=("multiple_accepted_overlap_complex_abstention",),
                    context_key_json=context_key_json,
                    location=location,
                    target_locus_key_json=target_locus_key_json,
                    evidence=evidence,
                    accepted_targets=accepted_targets,
                    conflicts=conflicts,
                    relationships=relationships,
                )
                drafts.append(_PlanDraft(
                    record,
                    location.row,
                    tuple(item.row for item in overlapping),
                    candidate_context.target_locus_key,
                ))
                continue

            conflict = conflicts[0]
            accepted_location = overlapping[0]
            if not conflict.source_authorized:
                action = ACTION_REJECT
                reasons = ("accepted_conflict_arbitration_rejected_source",)
                accepted_rows = ()
            elif conflict.source_intact_action == SOURCE_ATTACH_ALTERNATIVE:
                action = ACTION_ATTACH_ALTERNATIVE
                reasons = (conflict.decision,)
                accepted_rows = (accepted_location.row,)
            elif conflict.source_intact_action == SOURCE_ADMIT_POSSIBLE:
                if conflict.accepted_action == ACCEPTED_RETIRE:
                    if accepted_location.category_index != POSSIBLE_CATEGORY_INDEX:
                        action = ACTION_REJECT
                        reasons = (
                            "conflict_requested_non_possible_retirement",
                        )
                        accepted_rows = ()
                    else:
                        action = ACTION_REPLACE_POSSIBLE
                        reasons = (conflict.decision,)
                        accepted_rows = (accepted_location.row,)
                elif conflict.accepted_action == ACCEPTED_RETAIN:
                    action = ACTION_ADMIT_POSSIBLE
                    reasons = (conflict.decision,)
                    accepted_rows = ()
                else:
                    action = ACTION_REJECT
                    reasons = ("unsupported_accepted_conflict_action",)
                    accepted_rows = ()
            elif conflict.decision == "duplicate_source_geometry":
                action = ACTION_ATTACH_ALTERNATIVE
                reasons = ("identical_source_geometry_already_represented",)
                accepted_rows = (accepted_location.row,)
            else:
                action = ACTION_REJECT
                reasons = ("source_conflict_did_not_authorize_admission",)
                accepted_rows = ()
            record = _record(
                action=action,
                reasons=reasons,
                context_key_json=context_key_json,
                location=location,
                target_locus_key_json=target_locus_key_json,
                evidence=evidence,
                accepted_targets=accepted_targets,
                conflicts=conflicts,
                relationships=relationships,
            )
            drafts.append(_PlanDraft(
                record,
                location.row,
                accepted_rows,
                candidate_context.target_locus_key,
            ))

    drafts = tuple(sorted(
        drafts,
        key=lambda draft: (
            draft.record.context_key_json,
            draft.record.evidence.geometry.interval
            if draft.record.evidence.geometry is not None
            else (math.inf, math.inf),
            draft.record.evidence.geometry.geometry_id
            if draft.record.evidence.geometry is not None
            else "",
            draft.record.source_category_index,
            draft.record.source_locus_key_json,
            draft.record.source_row_ordinal,
        ),
    ))
    return SourceIntactRecoveryPlan(
        source_sha256=evaluation_context.source_sha256,
        drafts=drafts,
        cost_audit=SourceIntactRecoveryCostAudit(
            source_context_full_hashes=1,
            context_bound_strict_gate_evaluations=(
                cost_counts["strict_gate"]
            ),
            direct_tandem_recomputations=cost_counts["tandem"],
            accepted_conflict_arbitrations=cost_counts["conflict"],
            estimated_residual_full_source_hashes_during_planning=(
                cost_counts["tandem"] + 5 * cost_counts["conflict"]
            ),
            residual_cost_note=(
                "accepted-family tandem currently hashes the source once per "
                "recomputation; accepted/source conflict currently performs "
                "an estimated five full-source hashes per arbitration"
            ),
        ),
        _category_snapshot=category_snapshot,
    )


def _copy_categories(categories):
    return [
        {locus_key: list(rows) for locus_key, rows in category.items()}
        for category in categories
    ]


def _remove_row(new_categories, location, row):
    rows = new_categories[location.category_index].get(location.locus_key)
    if rows is None:
        raise RuntimeError("planned row locus disappeared during transaction")
    matches = [index for index, item in enumerate(rows) if item is row]
    if len(matches) != 1:
        raise RuntimeError("planned row identity changed during transaction")
    rows.pop(matches[0])
    if not rows:
        del new_categories[location.category_index][location.locus_key]


def _rebuilt_row(row, metadata):
    if isinstance(row, tuple):
        return tuple(row[:3]) + (metadata,)
    return list(row[:3]) + [metadata]


def _prior_basis(row):
    metadata = dict(row[3]) if len(row) == 4 else {}
    prior = metadata.get("classification_basis")
    if prior is not None:
        prior = _json_copy(prior, "prior classification_basis")
    return metadata, prior


def _admitted_row(row, record):
    metadata, prior = _prior_basis(row)
    metadata["classification_basis"] = {
        "method": CLASSIFICATION_METHOD,
        "prior_classification_basis": prior,
        "source_intact_recovery": _json_copy(
            record.as_dict(), "source-intact recovery record"
        ),
        "boundary_status": BOUNDARY_AMBIGUOUS,
        "exact_boundary_certified": False,
    }
    return _rebuilt_row(row, metadata)


def _retired_possible_row(row, record):
    metadata, prior = _prior_basis(row)
    metadata["classification_basis"] = {
        "method": "source_backed_intact_terminal_conflict_retirement",
        "prior_classification_basis": prior,
        "source_intact_recovery": _json_copy(
            record.as_dict(), "source-intact recovery record"
        ),
        "transition": {
            "source_category_index": POSSIBLE_CATEGORY_INDEX,
            "source_category_name": CATEGORY_NAMES[POSSIBLE_CATEGORY_INDEX],
            "target_category_index": DISCARDED_CATEGORY_INDEX,
            "target_category_name": CATEGORY_NAMES[DISCARDED_CATEGORY_INDEX],
        },
    }
    return _rebuilt_row(row, metadata)


def _alternative_attached_row(row, records):
    metadata = dict(row[3]) if len(row) == 4 else {}
    prior = metadata.get("source_intact_recovery_alternatives", [])
    if not isinstance(prior, (list, tuple)):
        raise ValueError(
            "source_intact_recovery_alternatives must be a sequence"
        )
    alternatives = list(_json_copy(
        list(prior), "prior source-intact alternatives"
    ))
    alternatives.extend(
        _json_copy(record.as_dict(), "source-intact alternative")
        for record in records
    )
    alternatives.sort(key=lambda item: (
        item["evidence"]["geometry"]["interval"],
        item["evidence"]["geometry"]["geometry_id"],
        _json_key(item["context_key"], "alternative context key"),
    ))
    metadata["source_intact_recovery_alternatives"] = alternatives
    return _rebuilt_row(row, metadata)


def _count_rows(categories, indexes=None):
    category_indexes = range(len(categories)) if indexes is None else indexes
    return sum(
        len(rows)
        for index in category_indexes
        for rows in categories[index].values()
    )


def apply_source_intact_recovery_plan(
    source_sequence, categories, plan
):
    """Apply a current plan to a copy of ``categories``.

    Returns ``(new_categories, records, stats)``.  No row is added or removed
    globally: admissions and retirements are category moves, while ambiguous
    alternatives rebuild only their accepted anchor metadata.
    """

    source = _validate_source(source_sequence)
    if not isinstance(plan, SourceIntactRecoveryPlan):
        raise TypeError("plan must be a SourceIntactRecoveryPlan")
    if _source_hash(source) != plan.source_sha256:
        raise ValueError("plan source hash does not match the current source")
    if not _snapshot_matches(categories, plan._category_snapshot):
        raise ValueError("categories or candidate geometry changed after planning")

    locations, _snapshot = _index_categories(categories)
    location_by_id = {id(location.row): location for location in locations}
    replacement_by_accepted = {}
    alternatives_by_accepted = {}
    for draft in plan.drafts:
        location = location_by_id.get(id(draft.recovery_row))
        if location is None or location.row is not draft.recovery_row:
            raise RuntimeError("planned recovery row is absent")
        if draft.record.action == ACTION_REPLACE_POSSIBLE:
            if len(draft.accepted_rows) != 1:
                raise RuntimeError("replacement needs exactly one accepted row")
            accepted_row = draft.accepted_rows[0]
            accepted_location = location_by_id.get(id(accepted_row))
            if (
                accepted_location is None
                or accepted_location.row is not accepted_row
                or accepted_location.category_index != POSSIBLE_CATEGORY_INDEX
            ):
                raise RuntimeError("replacement target is not a current Possible")
            if id(accepted_row) in replacement_by_accepted:
                raise RuntimeError("one Possible has multiple replacement plans")
            replacement_by_accepted[id(accepted_row)] = draft
        elif draft.record.action == ACTION_ATTACH_ALTERNATIVE:
            if not draft.accepted_rows:
                raise RuntimeError("alternative attachment needs an accepted row")
            for accepted_row in draft.accepted_rows:
                accepted_location = location_by_id.get(id(accepted_row))
                if (
                    accepted_location is None
                    or accepted_location.row is not accepted_row
                    or accepted_location.category_index
                    not in ACCEPTED_CATEGORY_INDEXES
                ):
                    raise RuntimeError("alternative anchor is not accepted")
                alternatives_by_accepted.setdefault(
                    id(accepted_row), []
                ).append(draft)
        elif draft.record.action in (
            ACTION_ADMIT_POSSIBLE,
            ACTION_ABSTAIN,
            ACTION_REJECT,
        ):
            if draft.accepted_rows:
                raise RuntimeError("non-anchor action retained accepted rows")
        else:
            raise RuntimeError("unknown source-intact recovery action")

    collision = set(replacement_by_accepted).intersection(
        alternatives_by_accepted
    )
    if collision:
        raise RuntimeError(
            "one accepted row cannot be replaced and annotated together"
        )

    new_categories = _copy_categories(categories)
    for draft in plan.drafts:
        record = draft.record
        if record.action not in (
            ACTION_ADMIT_POSSIBLE,
            ACTION_REPLACE_POSSIBLE,
        ):
            continue
        recovery_location = location_by_id[id(draft.recovery_row)]
        _remove_row(new_categories, recovery_location, draft.recovery_row)
        admitted = _admitted_row(draft.recovery_row, record)
        new_categories[POSSIBLE_CATEGORY_INDEX].setdefault(
            draft._target_locus_key, []
        ).append(admitted)

        if record.action == ACTION_REPLACE_POSSIBLE:
            accepted_row = draft.accepted_rows[0]
            accepted_location = location_by_id[id(accepted_row)]
            _remove_row(new_categories, accepted_location, accepted_row)
            retired = _retired_possible_row(accepted_row, record)
            new_categories[DISCARDED_CATEGORY_INDEX].setdefault(
                accepted_location.locus_key, []
            ).append(retired)

    for accepted_id, accepted_drafts in sorted(
        alternatives_by_accepted.items(),
        key=lambda item: (
            location_by_id[item[0]].category_index,
            location_by_id[item[0]].locus_key_json,
            location_by_id[item[0]].row_ordinal,
        ),
    ):
        accepted_location = location_by_id[accepted_id]
        accepted_row = accepted_location.row
        rows = new_categories[accepted_location.category_index][
            accepted_location.locus_key
        ]
        matches = [index for index, row in enumerate(rows) if row is accepted_row]
        if len(matches) != 1:
            raise RuntimeError("alternative anchor identity changed")
        sorted_drafts = sorted(
            accepted_drafts,
            key=lambda item: (
                item.record.evidence.geometry.interval,
                item.record.evidence.geometry.geometry_id,
                item.record.context_key_json,
            ),
        )
        rows[matches[0]] = _alternative_attached_row(
            accepted_row,
            tuple(draft.record for draft in sorted_drafts),
        )

    total_before = _count_rows(categories)
    total_after = _count_rows(new_categories)
    if total_before != total_after:
        raise RuntimeError("recovery transaction changed total row cardinality")
    accepted_before = _count_rows(categories, ACCEPTED_CATEGORY_INDEXES)
    accepted_after = _count_rows(new_categories, ACCEPTED_CATEGORY_INDEXES)
    admissions = sum(
        draft.record.action == ACTION_ADMIT_POSSIBLE
        for draft in plan.drafts
    )
    replacements = sum(
        draft.record.action == ACTION_REPLACE_POSSIBLE
        for draft in plan.drafts
    )
    if accepted_after != accepted_before + admissions:
        raise RuntimeError("accepted cardinality changed unexpectedly")
    stats = SourceIntactRecoveryStats(
        total_rows_before=total_before,
        total_rows_after=total_after,
        reviewed=len(plan.drafts),
        admitted_possible=admissions,
        replaced_possible=replacements,
        alternatives_attached=sum(
            draft.record.action == ACTION_ATTACH_ALTERNATIVE
            for draft in plan.drafts
        ),
        abstained=sum(
            draft.record.action == ACTION_ABSTAIN for draft in plan.drafts
        ),
        rejected=sum(
            draft.record.action == ACTION_REJECT for draft in plan.drafts
        ),
        accepted_rows_before=accepted_before,
        accepted_rows_after=accepted_after,
    )
    records = plan.records
    _json_copy(plan.as_dict(), "source-intact recovery plan")
    _json_copy(stats.as_dict(), "source-intact recovery stats")
    return new_categories, records, stats


__all__ = (
    "ACTION_ABSTAIN",
    "ACTION_ADMIT_POSSIBLE",
    "ACTION_ATTACH_ALTERNATIVE",
    "ACTION_REJECT",
    "ACTION_REPLACE_POSSIBLE",
    "SourceBackedRowFamilyEvidence",
    "SourceIntactAlternativeRelationship",
    "SourceIntactRecoveryCandidateContext",
    "SourceIntactRecoveryCostAudit",
    "SourceIntactRecoveryEvidence",
    "SourceIntactRecoveryGeometry",
    "SourceIntactRecoveryLocusContext",
    "SourceIntactRecoveryPlan",
    "SourceIntactRecoveryRecord",
    "SourceIntactRecoveryStats",
    "apply_source_intact_recovery_plan",
    "plan_source_intact_recoveries",
    "validate_source_intact_alternative_relationship",
)
