"""Pure arbitration between an accepted call and source-intact evidence.

The policy in this module is intentionally narrow.  It does not compare
classifier scores, prefer a longer path, mutate categories, or certify exact
boundaries.  Both candidates are reconstructed against the complete original
sequence before any interval relation is calculated.

Only one conflict topology can retire an accepted ``Possible`` call: a single
repeat/repeat collision confined to the two facing terminal repeats.  The
source-intact candidate must pass the strict existence gate, the configured
structural filter, a statistically assessed terminal-dominance test, and the
coding-tandem blacklist.  The accepted call must independently contradict both
spacer geometry and at least one other calibrated architecture domain.  Small
repeat count and insufficient terminal sample size never authorize retirement.

Containment, multi-repeat overlap, and stale evidence all fail closed.  An
eligible source-intact path always retains ``AMBIGUOUS`` boundary status.
"""

from dataclasses import asdict, dataclass
import hashlib
from typing import Optional, Tuple

from components.accepted_family_tandem_validator import (
    AcceptedFamilyTandemDecision,
)
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
    SourceIntactExistenceDecision,
    SourceIntactExistencePolicy,
    evaluate_source_intact_existence,
)


RULE_NAME = "accepted-vs-source-intact-conflict"
RULE_VERSION = "1.0.0"
COORDINATE_SYSTEM = "zero_based_half_open"

ACCEPTED_RETAIN = "retain"
ACCEPTED_RETIRE = "retire_conflicting_possible"
ACCEPTED_TANDEM_FIRST = "apply_tandem_demotion_first"
ACCEPTED_REVIEW_INVALID = "review_invalid_source_claim"

SOURCE_DO_NOT_ADMIT = "do_not_admit"
SOURCE_ADMIT_POSSIBLE = "admit_as_possible"
SOURCE_ATTACH_ALTERNATIVE = "attach_as_ambiguous_locus_alternative"


@dataclass(frozen=True)
class AcceptedSourceConflictPolicy:
    """Versioned, explicit authorization rule."""

    name: str = RULE_NAME
    version: str = RULE_VERSION
    replaceable_accepted_categories: Tuple[str, ...] = ("Possible",)
    required_primary_weakness_domain: str = "spacer_geometry"
    required_secondary_weakness_domains: Tuple[str, ...] = (
        "repeat_geometry",
        "repeat_coherence",
        "spacer_independence",
    )
    require_source_terminal_sample_assessed: bool = True

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("policy name must be a non-empty string")
        if not isinstance(self.version, str) or not self.version:
            raise ValueError("policy version must be a non-empty string")
        if (
                not isinstance(self.replaceable_accepted_categories, tuple)
                or not self.replaceable_accepted_categories
                or any(
                    not isinstance(value, str) or not value
                    for value in self.replaceable_accepted_categories
                )):
            raise ValueError(
                "replaceable_accepted_categories must be non-empty strings"
            )
        if (
                not isinstance(self.required_primary_weakness_domain, str)
                or not self.required_primary_weakness_domain):
            raise ValueError(
                "required_primary_weakness_domain must be non-empty"
            )
        if (
                not isinstance(self.required_secondary_weakness_domains, tuple)
                or not self.required_secondary_weakness_domains
                or any(
                    not isinstance(value, str) or not value
                    for value in self.required_secondary_weakness_domains
                )):
            raise ValueError(
                "required_secondary_weakness_domains must be non-empty strings"
            )
        if not isinstance(
                self.require_source_terminal_sample_assessed, bool):
            raise ValueError(
                "require_source_terminal_sample_assessed must be boolean"
            )

    def as_dict(self):
        return asdict(self)


ACCEPTED_SOURCE_CONFLICT_POLICY = AcceptedSourceConflictPolicy()


@dataclass(frozen=True)
class CandidateSourceSnapshot:
    """Immutable geometry re-sliced from the original sequence."""

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


@dataclass(frozen=True)
class ConflictArchitectureEvidence:
    """Source-derived strict and legacy structural evidence."""

    strict_existence_eligible: bool
    strict_existence_rejection_codes: Tuple[str, ...]
    weakness_domains: Tuple[str, ...]
    structural_filter_passed: bool
    structural_filter_terminal_status: Optional[str]
    strict_existence: SourceIntactExistenceDecision
    structural_filter: AdvancedFuzzySearchFilterEvidence

    def as_dict(self):
        return {
            "strict_existence_eligible": self.strict_existence_eligible,
            "strict_existence_rejection_codes": list(
                self.strict_existence_rejection_codes
            ),
            "weakness_domains": list(self.weakness_domains),
            "structural_filter_passed": self.structural_filter_passed,
            "structural_filter_terminal_status": (
                self.structural_filter_terminal_status
            ),
            "strict_existence": self.strict_existence.as_dict(),
            "structural_filter": self.structural_filter.as_dict(),
        }


@dataclass(frozen=True)
class ConflictRelationshipEvidence:
    """Topology calculated only after both paths reconstruct from source."""

    topology: str
    overlap_interval: Optional[Tuple[int, int]]
    overlap_bp: int
    accepted_overlap_fraction: float
    source_intact_overlap_fraction: float
    exact_shared_repeat_intervals: Tuple[Tuple[int, int], ...]
    repeat_overlap_pairs: Tuple[Tuple[int, int, int, int], ...]
    accepted_exact_repeat_subpath: bool
    source_intact_exact_repeat_subpath: bool
    single_facing_terminal_repeat_collision: bool

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class AcceptedSourceConflictDecision:
    """Non-mutating arbitration result and complete evidence record."""

    decision: str
    reason: str
    source_authorized: bool
    accepted_action: str
    source_intact_action: str
    boundary_status: str
    accepted_category: str
    source_sha256: Optional[str]
    accepted_snapshot: Optional[CandidateSourceSnapshot]
    source_intact_snapshot: Optional[CandidateSourceSnapshot]
    accepted_architecture: Optional[ConflictArchitectureEvidence]
    source_intact_architecture: Optional[ConflictArchitectureEvidence]
    relationship: Optional[ConflictRelationshipEvidence]
    accepted_tandem: Optional[AcceptedFamilyTandemDecision]
    source_intact_tandem: Optional[AcceptedFamilyTandemDecision]
    invalid_evidence_codes: Tuple[str, ...]
    policy: AcceptedSourceConflictPolicy

    def as_dict(self):
        return {
            "decision": self.decision,
            "reason": self.reason,
            "source_authorized": self.source_authorized,
            "accepted_action": self.accepted_action,
            "source_intact_action": self.source_intact_action,
            "boundary_status": self.boundary_status,
            "accepted_category": self.accepted_category,
            "source_sha256": self.source_sha256,
            "accepted_snapshot": (
                None
                if self.accepted_snapshot is None
                else self.accepted_snapshot.as_dict()
            ),
            "source_intact_snapshot": (
                None
                if self.source_intact_snapshot is None
                else self.source_intact_snapshot.as_dict()
            ),
            "accepted_architecture": (
                None
                if self.accepted_architecture is None
                else self.accepted_architecture.as_dict()
            ),
            "source_intact_architecture": (
                None
                if self.source_intact_architecture is None
                else self.source_intact_architecture.as_dict()
            ),
            "relationship": (
                None
                if self.relationship is None
                else self.relationship.as_dict()
            ),
            "accepted_tandem": (
                None
                if self.accepted_tandem is None
                else self.accepted_tandem.as_dict()
            ),
            "source_intact_tandem": (
                None
                if self.source_intact_tandem is None
                else self.source_intact_tandem.as_dict()
            ),
            "invalid_evidence_codes": list(self.invalid_evidence_codes),
            "policy": self.policy.as_dict(),
            "coordinate_system": COORDINATE_SYSTEM,
            "classifier_score_compared": False,
            "candidate_length_preference_used": False,
            "exact_boundary_certified": False,
        }


_WEAKNESS_DOMAIN_BY_CODE = {
    "invalid_consensus": "repeat_geometry",
    "consensus_length_out_of_range": "repeat_geometry",
    "repeat_length_out_of_range": "repeat_geometry",
    "repeat_length_cv_above_maximum": "repeat_geometry",
    "mean_repeat_identity_below_minimum": "repeat_coherence",
    "minimum_repeat_identity_below_minimum": "repeat_coherence",
    "spacer_length_out_of_range": "spacer_geometry",
    "spacer_length_cv_above_maximum": "spacer_geometry",
    "exact_spacer_duplicate_rule_failed": "spacer_independence",
    "distinct_spacer_identity_above_maximum": "spacer_independence",
    "prefix_kmer_dominance_above_maximum": "spacer_terminal_composition",
    "suffix_kmer_dominance_above_maximum": "spacer_terminal_composition",
    "repeat_count_below_minimum": "repeat_cardinality",
}


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _weakness_domains(existence, structural):
    domains = {
        _WEAKNESS_DOMAIN_BY_CODE[code]
        for code in existence.rejection_codes
        if code in _WEAKNESS_DOMAIN_BY_CODE
    }
    if not structural.passes:
        domains.add("legacy_structural_filter")
    return tuple(sorted(domains))


def _architecture_evidence(existence, structural):
    terminal = structural.terminal_dominance
    return ConflictArchitectureEvidence(
        strict_existence_eligible=existence.eligible,
        strict_existence_rejection_codes=existence.rejection_codes,
        weakness_domains=_weakness_domains(existence, structural),
        structural_filter_passed=structural.passes,
        structural_filter_terminal_status=(
            None if terminal is None else terminal.status
        ),
        strict_existence=existence,
        structural_filter=structural,
    )


def _snapshot(source, candidate, reconstruction):
    starts = tuple(reconstruction.genomic_repeat_starts)
    repeats = tuple(candidate.list_repeats)
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    spacer_intervals = tuple(
        (left_end, right_start)
        for (_left_start, left_end), (right_start, _right_end)
        in zip(repeat_intervals, repeat_intervals[1:])
    )
    interval = (repeat_intervals[0][0], repeat_intervals[-1][1])
    return CandidateSourceSnapshot(
        interval=interval,
        repeat_intervals=repeat_intervals,
        spacer_intervals=spacer_intervals,
        repeat_count=len(repeat_intervals),
        spacer_count=len(spacer_intervals),
        source_slice_sha256=_sha256(source[interval[0]:interval[1]]),
        repeat_source_slice_sha256s=tuple(
            _sha256(source[start:end]) for start, end in repeat_intervals
        ),
        spacer_source_slice_sha256s=tuple(
            _sha256(source[start:end]) for start, end in spacer_intervals
        ),
    )


def _interval_overlap(first, second):
    start = max(first[0], second[0])
    end = min(first[1], second[1])
    return None if end <= start else (start, end)


def _topology(accepted, source_intact):
    if accepted == source_intact:
        return "coincident"
    if accepted[1] <= source_intact[0]:
        return (
            "touching"
            if accepted[1] == source_intact[0]
            else "accepted_upstream_disjoint"
        )
    if source_intact[1] <= accepted[0]:
        return (
            "touching"
            if source_intact[1] == accepted[0]
            else "source_intact_upstream_disjoint"
        )
    if (
            source_intact[0] <= accepted[0]
            and accepted[1] <= source_intact[1]):
        return "source_intact_contains_accepted"
    if (
            accepted[0] <= source_intact[0]
            and source_intact[1] <= accepted[1]):
        return "accepted_contains_source_intact"
    if source_intact[0] < accepted[0]:
        return "source_intact_left_edge_overlap"
    return "accepted_left_edge_overlap"


def _relationship(accepted, source_intact):
    accepted_interval = accepted.interval
    source_interval = source_intact.interval
    overlap = _interval_overlap(accepted_interval, source_interval)
    overlap_bp = 0 if overlap is None else overlap[1] - overlap[0]

    accepted_repeats = set(accepted.repeat_intervals)
    source_repeats = set(source_intact.repeat_intervals)
    exact_shared = tuple(sorted(accepted_repeats.intersection(source_repeats)))
    repeat_overlap_pairs = tuple(
        (
            source_index,
            accepted_index,
            collision[0],
            collision[1],
        )
        for source_index, source_repeat in enumerate(
            source_intact.repeat_intervals
        )
        for accepted_index, accepted_repeat in enumerate(
            accepted.repeat_intervals
        )
        for collision in (_interval_overlap(source_repeat, accepted_repeat),)
        if collision is not None
    )
    source_is_left = source_interval[0] < accepted_interval[0]
    accepted_is_left = accepted_interval[0] < source_interval[0]
    single_terminal_collision = False
    if overlap is not None and len(repeat_overlap_pairs) == 1:
        source_index, accepted_index, start, end = repeat_overlap_pairs[0]
        collision_is_whole_envelope_overlap = (start, end) == overlap
        if source_is_left:
            terminal_indexes_match = (
                source_index == len(source_intact.repeat_intervals) - 1
                and accepted_index == 0
            )
        elif accepted_is_left:
            terminal_indexes_match = (
                accepted_index == len(accepted.repeat_intervals) - 1
                and source_index == 0
            )
        else:
            terminal_indexes_match = False
        single_terminal_collision = (
            collision_is_whole_envelope_overlap and terminal_indexes_match
        )

    return ConflictRelationshipEvidence(
        topology=_topology(accepted_interval, source_interval),
        overlap_interval=overlap,
        overlap_bp=overlap_bp,
        accepted_overlap_fraction=(
            overlap_bp / float(accepted_interval[1] - accepted_interval[0])
        ),
        source_intact_overlap_fraction=(
            overlap_bp / float(source_interval[1] - source_interval[0])
        ),
        exact_shared_repeat_intervals=exact_shared,
        repeat_overlap_pairs=repeat_overlap_pairs,
        accepted_exact_repeat_subpath=(
            bool(accepted_repeats)
            and accepted_repeats.issubset(source_repeats)
        ),
        source_intact_exact_repeat_subpath=(
            bool(source_repeats)
            and source_repeats.issubset(accepted_repeats)
        ),
        single_facing_terminal_repeat_collision=single_terminal_collision,
    )


def _tandem_binding_issue(source, snapshot, decision, label):
    if not isinstance(decision, AcceptedFamilyTandemDecision):
        return "{}_tandem_evidence_missing_or_wrong_type".format(label)
    if decision.accepted_interval != snapshot.interval:
        return "{}_tandem_interval_mismatch".format(label)
    if decision.source_sha256 != _sha256(source):
        return "{}_tandem_source_hash_mismatch".format(label)
    if decision.accepted_source_slice_sha256 != snapshot.source_slice_sha256:
        return "{}_tandem_slice_hash_mismatch".format(label)
    return None


def _base_decision(
        *,
        decision,
        reason,
        source_authorized,
        accepted_action,
        source_intact_action,
        accepted_category,
        source_sha256,
        accepted_snapshot,
        source_intact_snapshot,
        accepted_architecture,
        source_intact_architecture,
        relationship,
        accepted_tandem,
        source_intact_tandem,
        invalid_evidence_codes,
        policy):
    return AcceptedSourceConflictDecision(
        decision=decision,
        reason=reason,
        source_authorized=source_authorized,
        accepted_action=accepted_action,
        source_intact_action=source_intact_action,
        boundary_status=BOUNDARY_AMBIGUOUS,
        accepted_category=accepted_category,
        source_sha256=source_sha256,
        accepted_snapshot=accepted_snapshot,
        source_intact_snapshot=source_intact_snapshot,
        accepted_architecture=accepted_architecture,
        source_intact_architecture=source_intact_architecture,
        relationship=relationship,
        accepted_tandem=accepted_tandem,
        source_intact_tandem=source_intact_tandem,
        invalid_evidence_codes=tuple(invalid_evidence_codes),
        policy=policy,
    )


def arbitrate_accepted_source_intact_conflict(
    source_sequence,
    accepted_candidate,
    source_intact_candidate,
    *,
    accepted_category,
    structural_filter,
    accepted_tandem_decision,
    source_intact_tandem_decision,
    policy=ACCEPTED_SOURCE_CONFLICT_POLICY,
    existence_policy=SOURCE_INTACT_EXISTENCE_POLICY,
):
    """Return a source-bound, non-mutating accepted-conflict decision.

    The caller must provide direct tandem decisions for both paths.  Those
    records are accepted only when their interval and source hashes bind to
    the candidates supplied in this call.
    """

    if not isinstance(policy, AcceptedSourceConflictPolicy):
        raise TypeError("policy must be an AcceptedSourceConflictPolicy")
    if not isinstance(existence_policy, SourceIntactExistencePolicy):
        raise TypeError("existence_policy must be a SourceIntactExistencePolicy")
    if not isinstance(accepted_category, str) or not accepted_category:
        raise ValueError("accepted_category must be a non-empty string")
    if not isinstance(structural_filter, AdvancedFuzzySearchFilter):
        raise TypeError(
            "structural_filter must be an AdvancedFuzzySearchFilter"
        )

    source_valid = (
        isinstance(source_sequence, str)
        and bool(source_sequence)
        and not set(source_sequence).difference(IUPAC_DNA_ALPHABET)
    )
    source = source_sequence.upper() if source_valid else source_sequence
    source_hash = _sha256(source) if source_valid else None

    accepted_existence = evaluate_source_intact_existence(
        source_sequence,
        accepted_candidate,
        policy=existence_policy,
    )
    source_existence = evaluate_source_intact_existence(
        source_sequence,
        source_intact_candidate,
        policy=existence_policy,
    )
    invalid_codes = []
    if not source_valid:
        invalid_codes.append("invalid_source_sequence")
        return _base_decision(
            decision="invalid_evidence",
            reason="invalid_source_sequence",
            source_authorized=False,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_DO_NOT_ADMIT,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=None,
            source_intact_snapshot=None,
            accepted_architecture=None,
            source_intact_architecture=None,
            relationship=None,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=invalid_codes,
            policy=policy,
        )

    try:
        accepted_reconstruction = validate_candidate_reconstruction(
            source, accepted_candidate
        )
    except (AttributeError, TypeError, ValueError):
        accepted_reconstruction = None
    try:
        source_reconstruction = validate_candidate_reconstruction(
            source, source_intact_candidate
        )
    except (AttributeError, TypeError, ValueError):
        source_reconstruction = None

    if accepted_reconstruction is None or not accepted_reconstruction.valid:
        invalid_codes.append("accepted_source_reconstruction_failed")
    if source_reconstruction is None or not source_reconstruction.valid:
        invalid_codes.append("source_intact_source_reconstruction_failed")
    if invalid_codes:
        return _base_decision(
            decision="invalid_evidence",
            reason="candidate_source_reconstruction_failed",
            source_authorized=False,
            accepted_action=(
                ACCEPTED_REVIEW_INVALID
                if "accepted_source_reconstruction_failed" in invalid_codes
                else ACCEPTED_RETAIN
            ),
            source_intact_action=SOURCE_DO_NOT_ADMIT,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=None,
            source_intact_snapshot=None,
            accepted_architecture=None,
            source_intact_architecture=None,
            relationship=None,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=invalid_codes,
            policy=policy,
        )

    accepted_snapshot = _snapshot(
        source, accepted_candidate, accepted_reconstruction
    )
    source_snapshot = _snapshot(
        source, source_intact_candidate, source_reconstruction
    )
    accepted_structural = structural_filter.evidence(accepted_candidate)
    source_structural = structural_filter.evidence(source_intact_candidate)
    accepted_architecture = _architecture_evidence(
        accepted_existence, accepted_structural
    )
    source_architecture = _architecture_evidence(
        source_existence, source_structural
    )
    relationship = _relationship(accepted_snapshot, source_snapshot)

    for label, snapshot, tandem in (
        ("accepted", accepted_snapshot, accepted_tandem_decision),
        (
            "source_intact",
            source_snapshot,
            source_intact_tandem_decision,
        ),
    ):
        issue = _tandem_binding_issue(source, snapshot, tandem, label)
        if issue is not None:
            invalid_codes.append(issue)
    if invalid_codes:
        return _base_decision(
            decision="invalid_evidence",
            reason="tandem_evidence_not_bound_to_current_source_path",
            source_authorized=False,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_DO_NOT_ADMIT,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=invalid_codes,
            policy=policy,
        )

    if accepted_tandem_decision.supported:
        return _base_decision(
            decision="defer_to_accepted_tandem_demotion",
            reason="accepted_path_has_direct_coding_tandem_support",
            source_authorized=False,
            accepted_action=ACCEPTED_TANDEM_FIRST,
            source_intact_action=SOURCE_DO_NOT_ADMIT,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    terminal_assessed = (
        source_architecture.structural_filter_terminal_status == "passed"
    )
    source_authorized = (
        source_architecture.strict_existence_eligible
        and source_architecture.structural_filter_passed
        and not source_intact_tandem_decision.supported
        and (
            terminal_assessed
            or not policy.require_source_terminal_sample_assessed
        )
    )
    if not source_authorized:
        if source_intact_tandem_decision.supported:
            reason = "source_intact_path_has_direct_coding_tandem_support"
        elif not source_architecture.strict_existence_eligible:
            reason = "source_intact_existence_gate_rejected"
        elif not source_architecture.structural_filter_passed:
            reason = "source_intact_structural_filter_rejected"
        else:
            reason = "source_intact_terminal_sample_not_assessed"
        return _base_decision(
            decision="source_intact_not_authorized",
            reason=reason,
            source_authorized=False,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_DO_NOT_ADMIT,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    if relationship.overlap_bp == 0:
        return _base_decision(
            decision="independent_source_intact_array",
            reason="source_valid_intervals_are_disjoint_or_touching",
            source_authorized=True,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_ADMIT_POSSIBLE,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    if (
            accepted_snapshot == source_snapshot
            or (
                relationship.accepted_exact_repeat_subpath
                and relationship.source_intact_exact_repeat_subpath
            )):
        return _base_decision(
            decision="duplicate_source_geometry",
            reason="accepted_and_source_intact_paths_are_identical",
            source_authorized=True,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_DO_NOT_ADMIT,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    if (
            relationship.accepted_exact_repeat_subpath
            or relationship.source_intact_exact_repeat_subpath):
        return _base_decision(
            decision="same_array_boundary_alternative",
            reason="one_source_repeat_path_is_an_exact_subpath_of_the_other",
            source_authorized=True,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_ATTACH_ALTERNATIVE,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    weakness_domains = set(accepted_architecture.weakness_domains)
    primary_failed = (
        policy.required_primary_weakness_domain in weakness_domains
    )
    secondary_failed = bool(
        weakness_domains.intersection(
            policy.required_secondary_weakness_domains
        )
    )
    accepted_has_independent_contradictions = (
        not accepted_architecture.strict_existence_eligible
        and primary_failed
        and secondary_failed
    )
    replaceable_category = (
        accepted_category in policy.replaceable_accepted_categories
    )
    if (
            relationship.single_facing_terminal_repeat_collision
            and replaceable_category
            and accepted_has_independent_contradictions):
        return _base_decision(
            decision="source_intact_dominates_terminal_conflict",
            reason=(
                "single_terminal_repeat_collision_and_independent_accepted_"
                "architecture_contradictions"
            ),
            source_authorized=True,
            accepted_action=ACCEPTED_RETIRE,
            source_intact_action=SOURCE_ADMIT_POSSIBLE,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    if relationship.single_facing_terminal_repeat_collision:
        return _base_decision(
            decision="plausible_adjacent_arrays_boundary_ambiguous",
            reason=(
                "single_terminal_repeat_collision_without_retirement_"
                "authorization"
            ),
            source_authorized=True,
            accepted_action=ACCEPTED_RETAIN,
            source_intact_action=SOURCE_ADMIT_POSSIBLE,
            accepted_category=accepted_category,
            source_sha256=source_hash,
            accepted_snapshot=accepted_snapshot,
            source_intact_snapshot=source_snapshot,
            accepted_architecture=accepted_architecture,
            source_intact_architecture=source_architecture,
            relationship=relationship,
            accepted_tandem=accepted_tandem_decision,
            source_intact_tandem=source_intact_tandem_decision,
            invalid_evidence_codes=(),
            policy=policy,
        )

    return _base_decision(
        decision="complex_overlap_abstained",
        reason="containment_or_complex_overlap_is_not_safe_to_arbitrate",
        source_authorized=True,
        accepted_action=ACCEPTED_RETAIN,
        source_intact_action=SOURCE_ATTACH_ALTERNATIVE,
        accepted_category=accepted_category,
        source_sha256=source_hash,
        accepted_snapshot=accepted_snapshot,
        source_intact_snapshot=source_snapshot,
        accepted_architecture=accepted_architecture,
        source_intact_architecture=source_architecture,
        relationship=relationship,
        accepted_tandem=accepted_tandem_decision,
        source_intact_tandem=source_intact_tandem_decision,
        invalid_evidence_codes=(),
        policy=policy,
    )


__all__ = [
    "ACCEPTED_RETAIN",
    "ACCEPTED_RETIRE",
    "ACCEPTED_SOURCE_CONFLICT_POLICY",
    "ACCEPTED_TANDEM_FIRST",
    "ACCEPTED_REVIEW_INVALID",
    "AcceptedSourceConflictDecision",
    "AcceptedSourceConflictPolicy",
    "CandidateSourceSnapshot",
    "ConflictArchitectureEvidence",
    "ConflictRelationshipEvidence",
    "SOURCE_ADMIT_POSSIBLE",
    "SOURCE_ATTACH_ALTERNATIVE",
    "SOURCE_DO_NOT_ADMIT",
    "arbitrate_accepted_source_intact_conflict",
]
