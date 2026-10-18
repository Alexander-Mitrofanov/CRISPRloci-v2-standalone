"""Pure planning for source-backed completion of an accepted CRISPR core.

This module deliberately does not mutate score rows, categories, candidates,
or pipeline state.  It separates two questions that the historical pipeline
coupled:

* an already accepted, strictly adjudicated core authorizes that a CRISPR
  locus exists;
* source-only terminal evidence determines whether its reported boundaries
  should be extended.

Unaccepted candidates may be supplied as geometry witnesses.  They can remove
an otherwise over-conservative collision only when the proposed terminal
repeat is one of their exact source intervals and their source topology is
compatible with the accepted core.  A witness can never originate a plan.

Coordinates are zero-based and end-exclusive.
"""

from dataclasses import dataclass
import hashlib
import math
from numbers import Real
import statistics

from components.architecture_rescue_policy import ELIGIBLE
from components.components_evaluated_arrays_enhancement import (
    apply_source_backed_terminal_repeat_proposals,
    apply_source_backed_weak_terminal_repeat_proposal,
    build_source_backed_candidate_from_boundaries,
    source_backed_terminal_repeat_proposals,
    source_backed_weak_terminal_repeat_proposals,
    validate_candidate_source_reconstruction,
)


ACCEPTED_CATEGORIES = frozenset(("Bona-fide", "Possible"))
MAXIMUM_WITNESS_GAP = 140
MINIMUM_CONSENSUS_CONTAINMENT_IDENTITY = 0.80
MINIMUM_CONSENSUS_LENGTH_RATIO = 2.0 / 3.0
MINIMUM_PERIOD_TOLERANCE = 3.0
PERIOD_TOLERANCE_FRACTION = 0.10


def _finite_unit_interval(value):
    return (
        isinstance(value, Real)
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and 0.0 <= float(value) <= 1.0
    )


def _repeat_path(candidate):
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts,
            candidate.list_repeats,
        )
    )


def _path_interval(path):
    return path[0][0], path[-1][1]


def _intervals_overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _path_is_contiguous_subpath(shorter, longer):
    shorter = tuple(shorter)
    longer = tuple(longer)
    if not shorter or len(shorter) > len(longer):
        return False
    return any(
        longer[offset : offset + len(shorter)] == shorter
        for offset in range(len(longer) - len(shorter) + 1)
    )


def _normalized_consensus(candidate):
    consensus = getattr(candidate, "consensus", None)
    if not isinstance(consensus, str) or not consensus:
        raise ValueError("every boundary candidate needs a consensus")
    return consensus.upper().replace(" ", "").replace("-", "").replace("+", "")


def _containment_identity(first, second):
    if not first or not second:
        return 0.0
    shorter, longer = (first, second) if len(first) <= len(second) else (second, first)
    return max(
        sum(
            left == right
            for left, right in zip(
                shorter,
                longer[offset : offset + len(shorter)],
            )
        )
        for offset in range(len(longer) - len(shorter) + 1)
    ) / float(len(shorter))


def _median_period(path):
    starts = tuple(start for start, _end in path)
    if len(starts) < 2:
        return None
    return float(
        statistics.median(right - left for left, right in zip(starts, starts[1:]))
    )


def _interval_gap(first, second):
    if first[1] < second[0]:
        return second[0] - first[1]
    if second[1] < first[0]:
        return first[0] - second[1]
    return 0


def _witness_topology_is_compatible(candidate, witness):
    candidate_path = _repeat_path(candidate)
    witness_path = witness.geometry.repeat_intervals
    if (
        _interval_gap(_path_interval(candidate_path), _path_interval(witness_path))
        > MAXIMUM_WITNESS_GAP
    ):
        return False

    candidate_consensus = _normalized_consensus(candidate)
    witness_consensus = witness.consensus
    length_ratio = min(len(candidate_consensus), len(witness_consensus)) / float(
        max(len(candidate_consensus), len(witness_consensus))
    )
    if length_ratio < MINIMUM_CONSENSUS_LENGTH_RATIO:
        return False
    if (
        _containment_identity(
            candidate_consensus,
            witness_consensus,
        )
        < MINIMUM_CONSENSUS_CONTAINMENT_IDENTITY
    ):
        return False

    candidate_period = _median_period(candidate_path)
    witness_period = _median_period(witness_path)
    if candidate_period is None or witness_period is None:
        return False
    tolerance = max(
        MINIMUM_PERIOD_TOLERANCE,
        PERIOD_TOLERANCE_FRACTION * max(candidate_period, witness_period),
    )
    if abs(candidate_period - witness_period) > tolerance:
        return False

    candidate_starts = tuple(start for start, _end in candidate_path)
    witness_starts = tuple(start for start, _end in witness_path)
    return any(
        abs((source_start + source_period) - observed) <= tolerance
        for source_starts, source_period, observed_starts in (
            (candidate_starts, candidate_period, witness_starts),
            (witness_starts, witness_period, candidate_starts),
        )
        for source_start in source_starts
        for observed in observed_starts
    )


@dataclass(frozen=True)
class SourceGeometry:
    """An immutable source reconstruction of one array geometry."""

    repeat_intervals: tuple
    repeat_sequences: tuple
    spacer_sequences: tuple

    @property
    def source_interval(self):
        return _path_interval(self.repeat_intervals)

    @property
    def repeat_count(self):
        return len(self.repeat_intervals)

    def materialize(self, full_dna):
        """Build a fresh candidate after checking the exact source geometry."""

        candidate = build_source_backed_candidate_from_boundaries(
            full_dna,
            None,
            tuple(start for start, _end in self.repeat_intervals),
            tuple(end for _start, end in self.repeat_intervals),
        )
        observed = SourceGeometry.from_candidate(full_dna, candidate)
        if observed != self:
            raise ValueError("source sequence does not match the planned geometry")
        return candidate

    @classmethod
    def from_candidate(cls, full_dna, candidate):
        validate_candidate_source_reconstruction(
            full_dna,
            candidate.list_repeats,
            candidate.list_repeat_starts,
            candidate.list_spacers,
        )
        return cls(
            repeat_intervals=_repeat_path(candidate),
            repeat_sequences=tuple(candidate.list_repeats),
            spacer_sequences=tuple(candidate.list_spacers),
        )


@dataclass(frozen=True)
class CoreAuthorizationEvidence:
    """Explicit fail-closed evidence attached to the accepted core."""

    coding_context_assessed: bool
    tandem_evidence_assessed: bool
    coding_tandem_clear: bool

    def __post_init__(self):
        for field_name in (
            "coding_context_assessed",
            "tandem_evidence_assessed",
            "coding_tandem_clear",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ValueError("{} must be a boolean".format(field_name))


@dataclass(frozen=True)
class FinalBoundaryEvidence:
    """Fresh final-geometry evidence supplied by the integration caller."""

    fresh_score: float
    source_reconstruction_valid: bool
    structural_valid: bool
    phase_unambiguous: bool
    coding_context_assessed: bool
    tandem_evidence_assessed: bool
    coding_tandem_clear: bool
    reason_codes: tuple = ()

    def __post_init__(self):
        if not _finite_unit_interval(self.fresh_score):
            raise ValueError("fresh_score must be finite and in [0, 1]")
        for field_name in (
            "source_reconstruction_valid",
            "structural_valid",
            "phase_unambiguous",
            "coding_context_assessed",
            "tandem_evidence_assessed",
            "coding_tandem_clear",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ValueError("{} must be a boolean".format(field_name))
        if not isinstance(self.reason_codes, tuple) or not all(
            isinstance(reason, str) and reason for reason in self.reason_codes
        ):
            raise ValueError("reason_codes must be a tuple of non-empty strings")

    @property
    def passes(self):
        return all(
            (
                self.source_reconstruction_valid,
                self.structural_valid,
                self.phase_unambiguous,
                self.coding_context_assessed,
                self.tandem_evidence_assessed,
                self.coding_tandem_clear,
            )
        )


@dataclass(frozen=True)
class BoundaryGeometryWitness:
    """A source-family path that may corroborate but never authorize."""

    candidate: object
    path_id: str
    accepted: bool = False

    def __post_init__(self):
        if not isinstance(self.path_id, str) or not self.path_id:
            raise ValueError("a geometry witness needs a path_id")
        if not isinstance(self.accepted, bool):
            raise ValueError("witness accepted status must be boolean")


@dataclass(frozen=True)
class WitnessContextEvidence:
    """Attest that collision paths and their acceptance states are complete."""

    source_family_complete: bool
    accepted_statuses_complete: bool

    def __post_init__(self):
        for field_name in (
            "source_family_complete",
            "accepted_statuses_complete",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ValueError("{} must be a boolean".format(field_name))

    @property
    def complete(self):
        return self.source_family_complete and self.accepted_statuses_complete


@dataclass(frozen=True)
class BoundaryCompletionStep:
    """One source-authorized monotonic boundary addition."""

    tier: str
    side: str
    repeat_interval: tuple
    spacer_interval: tuple
    witness_path_ids: tuple
    evidence: tuple


@dataclass(frozen=True)
class AcceptedCoreBoundaryCompletionPlan:
    """Immutable output of boundary-completion planning."""

    action: str
    reason_codes: tuple
    core_category: str
    core_score: float
    final_score: object
    core_policy_name: str
    core_policy_version: str
    source_sha256: str
    core_geometry: SourceGeometry
    proposed_geometry: SourceGeometry
    steps: tuple
    final_evidence: object = None

    @property
    def changed(self):
        return self.core_geometry != self.proposed_geometry

    def materialize_proposed_candidate(self, full_dna):
        observed_hash = hashlib.sha256(full_dna.upper().encode("ascii")).hexdigest()
        if observed_hash != self.source_sha256:
            raise ValueError("source sequence changed after planning")
        return self.proposed_geometry.materialize(full_dna)


@dataclass(frozen=True)
class _NormalizedWitness:
    geometry: SourceGeometry
    consensus: str
    path_id: str
    accepted: bool


def _stable_witness_id(geometry):
    payload = repr(geometry.repeat_intervals).encode("ascii")
    return "geometry-" + hashlib.sha256(payload).hexdigest()[:16]


def _normalize_witness(full_dna, value):
    if isinstance(value, BoundaryGeometryWitness):
        candidate = value.candidate
        path_id = value.path_id
        accepted = value.accepted
    elif isinstance(value, dict):
        candidate = value.get("candidate")
        path_id = value.get("path_id")
        accepted = value.get("accepted", False)
        if not isinstance(path_id, str) or not path_id:
            raise ValueError("a lineage witness needs a path_id")
        if not isinstance(accepted, bool):
            raise ValueError("lineage accepted status must be boolean")
    else:
        candidate = value
        geometry = SourceGeometry.from_candidate(full_dna, candidate)
        path_id = _stable_witness_id(geometry)
        accepted = False
    geometry = SourceGeometry.from_candidate(full_dna, candidate)
    source_materialized = geometry.materialize(full_dna)
    return _NormalizedWitness(
        geometry=geometry,
        consensus=_normalized_consensus(source_materialized),
        path_id=path_id,
        accepted=accepted,
    )


def _proposal_key(proposal):
    if not isinstance(proposal, dict):
        raise ValueError("terminal proposal is malformed")
    interval = proposal.get("repeat_interval")
    if (
        proposal.get("side") not in ("left", "right")
        or not isinstance(interval, (tuple, list))
        or len(interval) != 2
        or any(
            not isinstance(value, int) or isinstance(value, bool) for value in interval
        )
        or interval[0] < 0
        or interval[1] <= interval[0]
    ):
        raise ValueError("terminal proposal is malformed")
    return proposal["side"], tuple(interval)


def _unique_side_proposals(proposals, closed_sides):
    by_side = {"left": {}, "right": {}}
    for proposal in proposals:
        side, interval = _proposal_key(proposal)
        if side in closed_sides:
            continue
        by_side[side].setdefault(interval, proposal)
    ambiguous = tuple(side for side in ("left", "right") if len(by_side[side]) > 1)
    selected = tuple(
        by_side[side][min(by_side[side])] for side in ("left", "right") if by_side[side]
    )
    return selected, ambiguous


def _added_span(candidate, proposal):
    path = _repeat_path(candidate)
    interval = tuple(proposal["repeat_interval"])
    if proposal["side"] == "left":
        return interval[0], path[0][0]
    return path[-1][1], interval[1]


def _proposal_spacer_interval(candidate, proposal):
    path = _repeat_path(candidate)
    interval = tuple(proposal["repeat_interval"])
    if proposal["side"] == "left":
        return interval[1], path[0][0]
    return path[-1][1], interval[0]


def _witness_has_next_exact_repeat(candidate, proposal, witness):
    current_path = _repeat_path(candidate)
    witness_path = witness.geometry.repeat_intervals
    proposal_interval = tuple(proposal["repeat_interval"])
    indexes = tuple(
        index
        for index, interval in enumerate(witness_path)
        if interval == proposal_interval
    )
    if len(indexes) != 1:
        return False
    index = indexes[0]
    if proposal["side"] == "right":
        inner_path = witness_path[:index]
        return not inner_path or (
            len(inner_path) <= len(current_path)
            and current_path[-len(inner_path) :] == inner_path
        )
    inner_path = witness_path[index + 1 :]
    return not inner_path or (
        len(inner_path) <= len(current_path)
        and current_path[: len(inner_path)] == inner_path
    )


def _proposal_collision_evidence(candidate, proposal, witnesses):
    added_span = _added_span(candidate, proposal)
    proposal_interval = tuple(proposal["repeat_interval"])
    new_path = (
        (proposal_interval,) + _repeat_path(candidate)
        if proposal["side"] == "left"
        else _repeat_path(candidate) + (proposal_interval,)
    )
    bridge_ids = []
    for witness in witnesses:
        witness_interval = witness.geometry.source_interval
        if not _intervals_overlap(added_span, witness_interval):
            continue
        if witness.accepted:
            if not _path_is_contiguous_subpath(
                witness.geometry.repeat_intervals, new_path
            ):
                return (), "overlapping_non_subpath_accepted_witness"
            continue
        if not (
            _witness_has_next_exact_repeat(candidate, proposal, witness)
            and _witness_topology_is_compatible(candidate, witness)
        ):
            return (), "incompatible_overlapping_geometry_witness"
        bridge_ids.append(witness.path_id)
    return tuple(sorted(set(bridge_ids))), None


def _step_evidence(proposal):
    excluded = frozenset(("side", "repeat_interval"))
    return tuple(
        sorted(
            (str(key), value)
            for key, value in proposal.items()
            if key not in excluded
            and isinstance(value, (str, int, float, bool, type(None)))
        )
    )


def _policy_identity(core_decision):
    policy = getattr(core_decision, "policy", None)
    name = getattr(policy, "name", "unspecified")
    version = getattr(policy, "version", "unspecified")
    if not isinstance(name, str) or not name:
        name = "unspecified"
    if not isinstance(version, str) or not version:
        version = "unspecified"
    return name, version


def _plan(
    action,
    reasons,
    category,
    core_score,
    core_decision,
    source_hash,
    core_geometry,
    proposed_geometry,
    steps=(),
    final_evidence=None,
):
    policy_name, policy_version = _policy_identity(core_decision)
    return AcceptedCoreBoundaryCompletionPlan(
        action=action,
        reason_codes=tuple(reasons),
        core_category=category,
        core_score=float(core_score),
        final_score=(
            None if final_evidence is None else float(final_evidence.fresh_score)
        ),
        core_policy_name=policy_name,
        core_policy_version=policy_version,
        source_sha256=source_hash,
        core_geometry=core_geometry,
        proposed_geometry=proposed_geometry,
        steps=tuple(steps),
        final_evidence=final_evidence,
    )


def plan_accepted_core_boundary_completion(
    full_dna,
    core_row,
    core_category,
    core_decision,
    core_authorization,
    witnesses=(),
    witness_context=None,
    final_evidence_provider=None,
    strong_proposal_provider=source_backed_terminal_repeat_proposals,
    weak_proposal_provider=source_backed_weak_terminal_repeat_proposals,
):
    """Return a mutation-free plan for maximal accepted-core completion.

    The final evidence provider receives a newly reconstructed candidate and
    must return :class:`FinalBoundaryEvidence`.  Its model score is recorded
    but is intentionally not an acceptance gate: array existence was already
    authorized by the strict core.  Every other final gate fails closed.
    """

    if not isinstance(full_dna, str):
        raise ValueError("full_dna must be a string")
    if (
        not isinstance(core_row, (tuple, list))
        or len(core_row) < 2
        or not _finite_unit_interval(core_row[0])
    ):
        raise ValueError("core_row must contain a valid score and candidate")
    if core_category not in ACCEPTED_CATEGORIES:
        raise ValueError("core_category must be an accepted category")
    if not isinstance(core_authorization, CoreAuthorizationEvidence):
        raise TypeError("core_authorization must be CoreAuthorizationEvidence")
    if isinstance(witnesses, (str, bytes)):
        raise ValueError("witnesses must be an iterable")
    if not callable(strong_proposal_provider):
        raise TypeError("strong_proposal_provider must be callable")
    if not callable(weak_proposal_provider):
        raise TypeError("weak_proposal_provider must be callable")
    if final_evidence_provider is not None and not callable(final_evidence_provider):
        raise TypeError("final_evidence_provider must be callable")

    core_candidate = core_row[1]
    core_geometry = SourceGeometry.from_candidate(full_dna, core_candidate)
    source_hash = hashlib.sha256(full_dna.upper().encode("ascii")).hexdigest()
    core_eligible = getattr(core_decision, "decision", None) == ELIGIBLE
    if not (
        core_eligible
        and core_authorization.coding_context_assessed
        and core_authorization.tandem_evidence_assessed
        and core_authorization.coding_tandem_clear
    ):
        return _plan(
            "abstain",
            ("core_not_strictly_eligible",),
            core_category,
            core_row[0],
            core_decision,
            source_hash,
            core_geometry,
            core_geometry,
        )
    if (
        not isinstance(witness_context, WitnessContextEvidence)
        or not witness_context.complete
    ):
        return _plan(
            "abstain",
            ("witness_context_incomplete",),
            core_category,
            core_row[0],
            core_decision,
            source_hash,
            core_geometry,
            core_geometry,
        )

    normalized_witnesses = tuple(
        sorted(
            (_normalize_witness(full_dna, witness) for witness in witnesses),
            key=lambda witness: (
                witness.geometry.repeat_intervals,
                witness.path_id,
                witness.accepted,
            ),
        )
    )
    if len({witness.path_id for witness in normalized_witnesses}) != len(
        normalized_witnesses
    ):
        raise ValueError("witness path IDs must be unique")

    current = core_geometry.materialize(full_dna)
    observed_geometries = {core_geometry.repeat_intervals}
    steps = []
    closed_sides = set()
    weak_used = False

    while True:
        strong, ambiguous = _unique_side_proposals(
            tuple(strong_proposal_provider(full_dna, current)),
            closed_sides,
        )
        if ambiguous:
            return _plan(
                "abstain",
                ("ambiguous_strong_terminal_{}".format("_".join(ambiguous)),),
                core_category,
                core_row[0],
                core_decision,
                source_hash,
                core_geometry,
                SourceGeometry.from_candidate(full_dna, current),
                steps,
            )
        if strong:
            bridge_by_proposal = {}
            for proposal in strong:
                bridge_ids, collision = _proposal_collision_evidence(
                    current, proposal, normalized_witnesses
                )
                if collision is not None:
                    return _plan(
                        "abstain",
                        (collision,),
                        core_category,
                        core_row[0],
                        core_decision,
                        source_hash,
                        core_geometry,
                        SourceGeometry.from_candidate(full_dna, current),
                        steps,
                    )
                bridge_by_proposal[_proposal_key(proposal)] = bridge_ids
            updated = apply_source_backed_terminal_repeat_proposals(
                full_dna, current, strong
            )
            for proposal in strong:
                steps.append(
                    BoundaryCompletionStep(
                        tier="strong",
                        side=proposal["side"],
                        repeat_interval=tuple(proposal["repeat_interval"]),
                        spacer_interval=_proposal_spacer_interval(current, proposal),
                        witness_path_ids=bridge_by_proposal[_proposal_key(proposal)],
                        evidence=_step_evidence(proposal),
                    )
                )
        elif not weak_used:
            weak, weak_ambiguous = _unique_side_proposals(
                tuple(weak_proposal_provider(full_dna, current)),
                closed_sides,
            )
            if weak_ambiguous or len(weak) > 1:
                return _plan(
                    "abstain",
                    ("ambiguous_weak_terminal",),
                    core_category,
                    core_row[0],
                    core_decision,
                    source_hash,
                    core_geometry,
                    SourceGeometry.from_candidate(full_dna, current),
                    steps,
                )
            if not weak:
                break
            proposal = weak[0]
            bridge_ids, collision = _proposal_collision_evidence(
                current, proposal, normalized_witnesses
            )
            if collision is not None:
                return _plan(
                    "abstain",
                    (collision,),
                    core_category,
                    core_row[0],
                    core_decision,
                    source_hash,
                    core_geometry,
                    SourceGeometry.from_candidate(full_dna, current),
                    steps,
                )
            updated = apply_source_backed_weak_terminal_repeat_proposal(
                full_dna, current, proposal
            )
            steps.append(
                BoundaryCompletionStep(
                    tier="weak",
                    side=proposal["side"],
                    repeat_interval=tuple(proposal["repeat_interval"]),
                    spacer_interval=_proposal_spacer_interval(current, proposal),
                    witness_path_ids=bridge_ids,
                    evidence=_step_evidence(proposal),
                )
            )
            weak_used = True
            closed_sides.add(proposal["side"])
        else:
            break

        updated_geometry = SourceGeometry.from_candidate(full_dna, updated)
        if (
            updated_geometry.repeat_intervals in observed_geometries
            or updated_geometry.repeat_count <= len(_repeat_path(current))
        ):
            return _plan(
                "abstain",
                ("non_monotonic_completion",),
                core_category,
                core_row[0],
                core_decision,
                source_hash,
                core_geometry,
                SourceGeometry.from_candidate(full_dna, current),
                steps,
            )
        observed_geometries.add(updated_geometry.repeat_intervals)
        current = updated

    proposed_geometry = SourceGeometry.from_candidate(full_dna, current)
    if proposed_geometry == core_geometry:
        return _plan(
            "unchanged",
            ("source_boundary_fixed_point",),
            core_category,
            core_row[0],
            core_decision,
            source_hash,
            core_geometry,
            proposed_geometry,
        )
    if final_evidence_provider is None:
        return _plan(
            "abstain",
            ("final_evidence_missing",),
            core_category,
            core_row[0],
            core_decision,
            source_hash,
            core_geometry,
            proposed_geometry,
            steps,
        )
    final_candidate = proposed_geometry.materialize(full_dna)
    final_evidence = final_evidence_provider(final_candidate)
    if not isinstance(final_evidence, FinalBoundaryEvidence):
        raise TypeError("final_evidence_provider must return FinalBoundaryEvidence")
    if not final_evidence.passes:
        failed = []
        for field_name in (
            "source_reconstruction_valid",
            "structural_valid",
            "phase_unambiguous",
            "coding_context_assessed",
            "tandem_evidence_assessed",
            "coding_tandem_clear",
        ):
            if not getattr(final_evidence, field_name):
                failed.append("final_{}_failed".format(field_name))
        failed.extend(final_evidence.reason_codes)
        return _plan(
            "abstain",
            tuple(failed),
            core_category,
            core_row[0],
            core_decision,
            source_hash,
            core_geometry,
            proposed_geometry,
            steps,
            final_evidence,
        )
    return _plan(
        "replace",
        ("accepted_core_boundary_completed",),
        core_category,
        core_row[0],
        core_decision,
        source_hash,
        core_geometry,
        proposed_geometry,
        steps,
        final_evidence,
    )


__all__ = [
    "ACCEPTED_CATEGORIES",
    "AcceptedCoreBoundaryCompletionPlan",
    "BoundaryCompletionStep",
    "BoundaryGeometryWitness",
    "CoreAuthorizationEvidence",
    "FinalBoundaryEvidence",
    "SourceGeometry",
    "WitnessContextEvidence",
    "plan_accepted_core_boundary_completion",
]
