"""Pipeline-facing evaluation of accepted arrays split at a circular origin.

The canonical detector still emits ordinary linear ``CrisprCandidate``
objects.  This module considers only accepted candidates close enough to the
two ends of an explicitly circular source to be separated by one legal
spacer.  It never mutates, removes, or merges those child objects.  Instead it
uses :mod:`components.circular_array_join` to build an ordered segmented
geometry and evaluates a freshly reconstructed linear view with the normal
structural filter, the source-CDS tandem veto, and the caller's normal model
scorer.

No topology is inferred here.  Unknown, incomplete, ambiguous, stale, or
linear topology evidence can produce audit decisions but cannot authorize a
join.  Every successful join can therefore be attached to the report as an
additional canonical-locus interpretation while the backward-compatible core
array records remain intact.
"""

from collections import Counter
from dataclasses import dataclass
import hashlib
import json
from numbers import Real
from typing import Callable, Optional, Tuple

from components.accepted_family_tandem_validator import (
    validate_accepted_family_tandem,
)
from components.architecture_rescue_coding_context import predict_source_cds
from components.circular_array_join import (
    GATE_STRUCTURAL,
    GATE_TANDEM,
    ORIENTATION_FORWARD,
    STATUS_JOINED,
    CircularArrayJoinDecision,
    CircularArrayJoinPolicy,
    CircularJoinCompatibilityEvidence,
    CircularJoinGateDecision,
    CircularJoinScoreEvidence,
    CircularTopologyEvidence,
    evaluate_circular_array_join,
    snapshot_circular_children,
    source_sha256,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
)
from components.module_evaluation import build_advanced_fuzzy_search_filter


RULE_NAME = "accepted-edge-source-replayed-circular-join"
RULE_VERSION = "1.0.0"
MAXIMUM_INDIVIDUAL_SPACER_LENGTH = 140
MAXIMUM_REPEAT_FAMILY_EDIT_RATE = 0.20


def _canonical_json(value):
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def _canonical_sha256(value):
    return hashlib.sha256(_canonical_json(value).encode("ascii")).hexdigest()


def _sequence_sha256(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def _levenshtein_distance(first, second):
    """Return exact unit-cost edit distance with bounded working memory."""

    if len(first) > len(second):
        first, second = second, first
    previous = list(range(len(first) + 1))
    for row_index, right_base in enumerate(second, 1):
        current = [row_index]
        for column_index, left_base in enumerate(first, 1):
            current.append(min(
                current[-1] + 1,
                previous[column_index] + 1,
                previous[column_index - 1] + (left_base != right_base),
            ))
        previous = current
    return previous[-1]


def _source_repeat_medoid(repeats):
    """Choose a deterministic exact source repeat representative.

    Exact multiplicity is deliberately retained.  Edge arrays are normally
    dominated by one repeat spelling, making this both cheaper and less
    dependent on an inherited gapped alignment than a legacy consensus.
    """

    counts = Counter(repeats)
    if not counts:
        raise ValueError("a circular child must contain repeats")
    return min(
        counts,
        key=lambda repeat: (-counts[repeat], -len(repeat), repeat),
    )


def _compatibility_metrics(first_candidate, second_candidate, parameters):
    first = _source_repeat_medoid(tuple(first_candidate.list_repeats))
    second = _source_repeat_medoid(tuple(second_candidate.list_repeats))
    distance = _levenshtein_distance(first, second)
    denominator = max(len(first), len(second))
    edit_rate = distance / float(denominator)
    maximum_edits = parameters.get("param_max_edit_distance", 6)
    if (
            not isinstance(maximum_edits, int)
            or isinstance(maximum_edits, bool)
            or maximum_edits < 0):
        raise ValueError("param_max_edit_distance must be non-negative")
    family_compatible = (
        distance <= maximum_edits
        and edit_rate <= MAXIMUM_REPEAT_FAMILY_EDIT_RATE
    )
    # A different nominal repeat length can be a real indel, but without a
    # joint boundary alignment it is not complete evidence that the two child
    # phases share one origin-crossing repeat lattice.  Abstain conservatively.
    phase_compatible = family_compatible and len(first) == len(second)
    return {
        "first_medoid_sha256": _sequence_sha256(first),
        "second_medoid_sha256": _sequence_sha256(second),
        "first_medoid_length": len(first),
        "second_medoid_length": len(second),
        "edit_distance": distance,
        "edit_rate": edit_rate,
        "maximum_edits": maximum_edits,
        "maximum_edit_rate": MAXIMUM_REPEAT_FAMILY_EDIT_RATE,
        "repeat_family_compatible": family_compatible,
        "repeat_phase_compatible": phase_compatible,
    }


def _linearized_content(view):
    elements = []
    for index, repeat in enumerate(view.repeats):
        elements.append(repeat)
        if index < len(view.spacers):
            elements.append(view.spacers[index])
    content = "".join(elements)
    for start, repeat in zip(view.linear_repeat_starts, view.repeats):
        if content[start:start + len(repeat)] != repeat:
            raise RuntimeError(
                "circular linearized repeat failed independent replay"
            )
    return content


def candidate_from_circular_view(view):
    """Rebuild a normal candidate solely from one circular linearized view."""

    content = _linearized_content(view)
    repeat_ends = tuple(
        start + len(repeat)
        for start, repeat in zip(view.linear_repeat_starts, view.repeats)
    )
    candidate = build_source_backed_candidate_from_boundaries(
        content,
        None,
        view.linear_repeat_starts,
        repeat_ends,
    )
    if tuple(candidate.list_repeats) != tuple(view.repeats):
        raise RuntimeError("circular candidate rebuilding changed repeats")
    if tuple(candidate.list_spacers) != tuple(view.spacers):
        raise RuntimeError("circular candidate rebuilding changed spacers")
    return content, candidate


@dataclass(frozen=True)
class CircularAcceptedJoinResult:
    """Immutable audit result for every eligible accepted edge pair."""

    source_sha256: str
    decisions: Tuple[CircularArrayJoinDecision, ...]

    @property
    def joined(self):
        return tuple(
            decision for decision in self.decisions
            if decision.status == STATUS_JOINED
        )

    def as_dict(self):
        return {
            "rule": {
                "name": RULE_NAME,
                "version": RULE_VERSION,
                "accepted_category_indexes": [0, 2],
                "child_calls_mutated": False,
                "topology_inferred": False,
                "maximum_individual_spacer_length": (
                    MAXIMUM_INDIVIDUAL_SPACER_LENGTH
                ),
            },
            "source_sha256": self.source_sha256,
            "attempt_count": len(self.decisions),
            "joined_count": len(self.joined),
            "decisions": [decision.as_dict() for decision in self.decisions],
        }


class _ViewEvaluator:
    """Cache one source-rebuilt view across all independent gates."""

    def __init__(
            self, parameters, score_candidate, source_cds_predictor):
        self.structural_filter = build_advanced_fuzzy_search_filter(parameters)
        self.score_candidate = score_candidate
        self.source_cds_predictor = source_cds_predictor
        self._cached_view_sha256 = None
        self._content = None
        self._candidate = None

    def _resolve(self, view):
        if self._cached_view_sha256 != view.view_sha256:
            self._content, self._candidate = candidate_from_circular_view(view)
            self._cached_view_sha256 = view.view_sha256
        return self._content, self._candidate

    def structural(self, view):
        _content, candidate = self._resolve(view)
        evidence = self.structural_filter.evidence(candidate)
        payload = {
            "subject_view_sha256": view.view_sha256,
            "filter_evidence": evidence.as_dict(),
        }
        return CircularJoinGateDecision(
            gate=GATE_STRUCTURAL,
            passed=evidence.passes,
            reason_code=(
                "fresh_linearized_structural_filter_passed"
                if evidence.passes
                else "fresh_linearized_structural_filter_rejected"
            ),
            subject_view_sha256=view.view_sha256,
            evidence_method=(
                "normal_advanced_fuzzy_search_filter_on_source_rebuilt_view"
            ),
            evidence_sha256=_canonical_sha256(payload),
        )

    def tandem(self, view):
        content, candidate = self._resolve(view)
        cds = self.source_cds_predictor(content)
        decision = validate_accepted_family_tandem(
            content,
            candidate,
            (candidate,),
            cds,
        )
        payload = {
            "subject_view_sha256": view.view_sha256,
            "tandem_decision": decision.as_dict(),
        }
        passed = not decision.supported
        return CircularJoinGateDecision(
            gate=GATE_TANDEM,
            passed=passed,
            reason_code=(
                "fresh_linearized_coding_tandem_not_supported"
                if passed
                else "fresh_linearized_coding_tandem_supported"
            ),
            subject_view_sha256=view.view_sha256,
            evidence_method=(
                "source_cds_tandem_validator_on_source_rebuilt_view"
            ),
            evidence_sha256=_canonical_sha256(payload),
        )

    def score(self, view):
        _content, candidate = self._resolve(view)
        score = self.score_candidate(candidate)
        if (
                not isinstance(score, Real)
                or isinstance(score, bool)
                or not 0.0 <= float(score) <= 1.0):
            raise ValueError("circular view scorer returned an invalid score")
        payload = {
            "subject_view_sha256": view.view_sha256,
            "score": float(score),
            "candidate_repeats": tuple(candidate.list_repeats),
            "candidate_spacers": tuple(candidate.list_spacers),
            "candidate_starts": tuple(candidate.list_repeat_starts),
        }
        return CircularJoinScoreEvidence(
            score=float(score),
            subject_view_sha256=view.view_sha256,
            evidence_method="normal_model_ensemble_on_source_rebuilt_view",
            evidence_sha256=_canonical_sha256(payload),
        )


def _accepted_candidates(categories):
    if not isinstance(categories, (tuple, list)) or len(categories) != 5:
        raise ValueError("categories must contain five mappings")
    accepted = []
    for category_index in (0, 2):
        category = categories[category_index]
        if not isinstance(category, dict):
            raise ValueError("accepted categories must be mappings")
        for interval in sorted(category):
            rows = category[interval]
            if isinstance(rows, (str, bytes)) or rows is None:
                raise ValueError("accepted category rows must be a sequence")
            for position, row in enumerate(rows):
                if not isinstance(row, (tuple, list)) or len(row) < 2:
                    raise ValueError("accepted category row is malformed")
                accepted.append((
                    category_index,
                    interval,
                    position,
                    row[1],
                ))
    return tuple(accepted)


def _edge_pairs(source, categories):
    accepted = _accepted_candidates(categories)
    geometries = snapshot_circular_children(
        source, (record[3] for record in accepted)
    )
    snapshots = tuple(zip(accepted, geometries))
    heads = tuple(
        item for item in snapshots
        if item[1].source_interval[0] <= MAXIMUM_INDIVIDUAL_SPACER_LENGTH
    )
    tails = tuple(
        item for item in snapshots
        if len(source) - item[1].source_interval[1]
        <= MAXIMUM_INDIVIDUAL_SPACER_LENGTH
    )
    pairs = []
    for head_record, head_geometry in heads:
        for tail_record, tail_geometry in tails:
            if head_geometry.geometry_sha256 == tail_geometry.geometry_sha256:
                continue
            junction_length = (
                len(source) - tail_geometry.source_interval[1]
                + head_geometry.source_interval[0]
            )
            if not 0 <= junction_length <= MAXIMUM_INDIVIDUAL_SPACER_LENGTH:
                continue
            pairs.append((
                head_record,
                head_geometry,
                tail_record,
                tail_geometry,
            ))
    return tuple(sorted(
        pairs,
        key=lambda pair: (
            pair[1].source_interval,
            pair[3].source_interval,
            pair[1].geometry_sha256,
            pair[3].geometry_sha256,
        ),
    ))


def _compatibility_evidence(
        source, pair, parameters, ambiguous_compatible_pair):
    head_record, head_geometry, tail_record, tail_geometry = pair
    metrics = _compatibility_metrics(
        head_record[3], tail_record[3], parameters
    )
    hashes = tuple(sorted((
        head_geometry.geometry_sha256,
        tail_geometry.geometry_sha256,
    )))
    payload = {
        "rule": RULE_NAME,
        "version": RULE_VERSION,
        "source_sha256": source_sha256(source),
        "child_geometry_sha256s": hashes,
        "metrics": metrics,
        "orientation": ORIENTATION_FORWARD,
        "assessment_ambiguous": ambiguous_compatible_pair,
    }
    return CircularJoinCompatibilityEvidence(
        child_geometry_sha256s=hashes,
        repeat_family_compatible=metrics["repeat_family_compatible"],
        repeat_phase_compatible=metrics["repeat_phase_compatible"],
        orientation_compatible=True,
        orientation=ORIENTATION_FORWARD,
        assessment_complete=True,
        assessment_ambiguous=ambiguous_compatible_pair,
        evidence_method=(
            "source_repeat_medoid_edit_distance_and_equal_length_phase"
        ),
        evidence_sha256=_canonical_sha256(payload),
    )


def evaluate_accepted_circular_edge_joins(
        source,
        categories,
        topology_evidence,
        parameters,
        score_candidate: Callable[[object], float],
        *,
        source_cds_predictor: Optional[Callable[[str], tuple]] = None):
    """Evaluate every accepted origin-edge pair without changing categories."""

    if not isinstance(source, str) or not source:
        raise ValueError("source must be one non-empty DNA string")
    if not isinstance(topology_evidence, CircularTopologyEvidence):
        raise TypeError("topology_evidence has the wrong type")
    if not isinstance(parameters, dict):
        raise TypeError("parameters must be a mapping")
    if not callable(score_candidate):
        raise TypeError("score_candidate must be callable")
    if source_cds_predictor is None:
        source_cds_predictor = predict_source_cds
    if not callable(source_cds_predictor):
        raise TypeError("source_cds_predictor must be callable")

    pairs = _edge_pairs(source, categories)
    preliminary = [
        _compatibility_metrics(pair[0][3], pair[2][3], parameters)
        for pair in pairs
    ]
    compatible_count = sum(
        metrics["repeat_family_compatible"]
        and metrics["repeat_phase_compatible"]
        for metrics in preliminary
    )
    evaluator = _ViewEvaluator(
        parameters,
        score_candidate,
        source_cds_predictor,
    )
    policy = CircularArrayJoinPolicy(
        minimum_junction_spacer_length=0,
        maximum_junction_spacer_length=MAXIMUM_INDIVIDUAL_SPACER_LENGTH,
        possible_score_threshold=0.5,
    )
    decisions = []
    for pair, metrics in zip(pairs, preliminary):
        is_compatible = (
            metrics["repeat_family_compatible"]
            and metrics["repeat_phase_compatible"]
        )
        compatibility = _compatibility_evidence(
            source,
            pair,
            parameters,
            ambiguous_compatible_pair=(
                is_compatible and compatible_count > 1
            ),
        )
        decisions.append(evaluate_circular_array_join(
            source,
            pair[0][3],
            pair[2][3],
            topology_evidence,
            compatibility,
            policy,
            structural_callback=evaluator.structural,
            tandem_callback=evaluator.tandem,
            scoring_callback=evaluator.score,
        ))

    result = CircularAcceptedJoinResult(
        source_sha256=source_sha256(source),
        decisions=tuple(decisions),
    )
    # Exercise the JSON boundary here so runtime objects or NaNs cannot be
    # carried into report construction unnoticed.
    _canonical_json(result.as_dict())
    return result


__all__ = (
    "CircularAcceptedJoinResult",
    "MAXIMUM_INDIVIDUAL_SPACER_LENGTH",
    "RULE_NAME",
    "RULE_VERSION",
    "candidate_from_circular_view",
    "evaluate_accepted_circular_edge_joins",
)
