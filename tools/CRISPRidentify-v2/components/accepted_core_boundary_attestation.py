"""Source-first phase attestation for accepted-core boundary completion.

Boundary completion and array existence are deliberately separate decisions.
An accepted core may establish that an array exists, but it does not prove
that a newly completed representation has exact nucleotide boundaries.  This
module supplies a pure, fail-closed attestation for that narrower claim.

The caller must provide an explicitly complete family and immutable lineage
entries for every family member and for the proposed completed candidate.
Every geometry is independently reconstructed from the original source.  The
proposed geometry is then challenged in three ways:

* a uniform +/-12-base source-only phase search;
* an immediate +/-1-base boundary competition, including sharp source-column
  boundary evidence; and
* every bounded same-occurrence phase represented in the complete family.

Longer same-occurrence family paths and incompatible overlapping paths remain
ambiguous.  Missing completeness attestations, malformed lineage, unknown
derivations, invalid source slices, and comparisons outside the calibrated
rule are unassessed.  Neither state authorizes ``phase_unambiguous=True``.

Coordinates are zero-based and end-exclusive.  No classifier score or
category is accepted by this API.
"""

from dataclasses import dataclass
import hashlib
import math
from typing import Optional, Tuple

from components.accepted_core_boundary_completion import (
    SourceGeometry,
)
from components.components_evaluated_arrays_enhancement import (
    SOURCE_CROSS_FAMILY_PHASE_RULE,
    build_source_backed_candidate_from_boundaries,
    compare_source_backed_phase_candidates,
    resolve_source_backed_global_phase,
)


BOUNDARY_EXACT = "EXACT"
BOUNDARY_AMBIGUOUS = "AMBIGUOUS"
BOUNDARY_UNASSESSED = "UNASSESSED"
COORDINATE_SYSTEM = "zero_based_half_open"

KNOWN_LINEAGE_DERIVATIONS = frozenset((
    "legacy_observed",
    "legacy_eligible",
    "complete_fuzzy_search",
    "local_chain_family_parent",
    "maximal_local_chain",
    "terminal_augmentation",
    "weak_terminal_augmentation",
    "exact_shared_repeat_union",
    "accepted_core_boundary_completion",
))
ROOT_LINEAGE_DERIVATIONS = frozenset((
    "legacy_observed",
    "legacy_eligible",
    "complete_fuzzy_search",
    "local_chain_family_parent",
    "maximal_local_chain",
))
DERIVED_LINEAGE_DERIVATIONS = (
    KNOWN_LINEAGE_DERIVATIONS.difference(ROOT_LINEAGE_DERIVATIONS)
)


def _source_hash(source):
    return hashlib.sha256(source.upper().encode("ascii")).hexdigest()


def _path(geometry):
    return geometry.repeat_intervals


def _path_interval(path):
    return path[0][0], path[-1][1]


def _intervals_overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _is_contiguous_subpath(shorter, longer):
    shorter = tuple(shorter)
    longer = tuple(longer)
    if not shorter or len(shorter) > len(longer):
        return False
    return any(
        longer[offset:offset + len(shorter)] == shorter
        for offset in range(len(longer) - len(shorter) + 1)
    )


def _bounded_path_equivalent(first, second, maximum_shift):
    return (
        len(first) == len(second)
        and all(
            abs(first_start - second_start) <= maximum_shift
            and abs(first_end - second_end) <= maximum_shift
            for (
                (first_start, first_end),
                (second_start, second_end),
            ) in zip(first, second)
        )
    )


def _bounded_alignments(shorter, longer, maximum_shift):
    if len(shorter) > len(longer):
        return ()
    return tuple(
        offset
        for offset in range(len(longer) - len(shorter) + 1)
        if _bounded_path_equivalent(
            shorter,
            longer[offset:offset + len(shorter)],
            maximum_shift,
        )
    )


def _geometry_as_dict(geometry):
    if geometry is None:
        return None
    return {
        "coordinate_system": COORDINATE_SYSTEM,
        "source_interval": list(geometry.source_interval),
        "repeat_count": geometry.repeat_count,
        "repeat_intervals": [
            list(interval) for interval in geometry.repeat_intervals
        ],
    }


@dataclass(frozen=True)
class BoundaryPhaseCompletenessEvidence:
    """Explicit caller attestations required before exactness is assessable."""

    source_family_complete: bool
    lineage_complete: bool
    same_occurrence_alternatives_complete: bool
    augmentation_converged: bool

    def __post_init__(self):
        for field_name in (
            "source_family_complete",
            "lineage_complete",
            "same_occurrence_alternatives_complete",
            "augmentation_converged",
        ):
            if not isinstance(getattr(self, field_name), bool):
                raise ValueError("{} must be a boolean".format(field_name))

    @property
    def complete(self):
        return all((
            self.source_family_complete,
            self.lineage_complete,
            self.same_occurrence_alternatives_complete,
            self.augmentation_converged,
        ))


@dataclass(frozen=True)
class BoundaryPhaseLineageEntry:
    """Immutable ancestry and geometry snapshot for one runtime candidate."""

    candidate: object
    source_geometry: SourceGeometry
    path_id: str
    family_id: str
    derivation: str
    parent_path_ids: Tuple[str, ...] = ()

    @classmethod
    def from_candidate(
            cls, source, candidate, path_id, family_id, derivation,
            parent_path_ids=()):
        """Create an immutable entry after exact source reconstruction."""

        return cls(
            candidate=candidate,
            source_geometry=SourceGeometry.from_candidate(source, candidate),
            path_id=path_id,
            family_id=family_id,
            derivation=derivation,
            parent_path_ids=tuple(parent_path_ids),
        )


@dataclass(frozen=True)
class BoundaryPhaseAlternative:
    """One source-valid geometry that challenged the proposed coordinates."""

    path_id: str
    relationship: str
    geometry: SourceGeometry

    def as_dict(self):
        return {
            "path_id": self.path_id,
            "relationship": self.relationship,
            "geometry": _geometry_as_dict(self.geometry),
        }


@dataclass(frozen=True)
class BoundaryPhaseCompetition:
    """Runtime-object-free source evidence for one bounded competition."""

    scope: str
    path_ids: Tuple[str, ...]
    status: str
    reason: str
    selected_path_id: Optional[str]
    preferred_path_id: Optional[str]
    winner_margin: Optional[float]
    source_scores: Tuple[Tuple[str, Optional[float]], ...]
    hypotheses_evaluated: int

    def as_dict(self):
        return {
            "scope": self.scope,
            "path_ids": list(self.path_ids),
            "status": self.status,
            "reason": self.reason,
            "selected_path_id": self.selected_path_id,
            "preferred_path_id": self.preferred_path_id,
            "winner_margin": self.winner_margin,
            "source_scores": [
                {"path_id": path_id, "source_score": source_score}
                for path_id, source_score in self.source_scores
            ],
            "hypotheses_evaluated": self.hypotheses_evaluated,
        }


@dataclass(frozen=True)
class SourcePhaseBoundaryAttestation:
    """Final source-only phase/boundary decision."""

    status: str
    reason_codes: Tuple[str, ...]
    source_sha256: str
    proposed_geometry: Optional[SourceGeometry]
    alternatives: Tuple[BoundaryPhaseAlternative, ...]
    competitions: Tuple[BoundaryPhaseCompetition, ...]

    @property
    def phase_unambiguous(self):
        return self.status == BOUNDARY_EXACT

    def as_dict(self):
        return {
            "status": self.status,
            "phase_unambiguous": self.phase_unambiguous,
            "reason_codes": list(self.reason_codes),
            "source_sha256": self.source_sha256,
            "proposed_geometry": _geometry_as_dict(self.proposed_geometry),
            "alternatives": [
                alternative.as_dict() for alternative in self.alternatives
            ],
            "competitions": [
                competition.as_dict() for competition in self.competitions
            ],
        }


def _result(
        status, reasons, source_sha256, proposed_geometry=None,
        alternatives=(), competitions=()):
    return SourcePhaseBoundaryAttestation(
        status=status,
        reason_codes=tuple(sorted(set(reasons))),
        source_sha256=source_sha256,
        proposed_geometry=proposed_geometry,
        alternatives=tuple(sorted(
            alternatives,
            key=lambda alternative: (
                alternative.geometry.repeat_intervals,
                alternative.path_id,
                alternative.relationship,
            ),
        )),
        competitions=tuple(competitions),
    )


def _validate_source(source):
    if not isinstance(source, str) or not source:
        return False
    try:
        source.encode("ascii")
    except UnicodeEncodeError:
        return False
    return True


def _lineage_cycle(path_ids, parents_by_path_id):
    visiting = set()
    visited = set()

    def visit(path_id):
        if path_id in visiting:
            return True
        if path_id in visited:
            return False
        visiting.add(path_id)
        if any(visit(parent) for parent in parents_by_path_id[path_id]):
            return True
        visiting.remove(path_id)
        visited.add(path_id)
        return False

    return any(visit(path_id) for path_id in path_ids)


def _validate_lineage(
        source, proposed_candidate, candidates, lineage_entries,
        observed_geometries):
    """Return ``(reasons, entries_by_candidate_id, entries_by_path_id)``."""

    reasons = []
    try:
        entries = tuple(lineage_entries)
    except TypeError:
        return (
            ["lineage_not_iterable"],
            {},
            {},
        )
    if not entries or any(
            not isinstance(entry, BoundaryPhaseLineageEntry)
            for entry in entries):
        return (
            ["lineage_entries_malformed"],
            {},
            {},
        )

    by_candidate_id = {}
    by_path_id = {}
    for entry in entries:
        candidate_id = id(entry.candidate)
        if candidate_id in by_candidate_id:
            reasons.append("lineage_candidate_duplicated")
        else:
            by_candidate_id[candidate_id] = entry
        if (
                not isinstance(entry.path_id, str)
                or not entry.path_id
                or entry.path_id in by_path_id):
            reasons.append("lineage_path_id_invalid_or_duplicated")
        else:
            by_path_id[entry.path_id] = entry
        if not isinstance(entry.family_id, str) or not entry.family_id:
            reasons.append("lineage_family_id_invalid")
        if entry.derivation not in KNOWN_LINEAGE_DERIVATIONS:
            reasons.append("unknown_lineage_derivation")
        if (
                not isinstance(entry.parent_path_ids, tuple)
                or any(
                    not isinstance(parent_id, str) or not parent_id
                    for parent_id in entry.parent_path_ids
                )
                or len(set(entry.parent_path_ids))
                != len(entry.parent_path_ids)):
            reasons.append("lineage_parent_ids_malformed")
        observed = observed_geometries.get(candidate_id)
        if (
                not isinstance(entry.source_geometry, SourceGeometry)
                or observed is None
                or entry.source_geometry != observed):
            reasons.append("lineage_source_geometry_mismatch")

    expected_ids = {id(candidate) for candidate in candidates}
    if set(by_candidate_id) != expected_ids:
        reasons.append("lineage_candidate_universe_incomplete")
    family_ids = {
        entry.family_id
        for entry in entries
        if isinstance(entry.family_id, str) and entry.family_id
    }
    if len(family_ids) != 1:
        reasons.append("lineage_family_identity_inconsistent")

    known_path_ids = set(by_path_id)
    parents_by_path_id = {}
    for path_id, entry in by_path_id.items():
        parent_ids = (
            entry.parent_path_ids
            if (
                isinstance(entry.parent_path_ids, tuple)
                and all(
                    isinstance(parent_id, str) and parent_id
                    for parent_id in entry.parent_path_ids
                )
                and len(set(entry.parent_path_ids))
                == len(entry.parent_path_ids)
            )
            else ()
        )
        if any(parent_id not in known_path_ids for parent_id in parent_ids):
            reasons.append("lineage_parent_unknown")
        if path_id in parent_ids:
            reasons.append("lineage_self_parent")
        if (
                entry.derivation in ROOT_LINEAGE_DERIVATIONS
                and parent_ids):
            reasons.append("lineage_root_has_parent")
        if (
                entry.derivation in DERIVED_LINEAGE_DERIVATIONS
                and not parent_ids):
            reasons.append("derived_lineage_parent_missing")
        parents_by_path_id[path_id] = tuple(
            parent_id for parent_id in parent_ids
            if parent_id in known_path_ids
        )
    if (
            len(parents_by_path_id) == len(by_path_id)
            and _lineage_cycle(tuple(by_path_id), parents_by_path_id)):
        reasons.append("lineage_cycle")

    proposed_entry = by_candidate_id.get(id(proposed_candidate))
    if proposed_entry is None:
        reasons.append("proposed_lineage_missing")
    elif (
            proposed_entry.derivation
            != "accepted_core_boundary_completion"):
        reasons.append("proposed_lineage_derivation_unknown")
    else:
        proposed_path = observed_geometries[id(proposed_candidate)].repeat_intervals
        valid_parent = any(
            parent_id in by_path_id
            and len(by_path_id[parent_id].source_geometry.repeat_intervals)
            < len(proposed_path)
            and _is_contiguous_subpath(
                by_path_id[parent_id].source_geometry.repeat_intervals,
                proposed_path,
            )
            for parent_id in parents_by_path_id.get(
                proposed_entry.path_id, ()
            )
        )
        if not valid_parent:
            reasons.append("proposed_lineage_not_terminal_extension")

    return reasons, by_candidate_id, by_path_id


def _geometry_within_rule(starts, ends, source_length, rule):
    if (
            not starts
            or len(starts) != len(ends)
            or len(starts) < rule["minimum_repeats"]):
        return False
    lengths = tuple(end - start for start, end in zip(starts, ends))
    spacers = tuple(
        next_start - end
        for end, next_start in zip(ends, starts[1:])
    )
    if (
            min(starts) < 0
            or max(ends) > source_length
            or min(lengths) < rule["minimum_repeat_length"]
            or max(lengths) > rule["maximum_repeat_length"]):
        return False
    mean_length = sum(lengths) / float(len(lengths))
    return (
        rule["minimum_mean_repeat_length"]
        <= mean_length
        <= rule["maximum_mean_repeat_length"]
        and bool(spacers)
        and all(
            rule["minimum_spacer_length"]
            <= spacer
            <= rule["maximum_spacer_length"]
            for spacer in spacers
        )
    )


def _comparison_record(scope, path_ids, result):
    evidence = result["candidate_evidence"]
    source_scores = tuple(
        (
            path_id,
            (
                None
                if item["source_score"] is None
                else float(item["source_score"])
            ),
        )
        for path_id, item in zip(path_ids, evidence)
    )
    ranked = sorted(
        (
            (score, path_id)
            for path_id, score in source_scores
            if score is not None
        ),
        key=lambda item: (-item[0], item[1]),
    )
    preferred = ranked[0][1] if ranked else None
    selected = (
        path_ids[result["winner_candidate_index"]]
        if result["status"] == "selected"
        else None
    )
    margin = result["winner_margin"]
    return BoundaryPhaseCompetition(
        scope=scope,
        path_ids=tuple(path_ids),
        status=result["status"],
        reason=result["reason"],
        selected_path_id=selected,
        preferred_path_id=preferred,
        winner_margin=None if margin is None else float(margin),
        source_scores=source_scores,
        hypotheses_evaluated=len(evidence),
    )


def _candidate_from_geometry(source, geometry):
    return build_source_backed_candidate_from_boundaries(
        source,
        None,
        tuple(start for start, _end in geometry.repeat_intervals),
        tuple(end for _start, end in geometry.repeat_intervals),
    )


def _uniform_phase_check(source, proposed_candidate, proposed_entry, rule):
    resolved, metadata = resolve_source_backed_global_phase(
        source,
        proposed_candidate,
        rule,
    )
    shift = (metadata["left_shift"], metadata["right_shift"])
    preferred_path_id = (
        proposed_entry.path_id
        if shift == (0, 0)
        else "uniform-phase:{:+d}:{:+d}".format(*shift)
    )
    source_scores = [(proposed_entry.path_id, metadata["current_score"])]
    alternative = None
    if shift != (0, 0):
        try:
            alternative_geometry = SourceGeometry.from_candidate(
                source, resolved
            )
            if alternative_geometry == proposed_entry.source_geometry:
                starts = tuple(
                    start + shift[0]
                    for start, _end
                    in proposed_entry.source_geometry.repeat_intervals
                )
                ends = tuple(
                    end + shift[1]
                    for _start, end
                    in proposed_entry.source_geometry.repeat_intervals
                )
                shifted = build_source_backed_candidate_from_boundaries(
                    source, None, starts, ends
                )
                alternative_geometry = SourceGeometry.from_candidate(
                    source, shifted
                )
            alternative = BoundaryPhaseAlternative(
                preferred_path_id,
                "uniform_source_phase_hypothesis",
                alternative_geometry,
            )
        except Exception:
            return (
                "unassessed",
                "uniform_phase_winner_geometry_invalid",
                None,
                None,
            )
        source_scores.append(
            (preferred_path_id, metadata["winner_score"])
        )

    margin = metadata["winner_margin"]
    current_wins = (
        metadata["reason"] == "current_geometry_wins"
        and shift == (0, 0)
        and isinstance(margin, (int, float))
        and math.isfinite(float(margin))
        and margin >= rule["minimum_winner_margin"]
    )
    status = "selected" if current_wins else "abstained"
    competition = BoundaryPhaseCompetition(
        scope="uniform_global_phase",
        path_ids=tuple(path_id for path_id, _score in source_scores),
        status=status,
        reason=metadata["reason"],
        selected_path_id=proposed_entry.path_id if current_wins else None,
        preferred_path_id=preferred_path_id,
        winner_margin=None if margin is None else float(margin),
        source_scores=tuple(
            (
                path_id,
                None if score is None else float(score),
            )
            for path_id, score in source_scores
        ),
        hypotheses_evaluated=int(metadata["hypotheses_scored"]),
    )
    if current_wins:
        return "selected", None, competition, None
    if metadata["reason"] in (
            "insufficient_repeats",
            "insufficient_valid_hypotheses",
            "geometry_outside_source_phase_rule"):
        return (
            "unassessed",
            "uniform_phase_{}".format(metadata["reason"]),
            competition,
            alternative,
        )
    if (
            metadata["reason"] == "current_geometry_wins"
            and margin is not None
            and margin < rule["minimum_winner_margin"]):
        reason = "uniform_phase_weak_source_score_margin"
    else:
        reason = "uniform_phase_{}".format(metadata["reason"])
    return "ambiguous", reason, competition, alternative


def _local_boundary_check(source, proposed_candidate, proposed_entry, rule):
    geometry = proposed_entry.source_geometry
    starts = tuple(start for start, _end in geometry.repeat_intervals)
    ends = tuple(end for _start, end in geometry.repeat_intervals)
    candidates = [proposed_candidate]
    path_ids = [proposed_entry.path_id]
    geometries = [geometry]
    shifts = [(0, 0)]
    for left_shift in (-1, 0, 1):
        for right_shift in (-1, 0, 1):
            if (left_shift, right_shift) == (0, 0):
                continue
            proposed_starts = tuple(
                start + left_shift for start in starts
            )
            proposed_ends = tuple(end + right_shift for end in ends)
            if not _geometry_within_rule(
                    proposed_starts, proposed_ends, len(source), rule):
                continue
            candidate = build_source_backed_candidate_from_boundaries(
                source, None, proposed_starts, proposed_ends
            )
            candidate_geometry = SourceGeometry.from_candidate(
                source, candidate
            )
            candidates.append(candidate)
            path_ids.append(
                "local-phase:{:+d}:{:+d}".format(
                    left_shift, right_shift
                )
            )
            geometries.append(candidate_geometry)
            shifts.append((left_shift, right_shift))

    left_tested = any(
        left_shift != 0 and right_shift == 0
        for left_shift, right_shift in shifts
    )
    right_tested = any(
        left_shift == 0 and right_shift != 0
        for left_shift, right_shift in shifts
    )
    if not left_tested or not right_tested:
        missing = "left" if not left_tested else "right"
        return (
            "unassessed",
            "{}_boundary_not_testable".format(missing),
            None,
            None,
        )
    try:
        comparison = compare_source_backed_phase_candidates(
            source, candidates, rule
        )
    except (TypeError, ValueError):
        return (
            "unassessed",
            "local_boundary_comparison_invalid",
            None,
            None,
        )
    competition = _comparison_record(
        "immediate_local_boundary", tuple(path_ids), comparison
    )
    if (
            comparison["status"] == "selected"
            and comparison["winner_candidate_index"] == 0):
        return "selected", None, competition, None

    evidence = comparison["candidate_evidence"]
    ranked = sorted(
        (
            (item["source_score"], index)
            for index, item in enumerate(evidence)
            if item["source_score"] is not None
        ),
        key=lambda item: (-item[0], item[1]),
    )
    alternative_index = (
        ranked[0][1]
        if ranked and ranked[0][1] != 0
        else (ranked[1][1] if len(ranked) > 1 else None)
    )
    alternative = (
        None
        if alternative_index is None
        else BoundaryPhaseAlternative(
            path_ids[alternative_index],
            "immediate_local_phase_hypothesis",
            geometries[alternative_index],
        )
    )
    return (
        "ambiguous",
        "local_boundary_{}".format(comparison["reason"]),
        competition,
        alternative,
    )


def attest_accepted_core_boundary_phase(
        source,
        proposed_candidate,
        complete_family_candidates,
        lineage_entries,
        completeness,
        rule=None):
    """Attest whether completed coordinates are source-exact.

    ``complete_family_candidates`` may include or omit ``proposed_candidate``;
    it is normalized into the audited universe exactly once.  ``lineage_entries``
    must contain exactly one entry for every normalized candidate.  The
    proposed entry must be an ``accepted_core_boundary_completion`` derived
    from a proper contiguous family subpath.

    The return value is always fail closed.  ``phase_unambiguous`` is true only
    for ``EXACT``.  ``AMBIGUOUS`` means source-valid alternatives remain;
    ``UNASSESSED`` means the supplied universe, lineage, source, or calibrated
    comparison was insufficient to make even that bounded claim.
    """

    if not _validate_source(source):
        return _result(
            BOUNDARY_UNASSESSED,
            ("source_invalid",),
            "",
        )
    source_sha256 = _source_hash(source)
    try:
        family = tuple(complete_family_candidates)
    except TypeError:
        return _result(
            BOUNDARY_UNASSESSED,
            ("complete_family_not_iterable",),
            source_sha256,
        )
    if not family:
        return _result(
            BOUNDARY_UNASSESSED,
            ("complete_family_empty",),
            source_sha256,
        )
    if any(candidate is None for candidate in family):
        return _result(
            BOUNDARY_UNASSESSED,
            ("complete_family_candidate_missing",),
            source_sha256,
        )
    family_ids = tuple(id(candidate) for candidate in family)
    if len(set(family_ids)) != len(family_ids):
        return _result(
            BOUNDARY_UNASSESSED,
            ("complete_family_candidates_duplicated",),
            source_sha256,
        )
    candidates = (
        family
        if any(candidate is proposed_candidate for candidate in family)
        else family + (proposed_candidate,)
    )

    observed_geometries = {}
    source_failures = []
    for candidate in candidates:
        try:
            observed_geometries[id(candidate)] = (
                SourceGeometry.from_candidate(source, candidate)
            )
        except Exception:
            source_failures.append(
                "proposed_source_reconstruction_failed"
                if candidate is proposed_candidate
                else "family_source_reconstruction_failed"
            )
    proposed_geometry = observed_geometries.get(id(proposed_candidate))
    if source_failures:
        return _result(
            BOUNDARY_UNASSESSED,
            tuple(source_failures),
            source_sha256,
            proposed_geometry,
        )

    reasons = []
    if not isinstance(completeness, BoundaryPhaseCompletenessEvidence):
        reasons.append("completeness_evidence_missing_or_malformed")
    else:
        for field_name in (
                "source_family_complete",
                "lineage_complete",
                "same_occurrence_alternatives_complete",
                "augmentation_converged"):
            if not getattr(completeness, field_name):
                reasons.append("{}_false".format(field_name))

    lineage_reasons, by_candidate_id, by_path_id = _validate_lineage(
        source,
        proposed_candidate,
        candidates,
        lineage_entries,
        observed_geometries,
    )
    reasons.extend(lineage_reasons)
    if reasons:
        return _result(
            BOUNDARY_UNASSESSED,
            reasons,
            source_sha256,
            proposed_geometry,
        )

    try:
        phase_rule = dict(
            SOURCE_CROSS_FAMILY_PHASE_RULE if rule is None else rule
        )
    except (TypeError, ValueError):
        return _result(
            BOUNDARY_UNASSESSED,
            ("phase_rule_malformed",),
            source_sha256,
            proposed_geometry,
        )
    # The public source-phase functions perform the complete rule validation.
    # This lightweight access check produces an ordinary fail-closed result
    # for malformed integration data instead of leaking a partial attestation.
    try:
        if set(phase_rule) != set(SOURCE_CROSS_FAMILY_PHASE_RULE):
            raise KeyError("phase rule field set changed")
        maximum_shift = phase_rule["maximum_boundary_shift"]
        minimum_margin = phase_rule["minimum_winner_margin"]
    except (KeyError, TypeError):
        return _result(
            BOUNDARY_UNASSESSED,
            ("phase_rule_malformed",),
            source_sha256,
            proposed_geometry,
        )
    if (
            not isinstance(maximum_shift, int)
            or isinstance(maximum_shift, bool)
            or maximum_shift < 1
            or isinstance(minimum_margin, bool)
            or not isinstance(minimum_margin, (int, float))
            or not math.isfinite(float(minimum_margin))):
        return _result(
            BOUNDARY_UNASSESSED,
            ("phase_rule_malformed",),
            source_sha256,
            proposed_geometry,
        )

    proposed_entry = by_candidate_id[id(proposed_candidate)]
    ambiguous_reasons = []
    unassessed_reasons = []
    alternatives = []
    competitions = []

    try:
        (
            uniform_state,
            uniform_reason,
            uniform_competition,
            uniform_alternative,
        ) = _uniform_phase_check(
            source, proposed_candidate, proposed_entry, phase_rule
        )
    except (KeyError, TypeError, ValueError):
        uniform_state = "unassessed"
        uniform_reason = "uniform_phase_rule_or_geometry_invalid"
        uniform_competition = None
        uniform_alternative = None
    if uniform_competition is not None:
        competitions.append(uniform_competition)
    if uniform_alternative is not None:
        alternatives.append(uniform_alternative)
    if uniform_state == "ambiguous":
        ambiguous_reasons.append(uniform_reason)
    elif uniform_state == "unassessed":
        unassessed_reasons.append(uniform_reason)

    try:
        (
            local_state,
            local_reason,
            local_competition,
            local_alternative,
        ) = _local_boundary_check(
            source, proposed_candidate, proposed_entry, phase_rule
        )
    except (KeyError, TypeError, ValueError):
        local_state = "unassessed"
        local_reason = "local_boundary_rule_or_geometry_invalid"
        local_competition = None
        local_alternative = None
    if local_competition is not None:
        competitions.append(local_competition)
    if local_alternative is not None:
        alternatives.append(local_alternative)
    if local_state == "ambiguous":
        ambiguous_reasons.append(local_reason)
    elif local_state == "unassessed":
        unassessed_reasons.append(local_reason)

    proposed_path = proposed_geometry.repeat_intervals
    phase_groups = {}
    ordered_entries = sorted(
        (
            entry for entry in by_path_id.values()
            if entry.candidate is not proposed_candidate
        ),
        key=lambda entry: (
            entry.source_geometry.repeat_intervals,
            entry.path_id,
        ),
    )
    for entry in ordered_entries:
        member_path = entry.source_geometry.repeat_intervals
        if member_path == proposed_path:
            continue
        if len(member_path) <= len(proposed_path):
            alignments = _bounded_alignments(
                member_path, proposed_path, maximum_shift
            )
            if not alignments:
                if _intervals_overlap(
                        _path_interval(member_path),
                        _path_interval(proposed_path)):
                    alternatives.append(BoundaryPhaseAlternative(
                        entry.path_id,
                        "overlapping_incompatible_family_geometry",
                        entry.source_geometry,
                    ))
                    ambiguous_reasons.append(
                        "overlapping_incompatible_family_geometry"
                    )
                continue
            if len(alignments) != 1:
                alternatives.append(BoundaryPhaseAlternative(
                    entry.path_id,
                    "multiply_aligned_family_geometry",
                    entry.source_geometry,
                ))
                ambiguous_reasons.append(
                    "same_occurrence_alignment_not_unique"
                )
                continue
            offset = alignments[0]
            reference_path = proposed_path[
                offset:offset + len(member_path)
            ]
            if member_path == reference_path:
                continue
            alternatives.append(BoundaryPhaseAlternative(
                entry.path_id,
                "same_occurrence_phase_alternative",
                entry.source_geometry,
            ))
            phase_groups.setdefault(
                (offset, len(member_path)), []
            ).append(entry)
        else:
            alignments = _bounded_alignments(
                proposed_path, member_path, maximum_shift
            )
            if alignments:
                alternatives.append(BoundaryPhaseAlternative(
                    entry.path_id,
                    "longer_same_occurrence_family_path",
                    entry.source_geometry,
                ))
                ambiguous_reasons.append(
                    "longer_same_occurrence_family_path"
                )
            elif _intervals_overlap(
                    _path_interval(member_path),
                    _path_interval(proposed_path)):
                alternatives.append(BoundaryPhaseAlternative(
                    entry.path_id,
                    "overlapping_incompatible_family_geometry",
                    entry.source_geometry,
                ))
                ambiguous_reasons.append(
                    "overlapping_incompatible_family_geometry"
                )

    for (offset, repeat_count), entries in sorted(phase_groups.items()):
        if repeat_count < phase_rule["minimum_repeats"]:
            unassessed_reasons.append(
                "family_phase_competitor_too_short"
            )
            continue
        reference_geometry = SourceGeometry(
            repeat_intervals=proposed_geometry.repeat_intervals[
                offset:offset + repeat_count
            ],
            repeat_sequences=proposed_geometry.repeat_sequences[
                offset:offset + repeat_count
            ],
            spacer_sequences=proposed_geometry.spacer_sequences[
                offset:offset + repeat_count - 1
            ],
        )
        competing_paths = [
            reference_geometry.repeat_intervals
        ] + [
            entry.source_geometry.repeat_intervals for entry in entries
        ]
        if any(
                not _bounded_path_equivalent(
                    first, second, maximum_shift
                )
                for first_index, first
                in enumerate(competing_paths[:-1])
                for second in competing_paths[first_index + 1:]):
            ambiguous_reasons.append(
                "wide_same_occurrence_phase_family"
            )
            continue
        reference = _candidate_from_geometry(source, reference_geometry)
        path_ids = (
            "{}:subpath:{}:{}".format(
                proposed_entry.path_id, offset, offset + repeat_count
            ),
        ) + tuple(entry.path_id for entry in entries)
        try:
            comparison = compare_source_backed_phase_candidates(
                source,
                (reference,) + tuple(
                    entry.candidate for entry in entries
                ),
                phase_rule,
            )
        except (TypeError, ValueError):
            unassessed_reasons.append(
                "family_phase_comparison_invalid"
            )
            continue
        competitions.append(_comparison_record(
            "complete_family_same_occurrence",
            path_ids,
            comparison,
        ))
        if not (
                comparison["status"] == "selected"
                and comparison["winner_candidate_index"] == 0):
            ambiguous_reasons.append(
                "family_phase_{}".format(comparison["reason"])
            )

    if unassessed_reasons:
        status = BOUNDARY_UNASSESSED
        final_reasons = unassessed_reasons + ambiguous_reasons
    elif ambiguous_reasons:
        status = BOUNDARY_AMBIGUOUS
        final_reasons = ambiguous_reasons
    else:
        status = BOUNDARY_EXACT
        final_reasons = ("source_phase_and_boundary_unambiguous",)
    return _result(
        status,
        final_reasons,
        source_sha256,
        proposed_geometry,
        alternatives,
        competitions,
    )


__all__ = [
    "BOUNDARY_AMBIGUOUS",
    "BOUNDARY_EXACT",
    "BOUNDARY_UNASSESSED",
    "BoundaryPhaseAlternative",
    "BoundaryPhaseCompetition",
    "BoundaryPhaseCompletenessEvidence",
    "BoundaryPhaseLineageEntry",
    "KNOWN_LINEAGE_DERIVATIONS",
    "SourcePhaseBoundaryAttestation",
    "attest_accepted_core_boundary_phase",
]
