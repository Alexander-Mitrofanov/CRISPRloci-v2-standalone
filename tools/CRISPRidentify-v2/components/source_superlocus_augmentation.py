"""Conservative source-superlocus and fixed-point augmentation primitive.

The primitive is deliberately isolated from scoring and pipeline categories.
It validates every legacy candidate against the complete source sequence,
freezes deterministic superloci from legacy geometry, and then applies the
existing one-step terminal, isolated weak-terminal, and exact-shared-repeat
primitives for at most four rounds. Originals that are not explicitly eligible
remain collision obstacles, but they never propose repeats, enter graph unions,
or appear in the returned candidates.

A superlocus that still changes in its bounded confirmation remains globally
nonconverged. One generated path may nevertheless remain output- and
architecture-eligible when it is the unique exact contiguous envelope of the
complete bounded and confirmation lineage and an exhaustive finite source-only
terminal closure proves that no descendant can leave that envelope. Shifted,
incomparable, disconnected, behavior-changing, and outward-growing paths all
retain the original superlocus-wide quarantine.

All coordinates are zero-based and end-exclusive.  Superlocus proximity does
not itself authorize an array union: only the shared-repeat primitive may
create a union, and that primitive requires exact repeat intervals and every
consecutive union edge to occur in an input path.
"""

from collections import Counter, defaultdict
import hashlib
from itertools import combinations, product
import math
import statistics

from components.components_detection import (
    FUZZY_SOURCE_PATH_SCHEMA,
    FuzzySearch,
    FuzzySourceOccurrence,
    FuzzySourcePathProvenance,
)
from components.components_detection_refinement import CrisprCandidate
from components.components_evaluated_arrays_enhancement import (
    apply_source_backed_terminal_repeat_proposals,
    apply_source_backed_weak_terminal_repeat_proposal,
    source_backed_terminal_repeat_proposals,
    source_backed_weak_terminal_repeat_proposals,
    validate_candidate_source_reconstruction,
)
from components.module_evaluated_arrays_enhancement import (
    build_source_backed_terminal_candidate_variant_records,
)
from components.shared_repeat_augmentation import (
    augment_source_backed_shared_repeat_candidates,
)


MAXIMUM_SUPERLOCUS_GAP = 140
MINIMUM_CONSENSUS_CONTAINMENT_IDENTITY = 0.80
MINIMUM_CONSENSUS_LENGTH_RATIO = 2.0 / 3.0
MINIMUM_PERIOD_TOLERANCE = 3.0
PERIOD_TOLERANCE_FRACTION = 0.10
MAXIMUM_ENHANCEMENT_ROUNDS = 4
PATH_LOCAL_CONVERGENCE_PROOF_SCHEMA = (
    "source-superlocus-path-local-convergence-proof"
)
PATH_LOCAL_CONVERGENCE_PROOF_VERSION = "1.0.0"


def _stable_value_key(value):
    return type(value).__name__, repr(value)


def _candidate_repeat_path(candidate):
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts,
            candidate.list_repeats,
        )
    )


def _candidate_geometry(candidate):
    return _candidate_repeat_path(candidate)


def _candidate_interval(candidate):
    path = _candidate_repeat_path(candidate)
    return path[0][0], path[-1][1]


def _candidate_sort_key(candidate):
    return (
        _candidate_geometry(candidate),
        str(getattr(candidate, "consensus", "")),
        tuple(candidate.list_repeats),
        tuple(candidate.list_spacers),
        tuple(getattr(candidate, "list_repeats_gaped", ())),
    )


def _source_geometry_key(candidate):
    return _candidate_geometry(candidate)


def _normalized_consensus(candidate):
    consensus = getattr(candidate, "consensus", None)
    if not isinstance(consensus, str) or not consensus:
        raise ValueError("every candidate must have a non-empty consensus")
    return (
        consensus.upper()
        .replace(" ", "")
        .replace("-", "")
        .replace("+", "")
    )


def _ungapped_containment_identity(first, second):
    if not first or not second:
        return 0.0
    if len(first) <= len(second):
        shorter, longer = first, second
    else:
        shorter, longer = second, first
    best_matches = max(
        sum(left == right for left, right in zip(
            shorter,
            longer[offset:offset + len(shorter)],
        ))
        for offset in range(len(longer) - len(shorter) + 1)
    )
    return best_matches / float(len(shorter))


def _median_period(candidate):
    starts = tuple(candidate.list_repeat_starts)
    if len(starts) < 2:
        return None
    return float(statistics.median(
        right - left for left, right in zip(starts, starts[1:])
    ))


def _interval_gap(first, second):
    if first[1] < second[0]:
        return second[0] - first[1]
    if second[1] < first[0]:
        return first[0] - second[1]
    return 0


def _period_tolerance(first_period, second_period):
    return max(
        MINIMUM_PERIOD_TOLERANCE,
        PERIOD_TOLERANCE_FRACTION * max(first_period, second_period),
    )


def _next_repeat_cadence_match(first, second, tolerance):
    matches = []
    first_period = _median_period(first)
    second_period = _median_period(second)
    for source_start in first.list_repeat_starts:
        expected = source_start + first_period
        for observed in second.list_repeat_starts:
            difference = abs(observed - expected)
            if difference <= tolerance:
                matches.append((
                    difference,
                    "first_to_second",
                    source_start,
                    expected,
                    observed,
                ))
    for source_start in second.list_repeat_starts:
        expected = source_start + second_period
        for observed in first.list_repeat_starts:
            difference = abs(observed - expected)
            if difference <= tolerance:
                matches.append((
                    difference,
                    "second_to_first",
                    source_start,
                    expected,
                    observed,
                ))
    return min(matches) if matches else None


def _candidate_pair_link_evidence(first, second):
    first_path = _candidate_repeat_path(first)
    second_path = _candidate_repeat_path(second)
    shared_intervals = tuple(sorted(set(first_path).intersection(second_path)))
    if shared_intervals:
        return {
            "mechanism": "exact_shared_repeat_interval",
            "shared_repeat_intervals": shared_intervals,
            "first_geometry": first_path,
            "second_geometry": second_path,
        }

    gap = _interval_gap(
        _candidate_interval(first),
        _candidate_interval(second),
    )
    if gap > MAXIMUM_SUPERLOCUS_GAP:
        return None

    first_consensus = _normalized_consensus(first)
    second_consensus = _normalized_consensus(second)
    length_ratio = min(len(first_consensus), len(second_consensus)) / float(
        max(len(first_consensus), len(second_consensus))
    )
    if length_ratio < MINIMUM_CONSENSUS_LENGTH_RATIO:
        return None
    containment_identity = _ungapped_containment_identity(
        first_consensus,
        second_consensus,
    )
    if containment_identity < MINIMUM_CONSENSUS_CONTAINMENT_IDENTITY:
        return None

    first_period = _median_period(first)
    second_period = _median_period(second)
    if first_period is None or second_period is None:
        return None
    tolerance = _period_tolerance(first_period, second_period)
    period_difference = abs(first_period - second_period)
    if period_difference > tolerance:
        return None
    cadence_match = _next_repeat_cadence_match(first, second, tolerance)
    if cadence_match is None:
        return None

    return {
        "mechanism": "compatible_source_topology",
        "gap": gap,
        "consensus_containment_identity": containment_identity,
        "consensus_length_ratio": length_ratio,
        "first_median_period": first_period,
        "second_median_period": second_period,
        "period_difference": period_difference,
        "period_tolerance": tolerance,
        "cadence_match": cadence_match,
        "first_geometry": first_path,
        "second_geometry": second_path,
    }


def _validated_source_path_family_keys(full_dna, candidate):
    """Return immutable fuzzy-parent identities after full source replay.

    Refined candidates no longer carry the ``FuzzySearch`` validation method,
    so downstream family construction must not trust copied namedtuples. A
    provenance record authorizes a family-context edge only after the complete
    parent path (including portions outside this local child), its spacers,
    deterministic family identifier, and the child's selected span all replay
    exactly against the original source sequence.
    """

    provenance = getattr(candidate, "source_path_provenance", None)
    if provenance is None:
        return ()
    if (
            not isinstance(provenance, tuple)
            or not provenance
            or any(
                not isinstance(record, FuzzySourcePathProvenance)
                for record in provenance
            )):
        raise ValueError("source path provenance is malformed")
    try:
        if len(set(provenance)) != len(provenance):
            raise ValueError("source path provenance is duplicated")
    except TypeError as error:
        raise ValueError("source path provenance is not immutable") from error

    candidate_starts = tuple(candidate.list_repeat_starts)
    candidate_repeats = tuple(candidate.list_repeats)
    candidate_spacers = tuple(candidate.list_spacers)
    family_keys = []
    for record in provenance:
        if record.schema != FUZZY_SOURCE_PATH_SCHEMA:
            raise ValueError("source path provenance schema is unsupported")
        if (
                not isinstance(record.repeat_candidate, str)
                or not record.repeat_candidate
                or not isinstance(record.parent_geometry, tuple)
                or not record.parent_geometry
                or not isinstance(record.parent_spacers, tuple)
                or len(record.parent_spacers)
                != len(record.parent_geometry) - 1):
            raise ValueError("source path provenance parent is malformed")

        previous_end = None
        parent_occurrence_ids = []
        for occurrence in record.parent_geometry:
            if (
                    not isinstance(occurrence, FuzzySourceOccurrence)
                    or not isinstance(occurrence.start, int)
                    or isinstance(occurrence.start, bool)
                    or not isinstance(occurrence.end, int)
                    or isinstance(occurrence.end, bool)
                    or occurrence.start < 0
                    or occurrence.end <= occurrence.start
                    or occurrence.end > len(full_dna)
                    or not isinstance(occurrence.sequence, str)
                    or not occurrence.sequence
                    or occurrence.end
                    != occurrence.start + len(occurrence.sequence)
                    or occurrence.occurrence_id
                    != FuzzySearch._source_occurrence_id(
                        occurrence.start, occurrence.end
                    )
                    or (
                        previous_end is not None
                        and occurrence.start < previous_end
                    )
                    or full_dna[occurrence.start:occurrence.end]
                    != occurrence.sequence):
                raise ValueError(
                    "source path provenance parent does not reconstruct "
                    "the original sequence"
                )
            parent_occurrence_ids.append(occurrence.occurrence_id)
            previous_end = occurrence.end
        if len(set(parent_occurrence_ids)) != len(parent_occurrence_ids):
            raise ValueError(
                "source path provenance occurrence IDs are duplicated"
            )

        for index, spacer in enumerate(record.parent_spacers):
            left = record.parent_geometry[index]
            right = record.parent_geometry[index + 1]
            if (
                    not isinstance(spacer, str)
                    or full_dna[left.end:right.start] != spacer):
                raise ValueError(
                    "source path provenance parent spacers do not "
                    "reconstruct the original sequence"
                )

        expected_family_id = FuzzySearch._source_family_id(
            record.repeat_candidate,
            record.parent_geometry,
            record.parent_spacers,
        )
        if record.family_id != expected_family_id:
            raise ValueError(
                "source path provenance family identity is inconsistent"
            )

        span = record.source_index_span
        if (
                not isinstance(span, tuple)
                or len(span) != 2
                or any(
                    not isinstance(index, int) or isinstance(index, bool)
                    for index in span
                )
                or span[0] < 0
                or span[1] > len(record.parent_geometry)
                or span[1] <= span[0]):
            raise ValueError("source path provenance span is invalid")
        selected_geometry = record.parent_geometry[span[0]:span[1]]
        selected_spacers = record.parent_spacers[
            span[0]:span[1] - 1
        ]
        if (
                tuple(record.occurrence_ids)
                != tuple(
                    occurrence.occurrence_id
                    for occurrence in selected_geometry
                )
                or candidate_starts
                != tuple(
                    occurrence.start for occurrence in selected_geometry
                )
                or candidate_repeats
                != tuple(
                    occurrence.sequence for occurrence in selected_geometry
                )
                or candidate_spacers != selected_spacers):
            raise ValueError(
                "source path provenance does not describe this candidate"
            )
        if record.derivation not in (
                "complete_fuzzy_search",
                "local_chain_family_parent",
                "maximal_local_chain"):
            raise ValueError("source path provenance derivation is unknown")
        if (
                record.derivation in (
                    "complete_fuzzy_search",
                    "local_chain_family_parent",
                )
                and span != (0, len(record.parent_geometry))):
            raise ValueError(
                "complete source path provenance must span its parent"
            )

        family_keys.append((
            record.family_id,
            record.parent_geometry,
            record.parent_spacers,
        ))
    return tuple(sorted(set(family_keys)))


def _source_path_family_link_evidence(
        first_candidate, first_families,
        second_candidate, second_families):
    """Return an edge only for one identical, source-validated parent."""

    shared_families = set(first_families).intersection(second_families)
    if not shared_families:
        return None
    family_id, parent_geometry, parent_spacers = min(shared_families)
    return {
        "mechanism": "validated_source_path_family",
        "family_id": family_id,
        "parent_geometry": parent_geometry,
        "parent_spacers": parent_spacers,
        "first_geometry": _candidate_repeat_path(first_candidate),
        "second_geometry": _candidate_repeat_path(second_candidate),
    }


def _link_evidence_sort_key(evidence):
    priority = {
        "exact_shared_repeat_interval": 0,
        "compatible_source_topology": 1,
        "validated_source_path_family": 2,
    }
    return (
        priority[evidence["mechanism"]],
        evidence["first_geometry"],
        evidence["second_geometry"],
    )


def _legacy_locus_link_evidence(first, second):
    evidence = []
    for first_candidate in first["eligible_candidates"]:
        for second_candidate in second["eligible_candidates"]:
            pair_evidence = _candidate_pair_link_evidence(
                first_candidate,
                second_candidate,
            )
            if pair_evidence is not None:
                evidence.append(pair_evidence)
    for first_candidate in first["candidates"]:
        first_families = first["source_path_families"][
            id(first_candidate)
        ]
        for second_candidate in second["candidates"]:
            pair_evidence = _source_path_family_link_evidence(
                first_candidate,
                first_families,
                second_candidate,
                second["source_path_families"][id(second_candidate)],
            )
            if pair_evidence is not None:
                evidence.append(pair_evidence)
    return min(evidence, key=_link_evidence_sort_key) if evidence else None


def _connected_components(adjacency):
    unseen = set(range(len(adjacency)))
    components = []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        stack = [first]
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            neighbours = sorted(
                adjacency[current].intersection(unseen),
                reverse=True,
            )
            for neighbour in neighbours:
                unseen.remove(neighbour)
                stack.append(neighbour)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _exact_shared_components(candidates):
    interval_owners = defaultdict(list)
    for candidate_index, candidate in enumerate(candidates):
        for interval in set(_candidate_repeat_path(candidate)):
            interval_owners[interval].append(candidate_index)
    adjacency = [set() for _candidate in candidates]
    for owners in interval_owners.values():
        if len(owners) < 2:
            continue
        for left_index, left in enumerate(owners[:-1]):
            for right in owners[left_index + 1:]:
                adjacency[left].add(right)
                adjacency[right].add(left)
    return _connected_components(adjacency)


def _path_identifier(superlocus_id, geometry):
    payload = repr((superlocus_id, geometry)).encode("utf-8")
    return "path-" + hashlib.sha256(payload).hexdigest()


class _MultipleWeakLineageError(RuntimeError):
    pass


def _derived_weak_interval_sets(
        state, geometry, stage, parent_path_ids, details):
    """Return viable per-derivation weak intervals for one geometry.

    An empty interval set is an ordinary/strong proof. A singleton is a proof
    containing one weakly added repeat and therefore closes further terminal
    completion and requires final source-coding assessment. Shared unions
    combine parent proofs, but an interval weak in one parent is discharged if
    another selected parent proof independently contains it as ordinary.
    """

    parent_by_id = {
        entry["path_id"]: entry for entry in state["paths"].values()
    }
    missing_parent_ids = set(parent_path_ids).difference(parent_by_id)
    if missing_parent_ids:
        raise RuntimeError("a lineage derivation references an unknown parent")
    parent_entries = tuple(
        parent_by_id[parent_id] for parent_id in parent_path_ids
    )
    if not parent_entries:
        return ((),)

    terminal_stage = stage in (
        "terminal_augmentation",
        "weak_terminal_augmentation",
    )
    parent_options = []
    for parent in parent_entries:
        options = tuple(parent["lineage_weak_interval_sets"])
        if terminal_stage:
            options = tuple(option for option in options if not option)
        if not options:
            return ()
        parent_options.append(options)

    weak_stage = stage == "weak_terminal_augmentation"
    weak_added_intervals = set()
    if weak_stage:
        proposals = tuple(
            proposal for proposal in details.get("proposals", ())
            if isinstance(proposal, dict)
            and proposal.get("tier") == "weak"
        )
        if len(proposals) != 1:
            raise RuntimeError(
                "a weak terminal derivation requires exactly one weak proposal"
            )
        weak_added_intervals.add(tuple(proposals[0]["repeat_interval"]))

    geometry_set = set(geometry)
    viable = set()
    rejected_multiple_weak = False
    for combination in product(*parent_options):
        weak_intervals = set()
        trusted_intervals = set()
        parent_geometry = set()
        for parent, parent_weak in zip(parent_entries, combination):
            parent_intervals = set(parent["source_geometry"])
            parent_weak = set(parent_weak)
            weak_intervals.update(parent_weak)
            trusted_intervals.update(parent_intervals.difference(parent_weak))
            parent_geometry.update(parent_intervals)

        added_intervals = geometry_set.difference(parent_geometry)
        if weak_stage:
            if not weak_added_intervals.issubset(added_intervals):
                raise RuntimeError(
                    "a weak proposal must add a new source interval"
                )
            weak_intervals.update(weak_added_intervals)
            trusted_intervals.update(
                added_intervals.difference(weak_added_intervals)
            )
        else:
            trusted_intervals.update(added_intervals)

        weak_intervals.difference_update(trusted_intervals)
        if not weak_intervals.issubset(geometry_set):
            raise RuntimeError(
                "lineage weak provenance is outside the derived geometry"
            )
        if len(weak_intervals) > 1:
            rejected_multiple_weak = True
            continue
        viable.add(tuple(sorted(weak_intervals)))

    if not viable and rejected_multiple_weak:
        raise _MultipleWeakLineageError(
            "a source path cannot contain two weakly added repeats"
        )
    return tuple(sorted(viable))


def _add_path(
    state,
    candidate,
    stage,
    round_number,
    parent_path_ids=(),
    details=None,
):
    geometry = _source_geometry_key(candidate)
    existing = state["paths"].get(geometry)
    details = {} if details is None else dict(details)
    parent_path_ids = tuple(sorted(set(parent_path_ids)))
    derived_weak_sets = _derived_weak_interval_sets(
        state,
        geometry,
        stage,
        parent_path_ids,
        details,
    )
    if not derived_weak_sets:
        raise RuntimeError("a terminal derivation has no open parent proof")
    if existing is not None:
        merged_parent_ids = tuple(sorted(
            set(existing["parent_path_ids"])
            .union(parent_path_ids)
            .difference((existing["path_id"],))
        ))
        existing["parent_path_ids"] = merged_parent_ids
        lineage_weak_interval_sets = tuple(sorted(set(
            existing["lineage_weak_interval_sets"]
        ).union(derived_weak_sets)))
        existing["lineage_weak_interval_sets"] = lineage_weak_interval_sets
        all_derivations_are_weak = all(lineage_weak_interval_sets)
        existing["terminal_completion_closed"] = all_derivations_are_weak
        existing["requires_source_coding_context"] = (
            all_derivations_are_weak
        )
        derivation = {
            "stage": stage,
            "round": round_number,
            "parent_path_ids": tuple(
                parent_id for parent_id in parent_path_ids
                if parent_id != existing["path_id"]
            ),
            "details": details,
        }
        if derivation not in existing["derivations"]:
            existing["derivations"] = (
                existing["derivations"] + (derivation,)
            )
        return existing, False
    validate_candidate_source_reconstruction(
        state["full_dna"],
        candidate.list_repeats,
        candidate.list_repeat_starts,
        candidate.list_spacers,
    )
    path_id = _path_identifier(state["superlocus_id"], geometry)
    normalized_parent_ids = tuple(
        parent_id for parent_id in parent_path_ids
        if parent_id != path_id
    )
    derivation = {
        "stage": stage,
        "round": round_number,
        "parent_path_ids": normalized_parent_ids,
        "details": details,
    }
    all_derivations_are_weak = all(derived_weak_sets)
    entry = {
        "path_id": path_id,
        "superlocus_id": state["superlocus_id"],
        "stage": stage,
        "round": round_number,
        "parent_path_ids": normalized_parent_ids,
        "source_geometry": geometry,
        "candidate": candidate,
        "details": details,
        "derivations": (derivation,),
        "lineage_weak_interval_sets": derived_weak_sets,
        "terminal_completion_closed": all_derivations_are_weak,
        "requires_source_coding_context": all_derivations_are_weak,
    }
    state["paths"][geometry] = entry
    state["candidate_path_ids"][id(candidate)] = path_id
    return entry, True


def _unique_sorted_candidates(candidates):
    by_runtime_id = {}
    for candidate in candidates:
        by_runtime_id.setdefault(id(candidate), candidate)
    return tuple(sorted(by_runtime_id.values(), key=_candidate_sort_key))


def _obstacle_mapping(states):
    mapping = {}
    eligible_runtime_ids = set()
    for state in states:
        path_candidates = [
            entry["candidate"] for entry in state["paths"].values()
        ]
        candidates = _unique_sorted_candidates(
            list(state["legacy_obstacles"]) + path_candidates
        )
        mapping[state["superlocus_id"]] = list(candidates)
        eligible_runtime_ids.update(
            id(entry["candidate"])
            for entry in state["paths"].values()
            if not entry.get("terminal_completion_closed", False)
        )
    return mapping, eligible_runtime_ids


def _sum_stats(target, source, prefix, fields):
    for field in fields:
        target[prefix + field] += source[field]


def _normalize_terminal_proposals(proposals):
    return tuple({
        "tier": proposal.get("tier", "strong"),
        "side": proposal["side"],
        "repeat_interval": tuple(proposal["repeat_interval"]),
        "spacer_length": proposal.get("spacer_length"),
        "substitutions": proposal.get("substitutions"),
        "repeat_identity": proposal.get("repeat_identity"),
        "exact_anchor_length": proposal.get("exact_anchor_length"),
        "source_phase_score": proposal.get("source_phase_score"),
        "source_phase_winner_margin": proposal.get(
            "source_phase_winner_margin"
        ),
        "translated_tandem_identity": proposal.get(
            "translated_tandem_identity"
        ),
        "requires_source_coding_context": proposal.get(
            "requires_source_coding_context", False
        ),
    } for proposal in proposals)


def _state_semantic_signature(state):
    """Return the deterministic state that can affect a later round.

    Derivation audit records deliberately do not participate: the same proof
    recorded in a later round has no construction consequence. Geometry,
    parent connectivity, viable weak-proof alternatives, and terminal-open
    status do participate because a duplicate derivation can change them
    without adding a geometry.
    """

    return tuple(
        (
            geometry,
            tuple(entry["parent_path_ids"]),
            tuple(entry["lineage_weak_interval_sets"]),
            bool(entry["terminal_completion_closed"]),
            bool(entry["requires_source_coding_context"]),
        )
        for geometry, entry in sorted(state["paths"].items())
    )


def _clone_state_for_confirmation(state):
    """Clone mutable lineage state while sharing immutable candidates/source."""

    cloned = dict(state)
    cloned["candidate_path_ids"] = dict(state["candidate_path_ids"])
    cloned["paths"] = {
        geometry: {
            **entry,
            "details": dict(entry["details"]),
            "derivations": tuple({
                **derivation,
                "details": dict(derivation["details"]),
            } for derivation in entry["derivations"]),
        }
        for geometry, entry in state["paths"].items()
    }
    return cloned


def _path_is_exact_contiguous_subpath(shorter, longer):
    """Return whether one exact repeat-interval path occurs in another."""

    shorter = tuple(shorter)
    longer = tuple(longer)
    if not shorter or len(shorter) > len(longer):
        return False
    return any(
        longer[offset:offset + len(shorter)] == shorter
        for offset in range(len(longer) - len(shorter) + 1)
    )


def _entry_behavior_signature(entry):
    """Return path state that can change its future construction behavior."""

    return (
        tuple(entry["lineage_weak_interval_sets"]),
        bool(entry["terminal_completion_closed"]),
        bool(entry["requires_source_coding_context"]),
    )


def _lineage_component_path_ids(paths, root_path_id):
    """Return one undirected lineage component from a path mapping."""

    by_path_id = {
        entry["path_id"]: entry for entry in paths.values()
    }
    if root_path_id not in by_path_id:
        raise RuntimeError("path-local convergence root is unknown")
    adjacency = {
        path_id: set() for path_id in by_path_id
    }
    for path_id, entry in by_path_id.items():
        for parent_path_id in entry["parent_path_ids"]:
            if parent_path_id not in by_path_id:
                raise RuntimeError(
                    "path-local convergence lineage has an unknown parent"
                )
            adjacency[path_id].add(parent_path_id)
            adjacency[parent_path_id].add(path_id)

    component = {root_path_id}
    stack = [root_path_id]
    while stack:
        current = stack.pop()
        for neighbour in sorted(adjacency[current].difference(component)):
            component.add(neighbour)
            stack.append(neighbour)
    return frozenset(component)


def _candidate_terminal_behavior_key(candidate):
    """Return all candidate state that can affect another terminal step."""

    return (
        _candidate_geometry(candidate),
        str(candidate.consensus),
        tuple(candidate.list_repeats_gaped),
    )


def _validated_terminal_proposals(full_dna, proposals, tier):
    proposals = tuple(proposals)
    by_side = {}
    for proposal in proposals:
        if not isinstance(proposal, dict):
            raise RuntimeError("terminal closure proposal is malformed")
        side = proposal.get("side")
        interval = proposal.get("repeat_interval")
        if (
                side not in ("left", "right")
                or side in by_side
                or not isinstance(interval, tuple)
                or len(interval) != 2
                or any(
                    not isinstance(value, int) or isinstance(value, bool)
                    for value in interval
                )
                or interval[0] < 0
                or interval[1] <= interval[0]
                or interval[1] > len(full_dna)):
            raise RuntimeError("terminal closure proposal is malformed")
        by_side[side] = interval
    if tier == "weak" and proposals and len(proposals) != 1:
        raise RuntimeError(
            "weak terminal closure produced multiple proposals"
        )
    return proposals


def _terminal_successors_without_collision_suppression(full_dna, candidate):
    """Return every candidate that collision filtering could construct.

    The production transaction computes strong proposals first and may retain
    any non-empty subset after peer-collision filtering.  Its strong builder
    deliberately carries the parent's alignment/consensus representation into
    the rebuilt candidate, while the weak builder uses its own source-only
    reconstruction.  Rebuilding one combined geometry canonically is therefore
    not an over-approximation of later behavior.  This proof helper calls the
    exact production builders for every admissible strong subset instead.
    """

    strong = _validated_terminal_proposals(
        full_dna,
        source_backed_terminal_repeat_proposals(full_dna, candidate),
        "strong",
    )
    successors = []
    if strong:
        for subset_size in range(1, len(strong) + 1):
            for subset in combinations(strong, subset_size):
                successor = apply_source_backed_terminal_repeat_proposals(
                    full_dna, candidate, subset
                )
                successors.append((successor, "strong"))
    else:
        weak = _validated_terminal_proposals(
            full_dna,
            source_backed_weak_terminal_repeat_proposals(
                full_dna, candidate
            ),
            "weak",
        )
        if weak:
            successor = apply_source_backed_weak_terminal_repeat_proposal(
                full_dna, candidate, weak[0]
            )
            successors.append((successor, "weak"))

    unique = {}
    for successor, tier in successors:
        validate_candidate_source_reconstruction(
            full_dna,
            successor.list_repeats,
            successor.list_repeat_starts,
            successor.list_spacers,
        )
        unique.setdefault(
            (_candidate_terminal_behavior_key(successor), tier),
            (successor, tier),
        )
    return tuple(
        unique[key] for key in sorted(unique)
    )


def _contained_terminal_closure(full_dna, confirmation_paths, envelope):
    """Prove all future terminal descendants remain inside one envelope.

    Every confirmation path is already required to be an exact contiguous
    subpath of ``envelope``. Strong descendants remain terminal-open; weak
    descendants are terminal-closed. Multiple alignment/consensus behavior
    states may share one geometry, so they are traversed separately. Every
    transition strictly adds a repeat, and an explicit finite state budget
    rejects rather than certifies an unexpectedly broad closure.
    """

    maximum_geometry_count = (
        len(envelope) * (len(envelope) + 1) // 2
    )
    behavior_state_budget = max(
        maximum_geometry_count,
        maximum_geometry_count * len(confirmation_paths),
    )
    candidates_by_state = {
        _candidate_terminal_behavior_key(entry["candidate"]): entry["candidate"]
        for entry in confirmation_paths.values()
    }
    initial_geometries = set(confirmation_paths)
    open_states = {
        _candidate_terminal_behavior_key(entry["candidate"])
        for entry in confirmation_paths.values()
        if not entry["terminal_completion_closed"]
    }
    pending = set(open_states)
    explored_states = set()
    explored_geometries = set()
    generated_geometries = set()
    generated_state_count = 0
    transition_count = 0
    weak_closed_states = set()

    while pending:
        state_key = min(pending)
        pending.remove(state_key)
        if state_key in explored_states:
            continue
        explored_states.add(state_key)
        geometry = state_key[0]
        explored_geometries.add(geometry)
        results = _terminal_successors_without_collision_suppression(
            full_dna,
            candidates_by_state[state_key],
        )
        for successor, tier in results:
            successor_geometry = _candidate_geometry(successor)
            transition_count += 1
            if (
                    len(successor_geometry) <= len(geometry)
                    or not _path_is_exact_contiguous_subpath(
                        geometry, successor_geometry
                    )):
                raise RuntimeError(
                    "terminal closure successor does not strictly extend parent"
                )
            if not _path_is_exact_contiguous_subpath(
                    successor_geometry, envelope):
                return {
                    "status": "rejected",
                    "reason": "terminal_successor_escaped_envelope",
                    "terminal_closure_explored_geometry_count": len(
                        explored_geometries
                    ),
                    "terminal_closure_explored_behavior_state_count": len(
                        explored_states
                    ),
                    "terminal_closure_generated_geometry_count": len(
                        generated_geometries
                    ),
                    "terminal_closure_generated_behavior_state_count": (
                        generated_state_count
                    ),
                    "terminal_closure_transition_count": transition_count,
                    "terminal_closure_weak_closed_geometry_count": len({
                        key[0] for key in weak_closed_states
                    }),
                    "terminal_closure_weak_closed_behavior_state_count": len(
                        weak_closed_states
                    ),
                    "terminal_closure_behavior_state_budget": (
                        behavior_state_budget
                    ),
                    "escaped_geometry": successor_geometry,
                }
            successor_key = _candidate_terminal_behavior_key(successor)
            if successor_key not in candidates_by_state:
                candidates_by_state[successor_key] = successor
                generated_state_count += 1
                if successor_geometry not in initial_geometries:
                    generated_geometries.add(successor_geometry)
                if len(candidates_by_state) > behavior_state_budget:
                    return {
                        "status": "rejected",
                        "reason": "terminal_closure_state_budget_exceeded",
                        "terminal_closure_explored_geometry_count": len(
                            explored_geometries
                        ),
                        "terminal_closure_explored_behavior_state_count": len(
                            explored_states
                        ),
                        "terminal_closure_generated_geometry_count": len(
                            generated_geometries
                        ),
                        "terminal_closure_generated_behavior_state_count": (
                            generated_state_count
                        ),
                        "terminal_closure_transition_count": transition_count,
                        "terminal_closure_weak_closed_geometry_count": len({
                            key[0] for key in weak_closed_states
                        }),
                        "terminal_closure_weak_closed_behavior_state_count": len(
                            weak_closed_states
                        ),
                        "terminal_closure_behavior_state_budget": (
                            behavior_state_budget
                        ),
                        "escaped_geometry": None,
                    }
            if tier == "strong":
                if successor_key not in open_states:
                    open_states.add(successor_key)
                    pending.add(successor_key)
            else:
                weak_closed_states.add(successor_key)

    return {
        "status": "converged",
        "reason": "finite_contained_terminal_closure",
        "terminal_closure_explored_geometry_count": len(explored_geometries),
        "terminal_closure_explored_behavior_state_count": len(
            explored_states
        ),
        "terminal_closure_generated_geometry_count": len(
            generated_geometries
        ),
        "terminal_closure_generated_behavior_state_count": (
            generated_state_count
        ),
        "terminal_closure_transition_count": transition_count,
        "terminal_closure_weak_closed_geometry_count": len({
            key[0] for key in weak_closed_states
        }),
        "terminal_closure_weak_closed_behavior_state_count": len(
            weak_closed_states
        ),
        "terminal_closure_behavior_state_budget": behavior_state_budget,
        "escaped_geometry": None,
    }


_PATH_LOCAL_WITNESS_KEYS = frozenset((
    "path_id",
    "source_geometry",
    "parent_path_ids",
    "list_repeats_gaped",
    "consensus",
    "lineage_weak_interval_sets",
    "terminal_completion_closed",
    "requires_source_coding_context",
))


def _path_local_source_sha256(full_dna):
    if not isinstance(full_dna, str):
        raise ValueError("path-local proof source must be a string")
    try:
        encoded = full_dna.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("path-local proof source must be ASCII DNA") from error
    return hashlib.sha256(encoded).hexdigest()


def _path_local_witness_sequence(value, label):
    if isinstance(value, (str, bytes)):
        raise ValueError("{} must be a sequence".format(label))
    try:
        return tuple(value)
    except TypeError as error:
        raise ValueError("{} must be a sequence".format(label)) from error


def _path_local_witness_geometry(full_dna, value, label):
    raw = _path_local_witness_sequence(value, label)
    if not raw:
        raise ValueError("{} must contain repeat intervals".format(label))
    geometry = []
    previous_end = None
    for index, raw_interval in enumerate(raw):
        interval = _path_local_witness_sequence(
            raw_interval, "{}[{}]".format(label, index)
        )
        if (
                len(interval) != 2
                or any(
                    not isinstance(item, int) or isinstance(item, bool)
                    for item in interval
                )
                or interval[0] < 0
                or interval[1] <= interval[0]
                or interval[1] > len(full_dna)
                or (
                    previous_end is not None
                    and interval[0] < previous_end
                )):
            raise ValueError("{} contains an invalid interval".format(label))
        geometry.append(interval)
        previous_end = interval[1]
    return tuple(geometry)


def _path_local_witness_weak_sets(value, geometry, label):
    raw_sets = _path_local_witness_sequence(value, label)
    if not raw_sets:
        raise ValueError("{} must contain a derivation state".format(label))
    geometry_set = set(geometry)
    weak_sets = []
    for option_index, raw_option in enumerate(raw_sets):
        raw_intervals = _path_local_witness_sequence(
            raw_option, "{}[{}]".format(label, option_index)
        )
        intervals = []
        for interval_index, raw_interval in enumerate(raw_intervals):
            interval = _path_local_witness_sequence(
                raw_interval,
                "{}[{}][{}]".format(
                    label, option_index, interval_index
                ),
            )
            if (
                    len(interval) != 2
                    or any(
                        not isinstance(item, int) or isinstance(item, bool)
                        for item in interval
                    )
                    or tuple(interval) not in geometry_set):
                raise ValueError(
                    "{} contains a weak interval outside its geometry".format(
                        label
                    )
                )
            intervals.append(tuple(interval))
        if len(intervals) > 1:
            raise ValueError(
                "{} contains multiple weak terminal intervals".format(label)
            )
        weak_sets.append(tuple(sorted(intervals)))
    normalized = tuple(sorted(set(weak_sets)))
    if tuple(weak_sets) != normalized:
        raise ValueError("{} is not canonical".format(label))
    return normalized


def _path_local_entry_from_witness(
        full_dna, superlocus_id, raw_record, label):
    if not isinstance(raw_record, dict):
        raise ValueError("{} must be a mapping".format(label))
    if set(raw_record) != _PATH_LOCAL_WITNESS_KEYS:
        raise ValueError("{} fields are incomplete or unknown".format(label))

    geometry = _path_local_witness_geometry(
        full_dna, raw_record["source_geometry"], label + ".source_geometry"
    )
    path_id = raw_record["path_id"]
    if (
            not isinstance(path_id, str)
            or not path_id
            or path_id != _path_identifier(superlocus_id, geometry)):
        raise ValueError("{} path ID is stale".format(label))

    parent_path_ids = _path_local_witness_sequence(
        raw_record["parent_path_ids"], label + ".parent_path_ids"
    )
    if (
            any(
                not isinstance(parent_id, str) or not parent_id
                for parent_id in parent_path_ids
            )
            or tuple(sorted(set(parent_path_ids))) != parent_path_ids
            or path_id in parent_path_ids):
        raise ValueError("{} parent path IDs are invalid".format(label))

    gapped_repeats = _path_local_witness_sequence(
        raw_record["list_repeats_gaped"],
        label + ".list_repeats_gaped",
    )
    if any(
            not isinstance(repeat, str) or not repeat
            for repeat in gapped_repeats):
        raise ValueError("{} gapped repeats are invalid".format(label))
    starts = tuple(start for start, _end in geometry)
    repeats = tuple(full_dna[start:end] for start, end in geometry)
    spacers = tuple(
        full_dna[left_end:right_start]
        for (_left_start, left_end), (right_start, _right_end)
        in zip(geometry, geometry[1:])
    )
    candidate = CrisprCandidate(
        list_repeats=list(repeats),
        list_repeats_gaped=list(gapped_repeats),
        list_spacers=list(spacers),
        list_repeat_starts=list(starts),
    )
    validate_candidate_source_reconstruction(
        full_dna,
        candidate.list_repeats,
        candidate.list_repeat_starts,
        candidate.list_spacers,
    )
    if tuple(candidate.list_repeats_gaped) != gapped_repeats:
        raise ValueError("{} gapped repeats are not canonical".format(label))
    if (
            not isinstance(raw_record["consensus"], str)
            or not raw_record["consensus"]
            or candidate.consensus != raw_record["consensus"]):
        raise ValueError("{} consensus is stale".format(label))

    weak_sets = _path_local_witness_weak_sets(
        raw_record["lineage_weak_interval_sets"],
        geometry,
        label + ".lineage_weak_interval_sets",
    )
    terminal_closed = raw_record["terminal_completion_closed"]
    requires_coding = raw_record["requires_source_coding_context"]
    expected_closed = all(weak_sets)
    if (
            not isinstance(terminal_closed, bool)
            or not isinstance(requires_coding, bool)
            or terminal_closed != expected_closed
            or requires_coding != expected_closed):
        raise ValueError("{} weak-state flags are inconsistent".format(label))

    return geometry, {
        "path_id": path_id,
        "parent_path_ids": parent_path_ids,
        "source_geometry": geometry,
        "candidate": candidate,
        "lineage_weak_interval_sets": weak_sets,
        "terminal_completion_closed": terminal_closed,
        "requires_source_coding_context": requires_coding,
    }


def _path_local_state_from_witness(
        full_dna, superlocus_id, raw_records, label):
    records = _path_local_witness_sequence(raw_records, label)
    if not records:
        raise ValueError("{} must contain paths".format(label))
    paths = {}
    path_ids = set()
    for index, raw_record in enumerate(records):
        geometry, entry = _path_local_entry_from_witness(
            full_dna,
            superlocus_id,
            raw_record,
            "{}[{}]".format(label, index),
        )
        if geometry in paths or entry["path_id"] in path_ids:
            raise ValueError("{} contains duplicate paths".format(label))
        paths[geometry] = entry
        path_ids.add(entry["path_id"])
    if any(
            not set(entry["parent_path_ids"]).issubset(path_ids)
            for entry in paths.values()):
        raise ValueError("{} contains an unknown parent path".format(label))
    return {"superlocus_id": superlocus_id, "paths": paths}


def _path_local_state_witness(full_dna, superlocus_id, paths):
    records = []
    for geometry, entry in sorted(paths.items()):
        candidate = entry["candidate"]
        validate_candidate_source_reconstruction(
            full_dna,
            candidate.list_repeats,
            candidate.list_repeat_starts,
            candidate.list_spacers,
        )
        if geometry != _candidate_geometry(candidate):
            raise RuntimeError(
                "path-local witness candidate geometry is inconsistent"
            )
        record = {
            "path_id": entry["path_id"],
            "source_geometry": tuple(tuple(interval) for interval in geometry),
            "parent_path_ids": tuple(entry["parent_path_ids"]),
            "list_repeats_gaped": tuple(candidate.list_repeats_gaped),
            "consensus": candidate.consensus,
            "lineage_weak_interval_sets": tuple(
                tuple(tuple(interval) for interval in option)
                for option in entry["lineage_weak_interval_sets"]
            ),
            "terminal_completion_closed": entry[
                "terminal_completion_closed"
            ],
            "requires_source_coding_context": entry[
                "requires_source_coding_context"
            ],
        }
        _path_local_entry_from_witness(
            full_dna, superlocus_id, record, "generated_path_witness"
        )
        records.append(record)
    return tuple(records)


def _path_local_proof_value(value):
    """Return a type-strict, tuple-normalized value for exact comparison."""

    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise ValueError("path-local proof keys must be strings")
        return (
            "mapping",
            tuple(sorted(
                (key, _path_local_proof_value(item))
                for key, item in value.items()
            )),
        )
    if isinstance(value, (list, tuple)):
        return (
            "sequence",
            tuple(_path_local_proof_value(item) for item in value),
        )
    if value is None:
        return ("none", None)
    if isinstance(value, bool):
        return ("boolean", value)
    if isinstance(value, int):
        return ("integer", value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("path-local proof numbers must be finite")
        return ("float", value)
    if isinstance(value, str):
        return ("string", value)
    raise ValueError("path-local proof contains an unsupported value")


def _path_local_convergence_result(full_dna, state, confirmation_state):
    """Return a unique source-envelope proof or a fail-closed rejection.

    A bounded path can survive a nonconverged superlocus only when it already
    contains every bounded and confirmation path as an exact contiguous
    subpath, spans the complete lineage component before and after
    confirmation, preserves its behavior-affecting weak state, and the exact
    finite terminal closure cannot leave its envelope. Because the envelope
    is itself present, every exact-shared-repeat union in that closure is the
    same envelope; no score or iteration order participates.
    """

    bounded_paths = state["paths"]
    confirmation_paths = confirmation_state["paths"]
    bounded_geometries = tuple(sorted(bounded_paths))
    confirmation_geometries = tuple(sorted(confirmation_paths))
    envelope_geometries = tuple(
        geometry for geometry in bounded_geometries
        if all(
            _path_is_exact_contiguous_subpath(other, geometry)
            for other in confirmation_geometries
        )
    )
    base_evidence = {
        "coordinate_system": "zero_based_half_open",
        "bounded_candidate_count": len(bounded_paths),
        "confirmation_candidate_count": len(confirmation_paths),
        "confirmation_added_geometry_count": len(
            set(confirmation_paths).difference(bounded_paths)
        ),
        "confirmation_changed_existing_geometry_count": sum(
            (
                tuple(bounded_paths[geometry]["parent_path_ids"]),
                _entry_behavior_signature(bounded_paths[geometry]),
            )
            != (
                tuple(confirmation_paths[geometry]["parent_path_ids"]),
                _entry_behavior_signature(confirmation_paths[geometry]),
            )
            for geometry in set(bounded_paths).intersection(
                confirmation_paths
            )
        ),
        "candidate_selection": "unique_exact_contiguous_envelope",
    }
    if len(envelope_geometries) != 1:
        return {
            **base_evidence,
            "status": "rejected",
            "reason": "no_unique_exact_contiguous_envelope",
            "envelope_geometry": None,
        }

    envelope = envelope_geometries[0]
    bounded_entry = bounded_paths[envelope]
    confirmation_entry = confirmation_paths[envelope]
    validate_candidate_source_reconstruction(
        full_dna,
        bounded_entry["candidate"].list_repeats,
        bounded_entry["candidate"].list_repeat_starts,
        bounded_entry["candidate"].list_spacers,
    )
    if (
            _entry_behavior_signature(bounded_entry)
            != _entry_behavior_signature(confirmation_entry)):
        return {
            **base_evidence,
            "status": "rejected",
            "reason": "envelope_behavior_changed_in_confirmation",
            "envelope_geometry": envelope,
        }

    bounded_path_ids = frozenset(
        entry["path_id"] for entry in bounded_paths.values()
    )
    confirmation_path_ids = frozenset(
        entry["path_id"] for entry in confirmation_paths.values()
    )
    envelope_path_id = bounded_entry["path_id"]
    if (
            _lineage_component_path_ids(
                bounded_paths, envelope_path_id
            ) != bounded_path_ids
            or _lineage_component_path_ids(
                confirmation_paths, envelope_path_id
            ) != confirmation_path_ids):
        return {
            **base_evidence,
            "status": "rejected",
            "reason": "envelope_does_not_span_complete_lineage",
            "envelope_geometry": envelope,
        }

    closure = _contained_terminal_closure(
        full_dna,
        confirmation_paths,
        envelope,
    )
    return {
        **base_evidence,
        **closure,
        "envelope_geometry": envelope,
        "envelope_path_id": envelope_path_id,
        "bounded_lineage_path_count": len(bounded_path_ids),
        "confirmation_lineage_path_count": len(confirmation_path_ids),
        "shared_union_closure": "envelope_absorbing",
    }


def _path_local_convergence_proof(full_dna, state, confirmation_state):
    """Bind one canonical convergence result to complete replay witnesses."""

    superlocus_id = state.get("superlocus_id")
    confirmation_superlocus_id = confirmation_state.get("superlocus_id")
    if (
            not isinstance(superlocus_id, str)
            or not superlocus_id
            or confirmation_superlocus_id != superlocus_id):
        raise RuntimeError(
            "path-local proof states must share one source superlocus ID"
        )
    bounded_witness = _path_local_state_witness(
        full_dna, superlocus_id, state["paths"]
    )
    confirmation_witness = _path_local_state_witness(
        full_dna, superlocus_id, confirmation_state["paths"]
    )
    result = _path_local_convergence_result(
        full_dna, state, confirmation_state
    )
    return {
        "proof_schema": {
            "name": PATH_LOCAL_CONVERGENCE_PROOF_SCHEMA,
            "version": PATH_LOCAL_CONVERGENCE_PROOF_VERSION,
        },
        "source_sha256": _path_local_source_sha256(full_dna),
        "source_length": len(full_dna),
        "superlocus_id": superlocus_id,
        "bounded_paths": bounded_witness,
        "confirmation_paths": confirmation_witness,
        **result,
    }


def replay_source_superlocus_path_local_convergence_proof(
        full_dna, proof, *, expected_candidate=None,
        expected_path_id=None, expected_superlocus_id=None):
    """Rebuild and exactly replay one untrusted path-local proof.

    Repeat and spacer sequences are sliced anew from ``full_dna`` using every
    claimed half-open interval. Gapped alignments reconstruct those repeats and
    regenerate their consensuses. The bounded and confirmation states are then
    rerun through the same exhaustive path-local proof. The supplied record is
    accepted only when its type-strict canonical value exactly equals that
    independently regenerated result. Optional expected objects bind the proof
    envelope to the row and lineage identifiers being reported downstream.
    """

    if not isinstance(proof, dict):
        raise ValueError("path-local convergence proof must be a mapping")
    schema = proof.get("proof_schema")
    if (
            not isinstance(schema, dict)
            or set(schema) != {"name", "version"}
            or schema.get("name") != PATH_LOCAL_CONVERGENCE_PROOF_SCHEMA
            or schema.get("version") != PATH_LOCAL_CONVERGENCE_PROOF_VERSION):
        raise ValueError("path-local convergence proof schema is unsupported")
    if proof.get("source_sha256") != _path_local_source_sha256(full_dna):
        raise ValueError("path-local convergence proof source hash is stale")
    if (
            not isinstance(proof.get("source_length"), int)
            or isinstance(proof.get("source_length"), bool)
            or proof["source_length"] != len(full_dna)):
        raise ValueError("path-local convergence proof source length is stale")
    superlocus_id = proof.get("superlocus_id")
    if not isinstance(superlocus_id, str) or not superlocus_id:
        raise ValueError("path-local convergence proof superlocus is invalid")
    if "bounded_paths" not in proof or "confirmation_paths" not in proof:
        raise ValueError("path-local convergence proof witnesses are missing")

    bounded_state = _path_local_state_from_witness(
        full_dna,
        superlocus_id,
        proof["bounded_paths"],
        "bounded_paths",
    )
    confirmation_state = _path_local_state_from_witness(
        full_dna,
        superlocus_id,
        proof["confirmation_paths"],
        "confirmation_paths",
    )
    bounded_paths = bounded_state["paths"]
    confirmation_paths = confirmation_state["paths"]
    if not set(bounded_paths).issubset(confirmation_paths):
        raise ValueError(
            "path-local confirmation omits a bounded source geometry"
        )
    for geometry, bounded_entry in bounded_paths.items():
        confirmation_entry = confirmation_paths[geometry]
        if (
                bounded_entry["path_id"] != confirmation_entry["path_id"]
                or _candidate_terminal_behavior_key(
                    bounded_entry["candidate"]
                ) != _candidate_terminal_behavior_key(
                    confirmation_entry["candidate"]
                )):
            raise ValueError(
                "path-local confirmation changed a bounded candidate"
            )

    regenerated = _path_local_convergence_proof(
        full_dna, bounded_state, confirmation_state
    )
    if _path_local_proof_value(regenerated) != _path_local_proof_value(proof):
        raise ValueError(
            "path-local convergence proof is stale or has been tampered with"
        )
    if (
            expected_superlocus_id is not None
            and expected_superlocus_id != regenerated["superlocus_id"]):
        raise ValueError(
            "path-local convergence proof names a different superlocus"
        )
    if (
            expected_path_id is not None
            and expected_path_id != regenerated.get("envelope_path_id")):
        raise ValueError(
            "path-local convergence proof names a different envelope path"
        )
    if expected_candidate is not None:
        try:
            validate_candidate_source_reconstruction(
                full_dna,
                expected_candidate.list_repeats,
                expected_candidate.list_repeat_starts,
                expected_candidate.list_spacers,
            )
            expected_behavior = _candidate_terminal_behavior_key(
                expected_candidate
            )
        except (AttributeError, IndexError, TypeError, ValueError) as error:
            raise ValueError(
                "reported path-local candidate is malformed or stale"
            ) from error
        envelope_geometry = regenerated.get("envelope_geometry")
        envelope_entry = bounded_paths.get(envelope_geometry)
        if (
                envelope_entry is None
                or expected_behavior != _candidate_terminal_behavior_key(
                    envelope_entry["candidate"]
                )):
            raise ValueError(
                "path-local convergence envelope differs from reported candidate"
            )
    return regenerated


def _apply_augmentation_round(full_dna, states, round_number):
    """Apply one terminal-plus-shared round and return isolated round stats."""

    states_by_id = {
        state["superlocus_id"]: state for state in states
    }
    obstacle_mapping, eligible_runtime_ids = _obstacle_mapping(states)
    terminal_records, terminal_stats = (
        build_source_backed_terminal_candidate_variant_records(
            full_dna,
            obstacle_mapping,
            eligible_candidate_ids=eligible_runtime_ids,
            include_weak=True,
        )
    )
    terminal_added = 0
    terminal_added_by_superlocus = Counter()
    geometry_duplicates_skipped = 0
    for record in sorted(
        terminal_records,
        key=lambda item: (
            item["locus_key"],
            _candidate_sort_key(item["variant"]),
        ),
    ):
        state = states_by_id[record["locus_key"]]
        parent_path_id = state["candidate_path_ids"].get(
            id(record["base_candidate"])
        )
        if parent_path_id is None:
            raise RuntimeError(
                "terminal variant was proposed by a non-eligible obstacle"
            )
        _entry, added = _add_path(
            state,
            record["variant"],
            stage=(
                "weak_terminal_augmentation"
                if record["terminal_tier"] == "weak"
                else "terminal_augmentation"
            ),
            round_number=round_number,
            parent_path_ids=(parent_path_id,),
            details={
                "terminal_tier": record["terminal_tier"],
                "proposals": _normalize_terminal_proposals(
                    record["safe_proposals"]
                ),
            },
        )
        if added:
            terminal_added += 1
            terminal_added_by_superlocus[state["superlocus_id"]] += 1
        else:
            geometry_duplicates_skipped += 1

    shared_round_stats = Counter()
    shared_added = 0
    shared_added_by_superlocus = Counter()
    for state in states:
        current_entries = sorted(
            state["paths"].values(),
            key=lambda entry: entry["source_geometry"],
        )
        current_candidates = tuple(
            entry["candidate"] for entry in current_entries
        )
        for component in _exact_shared_components(current_candidates):
            if len(component) < 2:
                continue
            component_candidates = [
                current_candidates[index] for index in component
            ]
            augmented, component_stats = (
                augment_source_backed_shared_repeat_candidates(
                    full_dna,
                    {state["superlocus_id"]: component_candidates},
                )
            )
            for field in (
                "multi_candidate_component_count",
                "variants_built",
                "duplicate_variants_skipped",
                "inconsistent_overlap_rejections",
                "spacer_out_of_bounds_rejections",
                "unsupported_consecutive_edge_rejections",
            ):
                shared_round_stats[field] += component_stats[field]
            variants = augmented[state["superlocus_id"]][
                len(component_candidates):
            ]
            parent_path_ids = tuple(
                state["candidate_path_ids"][id(candidate)]
                for candidate in component_candidates
            )
            component_geometries = tuple(
                _source_geometry_key(candidate)
                for candidate in component_candidates
            )
            union_details = {
                "component_geometries": component_geometries,
            }
            duplicate_union_count = component_stats[
                "duplicate_variants_skipped"
            ]
            if duplicate_union_count:
                if duplicate_union_count != 1:
                    raise RuntimeError(
                        "one shared component produced multiple duplicate "
                        "union geometries"
                    )
                union_geometry = tuple(sorted({
                    interval
                    for geometry in component_geometries
                    for interval in geometry
                }))
                existing_union = state["paths"].get(union_geometry)
                if existing_union is None:
                    raise RuntimeError(
                        "shared augmentation reported an unknown duplicate "
                        "union geometry"
                    )
                try:
                    _entry, added = _add_path(
                        state,
                        existing_union["candidate"],
                        stage="exact_shared_repeat_union",
                        round_number=round_number,
                        parent_path_ids=parent_path_ids,
                        details=union_details,
                    )
                except _MultipleWeakLineageError:
                    shared_round_stats[
                        "multiple_weak_lineage_rejections"
                    ] += 1
                    continue
                if added:
                    raise RuntimeError(
                        "duplicate shared union unexpectedly added a path"
                    )
            for variant in variants:
                try:
                    _entry, added = _add_path(
                        state,
                        variant,
                        stage="exact_shared_repeat_union",
                        round_number=round_number,
                        parent_path_ids=parent_path_ids,
                        details=union_details,
                    )
                except _MultipleWeakLineageError:
                    shared_round_stats[
                        "multiple_weak_lineage_rejections"
                    ] += 1
                    continue
                if added:
                    shared_added += 1
                    shared_added_by_superlocus[state["superlocus_id"]] += 1
                else:
                    geometry_duplicates_skipped += 1

    shared_stat_fields = (
        "multi_candidate_component_count",
        "variants_built",
        "duplicate_variants_skipped",
        "inconsistent_overlap_rejections",
        "spacer_out_of_bounds_rejections",
        "unsupported_consecutive_edge_rejections",
        "multiple_weak_lineage_rejections",
    )
    return {
        "terminal": dict(terminal_stats),
        "shared": {
            field: shared_round_stats[field]
            for field in shared_stat_fields
        },
        "terminal_geometries_added": terminal_added,
        "shared_geometries_added": shared_added,
        "terminal_added_by_superlocus": dict(
            terminal_added_by_superlocus
        ),
        "shared_added_by_superlocus": dict(shared_added_by_superlocus),
        "geometry_duplicates_skipped": geometry_duplicates_skipped,
    }


def _validate_inputs(full_dna, candidates_by_legacy_locus, eligible_candidate_ids):
    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    if not isinstance(candidates_by_legacy_locus, dict):
        raise ValueError("legacy candidates must be supplied as a mapping")
    try:
        eligible_ids = frozenset(eligible_candidate_ids)
    except TypeError as error:
        raise ValueError(
            "eligible candidate IDs must be an iterable of object IDs"
        ) from error
    if any(
        not isinstance(candidate_id, int)
        or isinstance(candidate_id, bool)
        or candidate_id < 0
        for candidate_id in eligible_ids
    ):
        raise ValueError("eligible candidate IDs must be non-negative integers")

    records = []
    observed_ids = set()
    for locus_key, candidates in candidates_by_legacy_locus.items():
        if not isinstance(candidates, (tuple, list)) or not candidates:
            raise ValueError(
                "each legacy locus must contain a non-empty candidate list"
            )
        validated = []
        source_path_families = {}
        for candidate in candidates:
            validate_candidate_source_reconstruction(
                full_dna,
                candidate.list_repeats,
                candidate.list_repeat_starts,
                candidate.list_spacers,
            )
            _normalized_consensus(candidate)
            source_path_families[id(candidate)] = (
                _validated_source_path_family_keys(full_dna, candidate)
            )
            validated.append(candidate)
            observed_ids.add(id(candidate))
        validated = _unique_sorted_candidates(validated)
        intervals = tuple(_candidate_interval(candidate) for candidate in validated)
        records.append({
            "locus_key": locus_key,
            "candidates": validated,
            "eligible_candidates": tuple(
                candidate for candidate in validated
                if id(candidate) in eligible_ids
            ),
            "source_path_families": {
                id(candidate): source_path_families[id(candidate)]
                for candidate in validated
            },
            "source_interval": (
                min(interval[0] for interval in intervals),
                max(interval[1] for interval in intervals),
            ),
        })
    records.sort(key=lambda record: (
        record["source_interval"],
        _stable_value_key(record["locus_key"]),
    ))
    return tuple(records), eligible_ids, observed_ids


def augment_source_superloci_to_fixed_point(
    full_dna,
    candidates_by_legacy_locus,
    eligible_candidate_ids,
    max_rounds=MAXIMUM_ENHANCEMENT_ROUNDS,
):
    """Build source superloci and iteratively add eligible path variants.

    Returns ``(candidates_by_superlocus, lineage, stats)``.  Output mappings
    contain eligible legacy candidates and their source-derived variants only;
    ineligible originals remain internal collision obstacles.  A round always
    runs one-step terminal augmentation first and exact-shared-repeat graph
    augmentation second.  Iteration stops at a geometry fixed point or after
    the requested cap, which must be between one and four inclusive.

    Reaching the cap still marks a changing superlocus nonconverged. A unique
    bounded envelope can be marked commit-eligible only by the independent
    path-local proof described in the module documentation; every other
    generated candidate stays quarantined and all downstream structural,
    tandem, coding, and classifier gates remain mandatory.
    """

    if (
        not isinstance(max_rounds, int)
        or isinstance(max_rounds, bool)
        or not 1 <= max_rounds <= MAXIMUM_ENHANCEMENT_ROUNDS
    ):
        raise ValueError("max_rounds must be an integer between one and four")

    legacy_records, eligible_ids, observed_ids = _validate_inputs(
        full_dna,
        candidates_by_legacy_locus,
        eligible_candidate_ids,
    )
    adjacency = [set() for _record in legacy_records]
    link_records = []
    for first_index, first in enumerate(legacy_records[:-1]):
        for second_index in range(first_index + 1, len(legacy_records)):
            second = legacy_records[second_index]
            evidence = _legacy_locus_link_evidence(first, second)
            if evidence is None:
                continue
            adjacency[first_index].add(second_index)
            adjacency[second_index].add(first_index)
            link_records.append({
                "first_legacy_locus": first["locus_key"],
                "second_legacy_locus": second["locus_key"],
                "evidence": evidence,
            })

    components = _connected_components(adjacency)
    states = []
    for component_index, component in enumerate(components, 1):
        component_records = tuple(
            legacy_records[index] for index in component
        )
        start = min(
            record["source_interval"][0] for record in component_records
        )
        end = max(
            record["source_interval"][1] for record in component_records
        )
        superlocus_id = "superlocus-{:04d}-{}-{}".format(
            component_index,
            start,
            end,
        )
        state = {
            "full_dna": full_dna,
            "superlocus_id": superlocus_id,
            "source_interval": (start, end),
            "legacy_locus_keys": tuple(
                record["locus_key"] for record in component_records
            ),
            "legacy_obstacles": _unique_sorted_candidates(
                candidate
                for record in component_records
                for candidate in record["candidates"]
            ),
            "paths": {},
            "candidate_path_ids": {},
        }
        eligible_originals = _unique_sorted_candidates(
            candidate
            for record in component_records
            for candidate in record["eligible_candidates"]
        )
        geometry_locus_keys = defaultdict(list)
        for record in component_records:
            for candidate in record["eligible_candidates"]:
                geometry_locus_keys[_source_geometry_key(candidate)].append(
                    record["locus_key"]
                )
        for candidate in eligible_originals:
            geometry = _source_geometry_key(candidate)
            entry, added = _add_path(
                state,
                candidate,
                stage="legacy_eligible",
                round_number=0,
                details={
                    "legacy_locus_keys": tuple(sorted(
                        geometry_locus_keys[geometry],
                        key=_stable_value_key,
                    )),
                    "runtime_candidate_id": id(candidate),
                },
            )
            state["candidate_path_ids"][id(candidate)] = entry["path_id"]
        states.append(state)

    stats = {
        "legacy_locus_count": len(legacy_records),
        "source_validated_candidate_count": sum(
            len(record["candidates"]) for record in legacy_records
        ),
        "eligible_legacy_candidate_count": sum(
            len(record["eligible_candidates"]) for record in legacy_records
        ),
        "ineligible_legacy_candidate_count": sum(
            len(record["candidates"]) - len(record["eligible_candidates"])
            for record in legacy_records
        ),
        "unmatched_eligible_candidate_ids": len(
            eligible_ids.difference(observed_ids)
        ),
        "superlocus_count": len(states),
        "superlocus_link_count": len(link_records),
        "exact_shared_repeat_link_count": sum(
            record["evidence"]["mechanism"]
            == "exact_shared_repeat_interval"
            for record in link_records
        ),
        "compatible_topology_link_count": sum(
            record["evidence"]["mechanism"]
            == "compatible_source_topology"
            for record in link_records
        ),
        "source_path_family_link_count": sum(
            record["evidence"]["mechanism"]
            == "validated_source_path_family"
            for record in link_records
        ),
        "source_path_provenance_candidate_count": sum(
            bool(record["source_path_families"][id(candidate)])
            for record in legacy_records
            for candidate in record["candidates"]
        ),
        "max_rounds": max_rounds,
        "rounds_run": 0,
        "fixed_point_reached": False,
        "max_rounds_reached": False,
        "confirmation_run": False,
        "confirmation_round": None,
        "confirmed_converged_superlocus_count": 0,
        "nonconverged_superlocus_count": 0,
        "commit_eligible_superlocus_count": 0,
        "path_local_confirmation_attempted_superlocus_count": 0,
        "path_local_converged_superlocus_count": 0,
        "path_local_rejected_superlocus_count": 0,
        "path_local_commit_eligible_candidate_count": 0,
        "path_local_terminal_closure_explored_geometry_count": 0,
        "path_local_terminal_closure_explored_behavior_state_count": 0,
        "path_local_terminal_closure_generated_geometry_count": 0,
        "path_local_terminal_closure_generated_behavior_state_count": 0,
        "path_local_terminal_closure_weak_closed_geometry_count": 0,
        "path_local_terminal_closure_weak_closed_behavior_state_count": 0,
        "path_local_terminal_closure_transition_count": 0,
        "path_local_terminal_closure_state_budget_rejection_count": 0,
        "quarantined_candidate_count": 0,
        "terminal_proposed_repeats": 0,
        "terminal_collision_skipped_repeats": 0,
        "terminal_candidate_peer_collision_skipped_repeats": 0,
        "terminal_proposal_collision_skipped_repeats": 0,
        "terminal_variants_built": 0,
        "terminal_variants_added": 0,
        "terminal_duplicate_variants_skipped": 0,
        "terminal_weak_proposed_repeats": 0,
        "terminal_weak_variants_built": 0,
        "terminal_weak_variants_added": 0,
        "terminal_weak_duplicate_variants_skipped": 0,
        "shared_multi_candidate_component_count": 0,
        "shared_variants_built": 0,
        "shared_variants_added": 0,
        "shared_duplicate_variants_skipped": 0,
        "shared_inconsistent_overlap_rejections": 0,
        "shared_spacer_out_of_bounds_rejections": 0,
        "shared_unsupported_consecutive_edge_rejections": 0,
        "shared_multiple_weak_lineage_rejections": 0,
        "geometry_duplicates_skipped": 0,
        "output_candidate_count": 0,
        "link_evidence": tuple(link_records),
        "round_stats": (),
        "superloci": (),
    }
    round_stats = []
    terminal_fields = (
        "proposed_repeats",
        "collision_skipped_repeats",
        "candidate_peer_collision_skipped_repeats",
        "proposal_collision_skipped_repeats",
        "variants_built",
        "duplicate_variants_skipped",
        "weak_proposed_repeats",
        "weak_variants_built",
        "weak_variants_added",
        "weak_duplicate_variants_skipped",
    )
    shared_fields = (
        "multi_candidate_component_count",
        "variants_built",
        "duplicate_variants_skipped",
        "inconsistent_overlap_rejections",
        "spacer_out_of_bounds_rejections",
        "unsupported_consecutive_edge_rejections",
    )

    for round_number in range(1, max_rounds + 1):
        before = {
            state["superlocus_id"]: _state_semantic_signature(state)
            for state in states
        }
        round_result = _apply_augmentation_round(
            full_dna, states, round_number
        )
        changed_superloci = tuple(
            state["superlocus_id"]
            for state in states
            if before[state["superlocus_id"]]
            != _state_semantic_signature(state)
        )

        terminal_added = round_result["terminal_geometries_added"]
        shared_added = round_result["shared_geometries_added"]
        stats["terminal_variants_added"] += terminal_added
        _sum_stats(
            stats, round_result["terminal"], "terminal_", terminal_fields
        )
        stats["shared_variants_added"] += shared_added
        _sum_stats(
            stats, round_result["shared"], "shared_", shared_fields
        )
        stats["shared_multiple_weak_lineage_rejections"] += (
            round_result["shared"].get(
                "multiple_weak_lineage_rejections", 0
            )
        )
        stats["geometry_duplicates_skipped"] += round_result[
            "geometry_duplicates_skipped"
        ]

        new_geometry_count = terminal_added + shared_added
        stats["rounds_run"] = round_number
        round_stats.append({
            "round": round_number,
            "terminal": round_result["terminal"],
            "shared": round_result["shared"],
            "terminal_geometries_added": terminal_added,
            "shared_geometries_added": shared_added,
            "new_geometry_count": new_geometry_count,
            "semantic_state_changed": bool(changed_superloci),
            "semantic_changed_superloci": changed_superloci,
        })
        if not changed_superloci:
            stats["fixed_point_reached"] = True
            break

    fixed_point_observed_during_mutation = stats["fixed_point_reached"]
    stats["max_rounds_reached"] = not fixed_point_observed_during_mutation
    stats["round_stats"] = tuple(round_stats)

    convergence = {}
    confirmation_states_by_id = {}
    if fixed_point_observed_during_mutation:
        for state in states:
            convergence[state["superlocus_id"]] = {
                "status": "converged",
                "reason": "semantic_fixed_point_observed",
                "confirmation_semantic_state_changed": False,
                "confirmation_terminal_geometries_added": 0,
                "confirmation_shared_geometries_added": 0,
            }
    else:
        stats["confirmation_run"] = True
        stats["confirmation_round"] = max_rounds + 1
        original_signatures = {
            state["superlocus_id"]: _state_semantic_signature(state)
            for state in states
        }
        confirmation_states = tuple(
            _clone_state_for_confirmation(state) for state in states
        )
        confirmation_states_by_id = {
            state["superlocus_id"]: state
            for state in confirmation_states
        }
        confirmation_result = _apply_augmentation_round(
            full_dna, confirmation_states, max_rounds + 1
        )
        for state in states:
            if (
                    _state_semantic_signature(state)
                    != original_signatures[state["superlocus_id"]]):
                raise RuntimeError(
                    "bounded confirmation mutated committed lineage state"
                )
        for confirmation_state in confirmation_states:
            superlocus_id = confirmation_state["superlocus_id"]
            semantic_changed = (
                _state_semantic_signature(confirmation_state)
                != original_signatures[superlocus_id]
            )
            terminal_added = confirmation_result[
                "terminal_added_by_superlocus"
            ].get(superlocus_id, 0)
            shared_added = confirmation_result[
                "shared_added_by_superlocus"
            ].get(superlocus_id, 0)
            if semantic_changed:
                mechanisms = []
                if terminal_added:
                    mechanisms.append("terminal_geometry")
                if shared_added:
                    mechanisms.append("shared_geometry")
                if not mechanisms:
                    mechanisms.append("lineage_semantics")
                status = "nonconverged"
                reason = "confirmation_changed_{}".format(
                    "_and_".join(mechanisms)
                )
            else:
                status = "converged"
                reason = "bounded_confirmation_fixed_point"
            convergence[superlocus_id] = {
                "status": status,
                "reason": reason,
                "confirmation_semantic_state_changed": semantic_changed,
                "confirmation_terminal_geometries_added": terminal_added,
                "confirmation_shared_geometries_added": shared_added,
            }

    stats["confirmed_converged_superlocus_count"] = sum(
        record["status"] == "converged"
        for record in convergence.values()
    )
    stats["nonconverged_superlocus_count"] = sum(
        record["status"] == "nonconverged"
        for record in convergence.values()
    )
    stats["commit_eligible_superlocus_count"] = stats[
        "confirmed_converged_superlocus_count"
    ]
    stats["fixed_point_reached"] = (
        stats["nonconverged_superlocus_count"] == 0
    )

    path_local_proofs = {}
    for state in states:
        superlocus_id = state["superlocus_id"]
        if convergence[superlocus_id]["status"] == "converged":
            continue
        proof = _path_local_convergence_proof(
            full_dna,
            state,
            confirmation_states_by_id[superlocus_id],
        )
        path_local_proofs[superlocus_id] = proof
    stats["path_local_confirmation_attempted_superlocus_count"] = len(
        path_local_proofs
    )
    stats["path_local_converged_superlocus_count"] = sum(
        proof["status"] == "converged"
        for proof in path_local_proofs.values()
    )
    stats["path_local_rejected_superlocus_count"] = sum(
        proof["status"] == "rejected"
        for proof in path_local_proofs.values()
    )
    stats["path_local_commit_eligible_candidate_count"] = stats[
        "path_local_converged_superlocus_count"
    ]
    for proof in path_local_proofs.values():
        stats[
            "path_local_terminal_closure_explored_geometry_count"
        ] += proof.get(
            "terminal_closure_explored_geometry_count", 0
        )
        stats[
            "path_local_terminal_closure_explored_behavior_state_count"
        ] += proof.get(
            "terminal_closure_explored_behavior_state_count", 0
        )
        stats[
            "path_local_terminal_closure_generated_geometry_count"
        ] += proof.get(
            "terminal_closure_generated_geometry_count", 0
        )
        stats[
            "path_local_terminal_closure_generated_behavior_state_count"
        ] += proof.get(
            "terminal_closure_generated_behavior_state_count", 0
        )
        stats[
            "path_local_terminal_closure_weak_closed_geometry_count"
        ] += proof.get(
            "terminal_closure_weak_closed_geometry_count", 0
        )
        stats[
            "path_local_terminal_closure_weak_closed_behavior_state_count"
        ] += proof.get(
            "terminal_closure_weak_closed_behavior_state_count", 0
        )
        stats[
            "path_local_terminal_closure_transition_count"
        ] += proof.get("terminal_closure_transition_count", 0)
        stats[
            "path_local_terminal_closure_state_budget_rejection_count"
        ] += proof.get("reason") == "terminal_closure_state_budget_exceeded"

    outputs = {}
    lineage = []
    superlocus_stats = []
    ineligible_runtime_ids = observed_ids.difference(eligible_ids)
    for state in states:
        entries = sorted(
            state["paths"].values(),
            key=lambda entry: entry["source_geometry"],
        )
        convergence_record = convergence[state["superlocus_id"]]
        converged = convergence_record["status"] == "converged"
        path_local_proof = path_local_proofs.get(state["superlocus_id"])
        path_local_envelope = (
            path_local_proof["envelope_geometry"]
            if path_local_proof is not None
            and path_local_proof["status"] == "converged"
            else None
        )
        for entry in entries:
            legacy_path = entry["stage"] == "legacy_eligible"
            path_local_converged = (
                entry["source_geometry"] == path_local_envelope
            )
            entry["superlocus_convergence_status"] = convergence_record[
                "status"
            ]
            entry["superlocus_convergence_reason"] = convergence_record[
                "reason"
            ]
            if converged:
                path_local_status = "not_required"
                path_local_reason = "superlocus_converged"
            elif path_local_converged:
                path_local_status = "converged"
                path_local_reason = path_local_proof["reason"]
            elif (
                    path_local_proof is not None
                    and path_local_proof["status"] == "converged"):
                path_local_status = "unresolved"
                path_local_reason = "dominated_contained_subpath"
            else:
                path_local_status = "rejected"
                path_local_reason = path_local_proof["reason"]
            entry["path_local_convergence_status"] = path_local_status
            entry["path_local_convergence_reason"] = path_local_reason
            entry["path_local_convergence_proof"] = (
                dict(path_local_proof)
                if path_local_converged else None
            )
            entry["augmentation_output_eligible"] = (
                converged or legacy_path or path_local_converged
            )
            entry["architecture_commit_eligible"] = (
                converged or path_local_converged
            )
        bounded_candidates = tuple(
            entry["candidate"] for entry in entries
        )
        candidates = tuple(
            entry["candidate"] for entry in entries
            if entry["augmentation_output_eligible"]
        )
        if any(
                id(candidate) in ineligible_runtime_ids
                for candidate in bounded_candidates):
            raise RuntimeError("an ineligible legacy candidate reached output")
        for candidate in bounded_candidates:
            validate_candidate_source_reconstruction(
                full_dna,
                candidate.list_repeats,
                candidate.list_repeat_starts,
                candidate.list_spacers,
            )
        outputs[state["superlocus_id"]] = candidates
        lineage.extend(entries)
        superlocus_stats.append({
            "superlocus_id": state["superlocus_id"],
            "source_interval": state["source_interval"],
            "legacy_locus_keys": state["legacy_locus_keys"],
            "collision_obstacle_count": len(state["legacy_obstacles"]),
            "bounded_candidate_count": len(bounded_candidates),
            "output_candidate_count": len(candidates),
            "quarantined_generated_candidate_count": (
                len(bounded_candidates) - len(candidates)
            ),
            "convergence_status": convergence_record["status"],
            "convergence_reason": convergence_record["reason"],
            "confirmation_semantic_state_changed": convergence_record[
                "confirmation_semantic_state_changed"
            ],
            "confirmation_terminal_geometries_added": convergence_record[
                "confirmation_terminal_geometries_added"
            ],
            "confirmation_shared_geometries_added": convergence_record[
                "confirmation_shared_geometries_added"
            ],
            "path_local_convergence_status": (
                "not_required" if converged
                else path_local_proof["status"]
            ),
            "path_local_convergence_reason": (
                "superlocus_converged" if converged
                else path_local_proof["reason"]
            ),
            "path_local_commit_eligible_candidate_count": int(
                path_local_envelope is not None
            ),
            "path_local_convergence_proof": (
                None if converged else dict(path_local_proof)
            ),
            "architecture_commit_eligible": converged,
        })

    lineage.sort(key=lambda entry: (
        entry["superlocus_id"],
        entry["source_geometry"],
    ))
    stats["output_candidate_count"] = sum(map(len, outputs.values()))
    stats["lineage_candidate_count"] = len(lineage)
    stats["quarantined_candidate_count"] = sum(
        record["quarantined_generated_candidate_count"]
        for record in superlocus_stats
    )
    stats["superloci"] = tuple(superlocus_stats)
    return outputs, tuple(lineage), stats


__all__ = [
    "MAXIMUM_ENHANCEMENT_ROUNDS",
    "MAXIMUM_SUPERLOCUS_GAP",
    "MINIMUM_CONSENSUS_CONTAINMENT_IDENTITY",
    "MINIMUM_CONSENSUS_LENGTH_RATIO",
    "PATH_LOCAL_CONVERGENCE_PROOF_SCHEMA",
    "PATH_LOCAL_CONVERGENCE_PROOF_VERSION",
    "augment_source_superloci_to_fixed_point",
    "replay_source_superlocus_path_local_convergence_proof",
]
