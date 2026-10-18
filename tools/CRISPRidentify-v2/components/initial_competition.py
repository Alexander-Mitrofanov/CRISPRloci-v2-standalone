"""Preserve ordinary competition; fence terminal-repair-origin scope bridges.

Repair-only connected groups touching several frozen ordinary components stay
noncanonical. Ties disclose a legacy representative choice, not unique phase.
"""
from dataclasses import dataclass
import hashlib
import json
import math

from components.candidate_origins import candidate_origin_kind
from components.fast_run_rescue import candidate_source_interval
from components.occurrence_protection import (
    occurrence_relation, occurrence_snapshot, validate_snapshot,
)

RULE = "repair_origin_scope_fence_v1"
OCCURRENCE_RULE = "repair_origin_occurrence_fence_v2"
METADATA_KEY = "initial_competition"
SCORE_TIE_TOLERANCE = 1e-12
REPRESENTATIVE_RULE = "maximum_score_then_legacy_stable_input_order"


def _geometry(candidate):
    candidate_source_interval(candidate)
    return tuple((start, start + len(repeat)) for start, repeat in zip(
        candidate.list_repeat_starts, candidate.list_repeats))


def _overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _score(row):
    value = float(row[0])
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("initial competition requires finite probability scores")
    return value


def _components(records):
    groups, maximum_end = [], -1
    for record in sorted(records, key=lambda item: (item[2], item[0])):
        if record[2][0] >= maximum_end:
            groups.append([])
            maximum_end = -1
        groups[-1].append(record)
        maximum_end = max(maximum_end, record[2][1])
    return [sorted(group, key=lambda item: item[0]) for group in groups]


def _descriptor(records):
    geometries = sorted({_geometry(row[1]) for _, row, _ in records})
    identity = hashlib.sha256(json.dumps(geometries, separators=(",", ":")).encode()).hexdigest()
    winner = max((row for _, row, _ in records), key=_score)
    return {"component_id": identity,
            "source_span": [min(x[2][0] for x in records), max(x[2][1] for x in records)],
            "geometry_count": len(geometries),
            "winner_geometry": [list(pair) for pair in _geometry(winner[1])],
            "winner_score": _score(winner)}


@dataclass(frozen=True)
class CompetitionGroup:
    rows: tuple
    evidence: dict


def _evidence(row, records, owners, repairs, conflict=False):
    tied = {}
    if not conflict:
        for _, other, _ in records:
            if 0 <= _score(row) - _score(other) <= SCORE_TIE_TOLERANCE:
                geometry = _geometry(other[1])
                tied[geometry] = max(_score(other), tied.get(geometry, 0.0))
    return {
        "rule": RULE, "decision": "repair_scope_conflict" if conflict else "ranked_scope",
        "coordinate_system": "zero_based_half_open",
        "candidate_geometry": [list(pair) for pair in _geometry(row[1])],
        "candidate_score": _score(row),
        "scope_span": [min(x[2][0] for x in records), max(x[2][1] for x in records)],
        "scope_kind": ("conflicting_repair_group" if conflict else
                       "ordinary_with_repairs" if owners and repairs else
                       "ordinary" if owners else "repair_only"),
        "ordinary_components": owners, "repair_components": repairs,
        "representative_rule": REPRESENTATIVE_RULE,
        "tie_tolerance": SCORE_TIE_TOLERANCE,
        "tied_candidates": [{"repeat_intervals": [list(pair) for pair in geometry], "score": tied[geometry]}
                            for geometry in sorted(tied)],
        "unique_phase": False if len(tied) > 1 else None,
    }


def _occurrence_evidence(row, records, owners, repairs, reference, observation,
                         conflicts, decision="ranked_scope", scope_span=None):
    # Noncanonical observations do not need a repeated traversal of the full
    # component to recompute ties or span. The complete scope was frozen once.
    local_records = (records if decision == "ranked_scope" else
                     [(0, row, candidate_source_interval(row[1]))])
    evidence = _evidence(row, local_records, owners, repairs)
    if scope_span is not None:
        evidence["scope_span"] = scope_span
    evidence.update(rule=OCCURRENCE_RULE, decision=decision)
    if decision != "ranked_scope":
        evidence.update(tied_candidates=[], unique_phase=None)
    # The reference and relation are shared immutable-value snapshots within
    # this categorization. Only the canonical row owns the conflict inventory.
    evidence["occurrence_protection"] = {
        "ordinary_reference": reference, "candidate": observation["candidate"],
        "relation": observation["relation"], "deferred_candidates": list(conflicts),
    }
    return evidence


def initial_competition_groups(rows, split_ordinary_components=True):
    """Return ranked groups and deferred rows, retaining each input once.

    Owners are computed against frozen ordinary components, never previously
    admitted repairs. A singly touching repair can be deferred if its repair
    component reaches another ordinary component through other repairs.
    """
    ordinary, repair = [], []
    for index, row in enumerate(rows):
        _score(row)
        record = (index, row, candidate_source_interval(row[1]))
        (ordinary if candidate_origin_kind(row[1]) == "ordinary" else repair).append(record)
    base = (_components(ordinary) if split_ordinary_components else
            [ordinary] if ordinary else [])
    descriptors = [_descriptor(group) for group in base]
    combined, attached = [list(group) for group in base], [[] for _ in base]
    orphans, conflicts = [], []
    for records in _components(repair):
        descriptor = _descriptor(records)
        owners = [i for i, original in enumerate(base)
                  if any(_overlap(record[2], other[2]) for record in records for other in original)]
        descriptor["direct_ordinary_component_ids"] = [descriptors[i]["component_id"] for i in owners]
        if len(owners) == 1:
            combined[owners[0]].extend(records)
            attached[owners[0]].append(descriptor)
        elif not owners:
            orphans.append((records, descriptor))
        else:
            for _, row, _ in records:
                conflicts.append(with_competition_evidence(row, _evidence(
                    row, records, [descriptors[i] for i in owners], [descriptor], True)))
    groups = []
    for i, records in enumerate(combined):
        records.sort(key=lambda item: item[0])
        protected = max((row for _, row, _ in base[i]), key=_score)
        reference = occurrence_snapshot(protected) if attached[i] else None
        admitted, deferred, memo = [], [], {}
        scope_span = [min(x[2][0] for x in records), max(x[2][1] for x in records)]
        ordinary_ids = {id(row) for _, row, _ in base[i]}
        for index, row, span in records:
            if reference is None or id(row) in ordinary_ids:
                admitted.append((index, row, span))
                continue
            geometry = _geometry(row[1])
            if geometry not in memo:
                memo[geometry] = occurrence_relation(reference["repeat_intervals"], geometry)
            relation = memo[geometry]
            observation = {"candidate": occurrence_snapshot(row), "relation": relation}
            evidence = _occurrence_evidence(row, records, [descriptors[i]], attached[i],
                reference, observation, (), decision=("repair_occurrence_preserved" if
                relation["status"] == "unique" else "repair_occurrence_conflict"), scope_span=scope_span)
            annotated = with_competition_evidence(row, evidence)
            if relation["status"] == "unique":
                admitted.append((index, annotated, span))
            else:
                deferred.append(observation)
                conflicts.append(annotated)
        winner = max((row for _, row, _ in admitted), key=_score)
        if reference is None:
            evidence = _evidence(winner, admitted, [descriptors[i]], attached[i])
        else:
            observation = {"candidate": occurrence_snapshot(winner),
                "relation": occurrence_relation(reference["repeat_intervals"], _geometry(winner[1]))}
            evidence = _occurrence_evidence(winner, admitted, [descriptors[i]], attached[i],
                reference, observation, deferred, scope_span=scope_span)
        groups.append(CompetitionGroup(tuple(row for _, row, _ in admitted), evidence))
    for records, descriptor in orphans:
        winner = max((row for _, row, _ in records), key=_score)
        groups.append(CompetitionGroup(tuple(row for _, row, _ in records),
            _evidence(winner, records, [], [descriptor])))
    return tuple(groups), tuple(conflicts)


def clear_competition_evidence(row):
    if len(row) < 4 or METADATA_KEY not in row[3]:
        return row
    metadata = dict(row[3])
    del metadata[METADATA_KEY]
    return list(row[:3]) + ([metadata] if metadata else [])


def with_competition_evidence(row, evidence):
    metadata = dict(row[3]) if len(row) > 3 else {}
    metadata[METADATA_KEY] = evidence
    return list(row[:3]) + [metadata]


def _checked_geometry(intervals, source_length):
    if not isinstance(intervals, list) or not intervals:
        raise ValueError("initial competition geometry is malformed")
    geometry = []
    for pair in intervals:
        if (not isinstance(pair, list) or len(pair) != 2
                or any(type(value) is not int for value in pair)
                or not 0 <= pair[0] < pair[1]
                or (geometry and pair[0] < geometry[-1][1])
                or (source_length is not None and pair[1] > source_length)):
            raise ValueError("initial competition intervals are invalid")
        geometry.append(tuple(pair))
    return tuple(geometry)


def competition_evidence(row, source_length=None, *, source=None):
    """Validate a bound observation, not the full historical scope inventory.

    Ties can be disjoint within a legacy transitive component: no assertion of
    direct winner overlap or independently resolved biological phase is made.
    """
    if len(row) < 4 or METADATA_KEY not in row[3]:
        return None
    evidence = row[3][METADATA_KEY]
    required = {"rule", "decision", "coordinate_system", "candidate_geometry", "candidate_score",
                "scope_span", "scope_kind", "ordinary_components", "repair_components",
                "representative_rule", "tie_tolerance", "tied_candidates", "unique_phase"}
    v2 = isinstance(evidence, dict) and evidence.get("rule") == OCCURRENCE_RULE
    if v2:
        required.add("occurrence_protection")
    if not isinstance(evidence, dict) or set(evidence) != required:
        raise ValueError("initial competition evidence is malformed")
    expected = {"rule": OCCURRENCE_RULE if v2 else RULE, "coordinate_system": "zero_based_half_open",
                "tie_tolerance": SCORE_TIE_TOLERANCE, "representative_rule": REPRESENTATIVE_RULE,
                "candidate_geometry": [list(pair) for pair in _geometry(row[1])],
                "candidate_score": _score(row)}
    if any(evidence.get(key) != value for key, value in expected.items()):
        raise ValueError("initial competition evidence differs from its candidate")
    span = _checked_geometry([evidence["scope_span"]], source_length)[0]
    candidate = _geometry(row[1])
    if not span[0] <= candidate[0][0] < candidate[-1][1] <= span[1]:
        raise ValueError("initial competition candidate is outside its scope")
    owners, repairs = evidence["ordinary_components"], evidence["repair_components"]
    if not isinstance(owners, list) or not isinstance(repairs, list):
        raise ValueError("initial competition components are malformed")
    for component in owners + repairs:
        own_span = _checked_geometry([component["source_span"]], source_length)[0]
        winner = _checked_geometry(component["winner_geometry"], source_length)
        if (not own_span[0] <= winner[0][0] < winner[-1][1] <= own_span[1]
                or type(component["geometry_count"]) is not int or component["geometry_count"] < 1
                or not isinstance(component["component_id"], str) or len(component["component_id"]) != 64):
            raise ValueError("initial competition component descriptor is invalid")
        _score([component["winner_score"]])
    owner_ids = [x["component_id"] for x in owners]
    if len(owner_ids) != len(set(owner_ids)):
        raise ValueError("initial competition repeats an ordinary component")
    for repair in repairs:
        if repair["direct_ordinary_component_ids"] != owner_ids:
            raise ValueError("initial competition repair owner union differs")
    if evidence["decision"] == "repair_scope_conflict":
        if (len(owners) < 2 or len(repairs) != 1 or evidence["tied_candidates"]
                or evidence["scope_kind"] != "conflicting_repair_group"
                or evidence["unique_phase"] is not None
                or any(not _overlap(span, owner["source_span"]) for owner in owners)):
            raise ValueError("initial competition conflict must touch multiple ordinary components")
    elif v2 and evidence["decision"] in ("repair_occurrence_conflict", "repair_occurrence_preserved"):
        if (len(owners) != 1 or not repairs or evidence["scope_kind"] != "ordinary_with_repairs"
                or evidence["tied_candidates"] or evidence["unique_phase"] is not None):
            raise ValueError("initial competition occurrence decision has invalid scope")
    elif evidence["decision"] == "ranked_scope":
        kind = "ordinary_with_repairs" if owners and repairs else "ordinary" if owners else "repair_only"
        if len(owners) > 1 or (not owners and len(repairs) != 1) or evidence["scope_kind"] != kind:
            raise ValueError("initial competition ranked scope has conflicting owners")
        tied = evidence["tied_candidates"]
        if not isinstance(tied, list) or not tied:
            raise ValueError("initial competition ranking lacks a representative")
        geometries, winner_seen = [], False
        for entry in tied:
            geometry = _checked_geometry(entry["repeat_intervals"], source_length)
            score = _score([entry["score"]])
            if (not span[0] <= geometry[0][0] < geometry[-1][1] <= span[1]
                    or not 0 <= _score(row) - score <= SCORE_TIE_TOLERANCE):
                raise ValueError("initial competition tie differs from its scope or score")
            geometries.append(geometry)
            winner_seen |= geometry == candidate and score == _score(row)
        if (not winner_seen or geometries != sorted(set(geometries))
                or evidence["unique_phase"] is not (False if len(geometries) > 1 else None)):
            raise ValueError("initial competition tied geometry observations are inconsistent")
    else:
        raise ValueError("initial competition decision is unknown")
    if v2:
        _validate_occurrence_evidence(row, evidence, source, source_length)
    return evidence


def _validate_occurrence_evidence(row, evidence, source, source_length):
    protection = evidence["occurrence_protection"]
    if (not isinstance(protection, dict) or set(protection) != {
            "ordinary_reference", "candidate", "relation", "deferred_candidates"}
            or len(evidence["ordinary_components"]) != 1 or not evidence["repair_components"]):
        raise ValueError("initial competition occurrence protection is malformed")
    reference = protection["ordinary_reference"]
    geometry = validate_snapshot(reference, source, source_length)
    owner = evidence["ordinary_components"][0]
    if (reference["repeat_intervals"] != owner["winner_geometry"]
            or reference["score"] != owner["winner_score"]):
        raise ValueError("initial competition ordinary witness differs from its frozen owner")
    current = protection["candidate"]
    validate_snapshot(current, source, source_length)
    if current != occurrence_snapshot(row):
        raise ValueError("initial competition occurrence candidate binding differs")
    relation = occurrence_relation(geometry, current["repeat_intervals"])
    if protection["relation"] != relation:
        raise ValueError("initial competition occurrence mapping differs")
    conflict = evidence["decision"] == "repair_occurrence_conflict"
    if conflict == (relation["status"] == "unique"):
        raise ValueError("initial competition occurrence disposition differs")
    deferred = protection["deferred_candidates"]
    if not isinstance(deferred, list) or (evidence["decision"] != "ranked_scope" and deferred):
        raise ValueError("initial competition deferred inventory is malformed")
    for item in deferred:
        if not isinstance(item, dict) or set(item) != {"candidate", "relation"}:
            raise ValueError("initial competition deferred candidate is malformed")
        proposed = validate_snapshot(item["candidate"], source, source_length)
        if (item["relation"] != occurrence_relation(geometry, proposed)
                or item["relation"]["status"] == "unique"):
            raise ValueError("initial competition deferred occurrence mapping differs")


def protected_ordinary_reference(row, *, source=None):
    """Return the original scope witness, never a newly selected candidate."""
    evidence = competition_evidence(row, source=source)
    if evidence is None or evidence["rule"] != OCCURRENCE_RULE:
        return None
    return evidence["occurrence_protection"]["ordinary_reference"]


def guard_reconstructed_occurrences(proposed_geometry, *original_rows, source=None):
    """Recheck each original observation against a final proposed geometry."""
    seen = set()
    for row in original_rows:
        reference = protected_ordinary_reference(row, source=source)
        if reference is None:
            continue
        key = (reference["representation_sha256"], reference["source_span_sha256"])
        if key in seen:
            continue
        seen.add(key)
        if occurrence_relation(reference["repeat_intervals"], proposed_geometry)["status"] != "unique":
            return False
    return True


def _initial_observations(row):
    if len(row) < 4:
        return
    yield row[3].get(METADATA_KEY)
    basis = row[3].get("classification_basis", {})
    if (basis.get("method") == "low_spacer_diversity_tandem_demotion"
            and isinstance(basis.get("prior_classification_basis"), dict)
            and basis["prior_classification_basis"].get("method") in ("source_backed_possible_reconstruction",
                "source_backed_localized_spacer_context")):
        basis = basis["prior_classification_basis"]
    if basis.get("method") in ("source_backed_candidate_reconstruction",
                               "source_backed_possible_reconstruction"):
        proof = basis["reconstruction"]
        for name in ("anchor_initial_competition", "donor_initial_competition"):
            yield proof.get(name)
    elif basis.get("method") == "source_backed_boundary_uncertain_reconstruction":
        proof = basis["boundary_family"]
        yield proof.get("anchor_initial_competition")
        for member in proof["members"] + proof["additional_peer_representations"]:
            yield member.get("initial_competition")

    elif basis.get("method") == "source_backed_localized_spacer_context":
        for observation in basis["contextual_spacer_context"]["observations"]:
            yield from observation["initial_competition"]
            yield observation["current_scope"].get("evidence")


def ordinary_reference_key(reference):
    return (tuple(map(tuple, reference["repeat_intervals"])), reference["score"],
            reference["representation_sha256"], reference["source_span_sha256"])


def validate_ordinary_witness_inventory(categories, source, *, observations=()):
    """Bind compact witnesses to preserved ordinary rows, once per inventory.

    This checks identities and source spans, not historical ranking or fresh
    feature extraction. Geometry/score indexing avoids hashing every candidate.
    """
    rows = [row for category in categories for values in category.values() for row in values]
    required, deferred = {}, {}
    def collect(evidence):
        if not isinstance(evidence, dict) or evidence.get("rule") != OCCURRENCE_RULE:
            return
        protection = evidence.get("occurrence_protection")
        if not isinstance(protection, dict) or "ordinary_reference" not in protection:
            raise ValueError("initial competition occurrence protection is malformed")
        reference = protection["ordinary_reference"]
        validate_snapshot(reference, source_length=len(source))
        key = ordinary_reference_key(reference)
        if key not in required:
            validate_snapshot(reference, source)
            required[key] = reference
        observations = protection.get("deferred_candidates")
        if not isinstance(observations, list):
            raise ValueError("initial competition deferred inventory is malformed")
        for observation in observations:
            if not isinstance(observation, dict) or "candidate" not in observation:
                raise ValueError("initial competition deferred candidate is malformed")
            snapshot = observation["candidate"]
            validate_snapshot(snapshot, source_length=len(source))
            key = ordinary_reference_key(snapshot)
            if key not in deferred:
                validate_snapshot(snapshot, source)
                deferred[key] = snapshot
    for row in rows:
        for evidence in _initial_observations(row):
            collect(evidence)
    for evidence in observations:
        collect(evidence)
    if not required:
        return frozenset()
    by_location = {}
    for key in required.keys() | deferred.keys():
        by_location.setdefault(key[:2], []).append(key)
    found, found_deferred = set(), set()
    for row in rows:
        keys = by_location.get((_geometry(row[1]), _score(row)), ())
        if keys:
            actual = occurrence_snapshot(row)
            kind = candidate_origin_kind(row[1])
            for key in keys:
                if kind == "ordinary" and actual == required.get(key):
                    found.add(key)
                if kind == "repair_only" and actual == deferred.get(key):
                    found_deferred.add(key)
    if found != set(required):
        raise ValueError("initial competition has no preserved ordinary occurrence witness")
    if found_deferred != set(deferred):
        raise ValueError("initial competition has no retained deferred repair snapshot")
    return frozenset(found)
