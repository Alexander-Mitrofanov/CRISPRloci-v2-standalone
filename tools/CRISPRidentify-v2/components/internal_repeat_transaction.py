"""Single-core, source-wide transaction for pending internal completions."""

from copy import deepcopy
import hashlib
import json

from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
)
from components.fast_run_rescue import candidate_score_fingerprint
from components.internal_local_core_bridge import (
    collect_bridge_inventory, compose_local_core, replay_local_core, _restore,
)
from components.internal_completion_authority import (
    authority_state, current_by_geometry, unrestricted, parent_current_allowed,
)
from components.internal_repeat_completion import (
    InternalSearchContext, _json, _restrictions, replay_internal_attestation,
)
from components.repeat_boundary_repair import classification_restriction, source_geometry
from components.terminal_repair_completion import (
    _geometry, _overlap, _rows, _score, _structural_from_rule,
)


METHOD = "source_backed_internal_repeat_completion"
TRANSACTION_POLICY = "single-core-internal-completion-v1"


def _ordered_subset(short, long):
    intervals = list(zip(*long))
    indices = [intervals.index(interval) if interval in intervals else -1 for interval in zip(*short)]
    return bool(indices) and indices[0] >= 0 and all(a < b for a, b in zip(indices, indices[1:]))


def _retained_anchors(record):
    parent = _geometry(record["parent_geometry"])
    closed = list(zip(*_geometry(record["gap_closed_geometry"])))
    final = list(zip(*_geometry(record["child_geometry"])))
    corrected = [final[closed.index(interval)] for interval in zip(*parent)]
    return tuple(zip(*corrected))


def _core_match(core, record):
    """Exactly one contiguous original-anchor match within attested trims."""
    parent = _geometry(record["parent_geometry"])
    corrected = _retained_anchors(record)
    origins = [parent]
    for authorization in record["parent_authorizations"]:
        if authorization["kind"] == "terminal_repair":
            origins.append(_geometry(authorization["boundary_attestation"]["parent_geometry"]))
    matches = set()
    count = len(core[0])
    for original in origins:
        for offset in range(len(parent[0]) - count + 1):
            if all(original[0][offset+i] <= start <= corrected[0][offset+i]
                   and corrected[1][offset+i] <= end <= original[1][offset+i]
                   for i, (start, end) in enumerate(zip(*core))):
                matches.add(offset)
    return next(iter(matches)) if len(matches) == 1 else None


def complete_internal_repeats(source, categories, pending, structural, score_candidates, search_context=None, boundary_records=None):
    """Freshly score, replay and atomically replace only one compatible core.

    Input rows and pending candidates are immutable. No candidate is added to
    ordinary evaluation or classified solely because its geometry reconstructs.
    """
    def authority_snapshot():
        rows = _rows(categories)
        score_state = tuple((category, key, index, id(row), _score(row[0]),
                             tuple(tuple(getattr(row[1], name)) for name in
                                   ("list_repeat_starts", "list_repeats", "list_repeats_gaped", "list_spacers")))
                            for category, key, index, row in rows)
        authority = ([authority_state(row[1]) for _, _, _, row in rows],
                     [(source_geometry(candidate), authority_state(candidate),
                       getattr(candidate, "internal_repair_record", None)) for candidate in pending],
                     boundary_records)
        return score_state, json.dumps(authority, sort_keys=True, default=repr)
    snapshot = authority_snapshot()
    context = search_context or InternalSearchContext(source)
    if context.source != source:
        raise ValueError("internal transaction search context source differs")
    audit = []
    current = current_by_geometry(categories)
    by_geometry = {}
    blocked = set()
    for proposal in pending:
        record = getattr(proposal, "internal_repair_record", None)
        geometry = source_geometry(proposal)
        try:
            valid_marker = (getattr(proposal, "source_path_classification_eligible", None) is False
                            and _restrictions(proposal) == ["pending_internal_completion"]
                            and not hasattr(proposal, "source_path_provenance"))
        except (ValueError, TypeError):
            valid_marker = False
        if not valid_marker:
            blocked.add(geometry)
            audit.append({"status": "abstained", "reason": "pending_restriction_changed",
                          "child_geometry": geometry, "authority": deepcopy(authority_state(proposal))})
            continue
        try:
            matching_geometry = isinstance(record, dict) and geometry == _geometry(record["child_geometry"])
        except (ValueError, TypeError, KeyError):
            matching_geometry = False
        if not matching_geometry:
            blocked.add(geometry)
            audit.append({"status": "abstained", "reason": "stale_pending_geometry"})
            continue
        by_geometry.setdefault(geometry, []).append(record)
    ready = {}
    for geometry, records in sorted(by_geometry.items()):
        if geometry in blocked or any(not unrestricted(source, candidate) for candidate in current.get(geometry, ())):
            audit.append({"status": "abstained", "reason": "restricted_child_derivation", "child_geometry": geometry,
                          "current_authorities": [deepcopy(authority_state(candidate)) for candidate in current.get(geometry, ())]})
            continue
        try:
            parents_allowed = all(parent_current_allowed(source, _geometry(record["parent_geometry"]), current,
                                                         structural, context.source_sha256) for record in records)
        except (ValueError, TypeError, KeyError):
            parents_allowed = False
        if not parents_allowed:
            audit.append({"status": "abstained", "reason": "restricted_parent_derivation", "child_geometry": geometry})
            continue
        canonical = None
        try:
            for record in records:
                rebuilt = replay_internal_attestation(source, record, structural, context)
                if source_geometry(rebuilt) != geometry:
                    raise ValueError("replay changed pending geometry")
                canonical = rebuilt
        except (ValueError, KeyError, TypeError, AttributeError) as error:
            audit.append({"status": "abstained", "reason": "invalid_attestation", "detail": str(error), "child_geometry": geometry})
            continue
        # No old provenance is copied to new source intervals. This explicit
        # marker describes the reattested new geometry only.
        canonical.source_path_classification_eligible = True
        ready[geometry] = (canonical, records)
    if not ready:
        return categories, audit
    key = "source-internal-completion-batch"
    candidates = [item[0] for item in ready.values()]
    fingerprints = [candidate_score_fingerprint(candidate) for candidate in candidates]
    scores = score_candidates({key: candidates})
    if authority_snapshot() != snapshot:
        raise ValueError("internal scoring mutated category snapshot")
    if set(scores) != {key} or len(scores[key]) != len(candidates):
        raise ValueError("internal scoring changed request cardinality")
    scored = {}
    for geometry, candidate, fingerprint, row in zip(ready, candidates, fingerprints, scores[key]):
        if row[1] is not candidate or candidate_score_fingerprint(candidate) != fingerprint:
            raise ValueError("internal scoring changed candidate identity or geometry")
        _score(row[0])
        if classification_restriction(candidate) or _restrictions(candidate) or not structural.evidence(candidate).passes:
            raise ValueError("internal scoring changed candidate eligibility or structure")
        scored[geometry] = row
    rows = _rows(categories)
    accepted = [item for item in rows if item[0] in (0, 2)]
    by_current_geometry = {}
    for item in rows:
        by_current_geometry.setdefault(source_geometry(item[3][1]), []).append(item[3][1])
    replacements = []
    phase_cache = {}
    bridge_proofs = {}
    for incumbent in accepted:
        category, location, _, old_row = incumbent
        core = source_geometry(old_row[1])
        threshold = 0.75 if category == 0 else 0.5
        if _score(old_row[0]) < threshold or any(
                not unrestricted(source, candidate)
                for candidate in by_current_geometry[core]):
            continue
        compatible = {}
        for geometry, (candidate, records) in ready.items():
            if len(geometry[0]) <= len(core[0]) or scored[geometry][0] < threshold:
                continue
            matches = []
            for record in records:
                if _core_match(core, record) is not None:
                    matches.append(record)
                elif boundary_records is not None and _overlap(core, geometry):
                    try:
                        inventory = collect_bridge_inventory(geometry, boundary_records, current, pending)
                        proof = compose_local_core(source, core, record, inventory, structural, context, phase_cache)
                    except (ValueError, TypeError, KeyError, AttributeError) as error:
                        audit.append({"status": "abstained", "reason": "local_core_composition_rejected",
                                      "detail": str(error), "core_geometry": core, "child_geometry": geometry})
                    else:
                        matches.append(record)
                        bridge_proofs[(core, _json(record))] = proof
            if matches:
                compatible[geometry] = matches
        if not compatible:
            continue
        maxima = [geometry for geometry in compatible if not any(
            geometry != other and _ordered_subset(geometry, other) for other in compatible)]
        if len(maxima) != 1:
            audit.append({"status": "abstained", "reason": "incomparable_completions", "core_geometry": core})
            continue
        geometry = maxima[0]
        if any(other[3] is not old_row and _overlap(geometry, source_geometry(other[3][1])) for other in accepted):
            audit.append({"status": "abstained", "reason": "other_accepted_overlap", "core_geometry": core})
            continue
        if any(other != geometry and _overlap(other, geometry)
               and not _ordered_subset(other, geometry) and not _ordered_subset(geometry, other)
               and scored[other][0] >= threshold for other in ready):
            audit.append({"status": "abstained", "reason": "incompatible_supported_path", "core_geometry": core})
            continue
        record = min(compatible[geometry], key=_json)
        evidence = {
            "policy": TRANSACTION_POLICY, "source_sha256": context.source_sha256,
            "coordinate_system": "zero_based_half_open",
            "core_geometry": core, "completed_geometry": geometry,
            "core_score": float(old_row[0]), "completed_score": float(scored[geometry][0]),
            "core_category_index": category,
            "completed_category_index": 0 if scored[geometry][0] >= 0.75 else 2,
            "repair_attestation": deepcopy(record),
            "compatible_geometries": tuple(sorted(compatible)),
            "accepted_event_geometries": tuple(source_geometry(item[3][1]) for item in accepted),
            "selection_reason": "unique_maximal_complete_internal_path",
        }
        local_proof = bridge_proofs.get((core, _json(record)))
        if local_proof is not None:
            evidence["local_core_composition"] = local_proof
        evidence["previous_completed_rows"] = [
            {"category_index": item[0], "score": float(item[3][0]),
             "authority": authority_state(item[3][1])}
            for item in rows if source_geometry(item[3][1]) == geometry]
        replacements.append((incumbent, scored[geometry], evidence))
    conflicting = {i for i, (_, _, evidence) in enumerate(replacements)
                   for j, (_, _, other) in enumerate(replacements) if i != j
                   and _overlap(evidence["completed_geometry"], other["completed_geometry"])}
    if authority_snapshot() != snapshot:
        raise ValueError("internal completion snapshot changed before commit")
    copied = [{location: list(items) for location, items in mapping.items()} for mapping in categories]
    committed = set()
    for index, (incumbent, new_row, evidence) in enumerate(replacements):
        if index in conflicting:
            audit.append({"status": "abstained", "reason": "planned_replacement_collision"})
            continue
        old_category, location, _, old_row = incumbent
        geometry = evidence["completed_geometry"]
        new_category = evidence["completed_category_index"]
        # Replay at commit, after scoring, against exactly the current scored
        # candidate. Full report validation repeats this check independently.
        replay_internal_completion_basis(source, new_row[1], evidence, new_row[0], new_category, context)
        metadata = dict(new_row[3]) if len(new_row) > 3 else {}
        metadata["classification_basis"] = {"method": METHOD, "internal_repeat_completion": evidence}
        replacement = list(new_row[:3]) + [metadata]
        for mapping in copied:
            for old_location in list(mapping):
                mapping[old_location] = [row for row in mapping[old_location]
                                         if row is not old_row and source_geometry(row[1]) != geometry]
                if not mapping[old_location]:
                    del mapping[old_location]
        copied[new_category].setdefault(location, []).append(replacement)
        lower = 1 if old_row[0] >= 0.75 or (new_category == 0 and old_row[0] >= 0.6) else 3
        retained = old_row
        if len(old_row) > 3 and "classification_basis" in old_row[3]:
            retained_metadata = dict(old_row[3])
            retained_metadata["classification_basis"] = {
                "method": METHOD, "retained_core": True,
                "internal_repeat_completion": evidence,
                "prior_classification_basis": deepcopy(old_row[3]["classification_basis"]),
            }
            retained = list(old_row[:3]) + [retained_metadata]
        copied[lower].setdefault(location, []).append(retained)
        committed.add(geometry)
        audit.append({"status": "committed", "evidence": evidence})
    if sum(len(items) for i in (0, 2) for items in copied[i].values()) != len(accepted):
        raise ValueError("internal completion changed accepted event count")
    for geometry in ready:
        if geometry not in committed:
            audit.append({"status": "pending", "reason": "no_qualified_single_core_transaction",
                          "child_geometry": geometry, "fresh_score": float(scored[geometry][0])})
    return copied, audit


def replay_internal_completion_basis(source, candidate, evidence, score, category, context=None):
    score = _score(score)
    if (evidence.get("policy") != TRANSACTION_POLICY
            or evidence.get("source_sha256") != hashlib.sha256(source.encode("ascii")).hexdigest()
            or evidence.get("selection_reason") != "unique_maximal_complete_internal_path"
            or _score(evidence.get("completed_score")) != score
            or type(category) is not int or category not in (0, 2)
            or type(evidence.get("completed_category_index")) is not int
            or evidence["completed_category_index"] != category):
        raise ValueError("internal completion source, score or category differs")
    record = evidence["repair_attestation"]
    structural = _structural_from_rule(record["structural_rule"])
    child = replay_internal_attestation(source, record, structural, context)
    core = _geometry(evidence["core_geometry"])
    geometry = _geometry(evidence["completed_geometry"])
    build(source, None, *core)
    if (source_geometry(candidate) != geometry or source_geometry(child) != geometry
            or candidate_score_fingerprint(candidate) != candidate_score_fingerprint(child)
            or len(geometry[0]) <= len(core[0]) or not structural.evidence(candidate).passes):
        raise ValueError("internal completion geometry or core relationship differs")
    direct_map = _core_match(core, record)
    local_proof = evidence.get("local_core_composition")
    if direct_map is None:
        if local_proof is None:
            raise ValueError("internal completion lacks an authorized core map")
        context = context or InternalSearchContext(source)
        replay_local_core(source, core, record, local_proof, structural, context)
    elif local_proof is not None:
        raise ValueError("unnecessary local-core authorization")
    for previous in evidence.get("previous_completed_rows", ()):
        previous_candidate = _restore(source, {"geometry": geometry, "authority": previous["authority"]})
        if not unrestricted(source, previous_candidate):
            raise ValueError("completed geometry had a contradictory current derivation")
        _score(previous["score"])
        if type(previous["category_index"]) is not int or previous["category_index"] not in range(5):
            raise ValueError("invalid previous completed category")
    old_category = evidence["core_category_index"]
    if type(old_category) is not int or old_category not in (0, 2):
        raise ValueError("invalid original category")
    threshold = 0.75 if old_category == 0 else 0.5
    if (_score(evidence["core_score"]) < threshold or score < threshold
            or category != (0 if score >= 0.75 else 2)):
        raise ValueError("internal completion changed score band")
    paths = [_geometry(path) for path in evidence["compatible_geometries"]]
    if geometry not in paths or any(not _ordered_subset(path, geometry) for path in paths):
        raise ValueError("internal completion is not a unique maximal path")
    accepted = [_geometry(path) for path in evidence["accepted_event_geometries"]]
    for path in accepted:
        build(source, None, *path)
    if accepted.count(core) != 1 or any(path != core and _overlap(path, geometry) for path in accepted):
        raise ValueError("internal completion overlaps another accepted event")


def replay_retained_internal_core(source, candidate, evidence, score, category, context=None):
    """Retain the old core's genuine band and geometry as lower evidence."""
    if source_geometry(candidate) != _geometry(evidence["core_geometry"]) or _score(score) != _score(evidence["core_score"]):
        raise ValueError("retained internal core geometry or score differs")
    completed_category = evidence["completed_category_index"]
    lower = 1 if score >= 0.75 or (completed_category == 0 and score >= 0.6) else 3
    if type(category) is not int or category != lower:
        raise ValueError("retained internal core category differs")
    completed = build(source, None, *_geometry(evidence["completed_geometry"]))
    replay_internal_completion_basis(source, completed, evidence, evidence["completed_score"], completed_category, context)
