"""Transactional exact-motif completion authorized by terminal repair.

Relative classifier ordering cannot decide how many exact source occurrences
belong to an already accepted core. This deliberately narrow transaction uses
only the bounded repair attestations, while retaining ordinary score bands.
"""

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
import math
from numbers import Real

from components.components_detection_refinement import AdvancedFuzzySearchFilter
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
)
from components.fast_run_rescue import candidate_score_fingerprint
from components.repeat_boundary_repair import (
    POLICY_VERSION, classification_restriction, propose_terminal_repairs,
    source_geometry,
)


METHOD = "source_backed_terminal_repair_completion"
TRANSACTION_POLICY = "attested-exact-core-completion-v1"
CATEGORY_NAMES = ("Bona-fide", "Alternative", "Possible", "Possible Discarded", "Low score")


def _geometry(value):
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError("terminal repair geometry must contain starts and ends")
    starts, ends = tuple(value[0]), tuple(value[1])
    if (not starts or len(starts) != len(ends)
            or any(type(x) is not int for x in starts + ends)
            or any(start < 0 or end <= start for start, end in zip(starts, ends))
            or any(start < end for start, end in zip(starts[1:], ends))):
        raise ValueError("terminal repair geometry has invalid source intervals")
    return starts, ends


def _score(value):
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(float(value)):
        raise ValueError("terminal completion needs a finite numeric model score")
    return float(value)


def _rows(categories):
    return [(category, key, index, row)
            for category, mapping in enumerate(categories)
            for key, rows in mapping.items() for index, row in enumerate(rows)]


def _snapshot(categories):
    return tuple((category, key, index, id(row), _score(row[0]),
                  candidate_score_fingerprint(row[1]), classification_restriction(row[1]))
                 for category, key, index, row in _rows(categories))


def _overlap(first, second):
    return first[0][0] < second[1][-1] and second[0][0] < first[1][-1]


def _subpath(short, long):
    count = len(short[0])
    return any(short == (long[0][offset:offset + count], long[1][offset:offset + count])
               for offset in range(len(long[0]) - count + 1))


def _corridor(core, parent, long):
    count = len(core[0])
    matches = []
    for offset in range(len(parent[0]) - count + 1):
        if all(parent[0][offset + index] <= start <= long[0][offset + index]
               and long[1][offset + index] <= end <= parent[1][offset + index]
               for index, (start, end) in enumerate(zip(*core))):
            matches.append(offset)
    return matches[0] if len(matches) == 1 else None


def _structural_from_rule(rule):
    return AdvancedFuzzySearchFilter(
        min_column_dominance_repeat=rule["column_dominance"],
        min_avg_spacer_length=rule["min_avg_spacer_length"],
        max_spacer_length=rule["max_spacer_length"],
        max_column_dominance_spacer=rule["max_column_dominance_spacer"],
        max_allowed_consecutive_spacers=rule["max_allowed_consecutive_spacers"],
        max_allowed_same_spacers=rule["max_allowed_same_spacers"],
        max_inconsistent_columns=rule["max_inconsistent_columns"],
        min_avg_repeat_length=rule["min_avg_repeat_length"],
        max_avg_repeat_length=rule["max_avg_repeat_length"],
        max_avg_spacer_length=rule["max_avg_spacer_length"],
        min_repeats=rule["min_number_repeats"],
    )


def _replay_attestation(source, record, structural, source_hash):
    if (record.get("policy") != POLICY_VERSION
            or record.get("source_sha256") != source_hash
            or record.get("parent_eligibility") is not True
            or record.get("structural_rule") != vars(structural)
            or record.get("reason") not in (
                "unique_minimum_supported_trim", "child_geometry_already_present")):
        raise ValueError("terminal repair attestation identity or policy differs")
    parent_geometry = _geometry(record["parent_geometry"])
    parent = build(source, None, *parent_geometry)
    children, replay, _ = propose_terminal_repairs(
        source, [parent], structural, source_sha256=source_hash,
    )
    if len(children) != 1 or len(replay) != 1:
        raise ValueError("terminal repair attestation cannot be replayed")
    replay = replay[0]
    for field in ("child_geometry", "left_trim_columns", "right_trim_columns",
                  "consumed_bases", "parent_failures", "resulting_checks", "boundary_evidence"):
        if json.dumps(record.get(field), sort_keys=True) != json.dumps(replay[field], sort_keys=True):
            raise ValueError("terminal repair attestation differs at " + field)
    return parent_geometry, children[0]


def _exact_path(candidate, structural):
    repeats, spacers = candidate.list_repeats, candidate.list_spacers
    motif = repeats[0]
    return (set(motif) <= set("ACGT") and len(set(repeats)) == 1
            and all(set(spacer) <= set("ACGT")
                    and structural.min_avg_spacer_length <= len(spacer) <= structural.max_avg_spacer_length
                    and motif not in spacer for spacer in spacers)
            and structural.evidence(candidate).passes)


@dataclass(frozen=True)
class CompletionPlan:
    snapshot: tuple
    replacements: tuple
    records: tuple


def plan_terminal_repair_completion(source, categories, records, structural):
    """Plan unique maximal exact extensions without changing rows or scores."""
    snapshot = _snapshot(categories)
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    rows = _rows(categories)
    by_geometry = {}
    for item in rows:
        by_geometry.setdefault(source_geometry(item[3][1]), []).append(item)
    accepted = [item for item in rows if item[0] in (0, 2)]
    supported = {}
    audit = []
    seen_attestations = set()
    for record in records:
        if "child_geometry" not in record:
            continue
        parent_geometry = _geometry(record["parent_geometry"])
        child_geometry = _geometry(record["child_geometry"])
        identity = parent_geometry, child_geometry
        if identity in seen_attestations:
            continue
        seen_attestations.add(identity)
        current = by_geometry.get(child_geometry, ())
        if not current:
            continue
        if any(classification_restriction(item[3][1]) for geometry in identity
               for item in by_geometry.get(geometry, ())):
            continue
        # Reject degenerate/long-gap rows cheaply before source-only replay.
        if not any(_exact_path(item[3][1], structural) for item in current):
            continue
        try:
            parent_geometry, canonical = _replay_attestation(source, record, structural, source_hash)
        except (ValueError, KeyError, TypeError):
            audit.append({"status": "abstained", "reason": "invalid_attestation", "child_geometry": child_geometry})
            continue
        if not _exact_path(canonical, structural):
            continue
        fingerprint = candidate_score_fingerprint(canonical)
        matching = [item for item in current
                    if candidate_score_fingerprint(item[3][1]) == fingerprint]
        if not matching or len({float(item[3][0]) for item in matching}) != 1:
            continue
        supported.setdefault(child_geometry, []).append((parent_geometry, matching[0], record))

    replacements = []
    for incumbent in accepted:
        category, key, _, row = incumbent
        core = source_geometry(row[1])
        if any(classification_restriction(item[3][1]) for item in by_geometry[core]):
            continue
        threshold = 0.75 if category == 0 else 0.5
        if row[0] < threshold:
            continue
        compatible = {}
        for long, attestations in supported.items():
            if len(long[0]) < len(core[0]) + 2:
                continue
            matches = [(parent, item, record) for parent, item, record in attestations
                       if float(item[3][0]) >= threshold
                       and _corridor(core, parent, long) is not None]
            if matches:
                compatible[long] = matches
        if not compatible:
            continue
        maxima = [long for long in compatible
                  if not any(long != other and _subpath(long, other) for other in compatible)]
        if len(maxima) != 1:
            audit.append({"status": "abstained", "reason": "incomparable_completions", "core_geometry": core})
            continue
        long = maxima[0]
        if any(other[3] is not row and _overlap(long, source_geometry(other[3][1]))
               for other in accepted):
            audit.append({"status": "abstained", "reason": "other_accepted_overlap", "core_geometry": core})
            continue
        # An attested, model-qualified overlapping phase outside the nesting
        # relation is unresolved evidence, even if it was not a corridor match.
        if any(other != long and _overlap(other, long)
               and not _subpath(other, long) and not _subpath(long, other)
               and any(float(item[3][0]) >= threshold for _, item, _ in members)
               for other, members in supported.items()):
            audit.append({"status": "abstained", "reason": "incompatible_attested_path", "core_geometry": core})
            continue
        parent, completed, attestation = min(compatible[long], key=lambda item: item[0])
        evidence = {
            "policy": TRANSACTION_POLICY, "source_sha256": source_hash,
            "coordinate_system": "zero_based_half_open",
            "core_geometry": core, "parent_geometry": parent, "completed_geometry": long,
            "core_score": float(row[0]), "completed_score": float(completed[3][0]),
            "core_category_index": category,
            "completed_category_index": 0 if completed[3][0] >= 0.75 else 2,
            "repair_attestation": deepcopy(attestation),
            "compatible_geometries": tuple(sorted(compatible)),
            "accepted_event_geometries": tuple(source_geometry(item[3][1]) for item in accepted),
            "eligibility_checked": True,
            "selection_reason": "unique_maximal_attested_exact_core_extension",
        }
        replacements.append((incumbent, completed, evidence))
    # Independently planned extensions must not create a collision together.
    conflicting = {i for i, (_, _, evidence) in enumerate(replacements)
                   for j, (_, _, other) in enumerate(replacements) if i != j
                   and _overlap(evidence["completed_geometry"], other["completed_geometry"])}
    replacements = tuple(item for index, item in enumerate(replacements) if index not in conflicting)
    return CompletionPlan(snapshot, replacements, tuple(audit))


def apply_terminal_repair_completion(source, categories, plan, structural, score_candidates):
    """Rescore selected geometry and commit only against the same snapshot."""
    if _snapshot(categories) != plan.snapshot:
        raise ValueError("terminal completion category snapshot changed")
    audit = list(plan.records)
    if not plan.replacements:
        return categories, audit
    requests = {}
    for incumbent, completed, _ in plan.replacements:
        requests.setdefault(incumbent[1], []).append(completed[3][1])
    fresh = score_candidates(requests)
    if _snapshot(categories) != plan.snapshot:
        raise ValueError("terminal completion scoring mutated category snapshot")
    if set(fresh) != set(requests) or any(len(fresh[key]) != len(rows) for key, rows in requests.items()):
        raise ValueError("terminal completion scoring changed cardinality")
    fresh_by_id = {}
    for key, candidates in requests.items():
        for candidate, row in zip(candidates, fresh[key]):
            _score(row[0])
            if row[1] is not candidate:
                raise ValueError("terminal completion scoring changed candidate identity")
            fresh_by_id[id(candidate)] = row
    copied = [{key: list(rows) for key, rows in category.items()} for category in categories]
    accepted_before = sum(len(rows) for index in (0, 2) for rows in categories[index].values())
    for incumbent, completed, evidence in plan.replacements:
        old_category, key, _, old_row = incumbent
        new_row = fresh_by_id[id(completed[3][1])]
        threshold = 0.75 if old_category == 0 else 0.5
        if new_row[0] < threshold or not _exact_path(new_row[1], structural):
            audit.append({"status": "abstained", "reason": "fresh_score_or_structure_rejected", "core_geometry": evidence["core_geometry"]})
            continue
        evidence = deepcopy(evidence)
        evidence["completed_score"] = float(new_row[0])
        new_category = 0 if new_row[0] >= 0.75 else 2
        evidence["completed_category_index"] = new_category
        metadata = dict(new_row[3]) if len(new_row) > 3 else {}
        metadata["classification_basis"] = {"method": METHOD, "terminal_repair_completion": evidence}
        replacement = list(new_row[:3]) + [metadata]
        # Remove all duplicate memberships of the selected scored source path.
        long = source_geometry(new_row[1])
        for category in copied:
            for location in list(category):
                category[location] = [row for row in category[location]
                                      if row is not old_row and source_geometry(row[1]) != long]
                if not category[location]:
                    del category[location]
        copied[new_category].setdefault(key, []).append(replacement)
        lower = 1 if old_row[0] >= 0.75 or (new_category == 0 and old_row[0] >= 0.6) else 3
        copied[lower].setdefault(key, []).append(old_row)
        audit.append({"status": "committed", "evidence": evidence})
    accepted_after = sum(len(rows) for index in (0, 2) for rows in copied[index].values())
    if accepted_after != accepted_before:
        raise ValueError("terminal completion changed accepted event count")
    return copied, audit


def replay_terminal_completion_basis(source, candidate, evidence, score, category):
    """Validate serialized source geometry and every reconstructible guard."""
    score = _score(score)
    completed_score = _score(evidence.get("completed_score"))
    core_score = _score(evidence.get("core_score"))
    if (evidence.get("policy") != TRANSACTION_POLICY
            or evidence.get("source_sha256") != hashlib.sha256(source.encode("ascii")).hexdigest()
            or evidence.get("eligibility_checked") is not True
            or evidence.get("selection_reason") != "unique_maximal_attested_exact_core_extension"
            or completed_score != score
            or evidence.get("completed_category_index") != category
            or type(category) is not int or category not in (0, 2)
            or type(evidence.get("completed_category_index")) is not int):
        raise ValueError("terminal completion source, score or category differs")
    core = _geometry(evidence["core_geometry"])
    long = _geometry(evidence["completed_geometry"])
    structural = _structural_from_rule(evidence["repair_attestation"]["structural_rule"])
    parent, canonical = _replay_attestation(source, evidence["repair_attestation"], structural,
                                          evidence["source_sha256"])
    build(source, None, *core)
    if (parent != _geometry(evidence["parent_geometry"])
            or source_geometry(candidate) != long or source_geometry(canonical) != long
            or _corridor(core, parent, long) is None
            or len(long[0]) < len(core[0]) + 2 or not _exact_path(candidate, structural)):
        raise ValueError("terminal completion source relationship differs")
    old_category = evidence["core_category_index"]
    threshold = 0.75 if old_category == 0 else 0.5
    if (type(old_category) is not int or old_category not in (0, 2)
            or core_score < threshold or score < threshold
            or category != (0 if score >= 0.75 else 2)):
        raise ValueError("terminal completion score band differs")
    paths = [_geometry(path) for path in evidence["compatible_geometries"]]
    if long not in paths or any(not _subpath(path, long) for path in paths):
        raise ValueError("terminal completion is not the unique maximal path")
    accepted = [_geometry(path) for path in evidence["accepted_event_geometries"]]
    for path in accepted:
        build(source, None, *path)
    if accepted.count(core) != 1 or any(path != core and _overlap(path, long) for path in accepted):
        raise ValueError("terminal completion overlaps another accepted event")
