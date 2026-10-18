"""Checked whole-path context for one localized duplicate-spacer block.

The original structural failure remains recorded. A separate transaction can
publish a Possible observation after independently checking both maximal
flanks and the complete generated path. It never promotes a cropped flank.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import math

from components.architecture_rescue_coding_context import SourceCDS
from components.candidate_origins import candidate_origin_evidence, inherit_candidate_origins, candidate_origin_kind
from components.components_detection_refinement import CrisprCandidate
from components.fast_run_rescue import compact_feature_values
from components.initial_competition import (
    _initial_observations, initial_competition_groups, ordinary_reference_key, competition_evidence,
    RULE as INITIAL_SCOPE_RULE, OCCURRENCE_RULE,
    validate_ordinary_witness_inventory,
)
from components.likely_tandem import evaluate_likely_tandem
from components.localized_spacer_context import (
    RULE, MINIMUM_OWN_SCORE, SCORE_TIE_TOLERANCE, PendingSpacerContext,
    candidate_capture, digest, encoded, localized_block, plain,
)
from components.occurrence_protection import occurrence_snapshot, validate_snapshot
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.reconstruction_boundary_family import (
    boundary_family_summary, describe_geometry_set, family_occurrence_evidence, geometry_id,
)
from components.reconstruction_transaction import (
    STRUCTURAL_PARAMETERS, _inventory_state, _restore, _snapshot,
)
from components.reconstruction_validation import (
    _coding_evidence, _known_negative_annotations, needs_coding_context,
)
from components.reconstruction_validator import validate_candidate_reconstruction


METHOD = "source_backed_localized_spacer_context"
EVIDENCE_FILENAME = "contextual_spacer_evidence.json"
OBSERVATION_KEY = "contextual_spacer_observation"
POLICY = {
    "rule": RULE, "minimum_own_score": MINIMUM_OWN_SCORE,
    "minimum_flank_repeats": 3, "confidence_ceiling": "Possible",
    "score_tie_tolerance": SCORE_TIE_TOLERANCE,
    "rank": "maximum_repeat_count_then_full_own_score",
    "flanks": "maximal_outside_context; correlated_checks_not_votes",
    "boundary_family": "complete_top_set_shared_starts_actual_member_view",
}


def _geometry(candidate):
    return [[start, start+len(repeat)] for start, repeat in
            zip(candidate.list_repeat_starts, candidate.list_repeats)]


def _flank(candidate, indexes):
    a, b = indexes
    child = CrisprCandidate(list(candidate.list_repeats[a:b]),
                            list(candidate.list_repeats_gaped[a:b]),
                            list(candidate.list_spacers[a:b-1]),
                            list(candidate.list_repeat_starts[a:b]))
    inherit_candidate_origins(child, (candidate,))
    return child


def _observation(row, source_hash):
    return {"snapshot": _snapshot(row), "metadata": plain(row[3] if len(row) == 4 else {}),
            "origins": candidate_origin_evidence(row[1], source_hash)}


def _validation(source, row, structural, source_cds, *, full=False, history=()):
    candidate = row[1]
    source_check = validate_candidate_reconstruction(source, candidate)
    alignment_valid = all(
        len(values) == len(candidate.list_repeats) and all(
            g.replace(" ", "").replace("-", "") == r for g, r in zip(values, candidate.list_repeats))
        for values in (candidate.list_repeats_gaped, candidate.list_gaped_repeats))
    structure = structural.evidence(candidate).as_dict()
    failures = [item["predicate"] for item in structure["predicate_results"]
                if item["status"] == "failed"]
    tandem = evaluate_likely_tandem(candidate.list_spacers)
    if tandem.likely_tandem:
        tandem_status = "failed"
    elif tandem.reason_code in ("spacer_count_below_policy_minimum", "spacer_count_above_policy_maximum"):
        tandem_status = "inapplicable"
    elif tandem.complete_link_threshold_satisfied is False:
        tandem_status = "passed"
    else:
        tandem_status = "not_run"
    evidence = {
        "source_reconstruction": source_check.as_dict(), "alignment_source_valid": alignment_valid, "structural_filter": structure,
        "candidate_eligibility": getattr(candidate, "source_path_classification_eligible", True),
        "known_negative_annotations": _known_negative_annotations(row, *history),
        "likely_tandem": {"status": tandem_status, "decision": tandem.as_dict()},
        "coding_tandem": _coding_evidence(source, candidate, source_cds, True),
        "canonical_sequence": all(set(s) <= set("ACGT") for s in
                                  candidate.list_repeats + candidate.list_spacers),
        "own_score": float(row[0]), "minimum_own_score": MINIMUM_OWN_SCORE,
        "scope": "enumerated_checks_only; not_biological_correctness",
        "original_structural_exception": "identical_spacer_burden" if full else None,
    }
    evidence["qualified"] = bool(
        source_check.valid and alignment_valid and evidence["candidate_eligibility"] is True
        and failures == (["identical_spacer_burden"] if full else [])
        and float(row[0]) >= MINIMUM_OWN_SCORE and evidence["canonical_sequence"]
        and all(evidence[key]["status"] in ("passed", "inapplicable") for key in
                ("known_negative_annotations", "likely_tandem", "coding_tandem")))
    return plain(evidence)


def _components(records):
    groups, end = [], -1
    for record in sorted(records, key=lambda record: (
            record["fresh"]["repeat_intervals"][0][0], record["ordinal"])):
        intervals = record["fresh"]["repeat_intervals"]
        if intervals[0][0] >= end:
            groups.append([])
            end = -1
        groups[-1].append(record)
        end = max(end, intervals[-1][1])
    return [sorted(group, key=lambda record: record["ordinal"]) for group in groups]


def _negative_signature(record):
    return tuple(record["validation"][name]["status"] for name in
                 ("known_negative_annotations", "likely_tandem", "coding_tandem"))


def _descriptor(source_hash, members, display_id):
    geometries = [member["fresh"]["repeat_intervals"] for member in members]
    neutral = describe_geometry_set(geometries)
    display = next((member for member in members if member["observation_id"] == display_id), None)
    if display is None:
        raise ValueError("context display member missing")
    scores = [member["fresh"]["score"] for member in members]
    return {
        **neutral, "rule": RULE, "coordinate_system": "zero_based_half_open",
        "status": "boundary_uncertain" if len(members) > 1 else "single_checked_geometry",
        "coordinate_semantics": "actual_admissible_member_view",
        "family_id": digest([RULE, source_hash, neutral["admissible_geometry_ids"]]),
        "confidence_category": "Possible", "accepted_core_present": False,
        "member_count": len(members), "repeat_count": neutral["supported_repeat_count"],
        "display_observation_id": display_id,
        "display_geometry_id": geometry_id(display["fresh"]["repeat_intervals"]),
        "display_own_score": display["fresh"]["score"],
        "member_score_min": min(scores), "member_score_max": max(scores),
        "own_score_range": [min(scores), max(scores)],
        "flank_checks": "diagnostic_correlated_context; no_score_borrowing",
    }


def _choose(source_hash, group, barriers):
    qualified = [record for record in group if record["qualified"]]
    if not qualified:
        return {"group_id": digest([r["observation_id"] for r in group]),
                "observation_ids": [r["observation_id"] for r in group], "top_peer_ids": [],
                "maximum_repeat_count": None, "maximum_own_score": None,
                "action": "abstained", "reason": "no_qualified_context"}, []
    count = max(len(record["fresh"]["repeat_intervals"]) for record in qualified)
    count_peers = [record for record in group if len(record["fresh"]["repeat_intervals"]) == count]
    score = max(record["fresh"]["score"] for record in qualified
                if len(record["fresh"]["repeat_intervals"]) == count)
    top = [record for record in count_peers
           if abs(record["fresh"]["score"]-score) <= SCORE_TIE_TOLERANCE]
    decision = {
        "group_id": digest([record["observation_id"] for record in group]),
        "observation_ids": [record["observation_id"] for record in group],
        "top_peer_ids": [record["observation_id"] for record in top],
        "maximum_repeat_count": count, "maximum_own_score": score,
        "action": "abstained", "reason": None,
    }
    if any(not record["qualified"] for record in top):
        return {**decision, "reason": "unsafe_same_rank_peer"}, []
    by_geometry = {}
    for record in top:
        by_geometry.setdefault(geometry_id(record["fresh"]["repeat_intervals"]), []).append(record)
    members = [min(by_geometry[key], key=lambda record: record["observation_id"]) for key in sorted(by_geometry)]
    try:
        descriptor = _descriptor(source_hash, members, members[0]["observation_id"])
    except ValueError:
        return {**decision, "reason": "different_start_tied_geometries"}, []
    if len({_negative_signature(record) for record in top}) != 1:
        return {**decision, "reason": "negative_applicability_disagreement"}, []
    start, end = descriptor["union_span"]
    collisions = [barrier for barrier in barriers
                  if start < barrier["interval"][1] and end > barrier["interval"][0]]
    if collisions:
        return {**decision, "reason": "accepted_envelope_collision", "collisions": collisions}, []
    references = [reference for record in top for reference in record["ordinary_references"]]
    guard = family_occurrence_evidence(references, [member["fresh"]["repeat_intervals"] for member in members])
    if guard["status"] != "preserved":
        return {**decision, "reason": "original_occurrence_conflict", "ordinary_witness_guard": guard}, []
    return {**decision, "action": "selected", "reason": "checked_localized_spacer_context",
            "descriptor": descriptor, "ordinary_witness_guard": guard,
            "member_observation_ids": [member["observation_id"] for member in members]}, members


def _barriers(categories, *, exclude=None):
    barriers = []
    for category in (0, 2):
        for key, rows in categories[category].items():
            for row in rows:
                if row is exclude:
                    continue
                summary = boundary_family_summary(row) or contextual_spacer_summary(row)
                geometry = _geometry(row[1])
                interval = summary["union_span"] if summary else [geometry[0][0], geometry[-1][1]]
                barriers.append({"category_index": category, "locus_key": plain(key),
                                 "interval": interval, "representation_sha256": digest(_snapshot(row))})
    return barriers


def _competitive_rows(categories, structural, pending_candidates):
    # Owners and repair bridges are qualified current observations, not fragments.
    spans = sorted((p.candidate.list_repeat_starts[0],
                    p.candidate.list_repeat_starts[-1] + len(p.candidate.list_repeats[-1]))
                   for p in pending_candidates)
    # Keep complete ordinary connected components touching a pending envelope.
    # Geometry-only sweep avoids running structural checks on unrelated rows.
    possible = [row for category in categories for rows in category.values() for row in rows
            if row[0] >= .50
            and getattr(row[1], "source_path_classification_eligible", True) is True
            and not (len(row) == 4 and OBSERVATION_KEY in row[3])]
    ordered = sorted(possible, key=lambda row: row[1].list_repeat_starts[0])
    components, end = [], -1
    for row in ordered:
        a, b = _geometry(row[1])[0][0], _geometry(row[1])[-1][1]
        if a >= end:
            components.append([])
            end = -1
        components[-1].append(row)
        end = max(end, b)
    return [row for group in components
            if any(group[0][1].list_repeat_starts[0] < b and
                   max(_geometry(r[1])[-1][1] for r in group) > a for a, b in spans)
            for row in group if structural(row[1]) is not None]


def _scope_observation(row, source_hash):
    return {"snapshot": occurrence_snapshot(row),
            "origins": candidate_origin_evidence(row[1], source_hash)}


def _scope_checks(ordinary, pending_rows, existing_repairs=()):
    repair = [row for row in pending_rows if candidate_origin_kind(row[1]) == "repair_only"]
    groups, conflicts = initial_competition_groups(ordinary + list(existing_repairs) + repair)
    annotations = {id(row[1]): row[3]["initial_competition"] for row in conflicts}
    for group in groups:
        for row in group.rows:
            if len(row) == 4 and "initial_competition" in row[3]:
                annotations[id(row[1])] = row[3]["initial_competition"]
    result = {}
    for row in pending_rows:
        evidence = annotations.get(id(row[1]))
        reference = (evidence.get("occurrence_protection", {}).get("ordinary_reference")
                     if evidence else None)
        result[id(row[1])] = {"origin_kind": candidate_origin_kind(row[1]),
            "evidence": evidence, "ordinary_references": [reference] if reference else [],
            "blocked": bool(evidence and evidence["decision"] in
                            ("repair_scope_conflict", "repair_occurrence_conflict"))}
    return result


def _history_evidence(prior, source=None):
    observations = []
    for row in prior:
        competition_evidence(row, source=source)
        for evidence in _initial_observations(row):
            if evidence is not None:
                if evidence.get("rule") not in (INITIAL_SCOPE_RULE, OCCURRENCE_RULE):
                    raise ValueError("context historical competition rule is unknown")
                observations.append(evidence)
    return observations


def _references(evidence, current_scope):
    references = list(current_scope["ordinary_references"])
    for item in evidence:
        reference = item.get("occurrence_protection", {}).get("ordinary_reference")
        if reference is not None:
            references.append(reference)
    return list({encoded(reference): reference for reference in references}.values())


def _qualified(record):
    history_blocked = any(item["decision"] in ("repair_scope_conflict", "repair_occurrence_conflict")
                          for item in record["initial_competition"])
    return (record["validation"]["qualified"] and not record["current_scope"]["blocked"]
            and not history_blocked and all(record["flanks"][role]["validation"]["qualified"]
                                           for role in ("left", "right")))


def rescue_contextual_spacer_arrays(source, categories, pending_candidates, parameters,
                                    score_candidates, get_source_cds):
    """Commit Possible views; retain every captured full observation as Low."""
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    structural_parameters = {key: parameters[key] for key in STRUCTURAL_PARAMETERS}
    structural = build_advanced_fuzzy_search_filter(structural_parameters)
    pending_candidates = tuple(pending_candidates)
    if not pending_candidates:
        return categories, {"schema": {"name": "CRISPRidentify-contextual-spacer-evidence", "version": "1.0.0"},
            "rule": RULE, "source_sha256": source_hash, "policy": dict(POLICY), "confidence_category": "Possible",
            "observations": [], "decisions": [], "actions": {}, "structural_parameters": structural_parameters,
            "ordinary_scope_rows": [], "repair_scope_rows": [], "accepted_barriers": [],
            "stats": {"pending_observations": 0, "fresh_candidates_scored": 0, "selections_committed": 0}}
    validate_ordinary_witness_inventory(categories, source)
    before = _inventory_state(categories)
    history = {}
    for category in categories:
        for rows in category.values():
            for row in rows:
                history.setdefault(geometry_id(_geometry(row[1])), []).append(row)
    pending_candidates = tuple(pending_candidates)
    if len({id(p) for p in pending_candidates}) != len(pending_candidates):
        raise ValueError("duplicate pending observation object; capture each observation separately")
    if len({id(p.candidate) for p in pending_candidates}) != len(pending_candidates):
        raise ValueError("pending observations must own distinct candidate copies")
    competitive = _competitive_rows(categories, structural, pending_candidates)
    ordinary = [row for row in competitive if candidate_origin_kind(row[1]) == "ordinary"]
    repairs = [row for row in competitive if candidate_origin_kind(row[1]) == "repair_only"]
    ordinary_scope_rows = [_scope_observation(row, source_hash) for row in ordinary]
    repair_scope_rows = [_scope_observation(row, source_hash) for row in repairs]
    candidates, plans = [], []
    for ordinal, pending in enumerate(pending_candidates):
        if not isinstance(pending, PendingSpacerContext):
            raise TypeError("context input must contain captured pending observations")
        pending.validate_current(source_hash)
        block = localized_block(pending.candidate, structural)
        if encoded(block) != pending.block_json:
            raise ValueError("pending context block or structural policy differs")
        candidate = pending.candidate
        roles = {"full": candidate,
                 "left": _flank(candidate, block["left_repeat_indexes"]),
                 "right": _flank(candidate, block["right_repeat_indexes"])}
        candidates.extend(roles.values())
        plans.append((ordinal, pending, block, roles))
    unique_candidates = list({id(candidate): candidate for candidate in candidates}.values())
    key = "localized-spacer-context-fresh-score-batch"
    scoring_capture = {id(c): encoded(candidate_capture(c, source_hash)) for c in unique_candidates}
    scored = score_candidates({key: unique_candidates}) if unique_candidates else {}
    if _inventory_state(categories) != before:
        raise ValueError("context scoring mutated the original categories")
    if unique_candidates and (set(scored) != {key} or len(scored[key]) != len(unique_candidates)):
        raise ValueError("context scoring changed candidate cardinality")
    fresh = {}
    for candidate, row in zip(unique_candidates, scored.get(key, ())):
        if row[1] is not candidate or type(row[0]) is bool or not math.isfinite(float(row[0])) or not 0 <= row[0] <= 1:
            raise ValueError("context scoring changed candidate identity/order or score")
        if encoded(candidate_capture(candidate, source_hash)) != scoring_capture[id(candidate)]:
            raise ValueError("context scoring changed full/flank representation or provenance")
        compact_feature_values(row)
        fresh[id(candidate)] = row
    cds = tuple(get_source_cds()) if any(needs_coding_context(c) for c in unique_candidates) else None
    scope = _scope_checks(ordinary, [fresh[id(plan[3]["full"])] for plan in plans], repairs)
    records, row_by_id = [], {}
    output = [{key: list(rows) for key, rows in category.items()} for category in categories]
    for ordinal, pending, block, roles in plans:
        pending.validate_current(source_hash)
        full_row = fresh[id(roles["full"])]
        prior = history.get(geometry_id(_geometry(roles["full"])), ())
        evidence = _history_evidence(prior, source)
        current_scope = scope[id(roles["full"])]
        references = _references(evidence, current_scope)
        record = {
            "ordinal": ordinal, "locus_key": plain(pending.locus_key),
            "capture": pending.observation, "block": block,
            "fresh": _snapshot(full_row), "fresh_metadata": plain(full_row[3] if len(full_row) == 4 else {}), "historical_rows": [_observation(row, source_hash) for row in prior],
            "initial_competition": evidence, "current_scope": current_scope,
            "ordinary_references": references,
            "validation": _validation(source, full_row, structural, cds, full=True, history=prior),
            "flanks": {role: {"fresh": _snapshot(fresh[id(roles[role])]),
                              "fresh_metadata": plain(fresh[id(roles[role])][3] if len(fresh[id(roles[role])]) == 4 else {}),
                              "validation": _validation(source, fresh[id(roles[role])], structural, cds)}
                       for role in ("left", "right")},
        }
        record["qualified"] = _qualified(record)
        record["observation_id"] = digest([source_hash, ordinal, record["locus_key"], record["capture"]])
        record = plain(record)
        records.append(record)
        metadata = deepcopy(full_row[3]) if len(full_row) == 4 else {}
        metadata[OBSERVATION_KEY] = {"observation_id": record["observation_id"], "capture": record["capture"]}
        metadata.setdefault("classification_basis", {"method": "model_score_band"})
        retained = list(full_row[:3]) + [metadata]
        output[4].setdefault(pending.locus_key, []).append(retained)
        row_by_id[record["observation_id"]] = retained
    barriers = _barriers(categories)
    planned = [_choose(source_hash, group, barriers) + (group,) for group in _components(records)]
    collisions = set()
    for index, (decision, _members, _group) in enumerate(planned):
        if decision["action"] != "selected":
            continue
        a, b = decision["descriptor"]["union_span"]
        for other, (second, _unused, _rest) in enumerate(planned[index+1:], index+1):
            if second["action"] == "selected":
                c, d = second["descriptor"]["union_span"]
                if a < d and b > c:
                    collisions.update((index, other))
    decisions = []
    for index, (decision, members, group) in enumerate(planned):
        if index in collisions:
            decision = {**decision, "action": "abstained", "reason": "symmetric_context_proposal_collision"}
        if decision["action"] == "selected":
            proof = plain({"rule": RULE, "source_sha256": source_hash, "policy": POLICY,
                           "structural_parameters": structural_parameters,
                           "confidence_category": "Possible", "accepted_core_present": False,
                           "observations": group, "members": members, "selection": decision,
                           "descriptor": decision["descriptor"], "accepted_barriers": barriers,
                           "ordinary_scope_rows": ordinary_scope_rows, "repair_scope_rows": repair_scope_rows})
            basis = {"method": METHOD, "contextual_spacer_context": proof, "evidence_sha256": digest(proof)}
            display_id = decision["descriptor"]["display_observation_id"]
            display = row_by_id[display_id]
            metadata = deepcopy(display[3])
            metadata["classification_basis"] = basis
            locus_key = tuple(next(member["locus_key"] for member in members if member["observation_id"] == display_id))
            output[2].setdefault(locus_key, []).append(list(display[:3]) + [metadata])
        decisions.append(decision)
    if _inventory_state(categories) != before:
        raise ValueError("context transaction mutated original categories")
    audit = {"schema": {"name": "CRISPRidentify-contextual-spacer-evidence", "version": "1.0.0"},
             "rule": RULE, "source_sha256": source_hash, "policy": POLICY,
             "confidence_category": "Possible", "observations": records, "decisions": decisions,
             "structural_parameters": structural_parameters, "ordinary_scope_rows": ordinary_scope_rows,
             "repair_scope_rows": repair_scope_rows,
             "accepted_barriers": barriers,
             "actions": dict(Counter(d["action"] for d in decisions)),
             "stats": {"pending_observations": len(records), "fresh_candidates_scored": len(unique_candidates),
                       "selections_committed": sum(d["action"] == "selected" for d in decisions)}}
    return output, plain(audit)


def _proof(basis):
    if not isinstance(basis, dict) or basis.get("method") != METHOD:
        return None
    if "accepted_family_tandem_review" in basis:
        if not isinstance(basis["accepted_family_tandem_review"], dict):
            raise ValueError("context tandem wrapper malformed")
        basis = {key: value for key, value in basis.items() if key != "accepted_family_tandem_review"}
    if set(basis) != {"method", "contextual_spacer_context", "evidence_sha256"}:
        raise ValueError("context basis keys differ")
    proof = basis["contextual_spacer_context"]
    if not isinstance(proof, dict) or set(proof) != {
            "rule", "source_sha256", "policy", "structural_parameters", "confidence_category",
            "accepted_core_present", "observations", "members", "selection", "descriptor",
            "accepted_barriers", "ordinary_scope_rows", "repair_scope_rows"}:
        raise ValueError("context proof keys differ")
    if (digest(proof) != basis["evidence_sha256"] or proof["rule"] != RULE
            or encoded(proof["policy"]) != encoded(POLICY)
            or proof["confidence_category"] != "Possible" or proof["accepted_core_present"] is not False):
        raise ValueError("context rule, confidence or proof digest differs")
    descriptor = _descriptor(proof["source_sha256"], proof["members"],
                             proof["descriptor"]["display_observation_id"])
    if encoded(descriptor) != encoded(proof["descriptor"]):
        raise ValueError("context descriptor differs")
    return proof


def contextual_spacer_summary(record):
    basis = record[3].get("classification_basis") if len(record) == 4 else None
    proof = _proof(basis)
    if proof is None:
        return None
    descriptor = proof["descriptor"]
    if (geometry_id(_geometry(record[1])) != descriptor["display_geometry_id"]
            or encoded(record[0]) != encoded(descriptor["display_own_score"])):
        raise ValueError("context displayed row differs")
    return {**descriptor, "method": METHOD, "proof_sha256": basis["evidence_sha256"]}


def report_contextual_spacer_summary(array):
    proof = _proof(array.get("classification_basis"))
    if proof is None:
        if "contextual_spacer_context" in array:
            raise ValueError("context summary lacks proof")
        return None
    result = {**proof["descriptor"], "method": METHOD,
              "proof_sha256": array["classification_basis"]["evidence_sha256"]}
    if any(type(r["source_interval"][key]) is not int for r in array["repeats"] for key in ("start", "end")):
        raise ValueError("context canonical coordinate types differ")
    geometry = [[r["source_interval"]["start"]-1, r["source_interval"]["end"]] for r in array["repeats"]]
    if (array["category"] != "Possible"
            or encoded(array["source_interval"]["start"]) != encoded(geometry[0][0] + 1)
            or encoded(array["source_interval"]["end"]) != encoded(geometry[-1][1])
            or encoded(array.get("contextual_spacer_context")) != encoded(result)
            or geometry_id(geometry) != result["display_geometry_id"]
            or encoded(array["certainty_score"]) != encoded(result["display_own_score"])
            or encoded(array["repeat_count"]) != encoded(result["repeat_count"])):
        raise ValueError("context canonical summary differs")
    return result


def contextual_spacer_member_intervals(record):
    """Return only complete actual admitted paths, never end combinations."""
    if contextual_spacer_summary(record) is None:
        return ()
    proof = _proof(record[3]["classification_basis"])
    return tuple(tuple(map(tuple, member["fresh"]["repeat_intervals"])) for member in proof["members"])


def _restored(source, snapshot):
    row = _restore(source, snapshot)
    if encoded(_snapshot(row)) != encoded(snapshot):
        raise ValueError("context snapshot typed round-trip differs")
    return row


def _cds(recorded):
    evidence = recorded["coding_tandem"]
    if type(evidence["source_cds_assessed"]) is not bool:
        raise ValueError("context CDS assessment type differs")
    return (tuple(SourceCDS(**item) for item in evidence["considered_source_cds"])
            if evidence["source_cds_assessed"] else None)


def validate_contextual_recorded_metrics(source, record, structural_parameters):
    """Replay source, stored-vector and enumerated checks in exported evidence.

    This artifact API cannot authenticate unexported origin birth objects or
    independently re-extract features. The caller separately joins complete
    retained observations, scope witnesses, and emitted decisions. Production
    publication uses the stronger live-inventory replay below.
    """
    structural = build_advanced_fuzzy_search_filter(structural_parameters)
    full = _restored(source, record["fresh"])
    full.append(deepcopy(record["fresh_metadata"]))
    capture = record["capture"]
    expected = candidate_capture(full[1], hashlib.sha256(source.encode("ascii")).hexdigest())
    for key in ("repeat_intervals", "repeats", "spacers", "gapped_repeats", "gapped_repeats_alias",
                "consensus", "eligible", "eligibility_reason", "score_status"):
        if encoded(capture[key]) != encoded(expected[key]):
            raise ValueError("context recorded capture/source representation differs")
    block = localized_block(full[1], structural)
    if block is None or encoded(block) != encoded(record["block"]):
        raise ValueError("context recorded localized block differs")
    history = []
    for observation in record["historical_rows"]:
        original = _restored(source, observation["snapshot"])
        original.append(deepcopy(observation["metadata"]))
        history.append(original)
    check = _validation(source, full, structural, _cds(record["validation"]), full=True, history=history)
    if encoded(check) != encoded(record["validation"]):
        raise ValueError("context recorded full metric replay differs")
    if set(record["flanks"]) != {"left", "right"}:
        raise ValueError("context recorded flank roles differ")
    for role in ("left", "right"):
        item = record["flanks"][role]
        if set(item) != {"fresh", "fresh_metadata", "validation"}:
            raise ValueError("context recorded flank keys differ")
        fresh = _restored(source, item["fresh"])
        candidate = _flank(full[1], block[role + "_repeat_indexes"])
        reconstructed = [fresh[0], candidate, fresh[2], deepcopy(item["fresh_metadata"])]
        if encoded(_snapshot(reconstructed)) != encoded(item["fresh"]):
            raise ValueError("context recorded flank source/recipe differs")
        check = _validation(source, reconstructed, structural, _cds(item["validation"]))
        if encoded(check) != encoded(item["validation"]):
            raise ValueError("context recorded flank metric replay differs")
    return {"source_and_stored_metrics": "passed", "feature_extraction_repeated": False,
            "raw_origin_authenticity": "requires_original_inventory"}


def _replay_records(source, categories, records, parameters, ordinary_scope_rows, repair_scope_rows):
    """Validate current retained observations once, without rerunning features.

    Scores are replayed from their own stored vectors. This does not certify
    independently extracted features or unrecorded upstream history.
    """
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    structural = build_advanced_fuzzy_search_filter(parameters)
    rows = [row for category in categories for values in category.values() for row in values]
    history, pending, scope_index = {}, {}, {}
    for row in rows:
        meta = row[3] if len(row) == 4 else {}
        if OBSERVATION_KEY in meta:
            # The accepted display is a second view; the Low observation is
            # the immutable retained owner of the original candidate payload.
            continue
        history.setdefault(geometry_id(_geometry(row[1])), []).append(row)
    for values in categories[4].values():
        for row in values:
            if len(row) == 4 and OBSERVATION_KEY in row[3]:
                oid = row[3][OBSERVATION_KEY]["observation_id"]
                if oid in pending:
                    raise ValueError("context retained observation ID is duplicated")
                pending[oid] = row
    scope_rows = ordinary_scope_rows + repair_scope_rows
    required_scope = {encoded(observation) for observation in scope_rows}
    for row in rows:
        if not (len(row) == 4 and OBSERVATION_KEY in row[3]):
            # Only candidate locations named by the compact witnesses require
            # source-span/representation hashing.
            if any(_geometry(row[1]) == item["snapshot"]["repeat_intervals"] and
                   row[0] == item["snapshot"]["score"] for item in scope_rows):
                key = encoded(_scope_observation(row, source_hash))
                if key in required_scope:
                    scope_index.setdefault(key, []).append(row)
    scope_actual = []
    for observation in scope_rows:
        validate_snapshot(observation["snapshot"], source)
        found = scope_index.get(encoded(observation), ())
        if not found:
            raise ValueError("context ordinary owner has no retained source-bound row")
        row = found[0]
        if row[0] < .50 or getattr(row[1], "source_path_classification_eligible", True) is not True or structural(row[1]) is None:
            raise ValueError("context ordinary owner is no longer qualified")
        scope_actual.append(row)
    ordinary, repairs = scope_actual[:len(ordinary_scope_rows)], scope_actual[len(ordinary_scope_rows):]
    if (any(candidate_origin_kind(row[1]) != "ordinary" for row in ordinary)
            or any(candidate_origin_kind(row[1]) != "repair_only" for row in repairs)):
        raise ValueError("context current scope origin kind differs")
    actual_full, actual_history = [], []
    seen = set()
    for record in records:
        if set(record) != {"ordinal", "locus_key", "capture", "block", "fresh", "fresh_metadata",
                           "historical_rows", "initial_competition", "current_scope", "ordinary_references",
                           "validation", "flanks", "qualified", "observation_id"}:
            raise ValueError("context observation keys differ")
        ordinal, oid = record["ordinal"], record["observation_id"]
        if type(ordinal) is not int or ordinal < 0 or oid in seen:
            raise ValueError("context observation ordinal/identity differs")
        seen.add(oid)
        if digest([source_hash, ordinal, record["locus_key"], record["capture"]]) != oid:
            raise ValueError("context capture identity differs")
        if oid not in pending:
            raise ValueError("context has no retained Low observation")
        retained = pending[oid]
        expected_metadata = deepcopy(record["fresh_metadata"])
        expected_metadata.setdefault("classification_basis", {"method": "model_score_band"})
        expected_metadata[OBSERVATION_KEY] = {"observation_id": oid, "capture": record["capture"]}
        if (encoded(retained[3]) != encoded(expected_metadata)
                or encoded(_snapshot(retained)) != encoded(record["fresh"])
                or encoded(candidate_capture(retained[1], source_hash)) != encoded(record["capture"])):
            raise ValueError("context retained observation/capture differs")
        _restored(source, record["fresh"])
        candidate = retained[1]
        block = localized_block(candidate, structural)
        if encoded(block) != encoded(record["block"]) or block is None:
            raise ValueError("context localized block differs")
        full = list(retained[:3]) + [deepcopy(record["fresh_metadata"])]
        prior = history.get(geometry_id(_geometry(candidate)), ())
        if encoded([_observation(row, source_hash) for row in prior]) != encoded(record["historical_rows"]):
            raise ValueError("context historical observations differ")
        if encoded(_history_evidence(prior, source)) != encoded(record["initial_competition"]):
            raise ValueError("context historical competition differs")
        actual_full.append(full)
        actual_history.append(prior)
        expected = _validation(source, full, structural, _cds(record["validation"]), full=True, history=prior)
        if encoded(expected) != encoded(record["validation"]):
            raise ValueError("context full validation differs")
        if set(record["flanks"]) != {"left", "right"}:
            raise ValueError("context flank roles differ")
        for role in ("left", "right"):
            item = record["flanks"][role]
            if set(item) != {"fresh", "fresh_metadata", "validation"}:
                raise ValueError("context flank keys differ")
            fresh = _restored(source, item["fresh"])
            flank = _flank(candidate, block[role + "_repeat_indexes"])
            reconstructed = [fresh[0], flank, fresh[2], deepcopy(item["fresh_metadata"])]
            if encoded(_snapshot(reconstructed)) != encoded(item["fresh"]):
                raise ValueError("context flank recipe differs")
            check = _validation(source, reconstructed, structural, _cds(item["validation"]))
            if encoded(check) != encoded(item["validation"]):
                raise ValueError("context flank validation differs")
    scope = _scope_checks(ordinary, actual_full, repairs)
    for record, row in zip(records, actual_full):
        current = scope[id(row[1])]
        if (encoded(current) != encoded(record["current_scope"])
                or encoded(_references(record["initial_competition"], current)) != encoded(record["ordinary_references"])
                or encoded(_qualified(record)) != encoded(record["qualified"])):
            raise ValueError("context scope/witness/qualification differs")
    validate_ordinary_witness_inventory(categories, source, observations=(
        evidence for record in records for evidence in
        [*record["initial_competition"], record["current_scope"]["evidence"]]))
    return actual_full


def _replay_decisions(source_hash, records, barriers):
    planned = [_choose(source_hash, group, barriers)[0] for group in _components(records)]
    collisions = set()
    for i, first in enumerate(planned):
        if first["action"] != "selected":
            continue
        a, b = first["descriptor"]["union_span"]
        for j, second in enumerate(planned[i+1:], i+1):
            if second["action"] == "selected":
                c, d = second["descriptor"]["union_span"]
                if a < d and b > c:
                    collisions.update((i, j))
    return [{**decision, "action": "abstained", "reason": "symmetric_context_proposal_collision"}
            if i in collisions else decision for i, decision in enumerate(planned)]


def _bind_barriers(categories, recorded):
    available = Counter(encoded(item) for item in _barriers(categories))
    expected = Counter(encoded(item) for item in recorded)
    if expected - available:
        raise ValueError("context original accepted barrier is missing")


def validate_contextual_spacer_basis(source, candidate, basis, row_score, category_index,
                                     current_row, categories, *, ordinary_witness_keys=None):
    proof = _proof(basis)
    if proof is None or type(category_index) is not int or category_index != 2:
        raise ValueError("context must retain the Possible confidence ceiling")
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    if source_hash != proof["source_sha256"] or current_row is None or current_row[1] is not candidate:
        raise ValueError("context source/current row differs")
    _replay_records(source, categories, proof["observations"], proof["structural_parameters"],
                    proof["ordinary_scope_rows"], proof["repair_scope_rows"])
    _bind_barriers(categories, proof["accepted_barriers"])
    decision, members = _choose(source_hash, proof["observations"], proof["accepted_barriers"])
    if (decision["action"] != "selected" or encoded(decision) != encoded(proof["selection"])
            or encoded(members) != encoded(proof["members"])):
        raise ValueError("context complete winning peer/member set differs")
    display = next(member for member in members if member["observation_id"] == proof["descriptor"]["display_observation_id"])
    if (encoded(_snapshot(current_row)) != encoded(display["fresh"])
            or encoded(row_score) != encoded(display["fresh"]["score"])
            or encoded(candidate_capture(candidate, source_hash)) != encoded(display["capture"])):
        raise ValueError("context actual-member display differs")
    a, b = proof["descriptor"]["union_span"]
    for barrier in _barriers(categories, exclude=current_row):
        # Report callers may pass an equivalent pre-quarantine row instead of
        # the same object; remove exactly one matching display, not siblings.
        if a < barrier["interval"][1] and b > barrier["interval"][0]:
            raise ValueError("context final accepted-envelope collision")
    if ordinary_witness_keys is not None:
        for member in proof["observations"]:
            for reference in member["ordinary_references"]:
                if ordinary_reference_key(reference) not in ordinary_witness_keys:
                    raise ValueError("context original ordinary witness missing from validation context")
    return basis


def _nested_proofs(basis):
    if not isinstance(basis, dict):
        return
    proof = _proof(basis)
    if proof is not None:
        yield proof
    elif "prior_classification_basis" in basis:
        yield from _nested_proofs(basis["prior_classification_basis"])


def validate_contextual_spacer_audit(report, audit, *, provenance_categories=None,
                                    likely_tandem_records=None, source_sequence=None):
    """Reconcile every recorded observation and decision with current rows.

    Quarantine dispositions are supplied only after their independent report
    extension validation. This is same-pass consistency, not authenticity of
    a coherently rewritten historical ledger.
    """
    published = {}
    for array in report["arrays"]:
        for proof in _nested_proofs(array.get("classification_basis")):
            key = digest(proof["selection"])
            if key in published:
                raise ValueError("context publication duplicates a selected proof")
            published[key] = proof
    if audit is None and not published:
        return
    if (not isinstance(audit, dict) or set(audit) != {"schema", "rule", "source_sha256", "policy",
            "confidence_category", "observations", "decisions", "actions", "stats",
            "structural_parameters", "ordinary_scope_rows", "repair_scope_rows", "accepted_barriers"}
            or encoded(audit["schema"]) != encoded({"name": "CRISPRidentify-contextual-spacer-evidence", "version": "1.0.0"})
            or audit["rule"] != RULE or audit["source_sha256"] != report["source"]["sha256"]
            or encoded(audit["policy"]) != encoded(POLICY) or audit["confidence_category"] != "Possible"):
        raise ValueError("context audit rule/source/policy differs")
    records = audit["observations"]
    if (not isinstance(records, list) or any(type(r["ordinal"]) is not int or r["ordinal"] != i for i, r in enumerate(records))
            or len({r["observation_id"] for r in records}) != len(records)):
        raise ValueError("context complete ordered observation accounting differs")
    if records:
        if source_sequence is None or provenance_categories is None or hashlib.sha256(source_sequence.encode("ascii")).hexdigest() != audit["source_sha256"]:
            raise ValueError("context complete audit replay requires source and full provenance inventory")
        actual = _replay_records(source_sequence, provenance_categories, records,
                                 audit["structural_parameters"], audit["ordinary_scope_rows"], audit["repair_scope_rows"])
        # Rebuild the qualified current owner universe, including components
        # bridged by pending repairs, without treating new Low/views as owners.
        class View:
            def __init__(self, candidate):
                self.candidate = candidate
        competitive = _competitive_rows(provenance_categories,
            build_advanced_fuzzy_search_filter(audit["structural_parameters"]), [View(row[1]) for row in actual])
        for kind, name in (("ordinary", "ordinary_scope_rows"), ("repair_only", "repair_scope_rows")):
            expected = [_scope_observation(row, audit["source_sha256"]) for row in competitive
                        if candidate_origin_kind(row[1]) == kind]
            if encoded(expected) != encoded(audit[name]):
                raise ValueError("context qualified owner/repair-bridge inventory is incomplete")
        _bind_barriers(provenance_categories, audit["accepted_barriers"])
    elif audit["ordinary_scope_rows"] or audit["repair_scope_rows"] or audit["accepted_barriers"]:
        raise ValueError("empty context audit has unexpected inventories")
    decisions = _replay_decisions(audit["source_sha256"], records, audit["accepted_barriers"])
    if encoded(decisions) != encoded(audit["decisions"]):
        raise ValueError("context full decision replay differs")
    selected = {digest(d): d for d in decisions if d["action"] == "selected"}
    expected_stats = {"pending_observations": len(records), "fresh_candidates_scored": 3*len(records),
                      "selections_committed": len(selected)}
    if (encoded(audit["stats"]) != encoded(expected_stats)
            or encoded(audit["actions"]) != encoded(dict(Counter(d["action"] for d in decisions)))):
        raise ValueError("context complete counters differ")
    if set(published).difference(selected):
        raise ValueError("context published selection is not in emitted audit")
    missing = set(selected).difference(published)
    observed = {}
    if provenance_categories is not None:
        for index, category in enumerate(provenance_categories):
            for locus, rows in category.items():
                for ordinal, row in enumerate(rows):
                    for proof in _nested_proofs(row[3].get("classification_basis") if len(row) == 4 else None):
                        observed[digest(proof["selection"])] = (index, locus, ordinal, row, proof)
    quarantine = [record.as_dict() if hasattr(record, "as_dict") else record for record in (likely_tandem_records or ())]
    for key in missing:
        if key not in observed:
            raise ValueError("context selected output is missing without retained disposition")
        index, locus, ordinal, row, proof = observed[key]
        valid = [record for record in quarantine if record["action"] == "quarantined"
                 and record["source_category_index"] == index == 2
                 and encoded(record["locus_key"]) == encoded(locus)
                 and encoded(record["row_ordinal"]) == encoded(ordinal)
                 and encoded(record["score"]) == encoded(row[0])
                 and encoded([[r["source_interval"]["start"], r["source_interval"]["end"]]
                              for r in record["source_architecture"]["repeats"]]) == encoded(_geometry(row[1]))]
        if len(valid) != 1:
            raise ValueError("context projected selection has no exact quarantine disposition")
        published[key] = proof
    lookup = {r["observation_id"]: r for r in records}
    for key, proof in published.items():
        if (proof["source_sha256"] != audit["source_sha256"]
                or encoded(proof["selection"]) != encoded(selected[key])
                or encoded(proof["observations"]) != encoded([lookup[oid] for oid in selected[key]["observation_ids"]])
                or encoded(proof["ordinary_scope_rows"]) != encoded(audit["ordinary_scope_rows"])
                or encoded(proof["repair_scope_rows"]) != encoded(audit["repair_scope_rows"])
                or encoded(proof["accepted_barriers"]) != encoded(audit["accepted_barriers"])
                or encoded(proof["structural_parameters"]) != encoded(audit["structural_parameters"])):
            raise ValueError("context proof and complete emitted audit differ")
