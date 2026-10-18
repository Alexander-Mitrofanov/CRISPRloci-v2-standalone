"""Observational terminal alternatives; never selection or confidence authority.

The supplied final inventory is the completeness boundary. No score, negative
predictor, generation history, or biological locus is inferred by this module.
Public replay checks recorded values; production replay can additionally bind
the entire extension to the actual source and current inventory.
"""
from bisect import bisect_left
from collections import Counter
import hashlib
import json
import math
from types import SimpleNamespace

from components.candidate_evidence import (
    _encoded, _identity, _intervals, _scope_key, candidate_observation,
)
from components.initial_competition import (
    OCCURRENCE_RULE, RULE as INITIAL_RULE, _initial_observations,
    competition_evidence, validate_ordinary_witness_inventory,
)
from components.likely_tandem import global_edit_distance
from components.likely_tandem_transaction import validated_row_tags
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.occurrence_protection import checked_geometry, occurrence_relation, validate_snapshot
from components.reconstruction_boundary_family import report_boundary_family_summary
from components.contextual_spacer_rescue import report_contextual_spacer_summary
from components.reconstruction_selection import ReconstructionPolicy
from components.reconstruction_validation import _known_negative_annotations
from components.report_model import _json_safe


KEY = "competing_terminal_extents"
SCHEMA = {"name": "CRISPRidentify-competing-terminal-extents", "version": "1.0.0"}
RULE = "exact_terminal_competing_extent_v1"
SEMANTICS = {
    "accepted_science": "unchanged; observations are not additional accepted arrays or votes",
    "scope": "same legacy search key is routing context, not biological locus identity",
    "scores": "existing_final_row; no independent feature or score refresh",
    "negatives": "not_recorded is unknown; absence of recorded failures is not clearance",
    "history": "current supplied observations only; no upstream generation replay",
}


def _plain(value):
    return json.loads(_encoded(value))


def _same(left, right):
    return _encoded(left) == _encoded(right)


def _sha(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _consensus(gapped):
    if not gapped or len({len(r) for r in gapped}) != 1:
        raise ValueError("competing extent alignment columns differ")
    columns = []
    for column in zip(*gapped):
        ranked = [c for c, _ in Counter(column).most_common() if c != "-"]
        if not ranked:
            raise ValueError("competing extent alignment has an all-deletion column")
        columns.append(ranked[0])
    return "".join(columns).replace(" ", "").replace("+", "")


def exact_terminal_relation(anchor, proposed):
    """Exact complete contiguous matching; no coordinate tolerance or cutoff."""
    anchor, proposed = checked_geometry(anchor), checked_geometry(proposed)
    if len(proposed) <= len(anchor):
        return None
    start = bisect_left(proposed, anchor[0])
    if proposed[start:start + len(anchor)] != anchor:
        return None
    end = start + len(anchor)
    return {"anchor_to_candidate_indexes": list(range(start, end)),
            "added_left_indexes": list(range(start)),
            "added_right_indexes": list(range(end, len(proposed)))}


def _geometry(array):
    return checked_geometry([[r["source_interval"]["start"] - 1,
                              r["source_interval"]["end"]] for r in array["repeats"]])


def _accepted(report):
    arrays = [a for a in report["arrays"] if a["category"] in ("Bona-fide", "Possible")]
    if len({a["id"] for a in arrays}) != len(arrays):
        raise ValueError("competing extent accepted identities are not unique")
    result = {}
    for array in arrays:
        geometry = _geometry(array)
        family = report_boundary_family_summary(array)
        contextual = report_contextual_spacer_summary(array)
        if contextual is not None:
            family = contextual
        result[array["id"]] = (array, geometry, family["union_span"] if family else
                               [geometry[0][0], geometry[-1][1]])
    return result


def _snapshot(row, category_index, position, key, source, source_hash, names,
              structural, parameters, canonical_ids):
    record = candidate_observation(
        row, category_index=category_index, category_position=position, key=key,
        source_sequence=source, source_hash=source_hash, feature_names=names,
        structural=structural, parameters=parameters, canonical_ids=canonical_ids)
    candidate = row[1]
    geometry = record["repeat_intervals"]
    if not candidate.consensus or candidate.consensus != _consensus(candidate.list_gaped_repeats):
        raise ValueError("competing extent accepted reference is stale")
    references = {}
    for observation in _initial_observations(row):
        if observation is not None and observation.get("rule") == OCCURRENCE_RULE:
            reference = observation["occurrence_protection"]["ordinary_reference"]
            validate_snapshot(reference, source)
            references[_identity(reference)] = reference
    record.update(
        repeat_sequences=list(candidate.list_repeats),
        spacer_sequences=list(candidate.list_spacers),
        consensus=candidate.consensus,
        row_metadata=_json_safe(row[3] if len(row) > 3 else {}),
        source_span_sha256=_sha(source[geometry[0][0]:geometry[-1][1]]),
        ordinary_witnesses=[references[k] for k in sorted(references)],
        recorded_negative_reasons=_known_negative_annotations(row)["recorded_failures"])
    record = _plain(record)
    record["snapshot_sha256"] = _identity(record)
    return record


def _snapshot_row(record):
    return [record["score"], SimpleNamespace(
        list_repeat_starts=[a for a, _ in record["repeat_intervals"]],
        list_repeats=record["repeat_sequences"], list_repeats_gaped=record["gapped_repeats"],
        list_spacers=record["spacer_sequences"]), [], record["row_metadata"]]


def _check_snapshot(record, source_metadata):
    """Validate the recorded snapshot without asserting an unseen inventory."""
    if not isinstance(record, dict) or "snapshot_sha256" not in record:
        raise ValueError("competing extent snapshot is missing")
    if record["snapshot_sha256"] != _identity({k: v for k, v in record.items() if k != "snapshot_sha256"}):
        raise ValueError("competing extent snapshot digest differs")
    geometry = checked_geometry(record["repeat_intervals"], source_metadata["length"])
    repeats, spacers, gapped = (record[k] for k in ("repeat_sequences", "spacer_sequences", "gapped_repeats"))
    if (type(record["score"]) not in (int, float) or not math.isfinite(record["score"])
            or not 0 <= record["score"] <= 1 or type(record["classification_eligible"]) is not bool
            or type(record["category_position"]) is not int or record["category_position"] < 0
            or type(record["category_index"]) is not int
            or type(record["source_valid"]) is not bool or record["source_valid"] is not True
            or type(record["accepted"]) is not bool
            or type(record["repeat_count"]) is not int or record["repeat_count"] != len(geometry)
            or len(repeats) != len(geometry) or len(gapped) != len(repeats)
            or len(spacers) != len(repeats) - 1):
        raise ValueError("competing extent snapshot cardinality/score differs")
    for i, ((a, b), repeat, aligned) in enumerate(zip(geometry, repeats, gapped)):
        if (not isinstance(repeat, str) or len(repeat) != b - a
                or not isinstance(aligned, str) or aligned.replace(" ", "").replace("-", "") != repeat
                or (i < len(spacers) and (not isinstance(spacers[i], str)
                    or len(spacers[i]) != geometry[i + 1][0] - b))):
            raise ValueError("competing extent source/alignment lengths differ")
    span_string = "".join(r + (spacers[i] if i < len(spacers) else "") for i, r in enumerate(repeats))
    representation = [[a for a, _ in geometry], repeats, gapped, spacers]
    origins = record["candidate_origins"]
    if (record["source_span"] != [geometry[0][0], geometry[-1][1]]
            or record["source_span_sha256"] != _sha(span_string)
            or record["source_geometry_id"] != _identity([source_metadata["sha256"], geometry])
            or record["alignment_id"] != _identity(gapped)
            or origins["current_representation_sha256"] != _identity(representation)
            or any(b["source_sha256"] != source_metadata["sha256"] for b in origins["repair_births"])):
        raise ValueError("competing extent source/representation binding differs")
    if (set(origins) != {"kind", "current_representation_sha256", "ordinary_birth_sha256", "repair_births"}
            or origins["kind"] not in ("ordinary", "repair_only")
            or (origins["kind"] == "ordinary") != (origins["ordinary_birth_sha256"] is not None)
            or (origins["kind"] == "repair_only" and not origins["repair_births"])):
        raise ValueError("competing extent origin kind differs")
    row = _snapshot_row(record)
    metadata = record["row_metadata"]
    basis = metadata.get("classification_basis", {})
    review = basis.get("accepted_family_tandem_review")
    expected_negative = {
        "tags": list(validated_row_tags(row)), "accepted_family_tandem_review": review,
        "accepted_family_tandem_review_status": "recorded" if review is not None else "not_recorded",
        "final_classification_method": basis.get("method")}
    if (not _same(record["negative_evidence"], expected_negative)
            or record["classification_basis_method"] != basis.get("method")
            or not _same(record.get("initial_competition"), metadata.get("initial_competition"))):
        raise ValueError("competing extent advertised metadata differs")
    competition_evidence(row, source_length=source_metadata["length"])
    references = {}
    for observation in _initial_observations(row):
        if observation is None:
            continue
        if observation.get("rule") not in (INITIAL_RULE, OCCURRENCE_RULE):
            raise ValueError("competing extent initial observation rule is unknown")
        if observation["rule"] == OCCURRENCE_RULE:
            protection = observation["occurrence_protection"]
            reference = protection["ordinary_reference"]
            references[_identity(reference)] = reference
            for check in [{"candidate": protection["candidate"], "relation": protection["relation"]}] + protection["deferred_candidates"]:
                if not _same(check["relation"], occurrence_relation(reference["repeat_intervals"], check["candidate"]["repeat_intervals"])):
                    raise ValueError("competing extent original occurrence relationship is not typed exact")
    if not _same(record["ordinary_witnesses"], [references[k] for k in sorted(references)]):
        raise ValueError("competing extent original witness set differs")
    if (record["consensus"] != _consensus(gapped)
            or record["recorded_negative_reasons"] != _known_negative_annotations(
                [record["score"], None, [], record["row_metadata"]])["recorded_failures"]):
        raise ValueError("competing extent reference/recorded negatives differ")
    expected_cid = _identity([record["source_geometry_id"], record["alignment_id"], record["score"],
        record["feature_sets"], record["classification_eligible"], record["category_index"], record["comparison_scope_id"]])
    if (record["candidate_id"] != expected_cid
            or record["observation_id"] != _identity([expected_cid, record["category_index"], record["category_position"]])
            or record["comparison_scope_id"] != _identity([source_metadata["sha256"], record["legacy_locus_key"]])):
        raise ValueError("competing extent observation identity differs")
    feature_sets = record["feature_sets"]
    feature_values = ({f["name"]: f["value"] for f in feature_sets[0]["features"]}
                      if len(feature_sets) == 1 else {})
    count_consistent = (feature_values["number_repeats"] == len(geometry)
                        if "number_repeats" in feature_values else None)
    if (not _same(record["feature_values"], feature_values)
            or not _same(record["feature_repeat_count_consistent"], count_consistent)
            or record["score_origin"] != "existing_final_row"):
        raise ValueError("competing extent advertised feature/score semantics differ")
    for reference in record["ordinary_witnesses"]:
        validate_snapshot(reference, source_length=source_metadata["length"])
    return geometry


def _diagnostics(anchor_id, anchor, candidate, relation, accepted):
    geometry = candidate["repeat_intervals"]
    span = (geometry[0][0], geometry[-1][1])
    collisions = [{"array_id": aid, "source_span": list(other_span)}
                  for aid, (_array, _geometry_, other_span) in sorted(accepted.items())
                  if aid != anchor_id and max(span[0], other_span[0]) < min(span[1], other_span[1])]
    references = {_identity(r): r for r in anchor["ordinary_witnesses"] + candidate["ordinary_witnesses"]}
    checks = [{"reference": references[key], "relation": occurrence_relation(
        references[key]["repeat_intervals"], geometry)} for key in sorted(references)]
    reasons = []
    if collisions:
        reasons.append("overlaps_other_accepted_array")
    if any(check["relation"]["status"] != "unique" for check in checks):
        reasons.append("original_occurrence_conflict")
    for role, record in (("anchor", anchor), ("candidate", candidate)):
        decision = record.get("initial_competition", {}).get("decision")
        if decision in ("repair_scope_conflict", "repair_occurrence_conflict"):
            reasons.append(role + "_" + decision)
        if record["recorded_negative_reasons"]:
            reasons.append(role + "_recorded_negative")
    reference = anchor["consensus"]
    identities = []
    for index in relation["added_left_indexes"] + relation["added_right_indexes"]:
        repeat = candidate["repeat_sequences"][index]
        distance = global_edit_distance(reference, repeat)
        denominator = max(len(reference), len(repeat))
        identities.append({"candidate_index": index, "edit_distance": distance,
                           "denominator": denominator, "identity": 1 - distance / denominator})
    return {"observation_id": candidate["observation_id"],
            "status": "blocked_observation" if reasons else "unaccepted_observation",
            "block_reasons": sorted(reasons), "collisions": collisions,
            "original_occurrences": {"status": ("no_recorded_witness" if not checks else
                "conflict" if any(c["relation"]["status"] != "unique" for c in checks) else "preserved"),
                "checks": checks},
            # These describe the BF/Alternative lane only. A remaining PD row
            # is not evidence of any particular Possible-lane decision; that
            # separate audit records actual fresh checks and abstentions.
            "routing_exclusions": ["anchor_not_bona_fide", "candidate_not_alternative"] +
                (["stored_score_below_current_donor_floor"] if candidate["score"] < ReconstructionPolicy().minimum_candidate_score else []),
            "identity_diagnostic": {"reference_rule": "accepted_consensus_diagnostic_only",
                "reference_sequence": reference, "admission_cutoff": None, "added_repeat_checks": identities}}


def build_competing_terminal_extents(report, *, source_sequence, accepted_categories,
                                    provenance_categories, parameters, feature_names):
    """Build complete observational links from the supplied inventory, or None.

    Legacy mode and APIs without a full inventory omit the optional extension.
    No search result or ranking is mutated, and no missing inventory means zero.
    """
    if parameters.get("param_reconstruction_policy") != "conservative" or provenance_categories is None:
        return None
    if len(accepted_categories) != 5 or len(provenance_categories) != 5:
        raise ValueError("competing extents require all five category slots")
    source = report["source"]
    source_metadata = {k: source[k] for k in ("id", "sha256", "length")}
    if _sha(source_sequence) != source["sha256"] or len(source_sequence) != source["length"]:
        raise ValueError("competing extent source differs from the canonical report")
    accepted = _accepted(report)
    validate_ordinary_witness_inventory(provenance_categories, source_sequence)
    structural = build_advanced_fuzzy_search_filter(parameters)
    canonical_ids = {}
    for aid, (array, geometry, _span) in accepted.items():
        canonical_ids.setdefault((array["category"], _identity([source["sha256"], geometry])), []).append(aid)
    by_position = {a["category_position"]: (aid, a) for aid, (a, _, _) in accepted.items() if a["category"] == "Possible"}
    current = {}
    position = 0
    for key in sorted(accepted_categories[2]):
        for row in accepted_categories[2][key]:
            if position not in by_position:
                raise ValueError("competing extent current Possible row lacks canonical anchor")
            aid, array = by_position[position]
            if (not _same(_scope_key(key), array["legacy_interval_key"])
                    or not _same(_intervals(row[1]), _geometry(array))
                    or not _same(float(row[0]), array["certainty_score"])):
                raise ValueError("competing extent canonical anchor differs from current row")
            current.setdefault(key, {}).setdefault(row[1].list_repeat_starts[0], []).append((aid, array, position, row))
            position += 1
    if position != len(by_position):
        raise ValueError("competing extent canonical Possible inventory differs")
    anchors, observations, grouped = {}, {}, {}
    position = 0
    for key in sorted(provenance_categories[3]):
        for row in provenance_categories[3][key]:
            row_position = position
            position += 1
            eligible = getattr(row[1], "source_path_classification_eligible", True)
            if type(eligible) is not bool:
                raise ValueError("competing extent eligibility must be boolean")
            if key not in current or not eligible:
                continue
            geometry = checked_geometry(_intervals(row[1]), len(source_sequence))
            pairs = [(aid, array, pos, anchor_row, relation)
                     for start, _ in geometry
                     for aid, array, pos, anchor_row in current[key].get(start, ())
                     if (relation := exact_terminal_relation(_intervals(anchor_row[1]), geometry)) is not None]
            if not pairs:
                continue
            record = _snapshot(row, 3, row_position, key, source_sequence, source["sha256"],
                               feature_names, structural, parameters, canonical_ids)
            if record["structural_filter"]["decision"] != "passed":
                continue
            observations[record["observation_id"]] = record
            for aid, array, pos, anchor_row, relation in pairs:
                if aid not in anchors:
                    anchors[aid] = {"array_sha256": _identity(array), "observation": _snapshot(
                        anchor_row, 2, pos, key, source_sequence, source["sha256"], feature_names,
                        structural, parameters, canonical_ids)}
                group = grouped.setdefault((aid, record["source_geometry_id"]), {
                    "anchor_array_id": aid, "source_geometry_id": record["source_geometry_id"],
                    "relation": relation, "observations": []})
                group["observations"].append(_diagnostics(aid, anchors[aid]["observation"], record, relation, accepted))
    links = []
    for key in sorted(grouped):
        link = grouped[key]
        link["observations"].sort(key=lambda d: d["observation_id"])
        link["link_id"] = _identity([RULE, source["sha256"], *key])
        links.append(link)
    result = {"schema": SCHEMA, "rule": RULE, "coordinate_system": "zero_based_half_open",
              "source": source_metadata, "semantics": SEMANTICS,
              "coverage": "complete_for_supplied_final_possible_and_possible_discarded_inventory",
              "anchors": anchors, "observations": observations, "links": links,
              "counts": {"accepted_anchors": len(anchors), "competing_geometries": len({k[1] for k in grouped}),
                         "candidate_observations": len(observations), "links": len(links),
                         "additional_accepted_arrays": 0, "additional_accepted_repeats": 0}}
    result = _plain(result)
    validate_competing_terminal_extents(report, result)
    return result


def validate_competing_terminal_extents(report, extension=None, **context):
    """Replay recorded links; with context also bind exact current inventory."""
    if extension is None:
        extension = report.get("extensions", {}).get(KEY)
    if extension is None:
        return None
    if context:
        expected = build_competing_terminal_extents(report, **context)
        if not _same(extension, expected):
            raise ValueError("competing extent differs from source/current inventory")
    if not isinstance(extension, dict):
        raise ValueError("competing extent must be a mapping")
    if (set(extension) != {"schema", "rule", "coordinate_system", "source", "semantics", "coverage",
                           "anchors", "observations", "links", "counts"}
            or not _same(extension["schema"], SCHEMA) or extension["rule"] != RULE
            or extension["coordinate_system"] != "zero_based_half_open"
            or extension["coverage"] != "complete_for_supplied_final_possible_and_possible_discarded_inventory"
            or not _same(extension["semantics"], SEMANTICS)
            or not _same(extension["source"], {k: report["source"][k] for k in ("id", "sha256", "length")})):
        raise ValueError("competing extent envelope differs")
    accepted = _accepted(report)
    for aid, anchor in extension["anchors"].items():
        if aid not in accepted or accepted[aid][0]["category"] != "Possible":
            raise ValueError("competing extent anchor is not an accepted Possible")
        array = accepted[aid][0]
        record = anchor["observation"]
        geometry = _check_snapshot(record, extension["source"])
        if (set(anchor) != {"array_sha256", "observation"} or anchor["array_sha256"] != _identity(array)
                or not _same(geometry, _geometry(array)) or not _same(record["score"], array["certainty_score"])
                or not _same(record["legacy_locus_key"], array["legacy_interval_key"])
                or not _same(record["feature_sets"], array["feature_sets"])
                or record["consensus"] != array["consensus_repeat"]
                or record["repeat_sequences"] != [r["sequence"] for r in array["repeats"]]
                or record["source_span_sha256"] != _sha(array["source_sequence"])
                or record["category_position"] != array["category_position"]
                or record["category"] != "Possible" or record["category_index"] != 2 or record["accepted"] is not True):
            raise ValueError("competing extent anchor binding differs")
    for oid, record in extension["observations"].items():
        _check_snapshot(record, extension["source"])
        if (oid != record["observation_id"] or record["category"] != "Possible Discarded"
                or record["category_index"] != 3 or record["classification_eligible"] is not True
                or record["structural_filter"]["decision"] != "passed" or record["accepted"] is not False):
            raise ValueError("competing extent candidate is outside the disclosure predicate")
    seen, used, anchors, gids = set(), set(), set(), set()
    for link in extension["links"]:
        aid, gid = link["anchor_array_id"], link["source_geometry_id"]
        if aid not in extension["anchors"] or (aid, gid) in seen or not link["observations"]:
            raise ValueError("competing extent link has absent/repeated anchor or geometry")
        anchor = extension["anchors"][aid]["observation"]
        expected_observations = []
        for diagnostic in link["observations"]:
            oid = diagnostic["observation_id"]
            if oid not in extension["observations"]:
                raise ValueError("competing extent link observation is missing")
            candidate = extension["observations"][oid]
            relation = exact_terminal_relation(anchor["repeat_intervals"], candidate["repeat_intervals"])
            if (candidate["source_geometry_id"] != gid or relation is None
                    or not _same(anchor["legacy_locus_key"], candidate["legacy_locus_key"])):
                raise ValueError("competing extent link relation differs")
            if not _same(relation, link["relation"]):
                raise ValueError("competing extent exact mapping differs")
            expected_observations.append(_diagnostics(aid, anchor, candidate, relation, accepted))
            used.add(oid)
        expected_observations.sort(key=lambda d: d["observation_id"])
        if (len({d["observation_id"] for d in expected_observations}) != len(expected_observations)
                or not _same(link, {"anchor_array_id": aid, "source_geometry_id": gid,
                    "relation": link["relation"], "observations": expected_observations,
                    "link_id": _identity([RULE, extension["source"]["sha256"], aid, gid])})):
            raise ValueError("competing extent diagnostics differ")
        seen.add((aid, gid))
        anchors.add(aid)
        gids.add(gid)
    counts = {"accepted_anchors": len(anchors), "competing_geometries": len(gids),
              "candidate_observations": len(used), "links": len(seen),
              "additional_accepted_arrays": 0, "additional_accepted_repeats": 0}
    if (used != set(extension["observations"]) or anchors != set(extension["anchors"])
            or not _same(counts, extension["counts"])):
        raise ValueError("competing extent registry/counts differ")
    return extension


def competing_terminal_extent_summaries(report):
    """Validated per-observation views; sequences reference the one registry."""
    extension = validate_competing_terminal_extents(report)
    result = {}
    if extension is None:
        return result
    for link in extension["links"]:
        aid = link["anchor_array_id"]
        anchor = extension["anchors"][aid]["observation"]
        for diagnostic in link["observations"]:
            record = extension["observations"][diagnostic["observation_id"]]
            result.setdefault(aid, []).append({
                "accepted_repeat_count": anchor["repeat_count"], "accepted_own_score": anchor["score"],
                "observation_id": record["observation_id"], "geometry_id": record["source_geometry_id"],
                "candidate_id": record["candidate_id"], "own_score": record["score"],
                "category": record["category"], "repeat_count": record["repeat_count"],
                "status": diagnostic["status"], "block_reasons": diagnostic["block_reasons"],
                "negative_review": record["negative_evidence"], "source_intervals": record["repeat_intervals"],
                "repeats": record["repeat_sequences"], "spacers": record["spacer_sequences"],
                "relation": link["relation"], "link_id": link["link_id"],
                "report_pointer": "report.json#/extensions/competing_terminal_extents/observations/" + record["observation_id"]})
    return result
