"""Atomic, freshly scored reconstruction of an already accepted array core.

Selection and candidate checks are separate pure modules. This adapter owns
batch scoring, source-CDS acquisition, category transitions and replayable
report evidence. It never changes candidate geometry or copies a core score.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import asdict
import hashlib
import json
import math

import numpy as np

from components.compact_model import predict_compact_scores
from components.components_detection_refinement import CrisprCandidate
from components.fast_run_rescue import candidate_score_fingerprint, compact_feature_values
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.initial_competition import (
    METADATA_KEY as COMPETITION_METADATA_KEY,
    competition_evidence,
    guard_reconstructed_occurrences,
    validate_ordinary_witness_inventory,
    ordinary_reference_key,
    protected_ordinary_reference,
    with_competition_evidence,
)
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.reconstruction_selection import (
    ReconstructionPolicy,
    candidate_geometry_id,
    reconstruction_candidate_locations,
    select_reconstructions,
)
from components.reconstruction_validation import (
    _known_negative_annotations,
    needs_coding_context,
    validate_reconstruction_candidate,
)
from components.architecture_rescue_coding_context import SourceCDS
from components.candidate_origins import candidate_origin_evidence
from components.reconstruction_boundary_family import (
    METHOD as FAMILY_METHOD, describe_family, boundary_family_summary, family_occurrence_evidence,
)
from components.reconstruction_selection import relationship
from components.reconstruction_reference import MAPPED_CORE_REFERENCE, recorded_reference_rule


METHOD = "source_backed_candidate_reconstruction"
EVIDENCE_FILENAME = "reconstruction_evidence.json"
STRUCTURAL_PARAMETERS = (
    "param_min_avg_repeat_length", "param_max_avg_repeat_length",
    "param_min_avg_spacer_length", "param_max_avg_spacer_length",
    "param_min_repeats", "param_max_identical_spacers",
    "param_max_identical_cluster_spacers",
)


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    if hasattr(value, "tolist"):
        return _plain(value.tolist())
    return value


def _encoded(value):
    return json.dumps(_plain(value), sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def _digest(value):
    return hashlib.sha256(_encoded(value)).hexdigest()


def _snapshot(row):
    candidate = row[1]
    return {
        "score": float(row[0]),
        "repeat_intervals": [[int(start), int(start) + len(repeat)]
                             for start, repeat in zip(candidate.list_repeat_starts,
                                                      candidate.list_repeats)],
        "gapped_repeats": list(candidate.list_gaped_repeats),
        "features": compact_feature_values(row),
        "eligible": getattr(candidate, "source_path_classification_eligible", True),
        "eligibility_reason": getattr(candidate, "source_path_classification_reason", None),
    }


def _inventory_state(categories):
    return tuple((index, repr(key), ordinal, id(row), float(row[0]),
                  candidate_score_fingerprint(row[1]), _digest(row[2:]))
                 for index, category in enumerate(categories)
                 for key, rows in category.items() for ordinal, row in enumerate(rows))


def reconstruct_accepted_arrays(source, categories, parameters, score_candidates,
                                get_source_cds, policy=None):
    """Return a new category graph and audit; leave all input rows untouched."""
    policy = policy or ReconstructionPolicy()
    structural_parameters = {k: parameters[k] for k in STRUCTURAL_PARAMETERS}
    structural = build_advanced_fuzzy_search_filter(structural_parameters)
    validate_ordinary_witness_inventory(categories, source)
    before = _inventory_state(categories)
    # One explicit recipe for both admission to fresh scoring and selection.
    reference_rule = MAPPED_CORE_REFERENCE
    donors = reconstruction_candidate_locations(categories, policy, reference_rule=reference_rule)
    rows = [location.row for location in donors]
    if donors:
        rows.extend(row for values in categories[0].values() for row in values)
    candidates = list({id(row[1]): row[1] for row in rows}.values())
    key = "retained-reconstruction-fresh-score-batch"
    scored = score_candidates({key: candidates}) if candidates else {}
    if _inventory_state(categories) != before:
        raise ValueError("reconstruction scoring mutated the input inventory")
    if candidates and (set(scored) != {key} or len(scored[key]) != len(candidates)):
        raise ValueError("reconstruction scoring changed request cardinality")
    fresh = {}
    for candidate, row in zip(candidates, scored.get(key, ())):
        if row[1] is not candidate:
            raise ValueError("reconstruction scoring changed candidate identity/order")
        if not math.isfinite(float(row[0])) or not 0 <= row[0] <= 1:
            raise ValueError("reconstruction scoring returned an invalid score")
        compact_feature_values(row)
        fresh[id(candidate)] = row
    source_cds = (tuple(get_source_cds())
                  if any(needs_coding_context(c) for c in candidates) else None)
    validated = {}

    def validate(row):
        if id(row) not in validated:
            candidate = row[1]
            if id(candidate) not in fresh:
                raise ValueError("selection requested an unscored reconstruction")
            validated[id(row)] = validate_reconstruction_candidate(
                source, row, fresh[id(candidate)], structural, source_cds=source_cds,
            )
        return validated[id(row)]

    result = select_reconstructions(categories, validate, policy, reference_rule=reference_rule)
    if _inventory_state(categories) != before:
        raise ValueError("reconstruction validation mutated the input inventory")
    output = [{key: list(rows) for key, rows in category.items()}
              for category in categories]
    selected = [decision for decision in result.decisions if decision.action == "selected"]
    # Only the chosen observation moves into the accepted slot. Other
    # alignments of its geometry retain their own provenance as Alternatives;
    # their diagnostic alias membership does not make them extra selections.
    # Locations, rather than object IDs, also preserve a repeated observation
    # when the same row object occurs in more than one input position.
    remove_locations = {(d.selected.locus_key, d.selected.row_index)
                        for d in selected if not d.family_members}
    for key, rows in list(output[1].items()):
        remaining = [row for index, row in enumerate(rows)
                     if (key, index) not in remove_locations]
        if remaining:
            output[1][key] = remaining
        else:
            del output[1][key]
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    for decision in selected:
        old = decision.anchor
        original = decision.selected.row
        fresh_row = decision.validated_row
        core_validation = validated[id(old.row)]
        proposed_validation = validated[id(original)]
        if decision.family_members:
            def member_record(location, relation, validation):
                return _plain({"location": location.as_dict(),
                    "original": _snapshot(location.row), "fresh": _snapshot(validation.row),
                    "initial_competition": competition_evidence(location.row, len(source), source=source),
                    "origins": candidate_origin_evidence(location.row[1], source_hash),
                    "relationship": relation, "validation": validation.as_dict()})
            members = [member_record(*item) for item in decision.family_members]
            member_locations = {_encoded(item["location"]) for item in members}
            peer_locations = {_encoded(item["candidate"]) for item in decision.evidence["same_rank_peers"]}
            additional = []
            for location in donors:
                encoded_location = _encoded(location.as_dict())
                if encoded_location in peer_locations and encoded_location not in member_locations:
                    relation = relationship(old.row[1], location.row[1], policy,
                                            reference_rule=reference_rule)
                    relation = {k: v for k, v in relation.items() if k not in ("qualified", "reason")}
                    additional.append(member_record(location, relation, validated[id(location.row)]))
            proof = _plain({
                "source_sha256": source_hash, "coordinate_system": "zero_based_half_open",
                "policy": asdict(policy), "structural_parameters": structural_parameters,
                "anchor_original": _snapshot(old.row), "anchor_fresh": _snapshot(core_validation.row),
                "anchor_initial_competition": competition_evidence(old.row, len(source), source=source),
                "anchor_validation": core_validation.evidence,
                "members": members, "additional_peer_representations": additional,
                "display_representation_id": decision.selected.representation_id,
                "selection": decision.as_dict(),
                "score_semantics": "own_fresh_member_model_scores; core_confidence_is_separate",
            })
            proof["descriptor"] = describe_family(source_hash, proof["anchor_original"],
                                                   members, decision.selected.representation_id)
            metadata = deepcopy(original[3]) if len(original) == 4 else {}
            metadata.pop(COMPETITION_METADATA_KEY, None)
            # Any prior boundary assessment was bound to a different ranking
            # decision. The report derives ambiguous ends from this full proof.
            metadata.pop("boundary_assessment", None)
            metadata["classification_basis"] = {"method": FAMILY_METHOD,
                "boundary_family": proof, "evidence_sha256": _digest(proof)}
            output[0][old.locus_key][old.row_index] = list(fresh_row[:3]) + [metadata]
            output[1].setdefault(old.locus_key, []).append(old.row)
            continue
        # Evidence dataclasses can contain tuples. Canonicalize before the
        # in-memory publication replay so it has the same representation as
        # a JSON-loaded report; otherwise an applicable CDS check can differ
        # solely because its intervals are tuples rather than lists.
        proof = _plain({
            "source_sha256": source_hash, "coordinate_system": "zero_based_half_open",
            "policy": asdict(policy), "structural_parameters": structural_parameters,
            "anchor_original": _snapshot(old.row),
            "anchor_fresh": _snapshot(core_validation.row),
            # Historical observations remain bound to the original score and
            # geometry. A fresh score or a larger donor does not rewrite an
            # earlier tie, and the observation is not selection authority.
            "anchor_initial_competition": competition_evidence(old.row, len(source), source=source),
            "donor_original": _snapshot(original),
            "donor_initial_competition": competition_evidence(original, len(source), source=source),
            "reconstruction": _snapshot(fresh_row),
            "anchor_validation": core_validation.evidence,
            "reconstruction_validation": proposed_validation.evidence,
            "selection": decision.as_dict(),
            "score_semantics": "own_fresh_model_score; core_confidence_is_separate",
        })
        metadata = deepcopy(original[3]) if len(original) == 4 else {}
        metadata.pop(COMPETITION_METADATA_KEY, None)
        metadata["classification_basis"] = {
            "method": METHOD, "reconstruction": proof,
            "evidence_sha256": _digest(proof),
        }
        replacement = list(fresh_row[:3]) + [metadata]
        output[0][old.locus_key][old.row_index] = replacement
        output[1].setdefault(old.locus_key, []).append(old.row)
    assert sum(len(v) for v in output[0].values()) == sum(len(v) for v in categories[0].values())
    audit = {
        "schema": {"name": "CRISPRidentify-reconstruction-evidence", "version": "1.0.0"},
        "source_sha256": source_hash,
        "coordinate_system": "zero_based_half_open", "policy": asdict(policy),
        "experimental": False, "biological_accuracy_established": False,
        "stats": dict(result.stats, fresh_candidates_scored=len(candidates),
                      candidates_validated=len(validated), selections_committed=len(selected)),
        "actions": dict(Counter(d.action for d in result.decisions)),
        "decisions": [d.as_dict() for d in result.decisions],
    }
    return output, _plain(audit)


def _restore(source, snapshot):
    intervals = snapshot["repeat_intervals"]
    starts = [a for a, _ in intervals]
    candidate = CrisprCandidate(
        list_repeats=[source[a:b] for a, b in intervals],
        list_repeats_gaped=snapshot["gapped_repeats"],
        list_spacers=[source[left[1]:right[0]] for left, right in zip(intervals, intervals[1:])],
        list_repeat_starts=starts,
    )
    candidate.source_path_classification_eligible = snapshot["eligible"]
    if snapshot["eligibility_reason"] is not None:
        candidate.source_path_classification_reason = snapshot["eligibility_reason"]
    row = [snapshot["score"], candidate,
           [np.array([[snapshot["features"][name] for name in COMPACT_FEATURE_SCHEMA]])]]
    if _snapshot(row) != snapshot:
        raise ValueError("reconstruction snapshot does not round-trip")
    score = float(predict_compact_scores(row[2][0])[0])
    if abs(score - row[0]) > 1e-12:
        raise ValueError("reconstruction score does not match its feature vector")
    return row


def validate_reconstruction_basis(source, candidate, basis, row_score,
                                  category_index, current_row, categories, *, ordinary_witness_keys=None,
                                  preserved_core_rows=None):
    """Replay source/score/negative evidence before report publication."""
    if set(basis) != {"method", "reconstruction", "evidence_sha256"}:
        raise ValueError("reconstruction classification basis keys differ")
    proof = basis["reconstruction"]
    if basis["evidence_sha256"] != _digest(proof):
        raise ValueError("reconstruction evidence digest differs")
    if proof["source_sha256"] != hashlib.sha256(source.encode("ascii")).hexdigest():
        raise ValueError("reconstruction source binding differs")
    if category_index != 0:
        raise ValueError("a core-confidence reconstruction must be Bona-fide")
    if current_row is None or current_row[1] is not candidate:
        raise ValueError("reconstruction report row is missing")
    if float(row_score) != proof["reconstruction"]["score"] or _snapshot(current_row) != proof["reconstruction"]:
        raise ValueError("reported reconstruction differs from the selected row")
    for before_key, after_key in (("anchor_original", "anchor_fresh"),
                                  ("donor_original", "reconstruction")):
        for field in ("repeat_intervals", "gapped_repeats", "eligible", "eligibility_reason"):
            if proof[before_key][field] != proof[after_key][field]:
                raise ValueError("fresh scoring changed reconstruction geometry/eligibility")
    if ordinary_witness_keys is None:
        ordinary_witness_keys = validate_ordinary_witness_inventory(categories, source,
            observations=(proof.get("anchor_initial_competition"), proof.get("donor_initial_competition")))
    policy = ReconstructionPolicy(**proof["policy"])
    if asdict(policy) != asdict(ReconstructionPolicy()):
        raise ValueError("unsupported reconstruction policy")
    anchor = _restore(source, proof["anchor_fresh"])
    proposal = _restore(source, proof["reconstruction"])
    matching = (preserved_core_rows if preserved_core_rows is not None else
                [row for values in categories[1].values() for row in values
                 if _snapshot(row) == proof["anchor_original"]])
    if any(_snapshot(row) != proof["anchor_original"] for row in matching):
        raise ValueError("preserved reconstruction core context differs")
    if not matching:
        raise ValueError("reconstruction has no preserved original core")
    historical_rows = []
    for name in ("anchor", "donor"):
        recorded = proof.get(name + "_initial_competition")
        if recorded is not None:
            historical = _restore(source, proof[name + "_original"])
            historical = with_competition_evidence(historical, recorded)
            competition_evidence(historical, len(source), source=source)
            reference = protected_ordinary_reference(historical, source=source)
            if reference is not None and ordinary_reference_key(reference) not in ordinary_witness_keys:
                raise ValueError("reconstruction has no preserved ordinary occurrence witness")
            historical_rows.append(historical)
    if not guard_reconstructed_occurrences(proof["reconstruction"]["repeat_intervals"],
                                          *historical_rows, source=source):
        raise ValueError("reconstruction violates protected ordinary occurrences")
    # Check the historical observation against the preserved original row.
    # The complete initial competitive inventory is not independently replayed.
    if not any(_plain(competition_evidence(row, len(source), source=source))
               == proof.get("anchor_initial_competition") for row in matching):
        raise ValueError("initial competition differs from the preserved original core")
    structural = build_advanced_fuzzy_search_filter(proof["structural_parameters"])
    validations = []
    for row, recorded in ((anchor, proof["anchor_validation"]),
                          (proposal, proof["reconstruction_validation"])):
        coding = recorded["coding_tandem"]
        cds = (tuple(SourceCDS(**raw) for raw in coding["considered_source_cds"])
               if coding["source_cds_assessed"] else None)
        checked = validate_reconstruction_candidate(source, row, row, structural, source_cds=cds)
        if not all((checked.source_valid, checked.structural_valid,
                    checked.negatives_passed, checked.score_fresh)):
            raise ValueError("reconstruction replay failed candidate-specific checks")
        if _plain(checked.evidence) != _plain(recorded):
            raise ValueError("reconstruction checks differ on source replay")
        validations.append(checked)
    # Replay local geometry/sequence eligibility and competing accepted spans.
    trial = [{(0, len(source)): [anchor]}, {(0, len(source)): [proposal]}, {}, {}, {}]
    for index in (0, 2):
        others = [row for rows in categories[index].values() for row in rows
                  if row[1] is not candidate]
        if others:
            trial[index].setdefault((0, len(source)), []).extend(others)
    lookups = {id(anchor): validations[0], id(proposal): validations[1]}
    def validate(row):
        if id(row) not in lookups:
            # Other accepted rows serve only as collision witnesses.
            raise ValueError("unexpected reconstruction replay validation")
        return lookups[id(row)]
    # Historical proofs without the fallback discriminator replay the old
    # anchor-only rule, independent of the current producer's default.
    reference_rule = recorded_reference_rule(proof["selection"]["relationship"])
    result = select_reconstructions(trial, validate, policy, reference_rule=reference_rule)
    core = next(d for d in result.decisions if d.anchor.row is anchor)
    if core.action != "selected":
        raise ValueError("reconstruction geometry/anti-fusion replay abstained")
    selection = proof["selection"]
    # New fallback evidence has an exact typed JSON contract. Python's ordinary
    # equality would accept True as index1 or an integral float as an index.
    # Retain the historical comparison semantics for legacy relationships.
    recorded_relation = selection["relationship"]
    relation_matches = (_encoded(recorded_relation) == _encoded(core.relation)
                        if "added_repeat_reference" in recorded_relation
                        else recorded_relation == core.relation)
    if selection["action"] != "selected" or not relation_matches:
        raise ValueError("recorded reconstruction relationship differs from source replay")
    for name, original, restored in (("anchor", proof["anchor_original"], anchor),
                                      ("selected", proof["donor_original"], proposal)):
        location = selection[name]
        if (location["repeat_intervals"] != original["repeat_intervals"]
                or location["geometry_id"] != candidate_geometry_id(restored[1])
                or location["score"] != original["score"]):
            raise ValueError("recorded reconstruction location differs from source replay")
    if selection["reconstruction_score"] != row_score:
        raise ValueError("recorded reconstruction score differs from final row")
    return basis


def validate_boundary_family_basis(source, candidate, basis, row_score,
                                   category_index, current_row, categories, *, ordinary_witness_keys=None):
    """Replay every coupled member against retained original observations.

    This checks the complete recorded same-rank peer set, not unrecorded
    upstream generation. Members cannot be substituted by source-valid hybrid
    ends: each original snapshot/representation/origin must actually exist.
    """
    if set(basis) != {"method", "boundary_family", "evidence_sha256"} or basis["method"] != FAMILY_METHOD:
        raise ValueError("boundary family classification basis differs")
    proof = basis["boundary_family"]
    if _digest(proof) != basis["evidence_sha256"]:
        raise ValueError("boundary family evidence digest differs")
    source_hash = hashlib.sha256(source.encode("ascii")).hexdigest()
    if proof["source_sha256"] != source_hash or proof["coordinate_system"] != "zero_based_half_open":
        raise ValueError("boundary family source binding differs")
    if category_index != 0 or current_row is None or current_row[1] is not candidate:
        raise ValueError("boundary family must have an actual Bona-fide view")
    if float(row_score) != float(current_row[0]):
        raise ValueError("boundary family report score differs")
    summary = boundary_family_summary(current_row)
    if summary is None or current_row[3]["classification_basis"].get("boundary_family") != proof:
        raise ValueError("boundary family current proof differs")
    policy = ReconstructionPolicy(**proof["policy"])
    if asdict(policy) != asdict(ReconstructionPolicy()):
        raise ValueError("unsupported boundary family policy")
    members = proof["members"]
    peer_records = members + proof["additional_peer_representations"]
    if len({_encoded(item["location"]) for item in peer_records}) != len(peer_records):
        raise ValueError("boundary family repeats a peer observation")
    observations = [proof.get("anchor_initial_competition")]
    observations.extend(item.get("initial_competition") for item in peer_records)
    if ordinary_witness_keys is None:
        ordinary_witness_keys = validate_ordinary_witness_inventory(categories, source, observations=observations)
    family_references = []
    originals = [(proof["anchor_original"], proof.get("anchor_initial_competition"))]
    originals.extend((item["original"], item.get("initial_competition")) for item in peer_records)
    for snapshot, observation in originals:
        if observation is not None:
            original = with_competition_evidence(_restore(source, snapshot), observation)
            family_references.append(protected_ordinary_reference(original, source=source))
    family_guard = family_occurrence_evidence(family_references,
        [member["fresh"]["repeat_intervals"] for member in members])
    if (family_guard["status"] != "preserved"
            or proof["selection"]["evidence"]["family_occurrence_guard"] != family_guard):
        raise ValueError("boundary family violates the union of original ordinary witnesses")
    retained = [r for values in categories[1].values() for r in values]
    # One local value-key pass, not one full retained-source scan per alias.
    # The index owns only existing row references during this family replay.
    required = {_digest(item["original"]) for item in peer_records}
    core_key = _digest(proof["anchor_original"])
    required.add(core_key)
    retained_index = {key: [] for key in required}
    for row in retained:
        key = _digest(_snapshot(row))
        if key in retained_index:
            retained_index[key].append(row)
    if any(_known_negative_annotations(row)["status"] != "passed"
           for row in retained_index[core_key]):
        raise ValueError("boundary family retained core has a recorded negative")
    selection = proof["selection"]
    if (selection["action"] != "selected" or not selection.get("boundary_uncertain")
            or selection["reason"] != "validated_shared_occurrence_extent"):
        raise ValueError("boundary family selection is not uncertain")
    peer_evidence = selection["evidence"]["same_rank_peers"]
    expected_peer_evidence = [{"candidate": item["location"],
        "fresh_score": item["fresh"]["score"], "vetoes": []} for item in peer_records]
    if sorted(map(_encoded, peer_evidence)) != sorted(map(_encoded, expected_peer_evidence)):
        raise ValueError("boundary family same-rank peer inventory differs")
    member_ids = proof["descriptor"]["admissible_geometry_ids"]
    peer_ids = sorted({item["location"]["geometry_id"] for item in peer_records})
    if member_ids != peer_ids or selection["tied_geometry_ids"] != peer_ids:
        raise ValueError("boundary family omitted or substituted a tied geometry")
    evidence = selection["evidence"]
    if (evidence["same_rank_peer_geometries"] != len(peer_ids)
            or evidence["direct_peer_validation_complete"] is not True
            or evidence["freshly_scored_peer_rows"] != evidence["direct_peer_rows"]
            or evidence["direct_peer_rows"] < len(peer_records)
            or evidence["direct_peer_geometries"] < len(peer_ids)):
        raise ValueError("boundary family peer completeness accounting differs")
    count = proof["descriptor"]["supported_repeat_count"]
    best_score = selection["evidence"]["best_reconstruction_score"]
    if (selection["evidence"]["best_repeat_count"] != count
            or best_score != max(item["fresh"]["score"] for item in peer_records)):
        raise ValueError("boundary family winning rank differs")
    from components.reconstruction_boundary_family import negative_signature
    signatures = set()
    for item in peer_records:
        location = item["location"]
        original = _restore(source, item["original"])
        proposed = _restore(source, item["fresh"])
        if (len(item["fresh"]["repeat_intervals"]) != count
                or abs(item["fresh"]["score"] - best_score) > policy.score_tie_tolerance):
            raise ValueError("boundary family peer is not equally ranked")
        found = [r for r in retained_index[_digest(item["original"])] if _snapshot(r) == item["original"]
                 and candidate_origin_evidence(r[1], source_hash) == item["origins"]
                 and _plain(competition_evidence(r, len(source), source=source)) == item.get("initial_competition")]
        if not found:
            raise ValueError("boundary family member has no retained original observation")
        if any(_known_negative_annotations(row)["status"] != "passed" for row in found):
            raise ValueError("boundary family retained member has a recorded negative")
        if (location["geometry_id"] != candidate_geometry_id(original[1])
                or location["score"] != original[0]):
            raise ValueError("boundary family original location differs")
        from components.reconstruction_selection import RowLocation
        expected_location = _plain(RowLocation(1, location["locus_key"], location["row_index"], original).as_dict())
        if expected_location != location:
            raise ValueError("boundary family original representation differs")
        validation = item["validation"]
        if not all(validation.get(k) is True for k in (
                "source_valid", "structural_valid", "negatives_passed", "score_fresh", "qualified")):
            raise ValueError("boundary family peer has an unresolved validation")
        signatures.add(negative_signature(validation["evidence"]))
        # Reuse the exact reconstruction source/model/negative/direct-map
        # validator for each actual member. Only the accepted view is replaced
        # temporarily; all other accepted collision witnesses stay intact.
        temporary = [{key: [proposed if r is current_row else r for r in rows]
                      for key, rows in categories[0].items()}] + list(categories[1:])
        member_selection = dict(selection, selected=location,
            relationship=item["relationship"], reconstruction_score=proposed[0])
        member_proof = {"source_sha256": source_hash, "coordinate_system": "zero_based_half_open",
            "policy": proof["policy"], "structural_parameters": proof["structural_parameters"],
            "anchor_original": proof["anchor_original"], "anchor_fresh": proof["anchor_fresh"],
            "anchor_initial_competition": proof.get("anchor_initial_competition"),
            "donor_original": item["original"], "donor_initial_competition": item.get("initial_competition"),
            "reconstruction": item["fresh"], "anchor_validation": proof["anchor_validation"],
            "reconstruction_validation": validation["evidence"], "selection": member_selection}
        member_basis = {"method": METHOD, "reconstruction": member_proof,
                        "evidence_sha256": _digest(member_proof)}
        validate_reconstruction_basis(source, proposed[1], member_basis, proposed[0], 0,
            proposed, temporary, ordinary_witness_keys=ordinary_witness_keys,
            preserved_core_rows=retained_index[core_key])
    if len(signatures) != 1:
        raise ValueError("boundary family negative applicability disagrees")
    for member in members:
        same_geometry = [item for item in peer_records
                         if item["location"]["geometry_id"] == member["location"]["geometry_id"]]
        best = min(same_geometry, key=lambda item: (-item["fresh"]["score"], item["location"]["representation_id"]))
        if member != best:
            raise ValueError("boundary family geometry alias selection differs")
    display = min(members, key=lambda item: (-item["fresh"]["score"],
                  item["location"]["geometry_id"], item["location"]["representation_id"]))
    if (display["location"]["representation_id"] != proof["display_representation_id"]
            or _snapshot(current_row) != display["fresh"]
            or selection["selected"] != display["location"]
            or selection["reconstruction_score"] != current_row[0]):
        raise ValueError("boundary family actual display member differs")
    # Every member span has the same start, so the conservative envelope equals
    # a longest member span. Check it explicitly rather than rely on the view.
    left, right = proof["descriptor"]["union_span"]
    for index in (0, 2):
        for rows in categories[index].values():
            for row in rows:
                if row is current_row:
                    continue
                intervals = _snapshot(row)["repeat_intervals"]
                if min(right, intervals[-1][1]) > max(left, intervals[0][0]):
                    raise ValueError("boundary family union span overlaps another accepted core")
    return basis
