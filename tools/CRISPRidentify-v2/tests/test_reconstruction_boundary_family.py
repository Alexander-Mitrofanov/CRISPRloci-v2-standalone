"""Coupled boundary-family selection, proof and export fixture."""
import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from components.components_detection_refinement import CrisprCandidate
from tests.test_reconstruction_transaction import fixture as exact_fixture
from tests.test_candidate_evidence import PARAMETERS
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.reconstruction_selection import CandidateValidation, select_reconstructions, RowLocation
from components.reconstruction_transaction import (
    reconstruct_accepted_arrays, validate_boundary_family_basis, _snapshot,
)
from components.reconstruction_boundary_family import (
    METHOD, boundary_family_summary, describe_family, digest,
    report_boundary_family_summary, validate_boundary_family_audit, family_occurrence_evidence,
)
from components.report_model import ReportModelBuilder


def fixture(*, end_trim=1):
    """Source-valid coupled ends and model-real synthetic vectors; no predictor.

    Six copies; only two whole donor geometries are supplied. The alternate
    shortens ends3/4 together. Mixed end choices are deliberately absent.
    """
    source, categories, scorer = exact_fixture()
    key = next(iter(categories[0]))
    donor = categories[1][key][0][1]
    starts = list(donor.list_repeat_starts)
    repeats = [r[:-end_trim] if i in (3, 4) else r for i, r in enumerate(donor.list_repeats)]
    gapped = [r + "-" * end_trim if i in (3, 4) else r for i, r in enumerate(repeats)]
    alternative = CrisprCandidate(
        list_repeat_starts=starts, list_repeats=repeats,
        list_repeats_gaped=gapped,
        list_spacers=[source[a + len(r):b] for a, r, b in zip(starts, repeats, starts[1:])])
    categories[1][key].append(scorer({key: [alternative]})[key][0])
    # Keep the fixture's own model-real equal donor scores. Its feature values
    # are synthetic, not independent predictor or biological evidence.
    assert np.array_equal(categories[1][key][0][2][0], categories[1][key][1][2][0])
    return source, categories, scorer


def qualified(row):
    return CandidateValidation(row, True, True, True, True, {
        "likely_tandem": {"status": "passed"},
        "coding_tandem": {"status": "inapplicable"},
        "known_negative_annotations": {"status": "passed"},
    })


class BoundaryFamilySelectionTests(unittest.TestCase):
    def test_coupled_geometry_aliases_and_permutations_never_vote(self):
        _source, categories, _score = fixture()
        key = next(iter(categories[1]))
        first = categories[1][key][0]
        baseline = select_reconstructions(categories, qualified).decisions[0]
        self.assertEqual(len(baseline.family_members), 2)
        for rows in (list(reversed(categories[1][key])), categories[1][key] + [first]):
            altered = list(categories)
            altered[1] = {key: rows}
            result = select_reconstructions(altered, qualified).decisions[0]
            self.assertEqual(result.tied_geometry_ids, baseline.tied_geometry_ids)
            self.assertEqual(result.selected.geometry_id, baseline.selected.geometry_id)
            self.assertEqual(len(result.family_members), 2)

    def test_same_rank_veto_does_not_manufacture_unique_survivor(self):
        _source, categories, _score = fixture()
        bad = next(iter(categories[1].values()))[1]
        for failed in (False, None):
            def validate(row):
                v = qualified(row)
                return replace(v, negatives_passed=failed) if row is bad else v
            result = select_reconstructions(categories, validate).decisions[0]
            self.assertEqual(result.action, "abstained")
            self.assertEqual(result.reason, "same_rank_peer_veto")
        def lower_failed(row):
            if row is bad:
                return replace(qualified([.61, row[1], row[2]]), negatives_passed=False)
            return qualified(row)
        result = select_reconstructions(categories, lower_failed).decisions[0]
        self.assertEqual(result.action, "selected")
        self.assertFalse(result.family_members)

    def test_negative_applicability_disagreement_abstains(self):
        _source, categories, _score = fixture()
        changed = next(iter(categories[1].values()))[1]
        def validate(row):
            v = qualified(row)
            if row is changed:
                v.evidence["coding_tandem"]["status"] = "passed"
            return v
        result = select_reconstructions(categories, validate).decisions[0]
        self.assertEqual(result.reason, "boundary_family_negative_disagreement")

    def test_existing_score_tolerance_retains_whole_family_only_within_limit(self):
        _source, categories, _score = fixture()
        changed = next(iter(categories[1].values()))[1]
        for difference, expected in ((.5e-12, 2), (2e-12, 0)):
            def validate(row):
                if row is changed:
                    return qualified([row[0] - difference, row[1], row[2]])
                return qualified(row)
            selected = select_reconstructions(categories, validate).decisions[0]
            self.assertEqual(selected.action, "selected")
            self.assertEqual(len(selected.family_members), expected)

    def test_different_start_tie_cannot_be_dropped_as_a_subgroup(self):
        from tests.test_reconstruction_selection import categories as make_categories, row
        a = row([60,100,140,180], .7)
        b = row([60,100,140,180], .7, repeats=["ACGTTGCAA"] * 4)
        c = row([100,140,180,220], .7)
        result = select_reconstructions(make_categories(donors=[a,b,c]), qualified).decisions[0]
        self.assertEqual(result.reason, "ambiguous_best_geometry")
        self.assertEqual(len(result.tied_geometry_ids), 3)

    def test_non_display_member_collision_still_blocks(self):
        from tests.test_reconstruction_selection import categories as make_categories, row, REPEAT
        a = row([60,100,140,180,220], .7)
        b = row([60,100,140,180,220], .7, repeats=[REPEAT] * 4 + [REPEAT + "A"])
        inputs = make_categories(donors=[a,b], possible=[row([230,270,310], .65)])
        decision = select_reconstructions(inputs, qualified).decisions[0]
        self.assertEqual(decision.reason, "same_rank_peer_veto")
        self.assertIn("overlaps_other_accepted_core", decision.evidence["same_rank_peers"][1]["vetoes"])

    def test_family_union_proposal_collision_is_symmetric(self):
        from tests.test_reconstruction_selection import row, REPEAT
        first, second = row([100,140,180], .9), row([500,540,580], .9)
        short = row([100,140,180,220,260,300,340], .7000000000005)
        long = row([100,140,180,220,260,300,340], .7,
                   repeats=[REPEAT]*6 + [REPEAT+'A'])
        other = row([350,390,430,470,500,540,580], .7)
        alone = select_reconstructions(
            [{"core":[first]}, {"donors":[short,long]}, {},{},{}], qualified).decisions[0]
        self.assertIs(alone.selected.row, short)
        self.assertEqual(alone.selected.span[1], 350)
        self.assertEqual(alone.collision_span[1], 351)
        for anchors in ([first,second], [second,first]):
            result = select_reconstructions(
                [{"cores":anchors}, {"donors":[short,long,other]}, {},{},{}], qualified)
            self.assertEqual(result.stats["selected"], 0)
            self.assertEqual([d.reason for d in result.decisions], ["proposal_collision"]*2)
            self.assertTrue(all(not d.family_members for d in result.decisions))

    def test_s1_witness_on_non_display_member_is_not_erased(self):
        from tests.test_initial_competition import row, evaluate, members
        ordinary = row(.8, [(10,16),(40,46),(70,76)])
        repaired = row(.9, [(13,19),(43,49),(73,79)], repair=True)
        inputs = evaluate([ordinary, repaired])
        core = members(inputs, 0)[0]
        a = row(.7, [(13,19),(43,49),(73,79),(103,109)])
        b = row(.7, [(13,20),(43,49),(73,79),(103,109)])
        for r in (core,a,b):
            r[1].consensus = "AAAAAA"
        inputs[1].setdefault((0,100000000), []).extend([a,b])
        decision = select_reconstructions(inputs, qualified).decisions[0]
        self.assertEqual(decision.reason, "same_rank_peer_veto")
        self.assertTrue(any("protected_ordinary_occurrence_conflict" in p["vetoes"]
                            for p in decision.evidence["same_rank_peers"]))

    def test_member_specific_witness_applies_to_every_family_geometry(self):
        from tests.test_initial_competition import row
        from components.initial_competition import initial_competition_groups
        ordinary = row(.8, [(10,16),(40,46),(70,76)])
        a = row(.7, [(13,19),(43,49),(73,79),(103,109)], repair=True)
        groups, conflicts = initial_competition_groups([ordinary,a])
        self.assertFalse(conflicts)
        a = next(r for r in groups[0].rows if r[1] is a[1])
        b = row(.7, [(13,20),(43,49),(73,79),(103,109)])
        core = row(.9, [(13,19),(43,49),(73,79)])
        for r in (core,a,b):
            r[1].consensus = "AAAAAA"
        decision = select_reconstructions(
            [{"core":[core]}, {"donors":[a,b,ordinary]}, {},{},{}], qualified).decisions[0]
        self.assertEqual(decision.reason, "boundary_family_occurrence_conflict")
        evidence = decision.evidence["family_occurrence_guard"]
        self.assertEqual(evidence["direct_checks"], 2)
        self.assertEqual(len(evidence["conflicts"]), 1)


class BoundaryFamilyTransactionTests(unittest.TestCase):
    def output(self):
        source, categories, scorer = fixture()
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        row = next(iter(output[0].values()))[0]
        return source, categories, output, audit, row

    def test_real_source_model_proof_and_all_original_members_survive(self):
        source, before, output, audit, row = self.output()
        basis = row[3]["classification_basis"]
        self.assertEqual(basis["method"], METHOD)
        self.assertEqual(audit["stats"]["selected"], 1)
        self.assertEqual(audit["stats"]["boundary_uncertain_selected"], 1)
        self.assertEqual(audit["stats"]["added_repeats"], 3)
        original = [r for c in before[:2] for rows in c.values() for r in rows]
        retained = [r for rows in output[1].values() for r in rows]
        self.assertTrue(all(any(r is old for r in retained) for old in original))
        core = next(iter(before[0].values()))[0]
        self.assertNotIn("classification_basis", core[3] if len(core) > 3 else {})
        self.assertEqual(basis, validate_boundary_family_basis(source, row[1], basis, row[0], 0, row, output))
        summary = boundary_family_summary(row)
        self.assertEqual(summary["variable_repeat_indexes"], [3,4])
        self.assertEqual(summary["display_own_score"], row[0])
        self.assertNotEqual(row[0], original[0][0])

    def test_alias_observations_json_replay_and_tandem_wrapper_summary(self):
        source, categories, scorer = fixture()
        key = next(iter(categories[1]))
        categories[1][key].append(copy.deepcopy(categories[1][key][0]))
        output, _ = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        row = next(iter(output[0].values()))[0]
        basis = row[3]["classification_basis"]
        proof = basis["boundary_family"]
        self.assertEqual(len(proof["members"]), 2)
        self.assertEqual(len(proof["additional_peer_representations"]), 1)
        restored = json.loads(json.dumps(basis))
        self.assertEqual(restored, validate_boundary_family_basis(
            source, row[1], restored, row[0], 0, row, output))
        wrapped = copy.deepcopy(row)
        wrapped[3]["classification_basis"]["accepted_family_tandem_review"] = {"supported": False}
        # Full tandem-wrapper source replay belongs to the report validator;
        # compact summaries must retain the family annotation through it.
        self.assertEqual(boundary_family_summary(wrapped), boundary_family_summary(row))

    def test_member_s1_witness_is_bound_after_accepted_metadata_is_removed(self):
        from components.candidate_origins import mark_repair_origin
        from components.initial_competition import initial_competition_groups
        source, categories, scorer = fixture()
        key = next(iter(categories[1]))
        core = categories[0][key][0]
        donor = categories[1][key][1]
        mark_repair_origin(donor[1], {
            "policy": "terminal-only-trim-v1", "status": "generated",
            "coordinate_system": "zero_based_half_open",
            "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry": [donor[1].list_repeat_starts,
                [a + len(r) for a, r in zip(donor[1].list_repeat_starts, donor[1].list_repeats)]],
        })
        groups, conflicts = initial_competition_groups([core, donor])
        self.assertFalse(conflicts)
        categories[1][key][1] = next(r for r in groups[0].rows if r[1] is donor[1])
        output, _ = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        row = output[0][key][0]
        proof = row[3]["classification_basis"]["boundary_family"]
        member = next(m for m in proof["members"] if m["initial_competition"] is not None)
        self.assertIsNone(proof["anchor_initial_competition"])
        self.assertNotIn("initial_competition", row[3])
        self.assertEqual(row[3]["classification_basis"], validate_boundary_family_basis(
            source, row[1], row[3]["classification_basis"], row[0], 0, row, output))
        for mutation in ("source_span_sha256", "representation_sha256", "score"):
            changed = copy.deepcopy(output)
            current = changed[0][key][0]
            basis = current[3]["classification_basis"]
            peer = next(m for m in basis["boundary_family"]["members"]
                        if m["location"] == member["location"])
            ref = peer["initial_competition"]["occurrence_protection"]["ordinary_reference"]
            ref[mutation] = .51 if mutation == "score" else "0"*64
            basis["evidence_sha256"] = digest(basis["boundary_family"])
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                validate_boundary_family_basis(source, current[1], basis, current[0], 0, current, changed)

    def test_source_valid_member_specific_witness_blocks_producer_and_replay(self):
        from components.candidate_origins import mark_repair_origin, candidate_origin_evidence
        from components.initial_competition import initial_competition_groups
        from components.reconstruction_validator import validate_candidate_reconstruction
        source, categories, scorer = fixture(end_trim=2)
        output, _ = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        key = next(iter(categories[1]))
        donor = categories[1][key][1]
        # A source-valid four-copy ordinary witness maps uniquely to the short
        # member. Its fourth interval intersects that 38-base repeat by19;
        # the same19 bases cover less than half of the other40-base member.
        starts = donor[1].list_repeat_starts
        intervals = [(a,a+len(r)) for a,r in zip(starts,donor[1].list_repeats)][:3]
        intervals.append((starts[3]-19, starts[3]+19))
        repeats = [source[a:b] for a,b in intervals]
        width = max(map(len,repeats))
        candidate = CrisprCandidate(list_repeat_starts=[a for a,b in intervals],
            list_repeats=repeats, list_repeats_gaped=[r+'-'*(width-len(r)) for r in repeats],
            list_spacers=[source[a[1]:b[0]] for a,b in zip(intervals,intervals[1:])])
        witness = scorer({key:[candidate]})[key][0]
        self.assertTrue(validate_candidate_reconstruction(source,candidate).valid)
        mark_repair_origin(donor[1], {
            "policy":"terminal-only-trim-v1", "status":"generated",
            "coordinate_system":"zero_based_half_open",
            "source_sha256":hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry":[starts,[a+len(r) for a,r in zip(starts,donor[1].list_repeats)]],
        })
        groups, conflicts = initial_competition_groups([witness,donor])
        self.assertFalse(conflicts)
        annotated = next(r for r in groups[0].rows if r[1] is donor[1])
        categories[1][key][1] = annotated
        categories[1][key].append(witness)
        result, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selected"], 0)
        self.assertEqual(audit["decisions"][0]["reason"], "boundary_family_occurrence_conflict")
        self.assertIs(result[0][key][0], categories[0][key][0])
        # Coherently attach this real retained historical observation to an
        # otherwise valid pre-witness family. Per-member self-checks pass;
        # publication must still reject the other geometry directly.
        current = output[0][key][0]
        basis = current[3]["classification_basis"]
        proof = basis["boundary_family"]
        member = next(m for m in proof["members"] if m["original"] == _snapshot(donor))
        member["initial_competition"] = annotated[3]["initial_competition"]
        member["origins"] = candidate_origin_evidence(donor[1], proof["source_sha256"])
        output[1][key] = [annotated if r[1] is donor[1] else r for r in output[1][key]] + [witness]
        basis["evidence_sha256"] = digest(proof)
        with self.assertRaisesRegex(ValueError, "union of original"):
            validate_boundary_family_basis(source,current[1],basis,current[0],0,current,output)

    def test_retention_off_canonical_report_contains_all_coupled_members(self):
        source, _before, output, _audit, row = self.output()
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/"source.fa"
            path.write_text(">source\n" + source + "\n")
            report = ReportModelBuilder(path, [output[0],{},{},{},{}], {}, [COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS, param_accepted_only=True), provenance_categories=output).build()
        array = report["arrays"][0]
        validate_boundary_family_audit(report, _audit)
        self.assertEqual(report_boundary_family_summary(array), boundary_family_summary(row))
        self.assertEqual(len(array["classification_basis"]["boundary_family"]["members"]), 2)
        dimensions = array["boundary_assessment"]["dimensions"]
        self.assertEqual(dimensions["repeat_phase"]["status"], "ambiguous")
        self.assertEqual(dimensions["occurrence_path"]["status"], "unassessed")
        for field in ("display_own_score", "family_id"):
            bad = copy.deepcopy(array)
            bad["boundary_family"][field] = 0 if field == "display_own_score" else "0"*64
            with self.assertRaises(ValueError):
                report_boundary_family_summary(bad)
        for audit_mutation in ("missing", "source", "member", "counts", "duplicate"):
            audit = copy.deepcopy(_audit)
            if audit_mutation == "missing":
                audit = None
            elif audit_mutation == "source":
                audit["source_sha256"] = "0" * 64
            elif audit_mutation == "member":
                audit["decisions"][0]["tied_geometry_ids"].pop()
            elif audit_mutation == "duplicate":
                audit["decisions"].append(copy.deepcopy(audit["decisions"][0]))
            else:
                audit["stats"]["boundary_uncertain_selected"] = 0
            with self.subTest(mutation=audit_mutation), self.assertRaises(ValueError):
                validate_boundary_family_audit(report, audit)

    def test_member_snapshot_and_membership_mutations_fail_with_new_digest(self):
        source, _before, output, _audit, row = self.output()
        for mutation in ("score", "origins", "initial", "omit", "core", "source"):
            bad_output = copy.deepcopy(output)
            bad = next(iter(bad_output[0].values()))[0]
            basis = bad[3]["classification_basis"]
            proof = basis["boundary_family"]
            member = proof["members"][0]
            if mutation == "score":
                member["fresh"]["score"] += .01
            elif mutation == "origins":
                member["origins"]["ordinary_birth_sha256"] = "0"*64
            elif mutation == "initial":
                member["initial_competition"] = {"rule":"bogus"}
            elif mutation == "omit":
                proof["members"].pop()
            elif mutation == "core":
                proof["anchor_original"]["repeat_intervals"][0][1] -= 1
            else:
                proof["source_sha256"] = "0"*64
            basis["evidence_sha256"] = digest(proof)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                validate_boundary_family_basis(source, bad[1], basis, bad[0], 0, bad, bad_output)

    def test_actual_retained_negative_and_other_accepted_collision_are_replayed(self):
        source, _before, output, _audit, _row = self.output()
        for corruption in ("negative", "collision", "source"):
            changed = copy.deepcopy(output)
            current = next(iter(changed[0].values()))[0]
            basis = current[3]["classification_basis"]
            actual_source = source
            if corruption == "negative":
                donor = next(iter(changed[1].values()))[0]
                if len(donor) == 3:
                    donor.append({})
                donor[3]["tags"] = ["tandem"]
            elif corruption == "collision":
                original_core = next(iter(changed[1].values()))[-1]
                changed[2]["other"] = [original_core]
            else:
                actual_source = "T" + source[1:]
            with self.subTest(corruption=corruption), self.assertRaises(ValueError):
                validate_boundary_family_basis(actual_source, current[1], basis, current[0], 0, current, changed)

    def test_source_valid_unadmitted_cartesian_hybrid_fails_membership(self):
        source, _before, output, _audit, _row = self.output()
        bad = copy.deepcopy(output)
        current = next(iter(bad[0].values()))[0]
        basis = current[3]["classification_basis"]
        proof = basis["boundary_family"]
        # Replace only a non-display member by an actual-source hybrid. Both
        # end values appear somewhere in the family; their combination never
        # existed as a retained observation. Recompute superficial IDs/digests.
        member = next(m for m in proof["members"]
                      if m["location"]["representation_id"] != proof["display_representation_id"])
        old_location = copy.deepcopy(member["location"])
        intervals = copy.deepcopy(member["original"]["repeat_intervals"])
        # Derive a mixed choice from the two actual member end sets, rather
        # than assume the deterministic display ordering.
        ends = proof["descriptor"]["observed_end_sets"]
        intervals[3][1] = ends[3][0]
        intervals[4][1] = ends[4][-1]
        repeats = [source[a:b] for a,b in intervals]
        width = max(map(len, repeats))
        c = CrisprCandidate(list_repeat_starts=[a for a,b in intervals],
            list_repeats=repeats, list_repeats_gaped=[r+'-'*(width-len(r)) for r in repeats],
            list_spacers=[source[a[1]:b[0]] for a,b in zip(intervals,intervals[1:])])
        vector = np.array([[member["original"]["features"][n] for n in COMPACT_FEATURE_SCHEMA]])
        hybrid = [member["original"]["score"], c, [vector]]
        member["original"] = member["fresh"] = _snapshot(hybrid)
        newloc = RowLocation(1, old_location["locus_key"], old_location["row_index"], hybrid).as_dict()
        member["location"] = newloc
        from components.candidate_origins import candidate_origin_evidence
        member["origins"] = candidate_origin_evidence(c, proof["source_sha256"])
        for peer in proof["selection"]["evidence"]["same_rank_peers"]:
            if peer["candidate"] == old_location:
                peer["candidate"] = newloc
        proof["members"].sort(key=lambda m:m["location"]["geometry_id"])
        proof["selection"]["tied_geometry_ids"] = sorted(m["location"]["geometry_id"] for m in proof["members"])
        proof["descriptor"] = describe_family(proof["source_sha256"], proof["anchor_original"],
                                               proof["members"], proof["display_representation_id"])
        proof["selection"]["evidence"]["family_occurrence_guard"] = family_occurrence_evidence(
            [], [member["fresh"]["repeat_intervals"] for member in proof["members"]])
        basis["evidence_sha256"] = digest(proof)
        with self.assertRaisesRegex(ValueError, "retained original"):
            validate_boundary_family_basis(source, current[1], basis, current[0], 0, current, bad)


if __name__ == "__main__":
    unittest.main()
