"""Possible-tier extents preserve original confidence and full observations."""

import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from components.components_detection_refinement import CrisprCandidate
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.compact_model import predict_compact_scores
from components.possible_reconstruction import (
    METHOD, POLICY, reconstruct_possible_arrays, validate_possible_reconstruction_basis,
    possible_reconstruction_summary, validate_possible_reconstruction_audit,
)
from components.reconstruction_selection import CandidateValidation, select_reconstructions
from components.reconstruction_transaction import _digest, reconstruct_accepted_arrays
from tests.test_candidate_evidence import PARAMETERS
from tests.test_reconstruction_transaction import fixture as bf_fixture
from tests.test_reconstruction_boundary_family import fixture as family_fixture


def fixture(*, equal_period=False, boundary_tie=False):
    """Source-valid paths and model-real synthetic vectors; no predictor.

    This models a qualified Possible core and a lower-scored PD extension.
    It does not claim the synthetic feature values were independently extracted.
    """
    if boundary_tie:
        source, categories, scorer = family_fixture()
    else:
        source, categories, scorer = bf_fixture(equal_period=equal_period)
    key = next(iter(categories[0]))

    def score_batch(requests):
        rows = scorer(requests)
        for values in rows.values():
            for row in values:
                row[2][0][0, 9] = 0.
                row[0] = float(predict_compact_scores(row[2][0])[0])
        return rows

    core = score_batch({key: [categories[0][key][0][1]]})[key]
    donors = score_batch({key: [row[1] for row in categories[1][key]]})[key]
    return source, [{}, {}, {key: core}, {key: donors}, {}], score_batch


def qualified(row):
    return CandidateValidation(row, True, True, True, True, {})


class PossibleReconstructionTests(unittest.TestCase):
    def run_fixture(self, **kwargs):
        source, categories, scorer = fixture(**kwargs)
        output, audit = reconstruct_possible_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        return source, categories, output, audit

    def test_actual_nc019_source_paths_with_explicit_fixture_only_scoring(self):
        data = json.loads((Path(__file__).parent / "fixtures/nc019892_possible_extent.json").read_text())
        source = data["source_window"]["sequence"]
        offset = data["source_window"]["start"]
        self.assertEqual(hashlib.sha256(source.encode()).hexdigest(), data["source_window"]["sha256"])
        categories = [{}, {}, {}, {}, {}]
        vectors = {}
        for observation in data["complete_observations"]:
            row = observation["observation"]
            intervals = [(a - offset, b - offset) for a, b in row["repeat_intervals"]]
            candidate = CrisprCandidate(
                list_repeat_starts=[a for a, _b in intervals],
                list_repeats=[source[a:b] for a, b in intervals],
                list_repeats_gaped=row["gapped_repeats"],
                list_spacers=[source[a[1]:b[0]] for a, b in zip(intervals, intervals[1:])])
            vector = np.array([[row["feature_values"][key] for key in COMPACT_FEATURE_SCHEMA]])
            self.assertEqual(float(predict_compact_scores(vector)[0]), row["score"])
            # Coordinates are synthetic local translations; constructor origins
            # are synthetic. Do not attach absolute historical birth attestations.
            vectors[id(candidate)] = vector
            categories[row["category_index"]].setdefault((0, len(source)), []).append(
                [row["score"], candidate, [vector]])
        def score(requests):
            # Stored model-real vectors only; no external extraction is claimed.
            return {key: [[float(predict_compact_scores(vectors[id(c)])[0]), c, [vectors[id(c)].copy()]]
                          for c in values] for key, values in requests.items()}
        output, audit = reconstruct_possible_arrays(source, categories, PARAMETERS, score, lambda: ())
        chosen = output[2][(0, len(source))][0]
        self.assertEqual(len(chosen[1].list_repeats), 189)
        self.assertEqual(chosen[0], .7226947397300213)
        relation = chosen[3]["classification_basis"]["reconstruction"]["selection"]["relationship"]
        self.assertEqual(relation["correspondence"]["anchor_to_candidate_indexes"], [105, 106, 107])
        self.assertNotIn("added_repeat_reference", relation)
        self.assertAlmostEqual(relation["minimum_added_repeat_identity"], .8378378378378378)
        self.assertEqual(audit["stats"]["added_repeats"], 186)
        self.assertEqual(sorted(len(row[1].list_repeats) for row in output[3][(0, len(source))]), [3, 164, 188])
        self.assertEqual(chosen[3]["classification_basis"]["method"], METHOD)

    def test_original_witness_survives_selected_proof_only_metadata(self):
        from components.candidate_origins import mark_repair_origin
        from components.initial_competition import initial_competition_groups, validate_ordinary_witness_inventory
        source, categories, scorer = fixture()
        key = next(iter(categories[2]))
        core, donor = categories[2][key][0], categories[3][key][0]
        mark_repair_origin(donor[1], {"policy": "terminal-only-trim-v1", "status": "generated",
            "coordinate_system": "zero_based_half_open", "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry": [donor[1].list_repeat_starts,
                [a + len(r) for a, r in zip(donor[1].list_repeat_starts, donor[1].list_repeats)]]})
        groups, conflicts = initial_competition_groups([core, donor])
        self.assertFalse(conflicts)
        categories[3][key] = [next(row for row in groups[0].rows if row[1] is donor[1])]
        out, _audit = reconstruct_possible_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        chosen = out[2][key][0]
        self.assertEqual(len(validate_ordinary_witness_inventory(out, source)), 1)
        self.assertNotIn("initial_competition", chosen[3])
        basis = chosen[3]["classification_basis"]
        validate_possible_reconstruction_basis(source, chosen[1], basis, chosen[0], 2, chosen, out)
        missing = list(out)
        missing[3] = {}
        with self.assertRaisesRegex(ValueError, "ordinary occurrence witness"):
            validate_ordinary_witness_inventory(missing, source)

    def test_s1_nontransitive_conflict_and_symmetric_proposal_collisions(self):
        from tests.test_initial_competition import row, evaluate, members
        ordinary = row(.68, [(10, 16), (40, 46), (70, 76)])
        repaired = row(.73, [(13, 19), (43, 49), (73, 79)], repair=True)
        inputs = evaluate([ordinary, repaired])
        core = members(inputs, 2)[0]
        first = row(.65, [(13, 19), (43, 49), (73, 79), (103, 109)])
        second = row(.65, [(13, 20), (43, 49), (73, 79), (103, 109)])
        for item in (core, first, second):
            item[1].consensus = "AAAAAA"
        inputs[3].setdefault((0, 100000000), []).extend([first, second])
        result = select_reconstructions(inputs, qualified, POLICY, _possible=True)
        self.assertEqual(result.decisions[0].reason, "same_rank_peer_veto")
        self.assertTrue(any("protected_ordinary_occurrence_conflict" in peer["vetoes"]
                            for peer in result.decisions[0].evidence["same_rank_peers"]))
        from tests.test_reconstruction_selection import row as simple_row
        left, right = simple_row([100, 140, 180], .7), simple_row([500, 540, 580], .7)
        a = simple_row([100, 140, 180, 220, 260, 300, 340], .65)
        b = simple_row([340, 390, 430, 470, 500, 540, 580], .65)
        for anchors in ([left, right], [right, left]):
            decision = select_reconstructions([{}, {}, {"core": anchors}, {"donors": [a, b]}, {}],
                                               qualified, POLICY, _possible=True)
            self.assertEqual([d.reason for d in decision.decisions], ["proposal_collision"] * 2)

    def test_recorded_multi_owner_scope_does_not_reenter_after_other_core_is_absent(self):
        from tests.test_initial_competition import row
        from components.initial_competition import initial_competition_groups
        a = row(.71, [(100, 110), (140, 150), (180, 190)])
        b = row(.70, [(300, 310), (340, 350), (380, 390)])
        bridge = row(.65, [(i, i + 10) for i in (100, 140, 180, 220, 260, 300, 340, 380)], repair=True)
        _groups, conflicts = initial_competition_groups([a, b, bridge])
        self.assertEqual(len(conflicts), 1)
        for item in (a, b, bridge):
            item[1].consensus = "A" * 10
        result = select_reconstructions([{}, {"prior": [b]}, {"core": [a]}, {"bridge": list(conflicts)}, {}],
                                        qualified, POLICY, _possible=True)
        self.assertEqual(result.stats["selected"], 0)
        self.assertEqual(result.decisions[0].evidence["geometry_rejections"]["repair_scope_conflict"], 1)

    def test_final_bf_actual_family_union_is_collision_barrier(self):
        from tests.test_competing_terminal_extents import source_candidate
        source, categories, scorer = bf_fixture()
        key = next(iter(categories[1]))
        donor = categories[1][key][0][1]
        repeats = donor.list_repeats[:]
        repeats[-1] = repeats[-1][:-6]
        short = CrisprCandidate(repeats, [r + ("-" * 6 if i == 5 else "") for i, r in enumerate(repeats)],
                                donor.list_spacers[:], donor.list_repeat_starts[:])
        categories[1][key].append(scorer({key: [short]})[key][0])
        motif = source[431:460]
        source += "ACGTACG"
        intervals = []
        for spacer in ("AC" * 16, "GT" * 16, "AG" * 16):
            start = len(source)
            source += motif
            intervals.append((start, len(source)))
            source += spacer
        categories[2][(431, len(source))] = [[.7, source_candidate(source, intervals), []]]
        categories[3][(431, len(source))] = [[.65, source_candidate(source, [(431, 460)] + intervals), []]]
        completed, _audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        family = completed[0][key][0]
        self.assertEqual(family[1].list_repeat_starts[-1] + len(family[1].list_repeats[-1]), 431)
        result = select_reconstructions(completed, qualified, POLICY, _possible=True)
        self.assertEqual(result.stats["selected"], 0)
        self.assertEqual(result.decisions[0].evidence["geometry_rejections"]["overlaps_other_accepted_core"], 1)

    def test_empty_and_nonselected_audits_are_validated_and_duplicate_proofs_rejected(self):
        source, categories, scorer = fixture(boundary_tie=True)
        _out, audit = reconstruct_possible_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        report = {"source": {"sha256": hashlib.sha256(source.encode()).hexdigest()}, "arrays": []}
        validate_possible_reconstruction_audit(report, audit)
        for mutation in ("source", "bool_count", "reason_count", "schema", "tier"):
            changed = copy.deepcopy(audit)
            if mutation == "source":
                changed["source_sha256"] = "0" * 64
            elif mutation == "bool_count":
                changed["stats"]["anchors"] = True
            elif mutation == "reason_count":
                changed["stats"]["reason_counts"] = {"unrelated": 1}
            elif mutation == "schema":
                changed["schema"]["version"] = "future"
            else:
                changed["confidence_category"] = "Bona-fide"
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                validate_possible_reconstruction_audit(report, changed)
        source, _before, after, audit = self.run_fixture()
        chosen = next(iter(after[2].values()))[0]
        report = {"source": {"sha256": hashlib.sha256(source.encode()).hexdigest()},
                  "arrays": [{"classification_basis": chosen[3]["classification_basis"]}] * 2}
        with self.assertRaisesRegex(ValueError, "duplicates"):
            validate_possible_reconstruction_audit(report, audit)

    def test_possible_tier_fresh_own_score_and_original_pd_core(self):
        source, before, after, audit = self.run_fixture()
        key = next(iter(after[2]))
        chosen = after[2][key][0]
        core = before[2][key][0]
        self.assertEqual(audit["stats"]["selected"], 1)
        self.assertIs(after[3][key][0], core)
        self.assertEqual(len(before[3][key]), 1)
        self.assertEqual(len(after[0]), 0)
        self.assertLess(chosen[0], core[0])
        summary = possible_reconstruction_summary(chosen)
        self.assertEqual(summary["confidence_category"], "Possible")
        self.assertEqual(summary["core_repeat_count"], 3)
        self.assertEqual(summary["repeat_count"], 6)
        proof = chosen[3]["classification_basis"]
        self.assertEqual(proof, validate_possible_reconstruction_basis(
            source, chosen[1], json.loads(json.dumps(proof)), chosen[0], 2, chosen, after))
        with self.assertRaisesRegex(ValueError, "Possible core category"):
            validate_possible_reconstruction_basis(source, chosen[1], proof, chosen[0], 0, chosen, after)

    def test_selected_location_only_alias_retention(self):
        source, categories, scorer = fixture()
        key = next(iter(categories[3]))
        donor = categories[3][key][0]
        other = copy.deepcopy(donor)
        other.append({"observation_note": "distinct original alias"})
        categories[3][key].extend([donor, other])
        out, audit = reconstruct_possible_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selected"], 1)
        self.assertEqual(len(out[3][key]), 3)
        self.assertIs(out[3][key][0], donor)
        self.assertIs(out[3][key][1], other)
        self.assertIs(out[3][key][2], categories[2][key][0])

    def test_only_direct_anchors_are_scored_and_original_origin_mutation_fails(self):
        from components.candidate_origins import inherit_candidate_origins
        source, categories, scorer = fixture()
        key = next(iter(categories[2]))
        core, donor = categories[2][key][0], categories[3][key][0]
        offset = len(source)
        source += source
        other = CrisprCandidate(core[1].list_repeats[:], core[1].list_gaped_repeats[:],
                                core[1].list_spacers[:], [a + offset for a in core[1].list_repeat_starts])
        untouched = [core[0], other, copy.deepcopy(core[2])]
        categories[2][(offset, len(source))] = [untouched]
        requested = []
        def score(requests):
            requested.extend(candidate for values in requests.values() for candidate in values)
            return scorer(requests)
        output, _audit = reconstruct_possible_arrays(source, categories, PARAMETERS, score, lambda: ())
        self.assertEqual({id(c) for c in requested}, {id(core[1]), id(donor[1])})
        self.assertIs(output[2][(offset, len(source))][0], untouched)
        def mutate(requests):
            result = scorer(requests)
            inherit_candidate_origins(donor[1], (core[1],))
            return result
        with self.assertRaisesRegex(ValueError, "mutated original provenance"):
            reconstruct_possible_arrays(source, categories, PARAMETERS, mutate, lambda: ())

    def test_floor_no_low_or_bf_reclassification(self):
        for score, selected in ((.60, 1), (.60 - 1e-10, 0), (.5540780107828767, 0)):
            _source, categories, _scorer = fixture()
            donor = next(iter(categories[3].values()))[0]
            donor[0] = score
            result = select_reconstructions(categories, qualified, POLICY, _possible=True)
            self.assertEqual(result.stats["selected"], selected)
        _source, categories, _scorer = fixture()
        categories[4] = categories[3]
        categories[3] = {}
        self.assertEqual(select_reconstructions(categories, qualified, POLICY, _possible=True).stats["selected"], 0)

    def test_fresh_floor_and_ineligible_or_negative_candidate_abstain(self):
        for kind in ("fresh_below", "lineage", "known_negative"):
            source, categories, scorer = fixture()
            donor = next(iter(categories[3].values()))[0]
            if kind == "lineage":
                donor[1].source_path_classification_eligible = False
                donor[1].source_path_classification_reason = "broken_local_edge_lineage_parent"
            if kind == "known_negative":
                donor.append({"tags": ["tandem"]})
            def score(requests):
                rows = scorer(requests)
                if kind == "fresh_below":
                    for values in rows.values():
                        for row in values:
                            if row[1] is donor[1]:
                                row[0] = .599
                return rows
            out, audit = reconstruct_possible_arrays(source, categories, PARAMETERS, score, lambda: ())
            self.assertEqual(audit["stats"]["selected"], 0)
            self.assertIs(next(iter(out[2].values()))[0], next(iter(categories[2].values()))[0])

    def test_tied_geometry_abstention_same_rank_veto_and_lower_rank_failure(self):
        _source, categories, _scorer = fixture(boundary_tie=True)
        bad = next(iter(categories[3].values()))[1]
        tied = select_reconstructions(categories, qualified, POLICY, _possible=True).decisions[0]
        self.assertEqual(tied.reason, "ambiguous_possible_geometry")
        self.assertEqual(len(tied.tied_geometry_ids), 2)
        for failed in (False, None):
            def validate(row):
                return replace(qualified(row), negatives_passed=failed) if row is bad else qualified(row)
            self.assertEqual(select_reconstructions(categories, validate, POLICY, _possible=True).decisions[0].reason,
                             "same_rank_peer_veto")
        def lower(row):
            return replace(qualified([.61, row[1], row[2]]), negatives_passed=False) if row is bad else qualified(row)
        self.assertEqual(select_reconstructions(categories, lower, POLICY, _possible=True).stats["selected"], 1)

    def test_bf_transaction_and_family_output_are_unchanged(self):
        for make in (bf_fixture, family_fixture):
            source, categories, scorer = make()
            baseline, bf_audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
            expected = copy.deepcopy(bf_audit)
            out, possible = reconstruct_possible_arrays(source, baseline, PARAMETERS, scorer, lambda: ())
            self.assertEqual(bf_audit, expected)
            self.assertEqual(possible["stats"]["selected"], 0)
            for key, rows in baseline[0].items():
                for original, current in zip(rows, out[0][key]):
                    self.assertIs(original, current)

    def test_typed_proof_source_feature_identity_and_original_core_tampering(self):
        source, _before, after, _audit = self.run_fixture()
        chosen = next(iter(after[2].values()))[0]
        original = chosen[3]["classification_basis"]
        for kind in ("rule", "source", "score", "feature", "bool_mapping", "geometry_id", "origins", "initial"):
            basis = copy.deepcopy(original)
            proof = basis["reconstruction"]
            if kind == "rule":
                proof["rule"] = "unknown"
            elif kind == "source":
                proof["source_sha256"] = "0" * 64
            elif kind == "score":
                proof["reconstruction"]["score"] += .01
            elif kind == "feature":
                proof["donor_original"]["features"]["repeat_len"] += 1
                proof["donor_observation"]["snapshot"] = copy.deepcopy(proof["donor_original"])
            elif kind == "bool_mapping":
                proof["selection"]["relationship"]["correspondence"]["anchor_to_candidate_indexes"][1] = True
            elif kind == "geometry_id":
                proof["selection"]["selected"]["geometry_id"] = "0" * 64
            elif kind == "origins":
                proof["donor_observation"]["origins"]["ordinary_birth_sha256"] = "0" * 64
            else:
                proof["anchor_initial_competition"] = {"rule": "unknown"}
            basis["evidence_sha256"] = _digest(proof)
            with self.subTest(kind=kind), self.assertRaises((ValueError, KeyError)):
                validate_possible_reconstruction_basis(source, chosen[1], basis, chosen[0], 2, chosen, after)
        missing = list(after)
        missing[3] = {}
        with self.assertRaisesRegex(ValueError, "retained exact original core"):
            validate_possible_reconstruction_basis(source, chosen[1], original, chosen[0], 2, chosen, missing)


if __name__ == "__main__":
    unittest.main()
