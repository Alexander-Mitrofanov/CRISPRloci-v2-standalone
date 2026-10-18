"""Adversarial contracts for core-versus-extent reconstruction."""

import copy
from dataclasses import replace
from types import SimpleNamespace
import unittest

from components.reconstruction_selection import (
    CandidateValidation,
    ReconstructionPolicy,
    candidate_geometry,
    candidate_geometry_id,
    monotone_repeat_mapping,
    reconstruction_candidate_locations,
    relationship,
    row_locations,
    select_reconstructions,
)


REPEAT = "ACGTTGCAAC"


def candidate(starts, repeats=None, alignment_suffix=""):
    repeats = list(repeats or [REPEAT] * len(starts))
    return SimpleNamespace(
        list_repeat_starts=list(starts), list_repeats=repeats,
        list_repeats_gaped=[repeat + alignment_suffix for repeat in repeats],
        list_spacers=["N" * (right - left - len(repeats[index]))
                      for index, (left, right) in enumerate(zip(starts, starts[1:]))],
        consensus=REPEAT, consensus_gaped=REPEAT + alignment_suffix,
    )


def row(starts, score=0.7, **kwargs):
    return [score, candidate(starts, **kwargs), [[float(len(starts))]], {"sentinel": "unchanged"}]


def categories(anchor=None, donors=(), possible=()):
    return [{"core": [anchor or row([100, 140, 180], 0.9)]},
            {"different_search_key": list(donors)} if donors else {},
            {"possible": list(possible)} if possible else {}, {}, {}]


def valid(scored_row):
    return CandidateValidation(scored_row, True, True, True, True,
                               {"source": "passed", "negative_checks": "passed"})


class GeometryMappingTests(unittest.TestCase):
    def test_internal_insertion_is_distinct_from_flank_extension(self):
        core = candidate([100, 180, 260])
        expanded = candidate([60, 100, 140, 180, 220, 260, 300])
        relation = relationship(core, expanded)
        self.assertTrue(relation["qualified"])
        self.assertEqual(relation["kind"], "internal_and_flank_extension")
        self.assertEqual(relation["correspondence"]["anchor_to_candidate_indexes"], [1, 3, 5])
        self.assertEqual(relation["internal_added_repeat_indexes"], [2, 4])
        self.assertEqual(relation["flank_added_repeat_indexes"], [0, 6])
        internal = relationship(core, candidate([100, 140, 180, 220, 260]))
        self.assertEqual(internal["kind"], "internal_repeat_insertion")

    def test_shift_tolerance_and_reciprocal_overlap_are_both_enforced(self):
        core = candidate_geometry(candidate([100, 140, 180]))
        shifted = candidate_geometry(candidate([99, 139, 179, 219]))
        relation = monotone_repeat_mapping(core, shifted)
        self.assertEqual(relation.status, "unique")
        self.assertEqual(relation.maximum_boundary_shift, 1)
        # Three-base shift satisfies the distance limit but not 80% of ten bp.
        self.assertEqual(monotone_repeat_mapping(core,
                         candidate_geometry(candidate([97, 137, 177, 217]))).status,
                         "no_mapping")
        strict = ReconstructionPolicy(boundary_tolerance=0)
        self.assertEqual(monotone_repeat_mapping(core, shifted, strict).status, "no_mapping")

    def test_ambiguous_mapping_is_not_hidden_by_path_convergence(self):
        policy = ReconstructionPolicy(boundary_tolerance=5, minimum_reciprocal_overlap=0.4)
        # The first anchor can map to either of two disjoint intervals; both
        # choices converge to the same second repeat and must remain ambiguous.
        result = monotone_repeat_mapping(((5, 15), (30, 40)),
                                        ((0, 10), (10, 20), (30, 40)), policy)
        self.assertEqual(result.status, "ambiguous_mapping")

    def test_incompatible_phase_does_not_become_internal_insertion(self):
        self.assertFalse(relationship(candidate([100, 140, 180]),
                                      candidate([115, 135, 155, 175]))["qualified"])

    def test_local_added_repeat_identity_is_required(self):
        core = candidate([100, 140, 180])
        donor = candidate([100, 140, 180, 220],
                          repeats=[REPEAT, REPEAT, REPEAT, "AAAAAAAAAA"])
        self.assertEqual(relationship(core, donor)["reason"], "added_repeat_identity_below_limit")
        donor.list_repeats[-1] = "ACGTTGCAAN"
        self.assertEqual(relationship(core, donor)["reason"], "added_repeat_not_canonical")

    def test_maximum_spacer_applies_to_each_gap(self):
        relation = relationship(candidate([100, 140, 180]), candidate([100, 140, 180, 340]))
        self.assertEqual(relation["reason"], "spacer_above_continuity_limit")

    def test_geometry_id_ignores_alignment_and_score(self):
        first = candidate([100, 140, 180])
        second = candidate([100, 140, 180], alignment_suffix="-")
        self.assertEqual(candidate_geometry_id(first), candidate_geometry_id(second))
        with self.assertRaisesRegex(ValueError, "nonoverlapping"):
            candidate_geometry(candidate([100, 105, 180]))


class ReconstructionSelectionTests(unittest.TestCase):
    def test_chooses_own_fresh_score_without_mutating_input(self):
        donor = row([60, 100, 140, 180, 220], 0.73)
        inputs = categories(donors=[donor])
        original = copy.deepcopy(inputs)

        def refresh(scored):
            fresh = list(scored)
            if scored is donor:
                fresh[0] = 0.61
                fresh[2] = [[5.0, 17.0]]
            return CandidateValidation(fresh, True, True, True, True, {"fresh": True})

        result = select_reconstructions(inputs, refresh)
        selected = result.decisions[0]
        self.assertEqual(selected.action, "selected")
        self.assertEqual(selected.validated_row[0], 0.61)
        self.assertEqual(selected.validated_row[2], [[5.0, 17.0]])
        self.assertEqual(selected.anchor.row[0], 0.9)
        self.assertEqual(inputs, original)
        self.assertEqual(result.stats["added_repeats"], 2)

    def test_internal_missed_repeats_are_actual_selection_candidates(self):
        inputs = categories(row([100, 180, 260], 0.9),
                            [row([100, 140, 180, 220, 260])])
        result = select_reconstructions(inputs, valid)
        self.assertEqual(result.stats["selected"], 1)
        self.assertEqual(result.decisions[0].relation["kind"], "internal_repeat_insertion")

    def test_floor_and_all_four_validation_gates_fail_closed(self):
        for attribute in ("source_valid", "structural_valid", "negatives_passed", "score_fresh"):
            for failure in (False, None):
                with self.subTest(attribute=attribute, failure=failure):
                    inputs = categories(donors=[row([60, 100, 140, 180])])
                    donor = inputs[1]["different_search_key"][0]

                    def validation(scored):
                        result = valid(scored)
                        return replace(result, **{attribute: failure}) if scored is donor else result

                    result = select_reconstructions(inputs, validation)
                    self.assertEqual(result.stats["selected"], 0)
                    self.assertEqual(result.decisions[0].reason, "no_validated_extension")
        inputs = categories(donors=[row([60, 100, 140, 180])])
        donor = inputs[1]["different_search_key"][0]

        def lowered(scored):
            return valid([0.59, scored[1], scored[2]]) if scored is donor else valid(scored)

        self.assertEqual(select_reconstructions(inputs, lowered).stats["selected"], 0)

    def test_unqualified_core_never_authorizes_reconstruction(self):
        inputs = categories(donors=[row([60, 100, 140, 180])])
        anchor = inputs[0]["core"][0]

        def failed_core(scored):
            return replace(valid(scored), negatives_passed=False) if scored is anchor else valid(scored)

        self.assertEqual(select_reconstructions(inputs, failed_core).decisions[0].reason,
                         "anchor_not_qualified")
        inputs[0]["core"][0][0] = 0.74
        self.assertFalse(reconstruction_candidate_locations(inputs))

    def test_equal_best_distinct_phases_abstain_in_both_input_orders(self):
        left = row([60, 100, 140, 180])
        right = row([100, 140, 180, 220])
        for donors in ([left, right], [right, left]):
            decision = select_reconstructions(categories(donors=donors), valid).decisions[0]
            self.assertEqual(decision.reason, "ambiguous_best_geometry")
            self.assertEqual(len(decision.tied_geometry_ids), 2)
            self.assertIsNone(decision.selected)

    def test_alignment_aliases_are_not_votes_and_all_locations_are_exposed(self):
        first = row([60, 100, 140, 180], 0.70)
        alias = row([60, 100, 140, 180], 0.72, alignment_suffix="-")
        inputs = categories(donors=[first, alias])
        inputs[1]["duplicate_scope"] = [first]
        decision = select_reconstructions(inputs, valid).decisions[0]
        self.assertEqual(decision.action, "selected")
        self.assertIs(decision.selected.row, alias)
        self.assertEqual(len(decision.selected_locations), 3)
        self.assertEqual(decision.evidence["qualified_extension_geometries"], 1)
        self.assertEqual(len(reconstruction_candidate_locations(inputs)), 2)

    def test_other_accepted_core_protects_against_fusion_even_when_compatible(self):
        donor = row([100, 140, 180, 220, 260, 300])
        inputs = categories(donors=[donor], possible=[row([220, 260, 300], 0.65)])
        self.assertFalse(reconstruction_candidate_locations(inputs))
        result = select_reconstructions(inputs, valid)
        self.assertEqual(result.stats["selected"], 0)
        self.assertEqual(result.decisions[0].evidence["geometry_rejections"],
                         {"overlaps_other_accepted_core": 1})

    def test_overlapping_proposals_abstain_symmetrically(self):
        first = row([100, 140, 180], 0.90)
        second = row([500, 540, 580], 0.90)
        donor_first = row([100, 140, 180, 220, 260, 300, 340])
        donor_second = row([300, 340, 380, 420, 460, 500, 540, 580])
        for anchors in ([first, second], [second, first]):
            inputs = [{"cores": anchors}, {"donors": [donor_first, donor_second]}, {}, {}, {}]
            result = select_reconstructions(inputs, valid)
            self.assertEqual(result.stats["selected"], 0)
            self.assertEqual([d.reason for d in result.decisions], ["proposal_collision"] * 2)
            self.assertTrue(all(not d.selected_locations for d in result.decisions))

    def test_no_transitive_merging_or_veto_from_unaccepted_witness(self):
        # An incompatible high-score Alternative cannot veto a qualified path.
        good = row([60, 100, 140, 180, 220], 0.65)
        unrelated_phase = row([105, 145, 185, 225], 0.99)
        disjoint = row([260, 300, 340, 380, 420, 460], 0.99)
        result = select_reconstructions(categories(donors=[good, unrelated_phase, disjoint]), valid)
        self.assertIs(result.decisions[0].selected.row, good)

    def test_only_alternatives_donate_not_possible_discarded_or_low_score(self):
        inputs = categories()
        inputs[3] = {"discarded": [row([60, 100, 140, 180], 0.70)]}
        inputs[4] = {"low": [row([60, 100, 140, 180, 220], 0.95)]}
        self.assertFalse(reconstruction_candidate_locations(inputs))
        self.assertEqual(len(row_locations(inputs)), 1)

    def test_validator_cannot_replace_or_mutate_candidate_object(self):
        inputs = categories(donors=[row([60, 100, 140, 180])])

        def replace_candidate(scored):
            changed = [scored[0], copy.deepcopy(scored[1]), scored[2]]
            return valid(changed)

        with self.assertRaisesRegex(ValueError, "exact candidate"):
            select_reconstructions(inputs, replace_candidate)

        def mutate(scored):
            scored[1].list_repeat_starts[0] += 1
            return valid(scored)

        with self.assertRaisesRegex(ValueError, "mutated"):
            select_reconstructions(inputs, mutate)

    def test_invalid_contracts_fail_before_output_commit(self):
        with self.assertRaisesRegex(ValueError, "five categories"):
            select_reconstructions([], valid)
        with self.assertRaisesRegex(ValueError, "candidate floor"):
            ReconstructionPolicy(minimum_candidate_score=0.9)
        with self.assertRaisesRegex(TypeError, "boolean"):
            CandidateValidation(row([100, 140, 180]), True, True, "passed", True)
        inputs = categories(donors=[row([60, 100, 140, 180])])
        with self.assertRaisesRegex(TypeError, "CandidateValidation"):
            select_reconstructions(inputs, lambda _row: True)


if __name__ == "__main__":
    unittest.main()
