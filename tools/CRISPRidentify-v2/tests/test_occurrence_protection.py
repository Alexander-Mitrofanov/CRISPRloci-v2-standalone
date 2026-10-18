"""S1 occurrence safety; no biological truth is inferred from this relation."""
import copy
import itertools
import unittest

from components.initial_competition import (
    OCCURRENCE_RULE, competition_evidence, guard_reconstructed_occurrences,
    initial_competition_groups, with_competition_evidence,
)
from components.occurrence_protection import occurrence_relation
from components.reconstruction_selection import select_reconstructions
from tests.test_initial_competition import row, evaluate, members
from tests.test_reconstruction_selection import valid


class OccurrenceRelationTests(unittest.TestCase):
    def test_half_boundary_and_shifts_are_explicit(self):
        relation = occurrence_relation([(10, 20)], [(15, 25)])
        self.assertEqual(relation["status"], "unique")
        self.assertEqual(relation["pairs"][0]["boundary_shifts"], [5, 5])
        self.assertEqual(occurrence_relation([(10, 20)], [(16, 26)])["status"], "no_complete_mapping")

    def test_complete_uniqueness_not_greedy_or_local_ambiguity(self):
        # First reference has two edges, but the second needs candidate1.
        relation = occurrence_relation([(5, 15), (15, 25)], [(0, 10), (10, 20)])
        self.assertEqual(relation["mapping"], [0, 1])
        self.assertEqual(relation["status"], "unique")
        relation = occurrence_relation([(5, 15), (30, 40)], [(0, 10), (10, 20), (30, 40)])
        self.assertEqual(relation["status"], "ambiguous_mapping")
        self.assertEqual(relation["mapping_count_capped_at_two"], 2)

    def test_global_injectivity_failure_is_not_a_spacer_omission(self):
        relation = occurrence_relation([(0, 10), (10, 20)], [(0, 20)])
        self.assertEqual(relation["status"], "no_complete_mapping")
        self.assertEqual(relation["reference_indexes_without_edge"], [])
        self.assertEqual(relation["reference_indexes_inside_candidate_spacers"], [])

    def test_dp_matches_exhaustive_small_path_oracle(self):
        intervals = [(a, b) for a in range(0, 7, 2) for b in range(a + 2, 9, 2)]
        geometries = [(x,) for x in intervals]
        geometries += [pair for pair in itertools.combinations(intervals, 2) if pair[0][1] <= pair[1][0]]
        for reference, proposed in itertools.product(geometries, repeat=2):
            paths = []
            for indexes in itertools.combinations(range(len(proposed)), len(reference)):
                if all(2 * (min(a[1], proposed[i][1]) - max(a[0], proposed[i][0]))
                       >= max(a[1] - a[0], proposed[i][1] - proposed[i][0])
                       for a, i in zip(reference, indexes)):
                    paths.append(indexes)
            actual = occurrence_relation(reference, proposed)
            self.assertEqual(actual["mapping_count_capped_at_two"], min(2, len(paths)))
            if len(paths) == 1:
                self.assertEqual(actual["mapping"], list(paths[0]))


class OccurrenceCompetitionTests(unittest.TestCase):
    def test_internal_deletion_and_displacement_preserve_ordinary_and_own_scores(self):
        ordinary_geometry = [(100 + 60 * i, 130 + 60 * i) for i in range(69)]
        displaced = ordinary_geometry[:]
        displaced[45] = (ordinary_geometry[45][0] - 22, ordinary_geometry[45][1] - 23)
        for score, proposed, missing, partial in (
                (.8894593126618535, displaced, [45], [45]),
                (.5734534, ordinary_geometry[:6] + ordinary_geometry[7:], [6], [])):
            ordinary = row(score, ordinary_geometry)
            repair = row(score + .00003, proposed, repair=True)
            for inputs in ([ordinary, repair], [repair, ordinary]):
                result = evaluate(inputs)
                accepted = members(result, 0 if score >= .75 else 2)[0]
                self.assertIs(accepted[1], ordinary[1])
                evidence = competition_evidence(accepted)
                self.assertEqual(evidence["rule"], OCCURRENCE_RULE)
                conflict = evidence["occurrence_protection"]["deferred_candidates"][0]
                self.assertEqual(conflict["candidate"]["score"], repair[0])
                self.assertEqual(conflict["relation"]["reference_indexes_without_edge"], missing)
                self.assertEqual(conflict["relation"]["reference_indexes_with_partial_overlap_only"], partial)
                retained = [r for cat in result for rows in cat.values() for r in rows if r[1] is repair[1]][0]
                self.assertEqual(competition_evidence(retained)["decision"], "repair_occurrence_conflict")
                self.assertTrue(retained[1].valid)
                self.assertTrue(retained[1].source_path_classification_eligible)
                self.assertEqual(retained[0], repair[0])

    def test_only_ordinary_winner_is_protected_and_next_supported_repair_can_win(self):
        ordinary = row(.7, [(10,20),(40,50),(70,80)])
        longer_loser = row(.6, [(10,20),(40,50),(70,80),(100,110)])
        unsupported = row(.95, [(10,20),(70,80)], repair=True)
        supported = row(.9, [(11,21),(41,51),(71,81)], repair=True)
        result = evaluate([unsupported, longer_loser, ordinary, supported])
        selected = members(result, 0)[0]
        self.assertIs(selected[1], supported[1])
        evidence = competition_evidence(selected)["occurrence_protection"]
        self.assertEqual(evidence["ordinary_reference"]["score"], .7)
        self.assertEqual(len(evidence["deferred_candidates"]), 1)
        self.assertEqual(sum(len(rows) for cat in result for rows in cat.values()), 4)

    def test_lower_scoring_conflict_is_still_bound_before_later_rescoring(self):
        ordinary = row(.9, [(10,20),(40,50),(70,80)])
        repair = row(.6, [(10,20),(70,80)], repair=True)
        result = evaluate([ordinary, repair])
        retained = members(result, 1)[0]
        self.assertEqual(competition_evidence(retained)["decision"], "repair_occurrence_conflict")
        self.assertEqual(retained[0], .6)

    def test_terminal_omission_is_explicit_unresolved_not_longest_selection(self):
        ordinary = row(.7, [(10, 20), (40, 50), (70, 80)])
        repair = row(.9, [(10, 20), (40, 50)], repair=True)
        result = evaluate([repair, ordinary])
        self.assertIs(members(result, 2)[0][1], ordinary[1])
        self.assertEqual(members(result, 1)[0][0], .9)
        relation = competition_evidence(members(result, 1)[0])["occurrence_protection"]["relation"]
        self.assertEqual(relation["reference_indexes_without_edge"], [2])
        self.assertEqual(relation["reference_indexes_inside_candidate_spacers"], [])

    def test_actual_conserved_core_boundary_correction_can_still_win(self):
        old = [(5314746,5314772),(5314824,5314848),(5314900,5314926),
               (5314940,5314966),(5314983,5315009)]
        new = [(5314744,5314766),(5314822,5314844),(5314898,5314920),
               (5314938,5314960),(5314982,5315003)]
        ordinary, repair = row(.7322381452640958, old), row(.7944477165426679, new, repair=True)
        result = evaluate([ordinary, repair])
        selected = members(result, 0)[0]
        self.assertIs(selected[1], repair[1])
        self.assertEqual(competition_evidence(selected)["occurrence_protection"]["relation"]["mapping"], list(range(5)))

    def test_aliases_remain_rows_not_ambiguity_votes_and_ties_keep_stable_order(self):
        ordinary = row(.8, [(10, 20), (40, 50)])
        repair = row(.8, [(15, 25), (45, 55)], repair=True)
        alias = copy.deepcopy(repair)
        for inputs in ([ordinary, repair, alias], [repair, alias, ordinary]):
            result = evaluate(inputs)
            accepted = members(result, 0)[0]
            self.assertIs(accepted[1], inputs[0][1])
            self.assertEqual(sum(len(rows) for cat in result for rows in cat.values()), 3)
            self.assertEqual(len(competition_evidence(accepted)["tied_candidates"]), 2)

    def test_final_mapping_must_preserve_original_not_transitive_repaired_core(self):
        ordinary = row(.8, [(10, 16), (40, 46), (70, 76)])
        repaired = row(.9, [(13, 19), (43, 49), (73, 79)], repair=True)
        result = evaluate([ordinary, repaired])
        core = members(result, 0)[0]
        donor = row(.7, [(14, 20), (44, 50), (74, 80), (104, 110)])
        core[1].consensus = donor[1].consensus = "AAAAAA"
        result[1].setdefault((0, 100000000), []).append(donor)
        decisions = select_reconstructions(result, valid).decisions
        self.assertEqual(decisions[0].action, "retained")
        self.assertFalse(guard_reconstructed_occurrences(
            [(14,20),(44,50),(74,80),(104,110)], core))
        donor[1].list_repeat_starts[:] = [13,43,73,103]
        donor[1].list_spacers[:] = ["N" * 24] * 3
        self.assertEqual(select_reconstructions(result, valid).decisions[0].action, "selected")

    def test_deferred_donor_carries_its_reference_after_category_movement(self):
        ordinary = row(.8, [(10,16),(40,46),(70,76)])
        repair = row(.9, [(14,20),(44,50),(74,80),(104,110)], repair=True)
        result = evaluate([ordinary, repair])
        deferred = members(result, 1)[0]
        # A separate immediate core could map at 80%; donor's original witness
        # independently prevents laundering its conflicting occurrence path.
        core = row(.95, [(13,19),(43,49),(73,79)])
        core[1].consensus = deferred[1].consensus = "AAAAAA"
        moved = [{"core": [core]}, {"other": [deferred]}, {}, {}, {}]
        self.assertEqual(select_reconstructions(moved, valid).decisions[0].action, "retained")

    def test_v2_candidate_mutations_are_rejected_and_v1_remains_readable(self):
        ordinary, repair = row(.8, [(10,20)]), row(.9, [(15,25)], repair=True)
        selected = members(evaluate([ordinary, repair]), 0)[0]
        for mutation in ("score", "geometry", "sequence", "span", "mapping"):
            altered = copy.deepcopy(selected)
            if mutation == "score":
                altered[0] -= .01
            elif mutation == "geometry":
                altered[1].list_repeat_starts[0] += 1
            elif mutation == "sequence":
                altered[1].list_repeats[0] = "C" * 10
            elif mutation == "span":
                altered[3]["initial_competition"]["occurrence_protection"]["candidate"]["source_span_sha256"] = "0" * 64
            else:
                altered[3]["initial_competition"]["occurrence_protection"]["relation"]["mapping"] = []
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                competition_evidence(altered)
        v1 = initial_competition_groups([ordinary])[0][0].evidence
        self.assertEqual(competition_evidence(with_competition_evidence(ordinary, v1)), v1)


class OccurrencePublicationTests(unittest.TestCase):
    def test_accepted_only_without_sidecar_discloses_conflict_and_binds_witness(self):
        import hashlib
        from pathlib import Path
        import tempfile
        from components.candidate_origins import mark_repair_origin
        from components.module_evaluation import ArrayEvaluation
        from components.report_model import ReportModelBuilder, ReportModelError
        from tests.test_candidate_evidence import _fixture, PARAMETERS
        source, original = _fixture()
        key = next(iter(original[0]))
        repair = original[0][key][0]  # three copies, .90
        ordinary = original[1][key][0]  # five copies, .70
        mark_repair_origin(repair[1], {
            "policy": "terminal-only-trim-v1", "status": "generated",
            "coordinate_system": "zero_based_half_open",
            "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry": [repair[1].list_repeat_starts,
                [a + len(r) for a, r in zip(repair[1].list_repeat_starts, repair[1].list_repeats)]],
        })
        params = dict(PARAMETERS, param_reconstruction_policy="conservative",
                      param_accepted_only=True, param_retain_candidate_evidence=False)
        categories = ArrayEvaluation({}, params, precomputed_scored_rows={key:[repair,ordinary]}).output()
        accepted = members(categories, 2)[0]
        self.assertIs(accepted[1], ordinary[1])
        projected = [categories[0], {}, categories[2], {}, {}]
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "source.fa"
            fasta.write_text(">source\n" + source + "\n")
            report = ReportModelBuilder(fasta, projected, {}, [], parameters=params,
                                         provenance_categories=categories).build()
            evidence = report["arrays"][0]["initial_competition"]
            self.assertEqual(len(evidence["occurrence_protection"]["deferred_candidates"]), 1)
            self.assertEqual(evidence["occurrence_protection"]["deferred_candidates"][0]["candidate"]["score"], .9)
            for field in ("score", "representation_sha256"):
                corrupt = copy.deepcopy(categories)
                snapshot = members(corrupt, 2)[0][3]["initial_competition"]["occurrence_protection"]["deferred_candidates"][0]["candidate"]
                snapshot[field] = .91 if field == "score" else "0" * 64
                with self.subTest(field=field), self.assertRaisesRegex(ReportModelError, "retained deferred repair"):
                    ReportModelBuilder(fasta, corrupt, {}, [], parameters=params).build()
            corrupt = copy.deepcopy(categories)
            from components.candidate_origins import mark_ordinary_birth
            mark_ordinary_birth(members(corrupt, 1)[0][1])
            with self.assertRaisesRegex(ReportModelError, "retained deferred repair"):
                ReportModelBuilder(fasta, corrupt, {}, [], parameters=params).build()
            # Updating a coherent metadata score pair cannot invent a new
            # score for the preserved original ordinary row.
            corrupt = copy.deepcopy(categories)
            current = members(corrupt, 2)[0]
            current[3]["initial_competition"]["ordinary_components"][0]["winner_score"] = .71
            current[3]["initial_competition"]["occurrence_protection"]["ordinary_reference"]["score"] = .71
            with self.assertRaisesRegex(ReportModelError, "ordinary occurrence witness"):
                ReportModelBuilder(fasta, corrupt, {}, [], parameters=params).build()
            corrupt = copy.deepcopy(categories)
            members(corrupt, 2)[0][3]["initial_competition"]["occurrence_protection"]["ordinary_reference"]["source_span_sha256"] = "0" * 64
            with self.assertRaisesRegex(ReportModelError, "source span"):
                ReportModelBuilder(fasta, corrupt, {}, [], parameters=params).build()


if __name__ == "__main__":
    unittest.main()
