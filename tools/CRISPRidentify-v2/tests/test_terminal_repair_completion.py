import copy
import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace

from components.components_evaluated_arrays_enhancement import build_source_backed_candidate_from_boundaries as build
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.report_model import ReportModelBuilder, ReportModelError
from components.repeat_boundary_repair import propose_terminal_repairs, source_geometry
from components.terminal_repair_completion import (
    METHOD, apply_terminal_repair_completion, plan_terminal_repair_completion,
    replay_terminal_completion_basis,
)
from tests.test_repeat_boundary_repair import CORE, PARAMETERS, synthetic


class TerminalRepairCompletionTests(unittest.TestCase):
    def setUp(self):
        self.structural = build_advanced_fuzzy_search_filter(PARAMETERS)

    def fixture(self, core_count=3, count=7, old_score=.79, new_score=.76,
                intermediate=False, left=False):
        dna, parent = synthetic(left=left, count=count)
        if intermediate:
            starts, ends = source_geometry(parent)
            bases = list(dna)
            for index, end in enumerate(ends):
                bases[end] = "AGCTACG"[index % 7]
            dna = "".join(bases)
            parent = build(dna, None, starts, [end + 1 for end in ends])
        children, records, _ = propose_terminal_repairs(dna, [parent], self.structural)
        long = children[0]
        starts, ends = source_geometry(parent)
        if intermediate:
            ends = tuple(end - 1 for end in ends)
        core = build(dna, None, starts[:core_count], ends[:core_count])
        # The incumbent's original geometry can be inside the trim corridor.
        if not self.structural(core):
            starts, ends = source_geometry(long)
            core = build(dna, None, starts[:core_count], ends[:core_count])
        category = 0 if old_score >= .75 else 2
        old = [old_score, core, []]
        new = [new_score, long, []]
        categories = [{}, {}, {}, {}, {}]
        categories[category][(0, len(dna))] = [old]
        categories[1 if category == 0 and new_score >= .6 else 3][(0, len(dna))] = [new]
        return dna, parent, long, records, categories, old, new

    def plan(self, fixture):
        dna, _, _, records, categories, _, _ = fixture
        return plan_terminal_repair_completion(dna, categories, records, self.structural)

    def apply(self, fixture, plan=None, scorer=None):
        dna, _, _, _, categories, _, new = fixture
        if scorer is None:
            def scorer(requested):
                return {key: [[new[0], c, []] for c in candidates]
                        for key, candidates in requested.items()}
        return apply_terminal_repair_completion(dna, categories, plan or self.plan(fixture),
                                               self.structural, scorer)

    def test_lower_relative_score_completes_core_and_preserves_original(self):
        for kwargs in ({}, {"old_score": .739, "new_score": .736},
                       {"old_score": .7, "new_score": .8},
                       {"intermediate": True}, {"left": True}):
            with self.subTest(kwargs=kwargs):
                fixture = self.fixture(**kwargs)
                dna, _, long, records, categories, old, new = fixture
                plan = self.plan(fixture)
                self.assertEqual(1, len(plan.replacements))
                result, audit = self.apply(fixture, plan)
                accepted = [row for index in (0, 2) for rows in result[index].values() for row in rows]
                self.assertEqual(1, len(accepted))
                self.assertIs(long, accepted[0][1])
                self.assertEqual(new[0], accepted[0][0])
                self.assertEqual(METHOD, accepted[0][3]["classification_basis"]["method"])
                self.assertTrue(any(row is old for index in (1, 3) for rows in result[index].values() for row in rows))
                repeated = plan_terminal_repair_completion(dna, result, records, self.structural)
                self.assertEqual((), repeated.replacements)
                category = 0 if new[0] >= .75 else 2
                evidence = audit[-1]["evidence"]
                replay_terminal_completion_basis(dna, long, json.loads(json.dumps(evidence)), new[0], category)

    def test_existing_geometry_attestation_also_authorizes_completion(self):
        fixture = list(self.fixture())
        dna, parent, long = fixture[:3]
        _, records, _ = propose_terminal_repairs(dna, [parent, long], self.structural)
        self.assertEqual("child_geometry_already_present", records[0]["reason"])
        fixture[3] = records
        self.assertEqual(1, len(self.plan(fixture).replacements))

    def test_wrong_hash_restrictions_threshold_and_one_extra_abstain(self):
        fixture = self.fixture()
        fixture[3][0]["source_sha256"] = "0" * 64
        self.assertEqual((), self.plan(fixture).replacements)
        for location in ("parent", "long", "core", "duplicate"):
            fixture = self.fixture()
            dna, parent, long, _, categories, old, _ = fixture
            restricted = {"parent": parent, "long": long, "core": old[1], "duplicate": copy.deepcopy(long)}[location]
            restricted.source_path_classification_eligible = False
            if location in ("parent", "duplicate"):
                categories[4].setdefault((0, len(dna)), []).append([.8, restricted, []])
            self.assertEqual((), self.plan(fixture).replacements)
        for old, new in ((.8, .74), (.7, .49)):
            self.assertEqual((), self.plan(self.fixture(old_score=old, new_score=new)).replacements)
        self.assertEqual((), self.plan(self.fixture(old_score=.4, new_score=.8)).replacements)
        self.assertEqual((), self.plan(self.fixture(count=6, core_count=5)).replacements)

    def test_other_accepted_overlap_across_keys_and_noncontiguous_core_abstain(self):
        fixture = self.fixture()
        dna, _, long, _, categories, old, _ = fixture
        starts, ends = source_geometry(long)
        other = build(dna, None, starts[-3:], ends[-3:])
        categories[2][(500, 999)] = [[.65, other, []]]
        self.assertEqual((), self.plan(fixture).replacements)
        fixture = self.fixture()
        dna, parent, _, _, _, old, _ = fixture
        starts, ends = source_geometry(parent)
        old[1] = build(dna, None, [starts[i] for i in (0, 2, 3)], [ends[i] for i in (0, 2, 3)])
        self.assertEqual((), self.plan(fixture).replacements)

    def test_nondegenerate_and_individual_spacer_guards(self):
        for reason in ("different_repeat", "oversized", "hidden_motif"):
            fixture = list(self.fixture())
            dna, parent = fixture[:2]
            starts, ends = source_geometry(parent)
            if reason == "different_repeat":
                position = starts[-1] + 8
                dna = dna[:position] + ("A" if dna[position] != "A" else "C") + dna[position + 1:]
            else:
                position = ends[3] + 10
                insert = "A" * 50 if reason == "oversized" else CORE
                dna = dna[:position] + insert + dna[position:]
                starts = tuple(s + len(insert) if s > position else s for s in starts)
                ends = tuple(e + len(insert) if e > position else e for e in ends)
            parent = build(dna, None, starts, ends)
            children, records, _ = propose_terminal_repairs(dna, [parent], self.structural)
            self.assertEqual(1, len(children))
            long = children[0]
            old = [.79, build(dna, None, starts[:3], ends[:3]), []]
            categories = [{(0, len(dna)): [old]}, {(0, len(dna)): [[.76, long, []]]}, {}, {}, {}]
            fixture = (dna, parent, long, records, categories, old, categories[1][(0, len(dna))][0])
            self.assertEqual((), self.plan(fixture).replacements, reason)

    def test_incomparable_extensions_abstain_regardless_of_order(self):
        dna, parent = synthetic(count=10)
        starts, ends = source_geometry(parent)
        parents = [build(dna, None, starts[:7], ends[:7]), build(dna, None, starts[3:], ends[3:])]
        children, records, _ = propose_terminal_repairs(dna, parents, self.structural)
        self.assertEqual(2, len(children))
        core = build(dna, None, starts[3:6], [e - 1 for e in ends[3:6]])
        categories = [{(0, len(dna)): [[.79, core, []]]},
                      {(0, len(dna)): [[.76, child, []] for child in children]}, {}, {}, {}]
        for sequence in (records, list(reversed(records))):
            plan = plan_terminal_repair_completion(dna, categories, sequence, self.structural)
            self.assertEqual((), plan.replacements)
            self.assertTrue(any(r["reason"] == "incomparable_completions" for r in plan.records))

    def test_stale_geometry_scores_eligibility_and_scoring_mutation_reject(self):
        for mutation in ("geometry", "score", "eligibility", "lineage"):
            fixture = self.fixture()
            plan = self.plan(fixture)
            new = fixture[-1]
            if mutation == "geometry":
                new[1].list_repeat_starts[0] += 1
            elif mutation == "score":
                new[0] = .77
            elif mutation == "eligibility":
                new[1].source_path_classification_eligible = False
            else:
                new[1].source_path_provenance = (SimpleNamespace(derivation="local_chain_family_parent"),)
            with self.assertRaises(ValueError):
                self.apply(fixture, plan)
        fixture = self.fixture()

        def mutate(requested):
            fixture[-1][1].source_path_classification_eligible = False
            return {key: [[.76, c, []] for c in candidates] for key, candidates in requested.items()}

        with self.assertRaises(ValueError):
            self.apply(fixture, scorer=mutate)

    def test_report_replays_completion_and_rejects_forged_evidence(self):
        fixture = self.fixture()
        dna = fixture[0]
        result, _ = self.apply(fixture)
        row = next(iter(result[0].values()))[0]
        basis = row[3]["classification_basis"]["terminal_repair_completion"]
        for field, value in (("core_score", True), ("completed_score", float("nan")),
                             ("eligibility_checked", 1), ("completed_category_index", False)):
            bad = copy.deepcopy(basis)
            bad[field] = value
            with self.assertRaises(ValueError):
                replay_terminal_completion_basis(dna, row[1], bad, row[0], 0)
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "source.fasta"
            fasta.write_text(">source\n" + dna + "\n")
            model = ReportModelBuilder(file_path=fasta, categories=result, non_array_data={},
                                       list_feature_names=[]).build()
            self.assertEqual(METHOD, model["arrays"][0]["classification_basis"]["method"])
            basis["source_sha256"] = "0" * 64
            with self.assertRaises(ReportModelError):
                ReportModelBuilder(file_path=fasta, categories=result, non_array_data={},
                                   list_feature_names=[]).build()


if __name__ == "__main__":
    unittest.main()
