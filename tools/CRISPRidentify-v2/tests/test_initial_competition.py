"""Regression contracts for repair-origin scope fencing and observable ties."""
import copy
import gzip
import itertools
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from components.candidate_origins import (
    candidate_origin_kind, candidate_origin_evidence, inherit_candidate_origins,
    mark_ordinary_birth, mark_repair_origin,
)
from components.initial_competition import (
    METADATA_KEY, competition_evidence,
)
from components.module_evaluation import ArrayEvaluation


def row(score, intervals, *, repair=False, eligible=True, valid=True):
    repeats = ["A" * (b - a) for a, b in intervals]
    candidate = SimpleNamespace(
        list_repeat_starts=[a for a, _ in intervals], list_repeats=repeats,
        list_repeats_gaped=repeats[:],
        list_spacers=["N" * (b[0] - a[1]) for a, b in zip(intervals, intervals[1:])],
        source_path_classification_eligible=eligible, valid=valid)
    mark_ordinary_birth(candidate)
    if repair:
        mark_repair_origin(candidate, {
            "policy": "terminal-only-trim-v1", "status": "generated",
            "coordinate_system": "zero_based_half_open", "source_sha256": "a" * 64,
            "child_geometry": [[a for a, _ in intervals], [b for _, b in intervals]],
        })
    return [score, candidate, []]


def evaluate(rows, policy="conservative", split=True):
    with patch("components.module_evaluation.build_advanced_fuzzy_search_filter",
               return_value=lambda candidate: candidate.valid):
        return ArrayEvaluation({}, {"param_reconstruction_policy": policy,
            "param_split_competitive_components": split},
            precomputed_scored_rows={(0, 100000000): rows}).output()


def members(categories, index):
    return [entry for rows in categories[index].values() for entry in rows]


def candidates(categories, index):
    return [id(entry[1]) for entry in members(categories, index)]


class InitialCompetitionTests(unittest.TestCase):
    def test_observed_repair_bridge_preserves_both_possible_cores(self):
        upstream = row(.7045212785746127, [(10250814, 10250839),
            (10250895, 10250920), (10251057, 10251082), (10251084, 10251109)])
        downstream = row(.6926388436032846, [(10251360, 10251386),
            (10251441, 10251467), (10251468, 10251494)])
        bridge = row(.6623743499125787, [(10250338, 10250392),
            (10250392, 10250445), (10250446, 10250500), (10250500, 10250553),
            (10250554, 10250608), (10250608, 10250661), (10250662, 10250716),
            (10250716, 10250770), (10250905, 10250959), (10250959, 10251013),
            (10251148, 10251202), (10251337, 10251391), (10251391, 10251446)], repair=True)
        result = evaluate([upstream, downstream, bridge])
        self.assertEqual(candidates(result, 2), [id(upstream[1]), id(downstream[1])])
        retained = members(result, 3)[0]
        self.assertEqual(retained[0], bridge[0])
        self.assertEqual(competition_evidence(retained)["decision"], "repair_scope_conflict")
        self.assertEqual(len(competition_evidence(retained)["ordinary_components"]), 2)
        self.assertTrue(retained[1].source_path_classification_eligible)
        self.assertTrue(retained[1].valid)
        self.assertEqual(candidates(evaluate([upstream, downstream, bridge], "legacy"), 2), [id(upstream[1])])

    def test_no_repair_parity_preserves_ordinary_long_array_competition(self):
        a, b, long = row(.9, [(10, 20)]), row(.89, [(100, 110)]), row(.88, [(5, 15), (50, 60), (105, 115)])
        disconnected = row(.7, [(200, 210)])
        for split in (False, True):
            for inputs in itertools.permutations([a, b, long, disconnected]):
                new, legacy = evaluate(inputs, split=split), evaluate(inputs, "legacy", split=split)
                for category in range(5):
                    self.assertEqual(candidates(new, category), candidates(legacy, category))
                self.assertEqual(candidates(new, 0), [id(a[1])])
                self.assertIn(id(long[1]), candidates(new, 1))

    def test_repair_chain_cannot_join_owners_even_if_each_child_touches_at_most_one(self):
        a, b = row(.8, [(0, 10)]), row(.7, [(100, 110)])
        repairs = [row(.99, [(5, 40)], repair=True), row(.98, [(35, 75)], repair=True),
                   row(.97, [(70, 105)], repair=True)]
        for inputs in itertools.permutations([a, b] + repairs):
            result = evaluate(inputs)
            self.assertEqual(candidates(result, 0), [id(a[1])])
            self.assertEqual(candidates(result, 2), [id(b[1])])
            self.assertEqual(set(candidates(result, 1)), {id(x[1]) for x in repairs})
            for entry in members(result, 1):
                evidence = competition_evidence(entry)
                self.assertEqual(evidence["decision"], "repair_scope_conflict")
                self.assertEqual(len(evidence["repair_components"][0]["direct_ordinary_component_ids"]), 2)

    def test_zero_owner_repair_group_competes_normally_and_single_owner_can_improve(self):
        core = row(.7, [(10, 20)])
        repair = row(.8, [(15, 25)], repair=True)
        orphan = row(.9, [(100, 120)], repair=True)
        orphan_alt = row(.85, [(110, 130)], repair=True)
        result = evaluate([core, repair, orphan, orphan_alt])
        self.assertEqual(set(candidates(result, 0)), {id(repair[1]), id(orphan[1])})
        self.assertEqual(set(candidates(result, 1)), {id(core[1]), id(orphan_alt[1])})
        kinds = {entry[3][METADATA_KEY]["scope_kind"] for entry in members(result, 0)}
        self.assertEqual(kinds, {"ordinary_with_repairs", "repair_only"})

    def test_high_score_multicore_repair_is_retained_alternative_even_with_possible_owners(self):
        a, b, repair = row(.7, [(10, 20)]), row(.6, [(40, 50)]), row(.99, [(15, 45)], repair=True)
        result = evaluate([repair, a, b])
        self.assertEqual(set(candidates(result, 2)), {id(a[1]), id(b[1])})
        self.assertEqual(candidates(result, 1), [id(repair[1])])
        self.assertEqual(members(result, 1)[0][0], .99)

    def test_dual_origin_is_ordinary_and_duplicate_repair_aliases_do_not_vote(self):
        a, b = row(.9, [(10, 20)]), row(.8, [(40, 50)])
        repair = row(.7, [(15, 45)], repair=True)
        duplicate = copy.deepcopy(repair)
        result = evaluate([a, b, repair, duplicate])
        self.assertEqual(set(candidates(result, 0)), {id(a[1]), id(b[1])})
        self.assertEqual(len(members(result, 1)), 2)
        natural = row(.7, [(15, 45)])
        inherit_candidate_origins(repair[1], (repair[1], natural[1]))
        self.assertEqual(candidate_origin_kind(repair[1]), "ordinary")
        result = evaluate([a, b, repair])
        self.assertEqual(candidates(result, 0), [id(a[1])])
        self.assertEqual(len(candidate_origin_evidence(repair[1])["repair_births"]), 1)

    def test_gates_and_full_inventory_are_preserved(self):
        a, b = row(.8, [(10, 20)]), row(.7, [(40, 50)])
        excluded = [row(.99, [(15, 45)], eligible=False),
                    row(.98, [(15, 45)], valid=False), row(.49, [(15, 45)])]
        result = evaluate(excluded + [a, b])
        self.assertEqual(set(candidates(result, 4)), {id(x[1]) for x in excluded})
        self.assertEqual(sum(len(x) for cat in result for x in cat.values()), 5)

    def test_unknown_and_stale_origin_fail_only_conservative_competition(self):
        candidate = row(.8, [(10, 20)])
        candidate[1].list_repeat_starts[0] += 1
        with self.assertRaisesRegex(ValueError, "binding is stale"):
            evaluate([candidate])
        self.assertEqual(candidates(evaluate([candidate], "legacy"), 0), [id(candidate[1])])
        del candidate[1].competition_origins
        with self.assertRaisesRegex(ValueError, "origin is missing"):
            evaluate([candidate])

    def test_disjoint_equal_score_representatives_in_one_ordinary_scope_are_explicit(self):
        a, b, connector = row(.8, [(10, 20)]), row(.8, [(40, 50)]), row(.6, [(15, 45)])
        observations = []
        for inputs in ([a, b, connector], [b, a, connector]):
            selected = members(evaluate(inputs), 0)[0]
            self.assertIs(selected[1], inputs[0][1])
            evidence = competition_evidence(selected)
            self.assertIs(evidence["unique_phase"], False)
            observations.append(evidence["tied_candidates"])
        self.assertEqual(observations[0], observations[1])

    def test_reclassification_refreshes_observation_without_mutating_inputs(self):
        a, b = row(.8, [(10, 20)]), row(.8, [(11, 21)])
        first = members(evaluate([a, b]), 0)[0]
        before = copy.deepcopy(first[3])
        updated = members(evaluate([first]), 0)[0]
        self.assertIsNone(competition_evidence(updated)["unique_phase"])
        self.assertEqual(first[3], before)
        corrupt = copy.deepcopy(updated)
        corrupt[0] = .9
        with self.assertRaisesRegex(ValueError, "competition"):
            competition_evidence(corrupt)


class InitialCompetitionReportingTests(unittest.TestCase):
    def test_report_and_sidecar_expose_same_tie_and_origin(self):
        from tests.test_candidate_evidence import _fixture, PARAMETERS
        from components.candidate_evidence import stage_candidate_evidence
        from components.output_transaction import OutputTransaction
        from components.report_model import ReportModelBuilder, ReportModelError
        source, original = _fixture()
        key = next(iter(original[0]))
        anchor, donor = original[0][key][0], original[1][key][0]
        donor[0] = anchor[0]
        params = dict(PARAMETERS, param_reconstruction_policy="conservative")
        categories = ArrayEvaluation({}, params, precomputed_scored_rows={key: [anchor, donor]}).output()
        expected = competition_evidence(members(categories, 0)[0])
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "source.fna"
            fasta.write_text(">source\n" + source + "\n")
            report = ReportModelBuilder(fasta, categories, {}, [], parameters=params).build()
            self.assertEqual(report["arrays"][0][METADATA_KEY], expected)
            path = Path(temporary) / "candidate_evidence.jsonl.gz"
            with OutputTransaction() as publication:
                stage_candidate_evidence(publication, path, source_header="source", source_sequence=source,
                    categories=categories, parameters=params, feature_names=[], canonical_report=report)
                publication.publish()
            records = [json.loads(line) for line in gzip.decompress(path.read_bytes()).splitlines()]
            accepted = [x for x in records if x["record_type"] == "candidate" and x["accepted"]][0]
            self.assertEqual(accepted[METADATA_KEY], expected)
            self.assertEqual(accepted["candidate_origins"]["kind"], "ordinary")
            members(categories, 0)[0][0] = .91
            with self.assertRaisesRegex(ReportModelError, "initial_competition"):
                ReportModelBuilder(fasta, categories, {}, [], parameters=params).build()


if __name__ == "__main__":
    unittest.main()
