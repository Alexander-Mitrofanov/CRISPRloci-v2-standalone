"""Low-score output competition must preserve rescue and audit inventories."""

import copy
import csv
import gzip
import itertools
import json
import tempfile
import unittest
from pathlib import Path

from components.evidence.initial_competition import select_low_score_representatives
from components.pipeline import Pipeline
from components.reporting.html_report import render_html_report
from components.reporting.report_extensions import ReportExtensionError, attach_report_extensions
from components.reporting.report_model import ReportModelBuilder
from components.reporting.stage import OutputMaker
from components.rescue.likely_tandem_transaction import apply_likely_tandem_quarantine
from tests.test_candidate_evidence import PARAMETERS, _fixture
from tests.test_initial_competition import evaluate, members, row


class LowScoreSelectionTests(unittest.TestCase):
    def test_matches_accepted_selection_for_scores_ties_and_disjoint_loci(self):
        for policy, split in itertools.product(("legacy", "conservative"), (False, True)):
            parameters = {"param_reconstruction_policy": policy,
                          "param_split_competitive_components": split}
            for inputs in itertools.permutations([
                    row(.3, [(10, 20)]), row(.4, [(11, 21)]),
                    row(.4, [(12, 22)]), row(.2, [(40, 50)])]):
                with self.subTest(policy=policy, split=split):
                    low = list(inputs)
                    categories = evaluate(low, policy, split)
                    self.assertEqual(len(members(categories, 4)), len(low))
                    selected = select_low_score_representatives(categories[4], parameters)
                    selected_candidates = [id(entry[1]) for entries in selected.values() for entry in entries]
                    for offset, accepted_index in ((.31, 2), (.55, 0)):
                        accepted = [[entry[0] + offset, entry[1], entry[2]] for entry in low]
                        expected = members(evaluate(accepted, policy, split), accepted_index)
                        self.assertEqual(selected_candidates, [id(entry[1]) for entry in expected])

    def test_conservative_repair_bridge_cannot_merge_low_score_loci(self):
        left, right = row(.3, [(10, 20)]), row(.2, [(40, 50)])
        bridge = row(.49, [(15, 45)], repair=True)
        category = {(0, 100): [bridge, left, right]}
        parameters = {"param_reconstruction_policy": "conservative",
                      "param_split_competitive_components": True}
        selected = select_low_score_representatives(category, parameters)
        self.assertEqual(selected[(0, 100)], [left, right])
        self.assertEqual(category[(0, 100)], [bridge, left, right])

    def test_filter_failures_remain_low_and_original_metadata_is_preserved(self):
        failed = row(.98, [(10, 20)], valid=False)
        lineage = row(.99, [(11, 21)], eligible=False)
        failed.append({"classification_basis": {"original": True}})
        categories = evaluate([row(.3, [(12, 22)]), failed, lineage])
        pipeline = Pipeline.__new__(Pipeline)
        pipeline.categories = categories
        pipeline.parameters = {"param_accepted_only": False}
        inventory = list(members(categories, 4))
        pipeline._drop_nonaccepted_categories()
        self.assertEqual(members(pipeline.categories, 4), [lineage])
        self.assertEqual(len(inventory), 3)
        self.assertFalse(members(pipeline.categories, 0))
        self.assertFalse(members(pipeline.categories, 2))
        pipeline.parameters["param_accepted_only"] = True
        pipeline._drop_nonaccepted_categories()
        self.assertEqual(pipeline.categories, [{}, {}, {}, {}, {}])

    def test_invalid_scores_or_stale_origins_fail_closed(self):
        for score in (float("nan"), float("inf"), -.1, 1.1, True, "0.4"):
            with self.subTest(score=score), self.assertRaises(ValueError):
                select_low_score_representatives({(0, 100): [row(score, [(10, 20)])]}, {})
        stale = row(.4, [(10, 20)])
        stale[1].list_repeat_starts[0] += 1
        with self.assertRaisesRegex(ValueError, "binding is stale"):
            select_low_score_representatives({(0, 100): [stale]},
                                            {"param_reconstruction_policy": "conservative"})


class LowScorePublicationTests(unittest.TestCase):
    def test_all_outputs_select_one_while_evidence_retains_every_candidate(self):
        source, categories = _fixture()
        key = next(iter(categories[4]))
        winner = copy.deepcopy(categories[4][key][0])
        winner[0] = .49
        categories[4][key].append(winner)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "source.fa"
            fasta.write_text(">source\n" + source + "\n")
            OutputMaker(
                file_path=fasta, parameters=PARAMETERS,
                flags={"output_profile": "standard", "flag_html_report": True,
                       "flag_fasta_report": True},
                result_path=str(root), pickle_result_path="", json_result_path="",
                categories=categories, non_array_data={"Strand": {
                    "Bona-fide": ["Forward"], "Possible": ["Forward"],
                }, "Leader": {0: [""]}, "Downstream": {0: [""]}, "IS": {}},
                list_features=[], header=">source\n", candidate_evidence_categories=categories,
            )
            report = json.loads((root / "report.json").read_text())
            low = [array for array in report["arrays"] if array["category"] == "Low score"]
            self.assertEqual([.49], [array["certainty_score"] for array in low])
            self.assertEqual(2, report["summary"]["accepted_array_count"])
            self.assertEqual(1, (root / "Low_Score_Candidates.txt").read_text().count("Low score CRISPR:"))
            with (root / "Summary.csv").open() as handle:
                summary = list(csv.DictReader(handle))
            self.assertEqual(2, len(summary))
            records = [json.loads(line) for line in gzip.decompress(
                (root / "candidate_evidence.jsonl.gz").read_bytes()).splitlines()]
            self.assertEqual(2, records[-1]["category_counts"]["Low score"])
            self.assertEqual(2, len(categories[4][key]))
            render_html_report(report)
            from tests.test_batch_runner import batch
            batch.validate_report(root / "report.json", {"source": source}, None, region_screen=False)

    def test_tandem_ledger_retains_demotions_omitted_by_selection(self):
        from components.refinement.candidates import CrisprCandidate
        repeat, spacer = "ACGTAC", "A" * 20
        source = repeat + spacer + repeat + spacer + repeat
        candidate = CrisprCandidate([repeat] * 3, [repeat] * 3, [spacer] * 2, [0, 26, 52])
        key = (0, len(source))
        categories = [{}, {}, {key: [[.55, candidate, []]]}, {},
                      {key: [[.9, copy.deepcopy(candidate), []]]}]
        projected, records, stats = apply_likely_tandem_quarantine(categories)
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "source.fa"
            fasta.write_text(">source\n" + source + "\n")
            core = ReportModelBuilder(fasta, projected, {}, [],
                                      parameters={"param_accepted_only": False},
                                      provenance_categories=categories).build()
            kwargs = dict(likely_tandem_records=records, likely_tandem_stats=stats,
                          likely_tandem_source_categories=categories, source_sequence=source)
            report = attach_report_extensions(core, **kwargs)
            ledger = report["extensions"]["likely_tandem"]
            self.assertEqual("omitted_by_low_score_selection", ledger["arrays"][0]["core_projection"])
            self.assertIsNone(ledger["arrays"][0]["reported_array_id"])
            self.assertEqual(2, ledger["low_score_projection"]["rows_before"])
            self.assertEqual(1, ledger["low_score_projection"]["rows_after"])
            render_html_report(report)
            from tests.test_batch_runner import batch
            report_path = Path(temporary) / "report.json"
            report_path.write_text(json.dumps(report))
            batch.validate_report(report_path, {"source": source}, None, region_screen=False)
            malformed = copy.deepcopy(report)
            malformed["extensions"]["likely_tandem"]["low_score_projection"]["rows_after"] = 2
            with self.assertRaisesRegex(ValueError, "Low score projection"):
                render_html_report(malformed)
            report_path.write_text(json.dumps(malformed))
            with self.assertRaisesRegex(batch.CompletionError, "Low score projection"):
                batch.validate_report(report_path, {"source": source}, None, region_screen=False)
            corrupted = copy.deepcopy(core)
            corrupted["arrays"][0]["certainty_score"] = .89
            with self.assertRaisesRegex(ReportExtensionError, "Low score selection does not replay"):
                attach_report_extensions(corrupted, **kwargs)


if __name__ == "__main__":
    unittest.main()
