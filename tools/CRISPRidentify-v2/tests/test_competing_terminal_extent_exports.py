"""Public terminal disclosure preserves accepted science and output lifecycle."""

from contextlib import nullcontext
import copy
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from components import module_output_maker, output_transaction
from components.components_detection_refinement import CrisprCandidate
from components.competing_terminal_extent_outputs import CSV_FILENAME
from components.competing_terminal_extents import competing_terminal_extent_summaries
from components.html_report import CanonicalReportMaker, render_html_report
from components.likely_tandem_transaction import apply_likely_tandem_quarantine
from tests.test_candidate_evidence import PARAMETERS, _fixture


def disclosure_fixture():
    source, originals = _fixture()
    repeat = next(iter(originals[0].values()))[0][1].list_repeats[0]
    starts = [i * (len(repeat) + 20) for i in range(8)]

    def candidate(indexes):
        positions = [starts[i] for i in indexes]
        return CrisprCandidate(
            list_repeats=[repeat] * len(positions),
            list_repeats_gaped=[repeat] * len(positions),
            list_repeat_starts=positions,
            list_spacers=[source[a + len(repeat):b]
                          for a, b in zip(positions, positions[1:])],
        )

    key = (0, len(source))
    categories = [{}, {}, {key: [[.61, candidate([2, 3, 4]), []]]},
                  {key: [[.55, candidate([1, 2, 3, 4, 5]), []]]}, {}]
    return source, categories


class CompetingTerminalExtentExportTests(unittest.TestCase):
    def publish(self, root, *, evidence=False, accepted_only=False,
                machine=False, baseline=False, categories=None):
        root.mkdir(exist_ok=True)
        source, default_categories = disclosure_fixture()
        categories = default_categories if categories is None else categories
        fasta = root / "source.fa"
        fasta.write_text(">source\n" + source + "\n")
        parameters = dict(PARAMETERS, param_reconstruction_policy="conservative",
                          param_accepted_only=accepted_only)
        projected, records, stats = apply_likely_tandem_quarantine(categories)
        if accepted_only:
            projected = [projected[i] if i in (0, 2) else {} for i in range(5)]
        annotations = {
            "Strand": {"Bona-fide": [], "Possible": ["Forward"], "Low score": []},
            "Leader": {0: []}, "Downstream": {0: []}, "IS": {},
        }
        attach = module_output_maker.attach_report_extensions

        def baseline_extensions(*args, **kwargs):
            kwargs.pop("competing_terminal_extent_context", None)
            return attach(*args, **kwargs)

        patch = mock.patch.object(module_output_maker, "attach_report_extensions",
                                  side_effect=baseline_extensions) if baseline else nullcontext()
        with patch:
            module_output_maker.OutputMaker(
                file_path=fasta, parameters=parameters,
                flags={"flag_fasta_report": True, "flag_html_report": True,
                       "output_profile": "machine" if machine else "standard"},
                result_path=str(root), pickle_result_path=str(root / "pickle"),
                json_result_path=str(root / "legacy-json"), categories=projected,
                non_array_data=annotations, list_features=[], header=">source\n",
                candidate_evidence_categories=categories if evidence else None,
                likely_tandem_records=records, likely_tandem_stats=stats,
                likely_tandem_source_categories=categories,
            )
        return json.loads((root / "report.json").read_text())

    def test_disclosure_changes_no_accepted_or_retained_scientific_content(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            old = self.publish(root / "old", baseline=True, evidence=True)
            new = self.publish(root / "new", evidence=True)
            self.assertEqual(old["arrays"], new["arrays"])
            self.assertEqual(old["summary"], new["summary"])
            self.assertEqual(1, new["summary"]["accepted_array_count"])
            for filename in ("Summary.csv", "Arrays.fasta", "Repeats.fasta", "Spacers.fasta",
                             "legacy-json/source.json", "pickle/source.pkl",
                             "candidate_evidence.jsonl.gz"):
                with self.subTest(filename=filename):
                    self.assertEqual((root / "old" / filename).read_bytes(),
                                     (root / "new" / filename).read_bytes())
            self.assertFalse((root / "old" / CSV_FILENAME).exists())
            summaries = competing_terminal_extent_summaries(new)
            item = next(iter(summaries.values()))[0]
            self.assertEqual((3, 5), (item["accepted_repeat_count"], item["repeat_count"]))
            self.assertEqual((.61, .55), (item["accepted_own_score"], item["own_score"]))
            text = (root / "new/Possible_Candidates.txt").read_text()
            html = (root / "new/report.html").read_text()
            for view in (text, html):
                self.assertIn("Unaccepted competing terminal extent", view)
                self.assertIn("not_recorded", view)
                self.assertIn(item["observation_id"], view)
            with (root / "new" / CSV_FILENAME).open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(1, len(rows))
            self.assertEqual("3", rows[0]["Accepted repeat count"])
            self.assertEqual("5", rows[0]["Observed repeat count"])
            self.assertEqual("Possible Discarded", rows[0]["Observed category"])
            sidecar = [json.loads(line) for line in gzip.decompress(
                (root / "new/candidate_evidence.jsonl.gz").read_bytes()).splitlines()]
            self.assertTrue(any(row.get("observation_id") == item["observation_id"]
                                for row in sidecar))

    def test_evidence_off_and_accepted_only_keep_same_disclosure(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reports = []
            for accepted_only in (False, True):
                for evidence in (False, True):
                    folder = root / (str(accepted_only) + str(evidence))
                    report = self.publish(folder, evidence=evidence, accepted_only=accepted_only)
                    reports.append(report)
                    self.assertEqual(evidence, (folder / "candidate_evidence.jsonl.gz").exists())
            for report in reports[1:]:
                self.assertEqual(reports[0]["arrays"], report["arrays"])
                self.assertEqual(reports[0]["extensions"]["competing_terminal_extents"],
                                 report["extensions"]["competing_terminal_extents"])

    def test_reused_folder_removes_stale_companion_and_notes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.publish(root)
            self.assertTrue((root / CSV_FILENAME).exists())
            _, categories = disclosure_fixture()
            categories[3] = {}
            self.publish(root, categories=categories)
            self.assertFalse((root / CSV_FILENAME).exists())
            self.assertNotIn("Unaccepted competing terminal extent",
                             (root / "Possible_Candidates.txt").read_text())
            self.publish(root)
            machine = self.publish(root, machine=True)
            self.assertTrue(competing_terminal_extent_summaries(machine))
            self.assertFalse((root / CSV_FILENAME).exists())
            self.assertFalse((root / "Possible_Candidates.txt").exists())

    def test_recorded_negative_stays_visible_as_blocked_in_every_view(self):
        _, categories = disclosure_fixture()
        next(iter(categories[3].values()))[0].append({"tags": ["tandem"]})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = self.publish(root, categories=categories)
            item = next(iter(competing_terminal_extent_summaries(report).values()))[0]
            self.assertEqual("blocked_observation", item["status"])
            self.assertIn("candidate_recorded_negative", item["block_reasons"])
            self.assertEqual(3, len(report["arrays"][0]["repeats"]))
            for filename in ("report.html", "Possible_Candidates.txt", CSV_FILENAME):
                with self.subTest(filename=filename):
                    text = (root / filename).read_text()
                    self.assertIn("blocked_observation", text)
                    self.assertIn("candidate_recorded_negative", text)
                    self.assertIn(item["observation_id"], text)
            self.assertIn(item["repeats"][-1], (root / "report.html").read_text())

    def test_companion_failure_cannot_publish_completed_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.publish(root)
            old_csv = (root / CSV_FILENAME).read_bytes()
            replace = output_transaction.os.replace

            def fail_csv(source, target):
                if (Path(source).name.startswith(".crispridentify-stage-")
                        and Path(target).name == CSV_FILENAME):
                    raise OSError("injected companion failure")
                return replace(source, target)

            with mock.patch.object(output_transaction.os, "replace", side_effect=fail_csv):
                with self.assertRaisesRegex(OSError, "injected companion failure"):
                    self.publish(root)
            self.assertFalse((root / "report.json").exists())
            self.assertEqual(old_csv, (root / CSV_FILENAME).read_bytes())
            self.assertFalse(list(root.glob(".crispridentify-*")))

    def test_old_reports_render_and_stale_anchor_cannot_hide_in_html(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = self.publish(Path(temporary))
        old = copy.deepcopy(report)
        old["extensions"].pop("competing_terminal_extents")
        self.assertNotIn("Unaccepted competing terminal extent", render_html_report(old))
        altered = copy.deepcopy(report)
        altered["arrays"][0]["certainty_score"] = .62
        with self.assertRaises(ValueError):
            render_html_report(altered)

    def test_standalone_report_writer_includes_available_inventory(self):
        source, categories = disclosure_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "source.fa"
            fasta.write_text(">source\n" + source + "\n")
            (root / CSV_FILENAME).write_text("stale companion\n")
            before = copy.deepcopy(categories)
            maker = CanonicalReportMaker(
                file_path=fasta, result_path=root, categories=categories,
                non_array_data={}, list_feature_names=[],
                parameters=dict(PARAMETERS, param_reconstruction_policy="conservative"),
                write_html=True,
            )
            self.assertTrue(competing_terminal_extent_summaries(maker.model))
            self.assertFalse((root / CSV_FILENAME).exists())
            for index in (2, 3):
                original = next(iter(before[index].values()))[0]
                after = next(iter(categories[index].values()))[0]
                self.assertEqual(original[:1] + original[2:], after[:1] + after[2:])
                self.assertEqual(original[1].__dict__, after[1].__dict__)


if __name__ == "__main__":
    unittest.main()
