"""Possible checked extents preserve confidence and publication provenance.

The shared synthetic fixture has real model outputs, not independently extracted
features or biological labels. These tests exercise output lifecycle, not native
prediction or the biological accuracy of the reconstructed path.
"""

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
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.html_report import CanonicalReportMaker, render_html_report
from components.likely_tandem_transaction import apply_likely_tandem_quarantine
from components.possible_reconstruction import (
    EVIDENCE_FILENAME, METHOD, reconstruct_possible_arrays,
)
from components.report_model import ReportModelBuilder, ReportModelError
from tests.test_candidate_evidence import PARAMETERS
from tests.test_possible_reconstruction import fixture
from tests.test_reconstruction_transaction import fixture as bf_fixture


def publication_fixture(*, tandem_sibling=False):
    """One checked Possible extent, one disjoint unchanged BF control."""
    source, categories, scorer = fixture()
    core = next(iter(categories[2].values()))[0]
    motif = core[1].list_repeats[0]

    def append_region(spacers):
        nonlocal source
        source += "ACGT" * 25
        starts = [len(source)]
        source += motif
        for spacer in spacers:
            source += spacer
            starts.append(len(source))
            source += motif
        candidate = CrisprCandidate(
            list_repeats=[motif] * len(starts),
            list_repeats_gaped=[motif] * len(starts),
            list_spacers=spacers, list_repeat_starts=starts,
        )
        return (starts[0], len(source)), candidate

    key, bf = append_region(list(core[1].list_spacers))
    _, original_bf_categories, _ = bf_fixture()
    original_bf = next(iter(original_bf_categories[0].values()))[0]
    categories[0][key] = [[original_bf[0], bf, copy.deepcopy(original_bf[2])]]
    if tandem_sibling:
        key, tandem = append_region(["A" * 31, "A" * 31])
        categories[2][key] = scorer({key: [tandem]})[key]
    output, audit = reconstruct_possible_arrays(
        source, categories, PARAMETERS, scorer, lambda: ())
    return source, categories, output, audit


class PossibleReconstructionExportTests(unittest.TestCase):
    def publish(self, root, state, *, evidence=False, accepted_only=False,
                audit_override=False, audit=None, omit_selected=False):
        root.mkdir(exist_ok=True)
        source, _before, provenance, recorded_audit = state
        fasta = root / "source.fa"
        fasta.write_text(">source\n" + source + "\n", encoding="ascii")
        projected, records, stats = apply_likely_tandem_quarantine(provenance)
        if accepted_only:
            projected = [projected[i] if i in (0, 2) else {} for i in range(5)]
        if omit_selected:
            projected[2] = {}
        annotations = {
            "Strand": {name: ["Forward"] * sum(map(len, projected[index].values()))
                       for index, name in ((0, "Bona-fide"), (2, "Possible"), (4, "Low score"))},
            "Leader": {0: [""] * sum(map(len, projected[0].values()))},
            "Downstream": {0: [""] * sum(map(len, projected[0].values()))},
            "IS": {},
        }
        module_output_maker.OutputMaker(
            file_path=fasta, result_path=str(root), categories=projected,
            parameters=dict(PARAMETERS, param_reconstruction_policy="conservative",
                            param_accepted_only=accepted_only),
            flags={"flag_fasta_report": True, "flag_html_report": True},
            pickle_result_path=str(root / "pickle"),
            json_result_path=str(root / "legacy-json"), non_array_data=annotations,
            list_features=[COMPACT_FEATURE_SCHEMA], header=">source\n",
            candidate_evidence_categories=provenance if evidence else None,
            likely_tandem_source_categories=provenance,
            likely_tandem_records=records, likely_tandem_stats=stats,
            possible_reconstruction_evidence=audit if audit_override else recorded_audit,
        )
        return json.loads((root / "report.json").read_text())

    def test_actual_outputs_preserve_possible_confidence_own_score_and_bf(self):
        state = publication_fixture()
        source, before, after, audit = state
        original = next(iter(before[2].values()))[0]
        chosen = next(iter(after[2].values()))[0]
        self.assertEqual(1, audit["stats"]["selections_committed"])
        self.assertNotEqual(original[0], chosen[0])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = self.publish(root, state, evidence=True)
            baseline = ReportModelBuilder(
                root / "source.fa", before, {}, [COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS, param_reconstruction_policy="conservative"),
            ).build()
            bf = [row for row in report["arrays"] if row["category"] == "Bona-fide"]
            self.assertEqual(bf, [row for row in baseline["arrays"] if row["category"] == "Bona-fide"])
            possible = next(row for row in report["arrays"] if row["category"] == "Possible")
            summary = possible["possible_reconstruction"]
            self.assertEqual(2, report["summary"]["accepted_array_count"])
            self.assertEqual((3, 6, 3), (summary["core_repeat_count"], summary["repeat_count"],
                                         summary["added_repeat_count"]))
            self.assertEqual("Possible", summary["confidence_category"])
            self.assertEqual(METHOD, possible["classification_basis"]["method"])
            self.assertEqual(original[0], summary["core_original_own_score"])
            self.assertEqual(chosen[0], possible["certainty_score"])
            self.assertEqual(chosen[0], summary["own_score"])
            for name in ("report.html", "Possible_Candidates.txt"):
                text = (root / name).read_text()
                self.assertIn("Checked extent of a Possible core", text)
                self.assertIn("remains Possible", text)
                self.assertIn("not calibrated probabilities", text)
            with (root / "Summary.csv").open(newline="") as handle:
                csv_rows = list(csv.reader(handle))
            possible_csv = next(row for row in csv_rows[1:] if "Possible" in row)
            self.assertIn(str(chosen[0]), possible_csv)
            self.assertNotIn(str(original[0]), possible_csv)
            legacy = json.loads((root / "legacy-json/source.json").read_text())
            self.assertEqual(summary, next(iter(legacy["Arrays"]["Possible"].values()))[
                "possible_reconstruction"])
            sidecar = [json.loads(line) for line in gzip.decompress(
                (root / "candidate_evidence.jsonl.gz").read_bytes()).splitlines()]
            current = next(row for row in sidecar if row.get("possible_reconstruction"))
            self.assertEqual(summary, current["possible_reconstruction"])
            self.assertTrue(any(row.get("category") == "Possible Discarded"
                                and row.get("repeat_count") == 3
                                and row.get("score") == original[0] for row in sidecar))
            self.assertEqual(audit, json.loads((root / EVIDENCE_FILENAME).read_text()))
            self.assertFalse(audit["biological_accuracy_established"])
            self.assertEqual(source, "".join((root / "source.fa").read_text().splitlines()[1:]))

    def test_retention_off_and_accepted_only_have_identical_accepted_proofs(self):
        state = publication_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reports = []
            for accepted_only in (False, True):
                for evidence in (False, True):
                    folder = root / f"{accepted_only}-{evidence}"
                    reports.append(self.publish(folder, state, evidence=evidence,
                                                accepted_only=accepted_only))
                    self.assertEqual(evidence, (folder / "candidate_evidence.jsonl.gz").exists())
                    self.assertTrue((folder / EVIDENCE_FILENAME).exists())
            for report in reports[1:]:
                self.assertEqual(reports[0]["arrays"], report["arrays"])
                self.assertEqual(reports[0]["summary"], report["summary"])
            # Retention-off publication still needs complete original-core provenance.
            source, _before, after, _audit = state
            fasta = root / "source.fa"
            fasta.write_text(">source\n" + source + "\n")
            with self.assertRaisesRegex(ReportModelError, "original core"):
                ReportModelBuilder(
                    fasta, [after[0], {}, after[2], {}, {}], {}, [COMPACT_FEATURE_SCHEMA],
                    parameters=dict(PARAMETERS, param_accepted_only=True),
                ).build()

    def test_real_unrelated_quarantine_does_not_hide_or_upgrade_selected_possible(self):
        state = publication_fixture(tandem_sibling=True)
        final, records, stats = apply_likely_tandem_quarantine(state[2])
        self.assertEqual(1, stats.possible_quarantined)
        self.assertEqual(1, sum(record.action == "quarantined" for record in records))
        self.assertEqual(1, sum(map(len, final[2].values())))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = self.publish(root, state, evidence=True, accepted_only=True)
            self.assertEqual(2, report["summary"]["accepted_array_count"])
            selected = [row for row in report["arrays"] if "possible_reconstruction" in row]
            self.assertEqual(1, len(selected))
            self.assertEqual("Possible", selected[0]["category"])
            self.assertEqual(state[3], json.loads((root / EVIDENCE_FILENAME).read_text()))

    def test_missing_selected_output_cannot_be_excused_by_retained_nonquarantine(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(ValueError):
                self.publish(root, publication_fixture(), accepted_only=True, omit_selected=True)
            self.assertFalse((root / "report.json").exists())
            self.assertFalse((root / "Summary.csv").exists())

    def test_corrupt_or_missing_audit_fails_before_legacy_publication(self):
        state = publication_fixture()
        audits = [None]
        for mutation in ("source", "count", "selection_score"):
            broken = copy.deepcopy(state[3])
            if mutation == "source":
                broken["source_sha256"] = "0" * 64
            elif mutation == "count":
                broken["stats"]["selections_committed"] = 0
            else:
                next(row for row in broken["decisions"] if row["action"] == "selected")[
                    "reconstruction_score"] += .01
            audits.append(broken)
        for index, audit in enumerate(audits):
            with self.subTest(index=index), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                self.publish(root, state, evidence=True)
                old_files = {name: (root / name).read_bytes() for name in (
                    "Summary.csv", "Possible_Candidates.txt", "candidate_evidence.jsonl.gz",
                    EVIDENCE_FILENAME)}
                with self.assertRaises(ValueError):
                    self.publish(root, state, evidence=True, audit_override=True, audit=audit)
                self.assertFalse((root / "report.json").exists())
                for name, value in old_files.items():
                    self.assertEqual(value, (root / name).read_bytes())
                self.assertFalse(list(root.glob(".crispridentify-*")))

    def test_corrupt_summary_and_wrong_confidence_cannot_render(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = self.publish(Path(temporary), publication_fixture())
        for field, value in (("confidence_category", "Bona-fide"), ("own_score", .99),
                             ("core_repeat_count", 6)):
            broken = copy.deepcopy(report)
            array = next(row for row in broken["arrays"] if "possible_reconstruction" in row)
            array["possible_reconstruction"][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                render_html_report(broken)

    def test_corrupt_own_score_or_missing_original_core_aborts_before_outputs(self):
        for mutation in ("score", "original_core"):
            state = publication_fixture()
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                self.publish(root, state)
                preserved_text = (root / "Possible_Candidates.txt").read_bytes()
                broken = copy.deepcopy(state)
                if mutation == "score":
                    next(iter(broken[2][2].values()))[0][0] += .01
                else:
                    broken[2][3].clear()
                with self.assertRaises(ReportModelError):
                    self.publish(root, broken)
                self.assertFalse((root / "report.json").exists())
                self.assertEqual(preserved_text, (root / "Possible_Candidates.txt").read_bytes())
                self.assertFalse(list(root.glob(".crispridentify-*")))

    def test_reused_folder_and_standalone_writer_retire_stale_possible_audit(self):
        state = publication_fixture()
        source, original, output, _audit = state
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.publish(root, state, evidence=True)
            plain = (source, original, original, None)
            report = self.publish(root, plain)
            self.assertFalse((root / EVIDENCE_FILENAME).exists())
            self.assertFalse((root / "candidate_evidence.jsonl.gz").exists())
            self.assertFalse(any("possible_reconstruction" in row for row in report["arrays"]))
            self.assertNotIn("Checked extent of a Possible core", (root / "report.html").read_text())
            self.publish(root, state)
            final, records, stats = apply_likely_tandem_quarantine(output)
            maker = CanonicalReportMaker(
                file_path=root / "source.fa", result_path=root, categories=final,
                non_array_data={}, list_feature_names=[COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS, param_reconstruction_policy="conservative"),
                write_html=True, likely_tandem_source_categories=output,
                likely_tandem_records=records, likely_tandem_stats=stats,
            )
            self.assertTrue(any("possible_reconstruction" in row for row in maker.model["arrays"]))
            self.assertFalse((root / EVIDENCE_FILENAME).exists())
            self.assertIn("Checked extent of a Possible core", (root / "report.html").read_text())

    def test_audit_publication_failure_never_leaves_completed_report(self):
        state = publication_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.publish(root, state)
            old_audit = (root / EVIDENCE_FILENAME).read_bytes()
            replace = output_transaction.os.replace

            def fail_audit(source, target):
                if (Path(source).name.startswith(".crispridentify-stage-")
                        and Path(target).name == EVIDENCE_FILENAME):
                    raise OSError("injected Possible audit publication failure")
                return replace(source, target)

            with mock.patch.object(output_transaction.os, "replace", side_effect=fail_audit):
                with self.assertRaisesRegex(OSError, "injected Possible audit"):
                    self.publish(root, state)
            self.assertFalse((root / "report.json").exists())
            self.assertEqual(old_audit, (root / EVIDENCE_FILENAME).read_bytes())
            self.assertFalse(list(root.glob(".crispridentify-*")))


if __name__ == "__main__":
    unittest.main()
