"""Public member views preserve coupled uncertainty with evidence retention off."""

import copy
import csv
import json
from pathlib import Path
import pickle
import tempfile
import unittest

from components.components_output_maker import (
    CompleteFolderSummaryMaker, FastaOutputArrayMaker, SummaryMakerCSV,
    _iter_fasta_records,
)
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.html_report import render_html_report
from components.module_output_maker import OutputMaker
from components.reconstruction_boundary_family import boundary_family_summary
from components.reconstruction_transaction import reconstruct_accepted_arrays
from components.sequence_utils import reverse_complement
from tests.test_candidate_evidence import PARAMETERS
from tests.test_reconstruction_boundary_family import fixture


def annotations(strand="Forward"):
    return {"Strand": {"Bona-fide": [strand], "Possible": [], "Low score": []},
            "Leader": {0: [""]}, "Downstream": {0: [""]}, "IS": {}}


class ExtentFamilyExportTests(unittest.TestCase):
    def setUp(self):
        self.source, original, scorer = fixture()
        self.categories, self.audit = reconstruct_accepted_arrays(
            self.source, original, PARAMETERS, scorer, lambda: ())
        self.row = next(iter(self.categories[0].values()))[0]
        self.family = boundary_family_summary(self.row)
        self.assertIsNotNone(self.family)

    def publish(self, root, strand="Forward"):
        fasta = root / "source.fa"
        fasta.write_text(">source\n" + self.source + "\n")
        output = root / "output"
        output.mkdir()
        OutputMaker(
            file_path=fasta, parameters=PARAMETERS,
            flags={"flag_fasta_report": True, "flag_html_report": True,
                   "flag_strand": True, "output_profile": "standard"},
            result_path=str(output), pickle_result_path=str(root / "pickle"),
            json_result_path=str(root / "legacy-json"), categories=self.categories,
            non_array_data=annotations(strand), list_features=[COMPACT_FEATURE_SCHEMA],
            header=">source\n", reconstruction_evidence=self.audit,
            candidate_evidence_categories=None,
        )
        return output

    def test_every_export_labels_one_actual_member_without_optional_sidecar(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = self.publish(root)
            report = json.loads((output / "report.json").read_text())
            array = report["arrays"][0]
            self.assertEqual(self.family, array["boundary_family"])
            self.assertEqual(self.family["member_count"],
                             len(array["classification_basis"]["boundary_family"]["members"]))
            self.assertFalse((output / "candidate_evidence.jsonl.gz").exists())
            text = (output / "Bona-Fide_Candidates.txt").read_text()
            self.assertIn("Boundary uncertainty:", text)
            self.assertIn(self.family["display_geometry_id"], text)
            html = (output / "report.html").read_text()
            self.assertIn("Shared extent, uncertain repeat ends", html)
            for member in self.family["admissible_geometry_ids"]:
                self.assertIn(member, html)
            with (output / "Summary.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(1, len(rows))
            self.assertEqual("boundary_uncertain", rows[0]["Geometry status"])
            self.assertEqual(self.family["display_geometry_id"], rows[0]["Display geometry ID"])
            self.assertEqual(self.row[1].list_repeat_starts[0] + 1, int(rows[0]["Start"]))
            self.assertEqual(self.family["supported_repeat_count"],
                             int(rows[0]["Supported repeat count"]))
            legacy = json.loads((root / "legacy-json/source.json").read_text())
            view = next(iter(legacy["Arrays"]["Bona-fide"].values()))
            self.assertEqual(self.family, view["boundary_family"])
            self.assertEqual(self.row[1].list_repeats, view["list_repeats"])
            with (root / "pickle/source.pkl").open("rb") as handle:
                restored = pickle.load(handle)
            self.assertEqual(self.family,
                             boundary_family_summary(next(iter(restored[0].values()))[0]))
            for name, sequences in (
                    ("Repeats.fasta", self.row[1].list_repeats),
                    ("Spacers.fasta", self.row[1].list_spacers),
                    ("Arrays.fasta", [array["source_sequence"]])):
                records = list(_iter_fasta_records(output / name))
                self.assertEqual(sequences, [sequence for _, sequence in records])
                for header, _ in records:
                    self.assertIn("boundary_uncertain=true", header)
                    self.assertIn("member=" + self.family["display_geometry_id"], header)

    def test_reverse_fasta_keeps_genomic_member_identity_and_exact_sequences(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            FastaOutputArrayMaker(str(output), self.categories, annotations("Reverse"))
            start = self.row[1].list_repeat_starts[0]
            end = self.row[1].list_repeat_starts[-1] + len(self.row[1].list_repeats[-1])
            for name, sequences in (
                    ("Repeats.fasta", self.row[1].list_repeats),
                    ("Spacers.fasta", self.row[1].list_spacers),
                    ("Arrays.fasta", [self.source[start:end]])):
                records = list(_iter_fasta_records(output / name))
                self.assertEqual([reverse_complement(r) for r in reversed(sequences)],
                                 [sequence for _, sequence in records])
                for header, _ in records:
                    self.assertIn(self.family["display_geometry_id"], header)

    def test_mixed_family_ordinary_empty_csvs_aggregate_with_one_schema(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ordinary = [dict(self.categories[0]), {}, {}, {}, {}]
            key = next(iter(ordinary[0]))
            ordinary[0][key] = [self.row[:3]]
            for name, categories in (("family", self.categories), ("ordinary", ordinary),
                                     ("empty", [{}, {}, {}, {}, {}])):
                output = root / name
                output.mkdir()
                SummaryMakerCSV(str(output), categories, annotations())
            CompleteFolderSummaryMaker(str(root), ["family", "ordinary", "empty"])
            with (root / "Complete_summary.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(2, len(rows))
            self.assertEqual({"boundary_uncertain", "reported_geometry"},
                             {r["Geometry status"] for r in rows})
            self.assertTrue(all(None not in row for row in rows))

    def test_html_rejects_missing_or_mismatched_uncertainty_disclosure(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = self.publish(Path(temporary))
            report = json.loads((output / "report.json").read_text())
        for change in ("missing", "member", "score"):
            altered = copy.deepcopy(report)
            array = altered["arrays"][0]
            if change == "missing":
                del array["boundary_family"]
            elif change == "member":
                array["boundary_family"]["display_geometry_id"] = "0" * 64
            else:
                array["certainty_score"] = .123
            with self.subTest(change=change), self.assertRaises(ValueError):
                render_html_report(altered)


if __name__ == "__main__":
    unittest.main()
