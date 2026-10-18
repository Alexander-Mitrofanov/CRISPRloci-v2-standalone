"""Contextual Possible calls keep their confidence and actual-member outputs.

The transaction fixture uses source-valid arrays and explicit fixture scoring;
these tests verify publication, not native extraction or biological accuracy.
"""
import copy
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from components.contextual_spacer_rescue import (
    EVIDENCE_FILENAME, METHOD, rescue_contextual_spacer_arrays,
)
from components.contextual_spacer_outputs import FASTA_FILENAMES
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.html_report import render_html_report
from components.likely_tandem_transaction import apply_likely_tandem_quarantine
from components.module_output_maker import OutputMaker
from tests.test_candidate_evidence import PARAMETERS
from tests.contextual_spacer_fixture import fixture


def publication_fixture():
    source, categories, pending, scorer = fixture(boundary_tie=True)
    output, audit = rescue_contextual_spacer_arrays(
        source, categories, pending, PARAMETERS, scorer, lambda: ())
    return source, output, audit


class ContextualSpacerExportTests(unittest.TestCase):
    def publish(self, root, state, *, evidence=False, accepted_only=False,
                machine=False, omit_selected=False, audit_override=False, audit=None):
        root.mkdir(exist_ok=True)
        source, provenance, recorded_audit = state
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
            "Downstream": {0: [""] * sum(map(len, projected[0].values()))}, "IS": {},
        }
        OutputMaker(
            file_path=fasta, result_path=str(root), categories=projected,
            parameters=dict(PARAMETERS, param_reconstruction_policy="conservative",
                            param_accepted_only=accepted_only),
            flags={"flag_fasta_report": True, "flag_html_report": True,
                   "output_profile": "machine" if machine else "standard"},
            pickle_result_path=str(root / "pickle"),
            json_result_path=str(root / "legacy-json"), non_array_data=annotations,
            list_features=[COMPACT_FEATURE_SCHEMA], header=">source\n",
            candidate_evidence_categories=provenance if evidence else None,
            likely_tandem_source_categories=provenance,
            likely_tandem_records=records, likely_tandem_stats=stats,
            contextual_spacer_evidence=audit if audit_override else recorded_audit,
        )
        return json.loads((root / "report.json").read_text())

    def test_all_outputs_keep_possible_ceiling_and_actual_member_uncertainty(self):
        state = publication_fixture()
        self.assertEqual(1, state[2]["stats"]["selections_committed"])
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            report = self.publish(root, state, evidence=True)
            selected = [row for row in report["arrays"] if "contextual_spacer_context" in row]
            self.assertEqual(1, len(selected))
            array = selected[0]
            summary = array["contextual_spacer_context"]
            self.assertEqual("Possible", array["category"])
            self.assertEqual(METHOD, array["classification_basis"]["method"])
            self.assertEqual(2, summary["member_count"])
            self.assertIs(False, summary["accepted_core_present"])
            self.assertEqual(array["certainty_score"], summary["display_own_score"])
            self.assertEqual("Possible", summary["confidence_category"])
            self.assertEqual(state[2]["stats"]["pending_observations"],
                             sum(map(len, state[1][4].values())))
            for filename in ("Possible_Candidates.txt", "report.html"):
                content = (root / filename).read_text()
                self.assertIn("Localized identical-spacer block", content)
                self.assertIn("remains Possible", content)
                self.assertIn("not calibrated probabilities", content)
                self.assertIn(summary["display_geometry_id"], content)
            for filename in FASTA_FILENAMES:
                content = (root / filename).read_text()
                self.assertIn("boundary_uncertain=true", content)
                self.assertIn("confidence_ceiling=Possible", content)
                self.assertIn(summary["display_geometry_id"], content)
            array_sequences = [line for line in (root / FASTA_FILENAMES[0]).read_text().splitlines()
                               if not line.startswith(">")]
            self.assertEqual([array["source_sequence"]], array_sequences)
            for filename, parts in zip(FASTA_FILENAMES[1:], ("repeats", "spacers")):
                self.assertEqual([part["sequence"] for part in array[parts]],
                                 [line for line in (root / filename).read_text().splitlines()
                                  if not line.startswith(">")])
            with (root / "Summary.csv").open(newline="") as handle:
                selected_csv = next(row for row in csv.DictReader(handle)
                                    if row["Category"] == "Possible")
            self.assertEqual("context_checked_boundary_family", selected_csv["Geometry status"])
            self.assertEqual("", selected_csv["Original core repeat count"])
            self.assertEqual(summary["display_geometry_id"], selected_csv["Display geometry ID"])
            legacy = json.loads((root / "legacy-json/source.json").read_text())
            self.assertEqual(summary, next(iter(legacy["Arrays"]["Possible"].values()))[
                "contextual_spacer_context"])
            retained = [json.loads(line) for line in gzip.decompress(
                (root / "candidate_evidence.jsonl.gz").read_bytes()).splitlines()]
            self.assertEqual(summary, next(row for row in retained
                                          if row.get("contextual_spacer_context"))[
                                              "contextual_spacer_context"])
            captures = [row["contextual_spacer_observation"]["observation_id"]
                        for row in retained if row.get("category") == "Low score"
                        and row.get("contextual_spacer_observation")]
            self.assertCountEqual([row["observation_id"] for row in state[2]["observations"]],
                                  captures)
            self.assertEqual(state[2], json.loads((root / EVIDENCE_FILENAME).read_text()))

    def test_retention_projection_and_machine_mode_preserve_accepted_proof(self):
        state = publication_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            accepted = []
            for index, (evidence, projection, machine) in enumerate((
                    (False, False, False), (True, False, False),
                    (False, True, False), (True, True, True))):
                folder = Path(temporary) / str(index)
                report = self.publish(folder, state, evidence=evidence,
                                      accepted_only=projection, machine=machine)
                accepted.append([row for row in report["arrays"]
                                 if row["category"] in ("Bona-fide", "Possible")])
                self.assertTrue((folder / EVIDENCE_FILENAME).exists())
                self.assertEqual(evidence, (folder / "candidate_evidence.jsonl.gz").exists())
                if machine:
                    self.assertFalse((folder / "Summary.csv").exists())
                    self.assertFalse(any((folder / name).exists() for name in FASTA_FILENAMES))
            for rows in accepted[1:]:
                self.assertEqual(accepted[0], rows)

    def test_missing_selected_or_corrupt_audit_aborts_before_publication(self):
        state = publication_fixture()
        broken = copy.deepcopy(state[2])
        broken["stats"]["selections_committed"] = 0
        for options in ({"omit_selected": True},
                        {"audit_override": True, "audit": None},
                        {"audit_override": True, "audit": broken}):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                with self.assertRaises(ValueError):
                    self.publish(root, state, **options)
                self.assertFalse((root / "report.json").exists())
                self.assertFalse((root / "Summary.csv").exists())

    def test_corrupt_canonical_confidence_member_and_score_cannot_render(self):
        with tempfile.TemporaryDirectory() as temporary:
            report = self.publish(Path(temporary), publication_fixture())
        for key, value in (("confidence_category", "Bona-fide"),
                           ("display_own_score", .001), ("member_count", 1)):
            broken = copy.deepcopy(report)
            array = next(row for row in broken["arrays"] if "contextual_spacer_context" in row)
            array["contextual_spacer_context"][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                render_html_report(broken)
        broken = copy.deepcopy(report)
        next(row for row in broken["arrays"] if "contextual_spacer_context" in row)["category"] = "Bona-fide"
        with self.assertRaises(ValueError):
            render_html_report(broken)

    def test_later_context_free_publication_removes_stale_audit(self):
        state = publication_fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.publish(root, state, evidence=True)
            self.assertTrue((root / EVIDENCE_FILENAME).exists())
            report = self.publish(root, (state[0], [{}, {}, {}, {}, {}], None))
            self.assertFalse((root / EVIDENCE_FILENAME).exists())
            self.assertFalse((root / "candidate_evidence.jsonl.gz").exists())
            self.assertFalse(any((root / name).exists() for name in FASTA_FILENAMES))
            self.assertEqual([], report["arrays"])


if __name__ == "__main__":
    unittest.main()
