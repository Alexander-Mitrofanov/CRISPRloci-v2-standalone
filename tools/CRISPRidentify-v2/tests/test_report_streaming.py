import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from components import html_report
from components import module_output_maker
from components.output_transaction import OutputTransaction
from components.report_model import ReportModelBuilder


class ReportStreamingTests(unittest.TestCase):
    def test_stream_matches_canonical_bytes_across_buffer_boundaries(self):
        model = {"unicode": "éλ🧬", "escaped": '"\\\n\t', "empty": [{}, []],
                 "values": [None, True, False, -0.0, 1e-30, 2 ** 70],
                 "rows": [{"z": index, "a": "ATGC" * 25} for index in range(2000)]}
        writes = []

        class RecordingStream(io.StringIO):
            def write(self, text):
                writes.append(len(text))
                return super().write(text)

        stream = RecordingStream()
        html_report._write_canonical_report_json(stream, model)
        self.assertEqual(html_report.canonical_report_json(model), stream.getvalue())
        self.assertGreater(len(writes), 2)
        self.assertLess(max(writes), 70 * 1024)

    def test_failed_stream_does_not_publish_or_leave_stage_files(self):
        for value in (float("nan"), float("inf"), object()):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "report.json").write_text("stale success")
                (root / "report.html").write_text("old html")
                model = {"a": ["ACGT" * 100] * 300, "z": value}
                with self.assertRaises((ValueError, TypeError)):
                    html_report._publish_report_model(model, root, False)
                self.assertFalse((root / "report.json").exists())
                self.assertEqual("old html", (root / "report.html").read_text())
                self.assertEqual([], list(root.glob(".crispridentify-*")))

    def test_stream_write_failure_discards_other_staged_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "report.json").write_text("stale success")
            with self.assertRaisesRegex(OSError, "disk full"), \
                    OutputTransaction() as transaction:
                transaction.stage_text(root / "other.csv", "new")
                with mock.patch.object(
                        html_report, "_write_canonical_report_json",
                        side_effect=OSError("disk full")):
                    html_report._publish_report_model({}, root, False, transaction)
            self.assertFalse((root / "report.json").exists())
            self.assertFalse((root / "other.csv").exists())
            self.assertEqual([], list(root.glob(".crispridentify-*")))

    def test_pipeline_reuses_validated_snapshot_after_legacy_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta = root / "source.fa"
            fasta.write_text(">source\nACGTACGT\n")
            parameters = {"param_accepted_only": False, "nested": {"value": [1]}}
            kwargs = dict(
                file_path=fasta, parameters=parameters,
                flags={"flag_fasta_report": False, "flag_html_report": False,
                       "output_profile": "standard"},
                result_path=root, pickle_result_path="", json_result_path="",
                categories=[{}, {}, {}, {}, {}], non_array_data={},
                list_features=[], header=">source\n",
            )
            original_build = ReportModelBuilder.build
            builds = []

            def tracked_build(builder):
                builds.append(builder)
                return original_build(builder)

            def legacy_mutation(**unused):
                parameters["nested"]["value"].append(2)

            with mock.patch.object(ReportModelBuilder, "build", tracked_build), \
                    mock.patch.object(module_output_maker, "SimpleOutputMaker"), \
                    mock.patch.object(module_output_maker, "SummaryOutputMaker"), \
                    mock.patch.object(module_output_maker, "SummaryMakerCSV",
                                      side_effect=legacy_mutation):
                module_output_maker.OutputMaker(**kwargs)
            report = json.loads((root / "report.json").read_text())
            self.assertEqual(1, len(builds))
            self.assertEqual([1, 2], parameters["nested"]["value"])
            self.assertEqual([1], report["run"]["parameters"]["nested"]["value"])


if __name__ == "__main__":
    unittest.main()
