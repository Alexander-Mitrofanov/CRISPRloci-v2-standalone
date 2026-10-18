import copy
import hashlib
import os
import tempfile
import unittest
from dataclasses import replace
from unittest import mock

from components.likely_tandem import evaluate_likely_tandem
from components.html_report import render_html_report
from components.components_detection_refinement import CrisprCandidate
from components.likely_tandem_transaction import (
    LikelyTandemQuarantineStats,
    apply_likely_tandem_quarantine,
)
from components.report_extensions import (
    ReportExtensionError,
    attach_report_extensions,
)
from components.module_output_maker import OutputMaker
from components.pipeline import Pipeline
from components.report_model import ReportModelBuilder, ReportModelError
from components.folder_scheduler import pipeline_result_key


class _Candidate:
    def __init__(self, spacers):
        repeat = "ACGTAC"
        self.list_spacers = list(spacers)
        self.list_repeats = [repeat] * (len(spacers) + 1)
        self.list_repeat_starts = [0]
        for spacer in spacers:
            self.list_repeat_starts.append(
                self.list_repeat_starts[-1] + len(repeat) + len(spacer)
            )


def _source(spacers, repeat="ACGTAC"):
    return repeat + repeat.join(spacers) + repeat


def _report(source_sequence):
    return {
        "schema": {
            "name": "CRISPRidentify-report",
            "version": "1.2.0",
        },
        "coordinate_system": {
            "indexing": 1,
            "end_inclusive": True,
            "reference": "input_forward_strand",
        },
        "source": {
            "id": "fixture",
            "header": "fixture",
            "file_name": "fixture.fa",
            "length": len(source_sequence),
            "sha256": hashlib.sha256(
                source_sequence.encode("ascii")
            ).hexdigest(),
            "sha256_scope": "uppercase_sequence",
        },
        "run": {
            "parameters": {
                "param_accepted_only": True,
                "param_min_avg_spacer_length": 18,
                "param_max_avg_spacer_length": 78,
            },
            "analyses": {},
        },
        "summary": {
            "array_count": 0,
            "accepted_array_count": 0,
            "category_counts": {
                "Bona-fide": 0,
                "Possible": 0,
                "Low score": 0,
            },
        },
        "validation": {
            "status": "valid",
            "method": "independent_original_sequence_slicing",
            "validated_array_count": 0,
            "invalid_array_count": 0,
        },
        "arrays": [],
    }


class LikelyTandemReleaseTests(unittest.TestCase):
    def test_pipeline_invalidates_stale_report_before_first_stage(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "source.fa")
            result_root = os.path.join(directory, "results")
            result_path = os.path.join(
                result_root, pipeline_result_key(fasta)
            )
            os.makedirs(result_path)
            completion_marker = os.path.join(result_path, "report.json")
            with open(completion_marker, "w", encoding="ascii") as handle:
                handle.write('{"stale": true}\n')

            with mock.patch.object(
                    Pipeline,
                    "_get_header",
                    side_effect=RuntimeError("early pipeline failure")):
                with self.assertRaisesRegex(
                        RuntimeError, "early pipeline failure"):
                    Pipeline(
                        result_root,
                        "",
                        "",
                        fasta,
                        {},
                        {},
                        False,
                    )

            self.assertFalse(os.path.exists(completion_marker))

    def test_complete_link_boundary_and_chain(self):
        at_boundary = evaluate_likely_tandem(
            ("A" * 20, "A" * 19 + "C")
        )
        chain = evaluate_likely_tandem((
            "A" * 20,
            "C" + "A" * 19,
            "CC" + "A" * 18,
        ))

        self.assertTrue(at_boundary.likely_tandem)
        self.assertFalse(chain.likely_tandem)

    def test_output_maker_requires_complete_audit_context(self):
        spacers = ("A" * 20, "A" * 20)
        source_sequence = _source(spacers)
        repeat = "ACGTAC"
        candidate = CrisprCandidate(
            list_repeats=[repeat, repeat, repeat],
            list_repeats_gaped=[repeat, repeat, repeat],
            list_spacers=list(spacers),
            list_repeat_starts=[0, 26, 52],
        )
        categories = [{}, {}, {(0, 58): [[
            0.55, candidate, [["features"]]
        ]]}, {}, {}]
        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "source.fa")
            with open(fasta, "w", encoding="ascii") as handle:
                handle.write(">source\n{}\n".format(source_sequence))
            marker = os.path.join(directory, "report.json")
            with open(marker, "w", encoding="ascii") as handle:
                handle.write('{"stale": true}\n')
            with self.assertRaisesRegex(
                    ReportExtensionError, "still qualifies"):
                OutputMaker(
                    file_path=fasta,
                    parameters={"param_accepted_only": False},
                    flags={
                        "flag_fasta_report": False,
                        "flag_html_report": False,
                        "output_profile": "standard",
                    },
                    result_path=directory,
                    pickle_result_path="",
                    json_result_path="",
                    categories=categories,
                    non_array_data={},
                    list_features=[],
                    header=">source\n",
                )
            self.assertFalse(os.path.exists(marker))

            with open(marker, "w", encoding="ascii") as handle:
                handle.write('{"stale": true}\n')
            with self.assertRaisesRegex(
                    ValueError, "accepted-only writers require explicit"):
                OutputMaker(
                    file_path=fasta,
                    parameters={"param_accepted_only": True},
                    flags={
                        "flag_fasta_report": False,
                        "flag_html_report": False,
                        "output_profile": "standard",
                    },
                    result_path=directory,
                    pickle_result_path="",
                    json_result_path="",
                    categories=[{}, {}, {}, {}, {}],
                    non_array_data={},
                    list_features=[],
                    header=">source\n",
                )
            self.assertFalse(os.path.exists(marker))

    def test_html_rejects_a_forged_empty_ledger(self):
        spacers = ("A" * 20, "A" * 20)
        source_sequence = _source(spacers)
        repeat = "ACGTAC"
        candidate = CrisprCandidate(
            list_repeats=[repeat, repeat, repeat],
            list_repeats_gaped=[repeat, repeat, repeat],
            list_spacers=list(spacers),
            list_repeat_starts=[0, 26, 52],
        )
        categories = [{}, {}, {(0, 58): [[
            0.55, candidate, [["features"]]
        ]]}, {}, {}]
        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "source.fa")
            with open(fasta, "w", encoding="ascii") as handle:
                handle.write(">source\n{}\n".format(source_sequence))
            report = ReportModelBuilder(
                file_path=fasta,
                categories=categories,
                provenance_categories=categories,
                non_array_data={},
                list_feature_names=[],
                parameters={"param_accepted_only": False},
                flags={},
            ).build()

        empty_categories = [{}, {}, {}, {}, {}]
        _result, records, stats = apply_likely_tandem_quarantine(
            empty_categories
        )
        scaffold = attach_report_extensions(
            _report(source_sequence),
            likely_tandem_records=records,
            likely_tandem_stats=stats,
            likely_tandem_source_categories=empty_categories,
            source_sequence=source_sequence,
        )
        extension = copy.deepcopy(
            scaffold["extensions"]["likely_tandem"]
        )
        extension["transaction_stats"].update({
            "total_rows_before": 1,
            "total_rows_after": 1,
            "accepted_rows_reviewed": 1,
            "accepted_rows_before": 1,
            "accepted_rows_after": 1,
            "retained_count": 1,
        })
        report["extensions"] = {"likely_tandem": extension}

        with self.assertRaisesRegex(ValueError, "still qualifies"):
            render_html_report(report)

        hidden_spacers = copy.deepcopy(report)
        for element in hidden_spacers["arrays"][0]["architecture"]:
            if element.get("type") == "spacer":
                element["type"] = "gap"
        with self.assertRaisesRegex(
                ValueError, "architecture.*tandem replay"):
            render_html_report(hidden_spacers)

    def test_transaction_and_report_extension_preserve_quarantine(self):
        spacers = ("A" * 20, "A" * 20)
        source_sequence = _source(spacers)
        candidate = _Candidate(spacers)
        row = [0.55, candidate, [["features"]]]
        categories = [{}, {}, {(0, 58): [row]}, {}, {}]

        result, records, stats = apply_likely_tandem_quarantine(categories)
        extension = attach_report_extensions(
            _report(source_sequence),
            likely_tandem_records=records,
            likely_tandem_stats=stats,
            likely_tandem_source_categories=categories,
            source_sequence=source_sequence,
        )["extensions"]["likely_tandem"]

        self.assertEqual({}, result[2])
        self.assertEqual(1, len(result[4][(0, 58)]))
        self.assertEqual(["tandem"], result[4][(0, 58)][0][3]["tags"])
        self.assertEqual("likely_tandem", extension["disposition"])
        self.assertEqual("tandem", extension["tag"])
        self.assertEqual("Low score", extension["reported_category"])
        self.assertEqual(
            "omitted_by_accepted_only",
            extension["arrays"][0]["core_projection"],
        )
        self.assertEqual(1, extension["array_count"])
        report = attach_report_extensions(
            _report(source_sequence),
            likely_tandem_records=records,
            likely_tandem_stats=stats,
            likely_tandem_source_categories=categories,
            source_sequence=source_sequence,
        )
        rendered = render_html_report(report)
        self.assertIn("Tandem demotion audit", rendered)
        self.assertIn("Omitted by accepted-only projection", rendered)

        forged_evidence = copy.deepcopy(report)
        evidence = forged_evidence["extensions"]["likely_tandem"][
            "arrays"
        ][0]["evidence"]
        evidence["minimum_pairwise_edit_identity"] = 0.0
        evidence["maximum_pairwise_edit_distance"] = 999
        with self.assertRaisesRegex(ValueError, "evidence replay"):
            render_html_report(forged_evidence)

        malformed_score = copy.deepcopy(report)
        malformed_score["extensions"]["likely_tandem"]["arrays"][0][
            "score"
        ] = "not-a-score"
        with self.assertRaisesRegex(ValueError, "ledger row is malformed"):
            render_html_report(malformed_score)

        out_of_source = copy.deepcopy(report)
        ledger = out_of_source["extensions"]["likely_tandem"]["arrays"][0]
        ledger["source_interval"]["start"] += 1000
        ledger["source_interval"]["end"] += 1000
        for field in ("repeats", "spacers"):
            for element in ledger["source_architecture"][field]:
                element["source_interval"]["start"] += 1000
                element["source_interval"]["end"] += 1000
        with self.assertRaisesRegex(ValueError, "ledger interval"):
            render_html_report(out_of_source)

        for field in ("accepted_rows_reviewed", "retained_count"):
            forged = replace(stats, **{field: 999})
            with self.subTest(field=field):
                with self.assertRaisesRegex(
                        ReportExtensionError, "internally inconsistent"):
                    attach_report_extensions(
                        _report(source_sequence),
                        likely_tandem_records=records,
                        likely_tandem_stats=forged,
                        likely_tandem_source_categories=categories,
                        source_sequence=source_sequence,
                    )

    def test_historical_two_field_row_projects_to_full_report(self):
        spacers = ("A" * 20, "A" * 20)
        source_sequence = _source(spacers)
        repeat = "ACGTAC"
        candidate = CrisprCandidate(
            list_repeats=[repeat, repeat, repeat],
            list_repeats_gaped=[repeat, repeat, repeat],
            list_spacers=list(spacers),
            list_repeat_starts=[0, 26, 52],
        )
        source_categories = [
            {}, {}, {(0, 58): [[0.55, candidate]]}, {}, {}
        ]
        result, records, stats = apply_likely_tandem_quarantine(
            source_categories
        )
        self.assertEqual(4, len(result[4][(0, 58)][0]))

        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "source.fa")
            with open(fasta, "w", encoding="ascii") as handle:
                handle.write(">source\n{}\n".format(source_sequence))
            core = ReportModelBuilder(
                file_path=fasta,
                categories=result,
                provenance_categories=source_categories,
                non_array_data={},
                list_feature_names=[],
                parameters={"param_accepted_only": False},
                flags={},
            ).build()
            report = attach_report_extensions(
                core,
                likely_tandem_records=records,
                likely_tandem_stats=stats,
                likely_tandem_source_categories=source_categories,
                source_sequence=source_sequence,
            )

        self.assertEqual("Low score", report["arrays"][0]["category"])
        self.assertEqual(["tandem"], report["arrays"][0]["tags"])
        render_html_report(report)

    def test_invalid_tandem_tag_fails_before_legacy_publication(self):
        source_sequence = _source(("A" * 20, "A" * 20))
        repeat = "ACGTAC"
        candidate = CrisprCandidate(
            list_repeats=[repeat, repeat, repeat],
            list_repeats_gaped=[repeat, repeat, repeat],
            list_spacers=["A" * 20, "A" * 20],
            list_repeat_starts=[0, 26, 52],
        )
        categories = [{}, {}, {}, {}, {
            (0, len(source_sequence)): [[
                0.4,
                candidate,
                [["features"]],
                {"tags": ["tandem"]},
            ]],
        }]
        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "source.fa")
            with open(fasta, "w", encoding="ascii") as handle:
                handle.write(">source\n{}\n".format(source_sequence))
            with self.assertRaisesRegex(
                    ReportModelError, "classification_basis|tandem tag"):
                OutputMaker(
                    file_path=fasta,
                    parameters={},
                    flags={
                        "flag_fasta_report": False,
                        "flag_html_report": False,
                        "output_profile": "standard",
                    },
                    result_path=directory,
                    pickle_result_path="",
                    json_result_path="",
                    categories=categories,
                    non_array_data={},
                    list_features=[],
                    header=">source\n",
                )
            self.assertFalse(os.path.exists(os.path.join(
                directory, "Low_Score_Candidates.txt"
            )))

    def test_forged_empty_ledger_fails_before_legacy_publication(self):
        spacers = ("A" * 20, "A" * 20)
        source_sequence = _source(spacers)
        repeat = "ACGTAC"
        candidate = CrisprCandidate(
            list_repeats=[repeat, repeat, repeat],
            list_repeats_gaped=[repeat, repeat, repeat],
            list_spacers=list(spacers),
            list_repeat_starts=[0, 26, 52],
        )
        categories = [{}, {}, {(0, 58): [[
            0.55, candidate, [["features"]]
        ]]}, {}, {}]
        forged_stats = LikelyTandemQuarantineStats(
            total_rows_before=1,
            total_rows_after=1,
            accepted_rows_reviewed=1,
            accepted_rows_before=1,
            accepted_rows_after=1,
            likely_tandem_count=0,
            retained_count=1,
            bona_fide_quarantined=0,
            possible_quarantined=0,
            diagnostic_rows_before=0,
            diagnostic_rows_after=0,
        )
        with tempfile.TemporaryDirectory() as directory:
            fasta = os.path.join(directory, "source.fa")
            with open(fasta, "w", encoding="ascii") as handle:
                handle.write(">source\n{}\n".format(source_sequence))
            with self.assertRaisesRegex(
                    ReportExtensionError, "still qualifies|incomplete"):
                OutputMaker(
                    file_path=fasta,
                    parameters={"param_accepted_only": False},
                    flags={
                        "flag_fasta_report": False,
                        "flag_html_report": True,
                        "output_profile": "standard",
                    },
                    result_path=directory,
                    pickle_result_path="",
                    json_result_path="",
                    categories=categories,
                    non_array_data={},
                    list_features=[],
                    header=">source\n",
                    likely_tandem_records=(),
                    likely_tandem_stats=forged_stats,
                    likely_tandem_source_categories=categories,
                )
            self.assertFalse(os.path.exists(os.path.join(
                directory, "Possible_Candidates.txt"
            )))
            self.assertFalse(os.path.exists(os.path.join(
                directory, "report.json"
            )))


if __name__ == "__main__":
    unittest.main()
