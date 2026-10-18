"""Research retention is complete evidence, never extra ordinary arrays."""

from contextlib import ExitStack
import copy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from components import candidate_evidence
from components import module_output_maker
from components import output_transaction
from components.components_detection_refinement import CrisprCandidate
from components.likely_tandem_transaction import apply_likely_tandem_quarantine
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.output_transaction import OutputTransaction
from components.pipeline import Pipeline


PARAMETERS = {
    "param_min_avg_repeat_length": 21,
    "param_max_avg_repeat_length": 55,
    "param_min_avg_spacer_length": 18,
    "param_max_avg_spacer_length": 78,
    "param_min_repeats": 3,
    "param_max_identical_spacers": 4,
    "param_max_identical_cluster_spacers": 3,
    "param_accepted_only": False,
}
CATEGORY_NAMES = (
    "Bona-fide", "Alternative", "Possible", "Possible Discarded", "Low score",
)


def _fixture():
    repeat = "GTTTCAGTATCCTCGCTGTTGGATCGTTGCCTAGACC"
    spacers = [
        "A" * 10 + "C" * 10,
        "C" * 10 + "G" * 10,
        "G" * 10 + "T" * 10,
        "T" * 10 + "A" * 10,
        "AG" * 10,
        "CT" * 10,
        "AC" * 10,
    ]
    source = repeat + "".join(spacer + repeat for spacer in spacers)
    starts = [index * (len(repeat) + 20) for index in range(8)]

    def candidate(indexes):
        positions = [starts[index] for index in indexes]
        return CrisprCandidate(
            list_repeats=[repeat] * len(positions),
            list_repeats_gaped=[repeat] * len(positions),
            list_spacers=[
                source[left + len(repeat):right]
                for left, right in zip(positions, positions[1:])
            ],
            list_repeat_starts=positions,
        )

    scope = (0, len(source))
    categories = [
        {scope: [[0.90, candidate([0, 1, 2]), []]]},
        {scope: [[0.70, candidate([0, 1, 2, 3, 4]), []]]},
        {scope: [[0.55, candidate([5, 6, 7]), []]]},
        {scope: [[0.52, candidate([1, 2, 3, 4]), []]]},
        {scope: [[0.30, candidate([2, 3, 4]), []]]},
    ]
    return source, categories


def _row(categories, index):
    return next(iter(categories[index].values()))[0]


def _canonical_links(categories):
    arrays = []
    for index in (0, 2):
        candidate = _row(categories, index)[1]
        arrays.append({
            "id": "accepted-{}".format(index),
            "category": CATEGORY_NAMES[index],
            "repeats": [
                {"source_interval": {"start": start + 1,
                                     "end": start + len(repeat)}}
                for start, repeat in zip(
                    candidate.list_repeat_starts, candidate.list_repeats)
            ],
        })
    return {"arrays": arrays}


def _records(path):
    return [json.loads(line) for line in gzip.decompress(path.read_bytes()).splitlines()]


def _hidden_files(root):
    return sorted(path.name for path in root.glob(".crispridentify-*"))


def _stage(publication, path, source, categories, *, feature_names=None):
    candidate_evidence.stage_candidate_evidence(
        publication, path,
        source_header="source fixture",
        source_sequence=source,
        categories=categories,
        parameters=dict(PARAMETERS),
        feature_names=[] if feature_names is None else feature_names,
        canonical_report=_canonical_links(categories),
    )


class CandidateEvidenceRecordTests(unittest.TestCase):
    def test_all_five_categories_are_source_valid_and_not_additional_arrays(self):
        source, categories = _fixture()
        before = copy.deepcopy(categories)
        for category in categories:
            for rows in category.values():
                for row in rows:
                    self.assertTrue(build_advanced_fuzzy_search_filter(PARAMETERS)(row[1]))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / candidate_evidence.EVIDENCE_FILENAME
            with OutputTransaction() as publication:
                _stage(publication, path, source, categories)
                self.assertFalse(path.exists())
                publication.publish()
            records = _records(path)

        header, footer = records[0], records[-1]
        candidates = [row for row in records if row["record_type"] == "candidate"]
        scopes = [row for row in records if row["record_type"] == "comparison_scope"]
        self.assertTrue(header["research_only"])
        self.assertTrue(header["not_additional_arrays"])
        self.assertEqual("final_before_projection", header["capture_stage"])
        self.assertEqual(hashlib.sha256(source.encode()).hexdigest(), header["source"]["sha256"])
        self.assertEqual(list(CATEGORY_NAMES), [row["category"] for row in candidates])
        self.assertEqual(2, sum(row["accepted"] for row in candidates))
        self.assertTrue(all(row["source_valid"] for row in candidates))
        self.assertTrue(all(row["structural_filter"]["structural_valid"] for row in candidates))
        self.assertEqual(["accepted-0"], candidates[0]["canonical_array_ids"])
        self.assertEqual(["accepted-2"], candidates[2]["canonical_array_ids"])
        self.assertEqual([], candidates[1]["canonical_array_ids"])
        self.assertEqual(5, candidates[1]["repeat_count"])
        self.assertEqual(0.70, candidates[1]["score"])
        self.assertEqual("not_recorded", candidates[1]["negative_evidence"]["accepted_family_tandem_review_status"])
        self.assertEqual(1, len(scopes))
        self.assertFalse(scopes[0]["is_biological_locus"])
        self.assertEqual(2, len(scopes[0]["accepted_candidate_ids"]))
        self.assertTrue(footer["complete"])
        self.assertEqual(5, footer["candidate_observations"])
        self.assertEqual(dict.fromkeys(CATEGORY_NAMES, 1), footer["category_counts"])
        for index in range(5):
            original, retained = _row(before, index), _row(categories, index)
            self.assertEqual(original[0], retained[0])
            self.assertEqual(original[1].__dict__, retained[1].__dict__)

    def test_gzip_is_deterministic_and_footer_binds_all_preceding_records(self):
        source, categories = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = [root / "first.gz", root / "second.gz"]
            for path in paths:
                with OutputTransaction() as publication:
                    _stage(publication, path, source, categories)
                    publication.publish()
            self.assertEqual(paths[0].read_bytes(), paths[1].read_bytes())
            self.assertEqual(b"\x00\x00\x00\x00", paths[0].read_bytes()[4:8])
            lines = gzip.decompress(paths[0].read_bytes()).splitlines(keepends=True)
        footer = json.loads(lines[-1])
        self.assertEqual(hashlib.sha256(b"".join(lines[:-1])).hexdigest(),
                         footer["uncompressed_records_sha256"])
        tampered = lines[:-1].copy()
        tampered[0] = tampered[0].replace(b'"research_only":true', b'"research_only":false')
        self.assertNotEqual(hashlib.sha256(b"".join(tampered)).hexdigest(),
                            footer["uncompressed_records_sha256"])

    def test_duplicate_observations_do_not_become_distinct_content_evidence(self):
        source, categories = _fixture()
        rows = next(iter(categories[1].values()))
        rows.append(copy.deepcopy(rows[0]))
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            with OutputTransaction() as publication:
                _stage(publication, path, source, categories)
                publication.publish()
            alternatives = [row for row in _records(path) if row.get("category") == "Alternative"]
        self.assertEqual(2, len(alternatives))
        self.assertEqual(alternatives[0]["candidate_id"], alternatives[1]["candidate_id"])
        self.assertEqual(alternatives[0]["source_geometry_id"], alternatives[1]["source_geometry_id"])
        self.assertNotEqual(alternatives[0]["observation_id"], alternatives[1]["observation_id"])

    def test_alignment_variants_and_internal_boundary_variants_have_distinct_identities(self):
        source, categories = _fixture()
        original = _row(categories, 1)[1]
        different_alignment = CrisprCandidate(
            list_repeats=list(original.list_repeats),
            list_repeats_gaped=[repeat[:index + 1] + " " + repeat[index + 1:]
                               for index, repeat in enumerate(original.list_repeats)],
            list_spacers=list(original.list_spacers),
            list_repeat_starts=list(original.list_repeat_starts),
        )
        # Array outer endpoints remain identical while the first repeat ends
        # one base earlier; that base becomes part of its following spacer.
        different_boundary = CrisprCandidate(
            list_repeats=[original.list_repeats[0][:-1]] + original.list_repeats[1:],
            list_repeats_gaped=[original.list_repeats[0][:-1] + "-"] + original.list_repeats[1:],
            list_spacers=[original.list_repeats[0][-1] + original.list_spacers[0]] + original.list_spacers[1:],
            list_repeat_starts=list(original.list_repeat_starts),
        )
        next(iter(categories[1].values())).extend([
            [0.70, different_alignment, []], [0.70, different_boundary, []],
        ])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            with OutputTransaction() as publication:
                _stage(publication, path, source, categories)
                publication.publish()
            alternatives = [row for row in _records(path) if row.get("category") == "Alternative"]
        baseline, alignment, boundary = alternatives
        self.assertEqual(baseline["source_geometry_id"], alignment["source_geometry_id"])
        self.assertNotEqual(baseline["alignment_id"], alignment["alignment_id"])
        self.assertNotEqual(baseline["candidate_id"], alignment["candidate_id"])
        self.assertEqual(baseline["source_span"], boundary["source_span"])
        self.assertEqual(baseline["repeat_count"], boundary["repeat_count"])
        self.assertNotEqual(baseline["source_geometry_id"], boundary["source_geometry_id"])

    def test_numpy_feature_values_are_preserved_without_rescoring(self):
        source, categories = _fixture()
        names = [["number_repeats", "repeat_len"]]
        for index in range(5):
            row = _row(categories, index)
            row[2] = [np.asarray([[len(row[1].list_repeats), len(row[1].list_repeats[0])]], dtype=float)]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            with OutputTransaction() as publication:
                _stage(publication, path, source, categories, feature_names=names)
                publication.publish()
            alternatives = [row for row in _records(path) if row.get("category") == "Alternative"]
        self.assertEqual({"number_repeats": 5.0, "repeat_len": 37.0}, alternatives[0]["feature_values"])
        self.assertEqual(0.70, alternatives[0]["score"])
        self.assertEqual("existing_final_row", alternatives[0]["score_origin"])
        self.assertTrue(alternatives[0]["feature_repeat_count_consistent"])

    def test_inconsistent_feature_count_is_explicit_without_rewriting_score(self):
        source, categories = _fixture()
        for index in range(5):
            row = _row(categories, index)
            row[2] = [np.asarray([[len(row[1].list_repeats), 37.0]], dtype=float)]
        _row(categories, 1)[2] = [np.asarray([[3.0, 37.0]])]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            with OutputTransaction() as publication:
                _stage(publication, path, source, categories,
                       feature_names=[["number_repeats", "repeat_len"]])
                publication.publish()
            alternative = next(row for row in _records(path)
                               if row.get("category") == "Alternative")
        self.assertEqual(5, alternative["repeat_count"])
        self.assertEqual(3.0, alternative["feature_values"]["number_repeats"])
        self.assertFalse(alternative["feature_repeat_count_consistent"])
        self.assertEqual(0.70, alternative["score"])

    def test_invalid_hidden_source_geometry_fails_before_any_artifact_is_staged(self):
        source, categories = _fixture()
        _row(categories, 1)[1].list_repeat_starts[0] += 1
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "evidence.gz"
            with OutputTransaction() as publication:
                with self.assertRaisesRegex(ValueError, "source reconstruction failed"):
                    _stage(publication, path, source, categories)
                self.assertEqual((), publication.targets)
            self.assertFalse(path.exists())
            self.assertEqual([], _hidden_files(root))

    def test_alignment_corruption_and_nonboolean_eligibility_fail_closed(self):
        for problem in ("alignment", "eligibility"):
            with self.subTest(problem=problem):
                source, categories = _fixture()
                candidate = _row(categories, 1)[1]
                if problem == "alignment":
                    candidate.list_gaped_repeats[0] = "A" + candidate.list_gaped_repeats[0][1:]
                else:
                    candidate.source_path_classification_eligible = "True"
                with tempfile.TemporaryDirectory() as temporary:
                    root = Path(temporary)
                    with OutputTransaction() as publication:
                        with self.assertRaisesRegex(ValueError, "alignment/source mismatch|eligibility must be boolean"):
                            _stage(publication, root / "evidence.gz", source, categories)
                    self.assertEqual([], _hidden_files(root))

    def test_missing_nonfinite_or_nonnumeric_scores_cannot_enter_research_inventory(self):
        for score in (None, float("nan"), float("inf"), "0.7", True):
            with self.subTest(score=score), tempfile.TemporaryDirectory() as temporary:
                source, categories = _fixture()
                _row(categories, 1)[0] = score
                root = Path(temporary)
                with OutputTransaction() as publication:
                    with self.assertRaises(ValueError):
                        _stage(publication, root / "evidence.gz", source, categories)
                self.assertFalse((root / "evidence.gz").exists())
                self.assertEqual([], _hidden_files(root))


class CandidateEvidencePublicationTests(unittest.TestCase):
    def _write(self, root, source, categories, evidence, *, accepted_only=False):
        fasta = root / "source.fa"
        fasta.write_text(">source fixture\n{}\n".format(source))
        parameters = dict(PARAMETERS, param_accepted_only=accepted_only)
        audit = {}
        ordinary = categories
        if accepted_only:
            projected, records, stats = apply_likely_tandem_quarantine(categories)
            ordinary = [projected[index] if index in (0, 2) else {} for index in range(5)]
            audit = dict(likely_tandem_records=records, likely_tandem_stats=stats,
                         likely_tandem_source_categories=categories)
        module_output_maker.OutputMaker(
            file_path=fasta, parameters=parameters,
            flags={"output_profile": "machine", "flag_html_report": False},
            result_path=root, pickle_result_path="", json_result_path="",
            categories=ordinary, non_array_data={}, list_features=[],
            header="source fixture", candidate_evidence_categories=evidence,
            **audit,
        )

    def test_enabled_retention_preserves_report_bytes_ids_and_accepted_count(self):
        source, categories = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(root, source, categories, None)
            ordinary_bytes = (root / "report.json").read_bytes()
            self._write(root, source, categories, categories)
            self.assertEqual(ordinary_bytes, (root / "report.json").read_bytes())
            report = json.loads(ordinary_bytes)
            records = _records(root / candidate_evidence.EVIDENCE_FILENAME)
            self.assertEqual(2, report["summary"]["accepted_array_count"])
            self.assertEqual(3, report["summary"]["array_count"])
            self.assertEqual({"Bona-fide", "Possible", "Low score"},
                             {row["category"] for row in report["arrays"]})
            report_accepted = {row["id"] for row in report["arrays"] if row["category"] in ("Bona-fide", "Possible")}
            linked = {array_id for row in records if row["record_type"] == "candidate"
                      for array_id in row["canonical_array_ids"]}
            self.assertEqual(report_accepted, linked)
            self.assertEqual([], _hidden_files(root))

    def test_accepted_only_output_still_retains_all_five_research_categories(self):
        source, categories = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(root, source, categories, None, accepted_only=True)
            ordinary_bytes = (root / "report.json").read_bytes()
            self._write(root, source, categories, categories, accepted_only=True)
            self.assertEqual(ordinary_bytes, (root / "report.json").read_bytes())
            report = json.loads(ordinary_bytes)
            self.assertEqual(2, len(report["arrays"]))
            self.assertEqual(2, report["summary"]["accepted_array_count"])
            candidates = [row for row in _records(root / candidate_evidence.EVIDENCE_FILENAME)
                          if row["record_type"] == "candidate"]
            self.assertEqual(set(CATEGORY_NAMES), {row["category"] for row in candidates})

    def test_disabling_retention_removes_stale_sidecar_with_fresh_report(self):
        source, categories = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            evidence = root / candidate_evidence.EVIDENCE_FILENAME
            self._write(root, source, categories, categories)
            ordinary_bytes = (root / "report.json").read_bytes()
            self.assertTrue(evidence.exists())
            self._write(root, source, categories, None)
            self.assertFalse(evidence.exists())
            self.assertEqual(ordinary_bytes, (root / "report.json").read_bytes())
            self.assertEqual([], _hidden_files(root))

    def test_invalid_unreported_alternative_prevents_completion_publication(self):
        source, categories = _fixture()
        invalid = copy.deepcopy(categories)
        _row(invalid, 1)[1].list_repeat_starts[0] += 1
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write(root, source, categories, categories)
            prior_evidence = (root / candidate_evidence.EVIDENCE_FILENAME).read_bytes()
            with self.assertRaisesRegex(ValueError, "source reconstruction failed"):
                self._write(root, source, categories, invalid)
            self.assertFalse((root / "report.json").exists())
            self.assertEqual(prior_evidence, (root / candidate_evidence.EVIDENCE_FILENAME).read_bytes())
            self.assertEqual([], _hidden_files(root))

    def test_binary_or_report_publication_failure_rolls_back_without_completion(self):
        source, categories = _fixture()
        for failed_target in (candidate_evidence.EVIDENCE_FILENAME, "report.json"):
            with self.subTest(failed_target=failed_target), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                self._write(root, source, categories, categories)
                sidecar = root / candidate_evidence.EVIDENCE_FILENAME
                old_evidence = sidecar.read_bytes()
                real_replace = output_transaction.os.replace

                def fail_target(source_path, target_path):
                    if (Path(source_path).name.startswith(".crispridentify-stage-")
                            and Path(target_path).name == failed_target):
                        raise OSError("injected research publication failure")
                    return real_replace(source_path, target_path)

                with mock.patch.object(output_transaction.os, "replace", side_effect=fail_target):
                    with self.assertRaisesRegex(OSError, "injected research publication failure"):
                        self._write(root, source, categories, categories)
                self.assertFalse((root / "report.json").exists())
                self.assertEqual(old_evidence, sidecar.read_bytes())
                self.assertEqual([], _hidden_files(root))

    def test_binary_writer_failure_discards_partial_stage_and_completion(self):
        source, categories = _fixture()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "report.json").write_text("stale successful run")
            with mock.patch.object(candidate_evidence.gzip.GzipFile, "write", side_effect=OSError("disk full")):
                with self.assertRaisesRegex(OSError, "disk full"):
                    self._write(root, source, categories, categories)
            self.assertFalse((root / "report.json").exists())
            self.assertFalse((root / candidate_evidence.EVIDENCE_FILENAME).exists())
            self.assertEqual([], _hidden_files(root))


class CandidateEvidencePipelineSeamTests(unittest.TestCase):
    def test_real_projection_clears_ordinary_categories_after_research_snapshot(self):
        # Execute the actual constructor ordering and projection, while expensive
        # detection/repair methods are replaced with a fixed evaluated inventory.
        source, categories = _fixture()
        stage_methods = (
            "_get_header", "_run_detection", "_run_detection_refinement",
            "_validate_refined_candidates", "_collapse_source_backed_phase_equivalent_candidates",
            "_run_fast_routed_start_end_rescue", "_results_enhancement",
            "_validate_candidate_categories", "_normalize_post_enhancement_source_phase",
            "_reevaluate_enhanced_candidates", "_run_fast_two_repeat_terminal_repair",
            "_run_fast_accepted_family_rescue", "_run_fast_bounded_phase_repair",
            "_complete_attested_terminal_repairs", "_complete_pending_internal_repeats",
            "_run_configured_source_backed_augmentation", "_annotate_model_score_classification_bases",
            "_demote_accepted_family_tandems_transactionally",
            "_complete_accepted_core_boundaries_transactionally",
            "_recover_source_intact_candidates_transactionally",
            "_repair_source_periodic_lattices_transactionally",
            "_recover_source_regular_subpaths_transactionally",
            "_quarantine_likely_tandems_transactionally", "_run_non_crispr_computation",
            "_validate_and_write_output",
        )

        def evaluated(pipeline):
            pipeline.categories = copy.deepcopy(categories)
            pipeline.genome_sequence = source

        def final_assessment(pipeline):
            # The extra observation proves capture follows the last assessment.
            rows = next(iter(pipeline.categories[1].values()))
            rows.append(copy.deepcopy(rows[0]))

        for enabled in (False, True):
            with self.subTest(enabled=enabled), tempfile.TemporaryDirectory() as temporary, ExitStack() as stack:
                for name in stage_methods:
                    stack.enter_context(mock.patch.object(Pipeline, name, autospec=True))
                stack.enter_context(mock.patch.object(Pipeline, "_run_evaluation", autospec=True, side_effect=evaluated))
                stack.enter_context(mock.patch.object(Pipeline, "_evaluate_circular_origin_joins", autospec=True, side_effect=final_assessment))
                pipeline = Pipeline(
                    temporary, "", "", str(Path(temporary) / "source.fa"),
                    dict(PARAMETERS, param_accepted_only=True), {}, False,
                    retain_candidate_evidence=enabled,
                )
                self.assertEqual([1, 0, 1, 0, 0],
                                 [sum(map(len, group.values())) for group in pipeline.categories])
                if enabled:
                    self.assertEqual([1, 2, 1, 1, 1], [sum(map(len, group.values()))
                                                      for group in pipeline._candidate_evidence_categories])
                else:
                    self.assertIsNone(pipeline._candidate_evidence_categories)


if __name__ == "__main__":
    unittest.main()
