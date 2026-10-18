import argparse
import ast
import pathlib
import unittest
from unittest import mock

from components import module_detection
from components.cli_utils import parse_boolean
from components.components_detection import FuzzySearch
from components.module_detection import Detection


V2_ROOT = pathlib.Path(__file__).resolve().parents[1]
CLI_PATH = V2_ROOT / "CRISPRidentify.py"
ERROR_LIMIT = "{i<=3,d<=3,s<=3,i+d+s<=6}"
REPEAT = "A" * 21
PHASE_FIXTURE = "C" + REPEAT + "T" * 20 + REPEAT
SOURCE_START = 137


def _cli_tree():
    return ast.parse(CLI_PATH.read_text(encoding="utf-8"))


def _argument_call(tree, option):
    matches = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "add_argument"
        and node.args
        and isinstance(node.args[0], ast.Constant)
        and node.args[0].value == option
    ]
    if len(matches) != 1:
        raise AssertionError("expected one CLI declaration for " + option)
    return matches[0]


def _flags_dictionary(tree):
    matches = [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and isinstance(node.value, ast.Dict)
        and any(
            isinstance(target, ast.Name) and target.id == "flags"
            for target in node.targets
        )
    ]
    if len(matches) != 1:
        raise AssertionError("expected one flags dictionary")
    return {
        ast.literal_eval(key): value
        for key, value in zip(matches[0].keys, matches[0].values)
    }


def _source():
    return "N" * SOURCE_START + PHASE_FIXTURE + "N" * 20


def _occurrence_signatures(search):
    return [
        (
            tuple(candidate.list_absolute_start),
            tuple(candidate.list_repeats),
            tuple(candidate.list_spacers),
            tuple(tuple(value) for value in candidate.list_fuzzy_counts),
            candidate.occurrence_path_id,
            candidate.occurrence_path_kind,
        )
        for candidate in search.complete_occurrence_path_family()
    ]


class BoundaryPathCliTests(unittest.TestCase):
    def test_cli_is_explicit_boolean_default_true_and_enters_flags(self):
        tree = _cli_tree()
        call = _argument_call(tree, "--fuzzy_boundary_paths")
        keywords = {keyword.arg: keyword.value for keyword in call.keywords}

        self.assertEqual("parse_boolean", keywords["type"].id)
        self.assertIs(ast.literal_eval(keywords["default"]), True)
        self.assertEqual((True, False), ast.literal_eval(keywords["choices"]))

        flag_value = _flags_dictionary(tree)["flag_fuzzy_boundary_paths"]
        self.assertIsInstance(flag_value, ast.Attribute)
        self.assertEqual("fuzzy_boundary_paths", flag_value.attr)
        self.assertEqual("args", flag_value.value.id)

    def test_cli_boolean_parser_rejects_ambiguous_values(self):
        self.assertIs(parse_boolean("True"), True)
        self.assertIs(parse_boolean("false"), False)
        for value in ("1", "yes", "", 0, 1, None):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    parse_boolean(value)


class BoundaryPathDetectionAndWorkerTests(unittest.TestCase):
    @staticmethod
    def _detection(flag=mock.DEFAULT):
        flags = {
            "flag_parallel": False,
            "flag_cpu": "1",
            "flag_fast_run": True,
            "flag_enhancement_min_max": True,
            "flag_enhancement_start_end": True,
        }
        if flag is not mock.DEFAULT:
            flags["flag_fuzzy_boundary_paths"] = flag
        parameters = {
            "param_fuzzy_search_mode": "native",
            "param_fuzzy_search_flank": 100,
            "param_fuzzy_family_padding": 15,
        }
        with mock.patch.object(Detection, "_get_complete_dna"), mock.patch.object(
            Detection, "_run_cluster_detection"
        ), mock.patch.object(
            Detection, "_extract_cluster_sequences"
        ), mock.patch.object(Detection, "_run_array_detection"):
            return Detection("unused.fna", flags, parameters, False)

    def test_detection_defaults_true_accepts_false_and_rejects_non_bool(self):
        self.assertIs(self._detection().flag_fuzzy_boundary_paths, True)
        self.assertIs(self._detection(False).flag_fuzzy_boundary_paths, False)
        for value in (0, 1, None, "False"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "must be boolean"):
                    self._detection(value)

    def test_worker_initializer_propagates_default_and_false(self):
        source = "NNNNNACGTNN"
        worker_input = ("AC", "ACGT", 5, ERROR_LIMIT)
        fake_search = mock.Mock(return_value=mock.sentinel.result)

        with mock.patch.object(module_detection, "FuzzySearch", fake_search):
            module_detection._initialize_fuzzy_source(source)
            self.assertIs(
                Detection._parallel_run_fuzzy_run(worker_input),
                mock.sentinel.result,
            )
            self.assertIs(
                fake_search.call_args.kwargs["_build_boundary_paths"], True
            )

            module_detection._initialize_fuzzy_source(source, False)
            Detection._parallel_run_fuzzy_run(worker_input)
            self.assertIs(
                fake_search.call_args.kwargs["_build_boundary_paths"], False
            )

        for value in (0, 1, None, "False"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(TypeError, "must be boolean"):
                    module_detection._initialize_fuzzy_source(source, value)


class BoundaryPathFuzzySearchTests(unittest.TestCase):
    def test_default_matches_explicit_true_and_calls_boundary_builder(self):
        original = FuzzySearch._initialize_boundary_path_alternatives
        calls = []

        def recording_builder(search, weighted_error, source_dna):
            calls.append(search)
            return original(search, weighted_error, source_dna)

        with mock.patch.object(
            FuzzySearch,
            "_initialize_boundary_path_alternatives",
            recording_builder,
        ):
            default = FuzzySearch(
                PHASE_FIXTURE,
                SOURCE_START,
                REPEAT,
                ERROR_LIMIT,
                source_dna=_source(),
            )
            explicit = FuzzySearch(
                PHASE_FIXTURE,
                SOURCE_START,
                REPEAT,
                ERROR_LIMIT,
                source_dna=_source(),
                _build_boundary_paths=True,
            )

        self.assertEqual(2, len(calls))
        self.assertEqual(
            _occurrence_signatures(default),
            _occurrence_signatures(explicit),
        )
        self.assertTrue(default._build_boundary_paths)
        self.assertTrue(all(
            alternative._build_boundary_paths
            for alternative in default.occurrence_path_alternatives
        ))

    def test_false_bypasses_only_boundary_builder_and_inherits_recursively(self):
        with mock.patch.object(
            FuzzySearch,
            "_initialize_boundary_path_alternatives",
        ) as boundary_builder:
            disabled = FuzzySearch(
                PHASE_FIXTURE,
                SOURCE_START,
                REPEAT,
                ERROR_LIMIT,
                source_dna=_source(),
                _build_boundary_paths=False,
            )

        boundary_builder.assert_not_called()
        self.assertIsNotNone(disabled.occurrence_path_lattice)
        self.assertEqual(1, len(disabled.occurrence_path_alternatives))
        self.assertEqual(2, len(disabled.complete_occurrence_path_family()))
        self.assertFalse(disabled._build_boundary_paths)
        self.assertTrue(all(
            alternative._build_boundary_paths is False
            for alternative in disabled.occurrence_path_alternatives
        ))
        self.assertIsNone(disabled.boundary_path_lattice)
        self.assertEqual((), disabled.boundary_path_alternatives)

    def test_constructor_rejects_non_boolean_boundary_flag(self):
        for value in (0, 1, None, "False"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(TypeError, "must be boolean"):
                    FuzzySearch(
                        PHASE_FIXTURE,
                        SOURCE_START,
                        REPEAT,
                        ERROR_LIMIT,
                        source_dna=_source(),
                        _build_boundary_paths=value,
                    )


if __name__ == "__main__":
    unittest.main()
