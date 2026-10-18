import argparse
import ast
import importlib.util
import sys
import types
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT
CLI_PATH = V2 / "CRISPRidentify.py"


def _load_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


CLI_UTILS = _load_file(
    "crispridentify_v2_timeout_cli_utils",
    V2 / "components" / "cli_utils.py",
)


def _load_non_array_module():
    components = types.ModuleType("components")
    components.__path__ = []

    computations = types.ModuleType(
        "components.components_non_array_computations"
    )
    computations.FullISElementSearch = object
    computations.FullLeaderSeqSearch = object
    computations.StrandComputation = object

    adapter = types.ModuleType("components.casandra_adapter")
    adapter.CasAndraPredictor = object

    reconstruction = types.ModuleType("components.reconstruction_validator")
    reconstruction.read_single_fasta = lambda _path: ("source", "ACGT")

    with mock.patch.dict(sys.modules, {
        "components": components,
        "components.components_non_array_computations": computations,
        "components.casandra_adapter": adapter,
        "components.reconstruction_validator": reconstruction,
    }):
        return _load_file(
            "crispridentify_v2_timeout_non_array",
            V2 / "components" / "module_non_array_computations.py",
        )


NON_ARRAY = _load_non_array_module()


def _cli_tree():
    return ast.parse(CLI_PATH.read_text(encoding="utf-8"), filename=str(CLI_PATH))


def _argument_call(tree, option):
    matches = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not node.args:
            continue
        function = node.func
        if not isinstance(function, ast.Attribute) or function.attr != "add_argument":
            continue
        if isinstance(node.args[0], ast.Constant) and node.args[0].value == option:
            matches.append(node)
    if len(matches) != 1:
        raise AssertionError(
            "expected one parser registration for {}, found {}".format(
                option, len(matches)
            )
        )
    return matches[0]


def _dict_assignment(tree, name):
    matches = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Dict):
            continue
        if any(
                isinstance(target, ast.Name) and target.id == name
                for target in node.targets):
            matches.append(node.value)
    if len(matches) != 1:
        raise AssertionError(
            "expected one dictionary assignment for {}, found {}".format(
                name, len(matches)
            )
        )
    return {
        ast.literal_eval(key): value
        for key, value in zip(matches[0].keys, matches[0].values)
    }


def _assert_args_timeout(test_case, node):
    test_case.assertIsInstance(node, ast.Attribute)
    test_case.assertEqual(node.attr, "cas_timeout_seconds")
    test_case.assertIsInstance(node.value, ast.Name)
    test_case.assertEqual(node.value.id, "args")


class PositiveIntegerTests(unittest.TestCase):
    def test_accepts_positive_integer_text(self):
        self.assertEqual(CLI_UTILS.positive_integer("1"), 1)
        self.assertEqual(CLI_UTILS.positive_integer("3600"), 3600)

    def test_rejects_nonpositive_and_noninteger_values(self):
        for value in ("0", "-1", "1.5", "many", "", 0, -1, 1.5, True):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    CLI_UTILS.positive_integer(value)


class TimeoutCliContractTests(unittest.TestCase):
    def test_cli_declares_positive_one_hour_default(self):
        call = _argument_call(_cli_tree(), "--cas_timeout_seconds")
        keywords = {keyword.arg: keyword.value for keyword in call.keywords}

        self.assertIsInstance(keywords["type"], ast.Name)
        self.assertEqual(keywords["type"].id, "positive_integer")
        self.assertEqual(ast.literal_eval(keywords["default"]), 3600)
        help_text = ast.literal_eval(keywords["help"])
        self.assertIn("each CasAndra genome annotation", help_text)
        self.assertIn("process group", help_text)

    def test_timeout_enters_runtime_flags_and_report_parameters(self):
        tree = _cli_tree()
        flags = _dict_assignment(tree, "flags")
        parameters = _dict_assignment(tree, "parameters")

        _assert_args_timeout(self, flags["cas_timeout_seconds"])
        _assert_args_timeout(self, parameters["param_cas_timeout_seconds"])

    def test_readme_documents_default_and_override(self):
        readme = (V2 / "README.md").read_text(encoding="utf-8")
        self.assertIn("--cas_timeout_seconds", readme)
        self.assertIn("one-hour", readme)


class PredictorTimeoutPropagationTests(unittest.TestCase):
    def test_predictor_receives_configured_timeout_and_default(self):
        cases = (
            ({"cas_timeout_seconds": 901}, 901),
            ({}, 3600),
        )
        for timeout_flags, expected in cases:
            with self.subTest(timeout_flags=timeout_flags):
                predictor_instance = mock.Mock()
                predictor_instance.run.return_value = {"status": "ok"}
                predictor = mock.Mock(return_value=predictor_instance)
                computations = NON_ARRAY.NonArrayComputations.__new__(
                    NON_ARRAY.NonArrayComputations
                )
                computations.file_path = "/inputs/genome.fasta"
                computations.source_header = "genome description"
                computations.dna = "ACGT"
                computations.flags_non_arrays_computations = {
                    "flag_cpu": "3",
                    "cas_gene_mode": "meta",
                    **timeout_flags,
                }

                with mock.patch.object(
                        NON_ARRAY, "CasAndraPredictor", predictor):
                    computations._calculate_cas()

                predictor.assert_called_once_with(
                    input_fasta="/inputs/genome.fasta",
                    source_id="genome",
                    source_dna="ACGT",
                    cpu_count=3,
                    gene_mode="meta",
                    timeout_seconds=expected,
                )
                self.assertEqual(computations.cas_result, {"status": "ok"})


if __name__ == "__main__":
    unittest.main()

