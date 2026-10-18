import importlib.util
from pathlib import Path
import sys
import types
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT
    / "vendor"
    / "casandra"
    / "src"
    / "casandra"
    / "classification_v4.py"
)
CASSETTE_MODULE_PATH = (
    ROOT
    / "vendor"
    / "casandra"
    / "src"
    / "casandra"
    / "cassette_hybrid_v2.py"
)


def load_module():
    spec = importlib.util.spec_from_file_location(
        "casandra_classification_v4_under_test", MODULE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


classification_v4 = load_module()


def load_cassette_module(architecture_result=("III-D", 0.75)):
    package = types.ModuleType("casandra")
    package.__path__ = []
    architecture = types.ModuleType("casandra.architecture_runtime")
    architecture.classify_architecture_local = mock.Mock(
        return_value=architecture_result
    )
    spec = importlib.util.spec_from_file_location(
        "casandra_cassette_hybrid_v2_under_test", CASSETTE_MODULE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(
        sys.modules,
        {
            "casandra": package,
            "casandra.architecture_runtime": architecture,
            "casandra.classification_v4": classification_v4,
            spec.name: module,
        },
    ):
        spec.loader.exec_module(module)
    return module, architecture.classify_architecture_local


class CasAndraClassificationV4Tests(unittest.TestCase):
    def test_compatible_hierarchy_is_unchanged(self):
        raw = {"class": "1", "type": "I", "subtype": "I-F1"}
        normalized, masked = classification_v4.normalize_classification(raw)
        self.assertEqual(raw, normalized)
        self.assertFalse(masked)

    def test_real_type_i_type_v_conflicts_are_masked(self):
        for conflicting_subtype in ("V-M", "V-F"):
            with self.subTest(subtype=conflicting_subtype):
                raw = {
                    "class": "1",
                    "type": "I",
                    "subtype": conflicting_subtype,
                }
                normalized, masked = (
                    classification_v4.normalize_classification(raw)
                )
                self.assertEqual(
                    {"class": "1", "type": "I", "subtype": None},
                    normalized,
                )
                self.assertTrue(masked)

    def test_type_determines_class(self):
        raw = {"class": "1", "type": "V", "subtype": "V-F"}
        normalized, masked = classification_v4.normalize_classification(raw)
        self.assertEqual(
            {"class": "2", "type": "V", "subtype": "V-F"},
            normalized,
        )
        self.assertTrue(masked)

    def test_negative_null_classification_is_stable(self):
        raw = {"class": None, "type": None, "subtype": None}
        normalized, masked = classification_v4.normalize_classification(raw)
        self.assertEqual(raw, normalized)
        self.assertFalse(masked)


class CasAndraCassetteHierarchyTests(unittest.TestCase):
    def test_real_mixed_class1_type_vi_route_is_normalized(self):
        module, architecture = load_cassette_module()
        direct = {
            "class": "1",
            "type": "VI",
            "subtype": "VI-B",
            "confidence": {
                "class_vote_fraction": 0.6489246084062577,
                "type_vote_fraction": 0.35107539159374235,
                "subtype_vote_fraction": 1.0,
            },
        }
        final = module.classify_hybrid_cassette_v2(
            {"classification": direct}, None, {}, {}
        )

        self.assertEqual(
            {"class": "2", "type": "VI", "subtype": "VI-B"},
            {key: final[key] for key in ("class", "type", "subtype")},
        )
        self.assertEqual(
            {
                "raw_classification": {
                    "class": "1", "type": "VI", "subtype": "VI-B",
                },
                "hierarchy_conflict_masked": True,
            },
            final["evidence"],
        )
        self.assertIs(final["direct_profile_result"], direct)
        self.assertEqual("direct_profile_aggregation", final["method"])
        self.assertEqual(1.0, final["confidence"])
        architecture.assert_not_called()

    def test_type_iii_architecture_route_normalizes_raw_class(self):
        module, architecture = load_cassette_module(("III-D", 0.625))
        direct = {
            "class": "2",
            "type": "III",
            "subtype": None,
            "confidence": {
                "class_vote_fraction": 0.51,
                "type_vote_fraction": 0.49,
                "subtype_vote_fraction": 0.0,
            },
        }
        cassette = {"classification": direct}
        genes = {"p1": object()}
        predictions = {"p1": object()}
        final = module.classify_hybrid_cassette_v2(
            cassette, "architecture", genes, predictions
        )

        architecture.assert_called_once_with(
            "architecture", cassette, genes, predictions
        )
        self.assertEqual(
            {"class": "1", "type": "III", "subtype": "III-D"},
            {key: final[key] for key in ("class", "type", "subtype")},
        )
        self.assertEqual(
            {
                "raw_classification": {
                    "class": "2", "type": "III", "subtype": "III-D",
                },
                "hierarchy_conflict_masked": True,
            },
            final["evidence"],
        )
        self.assertEqual(
            "type_iii_ordered_profile_architecture_extratrees",
            final["method"],
        )
        self.assertEqual(0.625, final["confidence"])


if __name__ == "__main__":
    unittest.main()
