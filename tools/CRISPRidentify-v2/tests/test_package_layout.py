"""Public import, serialization, worker-startup, and asset-location contracts."""

import ast
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from components._compat import LEGACY_MODULES
from components.cli.parser import build_parser

ROOT = Path(__file__).resolve().parents[1]


class PackageLayoutTests(unittest.TestCase):
    def run_probe(self, source, *arguments):
        # A fresh process prevents earlier tests from hiding import-order bugs.
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run(
                [sys.executable, "-c", textwrap.dedent(source), str(ROOT), *arguments],
                cwd=directory,
                env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
                capture_output=True,
                text=True,
                timeout=60,
            )
        self.assertEqual(0, result.returncode, result.stdout + result.stderr)

    def test_legacy_imports_share_identity_in_both_import_orders(self):
        probe = """
            import importlib, sys
            from unittest.mock import patch
            sys.path.insert(0, sys.argv[1])
            import components
            from components._compat import LEGACY_MODULES
            for old, new in LEGACY_MODULES.items():
                old_name, new_name = 'components.' + old, 'components.' + new
                names = (old_name, new_name)
                if sys.argv[2] == 'canonical-first':
                    names = names[::-1]
                first, second = map(importlib.import_module, names)
                assert first is second, (old_name, new_name)
                assert getattr(components, old) is first
                assert first.__name__ == new_name, first.__name__
                assert first.__spec__.name == new_name, first.__spec__
                assert first.__file__.endswith(new.replace('.', '/') + '.py')
                with patch(old_name + '._alias_probe', True, create=True):
                    assert second._alias_probe
                assert not hasattr(second, '_alias_probe')
        """
        for order in ("legacy-first", "canonical-first"):
            with self.subTest(order=order):
                self.run_probe(probe, order)

    def test_old_pickle_class_references_load_and_new_instances_round_trip(self):
        self.run_probe("""
            import pickle, sys
            sys.path.insert(0, sys.argv[1])
            # Protocol-0 GLOBAL references reproduce pre-refactor pickle names.
            candidate_type = pickle.loads(
                b'ccomponents.components_detection_refinement\\nCrisprCandidate\\n.'
            )
            from components.refinement.candidates import CrisprCandidate
            assert candidate_type is CrisprCandidate
            from components.pipeline import CandidateScoreCache as old_cache
            from components.scoring.cache import CandidateScoreCache
            assert old_cache is CandidateScoreCache
            candidate = CrisprCandidate(
                list_repeats=['ACGT', 'ACGT'],
                list_repeats_gaped=['ACGT', 'ACGT'],
                list_spacers=['TTTT'],
                list_repeat_starts=[0, 8],
            )
            restored = pickle.loads(pickle.dumps(candidate))
            assert type(restored) is CrisprCandidate
            assert restored == candidate
        """)

    def test_parser_and_package_import_do_not_load_scientific_runtime(self):
        self.run_probe("""
            import sys
            sys.path.insert(0, sys.argv[1])
            from components.cli.parser import build_parser
            assert build_parser().parse_args([]).preset == 'balanced'
            for name in ('numpy', 'components.pipeline', 'components.detection.stage'):
                assert name not in sys.modules, name
        """)

    def test_parser_instances_do_not_share_mutable_defaults(self):
        first, second = build_parser(), build_parser()
        first.set_defaults(preset="full")
        self.assertEqual("full", first.parse_args([]).preset)
        self.assertEqual("balanced", second.parse_args([]).preset)

    def test_assets_resolve_outside_repository_working_directory(self):
        self.run_probe("""
            from pathlib import Path
            import sys
            root = Path(sys.argv[1])
            sys.path.insert(0, str(root))
            from components.paths import PROJECT_ROOT, MODEL_ROOT, NATIVE_LIBRARY
            from components.scoring.compact_model import MANIFEST_PATH
            from components.scoring.compact_feature_extraction import BLAST_DATABASE
            from components.annotation.portable_strand import DEFAULT_MODEL_DIRECTORY
            from components.detection.literal_fuzzy import _LIBRARY_PATH
            from components.scoring.sequence_path_vectorizer import _NATIVE_LIBRARY_PATH
            assert PROJECT_ROOT == root
            assert MODEL_ROOT == root / 'trained_models'
            assert MANIFEST_PATH.is_file()
            assert Path(BLAST_DATABASE).is_file()
            assert Path(DEFAULT_MODEL_DIRECTORY).is_dir()
            assert NATIVE_LIBRARY.parent == root / 'components' / 'native'
            assert _LIBRARY_PATH == _NATIVE_LIBRARY_PATH == NATIVE_LIBRARY
        """)

    def test_first_party_implementation_uses_canonical_imports(self):
        paths = [ROOT / "CRISPRidentify.py"]
        paths.extend((ROOT / "components").rglob("*.py"))
        paths.extend((ROOT / "scripts").glob("*.py"))
        for path in paths:
            for node in ast.walk(ast.parse(path.read_text())):
                names = []
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.module:
                    names = [node.module]
                    if node.module == "components":
                        names.extend("components." + alias.name for alias in node.names)
                for name in names:
                    parts = name.split(".")
                    if len(parts) >= 2 and parts[0] == "components":
                        self.assertNotIn(
                            parts[1],
                            LEGACY_MODULES,
                            f"{path.relative_to(ROOT)}:{node.lineno}: {name}",
                        )


if __name__ == "__main__":
    unittest.main()
