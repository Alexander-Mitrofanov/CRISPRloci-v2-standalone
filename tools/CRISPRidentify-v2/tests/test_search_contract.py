import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from components.detection.search import VmatchRun
from components.detection.search_contract import effective_discovery
from components.refinement.candidates import AdvancedFuzzySearchFilter, CrisprCandidate
from components.sequence.validation import read_single_fasta
from tests.test_fuzzy_lattice_storage import PARAMETERS

ROOT = Path(__file__).resolve().parents[1]


class SearchContractTests(unittest.TestCase):
    def test_wider_structural_spacer_limit_does_not_expand_seed_query(self):
        source_path = ROOT / "tests/data/template_discovery/110nt_spacers.fa"
        _, source = read_single_fasta(source_path)
        starts = [81, 221, 361]
        repeats = [source[s:s+30] for s in starts]
        spacers = [source[s+30:t] for s, t in zip(starts, starts[1:])]
        candidate = CrisprCandidate(repeats, repeats, spacers, starts)
        self.assertFalse(AdvancedFuzzySearchFilter.from_parameters(PARAMETERS)(candidate))
        self.assertTrue(AdvancedFuzzySearchFilter.from_parameters(
            dict(PARAMETERS, param_max_avg_spacer_length=120))(candidate))
        for maximum in (78, 120):
            config = json.loads(subprocess.check_output([
                str(ROOT / "crispridentify"), "--print-effective-config", "--max_len_spacer", str(maximum)], text=True))
            self.assertEqual(maximum, config["effective_config"]["max_len_spacer"])
            self.assertEqual([23, -10, 100], config["effective_config"]["effective_discovery"]["vmatch"]["-l"])
        # Use the real Vmatch executable with an isolated working directory.
        code = """import sys
from components.detection.search import VmatchRun
assert not VmatchRun(sys.argv[1], True, 'balanced').output()
assert not VmatchRun(sys.argv[1], True, 'sensitive').output()
"""
        import os
        with tempfile.TemporaryDirectory() as directory:
            env = dict(os.environ, PYTHONPATH=str(ROOT),
                       PATH=str(ROOT / ".envs/crispridentify-core/bin") + os.pathsep + os.environ["PATH"])
            result = subprocess.run([sys.executable, "-c", code, str(source_path)], cwd=directory,
                                    env=env, capture_output=True, text=True)
            self.assertEqual(0, result.returncode, result.stderr)

    def test_captured_command_matches_report_for_all_modes(self):
        for fast in (True, False):
            for profile, fast_values in (("balanced", [23, -10, 100]),
                                         ("legacy", [23, 25, 60]), ("sensitive", [21, -10, 100])):
                expected = fast_values if fast else [21, -10, 100]
                run = VmatchRun.__new__(VmatchRun)
                run.flag_fast_run, run.fast_seed_profile = fast, profile
                run._get_run_settings()
                with mock.patch("builtins.open", mock.mock_open()), mock.patch("subprocess.run") as command:
                    run._vmatch_command()
                argv = command.call_args.args[0]
                self.assertEqual(["vmatch", "-l"] + list(map(str, expected)), argv[:5])
                self.assertEqual("1000001", argv[argv.index("-best") + 1])
                result = subprocess.run([sys.executable, str(ROOT / "CRISPRidentify.py"),
                                         "--print-effective-config", "--fast_run", str(fast),
                                         "--fast_run_seed_profile", profile, "--fuzzy_search_flank", "73"],
                                        capture_output=True, text=True, check=True)
                config = json.loads(result.stdout)["effective_config"]["effective_discovery"]
                self.assertEqual(expected, config["vmatch"]["-l"])
                self.assertEqual(config, effective_discovery(fast, profile, 73, "native"))
                self.assertEqual(55, config["seed_text_truncation"])
                self.assertEqual(140, config["local_chain_max_gap"])

    def test_contract_import_is_dependency_light(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run([sys.executable, "-c",
                "import sys; sys.path.insert(0, sys.argv[1]); "
                "import components.detection.search_contract; "
                "import components.detection.repeat_templates; "
                "assert not {'numpy', 'regex', 'components.pipeline'} & set(sys.modules)", str(ROOT)],
                cwd=directory, capture_output=True, text=True)
        self.assertEqual(0, result.returncode, result.stderr)
