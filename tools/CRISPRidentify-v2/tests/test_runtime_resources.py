import json
import os
from pathlib import Path
import subprocess
import sys
import unittest
from unittest import mock

from components import runtime_resources
from components.components_detection import FilterApproximation


ROOT = Path(__file__).resolve().parents[1]


class RuntimeResourcesTests(unittest.TestCase):
    def test_rejects_invalid_budget_without_changing_environment(self):
        before = dict(os.environ)
        for value in (0, -1, True, "2", 1.5, None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                runtime_resources.configure_native_threads(value)
        self.assertEqual(before, dict(os.environ))

    def test_cli_limits_loaded_blas_and_respects_lane_affinity(self):
        # Execute the actual CLI in a fresh interpreter; importing NumPy before
        # applying the policy would leave its native pool at the inherited 12.
        probe = '''
import contextlib, io, json, os, runpy, sys
sys.argv = ["CRISPRidentify.py", "--cpu", "2", "--print-effective-config"] + sys.argv[1:]
with contextlib.redirect_stdout(io.StringIO()):
    try:
        runpy.run_path("CRISPRidentify.py", run_name="__main__")
    except SystemExit as error:
        assert error.code in (0, None), error.code
from threadpoolctl import threadpool_info
from components.runtime_resources import external_thread_budget
print(json.dumps({"pools": threadpool_info(), "budget": external_thread_budget(),
                  "omp": os.environ["OMP_NUM_THREADS"]}))
'''
        allowed = sorted(os.sched_getaffinity(0))
        cases = [([], min(2, len(allowed)))]
        if hasattr(os, "sched_setaffinity"):
            cases.append((["--_folder_cpu_ids", str(allowed[-1])], 1))
        for options, expected in cases:
            environment = dict(os.environ, PYTHONHASHSEED="0",
                               OPENBLAS_NUM_THREADS="12", OMP_NUM_THREADS="12")
            result = subprocess.run(
                [sys.executable, "-c", probe, *options], cwd=ROOT,
                env=environment, check=True, capture_output=True, text=True,
            )
            state = json.loads(result.stdout)
            blas = [pool for pool in state["pools"] if pool["user_api"] == "blas"]
            self.assertTrue(blas)
            self.assertTrue(all(pool["num_threads"] == 1 for pool in blas))
            self.assertEqual(expected, state["budget"])
            self.assertEqual(str(expected), state["omp"])

    def test_clustal_receives_explicit_allocation(self):
        approximation = FilterApproximation.__new__(FilterApproximation)
        approximation.max_seq = "ACGT"
        approximation.min_seq = "ACG"
        with mock.patch(
                "components.components_detection.external_thread_budget",
                return_value=2), mock.patch(
                "components.components_detection.subprocess.run") as run, \
                mock.patch("builtins.open", mock.mock_open()):
            approximation._run_clustal_omega_repeats()
            approximation._run_clustal_omega_min_max()
        for call in run.call_args_list:
            argv = call.args[0]
            self.assertEqual("2", argv[argv.index("--threads") + 1])


if __name__ == "__main__":
    unittest.main()
