import importlib.util
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCHEDULER_PATH = (
    ROOT / "components" / "folder_scheduler.py"
)


def _load_v2_scheduler():
    spec = importlib.util.spec_from_file_location(
        "crispridentify_v2_folder_scheduler", SCHEDULER_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCHEDULER = _load_v2_scheduler()


def _option_value(arguments, option):
    index = arguments.index(option)
    return arguments[index + 1]


class CasAndraFolderLanePropagationTests(unittest.TestCase):
    def test_two_lane_launch_preserves_casandra_options_and_cpu_budgets(self):
        parent_arguments = (
            "--input_folder", "/inputs",
            "--result_folder", "/results",
            "--cas", "True",
            "--cas_gene_mode", "meta",
            "--cas_timeout_seconds", "901",
            "--cpu", "3",
            "--folder_jobs", "2",
            "--parallel", "False",
        )
        cpu_lanes = SCHEDULER.partition_cpu_lanes((2, 5, 8), 2)
        popen = mock.Mock()

        for lane_id, cpu_ids in enumerate(cpu_lanes, 1):
            child_arguments = SCHEDULER.rewrite_lane_argv(
                parent_arguments,
                input_folder="/queue/claimed/lane_{:02d}".format(lane_id),
                result_folder="/results",
                cpu=len(cpu_ids),
                folder_jobs=1,
            ) + (
                "--_folder_queue", "/queue",
                "--_folder_lane_id", str(lane_id),
                "--_folder_cpu_ids", ",".join(str(cpu) for cpu in cpu_ids),
            )
            popen(
                ("python", "CRISPRidentify.py") + child_arguments,
                env={"OMP_NUM_THREADS": str(len(cpu_ids))},
            )

        self.assertEqual(cpu_lanes, ((2, 5), (8,)))
        self.assertEqual(popen.call_count, 2)
        expected_lane_cpus = (("2", "2,5"), ("1", "8"))
        for call, (expected_cpu, expected_cpu_ids) in zip(
                popen.call_args_list, expected_lane_cpus):
            command = call.args[0]
            self.assertEqual(command.count("--cas"), 1)
            self.assertEqual(_option_value(command, "--cas"), "True")
            self.assertEqual(command.count("--cas_gene_mode"), 1)
            self.assertEqual(_option_value(command, "--cas_gene_mode"), "meta")
            self.assertEqual(command.count("--cas_timeout_seconds"), 1)
            self.assertEqual(
                _option_value(command, "--cas_timeout_seconds"), "901"
            )
            self.assertEqual(command.count("--cpu"), 1)
            self.assertEqual(_option_value(command, "--cpu"), expected_cpu)
            self.assertEqual(
                _option_value(command, "--_folder_cpu_ids"),
                expected_cpu_ids,
            )
            self.assertEqual(call.kwargs["env"]["OMP_NUM_THREADS"], expected_cpu)


if __name__ == "__main__":
    unittest.main()

