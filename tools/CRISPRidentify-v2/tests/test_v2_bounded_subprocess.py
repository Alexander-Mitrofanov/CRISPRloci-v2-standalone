import importlib.util
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "components" / "bounded_subprocess.py"
)


def _load_module():
    name = "crispridentify_v2_bounded_subprocess"
    spec = importlib.util.spec_from_file_location(name, MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BOUNDED = _load_module()


class BoundedSubprocessTests(unittest.TestCase):
    def test_success_returns_decoded_stdout_stderr_and_status(self):
        completed = BOUNDED.run_bounded(
            [
                sys.executable,
                "-c",
                (
                    "import sys; "
                    "sys.stdout.buffer.write('result ✓'.encode()); "
                    "sys.stderr.buffer.write(b'diagnostic')"
                ),
            ],
            timeout=5,
            stdout_limit=1024,
            stderr_limit=1024,
        )

        self.assertEqual(completed.returncode, 0)
        self.assertEqual(completed.stdout, "result ✓")
        self.assertEqual(completed.stderr, "diagnostic")
        self.assertFalse(completed.stdout_truncated)
        self.assertFalse(completed.stderr_truncated)

    def test_oversized_stdout_is_typed_error_while_stderr_keeps_tail(self):
        script = (
            "import sys\n"
            "for _ in range(32):\n"
            "    sys.stdout.buffer.write(b'O' * 8192)\n"
            "    sys.stderr.buffer.write(b'E' * 8192)\n"
            "sys.stderr.buffer.write(b'USEFUL-TAIL')\n"
        )
        with self.assertRaises(
                BOUNDED.BoundedSubprocessTruncationError) as raised:
            BOUNDED.run_bounded(
                [sys.executable, "-c", script],
                timeout=10,
                stdout_limit=4096,
                stderr_limit=64,
            )

        error = raised.exception
        self.assertEqual(error.returncode, 0)
        self.assertTrue(error.stdout_truncated)
        self.assertEqual(error.stdout, "O" * 4096)
        self.assertTrue(error.stderr_truncated)
        self.assertLessEqual(len(error.stderr.encode("utf-8")), 64)
        self.assertTrue(error.stderr.endswith("USEFUL-TAIL"))

    def test_oversized_stderr_alone_returns_its_useful_tail(self):
        completed = BOUNDED.run_bounded(
            [
                sys.executable,
                "-c",
                (
                    "import sys; "
                    "sys.stdout.write('ok'); "
                    "sys.stderr.write('old-' * 1000 + 'FINAL-DIAGNOSTIC')"
                ),
            ],
            timeout=5,
            stdout_limit=32,
            stderr_limit=48,
        )

        self.assertEqual(completed.stdout, "ok")
        self.assertTrue(completed.stderr_truncated)
        self.assertLessEqual(len(completed.stderr.encode("utf-8")), 48)
        self.assertTrue(completed.stderr.endswith("FINAL-DIAGNOSTIC"))

    @unittest.skipUnless(os.name == "posix", "process-group test requires POSIX")
    def test_timeout_kills_spawned_child_and_returns_without_hanging(self):
        with tempfile.TemporaryDirectory() as temporary_directory:
            sentinel = Path(temporary_directory) / "escaped-child.txt"
            child_script = (
                "import pathlib, sys, time; "
                "time.sleep(0.8); "
                "pathlib.Path(sys.argv[1]).write_text('escaped')"
            )
            parent_script = (
                "import subprocess, sys, time\n"
                "subprocess.Popen([sys.executable, '-c', sys.argv[1], sys.argv[2]])\n"
                "print('child-spawned', flush=True)\n"
                "time.sleep(60)\n"
            )
            started = time.monotonic()
            with self.assertRaises(
                    BOUNDED.BoundedSubprocessTimeoutError) as raised:
                BOUNDED.run_bounded(
                    [
                        sys.executable,
                        "-c",
                        parent_script,
                        child_script,
                        str(sentinel),
                    ],
                    timeout=0.25,
                    stdout_limit=1024,
                    stderr_limit=1024,
                )
            elapsed = time.monotonic() - started

            self.assertLess(elapsed, 5)
            self.assertIn("child-spawned", raised.exception.stdout)
            self.assertLess(raised.exception.returncode, 0)
            time.sleep(0.9)
            self.assertFalse(sentinel.exists())

    def test_rejects_shell_strings_and_invalid_limits_before_launch(self):
        with self.assertRaises(TypeError):
            BOUNDED.run_bounded("echo unsafe", timeout=1)
        for timeout in (0, -1, float("inf"), True):
            with self.subTest(timeout=timeout):
                with self.assertRaises(ValueError):
                    BOUNDED.run_bounded([sys.executable], timeout=timeout)
        with self.assertRaises(ValueError):
            BOUNDED.run_bounded(
                [sys.executable], timeout=1, stdout_limit=0
            )


if __name__ == "__main__":
    unittest.main()

