"""Containment tests use actual subprocesses, including forked descendants."""

from pathlib import Path
import os
import sys
import threading
import time

import pytest

from crisprleader.errors import ExecutionError
from crisprleader.runtime import run_command


pytestmark = pytest.mark.skipif(os.name != "posix", reason="POSIX process groups required")


def python(code: str, tmp_path: Path, **kwargs: object):
    return run_command([sys.executable, "-c", code], cwd=tmp_path, **kwargs)


def test_no_shell_and_private_cwd(tmp_path: Path):
    literal = "$(touch SHOULD_NOT_EXIST); * | &"
    result = run_command(
        [sys.executable, "-c", "import os,sys;print(os.getcwd());print(sys.argv[1])", literal],
        cwd=tmp_path,
    )
    assert result.stdout == f"{tmp_path}\n{literal}\n"
    assert result.stderr == "" and result.returncode == 0
    assert not (tmp_path / "SHOULD_NOT_EXIST").exists()


def test_environment_is_explicit_replacement(tmp_path: Path):
    result = python("import os; print(os.environ.get('ONLY_THIS'));print(os.environ.get('HOME'))",
                    tmp_path, env={"ONLY_THIS": "yes"})
    assert result.stdout == "yes\nNone\n"


def test_reads_both_pipes_without_deadlock(tmp_path: Path):
    result = python("import os;os.write(2,b'e'*150000);os.write(1,b'o'*150000)", tmp_path)
    assert result.stderr == "e" * 150000
    assert result.stdout == "o" * 150000


@pytest.mark.parametrize("code", ["import os;os.write(1,b'x'*101)",
                                  "import os;os.write(2,b'x'*101)",
                                  "import os;os.write(1,b'x'*60);os.write(2,b'x'*41)"])
def test_aggregate_byte_limit(tmp_path: Path, code: str):
    with pytest.raises(ExecutionError, match="combined stdout/stderr"):
        python(code, tmp_path, max_output_bytes=100)


def test_exact_limit_and_invalid_utf8(tmp_path: Path):
    result = python("import os;os.write(1,b'\\xff'*100)", tmp_path, max_output_bytes=100)
    assert result.stdout == "\ufffd" * 100
    assert python("pass", tmp_path, max_output_bytes=0).stdout == ""


@pytest.mark.parametrize("code", ["import time;time.sleep(20)",
                                  "import os,time;os.close(1);os.close(2);time.sleep(20)"])
def test_timeout_with_open_or_closed_pipes(tmp_path: Path, code: str):
    started = time.monotonic()
    with pytest.raises(ExecutionError, match="timed out"):
        python(code, tmp_path, timeout_seconds=0.15)
    assert time.monotonic() - started < 3


def test_cancellation(tmp_path: Path):
    cancelled = threading.Event()
    timer = threading.Timer(0.15, cancelled.set)
    timer.start()
    try:
        with pytest.raises(ExecutionError, match="cancelled"):
            python("import time;time.sleep(20)", tmp_path, cancel_event=cancelled)
    finally:
        timer.join()
    with pytest.raises(ExecutionError, match="before launch"):
        python("raise RuntimeError('should never launch')", tmp_path, cancel_event=cancelled)


def assert_not_running(pid: int):
    for _ in range(50):
        try:
            stat = Path(f"/proc/{pid}/stat").read_text()
        except FileNotFoundError:
            return
        if stat.split(") ", 1)[1].startswith("Z"):
            return  # Orphan awaiting OS/init reaping cannot execute.
        time.sleep(0.02)
    pytest.fail(f"descendant {pid} still executes after cleanup")


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="Linux /proc is used for independent check")
@pytest.mark.parametrize("close_pipes", [False, True])
def test_descendant_terminated_when_leader_exits(tmp_path: Path, close_pipes: bool):
    code = (
        "import os,time\n"
        "child=os.fork()\n"
        "if child:\n"
        " open('child.pid','w').write(str(child))\n"
        " os._exit(0)\n"
        + ("os.close(1);os.close(2)\n" if close_pipes else "")
        + "time.sleep(20)\n"
    )
    if close_pipes:
        assert python(code, tmp_path, timeout_seconds=0.5).returncode == 0
    else:
        with pytest.raises(ExecutionError, match="timed out"):
            python(code, tmp_path, timeout_seconds=0.15)
    assert_not_running(int((tmp_path / "child.pid").read_text()))


def test_failed_exit_is_actionable_and_bounded(tmp_path: Path):
    with pytest.raises(ExecutionError, match="status 7") as caught:
        python("import sys;sys.stderr.write('x'*10000);sys.exit(7)", tmp_path)
    assert len(str(caught.value)) < 4200


@pytest.mark.parametrize("kwargs", [{"timeout_seconds": 0}, {"timeout_seconds": float('nan')},
                                    {"timeout_seconds": True}, {"max_output_bytes": -1},
                                    {"max_output_bytes": True}])
def test_bad_bounds_rejected(tmp_path: Path, kwargs: dict):
    with pytest.raises(ExecutionError):
        python("pass", tmp_path, **kwargs)


def test_launch_errors(tmp_path: Path):
    with pytest.raises(ExecutionError, match="Could not start"):
        run_command(["/does-not-exist/crisprleader"], cwd=tmp_path)
    with pytest.raises(ExecutionError, match="argument sequence"):
        run_command("echo unsafe", cwd=tmp_path)
    with pytest.raises(ExecutionError, match="working directory"):
        python("pass", tmp_path / "missing")
