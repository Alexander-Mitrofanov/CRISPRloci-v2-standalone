"""Commands must not leave native descendants running after their launcher exits."""
import contextlib
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import pytest

from crisprtracrrna.errors import ToolError
from crisprtracrrna.process import Runner


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="Linux process-state control")
@pytest.mark.parametrize("exit_mode", ["success", "failure", "timeout"])
def test_native_descendants_are_stopped_on_every_exit(tmp_path, exit_mode):
    pid_file = tmp_path / "descendant.pid"
    script = (
        "import subprocess,pathlib,sys,time; "
        "child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)'], "
        "stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL); "
        f"pathlib.Path({str(pid_file)!r}).write_text(str(child.pid)); "
        + ("time.sleep(30)" if exit_mode == "timeout" else
           "sys.exit(7)" if exit_mode == "failure" else "sys.exit(0)")
    )
    try:
        if exit_mode == "success":
            Runner(tmp_path / "logs", timeout=0.5).run([sys.executable, "-c", script])
        else:
            with pytest.raises(ToolError):
                Runner(tmp_path / "logs", timeout=0.5).run([sys.executable, "-c", script])
        pid = int(pid_file.read_text())
        state = Path(f"/proc/{pid}/stat")
        for _ in range(40):
            if not state.exists() or state.read_text().split()[2] == "Z":
                break
            time.sleep(0.025)
        assert not state.exists() or state.read_text().split()[2] == "Z"
    finally:
        if pid_file.exists():
            with contextlib.suppress(ProcessLookupError):
                os.kill(int(pid_file.read_text()), signal.SIGKILL)


@pytest.mark.skipif(not Path("/proc").is_dir(), reason="Linux process-state control")
def test_sigterm_reaps_native_descendants(tmp_path):
    pid_file = tmp_path / "native.pid"
    native_script = (
        "import subprocess,pathlib,sys,time; "
        "child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(30)']); "
        f"pathlib.Path({str(pid_file)!r}).write_text(str(child.pid)); time.sleep(30)"
    )
    wrapper_script = (
        "from pathlib import Path; import sys; from crisprtracrrna.process import Runner; "
        f"command=[sys.executable,'-c',{native_script!r}]; "
        f"output=Path({str(tmp_path / 'wrapper-output')!r}); "
        "Runner(output, timeout=30).run(command)"
    )
    wrapper = subprocess.Popen([sys.executable, "-c", wrapper_script], start_new_session=True)
    try:
        deadline = time.monotonic() + 5
        while not pid_file.exists() and time.monotonic() < deadline:
            assert wrapper.poll() is None
            time.sleep(0.01)
        assert pid_file.exists()
        pid = int(pid_file.read_text())
        wrapper.send_signal(signal.SIGTERM)
        assert wrapper.wait(timeout=5) != 0
        state = Path(f"/proc/{pid}/stat")
        for _ in range(40):
            if not state.exists() or state.read_text().split()[2] == "Z":
                break
            time.sleep(0.025)
        assert not state.exists() or state.read_text().split()[2] == "Z"
    finally:
        with contextlib.suppress(ProcessLookupError):
            os.killpg(wrapper.pid, signal.SIGKILL)
        wrapper.wait()
        if pid_file.exists():
            with contextlib.suppress(ProcessLookupError):
                os.kill(int(pid_file.read_text()), signal.SIGKILL)


def test_native_command_is_usable_from_a_worker_thread(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    before = signal.getsignal(signal.SIGTERM)
    command = [sys.executable, "-c", "print('finished')"]
    with ThreadPoolExecutor(max_workers=1) as pool:
        pool.submit(
            lambda: Runner(tmp_path / "thread-logs", timeout=3).run(command)
        ).result(timeout=5)
    assert signal.getsignal(signal.SIGTERM) == before
