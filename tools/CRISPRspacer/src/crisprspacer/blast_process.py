"""Bounded, local-only BLAST subprocesses in private working directories."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time

from .io import InputError, LimitError, integer, regular


def blast_environment(cwd: Path):
    # Do not inherit .ncbirc, data loaders, telemetry or loader-injection variables.
    return {"PATH": "/usr/local/bin:/usr/bin:/bin", "LC_ALL": "C", "LANG": "C",
            "TMPDIR": str(cwd), "NCBI_DONT_USE_NCBIRC": "1",
            "NCBI_CONFIG__BLAST__DATA_LOADERS": "none", "BLAST_USAGE_REPORT": "false",
            "OMP_NUM_THREADS": "1"}


def _stop(process):
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    elif process.poll() is None:
        process.terminate()
    try:
        process.wait(timeout=1)
    except subprocess.TimeoutExpired:
        pass
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    elif process.poll() is None:
        process.kill()
    process.wait(timeout=5)


def run_blast(command: list[str], cwd: Path, stdout: Path, stderr: Path,
              timeout: int, max_output_bytes: int = 64 * 1024**2) -> None:
    integer(timeout, "native timeout", 1, 86_400)
    integer(max_output_bytes, "native output budget", 1, 20 * 1024**3)
    if not isinstance(command, list) or not command or any(not isinstance(x, str) or not x or "\0" in x for x in command):
        raise InputError("Invalid native command")
    if "-remote" in command:
        raise InputError("Remote BLAST searches are not supported")
    cwd = Path(cwd).absolute()
    if not cwd.is_dir() or any(p.is_symlink() for p in (cwd, *cwd.parents)):
        raise InputError("Native working directory must be a private regular directory")
    stdout, stderr = Path(stdout).absolute(), Path(stderr).absolute()
    if stdout == stderr:
        raise InputError("Native output streams must be distinct")
    for path in (stdout, stderr):
        if path.parent != cwd or path.exists() or path.is_symlink():
            raise InputError("Native outputs must be new files in the private working directory")
    executable = shutil.which(command[0])
    if executable is None:
        raise InputError(f"Required executable is unavailable: {Path(command[0]).name}")
    command = [str(regular(Path(executable).resolve())), *command[1:]]
    deadline = time.monotonic() + timeout
    process = None
    with stdout.open("xb") as output, stderr.open("xb") as error:
        try:
            process = subprocess.Popen(command, cwd=cwd, env=blast_environment(cwd), stdin=subprocess.DEVNULL,
                                       stdout=output, stderr=error, shell=False, start_new_session=os.name == "posix")
            while True:
                if stdout.stat().st_size + stderr.stat().st_size > max_output_bytes:
                    raise LimitError("Native BLAST output exceeds its byte budget")
                if time.monotonic() > deadline:
                    raise LimitError("Native BLAST deadline exceeded")
                code = process.poll()
                if code is not None:
                    if code != 0:
                        raise InputError(f"{Path(command[0]).name} exited unsuccessfully (status {code})")
                    return
                time.sleep(0.02)
        finally:
            if process is not None:
                # Kill descendants even if the launcher exited first.
                _stop(process)


def executable_identity(command: str) -> dict:
    location = shutil.which(command)
    if location is None:
        raise InputError(f"Required BLAST executable is unavailable: {Path(command).name}")
    path = regular(Path(location).resolve())
    if path.stat().st_size > 512 * 1024**2:
        raise LimitError("Executable exceeds identity byte budget")
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024**2), b""):
            sha.update(chunk)
    with tempfile.TemporaryDirectory(prefix="crisprspacer-version-") as tmp:
        work = Path(tmp)
        run_blast([str(path), "-version"], work, work / "version.txt", work / "stderr.txt", 10, 32_768)
        try:
            version = (work / "version.txt").read_text(encoding="ascii").strip()
        except UnicodeError as exc:
            raise InputError("Invalid BLAST version response") from exc
        if not version or len(version) > 4096:
            raise InputError("Invalid BLAST version response")
    return {"path": str(path), "sha256": sha.hexdigest(), "version": version}
