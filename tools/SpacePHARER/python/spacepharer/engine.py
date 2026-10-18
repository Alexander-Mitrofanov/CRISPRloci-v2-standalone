"""Bounded native execution; all descendants are terminated on interruption."""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import sys
import threading
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import BinaryIO

from . import __version__


class EngineError(RuntimeError):
    """Native engine failed or exceeded its deadline."""


def executable() -> str:
    configured = os.environ.get("SPACEPHARER_NATIVE")
    beside = Path(sys.executable).parent / "spacepharer-native"
    value = configured or (str(beside) if beside.is_file() else shutil.which("spacepharer-native"))
    if not value or not Path(value).is_file():
        raise EngineError("Install the pinned spacepharer-native binary")
    value = str(Path(value).resolve())
    try:
        with _process([value, "version"], subprocess.PIPE, subprocess.PIPE) as process:
            stdout, _ = process.communicate(timeout=10)
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode, process.args)
    except (OSError, subprocess.SubprocessError) as error:
        raise EngineError("Cannot verify native engine identity") from error
    if stdout.strip() != __version__.encode("ascii"):
        raise EngineError(f"Native engine version must be {__version__}")
    return value


@contextmanager
def _process(
    args: list[str], stdout: int | BinaryIO, stderr: int = subprocess.STDOUT
) -> Iterator[subprocess.Popen[bytes]]:
    """Contain searches and identity probes under the same process-group lifetime."""
    process = subprocess.Popen(
        args,
        stdout=stdout,
        stderr=stderr,
        start_new_session=True,
        env={**os.environ, "OMP_NUM_THREADS": "1"},
    )
    previous = {}

    def stop(signum, frame):
        raise EngineError("Native execution interrupted")

    try:
        # Python signal handlers belong to the main thread. Every caller still
        # owns the process group, deadline and unconditional descendant cleanup.
        if threading.current_thread() is threading.main_thread():
            for signum in (signal.SIGTERM, signal.SIGINT):
                previous[signum] = signal.signal(signum, stop)
        yield process
    finally:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
        for stream in (process.stdout, process.stderr):
            if stream is not None:
                stream.close()
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def run(args: list[str], log: Path, timeout: int = 900) -> None:
    with log.open("ab") as output, _process([executable(), *args], output) as process:
        try:
            if process.wait(timeout=timeout):
                raise EngineError("Native engine exited unsuccessfully; inspect the local log")
        except subprocess.TimeoutExpired as error:
            raise EngineError("Native deadline exceeded") from error
