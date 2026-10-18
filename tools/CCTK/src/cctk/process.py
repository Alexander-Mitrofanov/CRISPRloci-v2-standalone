"""Shell-free, bounded optional BLAST/MinCED execution and descendant containment."""

from __future__ import annotations
import os
from pathlib import Path
import sys
import shlex
import signal
import subprocess
import tempfile
import time
from collections.abc import Sequence
from .errors import CCTKError
from .io import MAX_BYTES


def _kill_group(child: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    child.wait()


def run(
    argv: Sequence[str | int | float | Path], *, input_text: str | None = None, timeout: float = 300
) -> subprocess.CompletedProcess[str]:
    if isinstance(argv, str) or not argv or not 0 < timeout <= 3600:
        raise CCTKError("external commands require argument vectors and a bounded deadline")
    command = [str(item) for item in argv]
    payload = b"" if input_text is None else input_text.encode()
    if len(payload) > MAX_BYTES:
        raise CCTKError("external input exceeds limit")
    with (
        tempfile.TemporaryFile() as inp,
        tempfile.TemporaryFile() as out,
        tempfile.TemporaryFile() as err,
    ):
        inp.write(payload)
        inp.seek(0)
        try:
            child = subprocess.Popen(
                [sys.executable, "-I", "-m", "cctk._exec", *command],
                stdin=inp,
                stdout=out,
                stderr=err,
                start_new_session=True,
            )
        except OSError as error:
            raise CCTKError(f"external tool unavailable: {Path(command[0]).name}") from error
        deadline = time.monotonic() + timeout
        try:
            while child.poll() is None:
                if time.monotonic() >= deadline:
                    raise CCTKError("external tool deadline exceeded")
                if any(os.fstat(stream.fileno()).st_size >= MAX_BYTES for stream in (out, err)):
                    raise CCTKError("external output exceeds limit")
                time.sleep(0.02)
        finally:
            # Descendants must not outlive successful leaders, errors or interruption.
            _kill_group(child)
        if any(os.fstat(stream.fileno()).st_size >= MAX_BYTES for stream in (out, err)):
            raise CCTKError("external output exceeds limit")
        out.seek(0)
        err.seek(0)
        try:
            stdout, stderr = out.read(MAX_BYTES + 1).decode(), err.read(MAX_BYTES + 1).decode()
        except UnicodeDecodeError as error:
            raise CCTKError("external tool emitted invalid UTF-8") from error
        if child.returncode:
            raise CCTKError(
                f"external tool failed: {Path(command[0]).name} (exit {child.returncode})"
            )
        return subprocess.CompletedProcess(command, child.returncode, stdout, stderr)


def blast_options(text: str | None) -> list[str]:
    """Extra BLAST arguments are parsed as arguments, never executable shell text."""
    args = shlex.split(text or "")
    if any(
        x in {"-out", "-outfmt", "-query", "-subject", "-db", "-export_search_strategy"}
        for x in args
    ):
        raise CCTKError("extra BLAST options may not replace input/output contracts")
    return args
