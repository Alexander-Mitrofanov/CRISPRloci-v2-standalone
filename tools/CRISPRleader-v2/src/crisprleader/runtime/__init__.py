"""Bounded POSIX subprocess execution for maintainer-selected scientific tools.

Containment is a process-group boundary, not an operating-system sandbox. A
malicious process can escape its group; callers must never run untrusted code.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
import os
from pathlib import Path
import selectors
import signal
import subprocess
import threading
import time

from ..errors import ExecutionError


@dataclass(frozen=True)
class CommandResult:
    """Successful command output, decoded as UTF-8 with replacement."""

    args: tuple[str, ...]
    stdout: str
    stderr: str
    returncode: int


def _kill_group(process: subprocess.Popen[bytes]) -> None:
    """Kill descendants even when the leader has exited, then reap the leader.

    A session leader's PID remains reserved until reaping. We deliberately do
    not poll/reap the leader while collecting pipes, so that this group ID
    cannot be recycled before cleanup. ESRCH means there are no group members.
    SIGKILL avoids an unbounded grace period and children ignoring SIGTERM.
    """
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_command(
    argv: Sequence[str | os.PathLike[str]],
    *,
    cwd: str | os.PathLike[str],
    timeout_seconds: float = 120.0,
    max_output_bytes: int = 1048576,
    cancel_event: threading.Event | None = None,
    env: Mapping[str, str] | None = None,
) -> CommandResult:
    """Run fixed argv without a shell and with one aggregate output bound.

    Timeout, cancellation, overflow, launch failure and nonzero exit raise
    :class:`ExecutionError`. ``env`` replaces the environment if supplied.
    Both pipes are consumed concurrently; the deadline includes descendants
    retaining an inherited output pipe after the command leader exits. Even on
    successful leader exit, remaining group members are terminated. ``cwd``
    must exist; this function does not create or reuse shared directories.
    """
    if os.name != "posix" or not hasattr(os, "waitid") or not hasattr(os, "WNOWAIT"):
        raise ExecutionError("Process-group containment requires POSIX waitid/WNOWAIT support")
    if isinstance(argv, (str, bytes)) or not argv:
        raise ExecutionError("Command must be a nonempty argument sequence, never a shell string")
    try:
        args = tuple(os.fspath(value) for value in argv)
    except TypeError as exc:
        raise ExecutionError("Every command argument must be text or a path") from exc
    if any(not isinstance(value, str) or "\0" in value for value in args) or not args[0]:
        raise ExecutionError("Command arguments must be NUL-free text; executable cannot be empty")
    if isinstance(timeout_seconds, bool) or not isinstance(timeout_seconds, (int, float)):
        raise ExecutionError("timeout_seconds must be a positive finite number")
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ExecutionError("timeout_seconds must be a positive finite number")
    if type(max_output_bytes) is not int or max_output_bytes < 0:
        raise ExecutionError("max_output_bytes must be a nonnegative integer")
    if cancel_event is not None and cancel_event.is_set():
        raise ExecutionError("Command cancelled before launch")
    if not Path(cwd).is_dir():
        raise ExecutionError("Command working directory does not exist or is not a directory")

    deadline = time.monotonic() + timeout_seconds
    try:
        process = subprocess.Popen(
            args,
            cwd=cwd,
            env=env,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False,
            close_fds=True,
            start_new_session=True,
        )
    except (OSError, ValueError, TypeError) as exc:
        raise ExecutionError(f"Could not start command: {type(exc).__name__}: {exc}") from exc

    output: dict[str, bytearray] = {"stdout": bytearray(), "stderr": bytearray()}
    consumed = 0
    try:
        assert process.stdout is not None and process.stderr is not None
        with selectors.DefaultSelector() as selector:
            for pipe, name in ((process.stdout, "stdout"), (process.stderr, "stderr")):
                os.set_blocking(pipe.fileno(), False)
                selector.register(pipe, selectors.EVENT_READ, name)
            while selector.get_map():
                if cancel_event is not None and cancel_event.is_set():
                    raise ExecutionError("Command cancelled; process group terminated")
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise ExecutionError("Command timed out; process group terminated")
                for key, _ in selector.select(min(remaining, 0.05)):
                    # Read at most one byte beyond the allowed aggregate size;
                    # overflow cannot allocate an attacker-controlled chunk.
                    chunk = os.read(key.fd, min(65536, max_output_bytes - consumed + 1))
                    if not chunk:
                        selector.unregister(key.fileobj)
                    else:
                        consumed += len(chunk)
                        if consumed > max_output_bytes:
                            raise ExecutionError(
                                "Command exceeded combined stdout/stderr byte limit; "
                                "process group terminated"
                            )
                        output[key.data].extend(chunk)
            # Pipes can close while the process remains alive. A timed wait
            # with cancellation checks is necessary after both EOFs too.
            while True:
                if cancel_event is not None and cancel_event.is_set():
                    raise ExecutionError("Command cancelled; process group terminated")
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise ExecutionError("Command timed out; process group terminated")
                # waitid(WNOWAIT) observes exit without reaping/reserving-ID loss.
                status = os.waitid(os.P_PID, process.pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
                if status is not None:
                    break
                time.sleep(min(remaining, 0.01))
        # Kill any descendants before reaping, including descendants that
        # deliberately closed their output pipes and outlived a successful parent.
        _kill_group(process)
        stdout = output["stdout"].decode("utf-8", errors="replace")
        stderr = output["stderr"].decode("utf-8", errors="replace")
        if process.returncode != 0:
            detail = stderr[-4096:]
            raise ExecutionError(f"Command exited with status {process.returncode}: {detail}")
        return CommandResult(args=args, stdout=stdout, stderr=stderr, returncode=process.returncode)
    except ExecutionError:
        if process.returncode is None:
            _kill_group(process)
        raise
    except BaseException:
        # Includes interruption by the calling worker, not only ordinary errors.
        if process.returncode is None:
            _kill_group(process)
        raise
    finally:
        if process.stdout is not None:
            process.stdout.close()
        if process.stderr is not None:
            process.stderr.close()
