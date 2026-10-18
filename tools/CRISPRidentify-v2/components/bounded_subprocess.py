"""Run an external command without unbounded pipe buffering.

The runner is intentionally small and POSIX-oriented: CRISPRidentify's
locked deployment target is Linux and process-group termination is required
to prevent a timed-out tool from leaving descendants behind.  Output pipes
are drained concurrently with ``selectors`` while only a bounded stdout head
and stderr tail are retained in memory.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import os
import selectors
import signal
import subprocess
import time
from typing import Mapping, Optional, Sequence, Union


Arg = Union[str, os.PathLike]
_READ_SIZE = 64 * 1024
_POLL_INTERVAL_SECONDS = 0.1
_KILL_DRAIN_SECONDS = 5.0


@dataclass(frozen=True)
class BoundedCompletedProcess:
    """Decoded result of a bounded subprocess invocation."""

    args: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    stdout_truncated: bool = False
    stderr_truncated: bool = False


class BoundedSubprocessError(RuntimeError):
    """Base class for failures raised after a child has been reaped."""

    def __init__(self, message: str, result: BoundedCompletedProcess):
        super().__init__(message)
        self.result = result
        self.args_list = result.args
        self.returncode = result.returncode
        self.stdout = result.stdout
        self.stderr = result.stderr
        self.stdout_truncated = result.stdout_truncated
        self.stderr_truncated = result.stderr_truncated


class BoundedSubprocessTimeoutError(BoundedSubprocessError, TimeoutError):
    """The command exceeded its deadline and its process group was killed."""

    def __init__(
            self,
            timeout: float,
            result: BoundedCompletedProcess,
    ):
        self.timeout = timeout
        super().__init__(
            "command exceeded its {:.3f} second timeout".format(timeout),
            result,
        )


class BoundedSubprocessTruncationError(BoundedSubprocessError):
    """Stdout exceeded its cap and cannot be treated as a complete payload."""

    def __init__(
            self,
            stdout_limit: int,
            result: BoundedCompletedProcess,
    ):
        self.stdout_limit = stdout_limit
        super().__init__(
            "command stdout exceeded its {} byte limit".format(stdout_limit),
            result,
        )


class _HeadBuffer:
    def __init__(self, limit: int):
        self.limit = limit
        self.value = bytearray()
        self.truncated = False

    def append(self, chunk: bytes) -> None:
        remaining = self.limit - len(self.value)
        if remaining > 0:
            self.value.extend(chunk[:remaining])
        if len(chunk) > remaining:
            self.truncated = True


class _TailBuffer:
    def __init__(self, limit: int):
        self.limit = limit
        self.value = bytearray()
        self.truncated = False

    def append(self, chunk: bytes) -> None:
        if len(chunk) >= self.limit:
            if self.value or len(chunk) > self.limit:
                self.truncated = True
            self.value[:] = chunk[-self.limit:]
            return
        overflow = len(self.value) + len(chunk) - self.limit
        if overflow > 0:
            del self.value[:overflow]
            self.truncated = True
        self.value.extend(chunk)


def _positive_timeout(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("timeout must be a positive finite number")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("timeout must be a positive finite number")
    return value


def _positive_limit(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError("{} must be a positive integer".format(name))
    return value


def _validated_argv(argv: Sequence[Arg]) -> tuple[str, ...]:
    if isinstance(argv, (str, bytes, os.PathLike)):
        raise TypeError("argv must be a non-empty sequence, not a shell string")
    try:
        values = tuple(os.fspath(value) for value in argv)
    except TypeError as error:
        raise TypeError("argv entries must be strings or path-like values") from error
    if not values:
        raise ValueError("argv must not be empty")
    if any(not isinstance(value, str) for value in values):
        raise TypeError("argv entries must resolve to text strings")
    return values


def _drain_ready(selector: selectors.BaseSelector, wait: float) -> None:
    for key, _events in selector.select(wait):
        stream = key.fileobj
        try:
            chunk = os.read(stream.fileno(), _READ_SIZE)
        except BlockingIOError:
            continue
        if chunk:
            key.data.append(chunk)
            continue
        selector.unregister(stream)
        stream.close()


def _kill_process_group(process: subprocess.Popen) -> None:
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
            return
        except ProcessLookupError:
            return
    try:
        process.kill()
    except ProcessLookupError:
        pass


def _reap_after_kill(
        process: subprocess.Popen,
        selector: selectors.BaseSelector,
) -> None:
    _kill_process_group(process)
    try:
        process.wait(timeout=_KILL_DRAIN_SECONDS)
    except subprocess.TimeoutExpired:
        process.kill()
        process.wait()

    drain_deadline = time.monotonic() + _KILL_DRAIN_SECONDS
    while selector.get_map() and time.monotonic() < drain_deadline:
        _drain_ready(
            selector,
            min(
                _POLL_INTERVAL_SECONDS,
                max(0.0, drain_deadline - time.monotonic()),
            ),
        )


def _close_streams(
        process: subprocess.Popen,
        selector: selectors.BaseSelector,
) -> None:
    for stream in (process.stdout, process.stderr):
        if stream is None or stream.closed:
            continue
        try:
            selector.unregister(stream)
        except KeyError:
            pass
        stream.close()
    selector.close()


def run_bounded(
        argv: Sequence[Arg],
        *,
        timeout: float,
        stdout_limit: int = 4 * 1024 * 1024,
        stderr_limit: int = 64 * 1024,
        cwd: Optional[Arg] = None,
        env: Optional[Mapping[str, str]] = None,
        encoding: str = "utf-8",
        errors: str = "replace",
) -> BoundedCompletedProcess:
    """Run ``argv`` with bounded output retention and a hard deadline.

    ``stdout`` retains its first ``stdout_limit`` bytes because callers may
    parse it as a complete protocol payload.  Exceeding that limit raises
    :class:`BoundedSubprocessTruncationError`.  ``stderr`` instead retains its
    most recent ``stderr_limit`` bytes so a failed command keeps its most
    useful diagnostic tail; truncation is exposed on the returned result.

    The command is never passed through a shell.  A new session is created so
    a timeout can kill the complete process group, including descendants.
    """

    command = _validated_argv(argv)
    timeout = _positive_timeout(timeout)
    stdout_limit = _positive_limit(stdout_limit, "stdout_limit")
    stderr_limit = _positive_limit(stderr_limit, "stderr_limit")
    stdout_buffer = _HeadBuffer(stdout_limit)
    stderr_buffer = _TailBuffer(stderr_limit)

    process = subprocess.Popen(
        command,
        cwd=os.fspath(cwd) if cwd is not None else None,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    selector = selectors.DefaultSelector()
    assert process.stdout is not None
    assert process.stderr is not None
    for stream, buffer in (
            (process.stdout, stdout_buffer),
            (process.stderr, stderr_buffer)):
        os.set_blocking(stream.fileno(), False)
        selector.register(stream, selectors.EVENT_READ, buffer)

    deadline = time.monotonic() + timeout
    timed_out = False
    try:
        while process.poll() is None or selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                timed_out = True
                break
            if selector.get_map():
                _drain_ready(
                    selector, min(_POLL_INTERVAL_SECONDS, remaining)
                )
            else:
                try:
                    process.wait(
                        timeout=min(_POLL_INTERVAL_SECONDS, remaining)
                    )
                except subprocess.TimeoutExpired:
                    pass
        if timed_out:
            _reap_after_kill(process, selector)
        else:
            process.wait()
    except BaseException:
        _reap_after_kill(process, selector)
        raise
    finally:
        _close_streams(process, selector)

    result = BoundedCompletedProcess(
        args=command,
        returncode=process.returncode,
        stdout=bytes(stdout_buffer.value).decode(encoding, errors),
        stderr=bytes(stderr_buffer.value).decode(encoding, errors),
        stdout_truncated=stdout_buffer.truncated,
        stderr_truncated=stderr_buffer.truncated,
    )
    if timed_out:
        raise BoundedSubprocessTimeoutError(timeout, result)
    if result.stdout_truncated:
        raise BoundedSubprocessTruncationError(stdout_limit, result)
    return result
