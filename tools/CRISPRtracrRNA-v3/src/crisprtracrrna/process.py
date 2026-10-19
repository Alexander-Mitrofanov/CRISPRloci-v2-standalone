"""Logged, shell-free commands with process-group cancellation."""

from __future__ import annotations

import contextlib
import hashlib
import json
import logging
import os
import shutil
import signal
import subprocess
import threading
import time
from pathlib import Path
from typing import NoReturn

from .errors import ToolError

LOG = logging.getLogger(__name__)


def _stderr_tail(path: Path, limit: int = 3000) -> str:
    # UTF-8 uses at most four bytes per character; retain the full log on disk.
    with path.open("rb") as handle:
        size = handle.seek(0, os.SEEK_END)
        handle.seek(max(0, size - 4 * limit))
        return handle.read(4 * limit).decode("utf-8", errors="replace")[-limit:]


def _stop_native(signum: int, frame: object) -> NoReturn:
    # Raising lets the command's finally block reap its isolated group.
    raise SystemExit(128 + signum)


class Runner:
    def __init__(
        self, directory: Path, timeout: float = 3600, environment: dict[str, str] | None = None
    ) -> None:
        self.directory = directory.resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.environment = {**os.environ, **(environment or {}), "PYTHONDONTWRITEBYTECODE": "1"}
        self.records: list[dict[str, object]] = []
        self._digests: dict[str, str] = {}

    def run(
        self,
        command: list[str] | tuple[str, ...],
        *,
        cwd: Path | None = None,
        stdout: Path | None = None,
        environment: dict[str, str] | None = None,
    ) -> Path:
        if (
            not isinstance(command, (list, tuple))
            or not command
            or any(not isinstance(part, str) for part in command)
            or not command[0]
        ):
            raise ToolError("Command must be an array of strings with a nonempty executable")
        ordinal = len(self.records) + 1
        output = stdout or self.directory / f"{ordinal:04d}.stdout.log"
        error = self.directory / f"{ordinal:04d}.stderr.log"
        env = {**self.environment, **(environment or {})}
        child_cwd = (cwd or self.directory).absolute()
        # Popen interprets explicit relative paths and relative PATH entries in
        # the child's cwd. Use that same context for provenance and execution.
        search_path = os.pathsep.join(
            str(child_cwd / entry) for entry in env.get("PATH", os.defpath).split(os.pathsep)
        )
        program = str(child_cwd / command[0]) if os.path.dirname(command[0]) else command[0]
        executable = shutil.which(program, path=search_path)
        if not executable:
            raise ToolError(f"Executable not found: {command[0]}")
        command = [executable, *command[1:]]
        if executable not in self._digests:
            with Path(executable).open("rb") as handle:
                self._digests[executable] = hashlib.file_digest(handle, "sha256").hexdigest()
        record: dict[str, object] = {
            "argv": command,
            "executable_sha256": self._digests[executable],
            "cwd": str(child_cwd),
            "stdout": str(output),
            "stderr": str(error),
            "timeout_seconds": self.timeout,
        }
        self.records.append(record)
        started = time.monotonic()
        completed = False
        try:
            with output.open("wb") as out, error.open("wb") as err:
                proc = subprocess.Popen(
                    command,
                    cwd=child_cwd,
                    env=env,
                    stdout=out,
                    stderr=err,
                    start_new_session=True,
                )
                previous = None
                if threading.current_thread() is threading.main_thread():
                    previous = signal.signal(signal.SIGTERM, _stop_native)
                try:
                    code = proc.wait(timeout=self.timeout)
                except BaseException:
                    # Reap descendants as well as the launcher; never leave inference running.
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(proc.pid, signal.SIGTERM)
                    with contextlib.suppress(subprocess.TimeoutExpired):
                        proc.wait(timeout=2)
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                    raise
                finally:
                    # The launcher can exit while descendants still hold output
                    # files open. Stop them before validation and publication.
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                    if previous is not None:
                        signal.signal(signal.SIGTERM, previous)
                record["returncode"] = code
                if code:
                    tail = _stderr_tail(error)
                    raise ToolError(f"{command[0]} exited {code}. See {error}\n{tail}")
            completed = True
        except subprocess.TimeoutExpired as exc:
            record["timed_out"] = True
            raise ToolError(f"{command[0]} timed out after {self.timeout}s; see {error}") from exc
        except OSError as exc:
            raise ToolError(f"Cannot execute {command[0]}: {exc}") from exc
        finally:
            record["seconds"] = round(time.monotonic() - started, 6)
            journal = self.directory / "commands.json"
            try:
                journal.write_text(json.dumps(self.records, indent=2) + "\n", encoding="utf-8")
            except OSError as exc:
                if completed:
                    raise ToolError(f"Cannot publish command records {journal}: {exc}") from exc
                # A secondary storage failure must preserve failure/timeout/cancellation.
                LOG.exception("Could not publish command records %s", journal)
        return output
