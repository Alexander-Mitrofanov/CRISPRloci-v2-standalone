"""Durable polling worker for shell-free Integration executions."""

from __future__ import annotations

import ctypes
import json
import logging
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import threading
import time
from typing import Any, Callable, Mapping
import uuid

from .artifacts import collect_artifacts
from .capacity import available_capacity, has_capacity
from .config import Settings
from .db import ClaimedJob, LeaseLostError, Store
from .execution.files import _capture_stream
from .execution.contracts import MODE_TO_TOOL
from .publication.annotations import write_annotation_exports
from .publication.repeats import write_repeat_exports
from .publication.protospacer import write_protospacer_exports
from .publication.viral import write_viral_exports
from .publication.repeat_type import write_repeat_type_exports
from .publication.association import write_association_exports
from .publication.array_compare import write_array_compare_exports
from .publication.locus_compare import write_locus_compare_exports
from .publication.prophage_context import write_prophage_context_exports
from .publication.insertion_context import write_insertion_context_exports
from .publication.defense_context import write_defense_context_exports
from .publication.rna_family import write_rna_family_exports
from .publication.repeat_map import write_repeat_map_exports
from .execution.input_integrity import _report_source_bindings, _verify_job_inputs
from .execution.provenance import _provenance
from .execution.storage import (
    StorageCleanupError,
    _accounted_storage_bytes,
    _job_tree_within_budget,
    _quarantine_and_purge_attempt,
)
from .execution.validation import (
    _adapter_group_names,
    _supported_schema_major,
    _validated_pipeline_status,
)
from .scientific_options import TOOL_OPTIONS_ENV, validate_tool_options
from .service import (
    CATEGORY_POLICY_CATEGORIES,
    JobService,
    require_unprivileged_process,
)
from .summary import build_summary, sanitize_public_summary
from .publication.sequence_ports import write_sequence_exports
from .publication.module_result import write_module_result


LOGGER = logging.getLogger("crispr_web.worker")
_CANCELLATION_CHECK_SECONDS = 0.2

__all__ = [
    "ProcessContainmentError",
    "StorageCleanupError",
    "Worker",
    "_accounted_storage_bytes",
    "_adapter_group_names",
    "_job_tree_within_budget",
    "_supported_schema_major",
    "main",
]


def _resolve_program(command: tuple[str, ...]) -> None:
    program = command[0]
    if "/" in program or "\\" in program:
        path = Path(program)
        if not path.is_absolute():
            raise ValueError("integration executable path must be absolute")
        if not path.is_file() or not os.access(path, os.X_OK):
            raise FileNotFoundError(f"integration executable is unavailable: {path}")
    elif shutil.which(program) is None:
        raise FileNotFoundError(f"integration executable is not on PATH: {program}")


_PR_SET_CHILD_SUBREAPER = 36
_PR_GET_CHILD_SUBREAPER = 37


def _enable_child_subreaper() -> None:
    """Make escaped grandchildren remain attributable to this one-job worker."""

    if os.name != "posix" or not Path("/proc/self").is_dir():
        raise RuntimeError("production worker containment requires Linux /proc")
    libc = ctypes.CDLL(None, use_errno=True)
    state = ctypes.c_int()
    if libc.prctl(_PR_GET_CHILD_SUBREAPER, ctypes.byref(state), 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    if state.value == 0 and libc.prctl(_PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))


def _proc_parent_and_group(pid: int) -> tuple[int, int] | None:
    try:
        text = Path(f"/proc/{pid}/stat").read_text(encoding="ascii")
        fields = text[text.rfind(")") + 2 :].split()
        return int(fields[1]), int(fields[2])
    except (OSError, ValueError, IndexError):
        return None


def _descendants(root_pid: int) -> dict[int, int]:
    children: dict[int, list[tuple[int, int]]] = {}
    try:
        entries = tuple(Path("/proc").iterdir())
    except OSError:
        return {}
    for entry in entries:
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        relationship = _proc_parent_and_group(pid)
        if relationship is not None:
            parent, group = relationship
            children.setdefault(parent, []).append((pid, group))
    descendants: dict[int, int] = {}
    pending = [root_pid]
    seen = {root_pid}
    while pending:
        parent = pending.pop()
        for child, group in children.get(parent, []):
            if child in seen:
                continue
            seen.add(child)
            pending.append(child)
            descendants[child] = group
    return descendants


def _signal_groups(groups: set[int], signum: signal.Signals) -> None:
    own_group = os.getpgrp()
    for group in sorted(groups, reverse=True):
        if group <= 1 or group == own_group:
            continue
        try:
            os.killpg(group, signum)
        except ProcessLookupError:
            pass
        except PermissionError:
            LOGGER.warning("could not signal process group %s", group)


def _signal_pids(pids: set[int], signum: signal.Signals) -> None:
    for pid in sorted(pids, reverse=True):
        if pid <= 1 or pid == os.getpid():
            continue
        try:
            os.kill(pid, signum)
        except ProcessLookupError:
            pass
        except PermissionError:
            LOGGER.warning("could not signal process %s", pid)


def _reap_adopted_children(pids: set[int], launcher_pid: int) -> None:
    for pid in tuple(pids):
        if pid == launcher_pid:
            continue
        try:
            observed, _status = os.waitpid(pid, os.WNOHANG)
        except (ChildProcessError, ProcessLookupError):
            continue
        if observed == pid:
            pids.discard(pid)


class ProcessContainmentError(RuntimeError):
    """A job-owned process survived the worker's mandatory kill sequence."""


def _terminate_tree(
    process: subprocess.Popen[bytes],
    known_groups: set[int],
    baseline_pids: set[int],
    grace_seconds: float = 10,
    kill_seconds: float = 2,
) -> None:
    """Terminate all attributed processes or fail the worker closed."""

    known_pids: set[int] = {process.pid}

    def refresh() -> None:
        for pid, group in _descendants(process.pid).items():
            known_pids.add(pid)
            known_groups.add(group)
        for pid, group in _descendants(os.getpid()).items():
            if pid not in baseline_pids:
                known_pids.add(pid)
                known_groups.add(group)
        try:
            known_groups.add(os.getpgid(process.pid))
        except ProcessLookupError:
            pass

    def live_processes() -> set[int]:
        process.poll()
        _reap_adopted_children(known_pids, process.pid)
        refresh()
        return {pid for pid in known_pids if Path(f"/proc/{pid}").exists()}

    refresh()
    deadline = time.monotonic() + max(0.0, grace_seconds)
    while True:
        _signal_groups(known_groups, signal.SIGTERM)
        _signal_pids(known_pids, signal.SIGTERM)
        live_pids = live_processes()
        if process.poll() is not None and not live_pids.difference({process.pid}):
            break
        if time.monotonic() >= deadline:
            break
        time.sleep(0.05)

    kill_deadline = time.monotonic() + max(0.0, kill_seconds)
    while True:
        refresh()
        _signal_groups(known_groups, signal.SIGKILL)
        _signal_pids(known_pids, signal.SIGKILL)
        live_pids = live_processes()
        if process.poll() is not None and not live_pids.difference({process.pid}):
            break
        if time.monotonic() >= kill_deadline:
            break
        time.sleep(0.05)
    if process.poll() is None:
        try:
            process.wait(timeout=max(0.05, min(1.0, kill_seconds)))
        except subprocess.TimeoutExpired:
            pass

    live_pids = live_processes()
    launcher_running = process.poll() is None
    residual_pids = (
        live_pids if launcher_running else live_pids.difference({process.pid})
    )
    if launcher_running or residual_pids:
        raise ProcessContainmentError(
            "job process tree survived SIGKILL; refusing further work"
        )


class _WorkerStopRequested(RuntimeError):
    pass


class _CancellationRequested(RuntimeError):
    pass


class _JobDeadlineExceeded(RuntimeError):
    pass


class _JobStorageExceeded(RuntimeError):
    pass



class Worker:
    def __init__(
        self,
        settings: Settings,
        *,
        store: Store | None = None,
        worker_id: str | None = None,
    ):
        self.settings = settings
        self.store = store or Store(settings)
        self.service = JobService(settings, self.store)
        self.worker_id = worker_id or (
            f"{socket.gethostname()}:{os.getpid()}:{uuid.uuid4().hex[:8]}"
        )
        self.stop_event = threading.Event()
        self._containment_failed = False
        self._storage_cleanup_failed = False

    def validate_runtime(self) -> None:
        _enable_child_subreaper()
        _resolve_program(self.settings.integration_command)
        self.settings.verify_scientific_release()
        runner = self.settings.runner_config
        if not runner.is_absolute() or not runner.is_file() or runner.is_symlink():
            raise FileNotFoundError(
                "operator runner config must be an existing non-symlink absolute file"
            )

    def command_for(self, claimed: ClaimedJob, output_root: Path) -> list[str]:
        input_root = self.service.job_root(claimed.job_id) / "input" / "records"
        categories = CATEGORY_POLICY_CATEGORIES[claimed.category_policy]
        command = [
            *self.settings.integration_command,
            "run",
            str(input_root),
            str(output_root),
            "--runner-config",
            str(self.settings.runner_config),
            "--tools",
            MODE_TO_TOOL[claimed.mode],
            "--categories",
            *categories,
            "--spacer-distance",
            str(claimed.spacer_distance),
            "--cpu",
            str(self.settings.worker_cpu),
            "--identify-folder-jobs",
            str(self.settings.identify_folder_jobs),
            "--stage-timeout",
            str(self.settings.stage_timeout_seconds),
        ]
        if claimed.mode == "array_compare":
            command.extend(("--array-map", str(self.service.job_root(claimed.job_id)/"array-map.json")))
        if claimed.mode == "spacer_association":
            if self.settings.association_reference_path is None or claimed.association_reference_sha256 is None or claimed.association_reference_sha256 != self.settings.association_reference_sha256:
                raise ValueError("Pinned association panel is unavailable or changed")
            command.extend(("--association-reference",str(self.settings.association_reference_path),"--association-reference-sha256",claimed.association_reference_sha256,
                "--association-groups",str(self.service.job_root(claimed.job_id)/"association-groups.json")))
        if claimed.mode == "repeat_map":
            if self.settings.repeat_reference_path is None or claimed.repeat_reference_sha256 is None or claimed.repeat_reference_sha256 != self.settings.repeat_reference_sha256:
                raise ValueError("Pinned repeat reference is unavailable or changed")
            command.extend(("--repeat-reference", str(self.settings.repeat_reference_path),
                            "--repeat-reference-sha256", claimed.repeat_reference_sha256))
        if claimed.mode == "viral_search":
            if self.settings.viral_reference_path is None or claimed.viral_reference_sha256 != self.settings.viral_reference_sha256 or claimed.viral_reference_sha256 is None:
                raise ValueError("Pinned viral reference is unavailable or changed")
            command.extend(("--viral-reference", str(self.settings.viral_reference_path),
                            "--viral-reference-sha256", claimed.viral_reference_sha256,
                            "--viral-max-mismatches", str(claimed.viral_max_mismatches)))
        if claimed.mode == "protospacer":
            if self.settings.spacer_database_path is None or claimed.spacer_database_sha256 != self.settings.spacer_database_sha256 or claimed.spacer_database_sha256 is None:
                raise ValueError("Pinned protospacer reference is unavailable or changed")
            command.extend(("--spacer-database", str(self.settings.spacer_database_path),
                            "--spacer-database-sha256", claimed.spacer_database_sha256))
        if claimed.mode in {"leader", "loci"}:
            command.extend(("--leader-flank-length", str(claimed.leader_flank_length)))
        if claimed.mode in {"repeats", "repeat_context", "rna_family"}:
            command.extend(("--molecule", claimed.molecule))
        if claimed.mode in {"tracrrna", "loci"}:
            command.extend(("--tracr-model-type", claimed.tracr_model_type))
        if not claimed.bias_corrections_effective:
            command.append("--disable-bias-corrections")
        return command

    def run_once(self) -> bool:
        _enable_child_subreaper()
        if self._containment_failed:
            raise ProcessContainmentError(
                "this worker is disabled after a process-containment failure"
            )
        if self._storage_cleanup_failed:
            raise StorageCleanupError(
                "this worker is disabled after an attempt-cleanup failure"
            )
        self.settings.verify_scientific_release()
        storage_available = has_capacity(
            self.settings.data_root,
            required_bytes=self.settings.worker_admission_free_bytes,
            required_inodes=self.settings.worker_admission_free_inodes,
        )
        if not storage_available:
            self.store.heartbeat_worker(self.worker_id, "storage_low")
            return False
        self.store.heartbeat_worker(self.worker_id)
        claimed = self.store.claim_next(self.worker_id)
        if claimed is None:
            return False
        try:
            self._execute(claimed)
        except LeaseLostError:
            LOGGER.warning("lease lost for job %s", claimed.job_id)
        except ProcessContainmentError:
            self._containment_failed = True
            LOGGER.critical("process containment failed for job %s", claimed.job_id)
            self._finish_containment_failure(claimed)
            raise
        except StorageCleanupError:
            self._storage_cleanup_failed = True
            LOGGER.critical("storage cleanup failed for job %s", claimed.job_id)
            self._finish_storage_cleanup_failure(claimed)
            raise
        except Exception:
            LOGGER.exception("worker failure for job %s", claimed.job_id)
            self._finish_unexpected_failure(claimed)
        return True

    def _execution_checkpoint(
        self,
        claimed: ClaimedJob,
        deadline: float,
        job_root: Path,
        log_limit_exceeded: threading.Event,
    ) -> Callable[[], None]:
        last_cancellation_check = 0.0
        last_database_check = 0.0
        last_disk_check = 0.0

        def checkpoint() -> None:
            nonlocal last_cancellation_check, last_database_check, last_disk_check
            now = time.monotonic()
            if now - last_cancellation_check >= _CANCELLATION_CHECK_SECONDS:
                if self.store.cancellation_requested(claimed.job_id, self.worker_id):
                    raise _CancellationRequested
                last_cancellation_check = now
            if log_limit_exceeded.is_set():
                raise _JobStorageExceeded
            if now - last_disk_check >= self.settings.job_disk_check_seconds:
                if not _job_tree_within_budget(
                    job_root,
                    self.settings.max_job_disk_bytes,
                    self.settings.max_job_entries,
                ):
                    raise _JobStorageExceeded
                try:
                    free_bytes, free_inodes = available_capacity(
                        self.settings.data_root
                    )
                except OSError as error:
                    raise _JobStorageExceeded from error
                if (
                    free_bytes < self.settings.min_free_bytes
                    or free_inodes < self.settings.min_free_inodes
                ):
                    raise _JobStorageExceeded
                last_disk_check = now
            if self.stop_event.is_set():
                raise _WorkerStopRequested
            if now >= deadline:
                raise _JobDeadlineExceeded
            interval = max(1.0, self.settings.worker_lease_seconds / 3)
            if now - last_database_check >= interval:
                self.store.heartbeat_job(claimed.job_id, self.worker_id)
                self.store.heartbeat_worker(self.worker_id)
                last_database_check = now

        return checkpoint

    def _execute(self, claimed: ClaimedJob) -> None:
        job_root = self.service.job_root(claimed.job_id)
        stored_job = self.store.get_job(claimed.job_id)
        if stored_job is None:
            raise LeaseLostError(claimed.job_id)
        if not _verify_job_inputs(self.settings, job_root, stored_job):
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="failed",
                summary=None,
                artifacts=[],
                error_code="input_integrity_failed",
                error_message="Stored normalized input failed its integrity check.",
            )
            return

        results_root = job_root / "results"
        results_root.mkdir(mode=0o700, exist_ok=True)
        output_root = results_root / f"attempt-{claimed.attempt:03d}"
        if output_root.exists():
            raise RuntimeError("fresh attempt output path already exists")

        launcher_logs = job_root / "worker-logs"
        launcher_logs.mkdir(mode=0o700, exist_ok=True)
        stdout_path = launcher_logs / f"attempt-{claimed.attempt:03d}.stdout.log"
        stderr_path = launcher_logs / f"attempt-{claimed.attempt:03d}.stderr.log"
        temporary_root = job_root / "tmp" / f"attempt-{claimed.attempt:03d}"
        temporary_root.mkdir(parents=True, mode=0o700, exist_ok=False)
        environment = os.environ.copy()
        environment[TOOL_OPTIONS_ENV] = json.dumps(
            validate_tool_options(claimed.tool_options, claimed.mode),
            separators=(",", ":"), allow_nan=False,
        )
        environment.update(
            {
                "TMPDIR": str(temporary_root),
                "TMP": str(temporary_root),
                "TEMP": str(temporary_root),
                "XDG_CACHE_HOME": str(temporary_root / "cache"),
            }
        )

        command = self.command_for(claimed, output_root)
        now_wall = time.time()
        if claimed.started_at <= 0 or claimed.started_at > now_wall + 5:
            remaining_seconds = 0.0
        else:
            elapsed_seconds = max(0.0, now_wall - claimed.started_at)
            remaining_seconds = max(
                0.0, self.settings.job_timeout_seconds - elapsed_seconds
            )
        deadline = time.monotonic() + remaining_seconds
        log_limit_exceeded = threading.Event()
        checkpoint = self._execution_checkpoint(
            claimed, deadline, job_root, log_limit_exceeded
        )
        known_groups: set[int] = set()
        baseline_pids = set(_descendants(os.getpid()))
        process: subprocess.Popen[bytes] | None = None
        capture_threads: list[threading.Thread] = []
        return_code = -1

        try:
            checkpoint()
            try:
                process = subprocess.Popen(
                    command,
                    cwd=job_root,
                    env=environment,
                    stdin=subprocess.DEVNULL,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    start_new_session=True,
                    close_fds=True,
                )
                if process.stdout is None or process.stderr is None:
                    raise RuntimeError("worker log pipes were not created")
                capture_threads = [
                    threading.Thread(
                        target=_capture_stream,
                        args=(
                            process.stdout,
                            stdout_path,
                            self.settings.max_launcher_log_bytes,
                            log_limit_exceeded,
                        ),
                        daemon=True,
                    ),
                    threading.Thread(
                        target=_capture_stream,
                        args=(
                            process.stderr,
                            stderr_path,
                            self.settings.max_launcher_log_bytes,
                            log_limit_exceeded,
                        ),
                        daemon=True,
                    ),
                ]
                for thread in capture_threads:
                    thread.start()
                known_groups.add(process.pid)
                while True:
                    known_groups.update(_descendants(process.pid).values())
                    observed_return_code = process.poll()
                    if observed_return_code is not None:
                        return_code = observed_return_code
                        break
                    checkpoint()
                    time.sleep(0.2)
            finally:
                containment_error: ProcessContainmentError | None = None
                if process is not None:
                    known_groups.update(_descendants(process.pid).values())
                    try:
                        _terminate_tree(process, known_groups, baseline_pids)
                    except ProcessContainmentError as error:
                        containment_error = error
                for thread in capture_threads:
                    thread.join(timeout=5)
                if containment_error is not None:
                    raise containment_error
                if any(thread.is_alive() for thread in capture_threads):
                    raise RuntimeError("bounded worker log capture did not terminate")

            checkpoint()
            if not _verify_job_inputs(self.settings, job_root, stored_job):
                self.store.finish_job(
                    job_id=claimed.job_id,
                    worker_id=self.worker_id,
                    status="failed",
                    summary=None,
                    artifacts=[],
                    error_code="input_integrity_failed",
                    error_message="A scientific subprocess modified its immutable input.",
                )
                return
            if return_code != 0:
                self.store.finish_job(
                    job_id=claimed.job_id,
                    worker_id=self.worker_id,
                    status="failed",
                    summary=None,
                    artifacts=[],
                    error_code="pipeline_failed",
                    error_message="The scientific workflow exited unsuccessfully.",
                )
                return

            expected_source_ids = tuple(stored_job.get("source_ids", ()))
            source_bindings = _report_source_bindings(
                job_root, expected_source_ids, max_bytes=self.settings.max_request_bytes
            )
            if source_bindings is None or len(expected_source_ids) != int(
                stored_job.get("record_count", -1)
            ):
                validated_status = None
            else:
                validated_status = _validated_pipeline_status(
                    claimed,
                    output_root,
                    source_bindings,
                    checkpoint=checkpoint,
                    input_root=job_root / "input" / "records",
                    spacer_reference=self.settings.spacer_database_path,
                    association_reference=self.settings.association_reference_path,
                    repeat_reference=self.settings.repeat_reference_path,
                    viral_reference=self.settings.viral_reference_path,
                )
            if validated_status is None:
                self.store.finish_job(
                    job_id=claimed.job_id,
                    worker_id=self.worker_id,
                    status="failed",
                    summary=None,
                    artifacts=[],
                    error_code="invalid_pipeline_result",
                    error_message=(
                        "The workflow outputs failed cross-file scientific validation."
                    ),
                )
                return

            summary = build_summary(
                output_root,
                max_arrays=self.settings.max_summary_arrays,
                max_rows=self.settings.max_summary_rows,
                max_files=self.settings.max_bundle_files,
                accepted_categories=CATEGORY_POLICY_CATEGORIES[claimed.category_policy],
                checkpoint=checkpoint,
            )
            checkpoint()
            self._finish(
                claimed,
                output_root,
                command,
                status=validated_status,
                summary=summary,
                checkpoint=checkpoint,
            )
        except _WorkerStopRequested:
            self.store.release_job(claimed.job_id, self.worker_id)
        except _CancellationRequested:
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="cancelled",
                summary=None,
                artifacts=[],
            )
        except _JobDeadlineExceeded:
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="failed",
                summary=None,
                artifacts=[],
                error_code="job_timeout",
                error_message=(
                    "The scientific workflow exceeded its fixed wall-time limit."
                ),
            )
        except _JobStorageExceeded:
            export_root = job_root / "exports" / output_root.name
            _quarantine_and_purge_attempt(
                job_root,
                claimed.attempt,
                (
                    ("results", output_root),
                    ("temporary", temporary_root),
                    ("exports", export_root),
                    ("stdout.log", stdout_path),
                    ("stderr.log", stderr_path),
                ),
            )
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="failed",
                summary=None,
                artifacts=[],
                error_code="job_storage_limit",
                error_message=(
                    "The scientific workflow exceeded its fixed storage or log limit."
                ),
            )

    def _private_output_roots(self, job_root: Path) -> tuple[Path, ...]:
        command_roots = tuple(
            Path(part).parent
            for part in self.settings.integration_command
            if Path(part).is_absolute()
        )
        return (
            self.settings.data_root,
            job_root,
            self.settings.runner_config,
            self.settings.runner_config.parent,
            *command_roots,
        )

    def _finish(
        self,
        claimed: ClaimedJob,
        output_root: Path,
        command: list[str],
        *,
        status: str,
        error_code: str | None = None,
        error_message: str | None = None,
        summary: dict[str, Any] | None = None,
        checkpoint: Callable[[], None] | None = None,
    ) -> None:
        if summary is None and output_root.is_dir():
            summary = build_summary(
                output_root,
                max_arrays=self.settings.max_summary_arrays,
                max_rows=self.settings.max_summary_rows,
                max_files=self.settings.max_bundle_files,
                accepted_categories=CATEGORY_POLICY_CATEGORIES[claimed.category_policy],
                checkpoint=checkpoint,
            )
        job_root = self.service.job_root(claimed.job_id)
        if summary is not None:
            identities = self.settings.scientific_release_identities
            summary["provenance"] = {
                "parameters": {"tool_options": claimed.tool_options, "leader_flank_length": claimed.leader_flank_length, "tracr_model_type": claimed.tracr_model_type, "category_policy": claimed.category_policy,
                               "spacer_distance": claimed.spacer_distance, "bias_corrections": claimed.bias_corrections_effective},
                "tool_versions": {entry["display_name"]: entry["version"] for entry in identities.values()},
                "repositories": [{key: entry[key] for key in ("display_name", "repository", "commit")}
                                 for entry in identities.values() if entry.get("repository")],
                "release_id": self.settings.release_id,
            }
            summary = sanitize_public_summary(
                summary,
                private_roots=self._private_output_roots(job_root),
            )
        if summary is not None and status == "completed":
            write_annotation_exports(output_root, summary, checkpoint)
            write_repeat_exports(output_root, summary, checkpoint)
            write_protospacer_exports(output_root, summary, checkpoint)
            write_viral_exports(output_root, summary, checkpoint)
            write_repeat_type_exports(output_root, summary, checkpoint)
            write_association_exports(output_root,summary,checkpoint)
            write_array_compare_exports(output_root,summary,checkpoint)
            write_locus_compare_exports(output_root, summary, checkpoint)
            write_prophage_context_exports(output_root, summary, checkpoint)
            write_insertion_context_exports(output_root, summary, checkpoint)
            write_defense_context_exports(output_root, summary, checkpoint)
            write_rna_family_exports(output_root, summary, checkpoint)
            write_repeat_map_exports(output_root, summary, checkpoint)
        if summary is not None and status in {"completed", "completed_no_eligible_groups"}:
            summary.setdefault("warnings", []).extend(write_sequence_exports(
                output_root, job_root, claimed, self.settings, checkpoint,
            ))
            summary["warnings"].extend(write_module_result(
                output_root, claimed.mode, status,
                max_bytes=self.settings.max_artifact_file_bytes,
                private_roots=self._private_output_roots(job_root), checkpoint=checkpoint,
            ))
        artifacts: list[Mapping[str, Any]] = []
        if output_root.is_dir():
            job = self.store.get_job(claimed.job_id)
            if job is None:
                raise LeaseLostError(claimed.job_id)
            effective_summary = summary or {}
            scientific_identities = self.settings.verify_scientific_release()
            provenance = _provenance(
                job=job,
                claimed=claimed,
                output_root=output_root,
                summary=effective_summary,
                command=command,
                backend_release_id=self.settings.release_id,
                backend_wheel_sha256=self.settings.backend_wheel_sha256,
                scientific_release_manifest_sha256=(
                    self.settings.scientific_release_manifest_sha256
                ),
                scientific_artifacts=scientific_identities,
            )
            handoff_path = job_root / "handoff.json"
            if handoff_path.exists():
                from .summaries.bounded_io import _json
                handoff = _json(job_root, handoff_path, checkpoint=checkpoint)
                provenance["input"]["handoff"] = {
                    "schema": handoff["schema"], "kind": handoff["kind"],
                    "submitted_sha256": handoff["submitted_sha256"],
                    "origin_evidence": "supplied_by_client",
                }
            artifacts, artifact_warnings = collect_artifacts(
                job_root,
                output_root,
                provenance=provenance,
                max_files=self.settings.max_bundle_files,
                max_file_bytes=self.settings.max_artifact_file_bytes,
                max_archive_bytes=self.settings.max_archive_bytes,
                private_roots=self._private_output_roots(job_root),
                checkpoint=checkpoint,
            )
            if summary is not None:
                if artifact_warnings:
                    summary.setdefault("warnings", []).extend(artifact_warnings)
                summary = sanitize_public_summary(
                    summary,
                    private_roots=self._private_output_roots(job_root),
                )
        self.store.finish_job(
            job_id=claimed.job_id,
            worker_id=self.worker_id,
            status=status,
            summary=summary,
            artifacts=artifacts,
            error_code=error_code,
            error_message=error_message,
        )

    def _finish_unexpected_failure(self, claimed: ClaimedJob) -> None:
        try:
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="failed",
                summary=None,
                artifacts=[],
                error_code="worker_error",
                error_message="The backend could not complete this job.",
            )
        except LeaseLostError:
            pass

    def _finish_containment_failure(self, claimed: ClaimedJob) -> None:
        try:
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="failed",
                summary=None,
                artifacts=[],
                error_code="process_containment_failed",
                error_message=(
                    "The scientific process tree could not be safely contained."
                ),
            )
        except Exception:
            LOGGER.exception(
                "could not persist containment failure for job %s", claimed.job_id
            )

    def _finish_storage_cleanup_failure(self, claimed: ClaimedJob) -> None:
        try:
            self.store.finish_job(
                job_id=claimed.job_id,
                worker_id=self.worker_id,
                status="failed",
                summary=None,
                artifacts=[],
                error_code="storage_cleanup_failed",
                error_message=(
                    "The failed attempt's storage could not be safely removed."
                ),
            )
        except Exception:
            LOGGER.exception(
                "could not persist storage cleanup failure for job %s",
                claimed.job_id,
            )

    def run_forever(self) -> None:
        self.validate_runtime()
        self.store.recover_stale()
        last_cleanup = 0.0
        while not self.stop_event.is_set():
            try:
                worked = self.run_once()
                now = time.monotonic()
                if now - last_cleanup >= 300:
                    self.service.cleanup_expired()
                    last_cleanup = now
                if not worked:
                    self.stop_event.wait(self.settings.worker_poll_seconds)
            except ProcessContainmentError:
                LOGGER.critical(
                    "worker process containment failed; terminating for systemd cleanup",
                    exc_info=True,
                )
                raise
            except StorageCleanupError:
                LOGGER.critical(
                    "worker storage cleanup failed; terminating without more claims",
                    exc_info=True,
                )
                raise
            except Exception:
                LOGGER.exception("worker loop iteration failed")
                self.stop_event.wait(self.settings.worker_poll_seconds)
        self.store.heartbeat_worker(self.worker_id, "stopped")


def main() -> None:
    logging.basicConfig(
        level=os.getenv("CRISPR_WEB_LOG_LEVEL", "INFO"),
        format="%(asctime)s %(levelname)s %(name)s %(message)s",
    )
    require_unprivileged_process()
    os.umask(0o077)
    worker = Worker(Settings.from_env())

    def stop(_signum: int, _frame: object) -> None:
        worker.stop_event.set()

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    worker.run_forever()


if __name__ == "__main__":
    main()
