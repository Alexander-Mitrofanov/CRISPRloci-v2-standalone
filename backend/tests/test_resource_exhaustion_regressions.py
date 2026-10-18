from __future__ import annotations

import json
import os
from pathlib import Path
import stat
import threading
from types import SimpleNamespace

from fastapi.testclient import TestClient
import pytest

import crispr_web.api as api_module
import crispr_web.capacity as capacity_module
import crispr_web.summary as summary_module
import crispr_web.worker as worker_module
from crispr_web.artifacts import _bundle_outputs, _safe_file
from crispr_web.db import SubmissionRateLimitError
from crispr_web.models import JobSubmission
from crispr_web.service import (
    JobService,
    StorageCapacityError,
    _remove_private_job_tree,
)
from crispr_web.worker import (
    Worker,
    _accounted_storage_bytes,
    _job_tree_within_budget,
)


def _submission(identifier: str) -> JobSubmission:
    return JobSubmission(sequence=f">{identifier}\nACGT\n", mode="detection")


def _unlimited_capacity(_path) -> tuple[int, int]:
    return 10**18, 10**18


@pytest.mark.parametrize("target", ["normalized", "record"])
def test_input_verification_rejects_oversized_files_before_opening(
    settings_factory, monkeypatch, target
) -> None:
    from crispr_web.execution.input_integrity import _verify_job_inputs

    service = JobService(settings_factory())
    _, queued = service.submit(_submission("bounded-input"), "192.0.2.45")
    root = service.job_root(queued["job_id"])
    job = service.store.get_job(queued["job_id"])
    assert _verify_job_inputs(service.settings, root, job)
    path = (root / "input/normalized.fasta" if target == "normalized"
            else next((root / "input/records").iterdir()))
    path.chmod(0o600)
    with path.open("wb") as stream:
        stream.truncate(service.settings.max_request_bytes + 1)
    original_open = Path.open

    def guarded_open(self, *args, **kwargs):
        assert self != path, "oversized input must be rejected before reading"
        return original_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded_open)
    assert not _verify_job_inputs(service.settings, root, job)


def test_source_bindings_reject_oversized_files_before_opening(settings_factory, monkeypatch):
    from crispr_web.execution.input_integrity import _report_source_bindings

    service = JobService(settings_factory())
    _, queued = service.submit(_submission("bounded-source"), "192.0.2.46")
    root = service.job_root(queued["job_id"])
    assert _report_source_bindings(root, ("bounded-source",))
    path = next((root / "input/records").iterdir())
    path.chmod(0o600)
    with path.open("wb") as stream:
        stream.truncate(20_000_001)
    original_open = Path.open

    def guarded_open(self, *args, **kwargs):
        assert self != path, "oversized source must be rejected before reading"
        return original_open(self, *args, **kwargs)

    monkeypatch.setattr(Path, "open", guarded_open)
    assert _report_source_bindings(root, ("bounded-source",)) is None


def test_global_submission_budget_survives_rotating_ip_cancel_churn_and_restart(
    settings_factory, monkeypatch
) -> None:
    monkeypatch.setattr(capacity_module, "available_capacity", _unlimited_capacity)
    settings = settings_factory(
        max_active_jobs=1,
        max_active_jobs_per_client=1,
        max_queued_jobs=1,
        rate_limit_submissions=20,
        global_rate_limit_submissions=3,
        global_rate_limit_window_seconds=3_600,
    )
    service = JobService(settings)

    for index in range(3):
        _token, queued = service.submit(
            _submission(f"churn-{index}"), f"198.51.100.{index + 10}"
        )
        root = service.job_root(queued["job_id"])
        assert root.is_dir()
        assert service.cancel(queued["job_id"])["status"] == "cancelled"
        assert not root.exists()

    restarted = JobService(settings)
    with pytest.raises(SubmissionRateLimitError):
        restarted.submit(_submission("rotated-again"), "203.0.113.200")

    assert list(settings.jobs_root.iterdir()) == []
    with restarted.store.connect() as connection:
        statuses = connection.execute(
            "SELECT status FROM jobs ORDER BY created_at"
        ).fetchall()
        global_events = connection.execute(
            "SELECT COUNT(*) FROM rate_events "
            "WHERE client_digest = 'global-submissions-v1'"
        ).fetchone()[0]
    assert [row["status"] for row in statuses] == ["cancelled"] * 3
    assert global_events == 3


def test_queued_cancel_removes_private_tree_but_keeps_authorized_db_tombstone(
    settings_factory, monkeypatch
) -> None:
    monkeypatch.setattr(capacity_module, "available_capacity", _unlimited_capacity)
    service = JobService(settings_factory())
    token, queued = service.submit(_submission("queued-purge"), "192.0.2.40")
    job_id = queued["job_id"]
    root = service.job_root(job_id)
    assert (root / "input" / "normalized.fasta").is_file()

    cancelled = service.cancel(job_id)

    assert cancelled["status"] == "cancelled"
    assert not root.exists()
    assert service.authorize(job_id, token) is True
    stored = service.store.get_job(job_id)
    assert stored is not None
    assert stored["status"] == "cancelled"
    assert bool(stored["cancel_requested"]) is True
    assert stored["summary"] is None
    assert stored["artifacts"] == []


def test_submission_admission_fails_closed_when_only_inode_reserve_is_low(
    settings_factory, monkeypatch
) -> None:
    settings = settings_factory(min_free_bytes=0, min_free_inodes=10)
    service = JobService(settings)
    monkeypatch.setattr(
        capacity_module,
        "available_capacity",
        lambda _path: (
            settings.submission_admission_free_bytes,
            settings.submission_admission_free_inodes - 1,
        ),
    )

    with pytest.raises(StorageCapacityError, match="byte or inode reserve"):
        service.submit(_submission("inode-admission"), "192.0.2.41")

    assert list(settings.jobs_root.iterdir()) == []


def test_readiness_is_low_when_only_inode_reserve_is_low(
    settings_factory, monkeypatch
) -> None:
    settings = settings_factory(
        readiness_cache_seconds=0,
        min_free_bytes=0,
        min_free_inodes=10,
    )
    app = api_module.create_app(settings)
    required_bytes = max(
        settings.submission_admission_free_bytes,
        settings.worker_admission_free_bytes,
    )
    required_inodes = max(
        settings.submission_admission_free_inodes,
        settings.worker_admission_free_inodes,
    )
    monkeypatch.setattr(app.state.store, "database_ok", lambda: True)
    monkeypatch.setattr(app.state.store, "worker_health", lambda: "ok")
    monkeypatch.setattr(
        api_module,
        "available_capacity",
        lambda _path: (required_bytes, required_inodes - 1),
    )

    with TestClient(app) as client:
        response = client.get("/api/v1/health")

    assert response.status_code == 503
    assert response.json()["storage"] == "low"


def test_low_inode_worker_and_store_guards_leave_job_unclaimed(
    settings_factory, monkeypatch
) -> None:
    monkeypatch.setattr(capacity_module, "available_capacity", _unlimited_capacity)
    settings = settings_factory(min_free_bytes=0, min_free_inodes=10)
    service = JobService(settings)
    _token, queued = service.submit(_submission("inode-worker"), "192.0.2.42")
    monkeypatch.setattr(
        capacity_module,
        "available_capacity",
        lambda _path: (
            settings.worker_admission_free_bytes,
            settings.worker_admission_free_inodes - 1,
        ),
    )

    assert service.store.claim_next("direct-inode-guard") is None
    worker = Worker(settings, store=service.store, worker_id="inode-guard-worker")
    assert worker.run_once() is False

    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "queued"
    assert stored["attempt"] == 0
    assert service.store.worker_health() == "storage_low"


def test_top_level_symlink_cleanup_never_follows_or_chmods_target(tmp_path) -> None:
    target = tmp_path / "operator-owned"
    target.mkdir()
    target.chmod(0o755)
    sentinel = target / "keep.txt"
    sentinel.write_text("keep\n", encoding="utf-8")
    link = tmp_path / "job-tree"
    try:
        link.symlink_to(target, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")

    with pytest.raises(OSError, match="not a directory"):
        _remove_private_job_tree(link)
    assert _remove_private_job_tree(link, ignore_errors=True) is False

    assert link.is_symlink()
    assert stat.S_IMODE(target.stat().st_mode) == 0o755
    assert sentinel.read_text(encoding="utf-8") == "keep\n"


@pytest.mark.parametrize(
    ("metadata", "expected"),
    [
        (
            SimpleNamespace(
                st_mode=stat.S_IFREG | 0o600,
                st_size=8_192,
                st_blocks=1,
            ),
            8_192,
        ),
        (
            SimpleNamespace(
                st_mode=stat.S_IFREG | 0o600,
                st_size=1,
                st_blocks=32,
            ),
            16_384,
        ),
        (
            SimpleNamespace(
                st_mode=stat.S_IFDIR | 0o700,
                st_size=99_999,
                st_blocks=8,
            ),
            4_096,
        ),
    ],
)
def test_storage_accounting_uses_larger_apparent_or_allocated_charge(
    metadata, expected
) -> None:
    assert _accounted_storage_bytes(metadata) == expected


def test_job_tree_budget_charges_root_nodes_and_entry_count() -> None:
    root_metadata = SimpleNamespace(
        st_mode=stat.S_IFDIR | 0o700,
        st_size=4_096,
        st_blocks=8,
        st_dev=1,
        st_ino=1,
    )
    child_metadata = (
        SimpleNamespace(st_mode=stat.S_IFDIR | 0o700, st_size=4_096, st_blocks=16),
        SimpleNamespace(st_mode=stat.S_IFREG | 0o600, st_size=10_000, st_blocks=1),
    )
    children = tuple(
        SimpleNamespace(lstat=lambda metadata=metadata: metadata)
        for metadata in child_metadata
    )
    root = SimpleNamespace(
        lstat=lambda: root_metadata,
        rglob=lambda pattern: iter(children),
    )
    total = _accounted_storage_bytes(root_metadata) + sum(
        _accounted_storage_bytes(metadata) for metadata in child_metadata
    )

    assert _job_tree_within_budget(root, total, len(children)) is True
    assert _job_tree_within_budget(root, total - 1, len(children)) is False
    assert _job_tree_within_budget(root, total, len(children) - 1) is False


def test_job_tree_budget_tolerates_removed_temporary_files(tmp_path) -> None:
    transient = tmp_path / "temporary"
    transient.write_bytes(b"temporary output")
    retained = tmp_path / "retained"
    retained.write_bytes(b"retained output")
    budget = sum(
        _accounted_storage_bytes(path.lstat()) for path in (tmp_path, retained)
    )

    def concurrent_cleanup(_pattern):
        transient.unlink(missing_ok=True)
        yield transient
        yield retained

    root = SimpleNamespace(lstat=tmp_path.lstat, rglob=concurrent_cleanup)
    assert _job_tree_within_budget(root, budget, 2) is True
    assert _job_tree_within_budget(root, budget - 1, 2) is False
    # Vanished entries cannot make traversal unbounded.
    assert _job_tree_within_budget(root, budget, 1) is False


def test_job_tree_budget_rejects_other_inspection_errors(tmp_path) -> None:
    def denied():
        raise PermissionError("cannot inspect scientific output")

    root = SimpleNamespace(
        lstat=tmp_path.lstat,
        rglob=lambda _pattern: iter((SimpleNamespace(lstat=denied),)),
    )
    assert _job_tree_within_budget(root, 1_000_000, 100) is False


def test_job_tree_budget_rejects_root_replacement_during_scan(tmp_path) -> None:
    root = tmp_path / "job"
    root.mkdir()

    def replaced(_pattern):
        root.rename(tmp_path / "moved")
        root.mkdir()
        yield from ()

    view = SimpleNamespace(lstat=root.lstat, rglob=replaced)
    assert _job_tree_within_budget(view, 1_000_000, 100) is False


def test_hard_linked_output_is_rejected_from_public_artifacts(tmp_path) -> None:
    output_root = tmp_path / "output"
    candidate_root = output_root / "spacerplacer"
    candidate_root.mkdir(parents=True)
    outside = tmp_path / "operator-owned.txt"
    outside.write_text("do not publish\n", encoding="utf-8")
    candidate = candidate_root / "summary.txt"
    try:
        os.link(outside, candidate)
    except OSError:
        pytest.skip("hard links unavailable on this filesystem")
    assert candidate.stat().st_nlink == 2

    assert _safe_file(output_root, candidate) is None
    files, _warnings = _bundle_outputs(
        output_root,
        max_file_bytes=1_000_000,
        max_files=100,
        max_total_bytes=1_000_000,
        markers=(),
        checkpoint=None,
    )
    assert files == []


def test_cancellation_database_checks_are_cadence_limited(
    settings_factory, monkeypatch
) -> None:
    monkeypatch.setattr(capacity_module, "available_capacity", _unlimited_capacity)
    settings = settings_factory(job_disk_check_seconds=60)
    service = JobService(settings)
    _token, queued = service.submit(_submission("cancel-cadence"), "192.0.2.43")
    worker_id = "cancel-cadence-worker"
    claimed = service.store.claim_next(worker_id)
    assert claimed is not None
    worker = Worker(settings, store=service.store, worker_id=worker_id)
    calls = 0

    def cancellation_requested(job_id, owner) -> bool:
        nonlocal calls
        calls += 1
        assert job_id == queued["job_id"]
        assert owner == worker_id
        return False

    now = [100.0]
    monkeypatch.setattr(service.store, "cancellation_requested", cancellation_requested)
    monkeypatch.setattr(worker_module.time, "monotonic", lambda: now[0])
    monkeypatch.setattr(worker_module, "_job_tree_within_budget", lambda *_args: True)
    monkeypatch.setattr(worker_module, "available_capacity", _unlimited_capacity)
    checkpoint = worker._execution_checkpoint(
        claimed,
        deadline=1_000.0,
        job_root=service.job_root(claimed.job_id),
        log_limit_exceeded=threading.Event(),
    )

    for _ in range(1_000):
        checkpoint()
    assert calls == 1

    now[0] = 100.19
    checkpoint()
    assert calls == 1
    now[0] = 100.21
    checkpoint()
    assert calls == 2


def test_final_stored_summary_is_recapped_after_artifact_warnings(
    settings_factory, monkeypatch
) -> None:
    monkeypatch.setattr(capacity_module, "available_capacity", _unlimited_capacity)
    settings = settings_factory()
    service = JobService(settings)
    _token, queued = service.submit(_submission("warning-cap"), "192.0.2.44")
    worker_id = "warning-cap-worker"
    claimed = service.store.claim_next(worker_id)
    assert claimed is not None
    worker = Worker(settings, store=service.store, worker_id=worker_id)
    output_root = (
        service.job_root(claimed.job_id) / "results" / f"attempt-{claimed.attempt:03d}"
    )
    output_root.mkdir(parents=True)
    artifact_warnings = ["x" * 10_000 for _ in range(600)]
    assert len(json.dumps(artifact_warnings).encode("utf-8")) > (
        summary_module._MAX_PUBLIC_SUMMARY_BYTES
    )
    monkeypatch.setattr(
        worker_module,
        "collect_artifacts",
        lambda *_args, **_kwargs: ([], artifact_warnings),
    )
    summary = {
        "schema_version": "1.0.0",
        "pipeline_status": "completed",
        "warnings": [],
    }

    worker._finish(
        claimed,
        output_root,
        worker.command_for(claimed, output_root),
        status="completed",
        summary=summary,
    )

    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "completed"
    encoded = json.dumps(
        stored["summary"], separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    assert len(encoded) <= summary_module._MAX_PUBLIC_SUMMARY_BYTES
    assert any(
        "exceeded the byte limit" in warning
        for warning in stored["summary"]["warnings"]
    )
