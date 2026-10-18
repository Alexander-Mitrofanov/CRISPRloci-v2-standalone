from __future__ import annotations

from types import SimpleNamespace
import os
import sqlite3
import stat
import threading
import time

import pytest

from crispr_web.db import QueueCapacityError, Store, SubmissionRateLimitError
import crispr_web.db as db_module
from crispr_web.models import JobSubmission
from crispr_web.service import JobService
import crispr_web.service as service_module


def _submission(identifier: str = "sample") -> JobSubmission:
    return JobSubmission(sequence=f">{identifier}\nACGT\n", mode="detection")


def test_expired_queued_job_is_never_authorized_cancelled_or_claimed(settings_factory):
    service = JobService(settings_factory())
    token, queued = service.submit(_submission("expired"), "192.0.2.60")
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET expires_at = ? WHERE job_id = ?",
            (time.time() - 1, queued["job_id"]),
        )

    assert service.authorize(queued["job_id"], token) is False
    with pytest.raises(KeyError):
        service.cancel(queued["job_id"])
    assert service.store.claim_next("expiry-worker") is None
    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "cancelled"
    assert stored["error_code"] == "job_expired"


@pytest.mark.parametrize("transition", ["release", "finish"])
def test_cancellation_wins_worker_completion_races(settings_factory, transition):
    service = JobService(settings_factory())
    _, queued = service.submit(_submission(transition), "192.0.2.61")
    worker_id = f"{transition}-worker"
    claim = service.store.claim_next(worker_id)
    assert claim is not None
    service.cancel(claim.job_id)

    if transition == "release":
        service.store.release_job(claim.job_id, worker_id)
    else:
        service.store.finish_job(
            job_id=claim.job_id,
            worker_id=worker_id,
            status="completed",
            summary={"private": "must not publish"},
            artifacts=[
                {
                    "artifact_id": "a" * 24,
                    "relative_path": "results.txt",
                    "name": "results.txt",
                    "size_bytes": 1,
                    "sha256": "b" * 64,
                    "media_type": "text/plain",
                }
            ],
        )

    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "cancelled"
    assert stored["summary"] is None
    assert stored["artifacts"] == []


def test_low_disk_worker_does_not_claim_job(settings_factory, monkeypatch):
    settings = settings_factory(min_free_bytes=1)
    service = JobService(settings)
    _, queued = service.submit(_submission("lowdisk"), "192.0.2.62")
    monkeypatch.setattr(
        "crispr_web.capacity.available_capacity", lambda _path: (0, 10**18)
    )

    assert service.store.claim_next("low-disk-worker") is None
    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "queued"
    assert stored["attempt"] == 0


def test_newer_database_schema_is_rejected_without_downgrade(settings_factory):
    settings = settings_factory()
    store = Store(settings)
    with store.connect() as connection:
        connection.execute(
            "UPDATE metadata SET value = '999' WHERE name = 'schema_version'"
        )

    with pytest.raises(RuntimeError, match="newer than this backend"):
        Store(settings)
    with store.connect() as connection:
        value = connection.execute(
            "SELECT value FROM metadata WHERE name = 'schema_version'"
        ).fetchone()[0]
    assert value == "999"


def test_database_ceiling_below_one_page_fails_initialization(
    settings_factory, monkeypatch
) -> None:
    monkeypatch.setattr(db_module, "_MAX_DATABASE_BYTES", 1)

    with pytest.raises(sqlite3.DatabaseError):
        Store(settings_factory())


def test_persisted_limits_reject_before_fasta_parsing(settings_factory, monkeypatch):
    rate_service = JobService(
        settings_factory(max_queued_jobs=4, rate_limit_submissions=1)
    )
    rate_service.submit(_submission("first"), "192.0.2.63")

    def parser_must_not_run(*_args, **_kwargs):
        raise AssertionError("FASTA parser ran before persisted preflight")

    monkeypatch.setattr("crispr_web.service.normalize_fasta", parser_must_not_run)
    with pytest.raises(SubmissionRateLimitError):
        rate_service.submit(_submission("rate"), "192.0.2.63")

    queue_service = JobService(
        settings_factory(max_queued_jobs=1, rate_limit_submissions=20)
    )
    with pytest.raises(QueueCapacityError):
        queue_service.submit(_submission("queue"), "192.0.2.64")


def test_materialized_inputs_are_read_only_and_private(settings_factory):
    service = JobService(settings_factory())
    _, queued = service.submit(_submission("permissions"), "192.0.2.65")
    input_root = service.job_root(queued["job_id"]) / "input"
    records_root = input_root / "records"
    files = [input_root / "normalized.fasta", *records_root.glob("*.fasta")]

    assert stat.S_IMODE(input_root.stat().st_mode) == 0o500
    assert stat.S_IMODE(records_root.stat().st_mode) == 0o500
    assert files
    assert all(stat.S_IMODE(path.stat().st_mode) == 0o400 for path in files)


def test_old_unregistered_job_tree_is_atomically_quarantined_and_removed(
    settings_factory,
):
    service = JobService(settings_factory())
    orphan = service.settings.jobs_root / ("a" * 32)
    orphan.mkdir(mode=0o700)
    (orphan / ".submission.lock").write_bytes(b"")
    (orphan / "partial-input.fasta").write_text(">partial\nACGT\n", encoding="ascii")
    old = time.time() - 7_200
    os.utime(orphan, (old, old))
    unrelated = service.settings.jobs_root / "operator-notes"
    unrelated.mkdir()
    os.utime(unrelated, (old, old))

    assert service.reconcile_orphaned_job_trees(now=time.time()) == 1
    assert not orphan.exists()
    assert unrelated.is_dir()
    quarantine = service.settings.data_root / ".orphan-quarantine"
    assert list(quarantine.iterdir()) == []


def test_orphan_reconciliation_rejects_symlinked_quarantine(settings_factory, tmp_path):
    service = JobService(settings_factory())
    outside = tmp_path / "outside-quarantine-target"
    outside.mkdir()
    sentinel = outside / "keep.txt"
    sentinel.write_text("keep\n", encoding="utf-8")
    quarantine = service.settings.data_root / ".orphan-quarantine"
    try:
        quarantine.symlink_to(outside, target_is_directory=True)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")

    assert service.reconcile_orphaned_job_trees(now=time.time() + 7_200) == 0
    assert sentinel.read_text(encoding="utf-8") == "keep\n"


def test_inflight_submission_lock_prevents_orphan_cleanup_race(
    settings_factory, monkeypatch
):
    service = JobService(settings_factory())
    entered = threading.Event()
    release = threading.Event()
    original_materialize = service_module.materialize_fasta
    result: list[object] = []

    def blocked_materialize(root, normalized):
        entered.set()
        assert release.wait(timeout=5)
        return original_materialize(root, normalized)

    def submit() -> None:
        try:
            result.append(service.submit(_submission("inflight"), "192.0.2.166"))
        except BaseException as error:
            result.append(error)

    monkeypatch.setattr(service_module, "materialize_fasta", blocked_materialize)
    thread = threading.Thread(target=submit)
    thread.start()
    assert entered.wait(timeout=3)
    roots = list(service.settings.jobs_root.iterdir())
    assert len(roots) == 1
    in_flight_root = roots[0]

    # Even if a cleanup clock makes the directory appear older than the grace
    # period, the held submission lock prevents quarantine.
    assert service.reconcile_orphaned_job_trees(now=time.time() + 7_200) == 0
    assert in_flight_root.is_dir()

    release.set()
    thread.join(timeout=5)
    assert not thread.is_alive()
    assert len(result) == 1 and not isinstance(result[0], BaseException)
    queued = result[0][1]
    assert service.store.get_job(queued["job_id"]) is not None
    assert not (in_flight_root / ".submission.lock").exists()


def test_database_failure_cleans_read_only_uncommitted_job_tree(
    settings_factory, monkeypatch
):
    service = JobService(settings_factory())

    def fail_create_job(**_kwargs):
        raise RuntimeError("simulated durable queue failure")

    monkeypatch.setattr(service.store, "create_job", fail_create_job)
    with pytest.raises(RuntimeError, match="durable queue failure"):
        service.submit(_submission("cleanup"), "192.0.2.66")
    assert list(service.settings.jobs_root.iterdir()) == []


def test_database_readiness_probe_is_transactional_and_full_sync(settings_factory):
    store = Store(settings_factory())
    assert store.database_ok() is True
    with store.connect() as connection:
        synchronous = connection.execute("PRAGMA synchronous").fetchone()[0]
        assert synchronous == 2
        assert connection.in_transaction is False


def test_total_execution_budget_persists_across_retries(settings_factory):
    settings = settings_factory(max_job_attempts=3)
    service = JobService(settings)
    _, queued = service.submit(_submission("budget"), "192.0.2.67")
    base = time.time()
    first = service.store.claim_next("budget-worker-1", now=base)
    assert first is not None
    service.store.release_job(first.job_id, "budget-worker-1")
    persisted = service.store.get_job(queued["job_id"])
    assert persisted is not None
    assert persisted["started_at"] == pytest.approx(first.started_at)

    second = service.store.claim_next("budget-worker-2", now=base + 1)
    assert second is not None
    assert second.started_at == pytest.approx(first.started_at)
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET started_at = ? WHERE job_id = ?",
            (time.time() - settings.job_timeout_seconds - 1, queued["job_id"]),
        )
    service.store.release_job(second.job_id, "budget-worker-2")
    exhausted = service.store.get_job(queued["job_id"])
    assert exhausted is not None
    assert exhausted["status"] == "failed"
    assert exhausted["error_code"] == "job_timeout"


def test_database_probe_validates_rate_salt_and_rolls_back_write(settings_factory):
    store = Store(settings_factory())
    assert store.database_ok() is True
    with store.connect() as connection:
        readiness_row = connection.execute(
            "SELECT value FROM service_state WHERE name = '__readiness_probe__'"
        ).fetchone()
        connection.execute(
            "UPDATE metadata SET value = 'invalid' WHERE name = 'rate_salt'"
        )
    assert readiness_row is None
    assert store.database_ok() is False


def test_concurrent_first_start_serializes_schema_initialization(settings_factory):
    settings = settings_factory()
    barrier = threading.Barrier(8)
    errors = []

    def initialize_store():
        try:
            barrier.wait(timeout=5)
            Store(settings)
        except BaseException as error:
            errors.append(error)

    threads = [threading.Thread(target=initialize_store) for _ in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=20)

    assert all(not thread.is_alive() for thread in threads)
    assert errors == []
    store = Store(settings)
    assert store.database_ok() is True
    with store.connect() as connection:
        assert connection.execute("PRAGMA journal_mode").fetchone()[0] == "wal"


class _WalConnection:
    def __init__(self, outcome):
        self.outcome = outcome
        self.closed = False

    def execute(self, statement):
        assert statement == "PRAGMA journal_mode = WAL"
        if isinstance(self.outcome, BaseException):
            raise self.outcome
        return SimpleNamespace(fetchone=lambda: (self.outcome,))

    def close(self):
        self.closed = True


def _sqlite_error(message: str, code: int) -> sqlite3.OperationalError:
    error = sqlite3.OperationalError(message)
    error.sqlite_errorcode = code
    return error


def test_wal_negotiation_retries_only_sqlite_busy(monkeypatch):
    store = object.__new__(Store)
    outcomes = iter(
        (
            _sqlite_error("database is locked", sqlite3.SQLITE_BUSY),
            _sqlite_error("database is locked", sqlite3.SQLITE_BUSY),
            "wal",
        )
    )
    connections = []

    def connect():
        connection = _WalConnection(next(outcomes))
        connections.append(connection)
        return connection

    delays = []
    monkeypatch.setattr(store, "connect", connect)
    monkeypatch.setattr(db_module.time, "sleep", delays.append)

    store._enable_wal()

    assert delays == [0.005, 0.01]
    assert len(connections) == 3
    assert all(connection.closed for connection in connections)


def test_wal_negotiation_propagates_non_busy_error(monkeypatch):
    store = object.__new__(Store)
    expected = _sqlite_error("disk I/O error", sqlite3.SQLITE_IOERR)
    connection = _WalConnection(expected)
    monkeypatch.setattr(store, "connect", lambda: connection)

    with pytest.raises(sqlite3.OperationalError) as captured:
        store._enable_wal()

    assert captured.value is expected
    assert connection.closed is True


def test_wal_negotiation_rejects_non_wal_result(monkeypatch):
    store = object.__new__(Store)
    connection = _WalConnection("delete")
    monkeypatch.setattr(store, "connect", lambda: connection)

    with pytest.raises(RuntimeError, match="refused WAL"):
        store._enable_wal()

    assert connection.closed is True
