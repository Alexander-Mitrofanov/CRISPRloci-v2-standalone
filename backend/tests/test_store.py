from __future__ import annotations

import time

import pytest

from crispr_web.db import QueueCapacityError, Store, SubmissionRateLimitError
from crispr_web.models import JobSubmission
from crispr_web.service import JobService


def test_token_digest_only_queue_capacity_and_rate_limit(settings_factory):
    settings = settings_factory(max_queued_jobs=1, rate_limit_submissions=1)
    service = JobService(settings)
    token, first = service.submit(
        JobSubmission(sequence=">first\nACGT\n", mode="detection"),
        "192.0.2.10",
    )
    with service.store.connect() as connection:
        digest = connection.execute("SELECT token_digest FROM jobs").fetchone()[0]
        assert token not in digest
        assert len(digest) == 64
    assert service.authorize(first["job_id"], token)
    assert not service.authorize(first["job_id"], "wrong")
    assert first["expires_at"] is not None

    with pytest.raises(SubmissionRateLimitError):
        service.submit(
            JobSubmission(sequence=">rate\nACGT\n", mode="detection"),
            "192.0.2.10",
        )
    with pytest.raises(QueueCapacityError):
        service.submit(
            JobSubmission(sequence=">queue\nACGT\n", mode="detection"),
            "192.0.2.11",
        )


def test_stale_jobs_retry_once_then_fail_and_stale_cancel_is_terminal(settings_factory):
    settings = settings_factory(max_job_attempts=2)
    base = time.time()
    service = JobService(settings)
    _, first = service.submit(
        JobSubmission(sequence=">retry\nACGT\n", mode="reconstruction"),
        "192.0.2.1",
    )
    claim = service.store.claim_next("worker-a", now=base)
    assert claim and claim.attempt == 1
    assert service.store.recover_stale(now=base + 11) == 1
    second = service.store.claim_next("worker-b", now=base + 12)
    assert second and second.attempt == 2
    assert service.store.recover_stale(now=base + 23) == 1
    exhausted = service.get(first["job_id"])
    assert exhausted["status"] == "failed"
    assert exhausted["error"]["code"] == "worker_restarts_exhausted"
    assert service.store.claim_next("worker-c", now=base + 24) is None

    _, cancel_job = service.submit(
        JobSubmission(sequence=">cancel\nACGT\n", mode="detection"),
        "192.0.2.2",
    )
    cancel_claim = service.store.claim_next("worker-d", now=base + 100)
    assert cancel_claim and cancel_claim.job_id == cancel_job["job_id"]
    service.cancel(cancel_job["job_id"])
    assert service.store.recover_stale(now=base + 111) == 1
    assert service.get(cancel_job["job_id"])["status"] == "cancelled"


def test_rate_event_expiry_index_is_created_idempotently(settings_factory):
    settings = settings_factory()
    Store(settings)
    Store(settings)
    with Store(settings).connect() as connection:
        indexes = {
            row["name"] for row in connection.execute("PRAGMA index_list(rate_events)")
        }
        assert "rate_events_lookup_idx" in indexes
        assert "rate_events_expiry_idx" in indexes
        plan = connection.execute(
            "EXPLAIN QUERY PLAN DELETE FROM rate_events WHERE submitted_at < ?",
            (0,),
        ).fetchall()
    assert any("rate_events_expiry_idx" in row["detail"] for row in plan)


def test_graceful_release_cannot_bypass_attempt_ceiling(settings_factory):
    service = JobService(settings_factory(max_job_attempts=2))
    _, queued = service.submit(
        JobSubmission(sequence=">bounded\nACGT\n", mode="detection"),
        "192.0.2.20",
    )
    base = time.time()
    first = service.store.claim_next("worker-a", now=base)
    assert first and first.attempt == 1
    service.store.release_job(first.job_id, "worker-a")
    second = service.store.claim_next("worker-b", now=base + 1)
    assert second and second.attempt == 2
    service.store.release_job(second.job_id, "worker-b")

    job = service.get(queued["job_id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "worker_restarts_exhausted"
    assert service.store.claim_next("worker-c", now=base + 2) is None


def test_expired_queued_job_and_private_directory_are_deleted(settings_factory):
    service = JobService(settings_factory())
    _, queued = service.submit(
        JobSubmission(sequence=">expires\nACGT\n", mode="detection"),
        "192.0.2.21",
    )
    root = service.job_root(queued["job_id"])
    assert root.is_dir()
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET expires_at = 0 WHERE job_id = ?",
            (queued["job_id"],),
        )

    assert service.cleanup_expired() == 1
    assert not root.exists()
    with pytest.raises(KeyError):
        service.get(queued["job_id"])


def test_cleanup_recovers_and_deletes_expired_stale_running_job(settings_factory):
    service = JobService(settings_factory())
    _, queued = service.submit(
        JobSubmission(sequence=">stale-running\nACGT\n", mode="detection"),
        "192.0.2.22",
    )
    claim = service.store.claim_next("unavailable-worker")
    assert claim is not None
    root = service.job_root(queued["job_id"])
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET lease_expires_at = 0, expires_at = 0 WHERE job_id = ?",
            (queued["job_id"],),
        )

    assert service.cleanup_expired() == 1
    assert not root.exists()
    assert service.store.get_job(queued["job_id"]) is None


def test_cleanup_preserves_running_job_with_live_lease(settings_factory):
    service = JobService(settings_factory())
    _, queued = service.submit(
        JobSubmission(sequence=">live-running\nACGT\n", mode="detection"),
        "192.0.2.23",
    )
    claim = service.store.claim_next("live-worker")
    assert claim is not None
    root = service.job_root(queued["job_id"])
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET lease_expires_at = ?, expires_at = 0 WHERE job_id = ?",
            (time.time() + 60, queued["job_id"]),
        )

    assert service.cleanup_expired() == 0
    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "running"
    assert root.is_dir()
