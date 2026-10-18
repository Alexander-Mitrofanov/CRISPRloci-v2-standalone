from __future__ import annotations

import time

import pytest

from crispr_web.db import Store
from crispr_web.models import JobSubmission
from crispr_web.service import JobService


def _submission(identifier: str) -> JobSubmission:
    return JobSubmission(sequence=f">{identifier}\nACGT\n", mode="detection")


def _identity_row(service: JobService, job_id: str) -> dict[str, object]:
    with service.store.connect() as connection:
        row = connection.execute(
            "SELECT status, client_digest, cancel_requested FROM jobs WHERE job_id = ?",
            (job_id,),
        ).fetchone()
    assert row is not None
    return dict(row)


def test_scheduled_cleanup_purges_idle_rate_events_at_exact_window_boundary(
    settings_factory, monkeypatch
) -> None:
    settings = settings_factory(rate_limit_window_seconds=60)
    service = JobService(settings)
    now = 10_000.0
    cutoff = now - settings.rate_limit_window_seconds
    with service.store.connect() as connection:
        connection.executemany(
            "INSERT INTO rate_events(client_digest, submitted_at) VALUES(?, ?)",
            [
                ("expired", cutoff - 1),
                ("boundary", cutoff),
                ("current", cutoff + 1),
            ],
        )

    monkeypatch.setattr("crispr_web.service.time.time", lambda: now)
    assert service.cleanup_expired() == 0

    with service.store.connect() as connection:
        rows = connection.execute(
            "SELECT client_digest, submitted_at FROM rate_events ORDER BY submitted_at"
        ).fetchall()
    assert [(row["client_digest"], row["submitted_at"]) for row in rows] == [
        ("boundary", cutoff),
        ("current", cutoff + 1),
    ]


@pytest.mark.parametrize(
    "terminal_status",
    ["completed", "completed_no_eligible_groups", "failed", "cancelled"],
)
def test_finish_erases_client_digest_for_every_terminal_status(
    settings_factory, terminal_status
) -> None:
    service = JobService(settings_factory())
    _token, job = service.submit(_submission(f"finish-{terminal_status}"), "192.0.2.80")
    claim = service.store.claim_next("finish-worker")
    assert claim is not None
    assert _identity_row(service, job["job_id"])["client_digest"] is not None

    service.store.finish_job(
        job_id=job["job_id"],
        worker_id="finish-worker",
        status=terminal_status,
        summary=None,
        artifacts=[],
    )

    row = _identity_row(service, job["job_id"])
    assert row["status"] == terminal_status
    assert row["client_digest"] is None


def test_cancel_retains_running_identity_until_release_then_erases_it(
    settings_factory,
) -> None:
    service = JobService(settings_factory())

    _token, queued = service.submit(_submission("queued-cancel"), "192.0.2.81")
    service.cancel(queued["job_id"])
    queued_row = _identity_row(service, queued["job_id"])
    assert queued_row["status"] == "cancelled"
    assert queued_row["client_digest"] is None

    _token, running = service.submit(_submission("running-cancel"), "192.0.2.82")
    claim = service.store.claim_next("cancel-worker")
    assert claim is not None
    service.cancel(running["job_id"])
    requested = _identity_row(service, running["job_id"])
    assert requested["status"] == "running"
    assert requested["cancel_requested"] == 1
    assert requested["client_digest"] is not None

    service.store.release_job(running["job_id"], "cancel-worker")
    released = _identity_row(service, running["job_id"])
    assert released["status"] == "cancelled"
    assert released["client_digest"] is None


@pytest.mark.parametrize(
    ("recovery_case", "expected_status"),
    [
        ("cancelled", "cancelled"),
        ("timed_out", "failed"),
        ("attempts_exhausted", "failed"),
    ],
)
def test_stale_recovery_erases_terminal_client_digest(
    settings_factory, recovery_case, expected_status
) -> None:
    service = JobService(settings_factory(max_job_attempts=1))
    _token, job = service.submit(_submission(f"recover-{recovery_case}"), "192.0.2.83")
    base = time.time()
    claim = service.store.claim_next("recovery-worker", now=base)
    assert claim is not None

    if recovery_case == "cancelled":
        service.cancel(job["job_id"])
    elif recovery_case == "timed_out":
        with service.store.connect() as connection:
            connection.execute(
                "UPDATE jobs SET started_at = ? WHERE job_id = ?",
                (base - 100, job["job_id"]),
            )

    assert service.store.recover_stale(now=base + 11) == 1
    row = _identity_row(service, job["job_id"])
    assert row["status"] == expected_status
    assert row["client_digest"] is None


@pytest.mark.parametrize(
    ("queue_case", "expected_status"),
    [
        ("expired", "cancelled"),
        ("timed_out", "failed"),
        ("attempts_exhausted", "failed"),
    ],
)
def test_claim_preflight_erases_terminal_queued_client_digest(
    settings_factory, queue_case, expected_status
) -> None:
    settings = settings_factory(max_job_attempts=1)
    service = JobService(settings)
    _token, job = service.submit(_submission(f"queue-{queue_case}"), "192.0.2.84")
    now = time.time()
    with service.store.connect() as connection:
        if queue_case == "expired":
            connection.execute(
                "UPDATE jobs SET expires_at = ? WHERE job_id = ?",
                (now - 1, job["job_id"]),
            )
        elif queue_case == "timed_out":
            connection.execute(
                "UPDATE jobs SET started_at = ? WHERE job_id = ?",
                (now - settings.job_timeout_seconds - 1, job["job_id"]),
            )
        else:
            connection.execute(
                "UPDATE jobs SET attempt = ? WHERE job_id = ?",
                (settings.max_job_attempts, job["job_id"]),
            )

    assert service.store.claim_next("queue-worker", now=now) is None
    row = _identity_row(service, job["job_id"])
    assert row["status"] == expected_status
    assert row["client_digest"] is None


@pytest.mark.parametrize("release_case", ["timed_out", "attempts_exhausted"])
def test_worker_release_erases_failed_client_digest(
    settings_factory, release_case
) -> None:
    max_attempts = 1 if release_case == "attempts_exhausted" else 2
    service = JobService(settings_factory(max_job_attempts=max_attempts))
    _token, job = service.submit(_submission(f"release-{release_case}"), "192.0.2.85")
    claim = service.store.claim_next("release-worker")
    assert claim is not None
    if release_case == "timed_out":
        with service.store.connect() as connection:
            connection.execute(
                "UPDATE jobs SET started_at = 0 WHERE job_id = ?",
                (job["job_id"],),
            )

    service.store.release_job(job["job_id"], "release-worker")
    row = _identity_row(service, job["job_id"])
    assert row["status"] == "failed"
    assert row["client_digest"] is None


def test_freezing_expired_queued_job_erases_client_digest(settings_factory) -> None:
    service = JobService(settings_factory())
    _token, job = service.submit(_submission("freeze-expired"), "192.0.2.86")
    now = time.time()
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET expires_at = ? WHERE job_id = ?",
            (now - 1, job["job_id"]),
        )

    assert service.store.freeze_expired_job(job["job_id"], now=now) is True
    row = _identity_row(service, job["job_id"])
    assert row["status"] == "cancelled"
    assert row["client_digest"] is None


def test_initialize_scrubs_existing_terminal_client_digests(settings_factory) -> None:
    settings = settings_factory()
    service = JobService(settings)
    _token, job = service.submit(_submission("upgrade-scrub"), "192.0.2.87")
    with service.store.connect() as connection:
        connection.execute(
            "UPDATE jobs SET status = 'failed' WHERE job_id = ?",
            (job["job_id"],),
        )
    assert _identity_row(service, job["job_id"])["client_digest"] is not None

    Store(settings)

    assert _identity_row(service, job["job_id"])["client_digest"] is None
