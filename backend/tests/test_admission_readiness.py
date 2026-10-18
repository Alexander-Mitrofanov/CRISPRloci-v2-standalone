from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import threading
import time

from fastapi.testclient import TestClient
import pytest

from crispr_web.api import create_app
from crispr_web.db import (
    ClientCapacityError,
    QueueCapacityError,
    Store,
    SubmissionRateLimitError,
)
from crispr_web.models import JobSubmission
from crispr_web.security import client_rate_identity
from crispr_web.service import JobService


def _submission(identifier: str) -> JobSubmission:
    return JobSubmission(sequence=f">{identifier}\nACGT\n", mode="detection")


def _create_direct(store: Store, index: int, client: str) -> None:
    store.create_job(
        job_id=f"{index:032x}",
        access_token_digest=f"{index + 100:064x}",
        mode="detection",
        filename=f"sample-{index}.fasta",
        record_count=1,
        base_count=4,
        input_sha256=f"{index + 200:064x}",
        source_ids=[f"sample-{index}"],
        category_policy="bona_fide_possible",
        spacer_distance=1,
        bias_corrections_requested=True,
        bias_corrections_effective=True,
        option_warnings=[],
        client_ip=client,
    )


def test_network_identity_normalizes_mapped_ipv4_and_ipv6_privacy_addresses() -> None:
    assert client_rate_identity("::ffff:192.0.2.9") == "192.0.2.9"
    assert client_rate_identity("192.0.2.9") == "192.0.2.9"
    assert client_rate_identity("2001:db8:1234:5678::1") == ("2001:db8:1234:5678::/64")
    assert client_rate_identity("2001:db8:1234:5678:ffff::abcd") == (
        "2001:db8:1234:5678::/64"
    )
    assert client_rate_identity("2001:db8:1234:5679::1") != (
        client_rate_identity("2001:db8:1234:5678::1")
    )


def test_job_persists_only_keyed_normalized_client_digest(settings_factory) -> None:
    settings = settings_factory(max_active_jobs_per_client=4)
    service = JobService(settings)
    _token, public = service.submit(
        _submission("privacy"), "2001:db8:abcd:42::dead:beef"
    )

    with service.store.connect() as connection:
        row = connection.execute(
            "SELECT client_digest FROM jobs WHERE job_id = ?", (public["job_id"],)
        ).fetchone()
        salt = connection.execute(
            "SELECT value FROM metadata WHERE name = 'rate_salt'"
        ).fetchone()["value"]
    assert row is not None
    digest = row["client_digest"]
    assert isinstance(digest, str) and len(digest) == 64
    assert "2001:db8" not in digest
    assert "client_digest" not in public

    import hashlib
    import hmac

    expected = hmac.new(
        bytes.fromhex(salt), b"2001:db8:abcd:42::/64", hashlib.sha256
    ).hexdigest()
    assert digest == expected


def test_ipv6_privacy_addresses_share_per_client_active_cap(settings_factory) -> None:
    service = JobService(
        settings_factory(
            max_active_jobs=4,
            max_active_jobs_per_client=1,
            max_queued_jobs=4,
            rate_limit_submissions=20,
        )
    )
    service.submit(_submission("first"), "2001:db8:1:2::1")
    with pytest.raises(ClientCapacityError):
        service.submit(_submission("same-network"), "2001:db8:1:2:ffff::abcd")
    service.submit(_submission("different-network"), "2001:db8:1:3::1")


def test_global_active_cap_counts_running_jobs(settings_factory) -> None:
    service = JobService(
        settings_factory(
            max_active_jobs=1,
            max_active_jobs_per_client=4,
            max_queued_jobs=4,
            rate_limit_submissions=20,
        )
    )
    service.submit(_submission("running"), "192.0.2.10")
    assert service.store.claim_next("worker") is not None
    with pytest.raises(QueueCapacityError):
        service.submit(_submission("blocked"), "192.0.2.11")


def test_per_client_cap_is_atomic_across_concurrent_creates(settings_factory) -> None:
    store = Store(
        settings_factory(
            max_active_jobs=4,
            max_active_jobs_per_client=1,
            max_queued_jobs=4,
            rate_limit_submissions=20,
        )
    )
    barrier = threading.Barrier(2)

    def create(index: int) -> object:
        barrier.wait(timeout=5)
        try:
            _create_direct(store, index, "198.51.100.50")
        except BaseException as error:
            return error
        return "created"

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(create, (1, 2)))
    assert results.count("created") == 1
    assert sum(isinstance(result, ClientCapacityError) for result in results) == 1
    with store.connect() as connection:
        active = connection.execute(
            "SELECT COUNT(*) FROM jobs WHERE status IN ('queued', 'running')"
        ).fetchone()[0]
    assert active == 1


def test_global_cap_is_atomic_during_distributed_submission_burst(
    settings_factory,
) -> None:
    store = Store(
        settings_factory(
            max_active_jobs=3,
            max_active_jobs_per_client=4,
            max_queued_jobs=10,
            rate_limit_submissions=20,
        )
    )
    barrier = threading.Barrier(12)

    def create(index: int) -> object:
        barrier.wait(timeout=5)
        try:
            _create_direct(store, index, f"198.51.100.{index + 1}")
        except BaseException as error:
            return error
        return "created"

    with ThreadPoolExecutor(max_workers=12) as executor:
        results = list(executor.map(create, range(12)))
    assert results.count("created") == 3
    assert sum(isinstance(result, QueueCapacityError) for result in results) == 9


def test_rate_limit_is_atomic_during_same_client_submission_burst(
    settings_factory,
) -> None:
    store = Store(
        settings_factory(
            max_active_jobs=32,
            max_active_jobs_per_client=32,
            max_queued_jobs=32,
            rate_limit_submissions=3,
        )
    )
    barrier = threading.Barrier(12)

    def create(index: int) -> object:
        barrier.wait(timeout=5)
        try:
            _create_direct(store, index, "203.0.113.44")
        except BaseException as error:
            return error
        return "created"

    with ThreadPoolExecutor(max_workers=12) as executor:
        results = list(executor.map(create, range(12)))
    assert results.count("created") == 3
    assert sum(isinstance(result, SubmissionRateLimitError) for result in results) == 9


def test_additive_digest_migration_preserves_legacy_insert_compatibility(
    settings_factory,
) -> None:
    settings = settings_factory(
        max_active_jobs=1,
        max_active_jobs_per_client=1,
        max_queued_jobs=4,
        rate_limit_submissions=20,
    )
    store = Store(settings)
    with store.connect() as connection:
        connection.execute("DROP INDEX jobs_client_active_idx")
        connection.execute("ALTER TABLE jobs DROP COLUMN client_digest")

    migrated = Store(settings)
    now = time.time()
    with migrated.connect() as connection:
        columns = {row["name"] for row in connection.execute("PRAGMA table_info(jobs)")}
        schema_version = connection.execute(
            "SELECT value FROM metadata WHERE name = 'schema_version'"
        ).fetchone()["value"]
        indexes = {row["name"] for row in connection.execute("PRAGMA index_list(jobs)")}
        connection.execute(
            """
            INSERT INTO jobs(
                job_id, token_digest, status, mode, created_at, updated_at,
                expires_at, filename, record_count, base_count, input_sha256,
                source_ids_json, category_policy, spacer_distance,
                bias_corrections_requested, bias_corrections_effective,
                option_warnings_json
            ) VALUES (
                :job_id, :token_digest, 'queued', 'detection', :now, :now,
                :expires_at, 'legacy.fasta', 1, 4, :input_sha256,
                '["legacy"]', 'bona_fide_possible', 1, 1, 1, '[]'
            )
            """,
            {
                "job_id": "f" * 32,
                "token_digest": "a" * 64,
                "now": now,
                "expires_at": now + 3600,
                "input_sha256": "b" * 64,
            },
        )
        legacy_digest = connection.execute(
            "SELECT client_digest FROM jobs WHERE job_id = ?", ("f" * 32,)
        ).fetchone()["client_digest"]
    assert "client_digest" in columns
    assert "jobs_client_active_idx" in indexes
    assert schema_version == "7"
    assert legacy_digest is None
    assert Store(settings).database_ok() is True
    with pytest.raises(QueueCapacityError):
        migrated.preflight_submission("203.0.113.25", now=now)


def test_public_health_is_single_flight_cached_and_fresh_is_local_only(
    settings_factory, monkeypatch
) -> None:
    settings = settings_factory(readiness_cache_seconds=60, min_free_bytes=1)
    app = create_app(settings)
    entered = threading.Event()
    release = threading.Event()
    calls = 0
    calls_lock = threading.Lock()

    def database_ok() -> bool:
        nonlocal calls
        with calls_lock:
            calls += 1
        entered.set()
        assert release.wait(timeout=5)
        return True

    monkeypatch.setattr(app.state.store, "database_ok", database_ok)
    monkeypatch.setattr(app.state.store, "worker_health", lambda: "ok")
    with TestClient(app, client=("127.0.0.1", 50000)) as client:
        with ThreadPoolExecutor(max_workers=8) as executor:
            futures = [executor.submit(client.get, "/api/v1/health") for _ in range(8)]
            assert entered.wait(timeout=3)
            release.set()
            responses = [future.result(timeout=5) for future in futures]
        assert all(response.status_code == 200 for response in responses)
        assert calls == 1
        assert client.get("/api/v1/health").status_code == 200
        assert calls == 1

        public_fresh = client.get(
            "/api/v1/health/fresh", headers={"X-Forwarded-For": "192.0.2.1"}
        )
        assert public_fresh.status_code == 404
        assert calls == 1
        assert client.get("/api/v1/health/fresh").status_code == 200
        assert calls == 2
        assert (
            "/api/v1/health/fresh"
            not in client.get("/api/v1/openapi.json").json()["paths"]
        )


def test_client_active_cap_has_actionable_public_response(settings_factory) -> None:
    app = create_app(
        settings_factory(
            max_active_jobs=4,
            max_active_jobs_per_client=1,
            max_queued_jobs=4,
            rate_limit_submissions=20,
        )
    )
    with TestClient(app, client=("192.0.2.55", 50000)) as client:
        first = client.post(
            "/api/v1/jobs", json={"sequence": ">first\nACGT\n", "mode": "detection"}
        )
        second = client.post(
            "/api/v1/jobs", json={"sequence": ">second\nACGT\n", "mode": "detection"}
        )
    assert first.status_code == 202
    assert second.status_code == 429
    assert second.headers["retry-after"] == "60"
    assert "maximum number of active jobs" in second.json()["detail"]
