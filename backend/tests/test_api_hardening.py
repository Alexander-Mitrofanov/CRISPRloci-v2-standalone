from __future__ import annotations

import json

from fastapi.testclient import TestClient
import pytest

from crispr_web.api import create_app
from crispr_web.db import Store
from crispr_web.service import JobService


@pytest.mark.parametrize(
    "field,value",
    [
        ("bias_corrections", "false"),
        ("bias_corrections", 0),
        ("bias_corrections", 1),
        ("spacer_distance", "1"),
        ("spacer_distance", 1.0),
        ("spacer_distance", True),
    ],
)
def test_public_options_reject_coerced_json_types(settings_factory, field, value):
    app = create_app(settings_factory())
    with TestClient(app) as client:
        response = client.post(
            "/api/v1/jobs",
            json={"sequence": ">sample\nACGT\n", "mode": "detection", field: value},
        )
    assert response.status_code == 422


def test_large_validation_error_never_echoes_raw_sequence(settings_factory):
    app = create_app(settings_factory(max_request_bytes=1_000_000))
    secret_sequence = "SENSITIVE_SEQUENCE_" * 25_000
    with TestClient(app) as client:
        response = client.post(
            "/api/v1/jobs",
            headers={"Origin": "https://frontend.example"},
            json={"sequence": [secret_sequence], "mode": "detection"},
        )
    assert response.status_code == 422
    assert len(response.content) < 10_000
    assert "SENSITIVE_SEQUENCE" not in response.text
    assert set(response.json()["detail"][0]) == {"loc", "msg", "type"}
    assert response.headers["access-control-allow-origin"] == "https://frontend.example"
    assert response.headers["x-content-type-options"] == "nosniff"


def test_preflight_and_unhandled_errors_keep_cors_and_security_headers(
    settings_factory, monkeypatch
):
    app = create_app(settings_factory())
    with TestClient(app, raise_server_exceptions=False) as client:
        preflight = client.options(
            "/api/v1/jobs",
            headers={
                "Origin": "https://frontend.example",
                "Access-Control-Request-Method": "POST",
                "Access-Control-Request-Headers": "content-type,authorization",
            },
        )
        assert preflight.status_code == 200
        assert (
            preflight.headers["access-control-allow-origin"]
            == "https://frontend.example"
        )
        assert preflight.headers["x-frame-options"] == "DENY"

        def explode(_self, _job_id, _token):
            raise RuntimeError("private backend failure")

        monkeypatch.setattr(JobService, "authorize", explode)
        failed = client.get(
            f"/api/v1/jobs/{'0' * 32}",
            headers={
                "Origin": "https://frontend.example",
                "Authorization": "Bearer test",
            },
        )
        assert failed.status_code == 500
        assert failed.json() == {
            "detail": "the backend could not complete this request"
        }
        assert (
            failed.headers["access-control-allow-origin"] == "https://frontend.example"
        )
        assert failed.headers["x-content-type-options"] == "nosniff"


def test_health_probe_failures_are_structured_degraded_responses(
    settings_factory, monkeypatch
):
    app = create_app(settings_factory())

    def fail_database(_self):
        raise RuntimeError("database path")

    def fail_worker(_self, now=None):
        raise RuntimeError("worker table")

    monkeypatch.setattr(Store, "database_ok", fail_database)
    monkeypatch.setattr(Store, "worker_health", fail_worker)
    with TestClient(app) as client:
        response = client.get("/api/v1/health")
    assert response.status_code == 503
    assert response.json()["database"] == "error"
    assert response.json()["worker"] == "error"
    assert len(json.dumps(response.json())) < 3_000


def test_manifest_tampering_degrades_readiness_without_exposing_paths(settings_factory):
    settings = settings_factory()
    settings.scientific_release_manifest.write_text("{}\n", encoding="utf-8")
    app = create_app(settings)
    with TestClient(app) as client:
        health = client.get("/api/v1/health")
        version = client.get("/api/v1/version")

    assert health.status_code == 503
    assert health.json()["scientific_release"] == "error"
    assert (
        health.json()["scientific_release_manifest_sha256"]
        == settings.scientific_release_manifest_sha256
    )
    identities = health.json()["scientific_artifacts"]
    assert set(identities) == {
        "crispridentify",
        "spacerplacer",
        "crispr_evor",
        "integration",
    }
    assert all("path" not in identity for identity in identities.values())
    assert version.json()["scientific_artifacts"] == identities


def test_low_storage_blocks_submission_and_degrades_health(settings_factory):
    settings = settings_factory(min_free_bytes=100_000_000_000)
    app = create_app(settings)
    with TestClient(app) as client:
        health = client.get("/api/v1/health")
        submission = client.post(
            "/api/v1/jobs",
            json={"sequence": ">sample\nACGT\n", "mode": "detection"},
        )
    assert health.status_code == 503
    assert health.json()["storage"] == "low"
    assert submission.status_code == 503
    assert submission.headers["retry-after"] == "300"
