from __future__ import annotations

import asyncio
import ipaddress

from fastapi import Request
from fastapi.testclient import TestClient

from crispr_web.api import RequestLimitMiddleware, create_app
from crispr_web.security import client_ip
from crispr_web.worker import Worker


def test_versioned_api_cors_security_and_capability_auth(settings_factory):
    settings = settings_factory(
        max_request_bytes=4_321,
        max_archive_bytes=87_654_321,
        max_header_characters=321,
    )
    app = create_app(settings)
    with TestClient(app) as client:
        health = client.get("/api/v1/health")
        assert health.status_code == 503
        assert health.json()["worker"] == "not_started"
        assert health.json()["release_id"] == settings.release_id
        assert health.json()["backend_wheel_sha256"] == settings.backend_wheel_sha256
        assert health.json()["scientific_release"] == "ok"
        assert (
            health.json()["scientific_release_manifest_sha256"]
            == settings.scientific_release_manifest_sha256
        )
        assert health.json()["scientific_artifacts"] == (
            settings.scientific_release_identities
        )
        version = client.get("/api/v1/version")
        assert version.status_code == 200
        assert version.json() == {
            "service": "crisprloci-v2",
            "version": settings.service_version,
            "api_version": "v1",
            "integration_contract": "crispr-tool-integration-run/1.0.0",
            "release_id": settings.release_id,
            "backend_wheel_sha256": settings.backend_wheel_sha256,
            "scientific_release_manifest_sha256": (
                settings.scientific_release_manifest_sha256
            ),
            "scientific_artifacts": settings.scientific_release_identities,
        }
        assert health.headers["x-content-type-options"] == "nosniff"
        assert health.headers["x-frame-options"] == "DENY"

        config = client.get(
            "/api/v1/config", headers={"Origin": "https://frontend.example"}
        )
        assert config.status_code == 200
        assert (
            config.headers["access-control-allow-origin"] == "https://frontend.example"
        )
        assert config.json()["category_policies"] == [
            "bona_fide_possible",
            "bona_fide_only",
        ]
        assert config.json()["identify_folder_jobs"] == 1
        assert config.json()["max_request_bytes"] == settings.max_request_bytes
        assert config.json()["max_archive_bytes"] == settings.max_archive_bytes
        assert config.json()["max_header_characters"] == settings.max_header_characters
        denied_origin = client.get(
            "/api/v1/config", headers={"Origin": "https://attacker.example"}
        )
        assert "access-control-allow-origin" not in denied_origin.headers

        created = client.post(
            "/api/v1/jobs",
            headers={"Origin": "https://frontend.example"},
            json={"sequence": ">sample\nACGTN\n", "mode": "detection"},
        )
        assert created.status_code == 202
        payload = created.json()
        job_id = payload["job"]["job_id"]
        token = payload["access_token"]
        assert len(token) == 43
        assert payload["job"]["expires_at"] is not None

        assert client.get(f"/api/v1/jobs/{job_id}").status_code == 401
        assert (
            client.get(
                f"/api/v1/jobs/{job_id}",
                headers={"Authorization": "Bearer wrong"},
            ).status_code
            == 404
        )
        authorized = client.get(
            f"/api/v1/jobs/{job_id}",
            headers={"Authorization": f"Bearer {token}"},
        )
        assert authorized.status_code == 200
        assert authorized.json()["input"]["source_ids"] == ["sample"]

        result = client.get(
            f"/api/v1/jobs/{job_id}/result",
            headers={"Authorization": f"Bearer {token}"},
        )
        assert result.status_code == 404

        cancelled = client.delete(
            f"/api/v1/jobs/{job_id}",
            headers={"Authorization": f"Bearer {token}"},
        )
        assert cancelled.status_code == 200
        assert cancelled.json()["job"]["status"] == "cancelled"

        openapi = client.get("/api/v1/openapi.json")
        assert openapi.status_code == 200
        assert all(path.startswith("/api/v1/") for path in openapi.json()["paths"])
        assert client.get("/api/v1/docs").status_code == 404


def test_submission_model_forbids_arbitrary_arguments(settings_factory):
    app = create_app(settings_factory())
    with TestClient(app) as client:
        response = client.post(
            "/api/v1/jobs",
            json={
                "sequence": ">sample\nACGT\n",
                "mode": "orientation",
                "command": ["/bin/sh", "-c", "anything"],
            },
        )
        assert response.status_code == 422


def test_authenticated_result_bundle_is_always_attachment(settings_factory):
    settings = settings_factory()
    app = create_app(settings)
    with TestClient(app) as client:
        created = client.post(
            "/api/v1/jobs",
            json={"sequence": ">a\nACGT\n>b\nACGA\n", "mode": "orientation"},
        ).json()
        job_id = created["job"]["job_id"]
        token = created["access_token"]
        Worker(settings, store=app.state.store, worker_id="api-worker").run_once()
        response = client.get(
            f"/api/v1/jobs/{job_id}/result",
            headers={"Authorization": f"Bearer {token}"},
        )
        assert response.status_code == 200
        assert response.headers["content-type"] == "application/zip"
        assert response.headers["content-disposition"].startswith("attachment;")
        assert response.content.startswith(b"PK")


def test_request_limit_counts_received_bytes_not_declared_length():
    downstream_called = False

    async def downstream(_scope, _receive, _send):
        nonlocal downstream_called
        downstream_called = True

    middleware = RequestLimitMiddleware(downstream, max_bytes=8)
    messages = iter(
        [
            {
                "type": "http.request",
                "body": b"0123456789",
                "more_body": False,
            }
        ]
    )
    sent = []

    async def receive():
        return next(messages)

    async def send(message):
        sent.append(message)

    scope = {
        "type": "http",
        "method": "POST",
        "path": "/api/v1/jobs",
        "headers": [(b"content-length", b"1")],
    }
    asyncio.run(middleware(scope, receive, send))
    assert downstream_called is False
    assert sent[0]["type"] == "http.response.start"
    assert sent[0]["status"] == 413


def test_request_limit_rejects_ambiguous_or_inconsistent_framing():
    async def exercise(headers, body=b"{}"):
        downstream_called = False

        async def downstream(_scope, _receive, _send):
            nonlocal downstream_called
            downstream_called = True

        middleware = RequestLimitMiddleware(downstream, max_bytes=32)
        messages = iter([{"type": "http.request", "body": body, "more_body": False}])
        sent = []

        async def receive():
            return next(messages)

        async def send(message):
            sent.append(message)

        scope = {
            "type": "http",
            "method": "POST",
            "path": "/api/v1/jobs",
            "headers": headers,
        }
        await middleware(scope, receive, send)
        return downstream_called, sent

    cases = (
        (
            [(b"content-length", b"2"), (b"content-length", b"2")],
            b"{}",
        ),
        (
            [(b"content-length", b"2"), (b"transfer-encoding", b"chunked")],
            b"{}",
        ),
        ([(b"content-length", b"3")], b"{}"),
    )
    for headers, body in cases:
        downstream_called, sent = asyncio.run(exercise(headers, body))
        assert downstream_called is False
        assert sent[0]["type"] == "http.response.start"
        assert sent[0]["status"] == 400


def test_request_limit_drops_disconnected_partial_submission():
    downstream_called = False

    async def downstream(_scope, _receive, _send):
        nonlocal downstream_called
        downstream_called = True

    middleware = RequestLimitMiddleware(downstream, max_bytes=32)
    messages = iter(
        [
            {"type": "http.request", "body": b"{", "more_body": True},
            {"type": "http.disconnect"},
        ]
    )
    sent = []

    async def receive():
        return next(messages)

    async def send(message):
        sent.append(message)

    scope = {
        "type": "http",
        "method": "POST",
        "path": "/api/v1/jobs",
        "headers": [(b"content-length", b"2")],
    }
    asyncio.run(middleware(scope, receive, send))
    assert downstream_called is False
    assert sent == []


def test_forwarded_chain_is_trusted_right_to_left():
    def request(peer: str, forwarded: str) -> Request:
        return Request(
            {
                "type": "http",
                "client": (peer, 12345),
                "headers": [(b"x-forwarded-for", forwarded.encode("ascii"))],
            }
        )

    loopback = ipaddress.ip_network("127.0.0.1/32")
    private_proxy = ipaddress.ip_network("10.0.0.0/8")
    chain = request("127.0.0.1", "198.51.100.9, 10.1.2.3")
    assert client_ip(chain, (loopback, private_proxy)) == "198.51.100.9"
    assert client_ip(chain, (loopback,)) == "10.1.2.3"

    untrusted = request("203.0.113.7", "198.51.100.9")
    assert client_ip(untrusted, (loopback, private_proxy)) == "203.0.113.7"
    invalid = request("127.0.0.1", "not-an-ip")
    assert client_ip(invalid, (loopback, private_proxy)) == "127.0.0.1"
