from __future__ import annotations

import json

from fastapi.testclient import TestClient

from crispr_web.api import create_app
from crispr_web.worker import Worker


def _submission(identifier: str, suffix: str) -> dict[str, object]:
    return {
        "sequence": f">{identifier}_A\nACGT{suffix}\n>{identifier}_B\nACGA{suffix}\n",
        "filename": f"{identifier}.fasta",
        "mode": "orientation",
        "category_policy": "bona_fide_possible",
        "spacer_distance": 1,
        "bias_corrections": True,
    }


def _bearer(token: str) -> dict[str, str]:
    return {"Authorization": f"Bearer {token}"}


def test_two_anonymous_users_have_independent_queued_jobs_results_and_artifacts(
    settings_factory,
) -> None:
    settings = settings_factory(
        max_queued_jobs=4,
        max_active_jobs=4,
        max_active_jobs_per_client=2,
        rate_limit_submissions=4,
    )
    app = create_app(settings)

    # Separate TestClient instances model independent browser sessions. The API
    # deliberately does not depend on cookies: possession of the per-job bearer
    # capability is the complete anonymous authorization boundary.
    with TestClient(app) as user_a, TestClient(app) as user_b:
        created_a = user_a.post("/api/v1/jobs", json=_submission("user_a", "T"))
        created_b = user_b.post("/api/v1/jobs", json=_submission("user_b", "C"))
        assert created_a.status_code == created_b.status_code == 202
        assert created_a.headers["cache-control"] == "no-store"
        assert created_b.headers["cache-control"] == "no-store"
        assert "set-cookie" not in created_a.headers
        assert "set-cookie" not in created_b.headers

        payload_a = created_a.json()
        payload_b = created_b.json()
        job_a = payload_a["job"]["job_id"]
        job_b = payload_b["job"]["job_id"]
        token_a = payload_a["access_token"]
        token_b = payload_b["access_token"]
        assert job_a != job_b
        assert token_a != token_b
        assert job_b not in json.dumps(payload_a)
        assert job_a not in json.dumps(payload_b)
        assert token_b not in json.dumps(payload_a)
        assert token_a not in json.dumps(payload_b)
        assert user_a.cookies == user_b.cookies == {}

        # No endpoint enumerates anonymous jobs, and a missing token is distinct
        # only from authentication failure—not from another user's locator.
        assert user_a.get("/api/v1/jobs").status_code == 405
        assert user_a.get(f"/api/v1/jobs/{job_a}").status_code == 401
        assert user_b.get(f"/api/v1/jobs/{job_b}").status_code == 401

        for client, foreign_job, foreign_token in (
            (user_a, job_b, token_a),
            (user_b, job_a, token_b),
        ):
            denied = client.get(
                f"/api/v1/jobs/{foreign_job}", headers=_bearer(foreign_token)
            )
            assert denied.status_code == 404
            assert denied.json() == {"detail": "job not found"}
            denied_cancel = client.delete(
                f"/api/v1/jobs/{foreign_job}", headers=_bearer(foreign_token)
            )
            assert denied_cancel.status_code == 404
            assert denied_cancel.json() == {"detail": "job not found"}

        # Both jobs coexist in distinct private roots and are independently
        # claimed and completed by the durable queue.
        assert app.state.service.job_root(job_a) != app.state.service.job_root(job_b)
        assert app.state.service.job_root(job_a).is_dir()
        assert app.state.service.job_root(job_b).is_dir()
        worker = Worker(settings, store=app.state.store, worker_id="two-user-worker")
        assert worker.run_once() is True
        assert worker.run_once() is True

        own_a = user_a.get(f"/api/v1/jobs/{job_a}", headers=_bearer(token_a))
        own_b = user_b.get(f"/api/v1/jobs/{job_b}", headers=_bearer(token_b))
        assert own_a.status_code == own_b.status_code == 200
        assert own_a.json()["status"] == own_b.json()["status"] == "completed"
        assert own_a.json()["input"]["source_ids"] == ["user_a_A", "user_a_B"]
        assert own_b.json()["input"]["source_ids"] == ["user_b_A", "user_b_B"]
        assert "user_b" not in json.dumps(own_a.json())
        assert "user_a" not in json.dumps(own_b.json())
        assert token_a not in own_a.text
        assert token_b not in own_b.text
        assert own_a.headers["cache-control"] == own_b.headers["cache-control"] == "no-store"

        artifact_a = own_a.json()["artifacts"][0]["artifact_id"]
        artifact_b = own_b.json()["artifacts"][0]["artifact_id"]
        assert user_a.get(
            f"/api/v1/jobs/{job_b}/artifacts/{artifact_b}", headers=_bearer(token_a)
        ).status_code == 404
        assert user_b.get(
            f"/api/v1/jobs/{job_a}/artifacts/{artifact_a}", headers=_bearer(token_b)
        ).status_code == 404
        own_artifact_a = user_a.get(
            f"/api/v1/jobs/{job_a}/artifacts/{artifact_a}", headers=_bearer(token_a)
        )
        own_artifact_b = user_b.get(
            f"/api/v1/jobs/{job_b}/artifacts/{artifact_b}", headers=_bearer(token_b)
        )
        assert own_artifact_a.status_code == own_artifact_b.status_code == 200
        assert own_artifact_a.headers["content-disposition"].startswith("attachment;")
        assert own_artifact_b.headers["content-disposition"].startswith("attachment;")

        assert user_a.get(
            f"/api/v1/jobs/{job_b}/result", headers=_bearer(token_a)
        ).status_code == 404
        assert user_b.get(
            f"/api/v1/jobs/{job_a}/result", headers=_bearer(token_b)
        ).status_code == 404
        own_bundle_a = user_a.get(
            f"/api/v1/jobs/{job_a}/result", headers=_bearer(token_a)
        )
        own_bundle_b = user_b.get(
            f"/api/v1/jobs/{job_b}/result", headers=_bearer(token_b)
        )
        assert own_bundle_a.status_code == own_bundle_b.status_code == 200
        assert own_bundle_a.content.startswith(b"PK")
        assert own_bundle_b.content.startswith(b"PK")

        with app.state.store.connect() as connection:
            rows = connection.execute(
                "SELECT job_id, token_digest FROM jobs ORDER BY job_id"
            ).fetchall()
        assert {row["job_id"] for row in rows} == {job_a, job_b}
        assert all(token_a != row["token_digest"] for row in rows)
        assert all(token_b != row["token_digest"] for row in rows)
