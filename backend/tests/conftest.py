from __future__ import annotations

from pathlib import Path
import hashlib
import json
import sys

import pytest

from crispr_web.config import Settings


def write_scientific_release_manifest(path: Path) -> tuple[str, dict[str, object]]:
    artifacts = {
        "crispridentify": {
            "name": "crispridentify",
            "display_name": "CRISPRidentify",
            "version": "2.0.0",
            "kind": "source_tree",
            "path": "/opt/crispr-suite/CRISPRidentify",
            "sha256": "3" * 64,
            "sha256_scope": "source_tree_v1",
        },
        "spacerplacer": {
            "name": "spacerplacer",
            "display_name": "SpacerPlacer",
            "version": "1.1.0",
            "kind": "wheel",
            "path": "/opt/crispr-suite/SpacerPlacer.whl",
            "sha256": "4" * 64,
            "sha256_scope": "file_bytes",
        },
        "crispr_evor": {
            "name": "crispr_evor",
            "display_name": "CRISPR-evOr",
            "version": "1.0.1.post1",
            "kind": "wheel",
            "path": "/opt/crispr-suite/CRISPR-evOr.whl",
            "sha256": "5" * 64,
            "sha256_scope": "file_bytes",
        },
        "integration": {
            "name": "integration",
            "display_name": "CRISPR tool integration",
            "version": "0.2.5",
            "kind": "wheel",
            "path": "/opt/crispr-suite/Integration.whl",
            "sha256": "6" * 64,
            "sha256_scope": "file_bytes",
        },
    }
    runtime_files = {
        key: {
            "path": f"/opt/crispr-suite/runtime/{key}",
            "sha256": "7" * 64,
            "sha256_scope": "file_bytes",
        }
        for key in (
            "runner_config",
            "run_crispridentify",
            "run_spacerplacer",
            "run_crispr_evor",
            "crispridentify_entrypoint",
            "spacerplacer_entrypoint",
            "crispr_evor_entrypoint",
            "integration_entrypoint",
        )
    }
    document: dict[str, object] = {
        "schema": "crispr-web-scientific-release",
        "schema_version": "1.0.0",
        "release_id": "test-release",
        "platform": {"os": "linux", "architecture": "x86_64", "python": "3.12"},
        "scientific_artifacts": artifacts,
        "runtime_files": runtime_files,
    }
    payload = (
        json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest(), document


@pytest.fixture
def scientific_manifest(tmp_path):
    path = tmp_path / "release-manifest.json"
    digest, document = write_scientific_release_manifest(path)
    return path, digest, document


@pytest.fixture
def settings_factory(tmp_path):
    def factory(**overrides):
        runner_config = tmp_path / "runners.json"
        runner_config.write_text("{}\n", encoding="utf-8")
        manifest_path = tmp_path / "release-manifest.json"
        manifest_digest, _document = write_scientific_release_manifest(manifest_path)
        fake = Path(__file__).with_name("fake_integration.py").resolve()
        values = {
            "data_root": tmp_path / "data",
            "database_path": tmp_path / "data" / "queue.sqlite3",
            "runner_config": runner_config,
            "integration_command": (sys.executable, str(fake)),
            "cors_origins": ("https://frontend.example",),
            "release_id": "1" * 64,
            "backend_wheel_sha256": "2" * 64,
            "scientific_release_manifest": manifest_path,
            "scientific_release_manifest_sha256": manifest_digest,
            "max_total_bases": 10_000,
            "max_record_bases": 8_000,
            "max_records": 10,
            "max_queued_jobs": 4,
            "rate_limit_submissions": 20,
            "rate_limit_window_seconds": 3_600,
            "stage_timeout_seconds": 30,
            "job_timeout_seconds": 60,
            "worker_lease_seconds": 10,
            "worker_stale_seconds": 20,
            "worker_poll_seconds": 1,
            "hsts": False,
            # Functional fixtures reserve resources for tiny synthetic jobs. Capacity
            # and containment tests override these values explicitly.
            "min_free_bytes": 1_000_000,
            "min_free_inodes": 100,
            "max_job_disk_bytes": 600_000_000,
            "max_job_entries": 5_000,
        }
        values.update(overrides)
        return Settings(**values)

    return factory
