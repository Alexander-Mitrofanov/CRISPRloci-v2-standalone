"""Scientific settings survive queuing and reach the isolated execution boundary."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import zipfile

import pytest
from pydantic import ValidationError

from crispr_web.db import Store
from crispr_web.models import JobSubmission
from crispr_web.scientific_options import TOOL_OPTIONS_ENV
from crispr_web.service import JobService
from crispr_web.worker import Worker


@pytest.mark.parametrize(("mode", "options"), [
    ("detection", {"crispridentify": {"--cpu": 999}}),
    ("detection", {"crispridentify": {"--file": "/tmp/private"}}),
    ("detection", {"spacerplacer": {"--seed": 1}}),
    ("tracrrna", {"crispr_tracrrna": {"--similarity": 0.8}}),
    ("repeats", {"crisprrepeat": {"--window": 100}}),
    ("detection", {"crispridentify": {"--min_repeats": True}}),
    ("detection", {"crispridentify": {"--fast_run": "False"}}),
    ("detection", {"crispridentify": {"--min_len_rep": 60}}),
    ("detection", {"crispridentify": {"--fast_run": False, "--fast_run_seed_profile": "legacy"}}),
    ("loci", {"crispr_tracrrna": {"--model-rescue": False, "--model-rescue-evalue": 0.01}}),
    ("repeat_context", {"crisprrepeat": {"--window": 50}}),
    ("reconstruction", {"spacerplacer": {"--insertion_rate": float("nan")}}),
    ("reconstruction", {"spacerplacer": {"--insertion_rate": 10**400}}),
    ("reconstruction", {"spacerplacer": {"--seed": "1; touch /tmp/no"}}),
    ("repeat_map", {"crisprmap": {"--max-distance": 4}}),
    ("protospacer", {"crisprspacer": {"--substitutions": 1}}),
])
def test_rejects_unsupported_unsafe_inert_and_invalid_options(mode, options):
    with pytest.raises(ValidationError):
        JobSubmission(sequence=">a\nACGT\n", mode=mode, tool_options=options)


def test_options_survive_store_reopen_and_worker_environment(settings_factory, monkeypatch):
    import crispr_web.worker as worker_module

    options = {"crispridentify": {"--fast_run": False, "--min_repeats": 4}}
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(JobSubmission(sequence=">sample\nACGT\n", mode="detection", tool_options=options), "198.51.100.90")
    reopened = Store(settings)
    assert reopened.get_job(queued["job_id"])["tool_options"] == options
    claim = reopened.claim_next("options-worker")
    assert claim.tool_options == options
    reopened.release_job(claim.job_id, "options-worker")

    seen = []
    real_popen = worker_module.subprocess.Popen

    def capturing_popen(*args, **kwargs):
        seen.append(json.loads(kwargs["env"][TOOL_OPTIONS_ENV]))
        return real_popen(*args, **kwargs)

    # An inherited override must never contaminate a queued job's selected flags.
    monkeypatch.setenv(TOOL_OPTIONS_ENV, '{"crispridentify":{"--cpu":999}}')
    monkeypatch.setattr(worker_module.subprocess, "Popen", capturing_popen)
    assert Worker(settings, store=reopened, worker_id="options-worker").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] == "completed"
    assert seen == [options]
    assert job["options"]["tool_options"] == options
    bundle_path, _ = service.result_bundle(job["job_id"])
    with zipfile.ZipFile(bundle_path) as bundle:
        provenance = json.loads(bundle.read("provenance.json"))
    assert provenance["options"]["tool_options"] == options
    assert provenance["execution"]["logical_invocation"]["tool_options"] == options

    # A stale runtime or changed launcher cannot silently drop or replace flags.
    from crispr_web.execution.validation_pipeline import _validated_pipeline_status
    output = service.job_root(claim.job_id) / "results" / "attempt-002"
    marker = output / "integration_result.json"
    report = json.loads(marker.read_text())
    sources = {"sample": (4, hashlib.sha256(b"ACGT").hexdigest())}
    assert _validated_pipeline_status(claim, output, sources) == "completed"
    for reported in ({}, {"crispridentify": {"--min_repeats": 3}}):
        report["tool_options"] = reported
        marker.write_text(json.dumps(report))
        assert _validated_pipeline_status(claim, output, sources) is None
    report.pop("tool_options")
    marker.write_text(json.dumps(report))
    assert _validated_pipeline_status(claim, output, sources) is None


def test_pre_options_database_migrates_atomically(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(JobSubmission(sequence=">sample\nACGT\n", mode="detection"), "198.51.100.91")
    with sqlite3.connect(settings.database_path) as connection:
        connection.execute("ALTER TABLE jobs DROP COLUMN tool_options_json")
        connection.execute("UPDATE metadata SET value = '6' WHERE name = 'schema_version'")
    migrated = Store(settings)
    assert migrated.get_job(queued["job_id"])["tool_options"] == {}
    assert migrated.claim_next("migration-worker").tool_options == {}
