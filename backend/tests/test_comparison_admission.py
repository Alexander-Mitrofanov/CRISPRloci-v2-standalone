"""Observed-array maps must fit the durable queue's admission envelope."""
from dataclasses import replace
from pathlib import Path

import pytest

from crispr_web.fasta import FastaValidationError
from crispr_web.models import JobSubmission
from crispr_web.service import JobService
from .test_array_compare import comparison_settings


def packet_text():
    return (Path(__file__).resolve().parents[2] / "deploy/fixtures/verification-array-sharing.json").read_text()


def test_inflight_comparison_map_is_included_in_inode_reservation(settings_factory, tmp_path, monkeypatch):
    settings = replace(comparison_settings(settings_factory, tmp_path), max_records=6)
    service = JobService(settings)
    create = service.store.create_job
    observed = []

    def measure_before_commit(**kwargs):
        root = service.job_root(kwargs["job_id"])
        assert (root / "array-map.json").is_file()
        count = 1 + sum(1 for _ in root.rglob("*"))
        observed.append(count)
        assert settings.submission_admission_free_inodes >= settings.min_free_inodes + count * settings.max_active_jobs
        return create(**kwargs)

    monkeypatch.setattr(service.store, "create_job", measure_before_commit)
    service.submit(JobSubmission(mode="array_compare", sequence=packet_text()), "192.0.2.24")
    assert observed == [12]  # Six FASTAs, three directories, canonical FASTA, map, lock.


def test_comparison_observations_honor_existing_operator_record_limit(settings_factory, tmp_path):
    settings = replace(comparison_settings(settings_factory, tmp_path), max_records=5)
    service = JobService(settings)
    with pytest.raises(FastaValidationError, match="5-record limit"):
        service.submit(JobSubmission(mode="array_compare", sequence=packet_text()), "192.0.2.25")
    assert list(settings.jobs_root.iterdir()) == []
