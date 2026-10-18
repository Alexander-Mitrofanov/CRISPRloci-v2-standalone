from __future__ import annotations

from dataclasses import replace

import json
from pathlib import Path
import os
import signal
import subprocess
import time

import pytest

import crispr_web.worker as worker_module
from crispr_web.models import JobSubmission
from crispr_web.service import JobService
from crispr_web.worker import ProcessContainmentError, StorageCleanupError, Worker


def _run_case(settings, marker: str, mode: str = "detection", *, record_count: int = 1):
    service = JobService(settings)
    sequence = (
        f">{marker}\nACGT\n"
        if record_count == 1
        else "".join(
            f">{marker}_{index}\nACGT\n" for index in range(1, record_count + 1)
        )
    )
    _, queued = service.submit(
        JobSubmission(sequence=sequence, mode=mode),
        f"198.51.100.{100 + sum(marker.encode()) % 100}",
    )
    assert Worker(
        settings, store=service.store, worker_id=f"{marker}-worker"
    ).run_once()
    return service, queued, service.get(queued["job_id"])


@pytest.mark.parametrize(
    ("marker", "mode"),
    [
        ("invalidadapter", "detection"),
        ("stalesource", "detection"),
        ("invalidarraycount", "detection"),
        ("headeronly", "reconstruction"),
        ("headeronly", "orientation"),
        ("invalidtree", "reconstruction"),
        ("invalidtree", "orientation"),
        ("invalidgroups", "reconstruction"),
        ("invalidgroups", "orientation"),
        ("invalidorientation", "orientation"),
    ],
)
def test_malformed_cross_file_scientific_outputs_fail_closed(
    settings_factory, marker, mode
):
    _service, _queued, job = _run_case(settings_factory(), marker, mode)
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"
    assert job["summary"] is None
    assert job["artifacts"] == []


@pytest.mark.parametrize(
    "marker",
    [
        "numericblank",
        "numerictext",
        "numericnan",
        "numericinf",
        "negativerate",
        "negativetree",
        "negativeruntime",
        "negativecount",
        "fractionalcount",
        "nanlikelihood",
        "invalidmodel",
        "negativeinsertions",
        "fractionalinsertions",
        "treelengthmismatch",
    ],
)
def test_spacerplacer_numeric_results_are_finite_and_domain_safe(
    settings_factory, marker
):
    _service, _queued, job = _run_case(settings_factory(), marker, "reconstruction")
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"
    assert job["summary"] is None
    assert job["artifacts"] == []


@pytest.mark.parametrize(
    "marker",
    [
        "nanlikelihood",
        "invalidmodel",
        "negativeinsertions",
        "fractionalinsertions",
        "treelengthmismatch",
    ],
)
def test_crispr_evor_reconstruction_statistics_fail_closed(settings_factory, marker):
    _service, _queued, job = _run_case(
        settings_factory(), marker, "orientation", record_count=2
    )
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"
    assert job["summary"] is None
    assert job["artifacts"] == []


@pytest.mark.parametrize(
    ("marker", "mode"),
    [
        ("treewrong", "reconstruction"),
        ("treemissing", "reconstruction"),
        ("treeextra", "reconstruction"),
        ("treeduplicate", "reconstruction"),
        ("treemissingbranch", "reconstruction"),
        ("treewrong", "orientation"),
        ("treemissing", "orientation"),
        ("treeextra", "orientation"),
        ("treeduplicate", "orientation"),
        ("treemissingbranch", "orientation"),
        ("forwardbad", "orientation"),
        ("reversebad", "orientation"),
        ("selectedbad", "orientation"),
        ("wrongselectedpolicy", "orientation"),
        ("forwardprediction_wrongselectedpolicy", "orientation"),
        ("ndprediction_wrongselectedpolicy", "orientation"),
    ],
)
def test_downstream_trees_bind_exact_unique_adapter_leaf_identities(
    settings_factory, marker, mode
):
    _service, _queued, job = _run_case(settings_factory(), marker, mode, record_count=2)
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"
    assert job["summary"] is None
    assert job["artifacts"] == []


@pytest.mark.parametrize("mode", ["reconstruction", "orientation"])
def test_tree_length_rounding_tolerance_accepts_valid_serialization(
    settings_factory, mode
):
    _service, _queued, job = _run_case(
        settings_factory(), "roundedtreelength", mode, record_count=2
    )
    assert job["status"] == "completed"


@pytest.mark.parametrize(
    ("mode", "section"),
    [("reconstruction", "reconstruction"), ("orientation", "orientation")],
)
def test_completed_summary_exposes_validated_selected_tree(
    settings_factory, mode, section
):
    _service, _queued, job = _run_case(
        settings_factory(), "treegraphic", mode, record_count=2
    )
    assert job["status"] == "completed"
    tree = job["summary"][section]["trees"][0]
    assert tree["group"] == "group_1"
    assert tree["selected_newick"].endswith(";")
    assert "treegraphic_1" in tree["selected_newick"]
    assert "treegraphic_2" in tree["selected_newick"]
    assert job["summary"][section]["trees_truncated"] is False
    if mode == "orientation":
        assert set(tree) == {
            "group",
            "forward_newick",
            "reverse_newick",
            "selected_newick",
        }


@pytest.mark.parametrize(
    ("marker", "expected_prediction", "expected_recommendation"),
    [
        ("forwardprediction", "Forward", False),
        ("ndprediction", "ND", None),
    ],
)
def test_supported_orientation_states_select_the_contract_tree(
    settings_factory, marker, expected_prediction, expected_recommendation
):
    _service, _queued, job = _run_case(
        settings_factory(), marker, "orientation", record_count=2
    )
    assert job["status"] == "completed"
    comparison = job["summary"]["orientation"]["comparisons"][0]
    assert comparison["prediction"] == expected_prediction
    assert comparison["recommended_reverse"] is expected_recommendation


@pytest.mark.parametrize("marker", ["tampernormalized", "tamperrecord"])
def test_subprocess_input_mutation_is_detected(settings_factory, marker):
    _service, _queued, job = _run_case(settings_factory(), marker)
    assert job["status"] == "failed"
    assert job["error"]["code"] == "input_integrity_failed"
    assert job["summary"] is None
    assert job["artifacts"] == []


def test_containment_failure_marks_current_job_and_disables_next_claim(
    settings_factory, monkeypatch
):
    settings = settings_factory()
    service = JobService(settings)
    _, first = service.submit(
        JobSubmission(sequence=">first\nACGT\n", mode="detection"),
        "198.51.100.211",
    )
    _, second = service.submit(
        JobSubmission(sequence=">second\nACGT\n", mode="detection"),
        "198.51.100.212",
    )
    worker = Worker(settings, store=service.store, worker_id="containment-worker")

    def fail_containment(*_args, **_kwargs):
        raise ProcessContainmentError("simulated residual process")

    monkeypatch.setattr(worker_module, "_terminate_tree", fail_containment)
    with pytest.raises(ProcessContainmentError):
        worker.run_once()
    first_job = service.get(first["job_id"])
    assert first_job["status"] == "failed"
    assert first_job["error"]["code"] == "process_containment_failed"
    assert service.get(second["job_id"])["status"] == "queued"

    with pytest.raises(ProcessContainmentError):
        worker.run_once()
    assert service.get(second["job_id"])["status"] == "queued"


def test_terminate_tree_detects_a_process_that_survives_signals(monkeypatch):
    baseline_pids = set(worker_module._descendants(os.getpid()))
    process = subprocess.Popen(
        ["/bin/sh", "-c", "exec sleep 30"],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
        close_fds=True,
    )
    monkeypatch.setattr(worker_module, "_signal_groups", lambda *_args: None)
    monkeypatch.setattr(worker_module, "_signal_pids", lambda *_args: None)
    try:
        with pytest.raises(ProcessContainmentError):
            worker_module._terminate_tree(
                process,
                {process.pid},
                baseline_pids,
                grace_seconds=0,
                kill_seconds=0,
            )
    finally:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=3)


@pytest.mark.parametrize(
    ("marker", "pid_name"),
    [("detach", "detached.pid"), ("doublefork", "doublefork.pid")],
)
def test_successful_job_reaps_detached_and_double_forked_children(
    settings_factory, marker, pid_name
):
    service, queued, job = _run_case(settings_factory(), marker)
    assert job["status"] == "completed"
    pid_path = service.job_root(queued["job_id"]) / "results" / "attempt-001" / pid_name
    child_pid = int(pid_path.read_text(encoding="ascii"))
    deadline = time.monotonic() + 2
    while Path(f"/proc/{child_pid}").exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert not Path(f"/proc/{child_pid}").exists()


def test_launcher_log_flood_purges_only_current_attempt(settings_factory):
    settings = settings_factory(max_launcher_log_bytes=4_096)
    service, queued, job = _run_case(settings, "logflood")
    assert job["status"] == "failed"
    assert job["error"]["code"] == "job_storage_limit"
    job_root = service.job_root(queued["job_id"])
    assert (job_root / "input" / "normalized.fasta").is_file()
    assert not (job_root / "results" / "attempt-001").exists()
    assert not (job_root / "tmp" / "attempt-001").exists()
    assert not (job_root / "exports" / "attempt-001").exists()
    assert not (job_root / "worker-logs" / "attempt-001.stdout.log").exists()
    assert not (job_root / "worker-logs" / "attempt-001.stderr.log").exists()
    assert not (job_root / ".storage-quarantine-attempt-001").exists()


def test_sparse_disk_flood_is_detected_from_apparent_size(settings_factory):
    settings = settings_factory(
        max_job_disk_bytes=10_000_000,
        max_archive_bytes=1_000_000,
        max_artifact_file_bytes=1_000_000,
        max_launcher_log_bytes=1_000_000,
        job_disk_check_seconds=1,
    )
    service, queued, job = _run_case(settings, "diskflood")
    assert job["status"] == "failed"
    assert job["error"]["code"] == "job_storage_limit"
    assert job["artifacts"] == []
    job_root = service.job_root(queued["job_id"])
    assert not (job_root / "results" / "attempt-001").exists()
    assert not (job_root / "tmp" / "attempt-001").exists()
    assert not (job_root / "worker-logs" / "attempt-001.stdout.log").exists()
    assert not (job_root / "worker-logs" / "attempt-001.stderr.log").exists()


def test_attempt_purge_failure_marks_job_and_disables_next_claim(
    settings_factory, monkeypatch
):
    settings = settings_factory()
    service = JobService(settings)
    _, first = service.submit(
        JobSubmission(sequence=">first-storage\nACGT\n", mode="detection"),
        "198.51.100.213",
    )
    _, second = service.submit(
        JobSubmission(sequence=">second-storage\nACGT\n", mode="detection"),
        "198.51.100.214",
    )
    worker = Worker(settings, store=service.store, worker_id="storage-cleanup-worker")
    monkeypatch.setattr(worker_module, "_job_tree_within_budget", lambda *_args: False)

    original_rmtree = worker_module.shutil.rmtree

    def fail_quarantine_removal(path, *args, **kwargs):
        if Path(path).name.startswith(".storage-quarantine-attempt-"):
            raise OSError("simulated quarantine removal failure")
        return original_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(worker_module.shutil, "rmtree", fail_quarantine_removal)
    with pytest.raises(StorageCleanupError):
        worker.run_once()

    first_job = service.get(first["job_id"])
    assert first_job["status"] == "failed"
    assert first_job["error"]["code"] == "storage_cleanup_failed"
    quarantine = service.job_root(first["job_id"]) / ".storage-quarantine-attempt-001"
    assert (quarantine / "temporary").is_dir()
    assert service.get(second["job_id"])["status"] == "queued"

    with pytest.raises(StorageCleanupError):
        worker.run_once()
    assert service.get(second["job_id"])["status"] == "queued"


def test_scientific_tool_temporary_paths_are_per_job(settings_factory):
    service, queued, job = _run_case(settings_factory(), "envcheck")
    assert job["status"] == "completed"
    job_root = service.job_root(queued["job_id"])
    environment = json.loads(
        (job_root / "results" / "attempt-001" / "environment.json").read_text(
            encoding="utf-8"
        )
    )
    temporary_root = job_root / "tmp" / "attempt-001"
    assert environment == {
        "TMPDIR": str(temporary_root),
        "TMP": str(temporary_root),
        "TEMP": str(temporary_root),
        "XDG_CACHE_HOME": str(temporary_root / "cache"),
    }


def test_additive_supported_major_schema_versions_are_accepted(settings_factory):
    _service, _queued, job = _run_case(settings_factory(), "compatibleminor")
    assert job["status"] == "completed"


@pytest.mark.parametrize("marker", ["unknownschema", "malformedschema"])
def test_unknown_major_and_malformed_schema_versions_fail_closed(
    settings_factory, marker
):
    _service, _queued, job = _run_case(settings_factory(), marker)
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"


def test_worker_refuses_to_claim_when_scientific_manifest_changes(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">manifesttamper\nACGT\n", mode="detection"),
        "198.51.100.198",
    )
    settings.scientific_release_manifest.write_text("{}\n", encoding="utf-8")
    worker = Worker(settings, store=service.store, worker_id="manifest-worker")

    with pytest.raises(ValueError, match="SHA-256"):
        worker.run_once()

    stored = service.store.get_job(queued["job_id"])
    assert stored is not None
    assert stored["status"] == "queued"
    assert stored["attempt"] == 0


def test_worker_refuses_launch_after_persisted_total_budget(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">budgetdeadline\nACGT\n", mode="detection"),
        "198.51.100.199",
    )
    worker = Worker(settings, store=service.store, worker_id="deadline-worker")
    claim = service.store.claim_next("deadline-worker")
    assert claim is not None
    expired_claim = replace(
        claim, started_at=time.time() - settings.job_timeout_seconds - 1
    )

    worker._execute(expired_claim)

    job = service.get(queued["job_id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "job_timeout"
    attempt_root = service.job_root(queued["job_id"]) / "results" / "attempt-001"
    assert not attempt_root.exists()
