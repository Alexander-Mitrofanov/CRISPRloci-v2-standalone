from __future__ import annotations

import json
from pathlib import Path
import threading
import time
import zipfile

import pytest

from crispr_web.artifacts import resolve_registered_file
from crispr_web.models import JobSubmission
from crispr_web.service import JobService
from crispr_web.worker import Worker


def test_orientation_worker_summary_options_and_reproducible_bundle(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    token, queued = service.submit(
        JobSubmission(
            sequence=">isolate_A\nACGTRYSWKMBDHVN\n>isolate_B\nACGTACGT\n",
            filename="../../cohort.fa",
            mode="orientation",
            spacer_distance=2,
        ),
        "198.51.100.1",
    )
    assert Worker(settings, store=service.store, worker_id="test-worker").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] == "completed"
    assert job["summary"]["adapter"]["unknown_strand_excluded_count"] == 1
    assert (
        job["summary"]["detection"]["arrays"][0]["model_score_is_probability"] is False
    )
    assert job["summary"]["detection"]["score_semantics"] == {
        "primary_interpretation": "category",
        "model_score_is_probability": False,
        "note": "Raw model scores are not calibrated probabilities.",
    }
    assert (
        job["summary"]["orientation"]["selected_reconstructions"][0]["name"]
        == "group_1"
    )
    assert job["options"]["spacer_distance"] == 2

    names = {item["name"] for item in job["artifacts"]}
    assert {"results.zip", "provenance.json", "adapter/manifest.json"} <= names
    assert "integration_result.json" not in names
    assert all(not name.endswith((".pkl", ".html")) for name in names)

    bundle_path, _ = service.result_bundle(job["job_id"])
    with zipfile.ZipFile(bundle_path) as bundle:
        assert all(
            info.file_size > 0
            for info in bundle.infolist()
            if not info.is_dir()
        )
        bundled = set(bundle.namelist())
        assert "input/normalized.fasta" in bundled
        assert "provenance.json" in bundled
        assert "results/adapter/manifest.json" in bundled
        assert "results/integration_result.json" not in bundled
        assert "results/crispr_evor/detailed_results/events.csv" in bundled
        assert "results/logs/integration.log" in bundled
        assert not any(name.startswith("logs/attempt-001") for name in bundled)
        assert not any(name.endswith((".pkl", ".html")) for name in bundled)
        assert "results/crispr_evor/external.log" not in bundled
        provenance = json.loads(bundle.read("provenance.json"))
        manifest = json.loads(bundle.read("results/adapter/manifest.json"))

    assert provenance["input"]["normalized_sha256"] == job["input"]["sha256"]
    from crispr_web import __version__
    assert provenance["software"]["backend_version"] == __version__
    assert provenance["software"]["backend_release_id"] == settings.release_id
    assert (
        provenance["software"]["backend_wheel_sha256"] == settings.backend_wheel_sha256
    )
    assert (
        provenance["software"]["scientific_release_manifest_sha256"]
        == settings.scientific_release_manifest_sha256
    )
    assert (
        provenance["software"]["scientific_artifacts"]
        == settings.scientific_release_identities
    )
    assert all(
        "path" not in identity
        for identity in provenance["software"]["scientific_artifacts"].values()
    )
    invocation = provenance["execution"]["logical_invocation"]
    assert invocation["analysis"] == "crispr-evor"
    assert invocation["spacer_distance"] == 2
    assert invocation["bias_corrections"] is True
    assert invocation["identify_folder_jobs"] == 1
    assert invocation["identify_folder_jobs_configured"] == 1
    assert "launcher_command" not in provenance["execution"]
    assert all(
        "command" not in stage for stage in provenance["execution"]["stages"].values()
    )
    assert {citation["doi"] for citation in provenance["citations"]} == {
        "10.1093/nar/gkaa1158",
        "10.1093/nar/gkae772",
        "10.1371/journal.pcbi.1013706",
    }
    assert manifest["sanitization"]["raw_internal_manifest_published"] is False
    assert manifest["sanitization"]["private_fields_removed"] >= 3
    assert any(
        "Worker log omitted because it was empty" in warning
        for warning in job["summary"]["warnings"]
    )
    published = json.dumps({"provenance": provenance, "manifest": manifest})
    assert token not in published
    for private_root in ("/srv/", "/opt/", "/etc/", str(settings.data_root)):
        assert private_root not in published


def test_mode_mapping_and_operator_bias_override_are_explicit(settings_factory):
    settings = settings_factory(
        disable_bias_corrections=True,
        worker_cpu=4,
        identify_folder_jobs=2,
    )
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(
            sequence=">detect\nACGT\n",
            mode="detection",
            category_policy="bona_fide_only",
            bias_corrections=True,
        ),
        "198.51.100.2",
    )
    view = service.get(queued["job_id"])
    assert view["options"]["bias_corrections_requested"] is True
    assert view["options"]["bias_corrections_effective"] is False
    assert view["options"]["warnings"]

    worker = Worker(settings, store=service.store, worker_id="mode-worker")
    claim = service.store.claim_next("mode-worker")
    assert claim is not None
    command = worker.command_for(claim, Path("/tmp/output"))
    assert command[command.index("--tools") + 1] == "identify"
    category_index = command.index("--categories")
    assert command[category_index + 1] == "Bona-fide"
    assert command[category_index + 2] == "--spacer-distance"
    assert "--disable-bias-corrections" in command
    assert command[command.index("--identify-folder-jobs") + 1] == "2"
    service.store.release_job(claim.job_id, "mode-worker")

    with pytest.raises(ValueError, match="cannot exceed worker_cpu"):
        settings_factory(worker_cpu=2, identify_folder_jobs=3)


def test_strict_category_policy_filters_primary_detection_list(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(
            sequence=">strict\nACGT\n",
            mode="detection",
            category_policy="bona_fide_only",
        ),
        "198.51.100.21",
    )
    Worker(settings, store=service.store, worker_id="strict-worker").run_once()
    detection = service.get(queued["job_id"])["summary"]["detection"]
    assert detection["array_count"] == 3
    assert detection["category_counts"] == {"Bona-fide": 1, "Possible": 2}
    assert detection["selected_array_count"] == 1
    assert detection["unselected_array_count"] == 2
    assert detection["selected_category_counts"] == {"Bona-fide": 1}
    assert detection["arrays_scope"] == "selected_for_analysis"
    assert {item["category"] for item in detection["arrays"]} == {"Bona-fide"}
    assert all(item["selected_for_analysis"] for item in detection["arrays"])


def test_no_eligible_groups_is_successful_scientific_outcome(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">noeligible\nACGT\n", mode="reconstruction"),
        "198.51.100.3",
    )
    Worker(settings, store=service.store, worker_id="no-result-worker").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] == "completed_no_eligible_groups"
    assert job["summary"]["adapter"]["emitted_group_count"] == 0


@pytest.mark.parametrize("mode", ["detection", "reconstruction", "orientation"])
def test_exit_zero_status_only_result_fails_mode_contract(settings_factory, mode):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">invalidcontract\nACGT\n", mode=mode),
        f"198.51.100.{30 + ('detection', 'reconstruction', 'orientation').index(mode)}",
    )
    Worker(settings, store=service.store, worker_id=f"invalid-{mode}").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"


def test_success_requires_one_identify_report_per_input_record(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">missingreport\nACGT\n", mode="detection"),
        "198.51.100.39",
    )
    Worker(settings, store=service.store, worker_id="missing-report-worker").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "invalid_pipeline_result"


def test_public_summary_recursively_redacts_tool_paths(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">pathleak\nACGT\n", mode="orientation"),
        "198.51.100.42",
    )
    Worker(settings, store=service.store, worker_id="path-redaction-worker").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] == "completed"
    rendered = json.dumps(job["summary"])
    for marker in ("/srv", "/opt", "/etc", "/home", "/tmp"):
        assert marker not in rendered
    assert any("summary value" in warning for warning in job["summary"]["warnings"])


def test_artifact_integrity_traversal_and_symlink_checks(settings_factory, tmp_path):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">integrity\nACGT\n", mode="reconstruction"),
        "198.51.100.40",
    )
    Worker(settings, store=service.store, worker_id="integrity-worker").run_once()
    job = service.get(queued["job_id"])
    artifact = next(
        item
        for item in job["artifacts"]
        if item["name"] == "spacerplacer/0_results.csv"
    )
    path, _ = service.artifact_file(job["job_id"], artifact["artifact_id"])
    original = path.read_bytes()
    path.write_bytes(bytes([original[0] ^ 1]) + original[1:])
    with pytest.raises(FileNotFoundError):
        service.artifact_file(job["job_id"], artifact["artifact_id"])

    root = service.job_root(job["job_id"])
    with pytest.raises(FileNotFoundError):
        resolve_registered_file(root, "../outside")
    target = tmp_path / "outside.txt"
    target.write_text("outside", encoding="utf-8")
    link = root / "linked-result.txt"
    try:
        link.symlink_to(target)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")
    with pytest.raises(FileNotFoundError):
        resolve_registered_file(root, "linked-result.txt")


def test_transient_store_error_terminates_live_process_tree(
    settings_factory, monkeypatch
):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">slow\nACGT\n", mode="detection"),
        "198.51.100.41",
    )
    pid_path = (
        service.job_root(queued["job_id"]) / "results" / "attempt-001" / "runner.pid"
    )
    observed_pids: list[int] = []
    original = service.store.cancellation_requested

    def fail_after_process_started(job_id: str, worker_id: str) -> bool:
        if pid_path.is_file():
            observed_pids.append(int(pid_path.read_text(encoding="ascii")))
            raise RuntimeError("transient sqlite read failure")
        return original(job_id, worker_id)

    monkeypatch.setattr(
        service.store, "cancellation_requested", fail_after_process_started
    )
    worker = Worker(settings, store=service.store, worker_id="fault-worker")
    thread = threading.Thread(target=worker.run_once)
    thread.start()
    thread.join(timeout=8)
    assert not thread.is_alive()
    assert observed_pids
    assert not Path(f"/proc/{observed_pids[0]}").exists()
    job = service.get(queued["job_id"])
    assert job["status"] == "failed"
    assert job["error"]["code"] == "worker_error"


def test_running_job_cancellation_terminates_worker_process(settings_factory):
    settings = settings_factory()
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=">slow\nACGT\n", mode="detection"),
        "198.51.100.4",
    )
    worker = Worker(settings, store=service.store, worker_id="cancel-worker")
    thread = threading.Thread(target=worker.run_once)
    thread.start()
    deadline = time.monotonic() + 5
    while service.get(queued["job_id"])["status"] != "running":
        assert time.monotonic() < deadline
        time.sleep(0.02)
    service.cancel(queued["job_id"])
    thread.join(timeout=8)
    assert not thread.is_alive()
    assert service.get(queued["job_id"])["status"] == "cancelled"
