"""Real isolated CRISPRmap pipeline plus independent tamper/admission gates."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import pytest
from crispr_integration.runner import RunnerConfig
from crispr_integration.repeat_map import run_repeat_map_pipeline
from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.publication.repeat_map import write_repeat_map_exports
from crispr_web.summaries.repeat_map import repeat_map_summary
from crispr_web.models import JobSubmission
from crispr_web.service import JobService, ReferenceUnavailableError
from crispr_web.fasta import FastaValidationError
from crispr_web.worker import Worker
from .test_viral_search import configured_settings

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tools/CRISPRmap-v2"
NATIVE = (sys.executable, "-c", f"import sys;sys.path.insert(0,{str(SOURCE / 'src')!r});from crisprmap.cli import main;sys.exit(main())")


def sha(b):
    return hashlib.sha256(b).hexdigest()


@pytest.fixture
def map_reference(tmp_path):
    reference = tmp_path / "reference.json"
    doc = json.loads((SOURCE / "examples/reference.json").read_text())
    # Undisplayed ties must still determine conflicts and missing annotations.
    first = doc["entries"][0]
    doc["entries"] = [{**first, "id": f"ref_{i:02}", "family": "CONFLICT" if i == 24 else first["family"],
                       "motif": None if i == 24 else first["motif"]} for i in range(25)] + [doc["entries"][1]]
    doc.pop("content_sha256")
    doc["content_sha256"] = sha((json.dumps(doc, sort_keys=True, separators=(",", ":")) + "\n").encode())
    reference.write_text(json.dumps(doc))
    return reference


@pytest.fixture
def map_pipeline(tmp_path, map_reference):
    source = tmp_path / "input"
    source.mkdir()
    records = {"exact": "ACGTCGATGCTAGCTAGTCGATCG", "near": "ACGTCGATGCTAGCTAGTCGATCA", "negative": "A" * 24,
               "ambiguous": "ACGTNNNN", "reverse": "CGATCGACTAGCTAGCATCGACGT"}
    for i, (key, sequence) in enumerate(records.items()):
        (source / f"{i:06d}.fasta").write_text(f">{key}\n{sequence}\n")
    digest = sha(map_reference.read_bytes())
    output = tmp_path / "output"
    run_repeat_map_pipeline(source, output, runner_config=RunnerConfig(("unused",), ("unused",), ("unused",), crisprmap=NATIVE),
                            reference_path=map_reference, reference_sha256=digest, timeout_seconds=60)
    claim = ClaimedJob("a" * 32, "repeat_map", 1, 0, "bona_fide_possible", 1, True, True, (), repeat_reference_sha256=digest)
    bindings = {k: (len(v), sha(v.encode())) for k, v in records.items()}
    return claim, output, bindings, source, map_reference


def validated(data):
    claim, output, bindings, source, reference = data
    return _validated_pipeline_status(claim, output, bindings, input_root=source, repeat_reference=reference)


def test_real_map_pipeline_and_complete_exports(map_pipeline):
    assert validated(map_pipeline) == "completed"
    out = map_pipeline[1]
    preview = repeat_map_summary(out, max_rows=1)
    assert preview["results_truncated"]
    row = preview["results"][0]
    assert row["best_hit_count"] == 25 and row["hits_truncated"]
    assert row["annotations"]["family"]["status"] == "conflicting"
    assert row["annotations"]["motif"]["status"] == "incomplete"
    write_repeat_map_exports(out, {"repeat_map": preview})
    public = json.loads((out / "repeat-mapping.json").read_text())
    assert len(public["results"]) == 5
    assert public["counts"] == {"query_count": 5, "unique_sequence_count": 5, "exact_count": 2, "near_count": 1, "no_hit_count": 1, "unsupported_count": 1}
    assert "/tmp/" not in json.dumps(public)
    assert len((out / "repeat-mapping.tsv").read_text().splitlines()) == 6


@pytest.mark.parametrize("mutation", ["source", "reference", "omitted_tie", "assignment", "alignment", "strand", "boolean", "status", "policy", "count", "stage", "symlink", "duplicate_key"])
def test_tampering_rejected(map_pipeline, mutation):
    claim, out, bindings, source, reference = map_pipeline
    path = out / "repeat-map/report.json"
    doc = json.loads(path.read_text())
    row = doc["results"][0]
    if mutation == "source":
        next(source.iterdir()).write_text(">forged\nACGT\n")
    elif mutation == "reference":
        reference.write_text(reference.read_text() + " ")
    elif mutation == "omitted_tie":
        row["best_hit_count"] -= 1
    elif mutation == "assignment":
        row["annotations"]["family"]["value"] = "DEMO-F1"
    elif mutation == "alignment":
        row["hits"][0]["reference_aligned"] = "T" * 24
    elif mutation == "strand":
        row["hits"][0]["orientation"] = "reverse"
    elif mutation == "boolean":
        doc["results"][1]["hits"][0]["distance"] = True
    elif mutation == "status":
        doc["results"][2]["status"] = "exact"
    elif mutation == "policy":
        doc["policy"]["max_distance"] = 4
    elif mutation == "count":
        doc["counts"]["near_count"] += 1
    elif mutation == "stage":
        p = out / "integration_result.json"
        stage = json.loads(p.read_text())
        stage["stages"]["crisprmap"]["reference_sha256"] = "a" * 64
        p.write_text(json.dumps(stage))
    path.write_text(json.dumps(doc))
    if mutation == "symlink":
        other = out / "outside.json"
        path.rename(other)
        path.symlink_to(other)
    elif mutation == "duplicate_key":
        path.write_text(path.read_text().replace('"execution_completed": true', '"execution_completed": false, "execution_completed": true'))
    assert validated(map_pipeline) is None


def map_settings(settings_factory, tmp_path, reference):
    settings = configured_settings(settings_factory, tmp_path)
    manifest = settings.scientific_release_manifest
    doc = json.loads(manifest.read_text())
    doc["schema_version"] = "6.0.0"
    doc["scientific_artifacts"]["crisprmap"] = {**doc["scientific_artifacts"]["crisprspacer"], "name": "crisprmap", "display_name": "CRISPRmap-v2", "version": "2.0.0a1"}
    manifest.write_text(json.dumps(doc))
    return replace(settings, scientific_release_manifest_sha256=sha(manifest.read_bytes()), repeat_reference_path=reference, repeat_reference_sha256=sha(reference.read_bytes()))


def test_admission_queue_and_command(settings_factory, tmp_path, map_reference):
    with pytest.raises(ReferenceUnavailableError):
        JobService(settings_factory()).submit(JobSubmission(mode="repeat_map", sequence=">r\nACGT\n"), "127.0.0.1")
    settings = map_settings(settings_factory, tmp_path, map_reference)
    service = JobService(settings)
    for sequence in (">r\n" + "A" * 201, "".join(f">r{i}\nACGT\n" for i in range(101))):
        with pytest.raises(FastaValidationError):
            service.submit(JobSubmission(mode="repeat_map", sequence=sequence), "127.0.0.1")
    _, job = service.submit(JobSubmission(mode="repeat_map", sequence=">r\nACGT\n"), "127.0.0.1")
    claim = service.store.claim_next("worker")
    assert claim.repeat_reference_sha256 == settings.repeat_reference_sha256
    assert job["mode"] == "repeat_map"
    worker = Worker(settings)
    command = worker.command_for(claim, tmp_path / "out")
    assert command[command.index("--repeat-reference-sha256") + 1] == settings.repeat_reference_sha256


def test_http_worker_download_recovery_and_reference_drift(settings_factory, tmp_path, map_reference):
    import io
    import zipfile
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app

    settings = map_settings(settings_factory, tmp_path, map_reference)
    settings = replace(settings, integration_command=(sys.executable, "-m", "crispr_integration"))
    settings.runner_config.write_text(json.dumps({"crispridentify": ["unused"], "spacerplacer": ["unused"],
                                                "crispr_evor": ["unused"], "crisprmap": list(NATIVE)}))
    worker = Worker(settings, worker_id="map-http-test")
    with TestClient(create_app(settings)) as client:
        assert "repeat_map" in client.get("/api/v1/config").json()["modes"]
        created = client.post("/api/v1/jobs", json={"mode": "repeat_map", "sequence": ">exact\nACGTCGATGCTAGCTAGTCGATCG\n"})
        assert created.status_code == 202, created.text
        data = created.json()
        job_id = data["job"]["job_id"]
        headers = {"Authorization": "Bearer " + data["access_token"]}
        assert worker.run_once()
        fetched = client.get(f"/api/v1/jobs/{job_id}", headers=headers)
        job = fetched.json()
        assert job["status"] == "completed", job
        assert job["summary"]["repeat_map"]["counts"]["exact_count"] == 1
        assert client.get(f"/api/v1/jobs/{job_id}").status_code in {401, 404}
        exports = {a["name"]: a for a in job["artifacts"]}
        for name in ("repeat-mapping.json", "repeat-mapping.tsv", "results.zip", "provenance.json"):
            response = client.get(exports[name]["download_url"], headers=headers)
            assert response.status_code == 200
            if name == "results.zip":
                with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
                    assert "results/repeat-mapping.json" in archive.namelist()
                    assert not any(n.startswith(("results/repeat-map/", "results/repeat-map-input/")) for n in archive.namelist())
                    assert json.loads(archive.read("provenance.json"))["options"]["repeat_reference_sha256"] == settings.repeat_reference_sha256
        # A queued job cannot silently move to a new reference.
        second = client.post("/api/v1/jobs", json={"mode": "repeat_map", "sequence": ">q\nACGT\n"}).json()
        map_reference.write_text(map_reference.read_text() + " ")
        assert worker.run_once()
        failed = client.get(f"/api/v1/jobs/{second['job']['job_id']}", headers={"Authorization": "Bearer " + second["access_token"]}).json()
        assert failed["status"] == "failed"
        assert not failed["artifacts"]


def test_detection_json_to_repeat_mapping_worker(settings_factory, tmp_path, map_reference):
    """Full published array port -> HTTP queue -> isolated native reference search."""
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app
    from .test_sequence_exchange import detection_ports

    packet = detection_ports(tmp_path)["repeats.json"]
    settings = replace(map_settings(settings_factory, tmp_path, map_reference),
                       integration_command=(sys.executable, "-m", "crispr_integration"))
    settings.runner_config.write_text(json.dumps({
        "crispridentify": ["unused"], "spacerplacer": ["unused"],
        "crispr_evor": ["unused"], "crisprmap": list(NATIVE),
    }))
    with TestClient(create_app(settings)) as client:
        response = client.post("/api/v1/jobs", json={"mode": "repeat_map", "sequence": json.dumps(packet)})
        assert response.status_code == 202, response.text
        submitted = response.json()
        assert Worker(settings, worker_id="handoff-map").run_once()
        headers = {"Authorization": "Bearer " + submitted["access_token"]}
        job = client.get(f"/api/v1/jobs/{submitted['job']['job_id']}", headers=headers).json()
        assert job["status"] == "completed", job
        assert job["summary"]["repeat_map"]["counts"]["query_count"] == 3
        files = {a["name"]: a for a in job["artifacts"]}
        assert client.get(files["input-sequences.json"]["download_url"]).status_code in {401, 404}
        exported = client.get(files["input-sequences.json"]["download_url"], headers=headers).json()
        assert exported == packet
        assert client.get(files["result.json"]["download_url"]).status_code in {401, 404}
        compact = client.get(files["result.json"]["download_url"], headers=headers).json()
        assert compact["complete"] is True and compact["module"] == "repeat_map"
        assert len(compact["data"]["results"]) == 3
        provenance = client.get(files["provenance.json"]["download_url"], headers=headers).json()
        assert provenance["input"]["handoff"]["submitted_sha256"] == sha(json.dumps(packet).encode())
