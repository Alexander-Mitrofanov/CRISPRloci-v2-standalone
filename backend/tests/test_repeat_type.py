"""Real isolated RepeatTyper execution and independent hostile-output validation."""

from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import pytest
from crispr_integration.runner import RunnerConfig
from crispr_integration.repeat_type import run_repeat_type_pipeline, validate_report
from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.publication.repeat_type import write_repeat_type_exports
from crispr_web.summaries.repeat_type import repeat_type_summary
from crispr_web.models import JobSubmission
from crispr_web.service import JobService, ReferenceUnavailableError
from crispr_web.fasta import FastaValidationError
from crispr_web.worker import Worker
from .test_viral_search import configured_settings

ROOT = Path(__file__).resolve().parents[2]


def sha(value):
    return hashlib.sha256(
        value.encode() if isinstance(value, str) else value
    ).hexdigest()


@pytest.fixture
def native():
    executable = os.environ.get("REPEAT_TYPER_EXECUTABLE") or shutil.which(
        "repeattyper"
    )
    if not executable:
        pytest.skip(
            "Install the pinned RepeatTyper wheel or set REPEAT_TYPER_EXECUTABLE for scientific execution"
        )
    return (executable,)


@pytest.fixture
def type_pipeline(tmp_path, native):
    source = tmp_path / "input"
    source.mkdir()
    records = {
        "first": "AACTCCCTAGGAGGGAGTAGAAG",
        "copy": "AACTCCCTAGGAGGGAGTAGAAG",
        "short": "ACGT",
        "ambiguous": "ACGTNNNN",
    }
    for i, (key, sequence) in enumerate(records.items()):
        (source / f"{i:06d}.fasta").write_text(f">{key}\n{sequence}\n")
    output = tmp_path / "output"
    run_repeat_type_pipeline(
        source,
        output,
        runner_config=RunnerConfig(
            ("unused",), ("unused",), ("unused",), repeattyper=native
        ),
        timeout_seconds=60,
    )
    claim = ClaimedJob(
        "a" * 32, "repeat_type", 1, 0, "bona_fide_possible", 1, True, True, ()
    )
    bindings = {k: (len(v), sha(v)) for k, v in records.items()}
    return claim, output, bindings, source


def validated(data):
    claim, output, bindings, source = data
    return _validated_pipeline_status(claim, output, bindings, input_root=source)


def test_native_report_and_complete_exports(type_pipeline):
    assert validated(type_pipeline) == "completed"
    output = type_pipeline[1]
    preview = repeat_type_summary(output, max_rows=1)
    assert preview["results_truncated"]
    write_repeat_type_exports(output, {"repeat_type": preview})
    public = json.loads((output / "repeat-typing.json").read_text())
    assert len(public["results"]) == 4
    assert public["counts"]["distinct_classified_sequences"] == 2
    assert public["results"][3]["native_score"] is None
    assert len(public["model"]["labels"]) == 37
    assert (
        public["model"]["upstream_commit"] == "a2d51831d4802533bb113bac23cf9df2853063b0"
    )
    assert "/tmp/" not in json.dumps(public)
    assert len((output / "repeat-typing.tsv").read_text().splitlines()) == 5


@pytest.mark.parametrize(
    "mutation",
    [
        "sequence",
        "id",
        "label",
        "score",
        "nan",
        "boolean",
        "status",
        "unknown",
        "model",
        "vocabulary",
        "threshold",
        "counts",
        "duplicate",
        "source",
        "stage",
        "symlink",
        "duplicate_key",
    ],
)
def test_independent_tamper_gate(type_pipeline, mutation):
    claim, output, bindings, source = type_pipeline
    path = output / "repeat-type/report.json"
    doc = json.loads(path.read_text())
    row = doc["results"][0]
    if mutation == "sequence":
        row["sequence"] = "AAAA"
    elif mutation == "id":
        row["id"] = "forged"
    elif mutation == "label":
        row["winner"] = "VII"
    elif mutation == "score":
        row["native_score"] = 1.5
    elif mutation == "nan":
        row["native_score"] = float("nan")
    elif mutation == "boolean":
        row["native_score"] = True
    elif mutation == "status":
        row["status"] = "unsupported"
    elif mutation == "unknown":
        doc["results"][3]["prediction"] = "I-A"
    elif mutation == "model":
        doc["model"]["files"]["xgb_repeats.json"] = "0" * 64
    elif mutation == "vocabulary":
        doc["model"]["labels"].append("VII")
    elif mutation == "threshold":
        doc["policy"]["threshold"] = 0.9
    elif mutation == "counts":
        doc["counts"]["unsupported"] += 1
    elif mutation == "duplicate":
        doc["results"][1] = dict(row)
    elif mutation == "source":
        next(source.iterdir()).write_text(">changed\nACGT\n")
    elif mutation == "stage":
        stage_path = output / "integration_result.json"
        stage = json.loads(stage_path.read_text())
        stage["stages"]["repeattyper"]["model_sha256"] = "0" * 64
        stage_path.write_text(json.dumps(stage))
    path.write_text(json.dumps(doc))
    if mutation == "symlink":
        target = output / "changed.json"
        path.rename(target)
        path.symlink_to(target)
    elif mutation == "duplicate_key":
        path.write_text(
            path.read_text().replace(
                '"execution_completed": true',
                '"execution_completed": false, "execution_completed": true',
            )
        )
    assert validated(type_pipeline) is None
    if mutation in {
        "sequence",
        "id",
        "label",
        "score",
        "nan",
        "boolean",
        "status",
        "unknown",
        "model",
        "vocabulary",
        "threshold",
        "counts",
        "duplicate",
    }:
        # Original immutable records; the adapter gate must reject independently.
        records = {
            "first": "AACTCCCTAGGAGGGAGTAGAAG",
            "copy": "AACTCCCTAGGAGGGAGTAGAAG",
            "short": "ACGT",
            "ambiguous": "ACGTNNNN",
        }
        with pytest.raises(ValueError):
            validate_report(
                doc, records, sha("".join(f">{k}\n{v}\n" for k, v in records.items()))
            )


def type_settings(settings_factory, tmp_path):
    settings = configured_settings(settings_factory, tmp_path)
    manifest = settings.scientific_release_manifest
    doc = json.loads(manifest.read_text())
    doc["schema_version"] = "7.0.0"
    for key, name, version in (
        ("crisprmap", "CRISPRmap-v2", "2.0.0a1"),
        ("repeattyper", "RepeatTyper", "2.0.0"),
    ):
        doc["scientific_artifacts"][key] = {
            **doc["scientific_artifacts"]["crisprspacer"],
            "name": key,
            "display_name": name,
            "version": version,
        }
    manifest.write_text(json.dumps(doc))
    return replace(
        settings, scientific_release_manifest_sha256=sha(manifest.read_bytes())
    )


def test_admission_requires_pinned_runtime_and_bounds(settings_factory, tmp_path):
    with pytest.raises(ReferenceUnavailableError):
        JobService(settings_factory()).submit(
            JobSubmission(mode="repeat_type", sequence=">r\nACGT\n"), "127.0.0.1"
        )
    service = JobService(type_settings(settings_factory, tmp_path))
    for sequence in (
        ">r\n" + "A" * 201,
        "".join(f">r{i}\nACGT\n" for i in range(1001)),
    ):
        with pytest.raises(FastaValidationError):
            service.submit(
                JobSubmission(mode="repeat_type", sequence=sequence), "127.0.0.1"
            )


def test_http_worker_handoff_download_and_uncertainty(
    settings_factory, tmp_path, native
):
    import io
    import zipfile
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app
    from .test_sequence_exchange import detection_ports

    packet = detection_ports(tmp_path)["repeats.json"]
    settings = replace(
        type_settings(settings_factory, tmp_path),
        integration_command=(sys.executable, "-m", "crispr_integration"),
    )
    settings.runner_config.write_text(
        json.dumps(
            {
                "crispridentify": ["unused"],
                "spacerplacer": ["unused"],
                "crispr_evor": ["unused"],
                "repeattyper": list(native),
            }
        )
    )
    with TestClient(create_app(settings)) as client:
        assert "repeat_type" in client.get("/api/v1/config").json()["modes"]
        created = client.post(
            "/api/v1/jobs",
            json={
                "mode": "repeat_type",
                "sequence": json.dumps(packet),
                "tool_options": {"repeattyper": {"--threshold": 1}},
            },
        )
        assert created.status_code == 202, created.text
        data = created.json()
        job_id = data["job"]["job_id"]
        headers = {"Authorization": "Bearer " + data["access_token"]}
        assert Worker(settings, worker_id="repeattyper-test").run_once()
        job = client.get(f"/api/v1/jobs/{job_id}", headers=headers).json()
        assert job["status"] == "completed", job
        assert job["summary"]["repeat_type"]["counts"]["uncertain"] == 3
        assert client.get(f"/api/v1/jobs/{job_id}").status_code in {401, 404}
        exports = {a["name"]: a for a in job["artifacts"]}
        for name in (
            "repeat-typing.json",
            "repeat-typing.tsv",
            "result.json",
            "input-sequences.json",
            "results.zip",
            "provenance.json",
        ):
            assert client.get(exports[name]["download_url"]).status_code in {401, 404}
            response = client.get(exports[name]["download_url"], headers=headers)
            assert response.status_code == 200
            if name == "input-sequences.json":
                assert response.json() == packet
            if name == "result.json":
                compact = response.json()
                assert (
                    compact["complete"] is True and compact["module"] == "repeat_type"
                )
                assert (
                    compact["data"]["results"][0]["origin"]
                    == packet["records"][0]["origin"]
                )
                assert compact["data"]["results"][0]["prediction"] is None
            if name == "results.zip":
                with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
                    assert "results/repeat-typing.json" in archive.namelist()
                    assert not any(
                        n.startswith(
                            ("results/repeat-type/", "results/repeat-type-input/")
                        )
                        for n in archive.namelist()
                    )


def test_release_schema_seven_does_not_require_unrelated_optional_tools(settings_factory, tmp_path):
    from crispr_web.config import _validated_scientific_manifest

    settings = type_settings(settings_factory, tmp_path)
    doc = json.loads(settings.scientific_release_manifest.read_text())
    for key in ("crisprmap", "crisprspacer", "crisprleader", "crisprrepeat"):
        doc["scientific_artifacts"].pop(key, None)
    assert "repeattyper" in _validated_scientific_manifest(doc)
    doc["scientific_artifacts"].pop("repeattyper")
    with pytest.raises(ValueError, match="artifact set"):
        _validated_scientific_manifest(doc)
