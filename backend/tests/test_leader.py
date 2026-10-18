from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import hashlib
import json
import shutil
import sqlite3
from contextlib import closing

import pytest
from pydantic import ValidationError

from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.models import JobSubmission
from crispr_web.publication.annotations import write_annotation_exports
from crispr_web.service import JobService
from crispr_web.summaries.leader import leader_summary
from crispr_web.worker import Worker

FIXTURE = Path(__file__).parent / "fixtures/leader"


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value) + "\n")


@pytest.fixture
def leader_output(tmp_path):
    root, inputs = tmp_path / "output", tmp_path / "input"
    inputs.mkdir()
    shutil.copyfile(FIXTURE / "source.fasta", inputs / "record.fasta")
    identify = json.loads((FIXTURE / "identify.json").read_text())
    native = json.loads((FIXTURE / "result.json").read_text())
    identity = identify["source"]["id"]
    write(root / "crispridentify/source/report.json", identify)
    index = dict(
        schema_version="1.0.0",
        analysis_mode="context_extraction",
        flank_length=37,
        records=[
            dict(key="record-000001", source_id=identity, native_source_id=identity)
        ],
        array_count=len(native["arrays"]),
        context_count=len(native["candidates"]),
    )
    pipeline = dict(
        schema=dict(name="crispr-tool-integration-run", version="1.0.0"),
        status="completed",
        workflow="leader",
        stages=dict(
            crispridentify=dict(status="completed"),
            crisprleader=dict(
                status="completed",
                analysis_mode="context_extraction",
                prediction_performed=False,
                record_count=1,
                array_count=index["array_count"],
                context_count=index["context_count"],
                flank_length=37,
            ),
        ),
    )

    def flush():
        write(root / "leader/index.json", index)
        write(root / "integration_result.json", pipeline)
        result = root / "leader/record-000001/result.json"
        write(result, native)
        write(
            result.with_name("manifest.json"),
            dict(
                schema_version="1.0.0",
                files={
                    "result.json": dict(
                        sha256=hashlib.sha256(result.read_bytes()).hexdigest(),
                        size_bytes=result.stat().st_size,
                    )
                },
            ),
        )

    flush()
    claimed = ClaimedJob(
        "a" * 32,
        "leader",
        1,
        0,
        "bona_fide_possible",
        1,
        True,
        True,
        (),
        leader_flank_length=37,
    )
    bindings = {identity: (identify["source"]["length"], identify["source"]["sha256"])}
    return root, inputs, native, index, pipeline, claimed, bindings, flush


def validate(data):
    root, inputs, _, _, _, claimed, bindings, _ = data
    return _validated_pipeline_status(claimed, root, bindings, input_root=inputs)


def test_context_reconstruction_full_export_and_preview(leader_output):
    assert validate(leader_output) == "completed"
    root, _, native, _, _, _, _, _ = leader_output
    preview = leader_summary(root, max_rows=1)
    assert preview["truncated"] is True and len(preview["contexts"]) == 1
    assert preview["prediction_performed"] is False
    write_annotation_exports(root, {"crisprleader": preview})
    full = json.loads((root / "leader-contexts.json").read_text())
    assert full["complete"] is True and full["truncated"] is False
    assert len(full["contexts"]) == len(native["candidates"]) == 2
    assert full["contexts"][1]["sequence"] == native["candidates"][1]["sequence"]
    assert "no_leader_prediction" in (root / "leader-contexts.gff3").read_text()
    gff = [
        line.split("\t")
        for line in (root / "leader-contexts.gff3").read_text().splitlines()
        if not line.startswith("#")
    ]
    bed = [
        line.split("\t")
        for line in (root / "leader-contexts.bed").read_text().splitlines()
        if not line.startswith("#")
    ]
    assert int(gff[0][3]) == int(bed[0][1]) + 1 and int(gff[0][4]) == int(bed[0][2])
    assert (root / "leader-contexts.fasta").read_text().count(">") == 2
    assert "identify_original" not in (root / "leader-contexts.json").read_text()


@pytest.mark.parametrize(
    "mutation",
    [
        "sequence",
        "coordinates",
        "strand",
        "side",
        "source",
        "anchor",
        "metadata",
        "predictions",
        "reference",
        "mode",
        "settings",
        "duplicate",
        "omitted",
        "stage_count",
        "index_source",
        "index_path",
        "manifest",
        "symlink",
        "native_array",
        "forged_repeat",
    ],
)
def test_context_evidence_fails_closed_even_with_recomputed_hashes(
    leader_output, mutation
):
    root, _, report, index, pipeline, _, _, flush = leader_output
    candidate = report["candidates"][0]
    if mutation == "sequence":
        candidate["sequence"] = "A" * len(candidate["sequence"])
    if mutation == "coordinates":
        candidate["segments"][0]["start"] += 1
    if mutation == "strand":
        candidate["strand"] = "-"
    if mutation == "side":
        candidate["side"] = "right"
    if mutation == "source":
        report["sources"][0]["sha256"] = "0" * 64
    if mutation == "anchor":
        candidate["anchor"]["sequence"] = "A"
    if mutation == "metadata":
        candidate["metadata"]["window_is_biological_boundary"] = True
    if mutation == "predictions":
        report["predictions"] = [{"status": "supported"}]
    if mutation == "reference":
        report["reference"] = {"bundle_id": "fake"}
    if mutation == "mode":
        report["analysis_mode"] = "prediction"
    if mutation == "settings":
        report["effective_config"]["flank_length"] = 500
    if mutation == "duplicate":
        report["candidates"].append(dict(candidate))
    if mutation == "omitted":
        report["candidates"].pop()
    if mutation == "stage_count":
        pipeline["stages"]["crisprleader"]["array_count"] += 1
    if mutation == "index_source":
        index["records"][0]["source_id"] = "wrong"
    if mutation == "index_path":
        index["records"][0]["key"] = "../record-000001"
    if mutation == "native_array":
        report["arrays"][0]["segments"][0]["end"] -= 1
    if mutation == "forged_repeat":
        identify_path = root / "crispridentify/source/report.json"
        identify = json.loads(identify_path.read_text())
        for row in identify["arrays"]:
            row["repeats"][0]["sequence"] = "A" * len(row["repeats"][0]["sequence"])
        write(identify_path, identify)
    flush()
    result = root / "leader/record-000001/result.json"
    if mutation == "manifest":
        result.write_text(result.read_text() + " ")
    if mutation == "symlink":
        destination = root.parent / "untrusted.json"
        result.rename(destination)
        result.symlink_to(destination)
    assert validate(leader_output) is None


@pytest.mark.parametrize("value", [0, 5001, -1, True, 1.5, "500", None])
def test_window_length_input_bounds(value):
    with pytest.raises(ValidationError):
        JobSubmission(mode="leader", sequence="ACGT", leader_flank_length=value)


def test_schema_three_migration_retains_existing_job_and_enables_leader(
    settings_factory, tmp_path
):
    settings = settings_factory()
    service = JobService(settings)
    token, job = service.submit(
        JobSubmission(mode="detection", sequence="ACGT"), "198.51.100.1"
    )
    assert Worker(settings, store=service.store).run_once()
    before = service.get(job["job_id"])
    with closing(service.store.connect()) as conn:
        conn.execute("ALTER TABLE jobs DROP COLUMN leader_flank_length")
        conn.execute("UPDATE metadata SET value='3' WHERE name='schema_version'")
        dump = "\n".join(conn.iterdump()).replace(", 'leader'", "")
    legacy = tmp_path / "schema-three.sqlite3"
    with closing(sqlite3.connect(legacy)) as conn:
        conn.executescript(dump)
    upgraded = JobService(replace(settings, database_path=legacy))
    assert upgraded.authorize(job["job_id"], token)
    assert upgraded.get(job["job_id"]) == before
    _, leader = upgraded.submit(
        JobSubmission(mode="leader", sequence="ACGT", leader_flank_length=5000),
        "198.51.100.2",
    )
    claimed = upgraded.store.claim_next("leader-worker")
    assert claimed.job_id == leader["job_id"] and claimed.leader_flank_length == 5000
    with closing(upgraded.store.connect()) as conn:
        assert conn.execute("PRAGMA foreign_key_check").fetchall() == []


def test_clipped_contexts_preserve_both_sides_and_unknown_topology(leader_output):
    root, inputs, report, index, pipeline, claimed, bindings, flush = leader_output
    report.clear()
    report.update(json.loads((FIXTURE / "clipped-result.json").read_text()))
    index["flank_length"] = 5000
    pipeline["stages"]["crisprleader"]["flank_length"] = 5000
    flush()
    assert (
        _validated_pipeline_status(
            replace(claimed, leader_flank_length=5000),
            root,
            bindings,
            input_root=inputs,
        )
        == "completed"
    )
    evidence = leader_summary(root, max_rows=100)
    assert {row["side"] for row in evidence["contexts"]} == {"left", "right"}
    assert all(
        row["truncated_reason"] == "unknown_topology_boundary"
        for row in evidence["contexts"]
    )
    assert all(row["available_length"] < 5000 for row in evidence["contexts"])
