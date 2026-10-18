"""Real native SpacePHARER pipeline, explicit grouping and independent publication gates."""

from dataclasses import replace
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest
from crispr_integration.runner import RunnerConfig
from crispr_integration.association import (
    run_association_pipeline,
    reference_data,
    validate_report,
)
from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.execution.validation_association import (
    reference_data as backend_reference,
)
from crispr_web.exchange import association_groups, sequence_set, sequence_input
from crispr_web.fasta import FastaValidationError
from crispr_web.models import JobSubmission
from crispr_web.service import JobService, ReferenceUnavailableError
from crispr_web.worker import Worker
from .test_repeat_type import type_settings

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tools/SpacePHARER"
NATIVE = os.environ.get("SPACEPHARER_EXECUTABLE") or shutil.which("spacepharer")
def sha(value):
    return hashlib.sha256(value.encode() if isinstance(value, str) else value).hexdigest()


@pytest.fixture(scope="module")
def native():
    if not NATIVE:
        pytest.skip(
            "Install pinned SpacePHARER CLI/native or set SPACEPHARER_EXECUTABLE"
        )
    return (NATIVE,)


@pytest.fixture(scope="module")
def panel(tmp_path_factory, native):
    root = tmp_path_factory.mktemp("association-panel")
    files = []
    for path in sorted((SOURCE / "examples").glob("*.gz")):
        target = root / path.stem
        target.write_bytes(gzip.decompress(path.read_bytes()))
        files += ["--genome", str(target)]
    subprocess.run(
        [
            *native,
            "prepare-reference",
            *files,
            "--name",
            "upstream-three-phages",
            "--output",
            str(root / "panel"),
        ],
        check=True,
        capture_output=True,
        timeout=120,
    )
    return root / "panel", sha((root / "panel/reference.json").read_bytes())


def upstream_records():
    lines = (SOURCE / "tests/fixtures/spacers.fasta").read_text().splitlines()
    return {lines[i][1:]: lines[i + 1] for i in range(0, len(lines), 2)}


@pytest.fixture(scope="module")
def original_pipeline(tmp_path_factory, native, panel):
    root = tmp_path_factory.mktemp("association-run")
    source = root / "input/records"
    source.mkdir(parents=True)
    records = upstream_records()
    for i, (key, seq) in enumerate(records.items()):
        (source / f"{i:06}.fasta").write_text(f">{key}\n{seq}\n")
    groups = root / "association-groups.json"
    groups.write_bytes((SOURCE / "tests/fixtures/groups.json").read_bytes())
    run_association_pipeline(
        source,
        root / "output",
        runner_config=RunnerConfig(
            ("unused",), ("unused",), ("unused",), spacepharer=native
        ),
        reference_path=panel[0],
        reference_sha256=panel[1],
        groups_path=groups,
        timeout_seconds=120,
    )
    return root


@pytest.fixture
def pipeline(tmp_path, original_pipeline, panel):
    root = tmp_path / "job"
    shutil.copytree(original_pipeline, root)
    source = root / "input/records"
    claim = ClaimedJob(
        "a" * 32,
        "spacer_association",
        1,
        0,
        "bona_fide_possible",
        1,
        True,
        True,
        (),
        association_reference_sha256=panel[1],
        association_grouping="per_source",
        association_groups_sha256=sha((root / "association-groups.json").read_bytes()),
    )
    records = upstream_records()
    return (
        claim,
        root / "output",
        {k: (len(v), sha(v)) for k, v in records.items()},
        source,
        panel[0],
    )


def validated(pipeline):
    claim, output, bindings, source, reference = pipeline
    return _validated_pipeline_status(
        claim, output, bindings, input_root=source, association_reference=reference
    )


def test_actual_native_grouped_support_and_exports(pipeline):
    from crispr_web.summaries.association import association_summary
    from crispr_web.publication.association import write_association_exports

    assert validated(pipeline) == "completed"
    out = pipeline[1]
    doc = json.loads((out / "association/report.json").read_text())
    assert [
        (r["group_id"], r["target_id"], r["combined_score"], r["native_hit_count"])
        for r in doc["associations"]
    ] == [
        ("CP003088", "AY369265.2", 12.24, 1),
        ("JNGQ01000029", "AJ242593.1", 74.56, 5),
    ]
    assert all(
        r["native_fdr"] is None and r["status"] == "uncalibrated_candidate"
        for r in doc["associations"]
    )
    assert any(
        h["query_start"] > h["query_end"]
        for r in doc["associations"]
        for h in r["hits"]
    )
    preview = association_summary(out, max_rows=1)
    assert preview["associations_truncated"]
    write_association_exports(out, {"spacer_association": preview})
    public = json.loads((out / "spacer-associations.json").read_text())
    assert len(public["associations"]) == 2 and public["counts"]["support_hits"] == 6
    assert "/tmp/" not in json.dumps(public)
    assert len((out / "spacer-associations.tsv").read_text().splitlines()) == 3


@pytest.mark.parametrize(
    "tamper",
    [
        "null_to_zero",
        "accepted",
        "score",
        "boolean",
        "nan",
        "group",
        "target",
        "sequence",
        "alignment",
        "double_gap",
        "reverse_coordinate",
        "pvalue",
        "occurrence",
        "missing_duplicate",
        "hitcount",
        "calibration",
        "receipt",
        "reference",
        "group_hash",
        "input_hash",
        "stage",
        "symlink",
        "duplicate_key",
    ],
)
def test_independent_output_mutations(pipeline, tamper):
    claim, out, bindings, source, reference = pipeline
    path = out / "association/report.json"
    doc = json.loads(path.read_text())
    row = doc["associations"][0]
    hit = row["hits"][0]
    if tamper == "null_to_zero":
        row["native_fdr"] = 0
    elif tamper == "accepted":
        row["status"] = "association"
    elif tamper == "score":
        row["combined_score"] = -1
    elif tamper == "boolean":
        row["native_hit_count"] = True
    elif tamper == "nan":
        hit["p_best_hit"] = float("nan")
    elif tamper == "group":
        row["group_id"] = "unrelated"
    elif tamper == "target":
        row["target_id"] = "unrelated"
    elif tamper == "sequence":
        hit["query_sequence"] = "ACGT"
    elif tamper == "alignment":
        hit["target_alignment"] = "AAAA"
    elif tamper == "double_gap":
        hit["query_alignment"] = "-" + hit["query_alignment"]
        hit["target_alignment"] = "-" + hit["target_alignment"]
    elif tamper == "reverse_coordinate":
        hit["query_start"] = hit["query_end"]
    elif tamper == "pvalue":
        hit["p_best_hit"] = 2
    elif tamper == "occurrence":
        hit["occurrence_ids"] = ["unrelated"]
    elif tamper == "missing_duplicate":
        hit["occurrence_ids"] = []
    elif tamper == "hitcount":
        row["native_hit_count"] += 1
    elif tamper == "calibration":
        doc["calibration_status"] = "estimated"
    elif tamper == "receipt":
        doc["calibration"] = None
    elif tamper == "reference":
        doc["reference"]["manifest_sha256"] = "0" * 64
    elif tamper == "group_hash":
        doc["groups_sha256"] = "0" * 64
    elif tamper == "input_hash":
        doc["input_sha256"] = "0" * 64
    elif tamper == "stage":
        stage = out / "integration_result.json"
        value = json.loads(stage.read_text())
        value["stages"]["spacepharer"]["counts"]["groups"] = True
        stage.write_text(json.dumps(value))
    path.write_text(json.dumps(doc))
    if tamper == "symlink":
        target = out / "moved.json"
        path.rename(target)
        path.symlink_to(target)
    elif tamper == "duplicate_key":
        path.write_text(
            path.read_text().replace(
                '"requested_fdr":', '"requested_fdr": 0.1, "requested_fdr":'
            )
        )
    assert validated(pipeline) is None
    if tamper not in {"stage", "symlink", "duplicate_key"}:
        ref, targets = reference_data(reference, claim.association_reference_sha256)
        with pytest.raises(ValueError):
            validate_report(
                doc,
                upstream_records(),
                json.loads(
                    (source.parent.parent / "association-groups.json").read_text()
                ),
                ref,
                targets,
                sha("".join(f">{k}\n{v}\n" for k, v in upstream_records().items())),
                claim.association_groups_sha256,
                0.05,
            )


@pytest.mark.parametrize("mutation", ["missing", "extra", "db_symlink", "digest"])
def test_reference_boundary_both_gates(panel, tmp_path, mutation):
    root = tmp_path / "panel"
    shutil.copytree(panel[0], root)
    if mutation == "missing":
        doc = json.loads((root / "reference.json").read_text())
        doc["files"].pop("db/target.index")
        (root / "db/target.index").unlink()
        (root / "reference.json").write_text(json.dumps(doc))
    elif mutation == "extra":
        (root / "db/extra").write_text("unlisted")
    elif mutation == "db_symlink":
        (root / "db").rename(tmp_path / "db")
        (root / "db").symlink_to(tmp_path / "db", target_is_directory=True)
    elif mutation == "digest":
        (root / "db/target").write_text("changed")
    digest = sha((root / "reference.json").read_bytes())
    for check in (reference_data, backend_reference):
        with pytest.raises(ValueError):
            check(root, digest)


def association_settings(settings_factory, tmp_path, panel):
    settings = type_settings(settings_factory, tmp_path)
    path = settings.scientific_release_manifest
    doc = json.loads(path.read_text())
    doc["schema_version"] = "8.0.0"
    doc["scientific_artifacts"]["spacepharer"] = {
        **doc["scientific_artifacts"]["repeattyper"],
        "name": "spacepharer",
        "display_name": "SpacePHARER",
        "version": "2.1.0",
        "repository": "https://github.com/Alexander-Mitrofanov/SpacePHARER.git",
        "commit": "3121a8eb3fd5b79996c7175552a19aa465486458",
    }
    path.write_text(json.dumps(doc))
    return replace(
        settings,
        scientific_release_manifest_sha256=sha(path.read_bytes()),
        association_reference_path=panel[0],
        association_reference_sha256=panel[1],
    )


def test_grouping_admission_and_missing_panel(settings_factory, tmp_path, panel):
    from pydantic import ValidationError

    with pytest.raises(ValidationError):
        JobSubmission(mode="spacer_association", sequence="ACGT")
    with pytest.raises(ReferenceUnavailableError):
        JobService(settings_factory()).submit(
            JobSubmission(
                mode="spacer_association",
                sequence="ACGT",
                association_grouping="single_set",
            ),
            "127.0.0.1",
        )
    service = JobService(association_settings(settings_factory, tmp_path, panel))
    for sequence, group in [
        ("A" * 201, "single_set"),
        ("ACGT", "per_source"),
        ("".join(f">r{i}\nACGT\n" for i in range(101)), "per_record"),
    ]:
        with pytest.raises(FastaValidationError):
            service.submit(
                JobSubmission(
                    mode="spacer_association",
                    sequence=sequence,
                    association_grouping=group,
                ),
                "127.0.0.1",
            )


def test_http_real_worker_source_groups_and_protected_downloads(
    settings_factory, tmp_path, panel, native, original_pipeline
):
    import io
    import zipfile
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app

    rows = []
    groups = json.loads((SOURCE / "tests/fixtures/groups.json").read_text())["groups"]
    records = upstream_records()
    for group in groups:
        for i, key in enumerate(group["members"]):
            seq = records[key]
            rows.append(
                {
                    "id": f"s{len(rows) + 1}",
                    "sequence": seq,
                    "origin": {
                        "source_id": group["id"],
                        "source_sha256": sha(group["id"]),
                        "array_id": "submitted-array",
                        "ordinal": i + 1,
                        "start": i * 200 + 1,
                        "end": i * 200 + len(seq),
                    },
                }
            )
    # Duplicate identity occurrence and unsupported record retain their own origins.
    support = json.loads(
        (original_pipeline / "output/association/report.json").read_text()
    )["associations"][0]["hits"][0]["query_sequence"]
    rows.append(
        {**next(row for row in rows if row["sequence"] == support), "id": "duplicate"}
    )
    rows.append({"id": "ambiguous", "sequence": "NNNN"})
    packet = sequence_set("spacer", "DNA", rows[:-1])
    settings = replace(
        association_settings(settings_factory, tmp_path, panel),
        integration_command=(sys.executable, "-m", "crispr_integration"),
        max_records=1000,
    )
    settings.runner_config.write_text(
        json.dumps(
            {
                "crispridentify": ["unused"],
                "spacerplacer": ["unused"],
                "crispr_evor": ["unused"],
                "spacepharer": list(native),
            }
        )
    )
    with TestClient(create_app(settings)) as client:
        assert "spacer_association" in client.get("/api/v1/config").json()["modes"]
        created = client.post(
            "/api/v1/jobs",
            json={
                "mode": "spacer_association",
                "sequence": json.dumps(packet),
                "association_grouping": "per_source",
            },
        )
        assert created.status_code == 202, created.text
        data = created.json()
        jid = data["job"]["job_id"]
        headers = {"Authorization": "Bearer " + data["access_token"]}
        assert Worker(settings, worker_id="spacepharer-test").run_once()
        job = client.get(f"/api/v1/jobs/{jid}", headers=headers).json()
        assert job["status"] == "completed", job
        evidence = job["summary"]["spacer_association"]
        assert (
            evidence["counts"]["groups"] == 2 and evidence["counts"]["candidates"] == 2
        )
        assert {r["group_id"] for r in evidence["associations"]} == {
            "source1",
            "source2",
        }
        assert any(
            "duplicate" in h["occurrence_ids"] and len(h["occurrence_ids"]) == 2
            for row in evidence["associations"]
            for h in row["hits"]
        )
        assert client.get(f"/api/v1/jobs/{jid}").status_code in {401, 404}
        exports = {a["name"]: a for a in job["artifacts"]}
        for name in (
            "spacer-associations.json",
            "spacer-associations.tsv",
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
                    compact["complete"] is True
                    and compact["module"] == "spacer_association"
                )
                assert compact["data"]["occurrences"] == packet["records"]
                assert compact["data"]["calibration_status"] == "unavailable"
            if name == "results.zip":
                with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
                    assert "results/spacer-associations.json" in archive.namelist()
                    assert not any(
                        n.startswith(
                            ("results/association/", "results/association-input/")
                        )
                        for n in archive.namelist()
                    )


def test_detection_spacer_handoff_never_silently_pools_sources(tmp_path):
    from .test_sequence_exchange import detection_ports

    packet = detection_ports(tmp_path)["spacers.json"]
    text, meta = sequence_input(
        json.dumps(packet), mode="spacer_association", molecule="DNA", max_bytes=100000
    )
    group = association_groups([r["id"] for r in packet["records"]], "per_source", meta)
    assert len(group["groups"]) == 1 and text.count(">") == len(packet["records"])
    # A second source ID never becomes part of the first group implicitly.
    meta["origins"][packet["records"][-1]["id"]]["source_id"] = "second-contig"
    assert (
        len(
            association_groups(
                [r["id"] for r in packet["records"]], "per_source", meta
            )["groups"]
        )
        == 2
    )


@pytest.mark.parametrize("case", ["detection", "unsupported", "negative"])
def test_real_http_controls_and_detection_occurrences(
    settings_factory, tmp_path, panel, native, case
):
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app
    from .test_sequence_exchange import detection_ports

    packet = (
        detection_ports(tmp_path)["spacers.json"]
        if case == "detection"
        else sequence_set(
            "spacer",
            "DNA",
            [
                {
                    "id": "control",
                    "sequence": "NNNN" if case == "unsupported" else "ACGT" * 9,
                }
            ],
        )
    )
    settings = replace(
        association_settings(settings_factory, tmp_path, panel),
        integration_command=(sys.executable, "-m", "crispr_integration"),
    )
    settings.runner_config.write_text(
        json.dumps(
            {
                "crispridentify": ["unused"],
                "spacerplacer": ["unused"],
                "crispr_evor": ["unused"],
                "spacepharer": list(native),
            }
        )
    )
    with TestClient(create_app(settings)) as client:
        created = client.post(
            "/api/v1/jobs",
            json={
                "mode": "spacer_association",
                "sequence": json.dumps(packet),
                "association_grouping": "per_source"
                if case == "detection"
                else "per_record",
            },
        )
        assert created.status_code == 202, created.text
        data = created.json()
        headers = {"Authorization": "Bearer " + data["access_token"]}
        assert Worker(settings, worker_id="association-controls").run_once()
        job = client.get(
            "/api/v1/jobs/" + data["job"]["job_id"], headers=headers
        ).json()
        assert job["status"] == "completed", job
        result = job["summary"]["spacer_association"]
        assert (
            result["counts"]["associations"] == 0
            and result["counts"]["candidates"] == 0
        )
        if case == "unsupported":
            assert (
                result["calibration_status"] == "not_run"
                and result["calibration"] is None
            )
            assert result["counts"]["unsupported"] == 1
        elif case == "negative":
            assert result["calibration_status"] == "unavailable"
        inputs = next(
            a for a in job["artifacts"] if a["name"] == "input-sequences.json"
        )
        assert client.get(inputs["download_url"], headers=headers).json() == packet


def test_reference_identity_participates_in_readiness(
    settings_factory, tmp_path, panel
):
    local = tmp_path / "panel"
    shutil.copytree(panel[0], local)
    settings = association_settings(settings_factory, tmp_path, (local, panel[1]))
    assert settings.scientific_release_ok()
    (local / "reference.json").write_text("{}")
    assert not settings.scientific_release_ok()


def test_schema_eight_has_no_unrelated_optional_requirements(
    settings_factory, tmp_path, panel
):
    from crispr_web.config import _validated_scientific_manifest

    settings = association_settings(settings_factory, tmp_path, panel)
    doc = json.loads(settings.scientific_release_manifest.read_text())
    for name in (
        "repeattyper",
        "crisprmap",
        "crisprspacer",
        "crisprleader",
        "crisprrepeat",
    ):
        doc["scientific_artifacts"].pop(name, None)
    assert "spacepharer" in _validated_scientific_manifest(doc)
    doc["scientific_artifacts"].pop("spacepharer")
    with pytest.raises(ValueError, match="artifact set"):
        _validated_scientific_manifest(doc)
