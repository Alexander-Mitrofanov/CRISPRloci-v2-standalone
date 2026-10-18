from __future__ import annotations

from dataclasses import replace
import hashlib
import json
import sqlite3
from contextlib import closing

import pytest
from fastapi.testclient import TestClient

from crispr_web.api import create_app
from crispr_web.db import ClaimedJob, Store
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.models import JobSubmission
from crispr_web.publication.annotations import write_annotation_exports
from crispr_web.service import JobService
from crispr_web.summaries.annotations import annotation_summary
from crispr_web.worker import Worker


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


def write_json(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document) + "\n")


@pytest.fixture
def annotation_output(tmp_path):
    root, inputs = tmp_path / "output", tmp_path / "input"
    root.mkdir()
    inputs.mkdir()
    raw = b">contig\nACGTACGT\n"
    (inputs / "record-1.fasta").write_bytes(raw)
    record = {
        "key": "f000001.s000001",
        "source": {
            "id": "contig",
            "length": 8,
            "sha256": sha(b"ACGTACGT"),
            "file_sha256": sha(raw),
        },
        "predictions": [
            {
                "export_id": "f000001.s000001.p000001",
                "interval": {"start": 2, "end": 4, "length": 3},
                "strand": "-",
                "sequence": "ACG",
                "ranking": {"score": 0.0},
            }
        ],
    }
    report = {
        "schema": {"name": "CRISPRtracrRNA-report", "version": "3.0.0"},
        "status": "complete",
        "configuration": {
            "mode": "model",
            "model_type": "V",
            "short_anti_repeat_rescue": False,
        },
        "coordinate_system": {
            "indexing": 1,
            "end_inclusive": True,
            "reference": "input_forward_strand",
        },
        "summary": {"sequence_count": 1, "prediction_count": 1},
        "records": [record],
    }

    def flush():
        retained = root / "tracrrna/work/f000001.s000001/record.json"
        write_json(retained, record)
        write_json(root / "tracrrna/report.json", report)
        write_json(
            root / "tracrrna/manifest.json",
            {
                "algorithm": "sha256",
                "files": {
                    "work/f000001.s000001/record.json": sha(retained.read_bytes())
                },
            },
        )

    flush()
    write_json(
        root / "integration_result.json",
        {
            "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
            "status": "completed",
            "workflow": "tracrrna",
            "stages": {
                "tracrrna": {
                    "status": "completed",
                    "mode": "model",
                    "model_type": "V",
                    "record_count": 1,
                    "prediction_count": 1,
                }
            },
        },
    )
    claimed = ClaimedJob(
        "a" * 32,
        "tracrrna",
        1,
        0,
        "bona_fide_possible",
        1,
        True,
        True,
        (),
        tracr_model_type="V",
    )
    return root, inputs, record, report, flush, claimed


def validated(fixture):
    root, inputs, _, _, _, claimed = fixture
    return _validated_pipeline_status(
        claimed, root, {"contig": (8, sha(b"ACGTACGT"))}, input_root=inputs
    )


def test_reverse_strand_annotation_binding_and_full_export(annotation_output):
    root, _, _, _, _, _ = annotation_output
    assert validated(annotation_output) == "completed"
    warnings = []
    cas, tracr, _ = annotation_summary(
        root, max_rows=0, accepted_categories=[], warnings=warnings
    )
    assert tracr["truncated"] is True
    write_annotation_exports(root, {"casandra": cas, "tracrrna": tracr})
    exported = json.loads((root / "annotations.json").read_text())
    assert exported["complete"] is True
    assert exported["tracrrna"]["candidates"][0]["sequence"] == "ACG"
    assert exported["tracrrna"]["candidates"][0]["score"] == 0.0
    assert exported["tracrrna"]["candidates"][0]["evalue"] is None
    assert exported["tracrrna"]["truncated"] is False
    assert "/work/" not in (root / "annotations.json").read_text()
    assert (
        root / "tracr-predictions.fasta"
    ).read_text() == ">f000001.s000001.p000001\nACG\n"


@pytest.mark.parametrize(
    "mutation",
    [
        "sequence",
        "coordinates",
        "identity",
        "input_hash",
        "duplicate",
        "model",
        "count",
        "nonfinite",
        "missing",
        "symlink",
    ],
)
def test_rejects_stale_or_malformed_annotation_evidence(annotation_output, mutation):
    root, inputs, record, report, flush, _ = annotation_output
    prediction = record["predictions"][0]
    if mutation == "sequence":
        prediction["sequence"] = "CGT"
    if mutation == "coordinates":
        prediction["interval"]["end"] = 9
    if mutation == "identity":
        record["source"]["id"] = "other"
    if mutation == "input_hash":
        record["source"]["file_sha256"] = "0" * 64
    if mutation == "duplicate":
        record["predictions"].append(dict(prediction))
        report["summary"]["prediction_count"] = 2
    if mutation == "model":
        report["configuration"]["model_type"] = "II"
    if mutation == "count":
        report["summary"]["prediction_count"] = 0
    if mutation == "nonfinite":
        prediction["ranking"]["score"] = float("inf")
    flush()
    retained = root / "tracrrna/work/f000001.s000001/record.json"
    if mutation == "missing":
        retained.unlink()
    if mutation == "symlink":
        retained.unlink()
        retained.symlink_to(inputs / "record-1.fasta")
    assert validated(annotation_output) is None


@pytest.mark.parametrize("mode", ["loci", "cas", "tracrrna", "leader"])
def test_annotation_options_queue_and_cli_contract(settings_factory, mode):
    settings = settings_factory()
    service = JobService(settings)
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/v1/jobs",
            json={"sequence": ">contig\nACGT\n", "mode": mode, "tracr_model_type": "V", "leader_flank_length": 37},
        )
        assert response.status_code == 202
        payload = response.json()
        job = payload["job"]
        assert job["options"]["tracr_model_type"] == "V"
        assert job["options"]["leader_flank_length"] == 37
        assert service.authorize(job["job_id"], payload["access_token"])
        claimed = service.store.claim_next("annotation-contract")
        assert claimed.mode == mode and claimed.tracr_model_type == "V"
        command = Worker(settings, store=service.store).command_for(
            claimed, service.job_root(claimed.job_id) / "output"
        )
        if mode in {"loci", "leader"}:
            assert command[command.index("--leader-flank-length") + 1] == "37"
        if mode in {"loci", "tracrrna"}:
            assert command[command.index("--tracr-model-type") + 1] == "V"
        assert (
            client.post(
                "/api/v1/jobs",
                json={"sequence": "ACGT", "mode": mode, "tracr_model_type": "VI"},
            ).status_code
            == 422
        )


def test_schema_two_migration_preserves_capabilities_jobs_artifacts_and_salt(
    settings_factory, tmp_path
):
    settings = settings_factory()
    service = JobService(settings)
    token, submitted = service.submit(
        JobSubmission(sequence=">a\nACGT\n", mode="detection"), "198.51.100.1"
    )
    assert Worker(settings, store=service.store).run_once()
    before = service.get(submitted["job_id"])
    assert before["status"] == "completed" and before["artifacts"]
    with closing(service.store.connect()) as connection:
        connection.execute("ALTER TABLE jobs DROP COLUMN tracr_model_type")
        connection.execute("ALTER TABLE jobs DROP COLUMN leader_flank_length")
        connection.execute("UPDATE metadata SET value='2' WHERE name='schema_version'")
        dump = "\n".join(connection.iterdump())
        salt = connection.execute(
            "SELECT value FROM metadata WHERE name='rate_salt'"
        ).fetchone()[0]
    dump = dump.replace(", 'cas', 'tracrrna', 'loci', 'leader'", "")
    legacy = tmp_path / "legacy.sqlite3"
    with closing(sqlite3.connect(legacy)) as connection:
        connection.executescript(dump)
    migrated_settings = replace(settings, database_path=legacy)
    migrated = JobService(migrated_settings)
    assert migrated.authorize(submitted["job_id"], token)
    assert migrated.get(submitted["job_id"]) == before
    with closing(migrated.store.connect()) as connection:
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert (
            connection.execute(
                "SELECT value FROM metadata WHERE name='schema_version'"
            ).fetchone()[0]
            == "7"
        )
        assert (
            connection.execute(
                "SELECT value FROM metadata WHERE name='rate_salt'"
            ).fetchone()[0]
            == salt
        )
    _, created = migrated.submit(
        JobSubmission(sequence=">b\nACGT\n", mode="loci", tracr_model_type="V"),
        "198.51.100.2",
    )
    claimed = migrated.store.claim_next("migrated")
    assert claimed.job_id == created["job_id"] and claimed.tracr_model_type == "V"
    assert Store(migrated_settings).database_ok()


@pytest.mark.parametrize(
    "corrupt", [None, "sequence", "checksum", "count", "confidence", "duplicate"]
)
def test_cas_native_evidence_and_public_summary(tmp_path, corrupt):
    root, inputs = tmp_path / "output", tmp_path / "input"
    inputs.mkdir()
    payload = b">contig\nACGTACGT\n"
    (inputs / "one.fasta").write_bytes(payload)
    native = root / "casandra/record-000001"
    native.mkdir(parents=True)
    gene = {
        "protein_id": "protein-1",
        "contig_id": "contig",
        "start_1based": 2,
        "end_1based_inclusive": 4,
        "strand": "-",
        "nucleotide_sequence": "ACG",
        "nucleotide_sha256": sha(b"ACG"),
        "prediction": {
            "is_cas": True,
            "classification": {"type": "II", "subtype": None},
            "score_margin": 0.0,
        },
    }
    cassette = {
        "cassette_id": "cassette-1",
        "contig_id": "contig",
        "start_1based": 2,
        "end_1based_inclusive": 4,
        "cas_gene_count": 1,
        "final_classification": {"type": "II", "subtype": None, "confidence": 0.0},
    }
    if corrupt == "sequence":
        gene["nucleotide_sequence"] = "CGT"
    if corrupt == "confidence":
        cassette["final_classification"]["confidence"] = "0.9"
    proteins = json.dumps(gene) + "\n"
    if corrupt == "duplicate":
        proteins += proteins
    (native / "proteins.jsonl").write_text(proteins)
    (native / "cassettes.jsonl").write_text(json.dumps(cassette) + "\n")
    write_json(
        native / "run.json",
        {
            "schema_version": 5,
            "program": "CasAndra",
            "inputs": [{"kind": "genome_fasta", "sha256": sha(payload)}],
            "genes": 2 if corrupt == "count" else 1,
            "cas_proteins": 1,
            "cas_cassettes": 1,
        },
    )
    write_json(
        native / "manifest.json",
        {
            "schema_version": 5,
            "files": {
                name: {"sha256": sha((native / name).read_bytes())}
                for name in ["proteins.jsonl", "cassettes.jsonl"]
            },
        },
    )
    if corrupt == "checksum":
        (native / "proteins.jsonl").write_text(proteins + "\n")
    write_json(
        root / "integration_result.json",
        {
            "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
            "status": "completed",
            "workflow": "casandra",
            "stages": {
                "casandra": {
                    "status": "completed",
                    "gene_count": 1,
                    "cas_gene_count": 1,
                    "cassette_count": 1,
                }
            },
        },
    )
    claimed = ClaimedJob("a" * 32, "cas", 1, 0, "bona_fide_possible", 1, True, True, ())
    result = _validated_pipeline_status(
        claimed, root, {"contig": (8, sha(b"ACGTACGT"))}, input_root=inputs
    )
    if corrupt:
        assert result is None
    else:
        assert result == "completed"
        warnings = []
        cas, _, _ = annotation_summary(
            root, max_rows=5, accepted_categories=[], warnings=warnings
        )
        assert warnings == [] and cas["cas_gene_count"] == 1
        assert cas["genes"][0]["score_margin"] == 0.0
        assert cas["cassettes"][0]["confidence"] == 0.0
        assert cas["cassettes"][0]["subtype"] is None
        write_annotation_exports(root, {"casandra": cas})
        from crispr_web.publication.module_result import module_result
        compact = module_result(root, "cas", "completed")["data"]
        assert compact["genes"][0]["score_margin"] == 0.0
        assert compact["cassettes"][0]["confidence"] == 0.0
        assert compact["cassettes"][0]["subtype"] is None


@pytest.mark.parametrize("field,value", [("identity", -0.1), ("identity", 1.1), ("coverage", True), ("coverage", "0.5"), ("energy_kcal_mol", "-7.5")])
def test_tracr_optional_evidence_rejected_independently(annotation_output, field, value):
    from crispr_integration.annotations import validate_tracr
    root, inputs, record, _, flush, _ = annotation_output
    key = "interaction" if field == "energy_kcal_mol" else "anti_repeat"
    record["predictions"][0][key] = {field: value}
    flush()
    assert validated(annotation_output) is None
    with pytest.raises(ValueError):
        validate_tracr(root, root / "tracrrna", inputs, mode="model", model_type="V")


@pytest.mark.parametrize("mutation", ["sequence", "duplicate", "coordinates", "experimental"])
def test_tracr_adapter_source_policy_binding(annotation_output, mutation):
    from crispr_integration.annotations import validate_tracr
    root, inputs, record, report, flush, _ = annotation_output
    assert validate_tracr(root, root / "tracrrna", inputs, mode="model", model_type="V")["prediction_count"] == 1
    if mutation == "sequence":
        record["predictions"][0]["sequence"] = "CGT"
    elif mutation == "duplicate":
        record["predictions"].append(dict(record["predictions"][0]))
        report["summary"]["prediction_count"] = 2
    elif mutation == "coordinates":
        report["coordinate_system"]["indexing"] = 0
    else:
        report["configuration"]["short_anti_repeat_rescue"] = True
    flush()
    with pytest.raises(ValueError):
        validate_tracr(root, root / "tracrrna", inputs, mode="model", model_type="V")


@pytest.mark.parametrize("energy", [None, 0, -7.5])
def test_tracr_native_energy_and_model_labels_preserved(annotation_output, energy):
    root, _, record, _, flush, _ = annotation_output
    record["predictions"][0]["interaction"] = {"energy_kcal_mol": energy}
    record["predictions"][0]["candidate_origin"] = None
    flush()
    assert validated(annotation_output) == "completed"
    cas, tracr, _ = annotation_summary(root, max_rows=10, accepted_categories=[], warnings=[])
    row = tracr["candidates"][0]
    assert row["origin"] is None
    assert row["evidence_class"] == "cm_only"
    assert row["boundary_kind"] == "merged_model_region"
    assert row["interaction_energy"] == energy
    write_annotation_exports(root, {"casandra": cas, "tracrrna": tracr})
    assert json.loads((root / "annotations.json").read_text())["tracrrna"]["candidates"][0]["interaction_energy"] == energy
    import csv
    with (root / "tracr-candidates.tsv").open() as table:
        assert next(csv.DictReader(table, delimiter="\t"))["interaction_energy"] == ("" if energy is None else str(energy))
    from crispr_web.publication.module_result import module_result
    assert module_result(root, "tracrrna", "completed")["data"]["candidates"][0]["interaction_energy"] == energy


@pytest.mark.parametrize("suffix", ["x", "\n", "\n>injected", ".extra"])
def test_tracr_export_identifiers_have_exact_native_shape(annotation_output, suffix):
    from crispr_integration.annotations import validate_tracr
    root, inputs, record, _, flush, _ = annotation_output
    record["predictions"][0]["export_id"] += suffix
    flush()
    assert validated(annotation_output) is None
    with pytest.raises(ValueError):
        validate_tracr(root, root / "tracrrna", inputs, mode="model", model_type="V")


def test_tracr_source_directory_ignores_subdirectories(annotation_output):
    from crispr_integration.annotations import validate_tracr
    root, inputs, _, _, _, _ = annotation_output
    (inputs / "unused-directory").mkdir()
    assert validate_tracr(root, root / "tracrrna", inputs, mode="model", model_type="V")["prediction_count"] == 1
