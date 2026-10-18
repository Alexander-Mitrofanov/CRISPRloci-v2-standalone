"""Independent protospacer validation, durable reference identity and publication."""

from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sqlite3
import sys
import zipfile

import pytest
from pydantic import ValidationError

from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.models import JobSubmission
from crispr_web.publication.protospacer import write_protospacer_exports
from crispr_web.service import JobService, ReferenceUnavailableError
from crispr_web.summaries.protospacer import protospacer_summary
from crispr_web.worker import Worker

FIXTURE = Path(__file__).parent / "fixtures/protospacer"


def configured_settings(settings_factory, tmp_path, **overrides):
    names = ("crispridentify", "spacerplacer", "crispr_evor", "casandra", "crispr_tracrrna", "crisprspacer")
    artifacts = {name: {"name": name, "display_name": name, "version": "0.2.0", "kind": "git_source",
                       "path": "/private/" + name, "sha256": "a" * 64, "sha256_scope": "source_inventory_v1",
                       "repository": "https://github.com/example/" + name + ".git", "commit": "b" * 40}
                 for name in names}
    artifacts["integration"] = {"name": "integration", "display_name": "Integration", "version": "0.6.0", "kind": "wheel",
                                "path": "/private/integration.whl", "sha256": "c" * 64, "sha256_scope": "file_bytes"}
    document = {"schema": "crispr-web-scientific-release", "schema_version": "5.0.0", "release_id": "d" * 64,
                "scientific_artifacts": artifacts}
    payload = json.dumps(document).encode()
    manifest = tmp_path / "spacer-release.json"
    manifest.write_bytes(payload)
    return settings_factory(scientific_release_manifest=manifest,
                            scientific_release_manifest_sha256=hashlib.sha256(payload).hexdigest(), **overrides)


def fixture(tmp_path):
    root, input_root = tmp_path / "output", tmp_path / "input"
    native = root / "protospacer"
    native.mkdir(parents=True)
    input_root.mkdir()
    (input_root / "000001.fasta").write_bytes((FIXTURE / "target.fasta").read_bytes())
    report = json.loads((FIXTURE / "analysis.json").read_text())
    inputs = json.loads((FIXTURE / "inputs.json").read_text())
    database = tmp_path / "database"
    database.mkdir()
    with sqlite3.connect(database / "spacers.sqlite3") as connection:
        connection.executescript("CREATE TABLE metadata(value TEXT); CREATE TABLE spacers(id TEXT PRIMARY KEY, sequence TEXT, status TEXT);")
        connection.execute("INSERT INTO metadata VALUES(?)", (json.dumps({k: v for k, v in report["reference"].items() if k != "database_sha256"}),))
        connection.executemany("INSERT INTO spacers VALUES(?,?,?)", [(row["spacer_id"], row["spacer_sequence"], "eligible") for row in report["matches"]])
    digest = hashlib.sha256((database / "spacers.sqlite3").read_bytes()).hexdigest()
    report["reference"]["database_sha256"] = digest
    inputs["reference"] = dict(report["reference"])
    (database / "reference.json").write_text(json.dumps(report["reference"]))
    pipeline = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
                "status": "completed", "workflow": "protospacer", "stages": {
                    "crisprspacer": {"status": "completed", "counts": dict(report["counts"]), "database_sha256": digest}}}

    def flush():
        (native / "analysis.json").write_text(json.dumps(report))
        (native / "inputs.json").write_text(json.dumps(inputs))
        (native / "matches.tsv").write_text("native table is private\n")
        manifest = {"schema_version": "1.0.0", "kind": "crisprspacer.database-search", "files": {}}
        for name in ("analysis.json", "inputs.json", "matches.tsv"):
            payload = (native / name).read_bytes()
            manifest["files"][name] = {"sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}
        (native / "manifest.json").write_text(json.dumps(manifest))
        (root / "integration_result.json").write_text(json.dumps(pipeline))

    flush()
    claim = ClaimedJob("a" * 32, "protospacer", 1, 0, "bona_fide_possible", 1, True, True, (), spacer_database_sha256=digest)
    binding = {row["id"]: (row["length"], row["sha256"]) for row in report["targets"]}
    return root, input_root, report, inputs, pipeline, claim, binding, flush


def validated(data):
    root, inputs, _, _, _, claim, binding, _ = data
    return _validated_pipeline_status(claim, root, binding, input_root=inputs, spacer_reference=root.parent / "database")


def test_source_bound_search_exports_all_rows_independently_of_preview(tmp_path):
    data = fixture(tmp_path)
    assert validated(data) == "completed"
    root = data[0]
    preview = protospacer_summary(root, max_rows=0)
    assert preview["available"] and preview["matches_truncated"] and preview["matches"] == []
    write_protospacer_exports(root, {"protospacer": preview})
    assert len(json.loads((root / "protospacer-analysis.json").read_text())["matches"]) == 1
    assert "spacer-1" in (root / "protospacer-matches.tsv").read_text()


@pytest.mark.parametrize("mutation", ["sequence", "start", "strand", "counts", "reference", "source", "policy", "pam", "duplicate", "retained", "stage", "symlink", "checksum"])
def test_corrupt_scientific_evidence_is_rejected_even_with_rehashed_manifest(tmp_path, mutation):
    data = fixture(tmp_path)
    root, _, report, inputs, pipeline, _, _, flush = data
    row = report["matches"][0]
    if mutation == "sequence":
        row["spacer_sequence"] = "T" * 16
    elif mutation == "start":
        row["target_start"] = True
    elif mutation == "strand":
        row["relative_strand"] = "-"
    elif mutation == "counts":
        report["counts"]["seed_windows"] += 1
    elif mutation == "reference":
        report["reference"]["database_sha256"] = "c" * 64
    elif mutation == "source":
        report["targets"][0]["sha256"] = "c" * 64
    elif mutation == "policy":
        report["policy"]["substitutions"] = 1
    elif mutation == "pam":
        row["pam_status"] = "compatible"
    elif mutation == "duplicate":
        report["matches"].append(dict(row))
    elif mutation == "retained":
        inputs["targets"][0]["sequence"] = "T" * 16
    elif mutation == "stage":
        pipeline["stages"]["crisprspacer"]["database_sha256"] = "c" * 64
    flush()
    if mutation == "symlink":
        file = root / "protospacer/analysis.json"
        target = tmp_path / "external.json"
        file.rename(target)
        file.symlink_to(target)
    elif mutation == "checksum":
        (root / "protospacer/matches.tsv").write_text("changed")
    assert validated(data) is None


def test_empty_reference_is_successful_no_assessable_windows(tmp_path):
    data = list(fixture(tmp_path))
    _, _, report, inputs, pipeline, claim, _, flush = data
    report["matches"] = []
    report["counts"].update(matches=0, matched_spacers=0, candidates_checked=0)
    report["reference"].update(eligible_record_count=0, skipped_length=1)
    report["outcome"] = "no_assessable_windows"
    database = tmp_path / "database"
    with sqlite3.connect(database / "spacers.sqlite3") as connection:
        connection.execute("UPDATE spacers SET sequence='A',status='skipped_length'")
        connection.execute("UPDATE metadata SET value=?", (json.dumps({k: v for k, v in report["reference"].items() if k != "database_sha256"}),))
    digest = hashlib.sha256((database / "spacers.sqlite3").read_bytes()).hexdigest()
    report["reference"]["database_sha256"] = digest
    (database / "reference.json").write_text(json.dumps(report["reference"]))
    inputs["reference"] = dict(report["reference"])
    pipeline["stages"]["crisprspacer"].update(counts=dict(report["counts"]), database_sha256=digest)
    data[5] = replace(claim, spacer_database_sha256=digest)
    flush()
    assert validated(data) == "completed"


def test_retained_hash_binds_wrapped_native_target_bytes(tmp_path):
    data = fixture(tmp_path)
    _, input_root, report, inputs, pipeline, _, binding, flush = data
    sequence = "AAAACCCCGGGGTTTA" + "N" * 84
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    # Worker input wrapping need not equal the private native wrapping. The
    # adapter emits 80-base lines, and backend reconstructs those bytes itself.
    (input_root / "000001.fasta").write_text(">target\n" + sequence + "\n")
    report["targets"][0].update(length=100, sha256=digest)
    report["matches"][0]["target_sha256"] = digest
    report["counts"].update(target_bases=100, seed_windows=85, ambiguous_seed_windows=84)
    pipeline["stages"]["crisprspacer"]["counts"] = dict(report["counts"])
    inputs["targets"][0].update(sequence=sequence, length=100, sha256=digest)
    native_bytes = (">target\n" + sequence[:80] + "\n" + sequence[80:] + "\n").encode()
    inputs["target_file_sha256"] = hashlib.sha256(native_bytes).hexdigest()
    binding["target"] = (100, digest)
    flush()
    assert validated(data) == "completed"


def test_reference_digest_is_persisted_and_cannot_drift(settings_factory, tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    settings = configured_settings(settings_factory, tmp_path, spacer_database_path=reference, spacer_database_sha256="a" * 64)
    service = JobService(settings)
    _, job = service.submit(JobSubmission(mode="protospacer", sequence=(FIXTURE / "target.fasta").read_text()), "127.0.0.1")
    service.store.initialize()
    claim = service.store.claim_next("spacer-worker")
    assert claim.spacer_database_sha256 == "a" * 64
    assert "spacer_database_sha256" not in job["options"]
    command = Worker(settings).command_for(claim, tmp_path / "out")
    assert command[command.index("--spacer-database-sha256") + 1] == "a" * 64
    with pytest.raises(ValueError, match="changed"):
        Worker(replace(settings, spacer_database_sha256="c" * 64)).command_for(claim, tmp_path / "out")
    service.store.request_cancel(claim.job_id)
    assert service.store.get_job(claim.job_id)["cancel_requested"]


def test_unconfigured_database_and_request_paths_are_rejected(settings_factory):
    service = JobService(settings_factory())
    with pytest.raises(ReferenceUnavailableError, match="unavailable"):
        service.submit(JobSubmission(mode="protospacer", sequence=">target\nACGT\n"), "127.0.0.1")
    with pytest.raises(ValidationError):
        JobSubmission(mode="protospacer", sequence="ACGT", spacer_database="/private")


def test_database_settings_require_digest_and_no_symlink(settings_factory, tmp_path):
    reference = tmp_path / "reference"
    reference.mkdir()
    with pytest.raises(ValueError, match="together"):
        settings_factory(spacer_database_path=reference)
    link = tmp_path / "link"
    link.symlink_to(reference, target_is_directory=True)
    with pytest.raises(ValueError, match="symlinks"):
        settings_factory(spacer_database_path=link, spacer_database_sha256="a" * 64)


def test_queued_worker_publication_preserves_private_boundary(settings_factory, tmp_path):
    data = fixture(tmp_path)
    reference = tmp_path / "database"
    digest = data[5].spacer_database_sha256
    settings = configured_settings(settings_factory, tmp_path, spacer_database_path=reference, spacer_database_sha256=digest,
                                integration_command=(sys.executable, str(Path(__file__).with_name("fake_protospacer.py"))))
    service = JobService(settings)
    token, job = service.submit(JobSubmission(mode="protospacer", sequence=(FIXTURE / "target.fasta").read_text()), "127.0.0.1")
    assert Worker(settings, store=service.store, worker_id="spacer-worker").run_once()
    result = service.get(job["job_id"])
    assert result["status"] == "completed"
    assert result["summary"]["protospacer"]["counts"]["matches"] == 1
    assert {row["name"] for row in result["artifacts"]} >= {"protospacer-analysis.json", "protospacer-matches.tsv", "results.zip"}
    assert service.authorize(job["job_id"], token)
    bundle, _ = service.result_bundle(job["job_id"])
    with zipfile.ZipFile(bundle) as archive:
        assert all("/protospacer/" not in name and "inputs.json" not in name for name in archive.namelist())
        provenance = json.loads(archive.read("provenance.json"))
        compact = json.loads(archive.read("results/result.json"))
        assert compact["complete"] is True and compact["module"] == "protospacer"
        assert len(compact["data"]["matches"]) == 1
    assert provenance["options"]["spacer_database_sha256"] == digest
    assert str(reference) not in json.dumps(result)


def test_protospacer_manifest_requires_version_five_and_tool_identity():
    from crispr_web.config import _validated_scientific_manifest
    names = ("crispridentify", "spacerplacer", "crispr_evor", "casandra", "crispr_tracrrna", "crisprspacer")
    artifacts = {name: {"name": name, "display_name": name, "version": "0.2.0", "kind": "git_source",
                       "path": "/private/" + name, "sha256": "a" * 64, "sha256_scope": "source_inventory_v1",
                       "repository": "https://github.com/example/" + name + ".git", "commit": "b" * 40}
                 for name in names}
    artifacts["integration"] = {"name": "integration", "display_name": "Integration", "version": "0.6.0", "kind": "wheel",
                                "path": "/private/integration.whl", "sha256": "c" * 64, "sha256_scope": "file_bytes"}
    document = {"schema": "crispr-web-scientific-release", "schema_version": "5.0.0", "release_id": "d" * 64,
                "scientific_artifacts": artifacts}
    assert _validated_scientific_manifest(document)["crisprspacer"]["version"] == "0.2.0"
    document["schema_version"] = "4.0.0"
    with pytest.raises(ValueError, match="artifact set"):
        _validated_scientific_manifest(document)


@pytest.mark.skipif(not (os.getenv("CRISPRSPACER_TEST_COMMAND") and os.getenv("CRISPR_INTEGRATION_TEST_COMMAND")),
                    reason="Explicit installed CRISPRspacer and Integration commands are required")
def test_real_installed_protospacer_lifecycle(settings_factory, tmp_path):
    """Opt-in proof across real process boundaries; no tool-package imports."""
    spacer_command = json.loads(os.environ["CRISPRSPACER_TEST_COMMAND"])
    integration_command = json.loads(os.environ["CRISPR_INTEGRATION_TEST_COMMAND"])
    for command in (spacer_command, integration_command):
        assert isinstance(command, list) and command and all(isinstance(part, str) and part for part in command)
        assert Path(command[0]).is_absolute()
    example = Path(os.getenv("CRISPRSPACER_TEST_FASTA", str(Path(__file__).resolve().parents[2] / "tools/CRISPRspacer/examples/ncbi-spacerome-three.fasta")))
    assert example.is_file(), "Supply the current pinned checkout's example FASTA explicitly"
    database = tmp_path / "database"
    subprocess.run([*spacer_command, "database-index", "--spacers", str(example), "--name", "ncbi-spacerome-three-acceptance",
                    "--source-url", "https://ftp.ncbi.nlm.nih.gov/pub/wolf/_suppl/spacerome/Supplementary_material_spacers.fna.gz",
                    "--source-release", "2017-04-25", "--output", str(database)],
                   check=True, capture_output=True, text=True, timeout=60, cwd=tmp_path)
    reference = json.loads((database / "reference.json").read_text())
    settings = configured_settings(settings_factory, tmp_path, spacer_database_path=database,
                                    spacer_database_sha256=reference["database_sha256"],
                                    integration_command=tuple(integration_command))
    settings.runner_config.write_text(json.dumps({"crispridentify": ["must-not-run"], "spacerplacer": ["must-not-run"],
                                                  "crispr_evor": ["must-not-run"], "crisprspacer": spacer_command}))
    # Parse only the first published example record; do not import scientific helpers.
    lines = example.read_text().splitlines()
    sequence = "".join(lines[1:next((i for i in range(1, len(lines)) if lines[i].startswith(">")), len(lines))])
    reverse = "".join({"A": "T", "T": "A", "G": "C", "C": "G"}[letter] for letter in sequence[::-1])
    service = JobService(settings)
    token, created = service.submit(JobSubmission(mode="protospacer", sequence=f">forward\n{sequence}\n>reverse\n{reverse}\n"), "127.0.0.1")
    assert Worker(settings, store=service.store, worker_id="real-spacer-worker").run_once()
    job = service.get(created["job_id"])
    assert job["status"] == "completed", job.get("error")
    rows = job["summary"]["protospacer"]["matches"]
    assert any(row["target_id"] == "forward" and row["relative_strand"] == "+" for row in rows)
    assert any(row["target_id"] == "reverse" and row["relative_strand"] == "-" for row in rows)
    assert job["summary"]["protospacer"]["reference"]["database_sha256"] == reference["database_sha256"]
    assert service.authorize(created["job_id"], token)
    bundle, _ = service.result_bundle(created["job_id"])
    with zipfile.ZipFile(bundle) as archive:
        assert "results/protospacer-analysis.json" in archive.namelist()
        assert "results/protospacer-matches.tsv" in archive.namelist()
        assert all("/protospacer/" not in name and "inputs.json" not in name for name in archive.namelist())
    _, empty = service.submit(JobSubmission(mode="protospacer", sequence=">nohit\n" + "A" * 80 + "\n"), "127.0.0.1")
    assert Worker(settings, store=service.store, worker_id="real-spacer-worker").run_once()
    nohit = service.get(empty["job_id"])
    assert nohit["status"] == "completed", nohit.get("error")
    assert nohit["summary"]["protospacer"]["outcome"] == "no_matches"
    assert nohit["summary"]["protospacer"]["counts"]["matches"] == 0


@pytest.mark.parametrize("mutation", ["id", "reference_label", "sequence", "eligibility", "bytes", "internal", "oversized_row", "missing_index"])
def test_pinned_catalogue_binding_at_both_gates(tmp_path, mutation):
    from crispr_integration.protospacer import validate_protospacer_report
    data = fixture(tmp_path)
    root, inputs, report, retained, pipeline, claim, binding, flush = data
    database = tmp_path / "database"
    if mutation == "id":
        report["matches"][0]["spacer_id"] = "absent"
    elif mutation == "reference_label":
        report["reference"]["source_release"] = "forged"
        retained["reference"] = report["reference"]
    elif mutation == "bytes":
        with (database / "spacers.sqlite3").open("ab") as stream:
            stream.write(b"changed")
    else:
        with sqlite3.connect(database / "spacers.sqlite3") as connection:
            if mutation == "sequence":
                connection.execute("UPDATE spacers SET sequence=?", ("A" * 16,))
            elif mutation == "eligibility":
                connection.execute("UPDATE spacers SET status='skipped_length'")
            elif mutation == "internal":
                connection.execute("UPDATE metadata SET value='{}'")
            elif mutation == "oversized_row":
                connection.execute("UPDATE spacers SET sequence=?", ("A" * 100_000,))
            elif mutation == "missing_index":
                connection.executescript("ALTER TABLE spacers RENAME TO old; CREATE TABLE spacers(id TEXT,sequence TEXT,status TEXT); INSERT INTO spacers SELECT * FROM old; DROP TABLE old;")
        digest = hashlib.sha256((database / "spacers.sqlite3").read_bytes()).hexdigest()
        report["reference"]["database_sha256"] = digest
        retained["reference"] = report["reference"]
        pipeline["stages"]["crisprspacer"]["database_sha256"] = digest
        claim = replace(claim, spacer_database_sha256=digest)
        (database / "reference.json").write_text(json.dumps(report["reference"]))
    flush()
    assert _validated_pipeline_status(claim, root, binding, input_root=inputs, spacer_reference=database) is None
    with pytest.raises(ValueError):
        validate_protospacer_report(report, {"target": "AAAACCCCGGGGTTTA"}, claim.spacer_database_sha256, database)


def test_catalogue_validation_preserves_checkpoint_interruption(tmp_path):
    from crispr_integration.protospacer_reference import bind_catalogue
    from crispr_web.execution.validation_protospacer_reference import validate_reference
    data = fixture(tmp_path)
    report = data[2]
    class LeaseLost(RuntimeError):
        pass
    def cancelled():
        raise LeaseLost("cancelled")
    for validator in (bind_catalogue, validate_reference):
        with pytest.raises(LeaseLost, match="cancelled"):
            validator(tmp_path / "database", data[5].spacer_database_sha256, report["reference"], {}, cancelled)


@pytest.mark.parametrize("kind", ["corrupt_database", "metadata_fifo", "database_fifo"])
def test_reference_setup_failures_are_bounded_value_errors(tmp_path, kind):
    from crispr_integration.protospacer_reference import bind_catalogue
    from crispr_web.execution.validation_protospacer_reference import validate_reference
    root = tmp_path / "catalogue"
    root.mkdir()
    database = root / "spacers.sqlite3"
    database.write_bytes(b"not a SQLite database")
    digest = hashlib.sha256(database.read_bytes()).hexdigest()
    reference = {"database_sha256": digest}
    (root / "reference.json").write_text(json.dumps(reference))
    if kind.endswith("fifo"):
        path = root / ("reference.json" if kind == "metadata_fifo" else "spacers.sqlite3")
        path.unlink()
        os.mkfifo(path)
    for bind in (bind_catalogue, validate_reference):
        with pytest.raises(ValueError):
            bind(root, digest, reference, {})
