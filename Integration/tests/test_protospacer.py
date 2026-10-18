"""The indexed search bypasses detection and binds tool output to source DNA."""

import hashlib
import json
from pathlib import Path
import sqlite3

import pytest

from crispr_integration.protospacer import run_protospacer_pipeline, validate_protospacer_report
from crispr_integration.runner import RunnerConfig, StageError

FIXTURE = Path(__file__).resolve().parents[2] / "backend/tests/fixtures/protospacer"


def test_operator_command_and_verified_manifest(tmp_path, monkeypatch):
    from crispr_integration import runner
    database = tmp_path / "database"
    database.mkdir()
    report = json.loads((FIXTURE / "analysis.json").read_text())
    retained = json.loads((FIXTURE / "inputs.json").read_text())
    with sqlite3.connect(database / "spacers.sqlite3") as connection:
        connection.executescript("CREATE TABLE metadata(value TEXT); CREATE TABLE spacers(id TEXT PRIMARY KEY, sequence TEXT, status TEXT);")
        connection.execute("INSERT INTO metadata VALUES(?)", (json.dumps({k: v for k, v in report["reference"].items() if k != "database_sha256"}),))
        connection.executemany("INSERT INTO spacers VALUES(?,?,?)", [(row["spacer_id"], row["spacer_sequence"], "eligible") for row in report["matches"]])
    digest = hashlib.sha256((database / "spacers.sqlite3").read_bytes()).hexdigest()
    report["reference"]["database_sha256"] = digest
    retained["reference"] = report["reference"]
    (database / "reference.json").write_text(json.dumps(report["reference"]))
    calls = []

    def execute(stage, command, logs, *, timeout_seconds):
        calls.append(command)
        assert stage == "crisprspacer" and command[:2] == ["/isolated/crisprspacer", "database-search"]
        assert command[command.index("--expected-database-sha256") + 1] == digest
        assert timeout_seconds == 60
        native = Path(command[command.index("--output") + 1])
        native.mkdir()
        (native / "analysis.json").write_text(json.dumps(report))
        (native / "inputs.json").write_text(json.dumps(retained))
        (native / "matches.tsv").write_text("native table\n")
        files = {}
        for name in ("analysis.json", "inputs.json", "matches.tsv"):
            data = (native / name).read_bytes()
            files[name] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
        (native / "manifest.json").write_text(json.dumps({"schema_version": "1.0.0", "kind": "crisprspacer.database-search", "files": files}))
        return {"status": "completed"}

    monkeypatch.setattr(runner, "_run_process", execute)
    config = RunnerConfig(("must-not-detect",), ("unused",), ("unused",), crisprspacer=("/isolated/crisprspacer",))
    result = run_protospacer_pipeline(FIXTURE / "target.fasta", tmp_path / "output", runner_config=config,
                                    database_path=database, database_sha256=digest, timeout_seconds=60)
    assert len(calls) == 1 and result["status"] == "completed"
    assert set(result["stages"]) == {"crisprspacer"}


@pytest.mark.parametrize("field,value", [("relative_strand", "-"), ("target_end", 17), ("substitutions", 1), ("pam_status", "compatible")])
def test_source_reconstruction_rejects_false_matches(field, value):
    report = json.loads((FIXTURE / "analysis.json").read_text())
    report["matches"][0][field] = value
    with pytest.raises(ValueError):
        validate_protospacer_report(report, {"target": "AAAACCCCGGGGTTTA"}, "a" * 64)


def test_large_target_is_wrapped_before_the_native_cli(tmp_path, monkeypatch):
    from crispr_integration import runner
    source = tmp_path / "large.fasta"
    sequence = "ACGT" * 250_001
    source.write_text(">large original description\n" + sequence + "\n")
    database = tmp_path / "database"
    database.mkdir()

    def inspect(stage, command, logs, *, timeout_seconds):
        target = Path(command[command.index("--targets") + 1])
        lines = target.read_text().splitlines()
        assert lines[0] == ">large"
        assert max(map(len, lines[1:])) == 80
        assert "".join(lines[1:]) == sequence
        assert command[command.index("--timeout") + 1] == "300"
        raise RuntimeError("inspection complete")

    monkeypatch.setattr(runner, "_run_process", inspect)
    with pytest.raises(StageError, match="inspection complete"):
        run_protospacer_pipeline(source, tmp_path / "out", runner_config=RunnerConfig(("unused",), ("unused",), ("unused",)),
                                database_path=database, database_sha256="a" * 64, timeout_seconds=7200)
