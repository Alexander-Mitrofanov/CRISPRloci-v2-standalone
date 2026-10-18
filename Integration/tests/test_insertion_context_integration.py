import csv
import io
import json
from pathlib import Path
import pytest
from crispr_integration.insertion_context import run_insertion_context_pipeline, validate_report, MODELS
from crispr_integration.insertion_evidence import document, FIELDS
from crispr_integration.runner import RunnerConfig

PAYLOAD = b">f000001.s000001\nAAAACCCCTTTT\n"


def table():
    values = [
        "seq000001",
        "IS1",
        "IS1_1",
        "1",
        "12",
        "12",
        "2",
        "1",
        "4",
        "9",
        "12",
        "8",
        "4",
        "4",
        "0",
        "0",
        "0",
        "",
        "1",
        "1e-20",
        "-1",
        "p",
        "0",
        "AAAA:AAAA",
    ]
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(FIELDS)
    writer.writerow(values)
    return buffer.getvalue().encode()


def test_source_bound_native_report():
    native = table()
    report = document(PAYLOAD, native, {"threads": 1, "include_partial": True}, MODELS)
    assert validate_report(report, PAYLOAD, native, 1) == {
        "sources": 1,
        "elements": 1,
        "complete": 0,
        "partial": 1,
        "with_terminal_inverted_repeats": 1,
    }
    assert report["elements"][0]["evalue_this_copy"] is None


@pytest.mark.parametrize("change", ["coordinate", "evalue", "tsd", "source", "model", "extra"])
def test_changed_evidence_is_rejected(change):
    native = table()
    report = document(PAYLOAD, native, {"threads": 1, "include_partial": True}, MODELS)
    if change == "coordinate":
        report["elements"][0]["end"] = 11
    if change == "evalue":
        report["elements"][0]["evalue_this_copy"] = 0
    if change == "tsd":
        report["elements"][0]["target_site_duplication"] = "AAAA"
    if change == "source":
        report["sources"][0]["sha256"] = "0" * 64
    if change == "model":
        report["models"] = {}
    if change == "extra":
        report["private_path"] = "/private"
    with pytest.raises(ValueError):
        validate_report(report, PAYLOAD, native, 1)


@pytest.mark.parametrize("worker_cpu,native_threads", [(1, 1), (8, 8), (14, 8), (256, 8)])
def test_worker_budget_is_capped_and_bound_to_native_evidence(
    tmp_path, monkeypatch, worker_cpu, native_threads
):
    source = tmp_path / "input.fasta"
    source.write_bytes(PAYLOAD)
    config = RunnerConfig(("identify",), ("spacerplacer",), ("evor",), isescan=("isescan-v2",))
    monkeypatch.setattr(
        "crispr_integration.context_annotations.run_cas_context",
        lambda *_: {"status": "completed", "gene_count": 0, "cas_gene_count": 0, "cassette_count": 0},
    )
    calls = []

    def native_run(stage, command, logs, *, timeout_seconds):
        assert stage == "isescan"
        threads = int(command[command.index("--threads") + 1])
        calls.append(threads)
        assert threads == native_threads
        output = Path(command[command.index("--output") + 1])
        output.mkdir()
        native = table()
        report = document(PAYLOAD, native, {"threads": threads, "include_partial": True}, MODELS)
        (output / "native.tsv").write_bytes(native)
        (output / "isescan.json").write_text(json.dumps(report))
        return {"status": "completed"}

    monkeypatch.setattr("crispr_integration.runner._run_process", native_run)
    result = run_insertion_context_pipeline(source, tmp_path / "output", runner_config=config, cpu=worker_cpu)
    assert calls == [native_threads]
    assert result["status"] == "completed"
    assert result["stages"]["isescan"]["counts"]["elements"] == 1


@pytest.mark.parametrize("cpu", [0, -1, True, 1.5, "14", None])
def test_invalid_worker_budget_fails_before_job_or_process_creation(tmp_path, monkeypatch, cpu):
    source = tmp_path / "input.fasta"
    source.write_bytes(PAYLOAD)
    output = tmp_path / "output"
    config = RunnerConfig(("identify",), ("spacerplacer",), ("evor",))

    def unexpected_process(*args, **kwargs):
        pytest.fail("Invalid budget dispatched a scientific process")

    monkeypatch.setattr("crispr_integration.runner._run_process", unexpected_process)
    with pytest.raises(ValueError, match="positive integer"):
        run_insertion_context_pipeline(source, output, runner_config=config, cpu=cpu)
    assert not output.exists()
