"""Run ISEScan in an isolated worker and validate native evidence against DNA."""

from pathlib import Path
import json
from .protospacer import _records, _payload, _json
from .insertion_evidence import document

MODELS = {
    "clusters.faa.hmm": "f5de1bf92059e39f66490340f1d4e585df3c4b63c79e92c74cf387d23b0007db",
    "clusters.single.faa": "da7b2b67f57394cca9c8952c21cfbf5077aaf3aabf38cefd285f0d702e234295",
}


def counts(report):
    elements = report["elements"]
    return {
        "sources": len(report["sources"]),
        "elements": len(elements),
        "complete": sum(e["native_type"] == "c" for e in elements),
        "partial": sum(e["native_type"] == "p" for e in elements),
        "with_terminal_inverted_repeats": sum(e["terminal_inverted_repeats"] is not None for e in elements),
    }


def validate_report(report, payload, native, cpu):
    computed = document(payload, native, {"threads": cpu, "include_partial": True}, MODELS)
    if json.dumps(report, sort_keys=True, allow_nan=False) != json.dumps(computed, sort_keys=True, allow_nan=False):
        raise ValueError("ISEScan source, model, interval or native evidence mismatch")
    return counts(computed)


def run_insertion_context_pipeline(input_path, output_path, *, runner_config, cpu=1, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    records = _records(Path(input_path).absolute())
    if not 1 <= len(records) <= 200 or sum(map(len, records.values())) > 2_000_000:
        raise ValueError("ISEScan source count or length exceeds budget")
    if type(cpu) is not int or cpu < 1:
        raise ValueError("ISEScan CPU budget must be a positive integer")
    # The worker budget can exceed an individual native engine's thread limit.
    native_threads = min(cpu, 8)
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "insertion-input").mkdir()
    payload = "".join(f">{key}\n{seq}\n" for key, seq in records.items()).encode("ascii")
    source = job / "insertion-input/source.fasta"
    source.write_bytes(payload)
    result = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "workflow": "insertion_context",
        "status": "running",
        "stages": {},
        "tool_options": runner_config.tool_options,
    }
    _atomic_json(job / "integration_result.json", result)
    try:
        from .context_annotations import run_cas_context
        result["stages"]["casandra"] = run_cas_context(job, records, runner_config, timeout_seconds)
        stage = _run_process(
            "isescan",
            [
                *runner_config.isescan,
                "--input",
                str(source),
                "--output",
                str(job / "isescan"),
                "--threads",
                str(native_threads),
                "--timeout-seconds",
                str(max(1, min(86400, int(timeout_seconds)))),
            ],
            job / "logs",
            timeout_seconds=timeout_seconds,
        )
        report = _json(_payload(job, job / "isescan/isescan.json"))
        stage.update(counts=validate_report(report, payload, _payload(job, job / "isescan/native.tsv"), native_threads))
        result["stages"]["isescan"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("isescan", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
