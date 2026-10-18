"""Identify → source spacer search → pinned Phigaro genomic context."""

from pathlib import Path

from .protospacer import _records, _payload, _json
from .prophage_evidence import validate_report
from .source_spacer_evidence import context_document


def run_prophage_context_pipeline(input_path, output_path, *, runner_config, cpu=1, timeout_seconds=3600):
    from .runner import run_pipeline, _run_process, _atomic_json, StageError

    records = _records(Path(input_path).absolute())
    if not 1 <= len(records) <= 20 or sum(map(len, records.values())) > 2_000_000:
        raise ValueError("Supply at most 20 sequences and 2 Mb")
    result = run_pipeline(
        input_path, output_path, runner_config=runner_config, tools=(), cpu=cpu, timeout_seconds=timeout_seconds
    )
    job = Path(output_path).absolute()
    result.update(workflow="prophage_context", status="running")
    (job / "phigaro-input").mkdir()
    payload = "".join(f">{key}\n{seq}\n" for key, seq in records.items()).encode("ascii")
    fasta = job / "phigaro-input/source.fasta"
    fasta.write_bytes(payload)
    try:
        stage = _run_process(
            "phigaro",
            [
                *runner_config.phigaro,
                "--input",
                str(fasta),
                "--output",
                str(job / "phigaro"),
                "--threads",
                str(min(cpu, 8)),
                "--timeout-seconds",
                str(max(1, min(14400, int(timeout_seconds)))),
            ],
            job / "logs",
            timeout_seconds=timeout_seconds,
        )
        report = _json(_payload(job, job / "phigaro/report.json"))
        proteins = (
            _payload(job, job / "phigaro/work/proteins.faa") if any(len(s) >= 20000 for s in records.values()) else b""
        )
        hits = _payload(job, job / "phigaro/work/hits.tsv") if report["genes"] else b""
        stage["counts"] = validate_report(report, records, payload, proteins, hits)
        result["stages"]["phigaro"] = stage
        paths = sorted((job / "crispridentify").rglob("report.json"))
        search = _run_process(
            "source-spacer",
            [
                *runner_config.crisprspacer,
                "self",
                "--source",
                str(fasta),
                "--identify",
                *map(str, paths),
                "--scope",
                "same-record",
                "--substitutions",
                "0",
                "--min-length",
                "16",
                "--max-hits",
                "10000",
                "--max-evidence-bytes",
                "20000000",
                "--timeout",
                str(max(1, min(86400, int(timeout_seconds)))),
                "--output",
                str(job / "source-spacer"),
            ],
            job / "logs",
            timeout_seconds=timeout_seconds,
        )
        native = _json(_payload(job, job / "source-spacer/analysis.json"))
        context = context_document(native, [_json(_payload(job, p)) for p in paths], records, report["regions"])
        search["counts"] = context["search_counts"]
        _atomic_json(job / "phigaro/context.json", context)
        result["stages"]["source_spacer"] = search
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("phigaro", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
