"""CRISPRleader's context-only contract, isolated from its Python environment."""

from __future__ import annotations

import hashlib
from pathlib import Path

from .annotations import read_object
from .parser import load_fasta_source_bindings, load_identified_arrays


def add_leader_context(job, source, result, config, timeout, flank_length=500):
    from .runner import _atomic_json, _run_process

    if type(flank_length) is not int or not 1 <= flank_length <= 5000:
        raise ValueError("Leader context length must be between 1 and 5000 bases")
    native = job / "leader"
    native.mkdir()
    records = []
    if result.get("workflow") == "loci":
        report = read_object(job, job / "tracrrna/report.json")
        for record in report["records"]:
            work = job / "tracrrna/work" / record["key"]
            records.append((record["source"]["id"], work / "input.fasta", work / "identify"))
    else:
        for fasta in sorted(source.iterdir()) if source.is_dir() else [source]:
            bindings = load_fasta_source_bindings(fasta)
            if len(bindings) != 1:
                raise ValueError("Leader server workflow requires one normalized record per file")
            records.append((bindings[0].source_id, fasta, job / "crispridentify"))
    index = {
        "schema_version": "1.0.0",
        "analysis_mode": "context_extraction",
        "flank_length": flank_length,
        "records": [],
        "array_count": 0,
        "context_count": 0,
    }
    for number, (source_id, fasta, identify) in enumerate(records, 1):
        key = f"record-{number:06d}"
        bound = load_fasta_source_bindings(fasta)[0]
        reports = [
            path
            for path in sorted(identify.rglob("report.json"))
            if read_object(job, path).get("source", {}).get("id") == bound.source_id
        ]
        if len(reports) != 1:
            raise ValueError("Leader requires exactly one canonical Identify report per source")
        load_identified_arrays(reports[0], expected_fasta_path=fasta)
        arrays = native / f"{key}.arrays.json"
        output = native / key
        commands = [
            (
                "adapt",
                ["adapt", "identify", "--fasta", str(fasta), "--reports", str(reports[0]), "--output", str(arrays)],
            ),
            (
                "extract",
                [
                    "extract",
                    "--fasta",
                    str(fasta),
                    "--arrays",
                    str(arrays),
                    "--flank-length",
                    str(flank_length),
                    "--output",
                    str(output),
                ],
            ),
            ("validate", ["validate", "--fasta", str(fasta), "--result", str(output / "result.json")]),
        ]
        for stage, args in commands:
            _run_process(f"leader-{key}-{stage}", [*config.crisprleader, *args], job / "logs", timeout_seconds=timeout)
        report = read_object(job, output / "result.json")
        manifest = read_object(job, output / "manifest.json")
        payload = (output / "result.json").read_bytes()
        if manifest.get("schema_version") != "1.0.0" or manifest["files"]["result.json"] != {
            "sha256": hashlib.sha256(payload).hexdigest(),
            "size_bytes": len(payload),
        }:
            raise ValueError("Leader result checksum mismatch")
        if (
            report.get("schema_version") != "1.0.0"
            or report.get("execution_status") != "completed"
            or report.get("analysis_mode") != "context_extraction"
            or report.get("predictions") != []
            or report.get("reference") is not None
            or report["effective_config"]["flank_length"] != flank_length
            or len(report["candidates"]) != 2 * len(report["arrays"])
        ):
            raise ValueError("Unsupported CRISPRleader context contract")
        index["records"].append({"key": key, "source_id": source_id, "native_source_id": bound.source_id})
        index["array_count"] += len(report["arrays"])
        index["context_count"] += len(report["candidates"])
    _atomic_json(native / "index.json", index)
    result["stages"]["crisprleader"] = {
        "status": "completed",
        "analysis_mode": "context_extraction",
        "prediction_performed": False,
        "record_count": len(records),
        "array_count": index["array_count"],
        "context_count": index["context_count"],
        "flank_length": flank_length,
    }


def run_leader_pipeline(input_path, output_path, *, runner_config, cpu=None, timeout_seconds=3600, flank_length=500):
    from .runner import StageError, _atomic_json, run_pipeline

    if type(flank_length) is not int or not 1 <= flank_length <= 5000:
        raise ValueError("Leader context length must be between 1 and 5000 bases")
    source = Path(input_path).resolve(strict=True)
    result = run_pipeline(
        source, output_path, runner_config=runner_config, tools=(), cpu=cpu, timeout_seconds=timeout_seconds
    )
    result["workflow"] = "leader"
    job = Path(output_path).resolve(strict=True)
    try:
        add_leader_context(job, source, result, runner_config, timeout_seconds, flank_length)
    except Exception as error:
        result.update(status="failed", error={"type": type(error).__name__, "message": str(error)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("crisprleader", str(error)) from error
    _atomic_json(job / "integration_result.json", result)
    return result
