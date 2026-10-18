"""Invoke standalone CRISPRrepeat; no scientific package imports in Integration."""

from __future__ import annotations

from pathlib import Path

from .annotations import read_object
from .repeat_validation import validate_repeat_report
from .scientific_options import tool_arguments


def run_repeat_pipeline(
    input_path,
    output_path,
    *,
    runner_config,
    context=False,
    molecule="DNA",
    categories=None,
    cpu=None,
    timeout_seconds=3600,
):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process, run_pipeline

    source = Path(input_path).resolve(strict=True)
    categories = categories or ["Bona-fide", "Possible"]
    if molecule not in {"DNA", "RNA"} or (context and molecule != "DNA"):
        raise ValueError("Unsupported repeat input molecule")
    if context:
        result = run_pipeline(
            source,
            output_path,
            runner_config=runner_config,
            tools=(),
            categories=categories,
            cpu=cpu,
            timeout_seconds=timeout_seconds,
        )
        job = Path(output_path).resolve(strict=True)
    else:
        job = _empty_job_root(output_path)
        (job / "logs").mkdir()
        result = {
            "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
            "status": "running",
            "stages": {},
        }
    result["workflow"] = "repeat_context" if context else "repeats"
    result["tool_options"] = runner_config.tool_options
    _atomic_json(job / "integration_result.json", result)
    native = job / "repeat"
    command = [
        *runner_config.crisprrepeat,
        "--input",
        str(source),
        "--output",
        str(native),
        "--molecule",
        molecule,
        "--orientation",
        "forward" if molecule == "RNA" else "both",
        "--timeout",
        str(max(1, min(3600, int(timeout_seconds)))),
        "--categories",
        *categories,
    ]
    if context:
        command += ["--reports", str(job / "crispridentify")]
    command += tool_arguments(runner_config, "crisprrepeat")
    try:
        stage = _run_process("crisprrepeat", command, job / "logs", timeout_seconds=timeout_seconds)
        report = read_object(job, native / "repeat-analysis.json")
        validate_repeat_report(report, source, job / "crispridentify" if context else None, molecule, categories,
                               tool_options=runner_config.tool_options.get("crisprrepeat", {}))
        stage.update(report["counts"])
        result["stages"]["crisprrepeat"] = stage
        result["status"] = "completed"
        _atomic_json(job / "integration_result.json", result)
        return result
    except Exception as error:
        result.update(status="failed", error={"type": type(error).__name__, "message": str(error)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("crisprrepeat", str(error)) from error
