"""Bounded isolated reverse spacer search with pinned reference and native replay."""

from pathlib import Path
import re

from .protospacer import _records, _payload, _json, _hash


def run_viral_pipeline(
    input_path,
    output_path,
    *,
    runner_config,
    reference_path,
    reference_sha256,
    max_mismatches=2,
    threads=2,
    timeout_seconds=3600,
):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    if (
        reference_path is None
        or not isinstance(reference_sha256, str)
        or not re.fullmatch("[a-f0-9]{64}", reference_sha256)
    ):
        raise ValueError("Viral search requires a pinned operator reference")
    reference = Path(reference_path).absolute()
    if not reference.is_dir() or any(p.is_symlink() for p in (reference, *reference.parents)):
        raise ValueError("Invalid viral reference directory")
    if len(runner_config.blastn) != 1 or type(max_mismatches) is not int or max_mismatches not in range(3):
        raise ValueError("Invalid viral-search runtime or mismatch policy")
    records = _records(Path(input_path).absolute())
    if any(len(s) > 80 for s in records.values()):
        raise ValueError("Supply spacer sequences at most 80 bases long")
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "viral-input").mkdir()
    fasta = job / "viral-input/spacers.fasta"
    fasta.write_text("".join(f">{key}\n{value}\n" for key, value in records.items()), encoding="ascii")
    native = job / "viral"
    result = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "status": "running",
        "workflow": "viral_search",
        "stages": {},
    }
    _atomic_json(job / "integration_result.json", result)
    try:
        command = [
            *runner_config.crisprspacer,
            "viral-search",
            "--spacers",
            str(fasta),
            "--reference",
            str(reference),
            "--expected-reference-sha256",
            reference_sha256,
            "--blastn",
            runner_config.blastn[0],
            "--max-mismatches",
            str(max_mismatches),
            "--threads",
            str(min(threads, 8)),
            "--timeout",
            str(max(1, min(300, int(timeout_seconds)))),
            "--output",
            str(native),
        ]
        stage = _run_process("crisprspacer", command, job / "logs", timeout_seconds=timeout_seconds)
        _run_process(
            "validate_viral",
            [
                *runner_config.crisprspacer,
                "validate-viral-result",
                str(native),
                "--reference",
                str(reference),
                "--expected-reference-sha256",
                reference_sha256,
            ],
            job / "logs",
            timeout_seconds=min(timeout_seconds, 300),
        )
        for name in ("validate_viral.stdout.log", "validate_viral.stderr.log"):
            (job / "logs" / name).unlink()
        report = _json(_payload(job, native / "result.json"))
        inputs = _json(_payload(job, native / "inputs.json"))
        schema = {"name": "crisprspacer.viral-search", "version": "1.0.0"}
        if (report.get("schema") != schema or inputs.get("schema") != schema
                or report.get("tool_version") not in {"0.3.0", "0.3.1"}
                or report.get("execution_completed") is not True
                or report.get("retrieval_completeness") != "heuristic"):
            raise ValueError("Unconfirmed viral search")
        if (
            report["reference"]["reference_sha256"] != reference_sha256
            or report["policy"]["max_substitutions"] != max_mismatches
        ):
            raise ValueError("Viral search reference/policy mismatch")
        if inputs["spacer_fasta_sha256"] != _hash(fasta.read_text()):
            raise ValueError("Viral query byte identity mismatch")
        if [(q["original_id"], q["sequence"], q["sha256"]) for q in inputs["queries"]] != [
            (k, v, _hash(v)) for k, v in records.items()
        ]:
            raise ValueError("Viral query source mismatch")
        stage.update(counts=report["counts"], reference_sha256=reference_sha256)
        result["stages"]["crisprspacer"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("crisprspacer", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
