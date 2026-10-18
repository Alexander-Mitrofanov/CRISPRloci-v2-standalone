"""Offline CRISPRmap reference mapping in an isolated, bounded worker process."""

from pathlib import Path
import hashlib
import re

from .protospacer import _records, _payload, _json
from .scientific_options import tool_arguments


def validate_report(report, records, reference_sha256, input_sha256, *, max_distance=3, top_k=20):
    """Check native completion, source identities, bounded hits and label semantics."""
    if (report.get("schema") != {"name": "crisprmap.search", "version": "1.0.0"}
            or report.get("execution_completed") is not True
            or report["reference"]["file_sha256"] != reference_sha256
            or report["input"] != {"sha256": input_sha256, "record_count": len(records)}
            or report["policy"]["max_distance"] != max_distance or report["policy"]["top_k"] != top_k
            or report["policy"]["max_comparisons"] != 100_000
            or len(report["results"]) != len(records)):
        raise ValueError("Repeat mapping completion, policy or identity mismatch")
    for row, (key, sequence) in zip(report["results"], records.items(), strict=True):
        if (row["id"] != key or row["sequence"] != sequence
                or row["sequence_sha256"] != hashlib.sha256(sequence.encode()).hexdigest()
                or row["status"] not in {"exact", "near", "no_hit", "unsupported"}
                or len(row["hits"]) > top_k or row["reported_hit_count"] != len(row["hits"])):
            raise ValueError("Repeat query or hit inventory mismatch")
        for annotation in row["annotations"].values():
            if annotation["value"] is not None and (row["status"] != "exact"
                    or annotation["status"] != "historical_association"
                    or annotation["values"] != [annotation["value"]]
                    or annotation["missing_reference_count"] != 0):
                raise ValueError("Unsupported repeat annotation assignment")
        for hit in row["hits"]:
            a, b = hit["query_aligned"], hit["reference_aligned"]
            if (a.replace("-", "") != sequence or len(a) != len(b)
                    or sum(x != y for x, y in zip(a, b, strict=True)) != hit["distance"]
                    or not 0 <= hit["distance"] <= max_distance):
                raise ValueError("Invalid repeat alignment")


def run_repeat_map_pipeline(input_path, output_path, *, runner_config,
                            reference_path, reference_sha256, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    if reference_path is None or not isinstance(reference_sha256, str) or not re.fullmatch(r"[a-f0-9]{64}", reference_sha256):
        raise ValueError("Repeat mapping requires a pinned operator reference")
    reference = Path(reference_path).absolute()
    if any(p.is_symlink() for p in (reference, *reference.parents)):
        raise ValueError("Unsafe repeat reference")
    raw = _payload(reference.parent, reference)
    if hashlib.sha256(raw).hexdigest() != reference_sha256:
        raise ValueError("Repeat reference checksum changed")
    records = _records(Path(input_path).absolute())
    if len(records) > 100 or any(len(s) > 200 for s in records.values()):
        raise ValueError("Supply at most 100 repeats, each at most 200 bases")
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "repeat-map-input").mkdir()
    fasta = job / "repeat-map-input/repeats.fasta"
    fasta.write_text("".join(f">{k}\n{v}\n" for k, v in records.items()), encoding="ascii")
    native = job / "repeat-map"
    result = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
              "workflow": "repeat_map", "status": "running", "stages": {},
              "tool_options": runner_config.tool_options}
    _atomic_json(job / "integration_result.json", result)
    try:
        stage = _run_process("crisprmap", [*runner_config.crisprmap, "classify",
            "--input", str(fasta), "--reference", str(reference), "--output", str(native),
            "--max-distance", "3", "--top-k", "20", "--max-comparisons", "100000",
            "--deadline-seconds", str(max(1, min(300, timeout_seconds))),
            *tool_arguments(runner_config, "crisprmap")],
            job / "logs", timeout_seconds=timeout_seconds)
        report = _json(_payload(job, native / "report.json"))
        options = runner_config.tool_options.get("crisprmap", {})
        validate_report(report, records, reference_sha256, hashlib.sha256(fasta.read_bytes()).hexdigest(),
                        max_distance=options.get("--max-distance", 3), top_k=options.get("--top-k", 20))
        stage.update(counts=report["counts"], reference_sha256=reference_sha256)
        result["stages"]["crisprmap"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("crisprmap", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
