"""Isolated RepeatTyper execution and independent adapter output validation."""
from pathlib import Path
import hashlib
import math
import json

from .protospacer import _records, _payload, _json
from .scientific_options import tool_arguments

MODEL = {'id': 'crisprcastyper-1.9.0-repeat-model',
 'upstream_repository': 'https://github.com/Russel88/CRISPRCasTyper',
 'upstream_version': '1.9.0',
 'upstream_commit': 'a2d51831d4802533bb113bac23cf9df2853063b0',
 'upstream_wheel_sha256': '9d1113ef7f6e375253986bc60740942a404b0439e1333da6ed2520ac0224cca7',
 'files': {'xgb_repeats.json': 'f9170b07805d1bb89fa1c02cfb861836e42498738e458f5e9e6db181dfcbe57a',
           'type_dict.tab': 'b01468938d868e038a14c7a42f3187bd234e2d63d49dfc4c04f92094cb4112c8'},
 'labels': ['I-A',
            'I-B',
            'I-C',
            'I-D',
            'I-E',
            'I-F',
            'I-F_T',
            'I-G',
            'II-A',
            'II-B',
            'II-C',
            'III-A',
            'III-B',
            'III-C',
            'III-D',
            'III-E',
            'III-F',
            'IV-A1',
            'IV-A2',
            'IV-A3',
            'IV-D',
            'IV-E',
            'V-A',
            'V-B1',
            'V-E',
            'V-F1',
            'V-F2',
            'V-F3',
            'V-G',
            'V-I',
            'V-J',
            'V-K',
            'VI-A',
            'VI-B1',
            'VI-B2',
            'VI-C',
            'VI-D'],
 'license': 'MIT',
 'feature_contract': 'canonical-4mer-counts-GC-fraction-length',
 'iteration_range': [0, 41]}


def _same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(right, sort_keys=True, allow_nan=False)


def validate_report(report, records, input_sha256, *, threshold=0.75):
    threshold = float(threshold)
    if (report.get("schema") != {"name": "repeattyper.classification", "version": "1.0.0"}
            or report.get("execution_completed") is not True
            or report["tool"] not in ({"name": "RepeatTyper", "version": "2.0.0"},
                                      {"name": "RepeatTyper", "version": "2.0.1"})
            or not _same(report["model"], MODEL)
            or not _same(report["input"], {"sha256": input_sha256, "record_count": len(records)})
            or not _same(report["policy"], {"threshold": threshold, "kmer": 4, "iteration_range": [0, 41],
                "score_semantics": "closed_set_model_score_not_calibrated_confidence",
                "orientation": "canonical_kmers_no_strand_inference"})
            or len(report["results"]) != len(records)):
        raise ValueError("RepeatTyper completion, policy or model identity mismatch")
    cache, canonical = {}, {}
    for row, (key, sequence) in zip(report["results"], records.items(), strict=True):
        score, winner = row["native_score"], row["winner"]
        unsupported = bool(set(sequence) - set("ACGT"))
        if (row["id"] != key or row["sequence"] != sequence
                or row["sequence_sha256"] != hashlib.sha256(sequence.encode()).hexdigest()
                or (unsupported and (score is not None or winner is not None))
                or (not unsupported and (type(score) not in (int, float) or not math.isfinite(score)
                    or not 0 <= score <= 1 or winner not in MODEL["labels"]))):
            raise ValueError("RepeatTyper occurrence or prediction mismatch")
        status = "unsupported" if unsupported else "predicted" if score >= threshold else "uncertain"
        if (row["status"] != status or row["prediction"] != (winner if status == "predicted" else None)
                or row["reason"] != ("ambiguous_iupac" if unsupported else "below_threshold" if status == "uncertain" else None)):
            raise ValueError("RepeatTyper uncertainty mismatch")
        if sequence in cache and cache[sequence] != (winner, score):
            raise ValueError("Repeated sequence has inconsistent evidence")
        cache[sequence] = winner, score
        if not unsupported:
            identity = min(sequence, sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1])
            if identity in canonical and canonical[identity] != (winner, score):
                raise ValueError("Reverse-complement equivalent sequence has inconsistent evidence")
            canonical[identity] = winner, score
    counts = {"records": len(records), "distinct_classified_sequences": sum(set(s) <= set("ACGT") for s in cache),
              **{status: sum(row["status"] == status for row in report["results"])
                 for status in ("predicted", "uncertain", "unsupported")}}
    if any(type(value) is not int for value in report["counts"].values()) or report["counts"] != counts:
        raise ValueError("RepeatTyper count mismatch")


def run_repeat_type_pipeline(input_path, output_path, *, runner_config, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    records = _records(Path(input_path).absolute())
    if not 1 <= len(records) <= 1000 or any(not 1 <= len(s) <= 200 for s in records.values()):
        raise ValueError("Supply at most 1000 repeats, each at most 200 bases")
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    (job / "repeat-type-input").mkdir()
    fasta = job / "repeat-type-input/repeats.fasta"
    fasta.write_text("".join(f">{key}\n{sequence}\n" for key, sequence in records.items()), encoding="ascii")
    result = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
              "workflow": "repeat_type", "status": "running", "stages": {}, "tool_options": runner_config.tool_options}
    _atomic_json(job / "integration_result.json", result)
    try:
        stage = _run_process("repeattyper", [*runner_config.repeattyper, "classify", "--input", str(fasta),
            "--output", str(job / "repeat-type"), *tool_arguments(runner_config, "repeattyper")],
            job / "logs", timeout_seconds=timeout_seconds)
        report = _json(_payload(job, job / "repeat-type/report.json"))
        validate_report(report, records, hashlib.sha256(fasta.read_bytes()).hexdigest(),
                        threshold=runner_config.tool_options.get("repeattyper", {}).get("--threshold", 0.75))
        stage.update(counts=report["counts"], model_id=MODEL["id"], model_sha256=MODEL["files"]["xgb_repeats.json"])
        result["stages"]["repeattyper"] = stage
        result["status"] = "completed"
    except Exception as exc:
        result.update(status="failed", error={"type": type(exc).__name__, "message": str(exc)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("repeattyper", str(exc)) from exc
    _atomic_json(job / "integration_result.json", result)
    return result
