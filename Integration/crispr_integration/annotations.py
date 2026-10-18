"""Isolated CasAndra and tracrRNA workflows using their versioned file contracts."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

from .parser import _read_fasta_records, load_fasta_source_bindings, load_identified_arrays
from .scientific_options import tool_arguments

MAX_REPORT_BYTES = 20_000_000


def _payload(root: Path, path: Path) -> bytes:
    relative = path.relative_to(root)
    current = root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise ValueError("Scientific outputs must not contain symbolic links")
    if not path.is_file() or path.stat().st_size > MAX_REPORT_BYTES:
        raise ValueError(f"Missing or oversized scientific output: {relative}")
    with path.open("rb") as stream:
        payload = stream.read(MAX_REPORT_BYTES + 1)
    if len(payload) > MAX_REPORT_BYTES:
        raise ValueError("Scientific output grew past the limit")
    return payload


def _strict_json(payload):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate JSON key in scientific output")
            result[key] = value
        return result

    def number(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError("Nonfinite JSON number in scientific output")
        return result

    def invalid(value):
        raise ValueError("Invalid JSON constant: " + value)

    return json.loads(payload, object_pairs_hook=pairs, parse_float=number, parse_constant=invalid)


def read_object(root: Path, path: Path) -> dict:
    value = _strict_json(_payload(root, path))
    if not isinstance(value, dict):
        raise ValueError(f"Expected object in {path.relative_to(root)}")
    return value


def checksum(path: Path) -> str:
    return hashlib.sha256(_payload(Path(path.anchor), path)).hexdigest()


def json_lines(root: Path, path: Path) -> list[dict]:
    rows = [_strict_json(line) for line in _payload(root, path).splitlines() if line.strip()]
    if any(not isinstance(row, dict) for row in rows):
        raise ValueError("JSONL records must be objects")
    return rows


def validate_casandra(root: Path, output: Path, fasta: Path) -> dict:
    run = read_object(root, output / "run.json")
    manifest = read_object(root, output / "manifest.json")
    if run.get("program") != "CasAndra" or run.get("schema_version") != 5 or manifest.get("schema_version") != 5:
        raise ValueError("Unsupported CasAndra genome output contract")
    if not any(row.get("kind") == "genome_fasta" and row.get("sha256") == checksum(fasta) for row in run["inputs"]):
        raise ValueError("CasAndra input checksum mismatch")
    sources = {row.source_id: row for row in load_fasta_source_bindings(fasta)}
    data = {}
    for filename, key in [("proteins.jsonl", "genes"), ("cassettes.jsonl", "cas_cassettes")]:
        path = output / filename
        if manifest["files"][filename]["sha256"] != checksum(path):
            raise ValueError("CasAndra output checksum mismatch")
        rows = json_lines(root, path)
        if type(run.get(key)) is not int or len(rows) != run[key]:
            raise ValueError("CasAndra output count mismatch")
        for row in rows:
            source = sources.get(row.get("contig_id"))
            start, end = row.get("start_1based"), row.get("end_1based_inclusive")
            if (
                source is None
                or type(start) is not int
                or type(end) is not int
                or not 1 <= start <= end <= source.length
            ):
                raise ValueError("CasAndra output does not match source coordinates")
            if key == "genes" and type(row["prediction"].get("is_cas")) is not bool:
                raise ValueError("Invalid CasAndra protein prediction")
            if key == "cas_cassettes" and not isinstance(row.get("final_classification"), dict):
                raise ValueError("Invalid CasAndra cassette classification")
        data[key] = rows
    count = sum(row["prediction"]["is_cas"] for row in data["genes"])
    if run.get("cas_proteins") != count:
        raise ValueError("CasAndra Cas-protein count mismatch")
    return {"gene_count": len(data["genes"]), "cas_gene_count": count, "cassette_count": len(data["cas_cassettes"])}


def validate_tracr(root: Path, output: Path, source: Path, *, mode: str, model_type: str, tool_options=None) -> dict:
    report = read_object(root, output / "report.json")
    if (
        report.get("schema") != {"name": "CRISPRtracrRNA-report", "version": "3.0.0"}
        or report.get("status") != "complete"
    ):
        raise ValueError("Unsupported or incomplete CRISPRtracrRNA report")
    config = report.get("configuration", {})
    if config.get("mode") != mode or config.get("model_type") != model_type:
        raise ValueError("tracrRNA effective settings differ from the requested workflow")
    for flag, value in (tool_options or {}).items():
        if config.get(flag[2:].replace("-", "_")) != value:
            raise ValueError("tracrRNA did not apply the requested scientific option")
    expected = {row.source_id: row for row in load_fasta_source_bindings(source)}
    manifest = read_object(root, output / "manifest.json")
    if manifest.get("algorithm") != "sha256":
        raise ValueError("Unsupported tracrRNA manifest")
    if config.get("short_anti_repeat_rescue") is not False:
        raise ValueError("Experimental tracrRNA policy was unexpectedly enabled")
    coordinates = report.get("coordinate_system", {})
    if (type(coordinates.get("indexing")) is not int or coordinates.get("indexing") != 1
            or coordinates.get("end_inclusive") is not True
            or coordinates.get("reference") != "input_forward_strand"):
        raise ValueError("Unsupported tracrRNA coordinates")
    files = sorted(path for path in source.iterdir() if path.is_file()) if source.is_dir() else [source]
    sequences = {header.split()[0]: sequence for path in files for header, sequence in _read_fasta_records(path)}
    if set(sequences) != set(expected) or any(
        (len(sequence), hashlib.sha256(sequence.encode("ascii")).hexdigest()) != (expected[key].length, expected[key].sha256)
        for key, sequence in sequences.items()
    ):
        raise ValueError("tracrRNA source changed during validation")
    observed, keys, prediction_ids = set(), set(), set()
    count = 0
    for record in report["records"]:
        identity = record["source"]
        bound = expected.get(identity.get("id"))
        key = record["key"]
        if not isinstance(key, str) or not re.fullmatch(r"f[0-9]{6}\.s[0-9]{6}", key) or key in keys:
            raise ValueError("Invalid tracrRNA record key")
        if (
            bound is None
            or identity["id"] in observed
            or (identity["length"], identity["sha256"]) != (bound.length, bound.sha256)
        ):
            raise ValueError("tracrRNA source binding mismatch")
        observed.add(identity["id"])
        keys.add(key)
        record_path = output / "work" / key / "record.json"
        if (
            manifest["files"][f"work/{key}/record.json"] != checksum(record_path)
            or read_object(root, record_path) != record
        ):
            raise ValueError("tracrRNA retained record mismatch")
        for prediction in record["predictions"]:
            start, end = prediction["interval"]["start"], prediction["interval"]["end"]
            if type(start) is not int or type(end) is not int or not 1 <= start <= end <= bound.length:
                raise ValueError("tracrRNA candidate outside source coordinates")
            if prediction["strand"] not in {"+", "-"} or len(prediction["sequence"]) != end - start + 1:
                raise ValueError("Invalid tracrRNA candidate sequence or strand")
            sequence = sequences[identity["id"]][start - 1:end]
            if prediction["strand"] == "-":
                sequence = sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]
            if prediction["sequence"] != sequence:
                raise ValueError("tracrRNA sequence does not reconstruct from source")
            export_id = prediction.get("export_id")
            if not isinstance(export_id, str) or not re.fullmatch(re.escape(key) + r"\.p[0-9]{6}", export_id) or export_id in prediction_ids:
                raise ValueError("Invalid tracrRNA export identity")
            prediction_ids.add(export_id)
            for name in ("identity", "coverage"):
                value = (prediction.get("anti_repeat") or {}).get(name)
                if value is not None and (type(value) not in {int, float} or not math.isfinite(value) or not 0 <= value <= 1):
                    raise ValueError("Invalid tracrRNA anti-repeat fraction")
            energy = (prediction.get("interaction") or {}).get("energy_kcal_mol")
            if energy is not None and (type(energy) not in {int, float} or not math.isfinite(energy)):
                raise ValueError("Invalid tracrRNA interaction energy")
            score = (prediction.get("ranking") or {}).get("score")
            if score is not None and (
                isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score)
            ):
                raise ValueError("Invalid tracrRNA ranking score")
        count += len(record["predictions"])
        if mode == "complete":
            work = output / "work" / key
            load_identified_arrays(work / "identify", expected_fasta_path=work / "input.fasta")
            validate_casandra(root, work / "casandra", work / "input.fasta")
    if observed != set(expected) or report["summary"] != {"sequence_count": len(expected), "prediction_count": count}:
        raise ValueError("tracrRNA record/count mismatch")
    return {"record_count": len(expected), "prediction_count": count, "mode": mode, "model_type": model_type}


def run_annotation_pipeline(
    input_path, output_path, *, runner_config, selection, cpu=None, timeout_seconds=3600, model_type="II", flank_length=500
):
    # Late import keeps the existing runner facade stable without a circular dependency.
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    source = Path(input_path).resolve(strict=True)
    load_fasta_source_bindings(source)
    if selection not in {"casandra", "tracrrna", "loci"} or model_type not in {"II", "V"}:
        raise ValueError("Unsupported annotation workflow or tracrRNA model")
    job = _empty_job_root(output_path)
    logs = job / "logs"
    logs.mkdir()
    result = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "status": "running",
        "workflow": selection,
        "tool_options": runner_config.tool_options,
        "stages": {},
    }
    output_marker = job / "integration_result.json"
    _atomic_json(output_marker, result)
    try:
        if selection == "casandra":
            files = sorted(source.iterdir()) if source.is_dir() else [source]
            native = job / "casandra"
            native.mkdir()
            counts = {"gene_count": 0, "cas_gene_count": 0, "cassette_count": 0}
            for index, fasta in enumerate(files, 1):
                output = native / f"record-{index:06d}"
                command = [
                    *runner_config.casandra,
                    "predict-genome",
                    "--genome",
                    str(fasta),
                    "--output",
                    str(output),
                    "--gene-mode",
                    "auto",
                    "--threads",
                    str(cpu or 1),
                ]
                command.extend(tool_arguments(runner_config, "casandra"))
                _run_process(f"casandra-{index:06d}", command, logs, timeout_seconds=timeout_seconds)
                for key, value in validate_casandra(job, output, fasta).items():
                    counts[key] += value
            result["stages"]["casandra"] = {"status": "completed", **counts}
        else:
            mode = "complete" if selection == "loci" else "model"
            config = job / "tracrrna-config.toml"
            # TOML string arrays and JSON string arrays share this scalar syntax.
            config.write_text(
                "identify_command = "
                + json.dumps(list(runner_config.crispridentify))
                + "\n"
                + "casandra_command = "
                + json.dumps(list(runner_config.casandra))
                + "\n"
            )
            output = job / "tracrrna"
            command = [
                *runner_config.crispr_tracrrna,
                "run",
                "--config",
                str(config),
                "--mode",
                mode,
                "--model-type",
                model_type,
                "--input",
                str(source),
                "--output",
                str(output),
                "--threads",
                str(cpu or 1),
                "--timeout",
                str(timeout_seconds),
            ]
            command.extend(tool_arguments(runner_config, "crispr_tracrrna"))
            stage = _run_process("tracrrna", command, logs, timeout_seconds=timeout_seconds)
            stage.update(validate_tracr(job, output, source, mode=mode, model_type=model_type, tool_options=runner_config.tool_options.get("crispr_tracrrna", {})))
            result["stages"]["tracrrna"] = stage
            if selection == "loci":
                result["stages"]["crispridentify"] = {"status": "completed", "provided_by": "tracrrna_complete"}
                result["stages"]["casandra"] = {"status": "completed", "provided_by": "tracrrna_complete"}
        if selection == "loci":
            from .leader import add_leader_context
            add_leader_context(job, source, result, runner_config, timeout_seconds, flank_length)
        result["status"] = "completed"
        _atomic_json(output_marker, result)
        return result
    except Exception as error:
        result["status"] = "failed"
        result["error"] = {"type": type(error).__name__, "message": str(error)}
        _atomic_json(output_marker, result)
        raise StageError(selection, str(error)) from error
