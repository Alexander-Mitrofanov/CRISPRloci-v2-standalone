"""Shell-free CRISPRspacer database search with source-bound output gates."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re

from .parser import _read_fasta_records

LIMIT = 20_000_000
POLICY = {
    "method": "full-length-exact", "substitutions": 0, "gaps": 0,
    "strands": "both", "topology": "linear", "array_overlap": "unknown", "pam": "not_evaluated",
}


def _hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _payload(root, path):
    current = root
    for part in path.relative_to(root).parts:
        current /= part
        if current.is_symlink():
            raise ValueError("Symlink in CRISPRspacer output")
    if not path.is_file() or path.stat().st_size > LIMIT:
        raise ValueError("Missing or oversized CRISPRspacer output")
    with path.open("rb") as stream:
        value = stream.read(LIMIT + 1)
    if len(value) > LIMIT:
        raise ValueError("CRISPRspacer output exceeded its byte budget")
    return value


def _json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate JSON key")
            result[key] = value
        return result
    def floating(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError("Nonfinite output number")
        return result
    def invalid(value):
        raise ValueError("Invalid JSON constant: " + value)
    return json.loads(data, object_pairs_hook=pairs, parse_float=floating, parse_constant=invalid)


def _records(source):
    records = {}
    paths = sorted(source.iterdir()) if source.is_dir() else [source]
    for path in paths:
        if path.is_symlink() or not path.is_file() or path.stat().st_size > LIMIT:
            raise ValueError("Protospacer input must contain regular FASTA files")
        for header, sequence in _read_fasta_records(path):
            identity = header.split()[0]
            if identity in records:
                raise ValueError("Duplicate protospacer source identity")
            records[identity] = sequence
            if len(records) > 1000 or sum(map(len, records.values())) > 5_000_000:
                raise ValueError("Protospacer target inventory exceeds its input budget")
    if not records:
        raise ValueError("No protospacer target records")
    return records


def validate_protospacer_report(report, records, database_sha256, database_path=None, checkpoint=None):
    """Check source reconstruction independently of the standalone search engine."""
    if (
        report.get("schema") != {"name": "crisprspacer.database-search", "version": "1.0.0"}
        or report.get("tool") not in ({"name": "CRISPRspacer", "version": "0.2.0"}, {"name": "CRISPRspacer", "version": "0.3.0"},
                                     {"name": "CRISPRspacer", "version": "0.3.1"})
        or report.get("status") != "completed" or report.get("search_complete") is not True
        or report.get("policy") != POLICY
    ):
        raise ValueError("Unsupported or incomplete protospacer output")
    reference = report["reference"]
    if not isinstance(reference, dict) or reference.get("database_sha256") != database_sha256:
        raise ValueError("Protospacer reference differs from pinned database")
    for name in ("source_file_sha256", "database_sha256"):
        if not isinstance(reference.get(name), str) or not re.fullmatch("[0-9a-f]{64}", reference[name]):
            raise ValueError("Invalid protospacer reference digest")
    for name in ("name", "source_url", "source_release"):
        value = reference.get(name)
        if not isinstance(value, str) or not value or len(value) > 2000 or any(ord(c) < 32 for c in value):
            raise ValueError("Invalid protospacer reference identity")
    for name in ("record_count", "eligible_record_count", "skipped_ambiguous", "skipped_length"):
        if type(reference.get(name)) is not int or reference[name] < 0:
            raise ValueError("Invalid protospacer reference count")
    if (
        reference.get("seed_length") != 16 or reference.get("min_length") != 16 or reference.get("max_length") != 1000
        or reference["record_count"] != sum(reference[x] for x in ("eligible_record_count", "skipped_ambiguous", "skipped_length"))
    ):
        raise ValueError("Unsupported protospacer reference policy/counts")
    targets = report["targets"]
    if not isinstance(targets, list) or len(targets) != len(records):
        raise ValueError("Protospacer source count mismatch")
    seen = set()
    for target in targets:
        key = target.get("id")
        if key in seen or key not in records or target != {"id": key, "length": len(records[key]), "sha256": _hash(records[key])}:
            raise ValueError("Protospacer target binding mismatch")
        seen.add(key)
    rows = report["matches"]
    if not isinstance(rows, list) or len(rows) > 10_000:
        raise ValueError("Protospacer match limit exceeded")
    unique, spacers = set(), {}
    hashes = {key: _hash(value) for key, value in records.items()}
    for row in rows:
        sequence, key, sid = row["spacer_sequence"], row["target_id"], row["spacer_id"]
        if not isinstance(sid, str) or not sid or len(sid) > 500 or any(ord(c) < 33 or ord(c) > 126 for c in sid):
            raise ValueError("Invalid spacer identity")
        if not isinstance(sequence, str) or not 16 <= len(sequence) <= 1000 or set(sequence) - set("ACGT"):
            raise ValueError("Invalid spacer sequence")
        if key not in records or row["spacer_sha256"] != _hash(sequence) or row["target_sha256"] != hashes[key]:
            raise ValueError("Protospacer sequence binding mismatch")
        if sid in spacers and spacers[sid] != sequence:
            raise ValueError("Spacer identity has inconsistent sequences")
        spacers[sid] = sequence
        start, end, strand = row["target_start"], row["target_end"], row["relative_strand"]
        if type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(records[key]) or end-start+1 != len(sequence):
            raise ValueError("Protospacer interval is invalid")
        reverse = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if strand not in {"+", "-", "both"} or (strand == "both") != (sequence == reverse):
            raise ValueError("Protospacer strand is inconsistent")
        expected = reverse if strand == "-" else sequence
        if records[key][start-1:end] != expected:
            raise ValueError("Protospacer match does not reconstruct")
        identity = sid, key, start, end, strand
        if identity in unique:
            raise ValueError("Duplicate protospacer match")
        unique.add(identity)
        for name in ("identities", "query_length"):
            if type(row.get(name)) is not int or row[name] != len(sequence):
                raise ValueError("Invalid protospacer alignment length")
        if any(type(row.get(name)) is not int or row[name] != 0 for name in ("substitutions", "gaps")):
            raise ValueError("Protospacer search must be exact")
        if row.get("query_coverage") != 1.0 or isinstance(row.get("query_coverage"), bool) or row.get("array_overlap") != "unknown" or row.get("pam_status") != "not_evaluated" or row.get("confidence") is not None or row.get("evalue") is not None:
            raise ValueError("Unsupported protospacer interpretation")
    counts = report["counts"]
    for name in ("matches", "matched_spacers", "targets", "target_bases", "seed_windows", "ambiguous_seed_windows", "candidates_checked"):
        if type(counts.get(name)) is not int or counts[name] < 0:
            raise ValueError("Invalid protospacer result count")
    windows = sum(max(0, len(value)-15) for value in records.values())
    ambiguous = sum(any(base not in "ACGT" for base in value[i:i+16]) for value in records.values() for i in range(max(0, len(value)-15)))
    if counts["matches"] != len(rows) or counts["matched_spacers"] != len(spacers) or counts["targets"] != len(records) or counts["target_bases"] != sum(map(len, records.values())) or counts["seed_windows"] != windows or counts["ambiguous_seed_windows"] != ambiguous or counts["candidates_checked"] < len(rows):
        raise ValueError("Protospacer result counts disagree")
    outcome = "no_assessable_windows" if windows == ambiguous or reference["eligible_record_count"] == 0 else "matches" if rows else "no_matches"
    if report.get("outcome") != outcome:
        raise ValueError("Protospacer outcome disagrees with evidence")
    from .protospacer_reference import bind_catalogue
    bind_catalogue(database_path, database_sha256, reference, spacers, checkpoint)
    return counts


def run_protospacer_pipeline(input_path, output_path, *, runner_config, database_path, database_sha256, timeout_seconds=3600):
    from .runner import StageError, _atomic_json, _empty_job_root, _run_process

    if database_path is None or not isinstance(database_sha256, str) or not re.fullmatch("[0-9a-f]{64}", database_sha256):
        raise ValueError("Protospacer requires an operator database and its SHA-256")
    database = Path(database_path).absolute()
    if not database.is_dir() or any(p.is_symlink() for p in (database, *database.parents)):
        raise ValueError("Protospacer database must be an existing directory without symlinks")
    records = _records(Path(input_path).absolute())
    job = _empty_job_root(output_path)
    (job / "logs").mkdir()
    private = job / "protospacer-input"
    private.mkdir()
    fasta = private / "targets.fasta"
    fasta.write_text("".join(f">{key}\n" + "".join(value[i:i+80] + "\n" for i in range(0, len(value), 80))
                             for key, value in records.items()), encoding="ascii")
    native = job / "protospacer"
    result = {"schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"}, "status": "running", "workflow": "protospacer", "stages": {}}
    _atomic_json(job / "integration_result.json", result)
    try:
        command = [*runner_config.crisprspacer, "database-search", "--database", str(database), "--targets", str(fasta), "--expected-database-sha256", database_sha256, "--output", str(native), "--timeout", str(max(1, min(300, int(timeout_seconds))))]
        stage = _run_process("crisprspacer", command, job / "logs", timeout_seconds=timeout_seconds)
        manifest = _json(_payload(job, native / "manifest.json"))
        if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != "crisprspacer.database-search" or set(manifest.get("files", {})) != {"analysis.json", "inputs.json", "matches.tsv"}:
            raise ValueError("Invalid protospacer manifest")
        payloads = {}
        for name in ("analysis.json", "inputs.json", "matches.tsv"):
            data = _payload(job, native / name)
            if manifest["files"][name] != {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}:
                raise ValueError("Protospacer artifact checksum mismatch")
            payloads[name] = data
        report = _json(payloads["analysis.json"])
        counts = validate_protospacer_report(report, records, database_sha256, database)
        retained = _json(payloads["inputs.json"])
        expected_targets = [{"id": key, "sequence": value, "sha256": _hash(value), "length": len(value), "assembly_id": None, "topology": "linear"} for key, value in records.items()]
        if retained != {"schema_version": "1.0.0", "targets": expected_targets, "target_file_sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(), "reference": report["reference"]}:
            raise ValueError("Protospacer retained input mismatch")
        stage.update(counts=counts, database_sha256=database_sha256)
        result["stages"]["crisprspacer"] = stage
        result["status"] = "completed"
    except Exception as error:
        result.update(status="failed", error={"type": type(error).__name__, "message": str(error)})
        _atomic_json(job / "integration_result.json", result)
        raise StageError("crisprspacer", str(error)) from error
    _atomic_json(job / "integration_result.json", result)
    return result
