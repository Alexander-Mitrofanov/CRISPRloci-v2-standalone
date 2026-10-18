"""Versioned bundles, tabular exports, and source-bound computational validation."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from pathlib import Path
import time

from .inputs import validate_inventory, validate_records
from . import __version__
from .io import (MAX_BASES, MAX_RECORDS, InputError, LimitError, atomic_directory, canonical, decode_json,
                 digest, identifier, integer, read_json, verify_manifest, write_json, write_manifest)
from .search import Policy, analyze, reverse_complement

RESULT_FILES = {"analysis.json", "inputs.json", "queries.tsv", "matches.tsv", "exclusions.tsv", "coverage.tsv"}
REFERENCE_FILES = {"reference.json"}


def cell(value):
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        value = canonical(value).decode().strip().replace("\n", "")
    elif isinstance(value, bool):
        return "true" if value else "false"
    text = str(value)
    # Preserve exact strings in JSON; protect spreadsheet opening of TSV strings.
    if isinstance(value, str) and text.startswith(("=", "+", "-", "@", "\t", "\r")):
        text = "'" + text
    return text


def write_table(path: Path, rows: list[dict], fields: list[str]):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in rows:
            writer.writerow([cell(row.get(key)) for key in fields])


QUERY_FIELDS = ["id", "kind", "status", "source_id", "array_id", "ordinal", "start", "end", "boundary",
                "sequence", "sha256", "category", "source_orientation"]
HIT_FIELDS = ["id", "query_id", "target_id", "target_start", "target_end", "target_segments", "wraps_origin",
              "relative_strand", "query_length", "query_coverage", "identities", "substitutions", "gaps",
              "same_source_record", "same_declared_assembly", "array_overlap", "overlapping_array_ids", "excluded",
              "target_sequence_forward", "target_sequence_query_frame", "left_flank_forward", "right_flank_forward",
              "pam_status", "evalue", "confidence", "query_sha256", "target_sha256"]
COVERAGE_FIELDS = ["query_sha256", "target_id", "possible_windows", "assessed_windows", "ambiguous_windows", "query_occurrence_count"]


def export_tables(directory: Path, report: dict):
    write_table(directory / "queries.tsv", report["queries"], QUERY_FIELDS)
    write_table(directory / "matches.tsv", [x for x in report["matches"] if not x["excluded"]], HIT_FIELDS)
    write_table(directory / "exclusions.tsv", [x for x in report["matches"] if x["excluded"]], HIT_FIELDS)
    write_table(directory / "coverage.tsv", report["coverage"], COVERAGE_FIELDS)


def write_result(destination: Path, inventory: dict, targets: list[dict], policy: Policy,
                 source_targets: bool = False, reference: dict | None = None) -> dict:
    # Search takes place inside a transaction directory; caps/errors never publish a completed bundle.
    with atomic_directory(destination) as work:
        deadline = time.monotonic() + policy.timeout
        validate_reference_identity(reference, targets)
        report = analyze(inventory, targets, policy, source_targets)
        report["reference"] = reference
        inputs = {"schema_version": "1.0.0", "inventory": inventory, "targets": targets,
                  "policy": policy.to_dict(), "source_targets": source_targets, "reference": reference}
        write_json(work / "inputs.json", inputs)
        write_json(work / "analysis.json", report)
        export_tables(work, report)
        write_manifest(work, "crisprspacer.result")
        if time.monotonic() > deadline:
            raise LimitError("Search/publication deadline exceeded; no completed result was published")
    return report


def prepare_bundle(destination: Path, inventory: dict):
    validate_inventory(inventory)
    with atomic_directory(destination) as work:
        write_json(work / "inventory.json", inventory)
        write_manifest(work, "crisprspacer.inventory")


def load_inventory(path: Path):
    if path.is_dir():
        _, payloads = verify_manifest(path, "crisprspacer.inventory", {"inventory.json"})
        return validate_inventory(decode_json(payloads["inventory.json"]))
    return validate_inventory(read_json(path))


def build_reference(destination: Path, targets: list[dict], name: str, source_digest: str,
                    source_label: str, release: str):
    validate_records(targets)
    identifier(name, "reference name")
    identifier(source_label, "source label")
    identifier(release, "source release")
    content = {"schema": {"name": "crisprspacer.reference", "version": "1.0.0"},
               "name": name, "source_label": source_label, "source_release": release,
               "packaged_at": datetime.now(timezone.utc).isoformat(),
               "source_file_sha256": source_digest, "inventory_sha256": digest(canonical(targets)),
               "record_count": len(targets), "base_count": sum(r["length"] for r in targets), "records": targets,
               "normalization": "uppercase DNA; unique first-token IDs; explicit assembly/topology metadata"}
    with atomic_directory(destination) as work:
        write_json(work / "reference.json", content)
        write_manifest(work, "crisprspacer.reference")


def load_reference(directory: Path) -> tuple[list[dict], dict]:
    manifest, payloads = verify_manifest(directory, "crisprspacer.reference", REFERENCE_FILES)
    content = decode_json(payloads["reference.json"])
    if not isinstance(content, dict) or content.get("schema") != {"name": "crisprspacer.reference", "version": "1.0.0"}:
        raise InputError("Unsupported reference schema")
    records = content.get("records")
    validate_records(records)
    integer(content.get("record_count"), "Reference record count", 1, MAX_RECORDS)
    integer(content.get("base_count"), "Reference base count", 1, MAX_BASES)
    if not records or content.get("inventory_sha256") != digest(canonical(records)) or content.get("record_count") != len(records) or content.get("base_count") != sum(r["length"] for r in records):
        raise InputError("Reference sequence inventory mismatch")
    for key in ("name", "source_label", "source_release"):
        identifier(content.get(key), key)
    return records, {"name": content["name"], "source_label": content["source_label"],
                     "source_release": content["source_release"], "inventory_sha256": content["inventory_sha256"],
                     "reference_sha256": manifest["files"]["reference.json"]["sha256"]}


def validate_reference_identity(reference: object, targets: list[dict]):
    if reference is None:
        return
    expected = {"name", "source_label", "source_release", "inventory_sha256", "reference_sha256"}
    if not isinstance(reference, dict) or set(reference) != expected:
        raise InputError("Invalid retained reference identity")
    for key in ("name", "source_label", "source_release"):
        identifier(reference[key], key)
    for key in ("inventory_sha256", "reference_sha256"):
        value = reference[key]
        if not isinstance(value, str) or len(value) != 64 or set(value) - set("0123456789abcdef"):
            raise InputError("Invalid reference digest")
    if reference["inventory_sha256"] != digest(canonical(targets)):
        raise InputError("Reference identity does not bind retained targets")


def validate_result(directory: Path) -> dict:
    _, payloads = verify_manifest(directory, "crisprspacer.result", RESULT_FILES)
    inputs, report = decode_json(payloads["inputs.json"]), decode_json(payloads["analysis.json"])
    if not isinstance(inputs, dict) or set(inputs) != {"schema_version", "inventory", "targets", "policy", "source_targets", "reference"} or inputs["schema_version"] != "1.0.0":
        raise InputError("Invalid retained input schema")
    if type(inputs["source_targets"]) is not bool:
        raise InputError("Source-target relation must be explicit boolean")
    inventory = validate_inventory(inputs["inventory"])
    target_map = validate_records(inputs["targets"])
    validate_reference_identity(inputs["reference"], inputs["targets"])
    query_map = {r["id"]: r for r in inventory["queries"]}
    if not isinstance(report, dict) or not isinstance(report.get("matches"), list):
        raise InputError("Invalid analysis report")
    # Direct source/metric checks do not depend on candidate discovery.
    identities = set()
    for hit in report["matches"]:
        if not isinstance(hit, dict):
            raise InputError("Invalid alignment row")
        query, target = query_map.get(hit.get("query_id")), target_map.get(hit.get("target_id"))
        if query is None or target is None or query["kind"] != "observed":
            raise InputError("Unknown alignment source")
        strand = hit.get("relative_strand")
        if strand not in {"+", "-", "both"}:
            raise InputError("Invalid alignment strand")
        parts = hit.get("target_segments")
        if not isinstance(parts, list) or not 1 <= len(parts) <= 2:
            raise InputError("Invalid alignment segments")
        sequence = ""
        for part in parts:
            if not isinstance(part, dict) or set(part) != {"start", "end"}:
                raise InputError("Invalid segment fields")
            start, end = part["start"], part["end"]
            if type(start) is not int or type(end) is not int or not 1 <= start <= end <= target["length"]:
                raise InputError("Segment outside target")
            sequence += target["sequence"][start-1:end]
        if len(parts) == 2 and (target["topology"] != "circular" or parts[0]["end"] != target["length"] or parts[1]["start"] != 1):
            raise InputError("Invalid circular segment geometry")
        oriented = reverse_complement(sequence) if strand == "-" else sequence
        qseq = query["sequence"]
        if len(sequence) != len(qseq) or len(sequence) > target["length"] or set(sequence) - set("ACGT"):
            raise InputError("Alignment is not full-length unambiguous DNA")
        distance = sum(a != b for a, b in zip(oriented, qseq))
        if hit.get("substitutions") != distance or hit.get("identities") != len(qseq)-distance or hit.get("gaps") != 0:
            raise InputError("Alignment metrics do not reconstruct")
        if hit.get("target_sequence_forward") != sequence or hit.get("target_sequence_query_frame") != oriented:
            raise InputError("Aligned sequence does not reconstruct")
        if strand == "both" and qseq != reverse_complement(qseq):
            raise InputError("Both-strand label requires a palindromic query")
        identity = (hit["query_id"], hit["target_id"], tuple((p["start"], p["end"]) for p in parts), strand)
        if identity in identities:
            raise InputError("Duplicate alignment occurrence")
        identities.add(identity)
    try:
        policy = Policy(**inputs["policy"])
    except TypeError as exc:
        raise InputError("Unsupported retained search policy") from exc
    expected = analyze(inventory, inputs["targets"], policy, inputs["source_targets"])
    expected["reference"] = inputs["reference"]
    if report.get("tool_version") not in {"0.3.0", __version__}:
        raise InputError("Unsupported analysis tool version")
    expected["tool_version"] = report["tool_version"]
    if canonical(expected) != canonical(report):
        raise InputError("Analysis differs from source-bound deterministic replay")
    # Compare export bytes to regenerated tables, not only their mutable checksums.
    import tempfile
    with tempfile.TemporaryDirectory(prefix="crisprspacer-validate-") as temp:
        root = Path(temp)
        export_tables(root, expected)
        for name in ("queries.tsv", "matches.tsv", "exclusions.tsv", "coverage.tsv"):
            if (root / name).read_bytes() != payloads[name]:
                raise InputError("Tabular artifact disagrees with report")
    return {"status": "valid", "checks": ["artifact_hashes", "source_binding", "alignment_metrics", "deterministic_replay", "table_consistency"],
            "biological_validation": False, "counts": report["counts"]}
