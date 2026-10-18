"""Versioned local reference inventories with explicit historical annotations."""

from __future__ import annotations

from pathlib import Path
import re

from .io import InputError, canonical, decode_json, digest, normalize_records, read_bytes

SCHEMA = {"name": "crisprmap.reference", "version": "1.0.0"}
NORMALIZATION = "uppercase_U_to_T"
LABELS = ("family", "motif", "accession")


def _label(value, field):
    if value is not None and (not isinstance(value, str) or not 1 <= len(value) <= 500
                              or any(ord(c) < 32 or ord(c) > 126 for c in value)):
        raise InputError(f"Reference {field} must be null or bounded printable text")
    return value


def _provenance(value):
    if not isinstance(value, dict) or not 3 <= len(value) <= 30:
        raise InputError("Reference provenance must be an object")
    for key in ("namespace", "source_url", "source_sha256"):
        if not isinstance(value.get(key), str) or not value[key]:
            raise InputError(f"Reference provenance requires {key}")
    if not re.fullmatch(r"[0-9a-f]{64}", value["source_sha256"]):
        raise InputError("Reference source_sha256 must be a SHA256 digest")
    for key, item in value.items():
        if not isinstance(key, str) or not 1 <= len(key) <= 100:
            raise InputError("Invalid provenance key")
        if item is not None and (not isinstance(item, str) or len(item) > 4000
                                 or any(ord(c) < 32 or ord(c) > 126 for c in item)):
            raise InputError("Provenance values must be null or bounded printable text")
    return dict(value)


def _content_digest(document):
    return digest(canonical({key: value for key, value in document.items() if key != "content_sha256"}))


def build_reference(records: list[dict], metadata_rows: dict, provenance: dict) -> dict:
    """Bind caller-owned labels to normalized records; missing labels stay null."""
    records = normalize_records(records)
    if not isinstance(metadata_rows, dict) or not set(metadata_rows) <= {r["id"] for r in records}:
        raise InputError("Metadata must map known reference identifiers to objects")
    entries = []
    for record in records:
        row = metadata_rows.get(record["id"], {})
        if not isinstance(row, dict) or not set(row) <= set(LABELS):
            raise InputError("Reference metadata supports family, motif and accession only")
        entry = {**record, "sequence_sha256": digest(record["sequence"])}
        entry.update({label: _label(row.get(label), label) for label in LABELS})
        entries.append(entry)
    result = {"schema": dict(SCHEMA), "provenance": _provenance(provenance),
              "normalization": NORMALIZATION, "entries": entries}
    result["content_sha256"] = _content_digest(result)
    return result


def validate_reference(document: object) -> dict:
    if not isinstance(document, dict) or set(document) != {
        "schema", "provenance", "normalization", "entries", "content_sha256"
    }:
        raise InputError("Reference fields do not match schema")
    if document["schema"] != SCHEMA or document["normalization"] != NORMALIZATION:
        raise InputError("Unsupported reference schema or normalization")
    entries = document["entries"]
    if not isinstance(entries, list):
        raise InputError("Reference entries must be a list")
    fields = {"id", "description", "sequence", "sequence_sha256", *LABELS}
    if any(not isinstance(e, dict) or set(e) != fields for e in entries):
        raise InputError("Reference entry fields do not match schema")
    records = normalize_records([{key: row[key] for key in ("id", "sequence", "description")} for row in entries])
    rebuilt = build_reference(
        records,
        {row["id"]: {key: row[key] for key in LABELS} for row in entries}, document["provenance"])
    if document != rebuilt:
        raise InputError("Reference digest, normalized sequence or annotation binding is inconsistent")
    return document


def load_reference(path: Path) -> tuple[dict, str]:
    raw = read_bytes(path)
    return validate_reference(decode_json(raw)), digest(raw)
