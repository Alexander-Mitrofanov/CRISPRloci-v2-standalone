"""Bounded FASTA and strict native annotations with independently bound sources."""
from __future__ import annotations

from dataclasses import replace
import io
import json
import math
import os
from pathlib import Path
import re
import stat
from typing import Any, BinaryIO

from crisprleader.coordinates import IUPAC_DNA, validate_array_geometry
from crisprleader.errors import InputError
from crisprleader.models import AnalysisConfig, ArrayAnnotation, RepeatAnchor, Segment, SourceRecord

MAX_JSON_BYTES = 64 * 1024 * 1024
MAX_JSON_DEPTH = 32
MAX_JSON_NODES = 1_000_000
MAX_HEADER_BYTES = 4096
MAX_SOURCES = 100_000
_ID = re.compile(r"^[!-~]{1,256}$")
_SHA = re.compile(r"^[0-9a-f]{64}$")


def checked_id(value: Any, label: str) -> str:
    if not isinstance(value, str) or not _ID.fullmatch(value) or ">" in value:
        raise InputError(f"{label} must be 1–256 printable ASCII characters without whitespace or >")
    return value


def _open_regular(path: str | Path) -> BinaryIO:
    try:
        descriptor = os.open(path, os.O_RDONLY | os.O_NONBLOCK)
        if not stat.S_ISREG(os.fstat(descriptor).st_mode):
            os.close(descriptor)
            raise InputError(f"input must be a regular file: {path}")
        return os.fdopen(descriptor, "rb")
    except OSError as error:
        raise InputError(f"cannot read input {path}: {error.strerror}") from error


def load_fasta(path: str | Path, max_bases: int = 50_000_000) -> dict[str, SourceRecord]:
    """Load plain multi-FASTA; sequence SHA256 covers upper-case IUPAC bases."""
    if type(max_bases) is not int or max_bases < 1:
        raise InputError("max_bases must be a positive integer")
    sources: dict[str, SourceRecord] = {}
    header: str | None = None
    identifier: str | None = None
    pieces = io.StringIO()
    total_bases = 0
    total_bytes = 0
    byte_bound = 2 * max_bases + 1024 * 1024

    def finish() -> None:
        if identifier is None:
            return
        if pieces.tell() == 0:
            raise InputError(f"FASTA record {identifier} has no sequence")
        sources[identifier] = SourceRecord(identifier, pieces.getvalue(), header=header)

    try:
        with _open_regular(path) as stream:
            while True:
                raw = stream.readline(min(max_bases + MAX_HEADER_BYTES + 1, byte_bound + 1))
                if not raw:
                    break
                total_bytes += len(raw)
                if total_bytes > byte_bound:
                    raise InputError(f"FASTA exceeds byte limit ({byte_bound}); reduce whitespace/headers")
                if raw.startswith(b">"):
                    if len(raw) > MAX_HEADER_BYTES:
                        raise InputError("FASTA header exceeds 4096 bytes")
                    finish()
                    header = raw[1:].rstrip(b"\r\n").decode("ascii")
                    if not header.strip() or any(ord(c) < 32 or ord(c) > 126 for c in header):
                        raise InputError("FASTA headers must be nonempty printable ASCII")
                    identifier = checked_id(header.split()[0], "FASTA source ID")
                    if identifier in sources:
                        raise InputError(f"duplicate FASTA source ID: {identifier}")
                    if len(sources) >= MAX_SOURCES:
                        raise InputError(f"FASTA exceeds {MAX_SOURCES} source records")
                    pieces = io.StringIO()
                else:
                    sequence = b"".join(raw.split()).decode("ascii").upper()
                    if not sequence:
                        continue
                    if identifier is None:
                        raise InputError("FASTA sequence appears before its >header")
                    if not set(sequence) <= IUPAC_DNA:
                        raise InputError(f"FASTA record {identifier} contains non-IUPAC DNA (no U/gaps)")
                    total_bases += len(sequence)
                    if total_bases > max_bases:
                        raise InputError(f"FASTA exceeds max_bases={max_bases}")
                    pieces.write(sequence)
            finish()
    except (UnicodeError, OSError) as error:
        raise InputError(f"invalid or unreadable ASCII FASTA: {error}") from error
    if not sources:
        raise InputError("FASTA must contain at least one nonempty sequence")
    return sources


def _object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise InputError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def validate_json_tree(value: Any) -> None:
    """Bound nested metadata too, and disallow non-finite or non-JSON values."""
    stack = [(value, 0)]
    count = 0
    while stack:
        node, depth = stack.pop()
        count += 1
        if depth > MAX_JSON_DEPTH or count > MAX_JSON_NODES:
            raise InputError("JSON exceeds nesting or element count limit")
        if isinstance(node, dict):
            if any(not isinstance(key, str) for key in node):
                raise InputError("JSON object keys must be strings")
            stack.extend((child, depth + 1) for child in node.values())
        elif isinstance(node, list):
            stack.extend((child, depth + 1) for child in node)
        elif isinstance(node, float):
            if not math.isfinite(node):
                raise InputError("JSON numeric values must be finite")
        elif node is not None and not isinstance(node, (str, int, bool)):
            raise InputError("annotations contain non-JSON data")


def read_json(path: str | Path, max_bytes: int = MAX_JSON_BYTES) -> Any:
    if type(max_bytes) is not int or max_bytes < 1:
        raise InputError("JSON byte budget must be a positive integer")
    try:
        with _open_regular(path) as stream:
            raw = stream.read(max_bytes + 1)
        if len(raw) > max_bytes:
            raise InputError(f"JSON input exceeds {max_bytes} bytes")
        value = json.loads(raw, object_pairs_hook=_object_pairs)
        validate_json_tree(value)
        return value
    except (ValueError, UnicodeError, RecursionError, OSError) as error:
        raise InputError(f"invalid or unreadable JSON in {path}: {error}") from error


def checked_object(
    value: Any, required: set[str], optional: set[str], label: str,
) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise InputError(f"{label} must be an object")
    missing, unknown = required - value.keys(), value.keys() - required - optional
    if missing or unknown:
        raise InputError(f"{label}: missing fields {sorted(missing)}, unknown fields {sorted(unknown)}")
    return value


def _metadata(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise InputError(f"{label} metadata must be an object")
    return value


def _dna(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value or not set(value.upper()) <= IUPAC_DNA:
        raise InputError(f"{label} must be nonempty IUPAC DNA")
    return value.upper()


def parse_segments(value: Any, label: str = "segments") -> tuple[Segment, ...]:
    if not isinstance(value, list) or len(value) not in (1, 2):
        raise InputError(f"{label} must be one interval or two origin-crossing intervals")
    result = []
    for item in value:
        segment = checked_object(item, {"start", "end"}, set(), label)
        try:
            result.append(Segment(segment["start"], segment["end"]))
        except ValueError as error:
            raise InputError(f"{label}: {error}") from error
    return tuple(result)


def load_annotation_document(
    document: Any, fasta_sources: dict[str, SourceRecord], config: AnalysisConfig | None = None,
) -> tuple[dict[str, SourceRecord], list[ArrayAnnotation]]:
    """Validate an in-memory native document, also used after format adaptation."""
    config = config or AnalysisConfig()
    for key in ("max_arrays", "max_candidates", "max_bases"):
        value = getattr(config, key)
        if type(value) is not int or value < 1:
            raise InputError(f"{key} must be a positive integer")
    if sum(source.length for source in fasta_sources.values()) > config.max_bases:
        raise InputError("source sequences exceed max_bases")
    validate_json_tree(document)
    doc = checked_object(document, {"schema_version", "coordinate_system", "sources", "arrays"},
                         set(), "annotations")
    if doc["schema_version"] != "1.0.0" or doc["coordinate_system"] != "zero_based_half_open":
        raise InputError("annotations require schema 1.0.0 and zero_based_half_open coordinates")
    if not isinstance(doc["sources"], list) or not 1 <= len(doc["sources"]) <= MAX_SOURCES:
        raise InputError("annotations sources must be a nonempty bounded list")
    sources: dict[str, SourceRecord] = {}
    for item in doc["sources"]:
        raw = checked_object(item, {"id", "length", "sha256", "topology"},
                             {"topology_provenance", "metadata"}, "source")
        identifier = checked_id(raw["id"], "source ID")
        if identifier in sources:
            raise InputError(f"duplicate annotation source ID: {identifier}")
        if identifier not in fasta_sources:
            raise InputError(f"annotation source ID {identifier} not present in FASTA")
        source = fasta_sources[identifier]
        if type(raw["length"]) is not int or raw["length"] != source.length:
            raise InputError(f"source {identifier} length does not match FASTA")
        if (not isinstance(raw["sha256"], str) or not _SHA.fullmatch(raw["sha256"])
                or raw["sha256"] != source.sha256):
            raise InputError(f"source {identifier} SHA256 does not match upper-case FASTA sequence")
        topology = raw["topology"]
        if topology not in ("linear", "circular", "unknown"):
            raise InputError(f"source {identifier} has invalid topology")
        provenance = raw.get("topology_provenance")
        if provenance is not None and (not isinstance(provenance, str) or not provenance.strip()):
            raise InputError("topology_provenance must be a nonempty string or null")
        if topology == "circular" and provenance is None:
            raise InputError("circular topology requires explicit topology_provenance")
        sources[identifier] = replace(source, topology=topology, topology_provenance=provenance,
                                      metadata=_metadata(raw.get("metadata", {}), "source"))
    if sources.keys() != fasta_sources.keys():
        raise InputError("annotation sources must bind every FASTA record exactly once")
    raw_arrays = doc["arrays"]
    if not isinstance(raw_arrays, list) or len(raw_arrays) > config.max_arrays:
        raise InputError(f"arrays must be a list with at most max_arrays={config.max_arrays} items")
    if len(raw_arrays) * 2 > config.max_candidates:
        raise InputError("two hypotheses per array exceed max_candidates")
    arrays: list[ArrayAnnotation] = []
    ids: set[str] = set()
    for item in raw_arrays:
        raw = checked_object(item, {"id", "source_id", "segments", "terminal_repeats",
                                   "consensus_repeat"},
                             {"strand", "strand_constraint", "category", "metadata"}, "array")
        identifier = checked_id(raw["id"], "array ID")
        if identifier in ids:
            raise InputError(f"duplicate array ID: {identifier}; IDs are document-global")
        ids.add(identifier)
        source_id = checked_id(raw["source_id"], "array source ID")
        if source_id not in sources:
            raise InputError(f"array {identifier} references unknown source {source_id}")
        raw_anchors = raw["terminal_repeats"]
        if not isinstance(raw_anchors, list) or len(raw_anchors) != 2:
            raise InputError(f"array {identifier} needs exactly two terminal repeats")
        anchors = []
        for anchor in raw_anchors:
            anchor = checked_object(anchor, {"segments", "sequence"}, set(), "terminal repeat")
            anchors.append(RepeatAnchor(parse_segments(anchor["segments"]),
                                        _dna(anchor["sequence"], "terminal-repeat sequence")))
        strand, constraint = raw.get("strand", "?"), raw.get("strand_constraint", False)
        category = raw.get("category")
        if category is not None and (not isinstance(category, str) or not category.strip()):
            raise InputError("array category must be a nonempty string or null")
        array = ArrayAnnotation(
            id=identifier, source_id=source_id, segments=parse_segments(raw["segments"]),
            terminal_repeats=(anchors[0], anchors[1]),
            consensus_repeat=_dna(raw["consensus_repeat"], "consensus_repeat"),
            strand=strand, strand_constraint=constraint, category=category,
            metadata=_metadata(raw.get("metadata", {}), "array"),
        )
        validate_array_geometry(sources[source_id], array)
        arrays.append(array)
    return sources, arrays


def load_inputs(
    fasta_path: str | Path, annotation_path: str | Path, config: AnalysisConfig | None = None,
) -> tuple[dict[str, SourceRecord], list[ArrayAnnotation]]:
    config = config or AnalysisConfig()
    sources = load_fasta(fasta_path, max_bases=config.max_bases)
    return load_annotation_document(read_json(annotation_path), sources, config)
