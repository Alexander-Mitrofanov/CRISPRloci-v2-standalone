"""Adapt canonical CRISPRidentify reports without importing the detector.

Report-v1 coordinates, repeat strings, and consensus are input-forward; an
orientation prediction never changes that frame. Only accepted categories
(Bona-fide/Possible) enter leader analysis. Original report fields and excluded
IDs/categories remain explicit provenance. Legacy presentation JSON is rejected.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Any

from crisprleader.coordinates import IUPAC_DNA
from crisprleader.errors import InputError
from crisprleader.models import ArrayAnnotation, RepeatAnchor, Segment, SourceRecord, json_record
from crisprleader.sequence_io import (
    MAX_JSON_BYTES, checked_id, load_annotation_document, load_fasta, read_json,
)

_CATEGORIES = frozenset(("Bona-fide", "Possible", "Alternative", "Possible Discarded", "Low score"))
_ACCEPTED = frozenset(("Bona-fide", "Possible"))
_VERSION = re.compile(r"^1\.[0-9]+\.[0-9]+$")
_HASH = re.compile(r"^[0-9a-f]{64}$")


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise InputError(f"Identify {label} must be an object")
    return value


def _interval(value: Any, source: SourceRecord, label: str) -> Segment:
    item = _mapping(value, label)
    start, end, length = item.get("start"), item.get("end"), item.get("length")
    if (type(start) is not int or type(end) is not int or type(length) is not int
            or not 1 <= start <= end <= source.length or length != end - start + 1):
        raise InputError(f"Identify {label} must be a source-bound 1-based inclusive interval")
    return Segment(start - 1, end)


def _sequence(value: Any, expected: str, label: str) -> str:
    if not isinstance(value, str) or value.upper() != expected:
        raise InputError(f"Identify {label} does not match its source-forward FASTA slice")
    return expected


def _orientation(value: Any) -> str:
    orientation = _mapping(value, "orientation")
    status, strand = orientation.get("status"), orientation.get("strand")
    if status == "computed" and strand in ("+", "-"):
        expected_labels = ("Forward",) if strand == "+" else ("Reverse", "Reversed")
        if orientation.get("label") not in expected_labels:
            raise InputError("Identify orientation label conflicts with computed strand")
        return strand
    if status in ("not_requested", "not_computed") and strand is None:
        return "?"
    raise InputError("Identify orientation has an unknown status or conflicting strand")


def _validation(value: Any, label: str) -> None:
    item = _mapping(value, label)
    if (item.get("status") != "valid"
            or item.get("method") != "independent_original_sequence_slicing"):
        raise InputError(f"Identify {label} is not a validated canonical result")


def _topology(report: dict[str, Any], source: SourceRecord) -> tuple[str, str | None]:
    extensions = _mapping(report.get("extensions", {}), "extensions")
    if "circular_origin_loci" in extensions:
        circular = _mapping(extensions["circular_origin_loci"], "circular_origin_loci")
        schema = _mapping(circular.get("schema"), "circular-origin schema")
        if (schema.get("name") != "CRISPRidentify-circular-origin-loci"
                or schema.get("version") != "1.0.0"
                or circular.get("source_sha256") != source.sha256):
            raise InputError("unsupported or unbound Identify circular-origin extension; use native annotations")
        joined, count = circular.get("joined_loci"), circular.get("joined_locus_count")
        if not isinstance(joined, list) or type(count) is not int or count != len(joined):
            raise InputError("Identify circular-origin extension has invalid joined-locus accounting")
        if joined:
            raise InputError(
                "Identify contains joined origin-crossing arrays without a canonical consensus "
                "contract supported by this adapter; provide native annotations with verified "
                "source-forward consensus and terminal-repeat anchors"
            )
    if "sequence_topology" not in extensions:
        return "unknown", None
    evidence = _mapping(extensions["sequence_topology"], "sequence_topology")
    if (evidence.get("schema_name") != "CRISPRidentify-sequence-topology-evidence"
            or evidence.get("schema_version") != "1.0.0"):
        raise InputError("unsupported Identify sequence-topology schema; use native annotations")
    binding = _mapping(evidence.get("source"), "topology source")
    if (binding.get("source_id") != source.id
            or type(binding.get("source_length")) is not int
            or binding.get("source_length") != source.length
            or binding.get("source_sha256") != source.sha256
            or binding.get("hash_normalization") != "uppercase_ascii_iupac_without_reformatting"):
        raise InputError("Identify topology evidence is not bound to the FASTA source")
    topology = evidence.get("topology")
    provenance = _mapping(evidence.get("provenance"), "topology provenance")
    if evidence.get("inference_used") is not False:
        raise InputError("Identify topology must have explicit, non-inferred evidence")
    if topology == "unknown":
        if (evidence.get("evidence_complete") is not False
                or evidence.get("evidence_ambiguous") is not True
                or provenance.get("method") != "unknown_default"
                or any(provenance.get(k) is not None for k in ("locator", "date", "evidence_sha256"))):
            raise InputError("Identify unknown topology evidence is contradictory")
        return "unknown", None
    if topology not in ("linear", "circular"):
        raise InputError("Identify has unknown topology value")
    if (evidence.get("evidence_complete") is not True
            or evidence.get("evidence_ambiguous") is not False
            or provenance.get("method") != "frozen_insdc_flatfile"
            or any(not isinstance(provenance.get(k), str) or not provenance[k].strip()
                   for k in ("accession", "locator", "date", "evidence_sha256"))
            or not _HASH.fullmatch(provenance["evidence_sha256"])):
        raise InputError("Identify known topology requires complete frozen source evidence")
    return topology, "CRISPRidentify topology assertion: " + json.dumps(provenance, sort_keys=True)


def _array(raw: Any, source: SourceRecord) -> ArrayAnnotation:
    raw = _mapping(raw, "array")
    identifier = checked_id(raw.get("id"), "Identify array ID")
    category = raw.get("category")
    if not isinstance(category, str) or category not in _CATEGORIES:
        raise InputError(f"Identify array {identifier} has unknown category {category!r}")
    _validation(raw.get("validation"), "array validation")
    if (raw["validation"].get("source_reconstructed") is not True
            or raw["validation"].get("issues") != []):
        raise InputError(f"Identify array {identifier} has invalid reconstruction evidence")
    segment = _interval(raw.get("source_interval"), source, "array interval")
    if "source_sequence" in raw:
        _sequence(raw["source_sequence"], source.sequence[segment.start:segment.end], "array sequence")
    repeats = raw.get("repeats")
    if (not isinstance(repeats, list) or len(repeats) < 2
            or type(raw.get("repeat_count")) is not int or raw["repeat_count"] != len(repeats)):
        raise InputError("Identify repeat_count must match at least two source-verified repeats")
    anchors: list[RepeatAnchor] = []
    for index, item in enumerate(repeats, 1):
        item = _mapping(item, "repeat")
        if type(item.get("index")) is not int or item["index"] != index:
            raise InputError("Identify repeats must be ordered with contiguous 1-based indices")
        repeat = _interval(item.get("source_interval"), source, "repeat interval")
        sequence = _sequence(item.get("sequence"), source.sequence[repeat.start:repeat.end], "repeat")
        if repeat.start < segment.start or repeat.end > segment.end:
            raise InputError("Identify repeat is outside its array")
        if anchors and anchors[-1].segments[0].end > repeat.start:
            raise InputError("Identify repeats overlap or are not source-forward ordered")
        anchors.append(RepeatAnchor((repeat,), sequence))
    if anchors[0].segments[0].start != segment.start or anchors[-1].segments[0].end != segment.end:
        raise InputError("Identify terminal repeats must match array endpoints")
    spacers = raw.get("spacers")
    if (not isinstance(spacers, list) or type(raw.get("spacer_count")) is not int
            or raw["spacer_count"] != len(repeats) - 1):
        raise InputError("Identify spacer_count must describe every inter-repeat gap")
    by_index: dict[int, Any] = {}
    for item in spacers:
        item = _mapping(item, "spacer")
        spacer_index = item.get("index")
        if type(spacer_index) is not int or not 1 <= spacer_index < len(repeats) or spacer_index in by_index:
            raise InputError("Identify spacer indices must be unique valid repeat-gap indices")
        by_index[spacer_index] = item
    for index, (left, right) in enumerate(zip(anchors, anchors[1:]), 1):
        start, end = left.segments[0].end, right.segments[0].start
        item = by_index.pop(index, None)
        if start == end:
            if item is not None:
                boundary = _mapping(item.get("source_boundary"), "empty-spacer boundary")
                if (item.get("sequence") != "" or "source_interval" in item
                        or type(boundary.get("position")) is not int or boundary["position"] != start
                        or boundary.get("coordinate_system") != "0-based interbase"
                        or boundary.get("reference") != "input_forward_strand"):
                    raise InputError("Identify empty spacer must be an exact source interbase boundary")
        else:
            if item is None or _interval(item.get("source_interval"), source, "spacer interval") != Segment(start, end):
                raise InputError("Identify spacer must exactly fill its inter-repeat source gap")
            _sequence(item.get("sequence"), source.sequence[start:end], "spacer")
    consensus = raw.get("consensus_repeat")
    if not isinstance(consensus, str) or not consensus or not set(consensus.upper()) <= IUPAC_DNA:
        raise InputError("Identify consensus_repeat must be nonempty source-forward IUPAC DNA")
    # Global IDs are independent of report filenames and colliding per-source display IDs.
    key = hashlib.sha256(json.dumps([source.id, identifier], separators=(",", ":")).encode()).hexdigest()
    return ArrayAnnotation(
        id="identify:" + key, source_id=source.id, segments=(segment,),
        terminal_repeats=(anchors[0], anchors[-1]), consensus_repeat=consensus.upper(),
        strand=_orientation(raw.get("orientation")), strand_constraint=False,  # type: ignore[arg-type]
        category=category, metadata={"identify_original": raw,
                                     "sequence_frame": "input_forward_strand",
                                     "orientation_is_supplied_evidence": True},
    )


def adapt_identify(fasta_path: str | Path, report_paths: list[Path]) -> dict[str, Any]:
    """Return validated native annotations for source-bound canonical report v1."""
    sources = load_fasta(fasta_path)
    if not isinstance(report_paths, (list, tuple)) or not report_paths:
        raise InputError("at least one canonical Identify report path is required")
    if len(report_paths) > len(sources):
        raise InputError("at most one Identify report is allowed per FASTA source")
    native_sources: dict[str, dict[str, Any]] = {}
    arrays: list[dict[str, Any]] = []
    remaining_bytes = MAX_JSON_BYTES
    for path in report_paths:
        report = _mapping(read_json(path, max_bytes=remaining_bytes), "report")
        semantic_bytes = json.dumps(
            report, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
        remaining_bytes -= len(semantic_bytes)
        if remaining_bytes < 0:
            raise InputError("combined Identify reports exceed bounded document size")
        report_digest = hashlib.sha256(semantic_bytes).hexdigest()
        del semantic_bytes
        schema = _mapping(report.get("schema"), "schema")
        version = schema.get("version")
        if (schema.get("name") != "CRISPRidentify-report"
                or not isinstance(version, str) or not _VERSION.fullmatch(version)):
            raise InputError("adapter requires canonical CRISPRidentify-report major version 1, not legacy JSON")
        frame = _mapping(report.get("coordinate_system"), "coordinate_system")
        if (type(frame.get("indexing")) is not int or frame["indexing"] != 1
                or frame.get("end_inclusive") is not True
                or frame.get("reference") != "input_forward_strand"
                or frame.get("orientation_is_separate") is not True):
            raise InputError("Identify coordinates must be declared 1-based inclusive input-forward")
        boundary = _mapping(frame.get("zero_width_boundaries"), "zero_width_boundaries")
        if (type(boundary.get("indexing")) is not int or boundary["indexing"] != 0
                or boundary.get("representation") != "interbase"
                or boundary.get("reference") != "input_forward_strand"):
            raise InputError("Identify zero-width boundary frame is unsupported")
        binding = _mapping(report.get("source"), "source")
        source_id = checked_id(binding.get("id"), "Identify source ID")
        if source_id not in sources:
            raise InputError(f"Identify source {source_id} not present in FASTA")
        if source_id in native_sources:
            raise InputError(f"duplicate Identify report for source {source_id}")
        source = sources[source_id]
        if (type(binding.get("length")) is not int or binding["length"] != source.length
                or binding.get("sha256") != source.sha256
                or binding.get("sha256_scope") != "uppercase_sequence"):
            raise InputError("Identify source length/SHA256/scope does not match FASTA")
        _validation(report.get("validation"), "report validation")
        raw_arrays = report.get("arrays")
        if not isinstance(raw_arrays, list) or len(raw_arrays) > 10000:
            raise InputError("Identify arrays must be a bounded list (maximum 10000 per source)")
        validation = report["validation"]
        if (type(validation.get("validated_array_count")) is not int
                or validation["validated_array_count"] != len(raw_arrays)
                or type(validation.get("invalid_array_count")) is not int
                or validation["invalid_array_count"] != 0):
            raise InputError("Identify report validation counts disagree with arrays")
        summary = _mapping(report.get("summary"), "summary")
        counts = {category: 0 for category in _CATEGORIES}
        excluded = []
        original_ids: set[str] = set()
        for raw_array in raw_arrays:
            array = _array(raw_array, source)
            original_id = raw_array["id"]
            if original_id in original_ids:
                raise InputError("duplicate Identify array ID within source report")
            original_ids.add(original_id)
            counts[array.category] += 1  # type: ignore[index]
            if array.category in _ACCEPTED:
                if len(arrays) >= 10000:
                    raise InputError("combined Identify reports exceed 10000 accepted arrays")
                arrays.append(json_record(array))
            else:
                excluded.append({"id": original_id, "category": array.category})
        if (type(summary.get("array_count")) is not int or summary["array_count"] != len(raw_arrays)
                or type(summary.get("accepted_array_count")) is not int
                or summary["accepted_array_count"] != sum(counts[k] for k in _ACCEPTED)):
            raise InputError("Identify summary array counts disagree with canonical rows")
        stated_counts = _mapping(summary.get("category_counts"), "category_counts")
        if (set(stated_counts) - _CATEGORIES
                or any(type(v) is not int or v != counts[k] for k, v in stated_counts.items())
                or any(counts[k] and k not in stated_counts for k in counts)):
            raise InputError("Identify category counts disagree with canonical rows")
        topology, provenance = _topology(report, source)
        native_sources[source_id] = {
            "id": source_id, "length": source.length, "sha256": source.sha256,
            "topology": topology, "topology_provenance": provenance,
            "metadata": {
                "identify_report_metadata": {k: v for k, v in report.items() if k != "arrays"},
                "adapter": "crisprleader.identify.v1",
                "accepted_categories": sorted(_ACCEPTED), "excluded_arrays": excluded,
                "report_semantic_sha256": report_digest,
                "adapter_validation_scope": (
                    "source identity, core array/repeat/spacer reconstruction, declared orientation "
                    "and topology; other upstream evidence is preserved without scientific revalidation"
                ),
            },
        }
    for source_id, source in sources.items():
        if source_id not in native_sources:
            native_sources[source_id] = {
                "id": source_id, "length": source.length, "sha256": source.sha256,
                "topology": "unknown", "topology_provenance": None,
                "metadata": {"identify_report_status": "not_supplied"},
            }
    document = {"schema_version": "1.0.0", "coordinate_system": "zero_based_half_open",
                "sources": [native_sources[key] for key in sources], "arrays": arrays}
    load_annotation_document(document, sources)
    return document
