"""Independent output checks; deliberately does not import coordinate helpers."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import stat
from typing import Any, Mapping

from .errors import InputError
from .models import SourceRecord

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_MAX_REPORT_BYTES = 128 * 1024 * 1024


def _unique_pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in items:
        if key in result:
            raise InputError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def read_json(path: str | Path, limit: int = _MAX_REPORT_BYTES) -> dict[str, Any]:
    """Bounded strict JSON object parsing, including duplicate and nonfinite checks."""
    if type(limit) is not int or limit < 1:
        raise InputError("JSON byte budget must be a positive integer")
    try:
        descriptor = os.open(Path(path), os.O_RDONLY | os.O_NONBLOCK)
        with os.fdopen(descriptor, "rb") as handle:
            metadata = os.fstat(handle.fileno())
            if not stat.S_ISREG(metadata.st_mode):
                raise InputError("JSON must be a regular file")
            if metadata.st_size > limit:
                raise InputError("JSON exceeds the configured byte limit")
            data = handle.read(limit + 1)
        if len(data) > limit:
            raise InputError("JSON exceeds the configured byte limit")
        result = json.loads(
            data,
            object_pairs_hook=_unique_pairs,
            parse_constant=lambda value: (_ for _ in ()).throw(
                InputError(f"Non-finite JSON value: {value}")
            ),
        )
    except (OSError, ValueError, RecursionError) as error:
        raise InputError(f"Cannot read JSON: {error}") from error
    if not isinstance(result, dict):
        raise InputError("JSON root must be an object")
    _finite_tree(result)
    return result


def _finite_tree(value: Any, depth: int = 0) -> None:
    if depth > 48:
        raise InputError("Scientific output exceeds maximum nesting")
    if isinstance(value, float) and not math.isfinite(value):
        raise InputError("Scientific output contains a non-finite number")
    if isinstance(value, dict):
        for item in value.values():
            _finite_tree(item, depth + 1)
    elif isinstance(value, list):
        for item in value:
            _finite_tree(item, depth + 1)


def _integer(value: Any, name: str) -> int:
    if type(value) is not int:
        raise InputError(f"{name} must be an integer")
    return value


def _identifier(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value or len(value) > 512 or any(
        ord(character) < 33 or ord(character) == 127 for character in value
    ):
        raise InputError(f"{name} must be a nonempty, bounded identifier without whitespace")
    return value


def _rows(value: Any, name: str, maximum: int) -> list[dict[str, Any]]:
    if not isinstance(value, list) or len(value) > maximum or any(
        not isinstance(row, dict) for row in value
    ):
        raise InputError(f"{name} must be a bounded list of objects")
    return value


def _reconstruct(sequence: str, segments: Any, strand: Any) -> str:
    if not isinstance(segments, list) or len(segments) > 2:
        raise InputError("Feature segments must be a list of at most two intervals")
    parts: list[str] = []
    intervals: list[tuple[int, int]] = []
    for segment in segments:
        if not isinstance(segment, dict) or set(segment) != {"start", "end"}:
            raise InputError("Malformed feature segment")
        start = _integer(segment["start"], "start")
        end = _integer(segment["end"], "end")
        if not 0 <= start < end <= len(sequence):
            raise InputError("Feature interval lies outside its source")
        if any(max(start, a) < min(end, b) for a, b in intervals):
            raise InputError("Feature segments overlap")
        intervals.append((start, end))
        parts.append(sequence[start:end])
    if len(intervals) == 2 and not (
        intervals[0][1] == len(sequence) and intervals[1][0] == 0
    ):
        raise InputError("Two segments must describe a contiguous origin-crossing arc")
    text = "".join(parts)
    if strand == "-":
        text = text.translate(_COMPLEMENT)[::-1]
    elif strand != "+":
        raise InputError("Coordinate-bound features require an explicit strand")
    return text


def _validate_feature(feature: Any, sources: Mapping[str, SourceRecord]) -> dict[str, Any]:
    if not isinstance(feature, dict):
        raise InputError("Feature must be an object")
    source = sources.get(str(feature.get("source_id", "")))
    if source is None:
        raise InputError("Feature refers to an unknown source")
    segments = feature.get("segments")
    if isinstance(segments, list) and len(segments) == 2 and source.topology != "circular":
        raise InputError("Wrapped feature lacks explicit circular topology")
    reconstructed = _reconstruct(source.sequence, segments, feature.get("strand"))
    if feature.get("sequence") != reconstructed:
        raise InputError("Feature sequence does not match its source intervals/strand")
    return feature


def _offset_segments(candidate: dict[str, Any], start: int, end: int) -> list[dict[str, int]]:
    """Independent projection of an oriented flank slice, without predictor helpers."""
    if not 0 <= start < end <= len(candidate["sequence"]):
        raise InputError("Core local interval lies outside the candidate flank")
    length = len(candidate["sequence"])
    lo, hi = (length - end, length - start) if candidate["strand"] == "-" else (start, end)
    cursor = 0
    output: list[dict[str, int]] = []
    for segment in candidate["segments"]:
        width = segment["end"] - segment["start"]
        left, right = max(lo, cursor), min(hi, cursor + width)
        if left < right:
            output.append({"start": segment["start"] + left - cursor,
                           "end": segment["start"] + right - cursor})
        cursor += width
    return output


def validate_result(document: dict[str, Any], sources: Mapping[str, SourceRecord]) -> None:
    """Verify identities, source reconstruction and evidence binding on a full result."""
    if not isinstance(document, dict):
        raise InputError("Result must be an object")
    _finite_tree(document)
    if document.get("schema_version") != "1.0.0":
        raise InputError("Unsupported result schema_version")
    if document.get("coordinate_system") != "zero_based_half_open":
        raise InputError("Unsupported result coordinate system")
    if document.get("execution_status") != "completed":
        raise InputError("Only completed scientific results may be exported")
    mode = document.get("analysis_mode")
    if mode not in ("context_extraction", "prediction"):
        raise InputError("Unknown analysis_mode")
    config = document.get("effective_config")
    if not isinstance(config, dict):
        raise InputError("Result effective_config must be an object")
    flank_length = _integer(config.get("flank_length"), "effective_config.flank_length")
    if not 1 <= flank_length <= 5000:
        raise InputError("effective_config.flank_length must be between 1 and 5000")
    source_rows = _rows(document.get("sources"), "sources", 100000)
    if len(source_rows) != len(sources):
        raise InputError("Result source count differs from supplied FASTA")
    seen: set[str] = set()
    for row in source_rows:
        identity = _identifier(row.get("id"), "source id")
        if identity not in sources or identity in seen:
            raise InputError("Duplicate or unknown result source")
        seen.add(row["id"])
        source = sources[row["id"]]
        digest = hashlib.sha256(source.sequence.upper().encode("ascii")).hexdigest()
        if _integer(row.get("length"), "source length") != len(source.sequence) or row.get("sha256") != digest:
            raise InputError("Result source checksum/length differs from supplied FASTA")
        if row.get("topology") != source.topology:
            raise InputError("Result topology differs from validated source metadata")
        if source.topology not in ("linear", "circular", "unknown"):
            raise InputError("Unknown source topology")
        if source.topology == "circular" and (
            not isinstance(row.get("topology_provenance"), str)
            or not row["topology_provenance"].strip()
        ):
            raise InputError("Circular source requires topology provenance")
        if row.get("topology_provenance") != source.topology_provenance:
            raise InputError("Result topology provenance differs from source metadata")
    arrays = _rows(document.get("arrays"), "arrays", 10000)
    candidates = _rows(document.get("candidates"), "candidates", 20000)
    predictions = _rows(document.get("predictions"), "predictions", 10000)
    if len(predictions) != (len(arrays) if mode == "prediction" else 0):
        raise InputError("Prediction count does not match the declared analysis mode")
    if mode == "context_extraction" and document.get("reference") is not None:
        raise InputError("Context extraction cannot assert a prediction reference")
    array_ids: set[tuple[str, str]] = set()
    arrays_by_id: dict[tuple[str, str], dict[str, Any]] = {}
    global_array_ids: set[str] = set()
    for array in arrays:
        key = (_identifier(array.get("source_id"), "array source"),
               _identifier(array.get("id"), "array id"))
        if key[1] in global_array_ids or key[0] not in sources:
            raise InputError("Duplicate or unknown array source identity")
        global_array_ids.add(key[1])
        array_ids.add(key)
        arrays_by_id[key] = array
        source = sources[key[0]]
        sequence = _reconstruct(source.sequence, array.get("segments"), "+")
        if not sequence or len(array["segments"]) == 2 and source.topology != "circular":
            raise InputError("Array is empty or wraps a non-circular source")
        repeats = _rows(array.get("terminal_repeats"), "terminal_repeats", 2)
        if len(repeats) != 2:
            raise InputError("Array must retain two observed terminal repeats")
        for repeat in repeats:
            if _reconstruct(source.sequence, repeat.get("segments"), "+") != repeat.get("sequence"):
                raise InputError("Array terminal repeat differs from original source")
            if not repeat.get("sequence"):
                raise InputError("Array terminal repeats must be nonempty")
        lengths = [len(repeat["sequence"]) for repeat in repeats]
        if sum(lengths) > len(sequence):
            raise InputError("Array terminal repeats overlap")
        frame = {"segments": array["segments"], "sequence": sequence, "strand": "+"}
        if (repeats[0]["segments"] != _offset_segments(frame, 0, lengths[0])
                or repeats[1]["segments"] != _offset_segments(frame, len(sequence) - lengths[1], len(sequence))):
            raise InputError("Terminal-repeat anchors must be flush with array endpoints")
        consensus = array.get("consensus_repeat")
        if not isinstance(consensus, str) or not consensus or set(consensus) - set("ACGTRYSWKMBDHVN"):
            raise InputError("Array consensus must be uppercase IUPAC DNA")
        if array.get("strand") not in ("+", "-", "?") or type(array.get("strand_constraint")) is not bool:
            raise InputError("Array strand metadata is invalid")
        if array["strand_constraint"] and array["strand"] == "?":
            raise InputError("Unknown strand cannot be a constraint")
    by_id: dict[str, dict[str, Any]] = {}
    sides: set[tuple[str, str, str]] = set()
    for candidate in candidates:
        _validate_feature(candidate, sources)
        key = (_identifier(candidate.get("source_id"), "candidate source"),
               _identifier(candidate.get("array_id"), "candidate array"))
        identity = _identifier(candidate.get("id"), "candidate id")
        if key not in array_ids or identity in by_id:
            raise InputError("Duplicate candidate or unbound array identity")
        if candidate.get("side") not in ("left", "right"):
            raise InputError("Candidate side is invalid")
        side_key = (*key, candidate["side"])
        if side_key in sides:
            raise InputError("Duplicate candidate side")
        sides.add(side_key)
        if candidate["strand"] != ("+" if candidate["side"] == "left" else "-"):
            raise InputError("Candidate side/strand mismatch")
        requested = _integer(candidate.get("requested_length"), "requested_length")
        if not 1 <= requested <= 5000 or len(candidate["sequence"]) > requested:
            raise InputError("Invalid candidate window length")
        if requested != flank_length:
            raise InputError("Candidate window length differs from effective_config.flank_length")
        anchor = candidate.get("anchor", {})
        expected_anchor = arrays_by_id[key]["terminal_repeats"][0 if candidate["side"] == "left" else 1]
        if anchor != expected_anchor:
            raise InputError("Candidate anchor is not the corresponding observed array terminal repeat")
        source = sources[candidate["source_id"]]
        metadata = candidate.get("metadata")
        if not isinstance(metadata, dict):
            raise InputError("Candidate metadata must be an object")
        if (
            metadata.get("source_sha256") != source.sha256
            or _integer(metadata.get("source_length"), "candidate metadata source_length") != source.length
            or metadata.get("topology") != source.topology
        ):
            raise InputError("Candidate source metadata differs from its validated source")
        if (
            metadata.get("extraction_policy") != "adjacent-window-excluding-focal-array"
            or metadata.get("window_is_biological_boundary") is not False
        ):
            raise InputError("Candidate context metadata contradicts its extraction policy")
        array = arrays_by_id[key]
        array_size = sum(segment["end"] - segment["start"] for segment in array["segments"])
        if source.topology == "circular":
            available = source.length - array_size
            reason = "circular_array_boundary"
        else:
            available = (array["segments"][0]["start"] if candidate["side"] == "left"
                         else source.length - array["segments"][-1]["end"])
            reason = "contig_boundary" if source.topology == "linear" else "unknown_topology_boundary"
        length = min(available, requested)
        start = ((array["segments"][0]["start"] - length if candidate["side"] == "left"
                  else array["segments"][-1]["end"]) % source.length)
        expected = [] if length == 0 else [{"start": start, "end": min(source.length, start + length)}]
        if start + length > source.length:
            expected.append({"start": 0, "end": start + length - source.length})
        if candidate["segments"] != expected:
            raise InputError("Candidate flank is not the requested adjacent window outside the array")
        if candidate.get("truncated_reason") != (reason if length < requested else None):
            raise InputError("Candidate clipping reason differs from its geometry")
        consensus = array["consensus_repeat"]
        if candidate["strand"] == "-":
            consensus = consensus.translate(_COMPLEMENT)[::-1]
        if candidate.get("consensus_repeat") != consensus:
            raise InputError("Candidate consensus is not in its declared strand frame")
        observed = _reconstruct(source.sequence, anchor.get("segments"), candidate["strand"])
        if candidate.get("observed_repeat") != observed:
            raise InputError("Observed terminal-repeat sequence mismatch")
        by_id[candidate["id"]] = candidate
    if len(candidates) != 2 * len(arrays):
        raise InputError("Every array must retain exactly two candidate flanks")
    predicted: set[tuple[str, str]] = set()
    for prediction in predictions:
        key = (_identifier(prediction.get("source_id"), "prediction source"),
               _identifier(prediction.get("array_id"), "prediction array"))
        if key not in array_ids or key in predicted:
            raise InputError("Prediction identity does not match exactly one array")
        predicted.add(key)
        if prediction.get("status") not in ("supported", "ambiguous", "no_supported_hit", "unsupported"):
            raise InputError("Unknown scientific decision status")
        evidence_rows = _rows(prediction.get("evidence"), "evidence", 100000)
        motif_rows = _rows(prediction.get("motifs"), "motifs", 1000000)
        for evidence in evidence_rows:
            evidence_candidate = by_id.get(_identifier(evidence.get("candidate_id"), "evidence candidate"))
            if evidence_candidate is None or (evidence_candidate["source_id"], evidence_candidate["array_id"]) != key:
                raise InputError("Evidence refers to an unbound candidate")
            start, end = evidence.get("start"), evidence.get("end")
            if (start is None) != (end is None):
                raise InputError("Evidence interval is incomplete")
            if start is not None:
                _offset_segments(evidence_candidate, _integer(start, "start"), _integer(end, "end"))
        core = prediction.get("core")
        if prediction["status"] == "supported":
            selected = by_id.get(_identifier(prediction.get("selected_candidate_id"), "selected candidate"))
            if selected is None or (selected["source_id"], selected["array_id"]) != key:
                raise InputError("Supported decision requires a bound selected candidate")
            core = _validate_feature(core, sources)
            if (core.get("source_id"), core.get("array_id")) != key:
                raise InputError("Core feature identity differs from its prediction")
            start = _integer(prediction.get("core_start"), "core_start")
            end = _integer(prediction.get("core_end"), "core_end")
            if core["segments"] != _offset_segments(selected, start, end):
                raise InputError("Core does not match selected candidate-local interval")
            if core["sequence"] != selected["sequence"][start:end] or core["strand"] != selected["strand"]:
                raise InputError("Core sequence/strand differs from selected candidate")
            if not any(item.get("passed") is True and item.get("candidate_id") == selected["id"]
                       and item.get("start") == start and item.get("end") == end
                       for item in prediction.get("evidence", [])):
                raise InputError("Supported core lacks matching passed raw evidence")
        elif any(prediction.get(field) is not None for field in ("core", "core_start", "core_end", "selected_candidate_id")):
            raise InputError("Unsupported/ambiguous decisions cannot select a supported core")
        for hit in motif_rows:
            motif_candidate = by_id.get(_identifier(hit.get("candidate_id"), "motif candidate"))
            if motif_candidate is None or (motif_candidate["source_id"], motif_candidate["array_id"]) != key:
                raise InputError("Motif refers to an unbound candidate")
            start, end = _integer(hit.get("start"), "motif start"), _integer(hit.get("end"), "motif end")
            junction = motif_candidate["sequence"] + motif_candidate["observed_repeat"]
            if not 0 <= start < end <= len(junction):
                raise InputError("Motif lies outside the observed junction context")
            boundary = len(motif_candidate["sequence"])
            if not ((hit.get("domain") == "leader" and end <= boundary)
                    or (hit.get("domain") == "repeat" and start >= boundary)
                    or (hit.get("domain") == "junction" and start < boundary < end)):
                raise InputError("Motif interval does not match its declared domain")
            feature = hit.get("feature")
            feature = _validate_feature(feature, sources)
            if (feature.get("source_id"), feature.get("array_id")) != key:
                raise InputError("Motif feature identity differs from its prediction")
            sequence = junction[start:end]
            if hit.get("strand") == "-":
                sequence = sequence.translate(_COMPLEMENT)[::-1]
            elif hit.get("strand") != "+":
                raise InputError("Motif orientation is invalid")
            expected_strand = motif_candidate["strand"]
            if hit["strand"] == "-":
                expected_strand = "-" if expected_strand == "+" else "+"
            if feature["strand"] != expected_strand:
                raise InputError("Motif genome strand differs from hypothesis and motif strands")
            combined = (motif_candidate["segments"] + motif_candidate["anchor"]["segments"]
                        if motif_candidate["strand"] == "+"
                        else motif_candidate["anchor"]["segments"] + motif_candidate["segments"])
            expected = _offset_segments({"segments": combined, "strand": motif_candidate["strand"],
                                         "sequence": junction}, start, end)
            merged: list[dict[str, int]] = []
            for segment in expected:
                if merged and merged[-1]["end"] == segment["start"]:
                    merged[-1]["end"] = segment["end"]
                else:
                    merged.append(dict(segment))
            if feature["segments"] != merged:
                raise InputError("Motif genomic coordinates differ from its local interval")
            if sequence != feature["sequence"]:
                raise InputError("Motif sequence differs from its observed local context")
