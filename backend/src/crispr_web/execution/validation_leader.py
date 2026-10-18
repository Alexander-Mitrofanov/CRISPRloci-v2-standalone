"""Reconstruct leader contexts independently of Integration and CRISPRleader.

The server accepts the pinned context-only, unknown-topology contract. Uploaded
FASTA does not establish circular topology. Scientific prediction stays absent.
"""

from __future__ import annotations

import hashlib
import json
import re

from .validation_annotations import _COMPLEMENT, _fasta, _json, _payload
from .validation_identify import _validate_identify_reports


def _hash(value):
    return hashlib.sha256(value.encode("ascii")).hexdigest()


def _segment(interval):
    return {"start": interval["start"] - 1, "end": interval["end"]}


def _array(raw, identity, sequence):
    interval = raw["source_interval"]
    start, end = interval["start"], interval["end"]
    if (
        type(start) is not int
        or type(end) is not int
        or not 1 <= start <= end <= len(sequence)
        or interval.get("length") != end - start + 1
    ):
        raise ValueError("Invalid source-bound array interval")
    repeats = raw["repeats"]
    if len(repeats) < 2 or len(repeats) != raw["repeat_count"]:
        raise ValueError("Invalid terminal repeats")
    previous = start - 1
    for number, repeat in enumerate(repeats, 1):
        span = repeat["source_interval"]
        low, high = span["start"], span["end"]
        if (
            type(low) is not int
            or type(high) is not int
            or not start <= low <= high <= end
            or low <= previous
            or span["length"] != high - low + 1
            or repeat["index"] != number
            or repeat["sequence"].upper() != sequence[low - 1 : high]
        ):
            raise ValueError("Repeat evidence differs from submitted sequence")
        previous = high
    if (
        repeats[0]["source_interval"]["start"] != start
        or repeats[-1]["source_interval"]["end"] != end
    ):
        raise ValueError("Terminal anchors do not match array endpoints")
    if (
        "source_sequence" in raw
        and raw["source_sequence"].upper() != sequence[start - 1 : end]
    ):
        raise ValueError("Array sequence differs from input")
    key = hashlib.sha256(
        json.dumps([identity, raw["id"]], separators=(",", ":")).encode()
    ).hexdigest()
    repeats = raw["repeats"]
    orientation = raw["orientation"]
    return {
        "id": "identify:" + key,
        "source_id": identity,
        "segments": [_segment(raw["source_interval"])],
        "terminal_repeats": [
            {
                "segments": [_segment(row["source_interval"])],
                "sequence": row["sequence"].upper(),
            }
            for row in (repeats[0], repeats[-1])
        ],
        "consensus_repeat": raw["consensus_repeat"].upper(),
        "strand": orientation["strand"] if orientation["status"] == "computed" else "?",
        "strand_constraint": False,
        "category": raw["category"],
        "metadata": {
            "identify_original": raw,
            "sequence_frame": "input_forward_strand",
            "orientation_is_supplied_evidence": True,
        },
    }


def _contexts(array, sequence, length):
    interval = array["segments"][0]
    rows = []
    for side, strand, anchor in [
        ("left", "+", array["terminal_repeats"][0]),
        ("right", "-", array["terminal_repeats"][1]),
    ]:
        start, end = (
            (max(0, interval["start"] - length), interval["start"])
            if side == "left"
            else (interval["end"], min(len(sequence), interval["end"] + length))
        )
        orient = (
            (lambda value: value)
            if strand == "+"
            else (lambda value: value.translate(_COMPLEMENT)[::-1])
        )
        rows.append(
            {
                "id": array["id"] + ":" + side,
                "source_id": array["source_id"],
                "array_id": array["id"],
                "side": side,
                "strand": strand,
                "segments": [{"start": start, "end": end}] if start < end else [],
                "sequence": orient(sequence[start:end]),
                "requested_length": length,
                "truncated_reason": "unknown_topology_boundary"
                if end - start < length
                else None,
                "anchor": anchor,
                "observed_repeat": orient(anchor["sequence"]),
                "consensus_repeat": orient(array["consensus_repeat"]),
                "metadata": {
                    "topology": "unknown",
                    "source_sha256": _hash(sequence),
                    "source_length": len(sequence),
                    "extraction_policy": "adjacent-window-excluding-focal-array",
                    "window_is_biological_boundary": False,
                },
            }
        )
    return rows


def validate_leader_pipeline(
    claimed, root, expected, integration, input_root, checkpoint=None
):
    if (
        input_root is None
        or integration.get("status") != "completed"
        or integration.get("workflow") != claimed.mode
    ):
        return False
    try:
        sources = {}
        for path in sorted(input_root.iterdir()):
            identity, sequence = _fasta(_payload(input_root, path, checkpoint))
            if identity in sources or expected.get(identity) != (
                len(sequence),
                _hash(sequence),
            ):
                return False
            sources[identity] = sequence
        if set(sources) != set(expected):
            return False
        index = _json(root, root / "leader/index.json", checkpoint)
        length = claimed.leader_flank_length
        if (
            type(length) is not int
            or not 1 <= length <= 5000
            or index.get("schema_version") != "1.0.0"
            or index.get("analysis_mode") != "context_extraction"
            or index.get("flank_length") != length
            or len(index["records"]) != len(sources)
        ):
            return False
        if claimed.mode == "leader":
            if integration["stages"]["crispridentify"].get("status") != "completed":
                return False
            if not _validate_identify_reports(root, expected, checkpoint):
                return False
            native_mapping = {identity: identity for identity in sources}
        else:
            tracr = _json(root, root / "tracrrna/report.json", checkpoint)
            native_mapping = {
                row["source"]["id"]: row["key"] for row in tracr["records"]
            }
        seen, keys = set(), set()
        total = 0
        for record in index["records"]:
            if checkpoint:
                checkpoint()
            identity, native_id, key = (
                record["source_id"],
                record["native_source_id"],
                record["key"],
            )
            if (
                identity not in sources
                or identity in seen
                or not re.fullmatch(r"record-[0-9]{6}", key)
                or key in keys
                or native_id != native_mapping.get(identity)
            ):
                return False
            seen.add(identity)
            keys.add(key)
            sequence = sources[identity]
            reports_root = root / "crispridentify"
            if claimed.mode == "loci":
                if not re.fullmatch(r"f[0-9]{6}\.s[0-9]{6}", native_id):
                    return False
                work = root / "tracrrna/work" / native_id
                if _fasta(_payload(root, work / "input.fasta", checkpoint)) != (
                    native_id,
                    sequence,
                ):
                    return False
                reports_root = work / "identify"
                if not _validate_identify_reports(
                    root,
                    {native_id: expected[identity]},
                    checkpoint,
                    reports_root=reports_root,
                ):
                    return False
            reports = []
            for path in reports_root.rglob("report.json"):
                document = _json(root, path, checkpoint)
                if document.get("source", {}).get("id") == native_id:
                    reports.append(document)
            if len(reports) != 1:
                return False
            raw = _payload(root, root / "leader" / key / "result.json", checkpoint)
            manifest = _json(root, root / "leader" / key / "manifest.json", checkpoint)
            if manifest.get("schema_version") != "1.0.0" or manifest["files"][
                "result.json"
            ] != {
                "sha256": hashlib.sha256(raw).hexdigest(),
                "size_bytes": len(raw),
            }:
                return False
            report = _json(root, root / "leader" / key / "result.json", checkpoint)
            if (
                report.get("schema_version") != "1.0.0"
                or report.get("execution_status") != "completed"
                or report.get("analysis_mode") != "context_extraction"
                or report.get("predictions") != []
                or report.get("reference") is not None
                or report.get("coordinate_system") != "zero_based_half_open"
                or report.get("effective_config", {}).get("flank_length") != length
            ):
                return False
            if report["sources"] != [
                {
                    "id": native_id,
                    "length": len(sequence),
                    "sha256": _hash(sequence),
                    "topology": "unknown",
                    "topology_provenance": None,
                }
            ]:
                return False
            arrays = [
                _array(row, native_id, sequence)
                for row in reports[0]["arrays"]
                if row["category"] in {"Bona-fide", "Possible"}
            ]
            if json.dumps(report["arrays"], sort_keys=True) != json.dumps(
                arrays, sort_keys=True
            ) or json.dumps(report["candidates"], sort_keys=True) != json.dumps(
                [row for array in arrays for row in _contexts(array, sequence, length)],
                sort_keys=True,
            ):
                return False
            total += len(arrays)
        stage = integration["stages"]["crisprleader"]
        return (
            stage.get("status") == "completed"
            and stage.get("analysis_mode") == "context_extraction"
            and stage.get("prediction_performed") is False
            and stage.get("flank_length") == length
            and all(
                type(stage.get(key)) is int and stage[key] == value
                for key, value in {
                    "record_count": len(sources),
                    "array_count": total,
                    "context_count": total * 2,
                }.items()
            )
            and index.get("array_count") == total
            and index.get("context_count") == total * 2
        )
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        IndexError,
        AttributeError,
        RecursionError,
    ):
        return False
