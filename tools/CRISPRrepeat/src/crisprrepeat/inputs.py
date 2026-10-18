"""Bounded, independently usable source reconstruction; no detector imports."""

from __future__ import annotations

import hashlib
from itertools import islice
import json
from pathlib import Path
import re

DNA = frozenset("ACGTRYSWKMBDHVN")
RNA = frozenset("ACGURYSWKMBDHVN")
MAX_FILE_BYTES = 100_000_000
MAX_REPORT_BYTES = 20_000_000
MAX_RECORDS = 1000
MAX_REPEAT_LENGTH = 200
MAX_CONTEXT_LENGTH = 5000
MAX_TOTAL_CONTEXT_BASES = 20000
MAX_INSTANCES = 1000
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def sha(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def bounded_paths(iterator, label):
    paths = list(islice(iterator, MAX_RECORDS + 1))
    if len(paths) > MAX_RECORDS:
        raise ValueError(f"Too many {label}")
    return sorted(paths)


def require_object(value, label):
    if not isinstance(value, dict):
        raise ValueError(f"Expected {label} object")
    return value


def bounded_bytes(path: Path, maximum=MAX_FILE_BYTES) -> bytes:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > maximum:
        raise ValueError("Missing, linked or oversized input file")
    with path.open("rb") as stream:
        data = stream.read(maximum + 1)
    if len(data) > maximum:
        raise ValueError("Input exceeds byte limit")
    return data


def bounded_text(path: Path, maximum=MAX_FILE_BYTES) -> str:
    return bounded_bytes(path, maximum).decode("utf-8")


def parse_json(text: str) -> dict:
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError("Duplicate JSON key")
            result[key] = value
        return result

    value = json.loads(
        text,
        object_pairs_hook=pairs,
        parse_constant=lambda _: (_ for _ in ()).throw(ValueError("Nonfinite JSON")),
    )
    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object")
    return value


def read_json(path: Path) -> dict:
    return parse_json(bounded_text(path, MAX_REPORT_BYTES))


def read_fasta(path: Path, molecule="DNA", *, repeats=False) -> dict[str, str]:
    if molecule not in {"DNA", "RNA"}:
        raise ValueError("Molecule must be DNA or RNA")
    paths = bounded_paths(path.iterdir(), "FASTA files") if path.is_dir() and not path.is_symlink() else [path]
    records = {}
    total = 0
    alphabet = DNA if molecule == "DNA" else RNA
    literal_alphabet = alphabet | {v.lower() for v in alphabet}
    for file in paths:
        name, parts, record_length = None, [], 0

        def finish():
            nonlocal total
            if name is None:
                return
            sequence = "".join(parts)
            if not sequence or set(sequence) - alphabet:
                raise ValueError("Empty sequence or invalid molecule alphabet; mixed T/U is not accepted")
            if repeats and len(sequence) > MAX_REPEAT_LENGTH:
                raise ValueError("Repeat exceeds 200 nucleotides")
            total += len(sequence)
            if total > MAX_FILE_BYTES or len(records) >= MAX_RECORDS:
                raise ValueError("FASTA input limit exceeded")
            records[name] = sequence

        for line in bounded_text(file).splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:]
                if len(header) > 500 or not header or any(ord(x) < 32 for x in header):
                    raise ValueError("Invalid FASTA header")
                name = header.split()[0]
                if not re.fullmatch(r"[A-Za-z0-9_.:|+-]{1,200}", name) or name in records:
                    raise ValueError("Unsafe or duplicate source identifier")
                parts = []
                record_length = 0
            else:
                if name is None:
                    raise ValueError("Input must be FASTA")
                literal = "".join(line.split())
                if set(literal) - literal_alphabet:
                    raise ValueError("Invalid nucleotide symbol")
                record_length += len(literal)
                if repeats and record_length > MAX_REPEAT_LENGTH:
                    raise ValueError("Repeat exceeds 200 nucleotides")
                if total + record_length > MAX_FILE_BYTES:
                    raise ValueError("FASTA input limit exceeded")
                parts.append(literal.upper())
        finish()
    if not records:
        raise ValueError("No FASTA records")
    return records


def interval(raw, sequence):
    require_object(raw, "detector interval")
    start, end = raw.get("start"), raw.get("end")
    if type(start) is not int or type(end) is not int or not 1 <= start <= end <= len(sequence):
        raise ValueError("Invalid detector interval (expected one-based inclusive)")
    length = raw.get("length", end - start + 1)
    if type(length) is not int or length != end - start + 1:
        raise ValueError("Interval length mismatch")
    return start - 1, end


def load_arrays(path: Path, sources: dict, categories: list[str]) -> list[dict]:
    paths = bounded_paths(path.rglob("report.json"), "reports") if path.is_dir() and not path.is_symlink() else [path]
    arrays, seen_sources, identifiers = [], set(), set()
    total_bases = instances = 0
    for report_path in paths:
        if any(parent.is_symlink() for parent in [report_path, *report_path.parents]):
            raise ValueError("Linked detector report")
        report = read_json(report_path)
        schema = require_object(report.get("schema"), "detector schema")
        if schema.get("name") != "CRISPRidentify-report" or not str(schema.get("version", "")).startswith(
            "1."
        ):
            raise ValueError("Expected canonical CRISPRidentify report schema 1.x")
        source = require_object(report.get("source"), "detector source")
        source_id = source.get("id")
        if not isinstance(source_id, str):
            raise ValueError("Expected detector source identifier")
        sequence = sources.get(source_id)
        if (
            sequence is None
            or source_id in seen_sources
            or type(source.get("length")) is not int
            or source.get("length") != len(sequence)
            or source.get("sha256") != sha(sequence)
            or source.get("sha256_scope") != "uppercase_sequence"
        ):
            raise ValueError("Detector source identity mismatch")
        seen_sources.add(source_id)
        raw_arrays = report.get("arrays")
        if not isinstance(raw_arrays, list):
            raise ValueError("Missing canonical arrays")
        for raw in raw_arrays:
            require_object(raw, "detector array")
            if raw.get("category") not in categories:
                continue
            identifier = raw.get("id")
            if not isinstance(identifier, str) or not identifier or len(identifier) > 200:
                raise ValueError("Missing array identifier")
            key = (source_id, identifier)
            if key in identifiers:
                raise ValueError("Duplicate array identity")
            identifiers.add(key)
            start, end = interval(raw.get("source_interval", {}), sequence)
            if raw.get("source_sequence") != sequence[start:end]:
                raise ValueError("Array does not reconstruct from submitted sequence")
            total_bases += end - start
            if end - start > MAX_CONTEXT_LENGTH or total_bases > MAX_TOTAL_CONTEXT_BASES:
                raise ValueError("Context length budget exceeded; context was not shortened")
            repeats, previous = [], start
            raw_repeats = raw.get("repeats")
            if not isinstance(raw_repeats, list):
                raise ValueError("Expected detector repeats list")
            for index, repeat in enumerate(raw_repeats, 1):
                require_object(repeat, "detector repeat")
                rs, re = interval(repeat.get("source_interval", {}), sequence)
                if (not start <= rs < re <= end or rs < previous
                        or type(repeat.get("index")) is not int or repeat.get("index") != index):
                    raise ValueError("Invalid repeat ordering or containment")
                if repeat.get("sequence") != sequence[rs:re] or re - rs > MAX_REPEAT_LENGTH:
                    raise ValueError("Repeat does not reconstruct or exceeds length limit")
                repeats.append({"index": index, "start": rs, "end": re, "sequence": sequence[rs:re]})
                previous = re
            instances += len(repeats)
            if (not repeats or type(raw.get("repeat_count")) is not int
                    or raw.get("repeat_count") != len(repeats) or instances > MAX_INSTANCES):
                raise ValueError("Invalid or excessive repeat count")
            consensus = raw.get("consensus_repeat")
            if (
                not isinstance(consensus, str)
                or not 1 <= len(consensus) <= MAX_REPEAT_LENGTH
                or set(consensus) - DNA
            ):
                raise ValueError("Invalid consensus")
            orientation = raw.get("orientation")
            if orientation is not None:
                require_object(orientation, "detector orientation")
            strand = (orientation or {}).get("strand")
            if strand is not None and not isinstance(strand, str):
                raise ValueError("Invalid detector orientation strand")
            arrays.append(
                {
                    "id": f"array_{len(arrays) + 1:06d}",
                    "native_id": identifier,
                    "source_id": source_id,
                    "start": start,
                    "end": end,
                    "category": raw["category"],
                    "source_sequence": sequence[start:end],
                    "consensus_repeat": consensus,
                    "repeats": repeats,
                    "detector_strand": strand if strand in {"+", "-"} else None,
                }
            )
    if seen_sources != set(sources):
        raise ValueError("Detector reports do not cover submitted records")
    return arrays
