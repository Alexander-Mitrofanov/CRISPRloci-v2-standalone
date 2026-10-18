"""Bounded input snapshots and deterministic, atomic artifact primitives."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import stat
import tempfile

MAX_BYTES = 64 * 1024 * 1024
MAX_RECORDS = 100_000
MAX_SEQUENCE = 200
MAX_BASES = 10_000_000
DNA = frozenset("ACGTRYSWKMBDHVN")


class InputError(ValueError):
    """Malformed input, unsupported schema, or inconsistent provenance."""


class LimitError(InputError):
    """A resource budget was exceeded; the analysis is not complete."""


def digest(data: bytes | str) -> str:
    return hashlib.sha256(data.encode("utf-8") if isinstance(data, str) else data).hexdigest()


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, ensure_ascii=True, allow_nan=False, separators=(",", ":")) + "\n").encode()


def identifier(value: object, label: str = "identifier") -> str:
    if not isinstance(value, str) or not 1 <= len(value) <= 500 or any(ord(c) < 33 or ord(c) > 126 for c in value):
        raise InputError(f"{label} must contain 1–500 printable ASCII characters without whitespace")
    return value


def integer(value: object, label: str, low: int, high: int) -> int:
    if type(value) is not int or not low <= value <= high:
        raise InputError(f"{label} must be an integer in [{low}, {high}]")
    return value


def regular(path: Path) -> Path:
    path = Path(path).absolute()
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise InputError("Input paths must not traverse symlinks")
    try:
        mode = path.stat().st_mode
    except OSError as exc:
        raise InputError(f"Cannot read input: {path.name}") from exc
    if not stat.S_ISREG(mode):
        raise InputError("Input must be a regular file")
    return path


def read_bytes(path: Path, limit: int = MAX_BYTES) -> bytes:
    path = regular(path)
    if path.stat().st_size > limit:
        raise LimitError("Input exceeds byte budget")
    with path.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise LimitError("Input exceeds byte budget")
    return data


def decode_json(data: bytes):
    def pairs(values):
        result = {}
        for key, value in values:
            if key in result:
                raise InputError("Duplicate JSON key")
            result[key] = value
        return result

    def constant(_):
        raise InputError("Nonfinite JSON number")

    def number(token, convert):
        if len(token) > 128:
            raise InputError("Oversized JSON numeric token")
        value = convert(token)
        if isinstance(value, float) and not math.isfinite(value):
            raise InputError("Nonfinite JSON number")
        return value

    try:
        return json.loads(data, object_pairs_hook=pairs, parse_constant=constant,
                          parse_int=lambda token: number(token, int),
                          parse_float=lambda token: number(token, float))
    except InputError:
        raise
    except (UnicodeError, ValueError, RecursionError) as exc:
        raise InputError("Malformed JSON") from exc


def read_json(path: Path):
    return decode_json(read_bytes(path))


def normalize_records(records: object) -> list[dict]:
    if not isinstance(records, list) or not 1 <= len(records) <= MAX_RECORDS:
        raise InputError(f"Expected 1–{MAX_RECORDS} sequence records")
    result, seen = [], set()
    total = 0
    for record in records:
        if not isinstance(record, dict):
            raise InputError("Sequence record must be an object")
        ident = identifier(record.get("id"))
        if ident in seen:
            raise InputError("Duplicate sequence identifier")
        seen.add(ident)
        sequence = record.get("sequence")
        if not isinstance(sequence, str):
            raise InputError("Sequence must be a string")
        sequence = sequence.upper().replace("U", "T")
        if not 1 <= len(sequence) <= MAX_SEQUENCE or not set(sequence) <= DNA:
            raise InputError(f"Sequences must contain 1–{MAX_SEQUENCE} IUPAC DNA/RNA letters")
        description = record.get("description", "")
        if (not isinstance(description, str) or len(description) > 2000
                or any(ord(c) < 32 or ord(c) > 126 for c in description)):
            raise InputError("Descriptions must be at most 2000 printable ASCII characters")
        total += len(sequence)
        if total > MAX_BASES:
            raise LimitError("Sequence collection exceeds total base budget")
        result.append({"id": ident, "sequence": sequence, "description": description})
    return result


def read_fasta(path: Path) -> tuple[list[dict], str]:
    """Read one bounded FASTA snapshot; normalize case and U to T explicitly."""
    raw = read_bytes(path)
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeError as exc:
        raise InputError("FASTA must contain ASCII text") from exc
    records, current, chunks = [], None, []
    sequence_length = 0
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current is not None:
                records.append({**current, "sequence": "".join(chunks)})
            parts = line[1:].split(maxsplit=1)
            if not parts:
                raise InputError("FASTA header is empty")
            current = {"id": parts[0], "description": parts[1] if len(parts) > 1 else ""}
            chunks = []
            sequence_length = 0
            if len(records) >= MAX_RECORDS:
                raise LimitError("Too many FASTA records")
        else:
            if current is None:
                raise InputError("Sequence precedes FASTA header")
            sequence_length += len(line)
            if sequence_length > MAX_SEQUENCE:
                raise LimitError("FASTA record exceeds sequence length budget")
            chunks.append(line)
    if current is not None:
        records.append({**current, "sequence": "".join(chunks)})
    return normalize_records(records), digest(raw)


def write_json(path: Path, value: object) -> None:
    data = canonical(value)
    if len(data) > MAX_BYTES:
        raise LimitError("JSON artifact exceeds byte budget")
    Path(path).write_bytes(data)


def atomic_json(path: Path, value: object) -> None:
    """Publish a complete JSON file with an exclusive hard-link operation."""
    path = Path(path).absolute()
    if any(p.is_symlink() for p in (path, *path.parents)):
        raise InputError("Output paths must not traverse symlinks")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=".crisprmap-", dir=path.parent)
    os.close(descriptor)
    try:
        write_json(Path(temporary), value)
        with open(temporary, "rb") as handle:
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as exc:
            raise InputError("Output already exists; refusing replacement") from exc
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        Path(temporary).unlink(missing_ok=True)
