"""Bounded, strict input and atomic output primitives; no network operations."""

from __future__ import annotations

from contextlib import contextmanager
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import stat
import tempfile
import zlib

DNA = frozenset("ACGTRYSWKMBDHVN")
MAX_BYTES = 256 * 1024 * 1024
MAX_BASES = 50_000_000
MAX_RECORDS = 100_000


class InputError(ValueError):
    """An input or evidence contract was violated."""


class LimitError(InputError):
    """A declared resource budget was exceeded; no complete result exists."""


def digest(value: str | bytes) -> str:
    return hashlib.sha256(value.encode("ascii") if isinstance(value, str) else value).hexdigest()


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, ensure_ascii=True, allow_nan=False, indent=2) + "\n").encode()


def identifier(value: object, label: str = "identifier") -> str:
    if not isinstance(value, str) or not value or len(value) > 500 or any(ord(c) < 33 or ord(c) > 126 for c in value):
        raise InputError(f"{label} must be 1–500 printable ASCII characters without whitespace")
    return value


def integer(value: object, label: str, lower: int, upper: int) -> int:
    if type(value) is not int or not lower <= value <= upper:
        raise InputError(f"{label} must be an integer in [{lower}, {upper}]")
    return value


def regular(path: Path) -> Path:
    path = Path(path).absolute()
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise InputError("Symlinks are not accepted as input paths")
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
        raise LimitError("Input exceeds its byte budget")
    with path.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise LimitError("Input exceeds its byte budget")
    return data


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise InputError(f"Duplicate JSON key: {key}")
        result[key] = value
    return result


def decode_json(data: bytes):
    def nonfinite(_):
        raise InputError("Non-finite JSON numbers are not supported")
    def finite_float(token):
        if len(token) > 128:
            raise InputError("Oversized JSON numeric token")
        number = float(token)
        if not math.isfinite(number):
            raise InputError("Non-finite JSON numbers are not supported")
        return number
    def bounded_int(token):
        if len(token) > 128:
            raise InputError("Oversized JSON numeric token")
        return int(token)
    try:
        return json.loads(data, object_pairs_hook=_pairs, parse_constant=nonfinite,
                          parse_float=finite_float, parse_int=bounded_int)
    except InputError:
        raise
    except (UnicodeError, ValueError, RecursionError) as exc:
        raise InputError("Malformed JSON input") from exc


def read_json(path: Path):
    return decode_json(read_bytes(path))


def read_fasta_snapshot(path: Path) -> tuple[list[dict], str]:
    """Read bounded plain/gzip DNA FASTA, preserving unique first-token IDs."""
    path = Path(path)
    raw_snapshot = read_bytes(path)
    raw_digest = digest(raw_snapshot)
    stream = gzip.GzipFile(fileobj=io.BytesIO(raw_snapshot)) if path.suffix == ".gz" else io.BytesIO(raw_snapshot)
    records, seen, chunks = [], set(), []
    current = None
    total = byte_count = 0

    def finish():
        if current is None:
            return
        sequence = "".join(chunks)
        if not sequence:
            raise InputError("Empty FASTA records are not supported")
        records.append({"id": current, "sequence": sequence, "sha256": digest(sequence),
                        "length": len(sequence), "assembly_id": None, "topology": "linear"})
        if len(records) > MAX_RECORDS:
            raise LimitError("Too many FASTA records")

    try:
        with stream as handle:
            while raw := handle.readline(1_000_001):
                byte_count += len(raw)
                if len(raw) > 1_000_000 or byte_count > MAX_BYTES:
                    raise LimitError("FASTA exceeds line or decompressed byte budget")
                line = raw.decode("ascii").strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish()
                    header = line[1:].split()
                    current = identifier(header[0] if header else None, "FASTA ID")
                    if current in seen:
                        raise InputError("Duplicate FASTA ID")
                    seen.add(current)
                    chunks = []
                else:
                    if current is None:
                        raise InputError("FASTA sequence precedes its header")
                    sequence = line.upper()
                    if not set(sequence) <= DNA:
                        raise InputError("Only DNA IUPAC letters are supported; gaps/RNA are not accepted")
                    total += len(sequence)
                    if total > MAX_BASES:
                        raise LimitError("FASTA exceeds total base budget")
                    chunks.append(sequence)
        finish()
    except (UnicodeError, EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise InputError("Malformed or corrupt FASTA") from exc
    if not records:
        raise InputError("FASTA contains no records")
    return records, raw_digest


def read_fasta(path: Path) -> list[dict]:
    return read_fasta_snapshot(path)[0]


def metadata(records: list[dict], path: Path | None) -> list[dict]:
    """Membership/topology are explicit declarations, never inferred from headers."""
    if path is None:
        return records
    data = read_json(path)
    if not isinstance(data, dict) or set(data) != {"schema_version", "records"} or data["schema_version"] != "1.0.0":
        raise InputError("Expected metadata schema_version 1.0.0 and records mapping")
    rows = data["records"]
    if not isinstance(rows, dict) or not set(rows) <= {r["id"] for r in records}:
        raise InputError("Metadata contains unknown record IDs")
    for record in records:
        row = rows.get(record["id"], {})
        if not isinstance(row, dict) or not set(row) <= {"assembly_id", "topology"}:
            raise InputError("Unknown record metadata fields")
        if row.get("assembly_id") is not None:
            record["assembly_id"] = identifier(row["assembly_id"], "assembly_id")
        if row.get("topology", "linear") not in {"linear", "circular"}:
            raise InputError("Topology must be linear or circular")
        record["topology"] = row.get("topology", "linear")
    return records


@contextmanager
def atomic_directory(destination: Path):
    """Never replace prior results; exclusive reservation protects concurrent writers."""
    destination = Path(destination).absolute()
    if any(p.is_symlink() for p in (destination, *destination.parents)):
        raise InputError("Output paths must not traverse symlinks")
    if destination.exists():
        raise InputError("Output directory already exists; choose a new destination")
    destination.parent.mkdir(parents=True, exist_ok=True)
    lock = destination.with_name(destination.name + ".partial")
    try:
        lock.mkdir(mode=0o700)
    except FileExistsError as exc:
        raise InputError("Another writer or interrupted reservation owns this output") from exc
    work = Path(tempfile.mkdtemp(prefix=".crisprspacer-", dir=destination.parent))
    try:
        yield work
        if destination.exists():
            raise InputError("Output appeared during execution; refusing replacement")
        for path in work.iterdir():
            if path.is_file():
                with path.open("rb") as handle:
                    os.fsync(handle.fileno())
        os.rename(work, destination)
        descriptor = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
    finally:
        if work.exists():
            shutil.rmtree(work)
        lock.rmdir()


def write_json(path: Path, value: object):
    encoder = json.JSONEncoder(sort_keys=True, ensure_ascii=True, allow_nan=False, indent=2)
    total = 0
    with path.open("wb") as handle:
        for chunk in encoder.iterencode(value):
            data = chunk.encode("ascii")
            total += len(data)
            if total + 1 > MAX_BYTES:
                raise LimitError("Serialized output exceeds byte budget")
            handle.write(data)
        handle.write(b"\n")


def write_manifest(directory: Path, kind: str):
    files = {p.name: {"sha256": digest(read_bytes(p)), "bytes": p.stat().st_size}
             for p in sorted(directory.iterdir()) if p.is_file()}
    write_json(directory / "manifest.json", {"schema_version": "1.0.0", "kind": kind, "files": files})


def verify_manifest(directory: Path, kind: str, expected: set[str]):
    manifest = read_json(directory / "manifest.json")
    if not isinstance(manifest, dict) or manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != kind:
        raise InputError("Unknown manifest schema/kind")
    files = manifest.get("files")
    if not isinstance(files, dict) or set(files) != expected:
        raise InputError("Manifest does not contain the exact required artifacts")
    if {p.name for p in directory.iterdir()} != expected | {"manifest.json"}:
        raise InputError("Unexpected files in bundle")
    payloads = {}
    total = 0
    for name in sorted(expected):
        identity = files[name]
        if not isinstance(identity, dict) or set(identity) != {"sha256", "bytes"}:
            raise InputError(f"Invalid artifact identity: {name}")
        integer(identity["bytes"], f"Artifact byte count: {name}", 0, MAX_BYTES)
        checksum = identity["sha256"]
        if not isinstance(checksum, str) or len(checksum) != 64 or set(checksum) - set("0123456789abcdef"):
            raise InputError(f"Invalid artifact digest: {name}")
        data = read_bytes(directory / name)
        total += len(data)
        if total > 512 * 1024 * 1024:
            raise LimitError("Bundle exceeds aggregate read budget")
        if identity != {"sha256": digest(data), "bytes": len(data)}:
            raise InputError(f"Artifact integrity check failed: {name}")
        payloads[name] = data
    return manifest, payloads
