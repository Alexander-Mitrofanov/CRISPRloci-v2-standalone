"""Bounded user inputs and atomic output transactions."""

from __future__ import annotations

import ctypes
import errno
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path


class InputError(ValueError):
    """Invalid bounded input or reference."""


ID = re.compile(r"[A-Za-z0-9_.:-]{1,160}\Z")
DNA = re.compile(r"[ACGTRYSWKMBDHVN]+\Z")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def reject_symlinks(path: Path) -> None:
    if any(part.is_symlink() for part in (path, *path.parents)):
        raise InputError("Symlink paths are not accepted")


def json_value(data: bytes) -> dict:
    def pairs(items):
        value = {}
        for key, item in items:
            if key in value:
                raise InputError("Duplicate JSON key")
            value[key] = item
        return value

    def invalid(value):
        raise InputError("Nonfinite JSON constant")

    value = json.loads(data, object_pairs_hook=pairs, parse_constant=invalid)
    if not isinstance(value, dict):
        raise InputError("JSON object required")
    return value


def file_digest(path: Path, limit: int = 512 * 1024 * 1024) -> str:
    reject_symlinks(path)
    if not path.is_file() or path.stat().st_size > limit:
        raise InputError("Invalid bounded reference file")
    result = hashlib.sha256()
    count = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            count += len(chunk)
            if count > limit:
                raise InputError("Reference grew beyond limit")
            result.update(chunk)
    return result.hexdigest()


def read(path: Path, limit: int = 2 * 1024 * 1024) -> bytes:
    reject_symlinks(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size > limit:
        raise InputError("Input must be a bounded regular file")
    with path.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise InputError("Input exceeded byte limit")
    return data


def fasta(path: Path, maximum: int = 1000, length: int = 200, limit: int = 2 * 1024 * 1024) -> list[dict]:
    return fasta_bytes(read(path, limit), maximum, length)


def fasta_bytes(data: bytes, maximum: int = 1000, length: int = 200) -> list[dict]:
    """Parse the exact bounded bytes used for the input identity receipt."""
    try:
        lines = data.decode("ascii").splitlines()
    except UnicodeError as error:
        raise InputError("FASTA must be ASCII") from error
    records: list[dict] = []
    seen: set[str] = set()
    chunks: list[str] = []
    bases = 0
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if records:
                records[-1]["sequence"] = "".join(chunks)
            chunks, bases = [], 0
            ident = line[1:].split()[0] if line[1:].strip() else ""
            if not ID.fullmatch(ident) or ident in seen or len(records) >= maximum:
                raise InputError("FASTA requires unique bounded identifiers and record count")
            seen.add(ident)
            records.append({"id": ident, "sequence": ""})
        elif records:
            bases += len(line)
            if bases > length:
                raise InputError("Sequence exceeded length limit")
            chunks.append(line.upper())
        else:
            raise InputError("FASTA sequence precedes header")
    if records:
        records[-1]["sequence"] = "".join(chunks)
    if not records or any(not DNA.fullmatch(row["sequence"]) for row in records):
        raise InputError("FASTA requires nonempty IUPAC DNA records")
    return records


def groups(path: Path, records: list[dict]) -> list[dict]:
    return groups_value(json_value(read(path)), records)


def groups_value(value: dict, records: list[dict]) -> list[dict]:
    """Validate the exact group-map value used for the input identity receipt."""
    if (
        not isinstance(value, dict)
        or set(value) != {"schema_version", "groups"}
        or value["schema_version"] != "spacepharer.groups/1.0.0"
    ):
        raise InputError("Unsupported group-map contract")
    rows = value["groups"]
    if not isinstance(rows, list) or not 1 <= len(rows) <= 100:
        raise InputError("Provide 1–100 explicit groups")
    members: list[str] = []
    identifiers: set[str] = set()
    for row in rows:
        if (
            not isinstance(row, dict)
            or set(row) != {"id", "members"}
            or not isinstance(row["id"], str)
            or not ID.fullmatch(row["id"])
            or row["id"] in identifiers
        ):
            raise InputError("Invalid group identifier")
        identifiers.add(row["id"])
        if (
            not isinstance(row["members"], list)
            or not row["members"]
            or any(not isinstance(x, str) for x in row["members"])
        ):
            raise InputError("Every group must name its occurrence members")
        members.extend(row["members"])
    if len(members) != len(set(members)) or set(members) != {row["id"] for row in records}:
        raise InputError("Groups must partition every input occurrence exactly once")
    return rows


def _rename_new_directory(source: Path, destination: Path) -> None:
    """Atomically publish without replacing a concurrently created directory."""
    if os.name == "nt":
        # Windows rename already refuses to replace an existing destination.
        source.rename(destination)
        return
    library = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith("linux") and hasattr(library, "renameat2"):
        rename = library.renameat2
        rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        # AT_FDCWD, RENAME_NOREPLACE: the existence check occurs inside rename.
        result = rename(-100, os.fsencode(source), -100, os.fsencode(destination), 1)
    elif sys.platform == "darwin" and hasattr(library, "renamex_np"):
        rename = library.renamex_np
        rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        result = rename(os.fsencode(source), os.fsencode(destination), 4)  # RENAME_EXCL
    else:
        raise InputError("This platform does not support atomic directory publication without replacement")
    if result:
        error = ctypes.get_errno()
        if error in (errno.EEXIST, errno.ENOTEMPTY):
            raise InputError("Output must be a new non-symlink directory")
        raise OSError(error, os.strerror(error), str(destination))



@contextmanager
def transaction(output: Path) -> Iterator[Path]:
    reject_symlinks(output)
    if output.exists() or output.is_symlink() or not output.parent.is_dir():
        raise InputError("Output must be a new directory under an existing parent")
    temporary = Path(tempfile.mkdtemp(prefix=".spacepharer-", dir=output.parent))
    try:
        yield temporary
        _rename_new_directory(temporary, output)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def write_json(path: Path, value: dict) -> None:
    data = json.dumps(value, ensure_ascii=True, allow_nan=False, sort_keys=True, indent=2).encode() + b"\n"
    if len(data) > 32 * 1024 * 1024:
        raise InputError("Report exceeds publication limit")
    path.write_bytes(data)
