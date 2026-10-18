"""Bounded regular-file I/O, shared by all standalone commands."""

from __future__ import annotations

import builtins
import ctypes
import errno
import json
import os
from pathlib import Path
import stat
from typing import Any

from .errors import CCTKError

MAX_BYTES = 32 * 1024 * 1024
MAX_ARRAYS = 200
MAX_OCCURRENCES = 10000
MAX_SPACERS_PER_ARRAY = 500


def publish_directory(staging: Path, target: Path) -> None:
    """Atomically publish a complete directory without replacing a competing output."""
    # Python's os.rename may overwrite an empty directory created after a prior
    # existence check. CCTK requires Linux; renameat2 enforces this in the kernel.
    library = ctypes.CDLL(None, use_errno=True)
    rename = getattr(library, "renameat2", None)
    if rename is None:
        raise CCTKError("atomic non-replacing output publication is unavailable")
    rename.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
    rename.restype = ctypes.c_int
    if rename(-100, os.fsencode(staging), -100, os.fsencode(target), 1) != 0:
        code = ctypes.get_errno()
        if code == errno.EEXIST:
            raise CCTKError("output directory appeared during analysis")
        raise OSError(code, os.strerror(code), str(target))


class BoundedFile:
    def __init__(self, stream: Any, limit: int = MAX_BYTES):
        self.stream = stream
        self.remaining = limit

    def __enter__(self) -> BoundedFile:
        return self

    def __exit__(self, *args: Any) -> None:
        self.stream.close()

    def __getattr__(self, name: str) -> Any:
        return getattr(self.stream, name)

    def _consume(self, value: Any) -> Any:
        size = len(value.encode("utf-8")) if isinstance(value, str) else len(value)
        self.remaining -= size
        if self.remaining < 0:
            raise CCTKError("file exceeds bounded I/O limit")
        return value

    def read(self, size: int = -1) -> Any:
        return self._consume(
            self.stream.read(min(size, self.remaining + 1) if size >= 0 else self.remaining + 1)
        )

    def readline(self, size: int = -1) -> Any:
        return self._consume(
            self.stream.readline(min(size, self.remaining + 1) if size >= 0 else self.remaining + 1)
        )

    def __iter__(self) -> BoundedFile:
        return self

    def __next__(self) -> Any:
        value = self.readline()
        if not value:
            raise StopIteration
        return value

    def readlines(self) -> list[Any]:
        return list(self)

    def write(self, value: Any) -> int:
        return self.stream.write(self._consume(value))

    def writelines(self, values: Any) -> None:
        for value in values:
            self.write(value)


def safe_open(
    path: str | Path, mode: str = "r", encoding: str | None = None, **kwargs: Any
) -> BoundedFile:
    """Refuse symlinks, devices and unbounded reads/writes, including growing files."""
    target = Path(path)
    if any(p.is_symlink() for p in (target, *target.parents)):
        raise CCTKError("symlinked paths are not supported")
    if set(mode) - set("rwaxb+t"):
        raise CCTKError("unsupported file mode")

    def opener(name: str, flags: int) -> int:
        descriptor = os.open(name, flags | os.O_NOFOLLOW | os.O_NONBLOCK, 0o600)
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode) or info.st_size > MAX_BYTES:
            os.close(descriptor)
            raise CCTKError("input must be a bounded regular file")
        return descriptor

    stream = builtins.open(
        target,
        mode,
        encoding=None if "b" in mode else (encoding or "utf-8"),
        opener=opener,
        **kwargs,
    )
    return BoundedFile(
        stream, MAX_BYTES - (os.fstat(stream.fileno()).st_size if "a" in mode else 0)
    )


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result = {}
    for key, value in pairs:
        if key in result:
            raise CCTKError("duplicate JSON field")
        result[key] = value
    return result


def load_json(path: str | Path) -> Any:
    with safe_open(path) as handle:
        try:
            return json.loads(
                handle.read(),
                object_pairs_hook=reject_duplicate_keys,
                parse_constant=lambda _: (_ for _ in ()).throw(CCTKError("non-finite JSON number")),
            )
        except (json.JSONDecodeError, RecursionError) as error:
            raise CCTKError("invalid JSON input") from error


def write_json(path: str | Path, value: Any) -> None:
    with safe_open(path, "w") as handle:
        json.dump(value, handle, allow_nan=False, indent=2, ensure_ascii=True)
        handle.write("\n")


def read_arrays(path: str | Path) -> dict[str, list[str]]:
    arrays: dict[str, list[str]] = {}
    occurrences = 0
    with safe_open(path) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            ident, *spacers = fields
            if ident in arrays or len(ident) > 256 or not spacers:
                raise CCTKError("array identifiers must be unique and arrays nonempty")
            if any(len(spacer) > 1000 for spacer in spacers):
                raise CCTKError("spacer symbol too long")
            occurrences += len(spacers)
            if (
                len(arrays) >= MAX_ARRAYS
                or len(spacers) > MAX_SPACERS_PER_ARRAY
                or occurrences > MAX_OCCURRENCES
            ):
                raise CCTKError("array input exceeds supported limits")
            arrays[ident] = spacers
    if not arrays:
        raise CCTKError("no arrays supplied")
    return arrays


def read_fasta(path: str | Path) -> dict[str, str]:
    fragments: dict[str, list[str]] = {}
    ident: str | None = None
    with safe_open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                ident = line[1:]
                if not ident or len(ident) > 1024 or ident in fragments:
                    raise CCTKError("FASTA identifiers must be nonempty and unique")
                fragments[ident] = []
            elif ident is None or not set(line.upper()) <= set("ACGTRYSWKMBDHVN"):
                raise CCTKError("invalid DNA FASTA")
            else:
                fragments[ident].append(line)
    sequences = {ident: "".join(parts) for ident, parts in fragments.items()}
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise CCTKError("FASTA contains no sequence or an empty record")
    return sequences
