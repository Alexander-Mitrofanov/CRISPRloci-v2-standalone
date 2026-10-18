"""Exactly capped streaming ZIP creation."""

from __future__ import annotations

from collections.abc import Callable
import os
from pathlib import Path
import zipfile

from .safety import _checkpoint


class _ArchiveLimitExceeded(RuntimeError):
    pass


class _CappedArchiveWriter:
    """File proxy that prevents ZIP output from crossing its exact byte cap."""

    def __init__(self, stream, max_bytes: int):
        self._stream = stream
        self._max_bytes = max_bytes
        self._high_water = 0

    def write(self, value):
        end = self._stream.tell() + len(value)
        if max(self._high_water, end) > self._max_bytes:
            raise _ArchiveLimitExceeded
        written = self._stream.write(value)
        self._high_water = max(self._high_water, self._stream.tell())
        return written

    def __getattr__(self, name):
        return getattr(self._stream, name)


def _write_zip(
    temporary: Path,
    bundle_files: list[tuple[str, Path]],
    *,
    max_bytes: int,
    checkpoint: Callable[[], None] | None,
) -> None:
    with temporary.open("xb") as raw_stream:
        capped_stream = _CappedArchiveWriter(raw_stream, max_bytes)
        with zipfile.ZipFile(
            capped_stream,
            "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=6,
        ) as bundle:
            for archive_name, path in bundle_files:
                _checkpoint(checkpoint)
                with (
                    path.open("rb") as source,
                    bundle.open(archive_name, "w") as target,
                ):
                    while block := source.read(1024 * 1024):
                        _checkpoint(checkpoint)
                        target.write(block)
        raw_stream.flush()
        os.fsync(raw_stream.fileno())
