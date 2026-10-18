"""Bounded discovery and privacy screening of bundle candidates."""

from __future__ import annotations

import codecs
from collections.abc import Callable
from pathlib import Path, PurePosixPath

from .contracts import _BUNDLE_SUFFIXES, _DENIED_PARTS, _INTERNAL_JSON
from .safety import (
    _checkpoint,
    _contains_private_path,
    _public_json_safe,
    _safe_file,
)


def _public_text_safe(
    path: Path,
    *,
    markers: tuple[str, ...],
    max_bytes: int,
    checkpoint: Callable[[], None] | None,
) -> bool:
    decoder = codecs.getincrementaldecoder("utf-8")(errors="strict")
    previous = ""
    observed_bytes = 0
    try:
        with path.open("rb") as stream:
            while block := stream.read(1024 * 1024):
                _checkpoint(checkpoint)
                observed_bytes += len(block)
                if observed_bytes > max_bytes:
                    return False
                combined = (previous + decoder.decode(block)).replace("\\", "/")
                if _contains_private_path(combined, markers):
                    return False
                previous = combined[-8192:]
            final = (previous + decoder.decode(b"", final=True)).replace("\\", "/")
    except (OSError, UnicodeError):
        return False
    _checkpoint(checkpoint)
    return not _contains_private_path(final, markers)


def _bundle_outputs(
    output_root: Path,
    *,
    max_file_bytes: int,
    max_files: int,
    max_total_bytes: int,
    markers: tuple[str, ...],
    checkpoint: Callable[[], None] | None,
) -> tuple[list[tuple[str, Path]], list[str]]:
    candidates: list[tuple[PurePosixPath, Path, int]] = []
    warnings: list[str] = []
    scanned = 0
    for path in output_root.rglob("*"):
        _checkpoint(checkpoint)
        scanned += 1
        if scanned > max_files:
            warnings.append(
                "Bundle filesystem-entry limit reached; remaining outputs were omitted."
            )
            break
        safe = _safe_file(output_root, path)
        if safe is None:
            continue
        relative = PurePosixPath(path.relative_to(output_root).as_posix())
        if (
            relative.suffix.lower() not in _BUNDLE_SUFFIXES
            or _DENIED_PARTS.intersection(relative.parts)
            or any(part.startswith(".") for part in relative.parts)
            or relative.as_posix() in _INTERNAL_JSON
        ):
            continue
        size = safe.stat().st_size
        if size == 0:
            warnings.append(f"Bundle file omitted because it is empty: {relative}")
            continue
        if size > max_file_bytes:
            warnings.append(
                f"Bundle file omitted because it exceeds the size limit: {relative}"
            )
            continue
        candidates.append((relative, safe, size))

    files: list[tuple[str, Path]] = []
    total_bytes = 0
    total_warning_added = False
    for relative, safe, size in sorted(candidates, key=lambda item: item[0].as_posix()):
        _checkpoint(checkpoint)
        if total_bytes + size > max_total_bytes:
            if not total_warning_added:
                warnings.append(
                    "Bundle byte limit reached; additional outputs were omitted."
                )
                total_warning_added = True
            continue
        suffix = relative.suffix.lower()
        if suffix == ".json":
            is_public_safe = _public_json_safe(
                safe,
                markers=markers,
                max_bytes=max_file_bytes,
                checkpoint=checkpoint,
            )
        else:
            try:
                is_public_safe = _public_text_safe(
                    safe,
                    markers=markers,
                    max_bytes=max_file_bytes,
                    checkpoint=checkpoint,
                )
            except OSError:
                is_public_safe = False
        if not is_public_safe:
            warnings.append(
                "Public artifact omitted because it was invalid or contained a private "
                f"path: {relative}"
            )
            continue
        files.append((f"results/{relative.as_posix()}", safe))
        total_bytes += size
    return files, warnings
