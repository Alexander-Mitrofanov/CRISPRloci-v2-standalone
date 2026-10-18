"""Private-path rejection, bounded JSON sanitization, and safe file identity."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
import hashlib
import json
import math
import mimetypes
import os
from pathlib import Path, PurePosixPath
import re
import stat

from .contracts import (
    _DENIED_PARTS,
    _EMBEDDED_ABSOLUTE_PATH,
    _IDENTIFY_FILENAMES,
    _INDIVIDUAL_EXACT,
    _MAX_PUBLIC_JSON_BYTES,
    _MEDIA_TYPES,
    _OMIT,
    _PRIVATE_MANIFEST_KEYS,
    _STATIC_PRIVATE_ROOTS,
)


def _checkpoint(callback: Callable[[], None] | None) -> None:
    if callback is not None:
        callback()


def _private_markers(private_roots: Iterable[str | Path]) -> tuple[str, ...]:
    markers = set(_STATIC_PRIVATE_ROOTS)
    for root in private_roots:
        rendered = Path(root).expanduser().resolve(strict=False).as_posix().rstrip("/")
        if rendered and rendered != "/":
            markers.add(rendered)
    return tuple(sorted(markers, key=len, reverse=True))


def _contains_private_path(value: object, markers: tuple[str, ...]) -> bool:
    if isinstance(value, str):
        normalized = value.replace("\\", "/")
        stripped = normalized.strip()
        if (
            stripped.startswith("/")
            or re.match(r"^[A-Za-z]:/", stripped)
            or _EMBEDDED_ABSOLUTE_PATH.search(normalized)
        ):
            return True
        return any(marker in normalized for marker in markers)
    if isinstance(value, Mapping):
        return any(
            _contains_private_path(key, markers)
            or _contains_private_path(item, markers)
            for key, item in value.items()
        )
    if isinstance(value, (list, tuple)):
        return any(_contains_private_path(item, markers) for item in value)
    return False


def _sanitize_manifest_value(
    value: object,
    markers: tuple[str, ...],
    redactions: list[int],
) -> object:
    if isinstance(value, Mapping):
        sanitized: dict[str, object] = {}
        for raw_key, raw_item in value.items():
            key = str(raw_key)
            normalized_key = key.lower()
            if normalized_key in _PRIVATE_MANIFEST_KEYS or normalized_key.endswith(
                ("_path", "_command")
            ):
                redactions[0] += 1
                continue
            item = _sanitize_manifest_value(raw_item, markers, redactions)
            if item is not _OMIT:
                sanitized[key] = item
        return sanitized
    if isinstance(value, (list, tuple)):
        sanitized_items = []
        for raw_item in value:
            item = _sanitize_manifest_value(raw_item, markers, redactions)
            if item is not _OMIT:
                sanitized_items.append(item)
        return sanitized_items
    if isinstance(value, str):
        if _contains_private_path(value, markers):
            redactions[0] += 1
            return _OMIT
        return value
    if value is None or isinstance(value, (bool, int, float)):
        if isinstance(value, float) and not math.isfinite(value):
            redactions[0] += 1
            return _OMIT
        return value
    redactions[0] += 1
    return _OMIT


def _read_json_bounded(
    path: Path,
    *,
    max_bytes: int,
    checkpoint: Callable[[], None] | None,
) -> object:
    if path.stat().st_size > max_bytes:
        raise ValueError("JSON exceeds the public artifact size limit")
    payload = bytearray()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            _checkpoint(checkpoint)
            if len(payload) + len(block) > max_bytes:
                raise ValueError("JSON exceeds the public artifact size limit")
            payload.extend(block)
    _checkpoint(checkpoint)
    return json.loads(payload)


def _write_sanitized_manifest(
    source: Path,
    target: Path,
    *,
    markers: tuple[str, ...],
    max_bytes: int,
    checkpoint: Callable[[], None] | None,
) -> tuple[Path | None, int]:
    try:
        raw = _read_json_bounded(
            source,
            max_bytes=min(max_bytes, _MAX_PUBLIC_JSON_BYTES),
            checkpoint=checkpoint,
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        return None, 0
    if not isinstance(raw, Mapping):
        return None, 0
    redactions = [0]
    sanitized = _sanitize_manifest_value(raw, markers, redactions)
    if not isinstance(sanitized, dict) or _contains_private_path(sanitized, markers):
        return None, redactions[0]
    sanitized["sanitization"] = {
        "private_fields_removed": redactions[0],
        "raw_internal_manifest_published": False,
    }
    with target.open("x", encoding="utf-8") as stream:
        json.dump(sanitized, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.chmod(target, 0o600)
    return target, redactions[0]


def _public_json_safe(
    path: Path,
    *,
    markers: tuple[str, ...],
    max_bytes: int,
    checkpoint: Callable[[], None] | None,
) -> bool:
    try:
        document = _read_json_bounded(
            path,
            max_bytes=min(max_bytes, _MAX_PUBLIC_JSON_BYTES),
            checkpoint=checkpoint,
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        return False
    return not _contains_private_path(document, markers)


def _allowed_individual(relative: PurePosixPath) -> bool:
    rendered = relative.as_posix()
    if rendered in _INDIVIDUAL_EXACT:
        return True
    parts = relative.parts
    if parts and parts[0] == "crispridentify" and relative.name in _IDENTIFY_FILENAMES:
        return True
    return bool(
        len(parts) >= 3
        and parts[0] in {"spacerplacer", "crispr_evor"}
        and parts[1] == "additional_data"
        and relative.suffix.lower() in {".json", ".nwk", ".tree", ".fa", ".fasta"}
        and not _DENIED_PARTS.intersection(parts)
    )


def _safe_file(root: Path, candidate: Path) -> Path | None:
    root = root.resolve()
    try:
        relative = candidate.relative_to(root)
    except ValueError:
        return None
    current = root
    for component in relative.parts:
        current = current / component
        try:
            metadata = current.lstat()
        except OSError:
            return None
        if stat.S_ISLNK(metadata.st_mode):
            return None
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError):
        return None
    try:
        metadata = resolved.lstat()
    except OSError:
        return None
    return resolved if stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1 else None


def resolve_registered_file(job_root: Path, relative_path: str) -> Path:
    relative = PurePosixPath(relative_path)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise FileNotFoundError(relative_path)
    candidate = job_root.joinpath(*relative.parts)
    safe = _safe_file(job_root, candidate)
    if safe is None:
        raise FileNotFoundError(relative_path)
    return safe


def _fsync_published(job_root: Path, paths: Iterable[Path]) -> None:
    directories: set[Path] = {job_root}
    for path in paths:
        descriptor = os.open(path, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        current = path.parent
        while True:
            current.relative_to(job_root)
            directories.add(current)
            if current == job_root:
                break
            current = current.parent
    for directory in sorted(directories, key=lambda value: len(value.parts), reverse=True):
        descriptor = os.open(
            directory,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0),
        )
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def file_sha256(path: Path, checkpoint: Callable[[], None] | None = None) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            _checkpoint(checkpoint)
            digest.update(block)
    return digest.hexdigest()


def _record(
    job_root: Path,
    path: Path,
    display_name: str,
    *,
    checkpoint: Callable[[], None] | None,
) -> dict[str, object]:
    relative = path.relative_to(job_root).as_posix()
    suffix = path.suffix.lower()
    size = path.stat().st_size
    if size <= 0:
        raise ValueError(f"refusing to register empty artifact: {display_name}")
    return {
        "artifact_id": hashlib.sha256(relative.encode("utf-8")).hexdigest()[:24],
        "relative_path": relative,
        "name": display_name,
        "size_bytes": size,
        "sha256": file_sha256(path, checkpoint),
        "media_type": _MEDIA_TYPES.get(
            suffix,
            mimetypes.guess_type(path.name)[0] or "application/octet-stream",
        ),
    }
