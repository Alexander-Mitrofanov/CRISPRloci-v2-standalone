"""Bounded filesystem and JSON primitives for public summaries."""

from __future__ import annotations

import json
from pathlib import Path
import stat
from typing import Any, Callable, Mapping

from .contracts import _MAX_RECONSTRUCTION_INPUT_BYTES, _MAX_SUMMARY_INPUT_BYTES


def _checkpoint(callback: Callable[[], None] | None) -> None:
    if callback is not None:
        callback()


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
    return resolved if resolved.is_file() else None


def _safe_directory(root: Path, candidate: Path) -> Path | None:
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
    return resolved if resolved.is_dir() else None


def _read_bounded(
    root: Path,
    path: Path,
    *,
    max_bytes: int,
    checkpoint: Callable[[], None] | None,
) -> bytes:
    safe = _safe_file(root, path)
    if safe is None:
        raise FileNotFoundError(path.name)
    if safe.stat().st_size > max_bytes:
        raise ValueError(f"{path.name} exceeds the summary parser size limit")
    payload = bytearray()
    with safe.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            _checkpoint(checkpoint)
            if len(payload) + len(block) > max_bytes:
                raise ValueError(f"{path.name} exceeds the summary parser size limit")
            payload.extend(block)
    _checkpoint(checkpoint)
    return bytes(payload)


def _json(
    root: Path,
    path: Path,
    *,
    checkpoint: Callable[[], None] | None,
    max_bytes: int = _MAX_SUMMARY_INPUT_BYTES,
) -> Any:
    return json.loads(
        _read_bounded(
            root,
            path,
            max_bytes=max_bytes,
            checkpoint=checkpoint,
        )
    )


def _read_reconstruction_bytes(
    root: Path,
    path: Path,
    *,
    budget: list[int],
    checkpoint: Callable[[], None] | None,
) -> bytes:
    safe = _safe_file(root, path)
    if safe is None:
        raise FileNotFoundError(path.name)
    size = safe.stat().st_size
    if size > _MAX_RECONSTRUCTION_INPUT_BYTES:
        raise ValueError("reconstruction input exceeds its per-file size limit")
    if size > budget[0]:
        raise ValueError("reconstruction inputs exceed their total size limit")
    payload = _read_bounded(
        root,
        safe,
        max_bytes=_MAX_RECONSTRUCTION_INPUT_BYTES,
        checkpoint=checkpoint,
    )
    budget[0] -= len(payload)
    return payload


def _reconstruction_json(
    root: Path,
    path: Path,
    *,
    budget: list[int],
    checkpoint: Callable[[], None] | None,
) -> Any:
    def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("reconstruction JSON contains duplicate keys")
            result[key] = value
        return result

    return json.loads(
        _read_reconstruction_bytes(
            root,
            path,
            budget=budget,
            checkpoint=checkpoint,
        ).decode("utf-8"),
        object_pairs_hook=unique_object,
    )


def _mapping(value: object) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


