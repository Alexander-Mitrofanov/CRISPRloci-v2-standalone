"""Safe bounded reads and log capture for worker attempts."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
import stat
import threading
from typing import Any, Callable, Mapping


LOGGER = logging.getLogger("crispr_web.worker")


def _safe_nonempty_output(root: Path, candidate: Path) -> bool:
    try:
        relative = candidate.relative_to(root)
    except ValueError:
        return False
    current = root
    for component in relative.parts:
        current = current / component
        try:
            metadata = current.lstat()
        except OSError:
            return False
        if stat.S_ISLNK(metadata.st_mode):
            return False
    try:
        resolved_root = root.resolve(strict=True)
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(resolved_root)
    except (OSError, ValueError):
        return False
    return resolved.is_file() and resolved.stat().st_size > 0



def _read_bounded_output(
    root: Path,
    path: Path,
    *,
    max_bytes: int = 20_000_000,
    checkpoint: Callable[[], None] | None = None,
) -> bytes | None:
    if not _safe_nonempty_output(root, path):
        return None
    try:
        if path.stat().st_size > max_bytes:
            return None
        payload = bytearray()
        with path.open("rb") as stream:
            while block := stream.read(1024 * 1024):
                if checkpoint is not None:
                    checkpoint()
                if len(payload) + len(block) > max_bytes:
                    return None
                payload.extend(block)
    except OSError:
        return None
    return bytes(payload)


def _read_json(
    root: Path,
    path: Path,
    *,
    max_bytes: int = 20_000_000,
    checkpoint: Callable[[], None] | None = None,
) -> Mapping[str, Any]:
    payload = _read_bounded_output(
        root,
        path,
        max_bytes=max_bytes,
        checkpoint=checkpoint,
    )
    if payload is None:
        return {}
    try:
        value = json.loads(payload)
    except (UnicodeError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, Mapping) else {}


def _capture_stream(
    stream,
    path: Path,
    max_bytes: int,
    exceeded: threading.Event,
) -> None:
    written = 0
    try:
        with path.open("xb") as target:
            os.chmod(path, 0o600)
            while block := stream.read(64 * 1024):
                remaining = max(0, max_bytes - written)
                if remaining:
                    portion = block[:remaining]
                    target.write(portion)
                    written += len(portion)
                if len(block) > remaining:
                    exceeded.set()
            target.flush()
            os.fsync(target.fileno())
    except OSError:
        LOGGER.exception("could not capture bounded worker log %s", path)
        exceeded.set()
    finally:
        try:
            stream.close()
        except OSError:
            pass


