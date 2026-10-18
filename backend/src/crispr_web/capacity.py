"""Fail-closed byte and inode capacity probes for private service storage."""

from __future__ import annotations

import os
from pathlib import Path
import shutil


def available_capacity(root: Path) -> tuple[int, int]:
    """Return bytes and unprivileged inodes available on the root filesystem."""

    free_bytes = int(shutil.disk_usage(root).free)
    filesystem = os.statvfs(root)
    free_inodes = int(filesystem.f_favail)
    if free_bytes < 0 or free_inodes < 0:
        raise OSError("storage availability returned a negative value")
    return free_bytes, free_inodes


def has_capacity(
    root: Path,
    *,
    required_bytes: int,
    required_inodes: int,
) -> bool:
    """Fail closed when either capacity dimension cannot satisfy its reserve."""

    try:
        free_bytes, free_inodes = available_capacity(root)
    except OSError:
        return False
    return free_bytes >= required_bytes and free_inodes >= required_inodes
