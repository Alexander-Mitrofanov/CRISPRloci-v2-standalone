"""Attempt storage accounting, quarantine, and fail-closed cleanup."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import stat

from ..fasta import fsync_directory


class StorageCleanupError(RuntimeError):
    """A storage-violating attempt could not be safely purged."""

    pass


def _accounted_storage_bytes(metadata: os.stat_result) -> int:
    """Count apparent sparse size and allocated blocks for every node type."""

    allocated = max(0, int(getattr(metadata, "st_blocks", 0))) * 512
    if stat.S_ISREG(metadata.st_mode):
        return max(int(metadata.st_size), allocated)
    return allocated


def _job_tree_within_budget(root: Path, max_bytes: int, max_entries: int) -> bool:
    total_bytes = 0
    entries = 0
    try:
        root_metadata = root.lstat()
        if not stat.S_ISDIR(root_metadata.st_mode):
            return False
        total_bytes = _accounted_storage_bytes(root_metadata)
        if total_bytes > max_bytes:
            return False
        for path in root.rglob("*"):
            entries += 1
            if entries > max_entries:
                return False
            try:
                metadata = path.lstat()
            except FileNotFoundError:
                # Scientific tools remove temporary files while this scan runs.
                # A vanished entry consumes no storage, but still counts toward
                # the bounded traversal. Other inspection errors fail closed.
                continue
            total_bytes += _accounted_storage_bytes(metadata)
            if total_bytes > max_bytes:
                return False
        final_metadata = root.lstat()
        if (
            not stat.S_ISDIR(final_metadata.st_mode)
            or (final_metadata.st_dev, final_metadata.st_ino)
            != (root_metadata.st_dev, root_metadata.st_ino)
        ):
            return False
    except OSError:
        return False
    return True


def _attempt_path_is_present(job_root: Path, candidate: Path) -> bool:
    """Validate fixed parents without following a tool-created symlink."""

    try:
        relative = candidate.relative_to(job_root)
    except ValueError as error:
        raise StorageCleanupError(
            "attempt cleanup path escaped its job root"
        ) from error
    if not relative.parts:
        raise StorageCleanupError("attempt cleanup cannot target the job root")

    current = job_root
    try:
        root_metadata = current.lstat()
    except OSError as error:
        raise StorageCleanupError(
            "job root could not be inspected for cleanup"
        ) from error
    if not stat.S_ISDIR(root_metadata.st_mode) or stat.S_ISLNK(root_metadata.st_mode):
        raise StorageCleanupError("job root is not a safe cleanup directory")

    for component in relative.parts[:-1]:
        current = current / component
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return False
        except OSError as error:
            raise StorageCleanupError(
                "attempt cleanup parent could not be inspected"
            ) from error
        if not stat.S_ISDIR(metadata.st_mode) or stat.S_ISLNK(metadata.st_mode):
            raise StorageCleanupError("attempt cleanup parent is not a safe directory")
    try:
        candidate.lstat()
    except FileNotFoundError:
        return False
    except OSError as error:
        raise StorageCleanupError(
            "attempt cleanup target could not be inspected"
        ) from error
    return True


def _make_private_tree_removable(root: Path) -> None:
    """Restore owner traversal on quarantined directories without following links."""

    for current, directories, _files in os.walk(root, followlinks=False):
        current_path = Path(current)
        current_path.chmod(0o700)
        for name in directories:
            child = current_path / name
            metadata = child.lstat()
            if stat.S_ISDIR(metadata.st_mode) and not stat.S_ISLNK(metadata.st_mode):
                child.chmod(0o700)


def _quarantine_and_purge_attempt(
    job_root: Path,
    attempt: int,
    candidates: tuple[tuple[str, Path], ...],
) -> None:
    """Atomically isolate and remove only one failed attempt's mutable paths."""

    present: list[tuple[str, Path]] = []
    for label, candidate in candidates:
        if not label or Path(label).name != label:
            raise StorageCleanupError("attempt cleanup label is invalid")
        if _attempt_path_is_present(job_root, candidate):
            present.append((label, candidate))
    if not present:
        return

    quarantine = job_root / f".storage-quarantine-attempt-{attempt:03d}"
    touched_parents: set[Path] = set()
    try:
        quarantine.mkdir(mode=0o700, exist_ok=False)
        for label, candidate in present:
            candidate.rename(quarantine / label)
            touched_parents.add(candidate.parent)
        for parent in touched_parents:
            fsync_directory(parent)
        fsync_directory(quarantine)
        _make_private_tree_removable(quarantine)
        shutil.rmtree(quarantine)
        fsync_directory(job_root)
    except StorageCleanupError:
        raise
    except (OSError, ValueError) as error:
        raise StorageCleanupError(
            "current attempt could not be safely quarantined and purged"
        ) from error
