#!/usr/bin/python3
"""Create or restore a durable SQLite release backup while services are stopped."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import stat
import time
from urllib.parse import quote


ROOT_UID = 0
ROOT_GID = 0
PRIVATE_DIRECTORY_MODE = 0o700
PRIVATE_FILE_MODE = 0o600
SQLITE_SIDECAR_SUFFIXES = ("-journal", "-wal", "-shm")
MAX_DATABASE_BACKUP_BYTES = 1_000_000_000
MIN_BACKUP_DESTINATION_FREE_BYTES = 2_000_000_000
MIN_BACKUP_DESTINATION_FREE_INODES = 1_024


def fail(message: str) -> None:
    raise SystemExit(f"database release backup failed: {message}")


def _absolute_path(path: Path, label: str) -> Path:
    path = Path(path)
    if (
        not path.is_absolute()
        or path.anchor != os.sep
        or path != Path(os.path.normpath(path))
    ):
        fail(f"{label} must be an absolute normalized path: {path}")
    if path == Path("/") or not path.name:
        fail(f"{label} must name a file or non-root directory: {path}")
    return path


def _open_directory_nofollow(path: Path, label: str) -> int:
    """Open every path component with O_NOFOLLOW and return the pinned directory."""
    path = _absolute_path(path, label)
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | os.O_NOFOLLOW
    descriptor = os.open("/", flags)
    try:
        for component in path.parts[1:]:
            next_descriptor = os.open(component, flags, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = next_descriptor
    except OSError as error:
        os.close(descriptor)
        fail(f"{label} is missing, unsafe, or contains a symlink: {error.errno}")
    return descriptor


def _require_private_root_directory(path: Path, label: str) -> int:
    descriptor = _open_directory_nofollow(path, label)
    metadata = os.fstat(descriptor)
    if (
        metadata.st_uid != ROOT_UID
        or metadata.st_gid != ROOT_GID
        or stat.S_IMODE(metadata.st_mode) != PRIVATE_DIRECTORY_MODE
    ):
        os.close(descriptor)
        fail(f"{label} must be exactly root:root mode 0700")
    return descriptor


def _open_parent(path: Path, label: str) -> tuple[Path, int, str]:
    path = _absolute_path(path, label)
    return path, _open_directory_nofollow(path.parent, f"{label} parent"), path.name


def _lstat_at(directory: int, name: str) -> os.stat_result | None:
    try:
        return os.stat(name, dir_fd=directory, follow_symlinks=False)
    except FileNotFoundError:
        return None


def _require_regular_at(
    directory: int,
    name: str,
    label: str,
    *,
    flags: int = os.O_RDONLY,
) -> tuple[int, os.stat_result]:
    try:
        descriptor = os.open(
            name,
            flags | os.O_CLOEXEC | os.O_NOFOLLOW,
            dir_fd=directory,
        )
    except OSError as error:
        fail(f"{label} is missing, unsafe, or not openable: {error.errno}")
    metadata = os.fstat(descriptor)
    if not stat.S_ISREG(metadata.st_mode):
        os.close(descriptor)
        fail(f"{label} is not a regular file")
    return descriptor, metadata


def _require_absent_at(directory: int, name: str, label: str) -> None:
    if _lstat_at(directory, name) is not None:
        fail(f"{label} already exists")


def _cleanup_regular_entries(directory: int, entries: tuple[str, ...]) -> None:
    """Remove known transient regular files, rejecting links and special files."""
    changed = False
    unsafe: list[str] = []
    try:
        for name in entries:
            metadata = _lstat_at(directory, name)
            if metadata is None:
                continue
            if not stat.S_ISREG(metadata.st_mode):
                unsafe.append(name)
                continue
            os.unlink(name, dir_fd=directory)
            changed = True
    finally:
        if changed:
            os.fsync(directory)
    if unsafe:
        fail("unsafe transient or SQLite sidecar path: " + ", ".join(unsafe))


def _reject_entries(directory: int, entries: tuple[str, ...], label: str) -> None:
    for name in entries:
        if _lstat_at(directory, name) is not None:
            fail(f"{label} must not exist: {name}")


def _require_safe_regular_or_absent_entries(
    directory: int, entries: tuple[str, ...], label: str
) -> None:
    unsafe: list[str] = []
    for name in entries:
        metadata = _lstat_at(directory, name)
        if metadata is not None and not stat.S_ISREG(metadata.st_mode):
            unsafe.append(name)
    if unsafe:
        fail(f"unsafe {label}: " + ", ".join(unsafe))


def _sidecars(name: str) -> tuple[str, ...]:
    return tuple(f"{name}{suffix}" for suffix in SQLITE_SIDECAR_SUFFIXES)


def _read_descriptor(descriptor: int, *, maximum: int | None = None) -> bytes:
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks: list[bytes] = []
    size = 0
    while block := os.read(descriptor, 1024 * 1024):
        size += len(block)
        if maximum is not None and size > maximum:
            fail("backup manifest exceeds the size limit")
        chunks.append(block)
    return b"".join(chunks)


def _sha256_descriptor(descriptor: int) -> str:
    digest = hashlib.sha256()
    os.lseek(descriptor, 0, os.SEEK_SET)
    while block := os.read(descriptor, 1024 * 1024):
        digest.update(block)
    return digest.hexdigest()


def _sha256_at(directory: int, name: str, label: str) -> str:
    descriptor, _ = _require_regular_at(directory, name, label)
    try:
        return _sha256_descriptor(descriptor)
    finally:
        os.close(descriptor)


def _sqlite_path_for_descriptor(descriptor: int, *, read_only: bool = False) -> str:
    mode = "ro" if read_only else "rw"
    immutable = "&immutable=1" if read_only else ""
    return f"file:/proc/self/fd/{descriptor}?mode={mode}{immutable}"


def _sqlite_path_in_directory(directory: int, name: str) -> str:
    return f"file:/proc/self/fd/{directory}/{quote(name, safe='')}?mode=rw"


def _transaction_name(
    operation: str, directory: Path, artifact: Path, binding: str
) -> str:
    digest = hashlib.sha256(
        b"\0".join(
            value.encode("utf-8")
            for value in (operation, str(directory), str(artifact), binding)
        )
    ).hexdigest()[:24]
    return f".{artifact.name}.{operation}-{digest}"


def _close_database(connection: sqlite3.Connection | None) -> None:
    if connection is None:
        return
    try:
        connection.close()
    except sqlite3.Error:
        # Cleanup of pinned transient paths must still run after a close failure.
        pass


def _validate_manifest(document: object) -> dict[str, object]:
    if not isinstance(document, dict):
        fail("backup manifest must be a JSON object")
    if document.get("schema") != {
        "name": "crispr-web-database-backup",
        "version": "1.0.0",
    }:
        fail("backup manifest schema is unsupported")
    release = document.get("release_before_switch")
    database_schema = document.get("database_schema_version")
    created_at = document.get("created_at_unix")
    size = document.get("size_bytes")
    digest = document.get("sha256")
    if not isinstance(release, str) or not release:
        fail("backup manifest release identity is invalid")
    if database_schema is not None and not isinstance(database_schema, str):
        fail("backup manifest database schema version is invalid")
    if (
        isinstance(created_at, bool)
        or not isinstance(created_at, (int, float))
        or not math.isfinite(created_at)
    ):
        fail("backup manifest creation time is invalid")
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        fail("backup manifest database size is invalid")
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        fail("backup manifest database digest is invalid")
    return document


def _checkpoint_existing_database(directory: int, name: str, label: str) -> None:
    descriptor, _ = _require_regular_at(directory, name, label, flags=os.O_RDWR)
    connection: sqlite3.Connection | None = None
    try:
        connection = sqlite3.connect(
            _sqlite_path_for_descriptor(descriptor), uri=True, timeout=30
        )
        connection.execute("PRAGMA busy_timeout = 30000")
        checkpoint = connection.execute("PRAGMA wal_checkpoint(TRUNCATE)").fetchone()
        if checkpoint is None or int(checkpoint[0]) != 0:
            fail("existing destination WAL checkpoint remained busy")
        integrity(connection)
        connection.close()
        connection = None
        os.fsync(descriptor)
        os.fsync(directory)
    finally:
        _close_database(connection)
        os.close(descriptor)


def regular_non_symlink(path: Path, *, must_exist: bool = True) -> Path:
    """Compatibility validator used by callers importing this control script."""
    path, parent, name = _open_parent(path, "path")
    try:
        metadata = _lstat_at(parent, name)
        if must_exist:
            if metadata is None or not stat.S_ISREG(metadata.st_mode):
                fail(f"path is not a regular file: {path}")
        elif metadata is not None and not stat.S_ISREG(metadata.st_mode):
            fail(f"path is not a regular file: {path}")
    finally:
        os.close(parent)
    return path


def fsync_directory(path: Path) -> None:
    descriptor = _open_directory_nofollow(path, "directory to fsync")
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def sha256(path: Path) -> str:
    path, parent, name = _open_parent(path, "file to hash")
    del path
    try:
        return _sha256_at(parent, name, "file to hash")
    finally:
        os.close(parent)


def integrity(connection: sqlite3.Connection) -> None:
    row = connection.execute("PRAGMA integrity_check").fetchone()
    if row is None or row[0] != "ok":
        fail("SQLite integrity_check did not return ok")


def schema_version(connection: sqlite3.Connection) -> str | None:
    try:
        row = connection.execute(
            "SELECT value FROM metadata WHERE name = 'schema_version'"
        ).fetchone()
    except sqlite3.Error:
        return None
    return None if row is None else str(row[0])


def require_backup_capacity(
    source: sqlite3.Connection,
    destination_directory: int,
) -> None:
    """Reject a backup that could consume the root filesystem reserve."""

    page_count = int(source.execute("PRAGMA page_count").fetchone()[0])
    page_size = int(source.execute("PRAGMA page_size").fetchone()[0])
    expected_bytes = page_count * page_size
    if expected_bytes < 0 or expected_bytes > MAX_DATABASE_BACKUP_BYTES:
        fail("source database exceeds the backup byte ceiling")
    capacity = os.fstatvfs(destination_directory)
    free_bytes = int(capacity.f_bavail) * int(capacity.f_frsize)
    free_inodes = int(capacity.f_favail)
    if free_bytes < expected_bytes + MIN_BACKUP_DESTINATION_FREE_BYTES:
        fail("backup destination is below its free-space reserve")
    if free_inodes < MIN_BACKUP_DESTINATION_FREE_INODES:
        fail("backup destination is below its free-inode reserve")


def create_backup(source: Path, destination: Path, release: str) -> None:
    destination = _absolute_path(destination, "backup destination")
    destination_parent = _require_private_root_directory(
        destination.parent, "backup destination parent"
    )
    try:
        source, source_parent, source_name = _open_parent(source, "backup source")
    except BaseException:
        os.close(destination_parent)
        raise
    destination_name = destination.name
    manifest_name = f"{destination_name}.json"
    temporary_name = _transaction_name(
        "backup", destination.parent, destination, release
    )
    temporary_manifest_name = f"{temporary_name}.json"
    temporary_entries = (
        temporary_name,
        temporary_manifest_name,
        *_sidecars(temporary_name),
    )
    destination_sidecars = _sidecars(destination_name)
    source_descriptor: int | None = None
    source_connection: sqlite3.Connection | None = None
    destination_connection: sqlite3.Connection | None = None
    published: list[str] = []
    manage_destination_sidecars = False
    succeeded = False
    try:
        _cleanup_regular_entries(destination_parent, temporary_entries)
        _require_absent_at(destination_parent, destination_name, "backup destination")
        _require_absent_at(destination_parent, manifest_name, "backup manifest")
        manage_destination_sidecars = True
        _cleanup_regular_entries(destination_parent, destination_sidecars)
        source_descriptor, source_metadata = _require_regular_at(
            source_parent,
            source_name,
            "backup source",
            flags=os.O_RDWR,
        )
        source_parent_metadata = os.fstat(source_parent)
        if (
            source_metadata.st_uid <= 0
            or source_metadata.st_uid != source_parent_metadata.st_uid
            or source_metadata.st_gid != source_parent_metadata.st_gid
            or stat.S_IMODE(source_metadata.st_mode) != PRIVATE_FILE_MODE
            or stat.S_IMODE(source_parent_metadata.st_mode) != PRIVATE_DIRECTORY_MODE
        ):
            fail(
                "backup source and parent must share a non-root owner/group "
                "with exact modes 0600/0700"
            )
        _require_safe_regular_or_absent_entries(
            source_parent, _sidecars(source_name), "backup source SQLite sidecar"
        )
        temporary_descriptor = os.open(
            temporary_name,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            PRIVATE_FILE_MODE,
            dir_fd=destination_parent,
        )
        os.close(temporary_descriptor)
        source_connection = sqlite3.connect(
            _sqlite_path_for_descriptor(source_descriptor),
            uri=True,
            timeout=30,
        )
        source_connection.execute("PRAGMA busy_timeout = 30000")
        checkpoint = source_connection.execute(
            "PRAGMA wal_checkpoint(TRUNCATE)"
        ).fetchone()
        if checkpoint is None or int(checkpoint[0]) != 0:
            fail("WAL checkpoint remained busy")
        integrity(source_connection)
        source_schema = schema_version(source_connection)
        require_backup_capacity(source_connection, destination_parent)
        destination_connection = sqlite3.connect(
            _sqlite_path_in_directory(destination_parent, temporary_name),
            uri=True,
            timeout=30,
        )
        source_connection.backup(destination_connection)
        integrity(destination_connection)
        if schema_version(destination_connection) != source_schema:
            fail("backup schema metadata differs from source")
        destination_connection.commit()
        destination_connection.close()
        destination_connection = None
        source_connection.close()
        source_connection = None
        os.fsync(source_parent)

        temporary_descriptor, temporary_metadata = _require_regular_at(
            destination_parent,
            temporary_name,
            "temporary backup",
            flags=os.O_RDWR,
        )
        try:
            os.fchown(temporary_descriptor, ROOT_UID, ROOT_GID)
            os.fchmod(temporary_descriptor, PRIVATE_FILE_MODE)
            os.fsync(temporary_descriptor)
            temporary_metadata = os.fstat(temporary_descriptor)
            backup_digest = _sha256_descriptor(temporary_descriptor)
        finally:
            os.close(temporary_descriptor)
        document = {
            "schema": {
                "name": "crispr-web-database-backup",
                "version": "1.0.0",
            },
            "release_before_switch": release,
            "database_schema_version": source_schema,
            "created_at_unix": time.time(),
            "size_bytes": temporary_metadata.st_size,
            "sha256": backup_digest,
        }
        document = _validate_manifest(document)
        manifest_descriptor = os.open(
            temporary_manifest_name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            PRIVATE_FILE_MODE,
            dir_fd=destination_parent,
        )
        try:
            payload = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode(
                "utf-8"
            )
            view = memoryview(payload)
            while view:
                written = os.write(manifest_descriptor, view)
                if written <= 0:
                    fail("backup manifest write did not make progress")
                view = view[written:]
            os.fchown(manifest_descriptor, ROOT_UID, ROOT_GID)
            os.fchmod(manifest_descriptor, PRIVATE_FILE_MODE)
            os.fsync(manifest_descriptor)
        finally:
            os.close(manifest_descriptor)

        os.rename(
            temporary_name,
            destination_name,
            src_dir_fd=destination_parent,
            dst_dir_fd=destination_parent,
        )
        published.append(destination_name)
        os.rename(
            temporary_manifest_name,
            manifest_name,
            src_dir_fd=destination_parent,
            dst_dir_fd=destination_parent,
        )
        published.append(manifest_name)
        _cleanup_regular_entries(
            destination_parent,
            (*_sidecars(temporary_name), *_sidecars(destination_name)),
        )
        os.fsync(destination_parent)
        succeeded = True
        print(json.dumps(document, sort_keys=True))
    finally:
        try:
            try:
                _close_database(destination_connection)
                _close_database(source_connection)
            finally:
                if source_descriptor is not None:
                    os.close(source_descriptor)
        finally:
            try:
                try:
                    if not succeeded and published:
                        _cleanup_regular_entries(destination_parent, tuple(published))
                finally:
                    try:
                        _cleanup_regular_entries(destination_parent, temporary_entries)
                    finally:
                        if manage_destination_sidecars:
                            _cleanup_regular_entries(
                                destination_parent, destination_sidecars
                            )
            finally:
                os.close(destination_parent)
                os.close(source_parent)


def _path_is_beneath(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _same_filesystem(first: os.stat_result, second: os.stat_result) -> bool:
    return first.st_dev == second.st_dev


def _restore_destination_metadata(
    directory: int,
    name: str,
    *,
    expected_uid: int,
    expected_gid: int,
    mode: int,
) -> os.stat_result | None:
    metadata = _lstat_at(directory, name)
    if metadata is not None:
        if not stat.S_ISREG(metadata.st_mode):
            fail("restore destination exists but is not a regular file")
        if (
            metadata.st_uid != expected_uid
            or metadata.st_gid != expected_gid
            or stat.S_IMODE(metadata.st_mode) != mode
        ):
            fail("restore destination ownership/mode differs from explicit policy")
    _require_safe_regular_or_absent_entries(
        directory, _sidecars(name), "restore destination SQLite sidecar"
    )
    return metadata


def restore_backup(
    source: Path,
    destination: Path,
    *,
    staging_directory: Path,
    expected_uid: int,
    expected_gid: int,
    mode: int,
) -> None:
    source = _absolute_path(source, "backup source")
    destination = _absolute_path(destination, "restore destination")
    staging_directory = _absolute_path(staging_directory, "restore staging directory")
    destination, destination_parent, destination_name = _open_parent(
        destination, "restore destination"
    )
    if _path_is_beneath(staging_directory, destination.parent):
        os.close(destination_parent)
        fail("restore staging directory must not be beneath the destination parent")
    try:
        staging = _require_private_root_directory(
            staging_directory, "restore staging directory"
        )
    except BaseException:
        os.close(destination_parent)
        raise
    staging_metadata = os.fstat(staging)
    destination_parent_metadata = os.fstat(destination_parent)
    if (
        staging_metadata.st_ino == destination_parent_metadata.st_ino
        and staging_metadata.st_dev == destination_parent_metadata.st_dev
    ):
        os.close(staging)
        os.close(destination_parent)
        fail("restore staging directory must differ from the destination parent")
    if not _same_filesystem(staging_metadata, destination_parent_metadata):
        os.close(staging)
        os.close(destination_parent)
        fail("restore staging directory must use the destination filesystem")
    if expected_uid <= 0 or expected_gid < 0 or mode != PRIVATE_FILE_MODE:
        os.close(staging)
        os.close(destination_parent)
        fail("restore ownership/mode policy is invalid")
    if (
        destination_parent_metadata.st_uid != expected_uid
        or destination_parent_metadata.st_gid != expected_gid
        or stat.S_IMODE(destination_parent_metadata.st_mode) != PRIVATE_DIRECTORY_MODE
    ):
        os.close(staging)
        os.close(destination_parent)
        fail(
            "restore destination parent must match the expected owner/group "
            "with exact mode 0700"
        )
    try:
        source, source_parent, source_name = _open_parent(source, "backup source")
    except BaseException:
        os.close(staging)
        os.close(destination_parent)
        raise

    source_descriptor: int | None = None
    manifest_descriptor: int | None = None
    source_connection: sqlite3.Connection | None = None
    destination_connection: sqlite3.Connection | None = None
    temporary_name = _transaction_name(
        "restore", staging_directory, destination, str(source)
    )
    manage_destination_sidecars = False
    succeeded = False
    try:
        _cleanup_regular_entries(staging, (temporary_name, *_sidecars(temporary_name)))
        initial_destination = _restore_destination_metadata(
            destination_parent,
            destination_name,
            expected_uid=expected_uid,
            expected_gid=expected_gid,
            mode=mode,
        )
        if initial_destination is None:
            manage_destination_sidecars = True
            _cleanup_regular_entries(destination_parent, _sidecars(destination_name))

        source_descriptor, source_metadata = _require_regular_at(
            source_parent, source_name, "backup source"
        )
        manifest_name = f"{source_name}.json"
        manifest_descriptor, manifest_metadata = _require_regular_at(
            source_parent, manifest_name, "backup manifest"
        )
        for label, metadata in (
            ("backup", source_metadata),
            ("backup manifest", manifest_metadata),
        ):
            if (
                metadata.st_uid != ROOT_UID
                or metadata.st_gid != ROOT_GID
                or stat.S_IMODE(metadata.st_mode) != PRIVATE_FILE_MODE
            ):
                fail(f"{label} must be exactly root:root mode 0600")
        source_parent_metadata = os.fstat(source_parent)
        if (
            source_parent_metadata.st_uid != ROOT_UID
            or source_parent_metadata.st_gid != ROOT_GID
            or stat.S_IMODE(source_parent_metadata.st_mode) != PRIVATE_DIRECTORY_MODE
        ):
            fail("backup parent must be exactly root:root mode 0700")
        _reject_entries(source_parent, _sidecars(source_name), "backup SQLite sidecar")
        try:
            document = json.loads(
                _read_descriptor(manifest_descriptor, maximum=64 * 1024).decode("utf-8")
            )
        except (UnicodeError, json.JSONDecodeError) as error:
            fail(f"backup manifest is malformed: {type(error).__name__}")
        document = _validate_manifest(document)
        source_digest = _sha256_descriptor(source_descriptor)
        if (
            document.get("size_bytes") != source_metadata.st_size
            or document.get("sha256") != source_digest
        ):
            fail("backup bytes differ from the adjacent manifest")

        temporary_entries = (temporary_name, *_sidecars(temporary_name))
        _cleanup_regular_entries(staging, temporary_entries)
        temporary_descriptor = os.open(
            temporary_name,
            os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC | os.O_NOFOLLOW,
            PRIVATE_FILE_MODE,
            dir_fd=staging,
        )
        os.close(temporary_descriptor)
        source_connection = sqlite3.connect(
            _sqlite_path_for_descriptor(source_descriptor, read_only=True),
            uri=True,
            timeout=30,
        )
        destination_connection = sqlite3.connect(
            _sqlite_path_in_directory(staging, temporary_name),
            uri=True,
            timeout=30,
        )
        integrity(source_connection)
        source_connection.backup(destination_connection)
        integrity(destination_connection)
        if schema_version(destination_connection) != document.get(
            "database_schema_version"
        ):
            fail("restored database schema differs from the backup manifest")
        destination_connection.commit()
        destination_connection.close()
        destination_connection = None
        source_connection.close()
        source_connection = None

        temporary_descriptor, _ = _require_regular_at(
            staging,
            temporary_name,
            "temporary restore database",
            flags=os.O_RDWR,
        )
        try:
            os.fchown(temporary_descriptor, expected_uid, expected_gid)
            os.fchmod(temporary_descriptor, mode)
            os.fsync(temporary_descriptor)
            restored_digest = _sha256_descriptor(temporary_descriptor)
            restored_metadata = os.fstat(temporary_descriptor)
        finally:
            os.close(temporary_descriptor)

        destination_metadata = _restore_destination_metadata(
            destination_parent,
            destination_name,
            expected_uid=expected_uid,
            expected_gid=expected_gid,
            mode=mode,
        )
        if destination_metadata is not None:
            _checkpoint_existing_database(
                destination_parent, destination_name, "existing restore destination"
            )
        manage_destination_sidecars = True
        _cleanup_regular_entries(destination_parent, _sidecars(destination_name))
        os.replace(
            temporary_name,
            destination_name,
            src_dir_fd=staging,
            dst_dir_fd=destination_parent,
        )
        os.fsync(staging)
        os.fsync(destination_parent)
        _cleanup_regular_entries(staging, _sidecars(temporary_name))
        final_descriptor, final_metadata = _require_regular_at(
            destination_parent, destination_name, "restored database"
        )
        try:
            if (
                final_metadata.st_uid != expected_uid
                or final_metadata.st_gid != expected_gid
                or stat.S_IMODE(final_metadata.st_mode) != mode
                or final_metadata.st_size != restored_metadata.st_size
                or _sha256_descriptor(final_descriptor) != restored_digest
            ):
                fail("restored database metadata or bytes changed during publication")
        finally:
            os.close(final_descriptor)
        succeeded = True
        print(
            json.dumps(
                {
                    "restored_sha256": restored_digest,
                    "restored_size_bytes": restored_metadata.st_size,
                },
                sort_keys=True,
            )
        )
    finally:
        try:
            try:
                _close_database(destination_connection)
                _close_database(source_connection)
            finally:
                if manifest_descriptor is not None:
                    os.close(manifest_descriptor)
                if source_descriptor is not None:
                    os.close(source_descriptor)
        finally:
            try:
                try:
                    _cleanup_regular_entries(
                        staging, (temporary_name, *_sidecars(temporary_name))
                    )
                finally:
                    if manage_destination_sidecars:
                        _cleanup_regular_entries(
                            destination_parent, _sidecars(destination_name)
                        )
                        if succeeded:
                            os.fsync(destination_parent)
            finally:
                os.close(staging)
                os.close(destination_parent)
                os.close(source_parent)


def main() -> int:
    if os.geteuid() != 0:
        fail("run as root while the API and worker are stopped")
    os.umask(0o077)
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    backup_parser = subparsers.add_parser("backup")
    backup_parser.add_argument("--source", type=Path, required=True)
    backup_parser.add_argument("--destination", type=Path, required=True)
    backup_parser.add_argument("--release", required=True)
    restore_parser = subparsers.add_parser("restore")
    restore_parser.add_argument("--source", type=Path, required=True)
    restore_parser.add_argument("--destination", type=Path, required=True)
    restore_parser.add_argument("--staging-directory", type=Path, required=True)
    restore_parser.add_argument("--uid", type=int, required=True)
    restore_parser.add_argument("--gid", type=int, required=True)
    restore_parser.add_argument(
        "--mode", type=lambda value: int(value, 8), required=True
    )
    args = parser.parse_args()
    if args.action == "backup":
        create_backup(args.source, args.destination, args.release)
    else:
        restore_backup(
            args.source,
            args.destination,
            staging_directory=args.staging_directory,
            expected_uid=args.uid,
            expected_gid=args.gid,
            mode=args.mode,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
