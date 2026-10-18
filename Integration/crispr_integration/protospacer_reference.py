"""Bounded membership binding for reported hits, not a completeness replay."""
from pathlib import Path
from contextlib import closing
import hashlib
import os
import sqlite3
import stat
import time

from .protospacer import _json

MAX_BYTES = 4 * 1024**3
MAX_SECONDS = 120
MAX_SQL_STEPS = 2_000_000


def _open_regular(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise ValueError("Spacer reference must be a regular file")
    return os.fdopen(descriptor, "rb")


def _bind_catalogue(directory, digest, reference, spacers, checkpoint=None):
    if directory is None:
        raise ValueError("Missing pinned spacer catalogue")
    directory = Path(directory).absolute()
    if any(path.is_symlink() for path in (directory, *directory.parents)):
        raise ValueError("Spacer catalogue contains symlinks")
    start = time.monotonic()

    def step():
        if checkpoint:
            checkpoint()
        if time.monotonic() - start > MAX_SECONDS:
            raise ValueError("Spacer catalogue validation deadline exceeded")

    def identity(stat):
        return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns

    metadata = directory / "reference.json"
    if metadata.is_symlink() or not metadata.is_file() or metadata.stat().st_size > 16_384:
        raise ValueError("Invalid spacer reference metadata")
    with _open_regular(metadata) as stream:
        raw = stream.read(16_385)
    if len(raw) > 16_384 or _json(raw) != reference or reference.get("database_sha256") != digest:
        raise ValueError("Spacer reference metadata differs from pinned result")
    path = directory / "spacers.sqlite3"
    if path.is_symlink() or not path.is_file():
        raise ValueError("Invalid spacer catalogue file")
    with _open_regular(path) as stream:
        before = os.fstat(stream.fileno())
        if not 0 < before.st_size <= MAX_BYTES:
            raise ValueError("Spacer catalogue exceeds byte budget")
        checksum, total = hashlib.sha256(), 0
        while block := stream.read(min(1024 * 1024, MAX_BYTES - total + 1)):
            step()
            total += len(block)
            if total > MAX_BYTES or total > before.st_size:
                raise ValueError("Spacer catalogue grew during hashing")
            checksum.update(block)
        if checksum.hexdigest() != digest or identity(before) != identity(os.fstat(stream.fileno())) or identity(before) != identity(path.stat()):
            raise ValueError("Spacer catalogue integrity mismatch")
        # Keep the hashed inode open while SQLite opens that same inode on Linux.
        pinned = Path(f"/proc/self/fd/{stream.fileno()}")
        if not pinned.exists():
            raise ValueError("Immutable file-descriptor SQLite binding unavailable")
        with closing(sqlite3.connect(pinned.as_uri() + "?mode=ro&immutable=1", uri=True)) as connection:
            try:
                connection.execute("PRAGMA trusted_schema=OFF")
                connection.execute("PRAGMA query_only=ON")
                connection.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 32_768)
                work, interrupted = 0, []

                def progress():
                    nonlocal work
                    work += 1000
                    try:
                        step()
                        if work > MAX_SQL_STEPS:
                            raise ValueError("Spacer catalogue SQL work budget exceeded")
                    except BaseException as error:
                        interrupted.append(error)
                        return 1
                    return 0

                connection.set_progress_handler(progress, 1000)

                def query(sql, args=()):
                    step()
                    try:
                        return connection.execute(sql, args).fetchall()
                    except sqlite3.Error as error:
                        if interrupted:
                            raise interrupted[0] from None
                        raise ValueError("Invalid spacer catalogue SQL structure") from error

                if query("SELECT type FROM sqlite_master WHERE name='metadata' LIMIT 2") != [("table",)] or query("SELECT type FROM sqlite_master WHERE name='spacers' LIMIT 2") != [("table",)]:
                    raise ValueError("Catalogue metadata and spacers must be tables")
                sizes = query("SELECT typeof(value),length(CAST(value AS BLOB)) FROM metadata LIMIT 2")
                if len(sizes) != 1 or sizes[0][0] != "text" or not 0 < sizes[0][1] <= 16_384:
                    raise ValueError("Malformed catalogue metadata")
                internal = _json(query("SELECT value FROM metadata LIMIT 2")[0][0])
                if internal != {key: value for key, value in reference.items() if key != "database_sha256"}:
                    raise ValueError("Internal catalogue reference differs")
                for key, sequence in spacers.items():
                    plan = query("EXPLAIN QUERY PLAN SELECT sequence,status FROM spacers WHERE id=? LIMIT 2", (key,))
                    if not plan or any("SEARCH spacers USING INDEX" not in row[3] and "SEARCH spacers USING COVERING INDEX" not in row[3] for row in plan):
                        raise ValueError("Spacer lookup requires an index")
                    shape = query("SELECT typeof(sequence),length(CAST(sequence AS BLOB)),typeof(status),length(CAST(status AS BLOB)) FROM spacers WHERE id=? LIMIT 2", (key,))
                    if len(shape) != 1 or shape[0][0] != "text" or not 16 <= shape[0][1] <= 1000 or shape[0][2] != "text" or not 1 <= shape[0][3] <= 30:
                        raise ValueError("Missing or oversized reported spacer")
                    if query("SELECT sequence,status FROM spacers WHERE id=? LIMIT 2", (key,)) != [(sequence, "eligible")]:
                        raise ValueError("Reported spacer differs from pinned catalogue")
                step()
                if identity(before) != identity(os.fstat(stream.fileno())) or identity(before) != identity(path.stat()):
                    raise ValueError("Spacer catalogue changed during validation")
            finally:
                connection.close()


def bind_catalogue(*args, **kwargs):
    try:
        return _bind_catalogue(*args, **kwargs)
    except sqlite3.Error as error:
        raise ValueError("Invalid spacer catalogue database") from error
