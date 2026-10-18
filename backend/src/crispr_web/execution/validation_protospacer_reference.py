"""Independently authenticate exact-search reported catalogue rows with fixed budgets.

Integrity hashing is bounded. Indexed membership validation cannot prove that the
scientific search omitted no hits. No Integration or scientific-package imports.
"""
import hashlib
import os
from pathlib import Path
import sqlite3
import stat
import time

from .validation_annotations import _strict_json

MAX_DATABASE_BYTES = 4 * 1024**3
VALIDATION_SECONDS = 120
SQL_WORK_LIMIT = 2_000_000


def _open_regular(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise ValueError("Spacer reference must be a regular file")
    return os.fdopen(descriptor, "rb")


def _validate_reference(directory, expected_digest, reported, hits, checkpoint=None):
    if directory is None:
        raise ValueError("Pinned spacer reference unavailable")
    root = Path(directory).absolute()
    for item in (root, *root.parents, root / "reference.json", root / "spacers.sqlite3"):
        if item.is_symlink():
            raise ValueError("Symlink in pinned spacer reference")
    deadline = time.monotonic() + VALIDATION_SECONDS

    def check():
        if checkpoint is not None:
            checkpoint()
        if time.monotonic() > deadline:
            raise ValueError("Spacer reference validation timed out")

    def signature(info):
        return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns

    with _open_regular(root / "reference.json") as metadata:
        payload = metadata.read(16_385)
    if len(payload) > 16_384:
        raise ValueError("Oversized spacer metadata")
    reference = _strict_json(payload)
    if reference != reported or reference.get("database_sha256") != expected_digest:
        raise ValueError("External spacer reference mismatch")
    filename = root / "spacers.sqlite3"
    with _open_regular(filename) as source:
        initial = os.fstat(source.fileno())
        if not 0 < initial.st_size <= MAX_DATABASE_BYTES:
            raise ValueError("Invalid spacer catalogue size")
        count = 0
        digest = hashlib.sha256()
        while chunk := source.read(min(262_144, MAX_DATABASE_BYTES - count + 1)):
            check()
            count += len(chunk)
            if count > MAX_DATABASE_BYTES or count > initial.st_size:
                raise ValueError("Spacer reference grew during read")
            digest.update(chunk)
        if digest.hexdigest() != expected_digest or signature(initial) != signature(os.fstat(source.fileno())) or signature(initial) != signature(filename.stat()):
            raise ValueError("Pinned spacer bytes changed")
        pinned = Path(f"/proc/self/fd/{source.fileno()}")
        if not pinned.exists():
            raise ValueError("Cannot bind SQLite to hashed file descriptor")
        database = sqlite3.connect(pinned.as_uri() + "?mode=ro&immutable=1", uri=True)
        caught, steps = [], 0
        try:
            database.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 32_768)
            database.execute("PRAGMA query_only=ON")
            database.execute("PRAGMA trusted_schema=OFF")

            def progress():
                nonlocal steps
                steps += 500
                try:
                    check()
                    if steps > SQL_WORK_LIMIT:
                        raise ValueError("Spacer reference SQL budget exhausted")
                except BaseException as error:
                    caught.append(error)
                    return 1
                return 0

            database.set_progress_handler(progress, 500)

            def read(sql, parameters=()):
                check()
                try:
                    return database.execute(sql, parameters).fetchall()
                except sqlite3.Error as error:
                    if caught:
                        raise caught[0]
                    raise ValueError("Malformed spacer reference database") from error

            if read("SELECT name,type FROM sqlite_master WHERE name IN ('metadata','spacers') ORDER BY name LIMIT 3") != [("metadata", "table"), ("spacers", "table")]:
                raise ValueError("Invalid spacer database tables")
            metadata_shape = read("SELECT typeof(value),length(CAST(value AS BLOB)) FROM metadata LIMIT 2")
            if len(metadata_shape) != 1 or metadata_shape[0][0] != "text" or not 1 <= metadata_shape[0][1] <= 16_384:
                raise ValueError("Invalid internal spacer metadata size/type")
            internal = _strict_json(read("SELECT value FROM metadata LIMIT 2")[0][0])
            if internal != {key: val for key, val in reference.items() if key != "database_sha256"}:
                raise ValueError("Internal/external spacer metadata mismatch")
            for spacer_id, dna in hits.items():
                plan = read("EXPLAIN QUERY PLAN SELECT sequence,status FROM spacers WHERE id=? LIMIT 2", (spacer_id,))
                if len(plan) != 1 or not plan[0][3].startswith(("SEARCH spacers USING INDEX", "SEARCH spacers USING COVERING INDEX")):
                    raise ValueError("Spacer membership lookup is not indexed")
                lengths = read("SELECT typeof(sequence),length(CAST(sequence AS BLOB)),typeof(status),length(CAST(status AS BLOB)) FROM spacers WHERE id=? LIMIT 2", (spacer_id,))
                if len(lengths) != 1 or lengths[0][0] != "text" or not 16 <= lengths[0][1] <= 1000 or lengths[0][2] != "text" or not 0 < lengths[0][3] <= 30:
                    raise ValueError("Invalid reported catalogue row")
                if read("SELECT sequence,status FROM spacers WHERE id=? LIMIT 2", (spacer_id,)) != [(dna, "eligible")]:
                    raise ValueError("Reported match does not belong to pinned catalogue")
            check()
            if signature(initial) != signature(os.fstat(source.fileno())) or signature(initial) != signature(filename.stat()):
                raise ValueError("Spacer reference changed during validation")
        finally:
            database.close()


def validate_reference(*args, **kwargs):
    try:
        return _validate_reference(*args, **kwargs)
    except sqlite3.Error as error:
        raise ValueError("Invalid spacer catalogue database") from error
