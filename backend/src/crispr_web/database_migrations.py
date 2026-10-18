"""Transactional migration of the historical three-mode queue to CRISPRloci v2."""

from __future__ import annotations

import re
import sqlite3

from .models import ScientificMode


def expand_analysis_modes(connection: sqlite3.Connection) -> None:
    """Rebuild SQLite's CHECK constraint while retaining jobs and artifact FKs.

    Called inside Store.initialize's IMMEDIATE transaction. No writes escape the
    transaction; queued/running jobs, token digests and registered downloads survive.
    """
    sql = connection.execute(
        "SELECT sql FROM sqlite_master WHERE name='jobs'"
    ).fetchone()[0]
    constraint = (
        "CHECK (mode IN ("
        + ",".join(repr(mode.value) for mode in ScientificMode)
        + "))"
    )
    pattern = r"CHECK\s*\(mode\s+IN\s*\([^)]*\)\)"
    existing = re.search(pattern, sql, re.IGNORECASE)
    if existing is not None and set(re.findall(r"'([^']+)'", existing.group())) == {
        mode.value for mode in ScientificMode
    }:
        return
    revised, count = re.subn(pattern, constraint, sql, flags=re.IGNORECASE)
    if count != 1:
        raise RuntimeError(
            "Unsupported historical jobs mode constraint; database was not migrated"
        )
    revised = re.sub(
        r'^CREATE TABLE\s+(?:"jobs"|jobs)',
        "CREATE TABLE jobs_v4",
        revised,
        count=1,
        flags=re.IGNORECASE,
    )
    artifact_sql = connection.execute(
        "SELECT sql FROM sqlite_master WHERE name='artifacts'"
    ).fetchone()[0]
    indexes = [
        row[0]
        for row in connection.execute(
            "SELECT sql FROM sqlite_master WHERE type='index' AND tbl_name IN ('jobs','artifacts') AND sql IS NOT NULL"
        )
    ]
    connection.execute(revised)
    connection.execute("INSERT INTO jobs_v4 SELECT * FROM jobs")
    connection.execute(
        "CREATE TEMP TABLE migration_artifacts AS SELECT * FROM artifacts"
    )
    connection.execute("DROP TABLE artifacts")
    connection.execute("DROP TABLE jobs")
    connection.execute("ALTER TABLE jobs_v4 RENAME TO jobs")
    connection.execute(artifact_sql)
    connection.execute("INSERT INTO artifacts SELECT * FROM migration_artifacts")
    connection.execute("DROP TABLE migration_artifacts")
    for statement in indexes:
        connection.execute(statement)
    if connection.execute("PRAGMA foreign_key_check").fetchone() is not None:
        raise RuntimeError("Queue migration failed foreign-key verification")
