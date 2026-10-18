"""SQLite-backed durable queue and capability-protected job metadata."""

from __future__ import annotations

from contextlib import closing, contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import hmac
import json
import os
import secrets
import sqlite3
import time
from typing import Any, Iterator, Mapping, Sequence

from .capacity import has_capacity
from .config import Settings
from .database_migrations import expand_analysis_modes
from .security import client_rate_identity, token_matches
from .scientific_options import validate_tool_options


_SCHEMA_VERSION = 7
_GLOBAL_RATE_IDENTITY = "global-submissions-v1"
_WAL_BUSY_RETRY_SECONDS = 15.0
_WAL_BUSY_INITIAL_DELAY_SECONDS = 0.005
_WAL_BUSY_MAX_DELAY_SECONDS = 0.1
_MAX_DATABASE_BYTES = 1_000_000_000
_WAL_JOURNAL_SIZE_LIMIT_BYTES = 64_000_000

TERMINAL_STATUSES = frozenset(
    {"completed", "completed_no_eligible_groups", "failed", "cancelled"}
)


class QueueCapacityError(RuntimeError):
    pass


class ClientCapacityError(RuntimeError):
    pass


class SubmissionRateLimitError(RuntimeError):
    def __init__(self, retry_after_seconds: int):
        self.retry_after_seconds = retry_after_seconds
        super().__init__("submission rate limit exceeded")


class LeaseLostError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class ClaimedJob:
    job_id: str
    mode: str
    attempt: int
    started_at: float
    category_policy: str
    spacer_distance: int
    bias_corrections_requested: bool
    bias_corrections_effective: bool
    option_warnings: tuple[str, ...]
    leader_flank_length: int = 500
    tracr_model_type: str = "II"
    molecule: str = "DNA"
    spacer_database_sha256: str | None = None
    association_reference_sha256: str | None = None
    association_grouping: str | None = None
    array_map_sha256: str | None = None
    association_groups_sha256: str | None = None
    repeat_reference_sha256: str | None = None
    viral_reference_sha256: str | None = None
    viral_max_mismatches: int = 2
    tool_options: dict = field(default_factory=dict)


def utc_iso(epoch: float | None) -> str | None:
    if epoch is None:
        return None
    return (
        datetime.fromtimestamp(epoch, timezone.utc).isoformat().replace("+00:00", "Z")
    )


class Store:
    def __init__(self, settings: Settings):
        self.settings = settings
        self._prepare_paths()
        self.initialize()

    def _prepare_paths(self) -> None:
        self.settings.data_root.mkdir(parents=True, exist_ok=True, mode=0o700)
        self.settings.jobs_root.mkdir(parents=True, exist_ok=True, mode=0o700)
        os.chmod(self.settings.data_root, 0o700)
        os.chmod(self.settings.jobs_root, 0o700)
        self.settings.database_path.parent.mkdir(
            parents=True, exist_ok=True, mode=0o700
        )

    def connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(
            self.settings.database_path,
            timeout=15,
            isolation_level=None,
        )
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        connection.execute("PRAGMA busy_timeout = 15000")
        connection.execute("PRAGMA synchronous = FULL")
        page_size = int(connection.execute("PRAGMA page_size").fetchone()[0])
        maximum_pages = max(1, _MAX_DATABASE_BYTES // page_size)
        configured_pages = int(
            connection.execute(f"PRAGMA max_page_count = {maximum_pages}").fetchone()[0]
        )
        if configured_pages > maximum_pages:
            connection.close()
            raise RuntimeError("database already exceeds its fixed byte ceiling")
        connection.execute(
            f"PRAGMA journal_size_limit = {_WAL_JOURNAL_SIZE_LIMIT_BYTES}"
        )
        connection.execute("PRAGMA wal_autocheckpoint = 1000")
        return connection

    @contextmanager
    def transaction(self, *, immediate: bool = False) -> Iterator[sqlite3.Connection]:
        connection = self.connect()
        try:
            connection.execute("BEGIN IMMEDIATE" if immediate else "BEGIN")
            yield connection
            connection.commit()
        except BaseException:
            connection.rollback()
            raise
        finally:
            connection.close()

    def initialize(self) -> None:
        self._enable_wal()
        with self.transaction(immediate=True) as connection:
            metadata_exists = connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'metadata'"
            ).fetchone()
            current_schema_version = 0
            if metadata_exists is not None:
                version_row = connection.execute(
                    "SELECT value FROM metadata WHERE name = 'schema_version'"
                ).fetchone()
                if version_row is not None:
                    try:
                        current_schema_version = int(version_row["value"])
                    except (TypeError, ValueError) as error:
                        raise RuntimeError(
                            "database schema version is invalid"
                        ) from error
                    if current_schema_version > _SCHEMA_VERSION:
                        raise RuntimeError(
                            "database schema is newer than this backend release"
                        )
            for statement in (
                """
                CREATE TABLE IF NOT EXISTS jobs (
                    job_id TEXT PRIMARY KEY,
                    token_digest TEXT NOT NULL,
                    client_digest TEXT,
                    status TEXT NOT NULL CHECK (status IN (
                        'queued', 'running', 'completed',
                        'completed_no_eligible_groups', 'failed', 'cancelled'
                    )),
                    mode TEXT NOT NULL CHECK (mode IN (
                        'detection', 'orientation', 'reconstruction', 'cas', 'tracrrna', 'loci', 'leader', 'repeats', 'repeat_context', 'protospacer'
                    )),
                    created_at REAL NOT NULL,
                    updated_at REAL NOT NULL,
                    started_at REAL,
                    finished_at REAL,
                    expires_at REAL,
                    filename TEXT NOT NULL,
                    record_count INTEGER NOT NULL,
                    base_count INTEGER NOT NULL,
                    input_sha256 TEXT NOT NULL,
                    source_ids_json TEXT NOT NULL,
                    category_policy TEXT NOT NULL,
                    spacer_distance INTEGER NOT NULL,
                    bias_corrections_requested INTEGER NOT NULL,
                    bias_corrections_effective INTEGER NOT NULL,
                    option_warnings_json TEXT NOT NULL,
                    attempt INTEGER NOT NULL DEFAULT 0,
                    worker_id TEXT,
                    lease_expires_at REAL,
                    cancel_requested INTEGER NOT NULL DEFAULT 0,
                    summary_json TEXT,
                    error_code TEXT,
                    error_message TEXT
                );
                CREATE INDEX IF NOT EXISTS jobs_queue_idx
                    ON jobs(status, created_at);
                CREATE INDEX IF NOT EXISTS jobs_expiry_idx
                    ON jobs(expires_at);

                CREATE TABLE IF NOT EXISTS artifacts (
                    job_id TEXT NOT NULL REFERENCES jobs(job_id) ON DELETE CASCADE,
                    artifact_id TEXT NOT NULL,
                    relative_path TEXT NOT NULL,
                    name TEXT NOT NULL,
                    size_bytes INTEGER NOT NULL,
                    sha256 TEXT NOT NULL,
                    media_type TEXT NOT NULL,
                    PRIMARY KEY(job_id, artifact_id),
                    UNIQUE(job_id, relative_path)
                );

                CREATE TABLE IF NOT EXISTS rate_events (
                    client_digest TEXT NOT NULL,
                    submitted_at REAL NOT NULL
                );
                CREATE INDEX IF NOT EXISTS rate_events_lookup_idx
                    ON rate_events(client_digest, submitted_at);
                CREATE INDEX IF NOT EXISTS rate_events_expiry_idx
                    ON rate_events(submitted_at);

                CREATE TABLE IF NOT EXISTS service_state (
                    name TEXT PRIMARY KEY,
                    value TEXT NOT NULL,
                    updated_at REAL NOT NULL
                );

                CREATE TABLE IF NOT EXISTS metadata (
                    name TEXT PRIMARY KEY,
                    value TEXT NOT NULL
                );
                """
            ).split(";"):
                if statement.strip():
                    connection.execute(statement)
            columns = {
                row["name"] for row in connection.execute("PRAGMA table_info(jobs)")
            }
            migrations = {
                "tool_options_json": "TEXT NOT NULL DEFAULT '{}'",
                "spacer_database_sha256": "TEXT",
                "association_reference_sha256": "TEXT",
                "association_grouping": "TEXT",
                "association_groups_sha256": "TEXT",
                "array_map_sha256": "TEXT",
                "repeat_reference_sha256": "TEXT",
                "viral_reference_sha256": "TEXT",
                "viral_max_mismatches": "INTEGER NOT NULL DEFAULT 2",
                "leader_flank_length": "INTEGER NOT NULL DEFAULT 500",
                "molecule": "TEXT NOT NULL DEFAULT 'DNA'",
                "tracr_model_type": "TEXT NOT NULL DEFAULT 'II'",
                "client_digest": "TEXT",
                "input_sha256": "TEXT NOT NULL DEFAULT ''",
                "source_ids_json": "TEXT NOT NULL DEFAULT '[]'",
                "category_policy": "TEXT NOT NULL DEFAULT 'bona_fide_possible'",
                "spacer_distance": "INTEGER NOT NULL DEFAULT 1",
                "bias_corrections_requested": "INTEGER NOT NULL DEFAULT 1",
                "bias_corrections_effective": "INTEGER NOT NULL DEFAULT 1",
                "option_warnings_json": "TEXT NOT NULL DEFAULT '[]'",
            }
            for name, declaration in migrations.items():
                if name not in columns:
                    connection.execute(
                        f"ALTER TABLE jobs ADD COLUMN {name} {declaration}"
                    )
            expand_analysis_modes(connection)
            connection.execute(
                "CREATE INDEX IF NOT EXISTS jobs_client_active_idx "
                "ON jobs(client_digest, status) WHERE client_digest IS NOT NULL"
            )
            connection.execute(
                "UPDATE jobs SET client_digest = NULL "
                "WHERE client_digest IS NOT NULL AND status IN "
                "('completed', 'completed_no_eligible_groups', 'failed', 'cancelled')"
            )
            if current_schema_version < _SCHEMA_VERSION:
                connection.execute(
                    "INSERT INTO metadata(name, value) VALUES('schema_version', ?) "
                    "ON CONFLICT(name) DO UPDATE SET value = excluded.value",
                    (str(_SCHEMA_VERSION),),
                )
            connection.execute(
                "INSERT OR IGNORE INTO metadata(name, value) VALUES('rate_salt', ?)",
                (secrets.token_hex(32),),
            )
        os.chmod(self.settings.database_path, 0o600)

    def _enable_wal(self) -> None:
        deadline = time.monotonic() + _WAL_BUSY_RETRY_SECONDS
        delay = _WAL_BUSY_INITIAL_DELAY_SECONDS
        while True:
            try:
                with closing(self.connect()) as connection:
                    row = connection.execute("PRAGMA journal_mode = WAL").fetchone()
            except sqlite3.OperationalError as error:
                error_code = getattr(error, "sqlite_errorcode", None)
                primary_error_code = (
                    error_code & 0xFF if isinstance(error_code, int) else None
                )
                if primary_error_code != sqlite3.SQLITE_BUSY:
                    raise
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise
                time.sleep(min(delay, remaining))
                delay = min(delay * 2, _WAL_BUSY_MAX_DELAY_SECONDS)
                continue
            if row is None or row[0] != "wal":
                raise RuntimeError("database refused WAL journal mode")
            return

    def _client_digest(self, connection: sqlite3.Connection, client_ip: str) -> str:
        row = connection.execute(
            "SELECT value FROM metadata WHERE name = 'rate_salt'"
        ).fetchone()
        if row is None:
            raise RuntimeError("database rate salt is missing")
        key = bytes.fromhex(row["value"])
        identity = client_rate_identity(client_ip)
        return hmac.new(key, identity.encode("ascii"), hashlib.sha256).hexdigest()

    @staticmethod
    def _enforce_submission_rate(
        connection: sqlite3.Connection,
        identity: str,
        *,
        cutoff: float,
        limit: int,
        window_seconds: int,
        timestamp: float,
    ) -> None:
        row = connection.execute(
            "SELECT COUNT(*) AS count, MIN(submitted_at) AS oldest "
            "FROM rate_events WHERE client_digest = ? AND submitted_at >= ?",
            (identity, cutoff),
        ).fetchone()
        if row["count"] < limit:
            return
        retry = max(
            1,
            int(row["oldest"] + window_seconds - timestamp + 0.999),
        )
        raise SubmissionRateLimitError(retry)

    def purge_expired_rate_events(self, now: float | None = None) -> int:
        timestamp = time.time() if now is None else now
        cutoff = timestamp - max(
            self.settings.rate_limit_window_seconds,
            self.settings.global_rate_limit_window_seconds,
        )
        with self.transaction(immediate=True) as connection:
            return connection.execute(
                "DELETE FROM rate_events WHERE submitted_at < ?", (cutoff,)
            ).rowcount

    def preflight_submission(self, client_ip: str, now: float | None = None) -> None:
        """Reject obvious persisted limits before parsing or materializing FASTA."""

        timestamp = time.time() if now is None else now
        client_cutoff = timestamp - self.settings.rate_limit_window_seconds
        global_cutoff = timestamp - self.settings.global_rate_limit_window_seconds
        with self.connect() as connection:
            client_digest = self._client_digest(connection, client_ip)
            self._enforce_submission_rate(
                connection,
                client_digest,
                cutoff=client_cutoff,
                limit=self.settings.rate_limit_submissions,
                window_seconds=self.settings.rate_limit_window_seconds,
                timestamp=timestamp,
            )
            self._enforce_submission_rate(
                connection,
                _GLOBAL_RATE_IDENTITY,
                cutoff=global_cutoff,
                limit=self.settings.global_rate_limit_submissions,
                window_seconds=self.settings.global_rate_limit_window_seconds,
                timestamp=timestamp,
            )
            client_active = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs WHERE client_digest = ? AND "
                "(status = 'running' OR (status = 'queued' AND expires_at > ?))",
                (client_digest, timestamp),
            ).fetchone()["count"]
            if client_active >= self.settings.max_active_jobs_per_client:
                raise ClientCapacityError(
                    "this client already has the maximum active jobs"
                )
            active = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs WHERE status = 'running' "
                "OR (status = 'queued' AND expires_at > ?)",
                (timestamp,),
            ).fetchone()["count"]
            if active >= self.settings.max_active_jobs:
                raise QueueCapacityError(
                    "the compute service is at active-job capacity"
                )
            queued = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs "
                "WHERE status = 'queued' AND expires_at > ?",
                (timestamp,),
            ).fetchone()["count"]
            if queued >= self.settings.max_queued_jobs:
                raise QueueCapacityError("the compute queue is full")

    def create_job(
        self,
        *,
        job_id: str,
        access_token_digest: str,
        mode: str,
        filename: str,
        record_count: int,
        base_count: int,
        input_sha256: str,
        source_ids: Sequence[str],
        category_policy: str,
        spacer_distance: int,
        bias_corrections_requested: bool,
        bias_corrections_effective: bool,
        option_warnings: Sequence[str],
        client_ip: str,
        leader_flank_length: int = 500,
        tracr_model_type: str = "II",
        molecule: str = "DNA",
        spacer_database_sha256: str | None = None,
        association_reference_sha256: str | None = None,
        association_grouping: str | None = None,
        association_groups_sha256: str | None = None,
        array_map_sha256: str | None = None,
        repeat_reference_sha256: str | None = None,
        viral_reference_sha256: str | None = None,
        viral_max_mismatches: int = 2,
        tool_options: dict | None = None,
        now: float | None = None,
    ) -> None:
        selected_tool_options = validate_tool_options(tool_options or {}, mode)
        timestamp = time.time() if now is None else now
        client_cutoff = timestamp - self.settings.rate_limit_window_seconds
        global_cutoff = timestamp - self.settings.global_rate_limit_window_seconds
        purge_cutoff = min(client_cutoff, global_cutoff)
        with self.transaction(immediate=True) as connection:
            connection.execute(
                "DELETE FROM rate_events WHERE submitted_at < ?", (purge_cutoff,)
            )
            client_digest = self._client_digest(connection, client_ip)
            self._enforce_submission_rate(
                connection,
                client_digest,
                cutoff=client_cutoff,
                limit=self.settings.rate_limit_submissions,
                window_seconds=self.settings.rate_limit_window_seconds,
                timestamp=timestamp,
            )
            self._enforce_submission_rate(
                connection,
                _GLOBAL_RATE_IDENTITY,
                cutoff=global_cutoff,
                limit=self.settings.global_rate_limit_submissions,
                window_seconds=self.settings.global_rate_limit_window_seconds,
                timestamp=timestamp,
            )
            self._expire_queued_in_transaction(connection, timestamp)
            client_active = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs "
                "WHERE client_digest = ? AND status IN ('queued', 'running')",
                (client_digest,),
            ).fetchone()["count"]
            if client_active >= self.settings.max_active_jobs_per_client:
                raise ClientCapacityError(
                    "this client already has the maximum active jobs"
                )
            active = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs "
                "WHERE status IN ('queued', 'running')"
            ).fetchone()["count"]
            if active >= self.settings.max_active_jobs:
                raise QueueCapacityError(
                    "the compute service is at active-job capacity"
                )
            queued = connection.execute(
                "SELECT COUNT(*) AS count FROM jobs WHERE status = 'queued'"
            ).fetchone()["count"]
            if queued >= self.settings.max_queued_jobs:
                raise QueueCapacityError("the compute queue is full")
            connection.execute(
                """
                INSERT INTO jobs(
                    job_id, token_digest, client_digest, status, mode, created_at,
                    updated_at, expires_at, filename, record_count, base_count,
                    input_sha256, source_ids_json, category_policy, spacer_distance,
                    bias_corrections_requested, bias_corrections_effective,
                    option_warnings_json, tracr_model_type, leader_flank_length, molecule, spacer_database_sha256, viral_reference_sha256, viral_max_mismatches, repeat_reference_sha256, association_reference_sha256, association_grouping, association_groups_sha256, array_map_sha256, tool_options_json
                ) VALUES (?, ?, ?, 'queued', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    job_id,
                    access_token_digest,
                    client_digest,
                    mode,
                    timestamp,
                    timestamp,
                    timestamp
                    + self.settings.retention_seconds
                    + (self.settings.max_queued_jobs + 1)
                    * self.settings.job_timeout_seconds,
                    filename,
                    record_count,
                    base_count,
                    input_sha256,
                    json.dumps(list(source_ids), separators=(",", ":")),
                    category_policy,
                    spacer_distance,
                    int(bias_corrections_requested),
                    int(bias_corrections_effective),
                    json.dumps(list(option_warnings), separators=(",", ":")),
                    tracr_model_type,
                    leader_flank_length,
                    molecule,
                    spacer_database_sha256, viral_reference_sha256, viral_max_mismatches, repeat_reference_sha256, association_reference_sha256, association_grouping, association_groups_sha256, array_map_sha256,
                    json.dumps(selected_tool_options, separators=(",", ":"), allow_nan=False),
                ),
            )
            connection.executemany(
                "INSERT INTO rate_events(client_digest, submitted_at) VALUES(?, ?)",
                (
                    (client_digest, timestamp),
                    (_GLOBAL_RATE_IDENTITY, timestamp),
                ),
            )

    def registered_job_ids(self) -> set[str]:
        """Load lightweight queue membership once for filesystem reconciliation."""

        with self.connect() as connection:
            return {
                str(row["job_id"])
                for row in connection.execute("SELECT job_id FROM jobs")
            }

    def get_job(self, job_id: str) -> dict[str, Any] | None:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT * FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
            if row is None:
                return None
            job = dict(row)
            job["tool_options"] = validate_tool_options(json.loads(job["tool_options_json"]), job["mode"])
            try:
                job["source_ids"] = json.loads(job["source_ids_json"])
                job["option_warnings"] = json.loads(job["option_warnings_json"])
            except (TypeError, json.JSONDecodeError):
                job["source_ids"] = []
                job["option_warnings"] = [
                    "Stored provenance metadata could not be decoded."
                ]
            artifacts = connection.execute(
                "SELECT artifact_id, relative_path, name, size_bytes, sha256, media_type "
                "FROM artifacts WHERE job_id = ? ORDER BY name, artifact_id",
                (job_id,),
            ).fetchall()
            job["artifacts"] = [dict(item) for item in artifacts]
            if job["summary_json"]:
                try:
                    job["summary"] = json.loads(job["summary_json"])
                except json.JSONDecodeError:
                    job["summary"] = None
            else:
                job["summary"] = None
            if job["status"] == "queued":
                job["queue_position"] = connection.execute(
                    "SELECT COUNT(*) AS count FROM jobs "
                    "WHERE status = 'queued' AND "
                    "(created_at < ? OR (created_at = ? AND job_id <= ?))",
                    (job["created_at"], job["created_at"], job_id),
                ).fetchone()["count"]
            else:
                job["queue_position"] = None
            return job

    def authorize(self, job_id: str, token: str, now: float | None = None) -> bool:
        timestamp = time.time() if now is None else now
        with self.connect() as connection:
            row = connection.execute(
                "SELECT token_digest, expires_at FROM jobs WHERE job_id = ?",
                (job_id,),
            ).fetchone()
        return bool(
            row is not None
            and row["expires_at"] is not None
            and float(row["expires_at"]) > timestamp
            and token_matches(token, row["token_digest"])
        )

    def request_cancel(self, job_id: str) -> bool:
        now = time.time()
        expiry = now + self.settings.retention_seconds
        with self.transaction(immediate=True) as connection:
            row = connection.execute(
                "SELECT status, expires_at FROM jobs WHERE job_id = ?", (job_id,)
            ).fetchone()
            if row is None:
                raise KeyError(job_id)
            if row["expires_at"] is None or float(row["expires_at"]) <= now:
                raise KeyError(job_id)
            if row["status"] == "queued":
                cursor = connection.execute(
                    "UPDATE jobs SET status = 'cancelled', client_digest = NULL, cancel_requested = 1, "
                    "updated_at = ?, finished_at = ?, expires_at = ? WHERE job_id = ?",
                    (now, now, expiry, job_id),
                )
                return cursor.rowcount == 1
            if row["status"] == "running":
                connection.execute(
                    "UPDATE jobs SET cancel_requested = 1, updated_at = ? "
                    "WHERE job_id = ? AND cancel_requested = 0",
                    (now, job_id),
                )
            return False

    def _expire_queued_in_transaction(
        self, connection: sqlite3.Connection, timestamp: float
    ) -> int:
        return connection.execute(
            """
            UPDATE jobs
            SET status = 'cancelled', client_digest = NULL, cancel_requested = 1,
                updated_at = ?, finished_at = ?,
                error_code = 'job_expired',
                error_message = 'The queued job expired before execution.'
            WHERE status = 'queued' AND expires_at IS NOT NULL
                AND expires_at <= ?
            """,
            (timestamp, timestamp, timestamp),
        ).rowcount

    def _recover_stale_in_transaction(
        self, connection: sqlite3.Connection, timestamp: float
    ) -> int:
        cancelled = connection.execute(
            """
            UPDATE jobs
            SET status = 'cancelled', client_digest = NULL, worker_id = NULL, lease_expires_at = NULL,
                updated_at = ?, finished_at = ?, expires_at = ?
            WHERE status = 'running' AND lease_expires_at <= ?
                AND cancel_requested = 1
            """,
            (
                timestamp,
                timestamp,
                timestamp + self.settings.retention_seconds,
                timestamp,
            ),
        ).rowcount
        timed_out = connection.execute(
            """
            UPDATE jobs
            SET status = 'failed', client_digest = NULL, worker_id = NULL, lease_expires_at = NULL,
                updated_at = ?, finished_at = ?, expires_at = ?,
                error_code = 'job_timeout',
                error_message = 'The scientific workflow exhausted its total execution budget.'
            WHERE status = 'running' AND lease_expires_at <= ?
                AND cancel_requested = 0 AND started_at IS NOT NULL
                AND (started_at > ? OR started_at + ? <= ?)
            """,
            (
                timestamp,
                timestamp,
                timestamp + self.settings.retention_seconds,
                timestamp,
                timestamp + 5,
                self.settings.job_timeout_seconds,
                timestamp,
            ),
        ).rowcount
        exhausted = connection.execute(
            """
            UPDATE jobs
            SET status = 'failed', client_digest = NULL, worker_id = NULL, lease_expires_at = NULL,
                updated_at = ?, finished_at = ?, expires_at = ?,
                error_code = 'worker_restarts_exhausted',
                error_message = 'The backend restarted repeatedly while this job was running.'
            WHERE status = 'running' AND lease_expires_at <= ?
                AND cancel_requested = 0 AND attempt >= ?
            """,
            (
                timestamp,
                timestamp,
                timestamp + self.settings.retention_seconds,
                timestamp,
                self.settings.max_job_attempts,
            ),
        ).rowcount
        queued = connection.execute(
            """
            UPDATE jobs
            SET status = 'queued', worker_id = NULL, lease_expires_at = NULL,
                updated_at = ?
            WHERE status = 'running' AND lease_expires_at <= ?
                AND cancel_requested = 0 AND attempt < ?
            """,
            (timestamp, timestamp, self.settings.max_job_attempts),
        ).rowcount
        return cancelled + timed_out + exhausted + queued

    def recover_stale(self, now: float | None = None) -> int:
        timestamp = time.time() if now is None else now
        with self.transaction(immediate=True) as connection:
            return self._recover_stale_in_transaction(connection, timestamp)

    def claim_next(self, worker_id: str, now: float | None = None) -> ClaimedJob | None:
        timestamp = time.time() if now is None else now
        if not has_capacity(
            self.settings.data_root,
            required_bytes=self.settings.worker_admission_free_bytes,
            required_inodes=self.settings.worker_admission_free_inodes,
        ):
            return None
        with self.transaction(immediate=True) as connection:
            self._recover_stale_in_transaction(connection, timestamp)
            self._expire_queued_in_transaction(connection, timestamp)
            connection.execute(
                """
                UPDATE jobs
                SET status = 'failed', client_digest = NULL, updated_at = ?, finished_at = ?, expires_at = ?,
                    error_code = 'job_timeout',
                    error_message = 'The scientific workflow exhausted its total execution budget.'
                WHERE status = 'queued' AND started_at IS NOT NULL
                    AND (started_at > ? OR started_at + ? <= ?)
                """,
                (
                    timestamp,
                    timestamp,
                    timestamp + self.settings.retention_seconds,
                    timestamp + 5,
                    self.settings.job_timeout_seconds,
                    timestamp,
                ),
            )
            connection.execute(
                """
                UPDATE jobs
                SET status = 'failed', client_digest = NULL, updated_at = ?, finished_at = ?, expires_at = ?,
                    error_code = 'worker_restarts_exhausted',
                    error_message = 'The backend restarted repeatedly while this job was running.'
                WHERE status = 'queued' AND attempt >= ?
                """,
                (
                    timestamp,
                    timestamp,
                    timestamp + self.settings.retention_seconds,
                    self.settings.max_job_attempts,
                ),
            )
            row = connection.execute(
                "SELECT job_id, mode, attempt, category_policy, spacer_distance, "
                "bias_corrections_requested, bias_corrections_effective, "
                "option_warnings_json, started_at, tracr_model_type, leader_flank_length, molecule, spacer_database_sha256, viral_reference_sha256, viral_max_mismatches, repeat_reference_sha256, association_reference_sha256, association_grouping, association_groups_sha256, array_map_sha256, tool_options_json FROM jobs "
                "WHERE status = 'queued' AND attempt < ? AND expires_at > ? "
                "ORDER BY created_at, job_id LIMIT 1",
                (self.settings.max_job_attempts, timestamp),
            ).fetchone()
            if row is None:
                return None
            attempt = int(row["attempt"]) + 1
            execution_started_at = (
                timestamp if row["started_at"] is None else float(row["started_at"])
            )
            cursor = connection.execute(
                """
                UPDATE jobs
                SET status = 'running', attempt = ?, worker_id = ?,
                    started_at = ?, updated_at = ?, lease_expires_at = ?,
                    cancel_requested = 0, error_code = NULL, error_message = NULL
                WHERE job_id = ? AND status = 'queued'
                """,
                (
                    attempt,
                    worker_id,
                    execution_started_at,
                    timestamp,
                    timestamp + self.settings.worker_lease_seconds,
                    row["job_id"],
                ),
            )
            if cursor.rowcount != 1:
                return None
            return ClaimedJob(
                job_id=row["job_id"],
                mode=row["mode"],
                attempt=attempt,
                started_at=execution_started_at,
                category_policy=row["category_policy"],
                spacer_distance=int(row["spacer_distance"]),
                bias_corrections_requested=bool(row["bias_corrections_requested"]),
                bias_corrections_effective=bool(row["bias_corrections_effective"]),
                option_warnings=tuple(json.loads(row["option_warnings_json"])),
                tracr_model_type=row["tracr_model_type"],
                leader_flank_length=row["leader_flank_length"],
                molecule=row["molecule"],
                spacer_database_sha256=row["spacer_database_sha256"],
                association_reference_sha256=row["association_reference_sha256"],
                association_grouping=row["association_grouping"],
                association_groups_sha256=row["association_groups_sha256"],
                array_map_sha256=row["array_map_sha256"],
                repeat_reference_sha256=row["repeat_reference_sha256"],
                viral_reference_sha256=row["viral_reference_sha256"],
                viral_max_mismatches=row["viral_max_mismatches"],
                tool_options=validate_tool_options(json.loads(row["tool_options_json"]), row["mode"]),
            )

    def release_job(self, job_id: str, worker_id: str) -> None:
        """Release an owned job while preserving cancellation and retry ceilings."""

        now = time.time()
        with self.transaction(immediate=True) as connection:
            row = connection.execute(
                "SELECT attempt, cancel_requested, started_at FROM jobs "
                "WHERE job_id = ? AND status = 'running' AND worker_id = ?",
                (job_id, worker_id),
            ).fetchone()
            if row is None:
                raise LeaseLostError(job_id)
            started_at = row["started_at"]
            timed_out = (
                started_at is None
                or float(started_at) <= 0
                or float(started_at) > now + 5
                or float(started_at) + self.settings.job_timeout_seconds <= now
            )
            if bool(row["cancel_requested"]):
                connection.execute(
                    """
                    UPDATE jobs
                    SET status = 'cancelled', client_digest = NULL, worker_id = NULL,
                        lease_expires_at = NULL, updated_at = ?, finished_at = ?,
                        expires_at = ?, error_code = NULL, error_message = NULL
                    WHERE job_id = ? AND status = 'running' AND worker_id = ?
                    """,
                    (
                        now,
                        now,
                        now + self.settings.retention_seconds,
                        job_id,
                        worker_id,
                    ),
                )
            elif timed_out:
                connection.execute(
                    """
                    UPDATE jobs
                    SET status = 'failed', client_digest = NULL, worker_id = NULL, lease_expires_at = NULL,
                        updated_at = ?, finished_at = ?, expires_at = ?,
                        error_code = 'job_timeout',
                        error_message = 'The scientific workflow exhausted its total execution budget.'
                    WHERE job_id = ? AND status = 'running' AND worker_id = ?
                    """,
                    (
                        now,
                        now,
                        now + self.settings.retention_seconds,
                        job_id,
                        worker_id,
                    ),
                )
            elif int(row["attempt"]) >= self.settings.max_job_attempts:
                connection.execute(
                    """
                    UPDATE jobs
                    SET status = 'failed', client_digest = NULL, worker_id = NULL, lease_expires_at = NULL,
                        updated_at = ?, finished_at = ?, expires_at = ?,
                        error_code = 'worker_restarts_exhausted',
                        error_message = 'The backend restarted repeatedly while this job was running.'
                    WHERE job_id = ? AND status = 'running' AND worker_id = ?
                    """,
                    (
                        now,
                        now,
                        now + self.settings.retention_seconds,
                        job_id,
                        worker_id,
                    ),
                )
            else:
                connection.execute(
                    """
                    UPDATE jobs
                    SET status = 'queued', worker_id = NULL, lease_expires_at = NULL,
                        updated_at = ?
                    WHERE job_id = ? AND status = 'running' AND worker_id = ?
                    """,
                    (now, job_id, worker_id),
                )

    def heartbeat_job(self, job_id: str, worker_id: str) -> None:
        now = time.time()
        with self.connect() as connection:
            cursor = connection.execute(
                "UPDATE jobs SET lease_expires_at = ?, updated_at = ? "
                "WHERE job_id = ? AND status = 'running' AND worker_id = ?",
                (
                    now + self.settings.worker_lease_seconds,
                    now,
                    job_id,
                    worker_id,
                ),
            )
        if cursor.rowcount != 1:
            raise LeaseLostError(job_id)

    def cancellation_requested(self, job_id: str, worker_id: str) -> bool:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT cancel_requested FROM jobs "
                "WHERE job_id = ? AND status = 'running' AND worker_id = ?",
                (job_id, worker_id),
            ).fetchone()
        if row is None:
            raise LeaseLostError(job_id)
        return bool(row["cancel_requested"])

    def finish_job(
        self,
        *,
        job_id: str,
        worker_id: str,
        status: str,
        summary: Mapping[str, Any] | None,
        artifacts: Sequence[Mapping[str, Any]],
        error_code: str | None = None,
        error_message: str | None = None,
    ) -> None:
        if status not in TERMINAL_STATUSES:
            raise ValueError(f"invalid terminal status: {status}")
        now = time.time()
        with self.transaction(immediate=True) as connection:
            row = connection.execute(
                "SELECT cancel_requested FROM jobs "
                "WHERE job_id = ? AND status = 'running' AND worker_id = ?",
                (job_id, worker_id),
            ).fetchone()
            if row is None:
                raise LeaseLostError(job_id)
            if bool(row["cancel_requested"]):
                status = "cancelled"
                summary = None
                artifacts = ()
                error_code = None
                error_message = None
            serialized = (
                json.dumps(summary, separators=(",", ":"), allow_nan=False)
                if summary is not None
                else None
            )
            cursor = connection.execute(
                """
                UPDATE jobs
                SET status = ?, client_digest = NULL, updated_at = ?, finished_at = ?, expires_at = ?,
                    lease_expires_at = NULL, worker_id = NULL, summary_json = ?,
                    error_code = ?, error_message = ?
                WHERE job_id = ? AND status = 'running' AND worker_id = ?
                """,
                (
                    status,
                    now,
                    now,
                    now + self.settings.retention_seconds,
                    serialized,
                    error_code,
                    error_message,
                    job_id,
                    worker_id,
                ),
            )
            if cursor.rowcount != 1:
                raise LeaseLostError(job_id)
            connection.execute("DELETE FROM artifacts WHERE job_id = ?", (job_id,))
            for artifact in artifacts:
                connection.execute(
                    """
                    INSERT INTO artifacts(
                        job_id, artifact_id, relative_path, name,
                        size_bytes, sha256, media_type
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        job_id,
                        artifact["artifact_id"],
                        artifact["relative_path"],
                        artifact["name"],
                        artifact["size_bytes"],
                        artifact["sha256"],
                        artifact["media_type"],
                    ),
                )

    def heartbeat_worker(self, worker_id: str, state: str = "running") -> None:
        with self.connect() as connection:
            connection.execute(
                """
                INSERT INTO service_state(name, value, updated_at)
                VALUES(?, ?, ?)
                ON CONFLICT(name) DO UPDATE SET
                    value = excluded.value, updated_at = excluded.updated_at
                """,
                (f"worker:{worker_id}", state, time.time()),
            )

    def worker_health(self, now: float | None = None) -> str:
        timestamp = time.time() if now is None else now
        with self.connect() as connection:
            row = connection.execute(
                "SELECT value, updated_at FROM service_state "
                "WHERE name LIKE 'worker:%' ORDER BY updated_at DESC LIMIT 1"
            ).fetchone()
        if row is None:
            return "not_started"
        if timestamp - float(row["updated_at"]) > self.settings.worker_stale_seconds:
            return "stale"
        if row["value"] == "running":
            return "ok"
        if row["value"] == "storage_low":
            return "storage_low"
        return "stale"

    def database_ok(self) -> bool:
        connection: sqlite3.Connection | None = None
        try:
            connection = self.connect()
            connection.execute("PRAGMA busy_timeout = 1000")
            connection.execute("BEGIN IMMEDIATE")
            rows = {
                row["name"]: row["value"]
                for row in connection.execute(
                    "SELECT name, value FROM metadata "
                    "WHERE name IN ('schema_version', 'rate_salt')"
                )
            }
            schema_ok = int(rows.get("schema_version", "")) == _SCHEMA_VERSION
            columns = {
                row["name"] for row in connection.execute("PRAGMA table_info(jobs)")
            }
            schema_ok = schema_ok and "client_digest" in columns
            salt = rows.get("rate_salt")
            salt_ok = (
                isinstance(salt, str)
                and len(salt) == 64
                and len(bytes.fromhex(salt)) == 32
            )
            connection.execute("SELECT job_id FROM jobs LIMIT 1").fetchone()
            connection.execute(
                "INSERT INTO service_state(name, value, updated_at) VALUES(?, ?, ?) "
                "ON CONFLICT(name) DO UPDATE SET value = excluded.value, "
                "updated_at = excluded.updated_at",
                ("__readiness_probe__", "rollback", time.time()),
            )
            connection.rollback()
            return schema_ok and salt_ok
        except (OSError, TypeError, ValueError, sqlite3.Error):
            if connection is not None and connection.in_transaction:
                connection.rollback()
            return False
        finally:
            if connection is not None:
                connection.close()

    def expired_job_ids(self, now: float | None = None) -> list[str]:
        timestamp = time.time() if now is None else now
        with self.connect() as connection:
            return [
                row["job_id"]
                for row in connection.execute(
                    "SELECT job_id FROM jobs WHERE expires_at IS NOT NULL "
                    "AND expires_at <= ? AND status IN "
                    "('queued', 'completed', 'completed_no_eligible_groups', 'failed', 'cancelled')",
                    (timestamp,),
                ).fetchall()
            ]

    def freeze_expired_job(self, job_id: str, now: float | None = None) -> bool:
        """Atomically prevent an expired queued job from being claimed."""

        timestamp = time.time() if now is None else now
        with self.transaction(immediate=True) as connection:
            connection.execute(
                """
                UPDATE jobs
                SET status = 'cancelled', client_digest = NULL, cancel_requested = 1,
                    updated_at = ?, finished_at = ?
                WHERE job_id = ? AND status = 'queued'
                    AND expires_at IS NOT NULL AND expires_at <= ?
                """,
                (timestamp, timestamp, job_id, timestamp),
            )
            row = connection.execute(
                "SELECT status, expires_at FROM jobs WHERE job_id = ?",
                (job_id,),
            ).fetchone()
        return bool(
            row is not None
            and row["status"] in TERMINAL_STATUSES
            and row["expires_at"] is not None
            and float(row["expires_at"]) <= timestamp
        )

    def delete_expired_job(self, job_id: str, now: float | None = None) -> bool:
        timestamp = time.time() if now is None else now
        with self.connect() as connection:
            cursor = connection.execute(
                "DELETE FROM jobs WHERE job_id = ? AND expires_at <= ? AND status IN "
                "('completed', 'completed_no_eligible_groups', 'failed', 'cancelled')",
                (job_id, timestamp),
            )
            return cursor.rowcount == 1

    def artifact(self, job_id: str, artifact_id: str) -> dict[str, Any] | None:
        with self.connect() as connection:
            row = connection.execute(
                "SELECT artifact_id, relative_path, name, size_bytes, sha256, media_type "
                "FROM artifacts WHERE job_id = ? AND artifact_id = ?",
                (job_id, artifact_id),
            ).fetchone()
            return None if row is None else dict(row)
