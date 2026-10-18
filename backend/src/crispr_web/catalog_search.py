"""Bounded lookup in an offline exact-search index pinned to one snapshot."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
from pathlib import Path
import sqlite3
import stat
import threading
import time
from typing import Iterator

from .catalog_terms import KINDS, TERM_VERSION, query_hash
from .catalog_similarity import neighbours


class SearchIndexUnavailable(Exception):
    """Global exact search cannot safely use the configured sidecar."""


def _safe_file(path: Path) -> None:
    for item in (*reversed(path.parents), path):
        if stat.S_ISLNK(item.lstat().st_mode):
            raise SearchIndexUnavailable("search index path is a symlink")
    if not stat.S_ISREG(path.lstat().st_mode):
        raise SearchIndexUnavailable("search index is not a regular file")


class CatalogSearchIndex:
    def __init__(self, registry: Path | None, path: Path | None = None):
        self.registry = registry
        self.path = path or (registry.with_name("exact-search.sqlite3") if registry else None)
        self._digest_lock = threading.Lock()
        self._registry_identity: tuple[tuple[int, ...], str] | None = None

    def _registry_sha256(self) -> str:
        if self.registry is None:
            raise SearchIndexUnavailable("catalogue registry is not configured")
        _safe_file(self.registry)
        metadata = self.registry.stat()
        identity = (metadata.st_dev, metadata.st_ino, metadata.st_size, metadata.st_mtime_ns, metadata.st_ctime_ns)
        with self._digest_lock:
            if self._registry_identity is not None and self._registry_identity[0] == identity:
                return self._registry_identity[1]
            with self.registry.open("rb") as source:
                digest = hashlib.file_digest(source, "sha256").hexdigest()
            # The publication is immutable. An operator replacing it during the
            # first digest must retry against a consistently pinned registry.
            after = self.registry.stat()
            if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns):
                raise SearchIndexUnavailable("catalogue registry changed during verification")
            self._registry_identity = (identity, digest)
            return digest

    @contextmanager
    def _verified(self, *, snapshot_id: str, shards_sha256: str) -> Iterator[sqlite3.Connection]:
        if self.path is None:
            raise SearchIndexUnavailable("global exact search is not configured")
        connection = None
        try:
            _safe_file(self.path)
            registry_digest = self._registry_sha256()
            connection = sqlite3.connect(self.path.as_uri() + "?mode=ro&immutable=1", uri=True, timeout=0.2)
            connection.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 131_072)
            connection.execute("PRAGMA query_only=ON")
            connection.execute("PRAGMA trusted_schema=OFF")
            connection.execute("PRAGMA cache_size=-1024")
            deadline = time.monotonic() + 0.5
            operations = 0

            def progress():
                nonlocal operations
                operations += 1000
                return int(operations > 100_000 or time.monotonic() > deadline)

            connection.set_progress_handler(progress, 1000)
            tables = dict(connection.execute(
                "SELECT name,type FROM sqlite_master WHERE name IN ('metadata','hits')"
            ))
            if tables != {"metadata": "table", "hits": "table"}:
                raise SearchIndexUnavailable("exact-search index schema is invalid")
            connection.execute("SELECT kind,key_hash,shard_id FROM hits LIMIT 0")
            expected = {
                "schema": "crisprloci.exact-search.v1", "term_version": TERM_VERSION,
                "status": "complete", "snapshot_id": snapshot_id,
                "registry_sha256": registry_digest, "shards_sha256": shards_sha256,
            }
            actual = dict(connection.execute(
                "SELECT key,value FROM metadata WHERE key IN (?,?,?,?,?,?)",
                tuple(expected),
            ))
            if actual != expected or not shards_sha256:
                raise SearchIndexUnavailable("exact-search index does not match the published snapshot")
            yield connection
        except (OSError, sqlite3.Error, ValueError, TypeError) as error:
            raise SearchIndexUnavailable("global exact search is temporarily unavailable") from error
        finally:
            if connection is not None:
                connection.close()

    def available(self, *, snapshot_id: str, shards_sha256: str) -> bool:
        try:
            with self._verified(snapshot_id=snapshot_id, shards_sha256=shards_sha256):
                return True
        except SearchIndexUnavailable:
            return False

    def candidates(self, entity: str, query: str, first_shard: int, *, snapshot_id: str, shards_sha256: str,
                   similar: bool = False) -> list[int]:
        with self._verified(snapshot_id=snapshot_id, shards_sha256=shards_sha256) as connection:
            if similar:
                # Probe the existing immutable exact index, never scan its hits.
                # Each probe returns at most five shards. Their five smallest
                # distinct IDs are exactly the next candidate-shard page.
                deadline = time.monotonic() + 2
                operations = 0

                def progress():
                    nonlocal operations
                    operations += 1000
                    return int(operations > 2_000_000 or time.monotonic() > deadline)

                connection.set_progress_handler(progress, 1000)
                candidates: set[int] = set()
                for term in neighbours(query):
                    if time.monotonic() > deadline:
                        raise SearchIndexUnavailable("sequence search exceeded its read budget")
                    rows = connection.execute(
                        "SELECT shard_id FROM hits WHERE kind=? AND key_hash=? AND shard_id>=? "
                        "AND shard_id<=? ORDER BY shard_id LIMIT 5",
                        (KINDS[entity], query_hash(KINDS[entity], term), first_shard,
                         max(candidates) if len(candidates) == 5 else 999_999),
                    ).fetchall()
                    if any(type(row[0]) is not int or not 0 < row[0] < 1_000_000 for row in rows):
                        raise SearchIndexUnavailable("search index contains invalid shard references")
                    candidates = set(sorted(candidates | {row[0] for row in rows})[:5])
                return sorted(candidates)
            rows = connection.execute(
                "SELECT shard_id FROM hits WHERE kind=? AND key_hash=? AND shard_id>=? "
                "ORDER BY shard_id LIMIT 5",
                (KINDS[entity], query_hash(KINDS[entity], query), first_shard),
            ).fetchall()
            if any(type(row[0]) is not int or not 0 < row[0] < 1_000_000 for row in rows):
                raise SearchIndexUnavailable("exact-search index contains invalid shard references")
            return [row[0] for row in rows]
