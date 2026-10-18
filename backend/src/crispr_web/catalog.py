"""Bounded reads from an operator-imported, immutable public catalogue.

This module never opens the job database, runs a scientific tool, or downloads
data. A catalogue is replaced atomically by the offline importer.
"""

from __future__ import annotations

import base64
from collections import OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import stat
import threading
import time
from typing import Iterator
import zlib

from pydantic import ValidationError

from .catalog_models import (
    CatalogArray,
    CatalogAssembly,
    CatalogCasGene,
    CatalogCasSystem,
    CatalogModel,
    CatalogRepeat,
    CatalogSpacer,
    CatalogSummary,
)
from .catalog_search import CatalogSearchIndex, SearchIndexUnavailable
from .catalog_similarity import sequence_matcher
from .catalog_terms import gene_matches


class CatalogUnavailable(Exception):
    """The configured catalogue cannot safely answer a request."""


class CatalogCursorError(ValueError):
    """A cursor is invalid or belongs to different data or filters."""


class CatalogQueryError(ValueError):
    """A search needs a more specific public filter."""


@dataclass(frozen=True)
class Entity:
    table: str
    primary_key: str
    model: type[CatalogModel]


ENTITIES = {
    "assemblies": Entity("assemblies", "assembly_pk", CatalogAssembly),
    "arrays": Entity("arrays", "array_pk", CatalogArray),
    "repeats": Entity("repeats", "repeat_pk", CatalogRepeat),
    "spacers": Entity("spacers", "spacer_pk", CatalogSpacer),
    "cas-systems": Entity("cas_systems", "cas_system_pk", CatalogCasSystem),
    "cas-genes": Entity("cas_genes", "cas_gene_pk", CatalogCasGene),
}
_DNA = re.compile(r"[ACGTRYSWKMBDHVN]+", re.IGNORECASE)
_MAX_RESPONSE_BYTES = 2_000_000
_SHARD_ID_BASE = 1_000_000_000
_MAX_SHARD_BYTES = 256 * 1024**2
_MAX_COMPRESSED_BYTES = 64 * 1024**2
_MAX_ACTIVE_SHARD_BYTES = 544 * 1024**2
_MAX_CACHED_SHARD_BYTES = 128 * 1024**2
_ASSEMBLY = re.compile(r"GC[AF]_\d+\.\d+")


def _configure(connection: sqlite3.Connection) -> None:
    connection.row_factory = sqlite3.Row
    connection.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 262_144)
    connection.execute("PRAGMA query_only = ON")
    connection.execute("PRAGMA trusted_schema = OFF")
    connection.execute("PRAGMA cache_size = -1024")
    connection.create_function("catalog_gene_match", 3, gene_matches, deterministic=True)
    deadline = time.monotonic() + 0.5
    operations = 0

    def progress() -> int:
        nonlocal operations
        operations += 1000
        return int(operations > 2_000_000 or time.monotonic() > deadline)

    connection.set_progress_handler(progress, 1000)


def _bounded_rows(cursor: sqlite3.Cursor, budget: int = _MAX_RESPONSE_BYTES) -> tuple[list[sqlite3.Row], int]:
    rows = []
    used = 0
    for row in cursor:
        # Account incrementally, before an oversized result is materialized.
        # Allow room for public boolean normalization and JSON separators.
        used += len(json.dumps(dict(row)).encode()) + 128
        if used > budget:
            raise CatalogUnavailable("catalogue page exceeds the response limit")
        rows.append(row)
    return rows, used


def _query_parts(entity: Entity) -> tuple[str, str]:
    columns = []
    for field in entity.model.model_fields:
        if field == "id":
            source = f"r.{entity.primary_key}"
        elif field == "array_id":
            source = "r.array_pk"
        elif field == "cas_system_id":
            source = "r.cas_system_pk"
        elif field == "accession" and entity.table != "assemblies":
            source = "a.accession"
        elif field == "sequence_id":
            source = "s.sequence_id"
        else:
            source = f"r.{field}"
        columns.append(f"{source} AS {field}")
    source = f"{entity.table} AS r"
    if entity.table in {"repeats", "spacers"}:
        source += " JOIN arrays AS ar ON ar.array_pk = r.array_pk"
        source += " JOIN sequences AS s ON s.sequence_pk = ar.sequence_pk"
    elif entity.table != "assemblies":
        source += " JOIN sequences AS s ON s.sequence_pk = r.sequence_pk"
    if entity.table != "assemblies":
        source += " JOIN assemblies AS a ON a.assembly_pk = s.assembly_pk"
    return ", ".join(columns), source


class Catalog:
    def __init__(self, path: Path | None, search_index_path: Path | None = None):
        self.path = path
        self.search_index = CatalogSearchIndex(path, search_index_path)
        self._shard_cache: OrderedDict[tuple[str, str], bytearray] = OrderedDict()
        self._shard_lock = threading.Lock()
        self._shard_slots = threading.BoundedSemaphore(2)
        self._shard_budget = threading.Condition()
        self._active_shard_bytes = 0

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        if self.path is None:
            raise CatalogUnavailable("catalogue is not configured")
        connection = None
        try:
            # Do not resolve symlinks: an operator path must itself be safe.
            for part in (*reversed(self.path.parents), self.path):
                if stat.S_ISLNK(part.lstat().st_mode):
                    raise CatalogUnavailable("catalogue path is a symlink")
            if not stat.S_ISREG(self.path.lstat().st_mode):
                raise CatalogUnavailable("catalogue is not a regular file")
            connection = sqlite3.connect(
                self.path.as_uri() + "?mode=ro&immutable=1",
                uri=True,
                timeout=0.2,
            )
            _configure(connection)
            yield connection
        except (OSError, sqlite3.Error, ValidationError, TypeError, KeyError, EOFError, zlib.error) as error:
            raise CatalogUnavailable("catalogue is unavailable") from error
        finally:
            if connection is not None:
                connection.close()

    @contextmanager
    def _admit_shard(self, size: int) -> Iterator[None]:
        if type(size) is not int or not 0 < size <= _MAX_SHARD_BYTES:
            raise CatalogUnavailable("catalogue shard expansion is too large")
        # Each active reader retains an exact-size input image and SQLite's
        # copied image. The allowance also covers bounded decompression buffers
        # and SQLite's small pager cache. Admission precedes either allocation.
        weight = 2 * size + 8 * 1024**2
        deadline = time.monotonic() + 2
        if not self._shard_slots.acquire(timeout=2):
            raise CatalogUnavailable("catalogue read capacity is busy")
        reserved = False
        try:
            with self._shard_budget:
                while self._active_shard_bytes + weight > _MAX_ACTIVE_SHARD_BYTES:
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        raise CatalogUnavailable("catalogue memory capacity is busy")
                    self._shard_budget.wait(remaining)
                self._active_shard_bytes += weight
                reserved = True
            yield
        finally:
            if reserved:
                with self._shard_budget:
                    self._active_shard_bytes -= weight
                    self._shard_budget.notify_all()
            self._shard_slots.release()

    @contextmanager
    def _shard(self, row: sqlite3.Row) -> Iterator[sqlite3.Connection]:
        with self._admit_shard(row["uncompressed_bytes"]):
            with self._shard_connection(row) as connection:
                yield connection

    @contextmanager
    def _shard_connection(self, row: sqlite3.Row) -> Iterator[sqlite3.Connection]:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise CatalogUnavailable("catalogue shard path is invalid")
        if self.path is None:
            raise CatalogUnavailable("catalogue is not configured")
        path = self.path.parent / relative
        identity = (str(path), row["sha256"])
        with self._shard_lock:
            payload = self._shard_cache.get(identity)
            if payload is None:
                for part in (*reversed(path.parents), path):
                    if stat.S_ISLNK(part.lstat().st_mode):
                        raise CatalogUnavailable("catalogue shard is a symlink")
                metadata = path.stat()
                if not stat.S_ISREG(metadata.st_mode) or metadata.st_size > _MAX_COMPRESSED_BYTES:
                    raise CatalogUnavailable("catalogue shard size is invalid")
                with path.open("rb") as handle:
                    digest = hashlib.sha256()
                    total = 0
                    for block in iter(lambda: handle.read(1024**2), b""):
                        total += len(block)
                        if total > _MAX_COMPRESSED_BYTES:
                            raise CatalogUnavailable("catalogue compressed shard is too large")
                        digest.update(block)
                    if digest.hexdigest() != row["sha256"]:
                        raise CatalogUnavailable("catalogue shard checksum is invalid")
                    handle.seek(0)
                    payload = bytearray(row["uncompressed_bytes"])
                    view = memoryview(payload)
                    with gzip.GzipFile(fileobj=handle) as archive:
                        offset = 0
                        while offset < len(payload):
                            read = archive.readinto(view[offset:offset + 1024**2])
                            if not read:
                                raise CatalogUnavailable("catalogue shard expansion is truncated")
                            offset += read
                        if archive.read(1):
                            raise CatalogUnavailable("catalogue shard expansion is oversized")
                    view.release()
                if len(payload) <= _MAX_CACHED_SHARD_BYTES:
                    self._shard_cache[identity] = payload
                    while len(self._shard_cache) > 2 or sum(map(len, self._shard_cache.values())) > _MAX_CACHED_SHARD_BYTES:
                        self._shard_cache.popitem(last=False)
            if identity in self._shard_cache:
                self._shard_cache.move_to_end(identity)
        connection = sqlite3.connect(":memory:")
        try:
            connection.deserialize(payload)
            _configure(connection)
            yield connection
        finally:
            connection.close()

    @staticmethod
    def _summary(connection: sqlite3.Connection) -> CatalogSummary:
        row = connection.execute(
            "SELECT value FROM export_metadata WHERE key = 'summary'"
        ).fetchone()
        if row is None or not isinstance(row[0], str) or len(row[0]) > 131_072:
            raise CatalogUnavailable("catalogue metadata is invalid")
        summary = CatalogSummary.model_validate_json(row[0])
        if not summary.available or not summary.snapshot_id or summary.counts is None:
            raise CatalogUnavailable("catalogue metadata is incomplete")
        return summary

    def summary(self) -> CatalogSummary:
        if self.path is None:
            return CatalogSummary(
                available=False,
                limitations=["A reference catalogue has not been configured on this server."],
            )
        with self._connect() as connection:
            summary = self._summary(connection)
            if self.search_index.available(
                snapshot_id=summary.snapshot_id,
                shards_sha256=summary.provenance.get("shards_sha256", ""),
            ):
                previous = "Exact repeat, spacer and Cas profile searches require a selected genome. No sequence alignment or similarity search is performed."
                replacement = (
                    "Repeat and spacer searches match exact sequences across the snapshot; "
                    "Cas searches match source annotations, profiles and explicit gene-name tokens, ignoring letter case. "
                    "Search in DB also offers full-length repeat/spacer candidates with up to one "
                    "literal nucleotide substitution on either strand. No gapped alignment or "
                    "inferred Cas-family expansion is performed."
                )
                summary = summary.model_copy(update={
                    "limitations": [replacement if item == previous else item for item in summary.limitations]
                })
            return summary

    @staticmethod
    def _cursor(
        value: str | None, snapshot: str, entity: str, filters: str
    ) -> int:
        if value is None:
            return 0
        try:
            data = json.loads(base64.b64decode(value, altchars=b"-_", validate=True))
        except (ValueError, UnicodeError) as error:
            raise CatalogCursorError("invalid catalogue cursor") from error
        if (
            not isinstance(data, list)
            or len(data) != 4
            or data[:3] != [snapshot, entity, filters]
            or type(data[3]) is not int
            or not 0 <= data[3] < 2**63 - 1
        ):
            raise CatalogCursorError("catalogue cursor does not match this query or snapshot")
        return data[3]

    def page(
        self,
        name: str,
        *,
        q: str | None = None,
        assembly_accession: str | None = None,
        array_id: int | None = None,
        cursor: str | None = None,
        limit: int = 25,
        match: str = "exact",
    ) -> dict[str, object]:
        entity = ENTITIES[name]
        if match not in {"exact", "similar"}:
            raise CatalogQueryError("unknown sequence match mode")
        similar = match == "similar"
        if similar:
            if name not in {"repeats", "spacers"} or not q or len(q) > 1000 or not _DNA.fullmatch(q):
                raise CatalogQueryError("similar search requires a repeat or spacer DNA sequence of 1–1000 bases")
            q = q.upper()
        matcher = sequence_matcher(q) if similar else None
        filters = hashlib.sha256(
            json.dumps([q, assembly_accession, array_id, match] if similar else [q, assembly_accession, array_id]).encode()
        ).hexdigest()[:16]
        with self._connect() as connection:
            summary = self._summary(connection)
            last = self._cursor(cursor, summary.snapshot_id, name, filters)
            sharded = name != "assemblies" and connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'shards'"
            ).fetchone() is not None
            candidate_ids = None
            if sharded and q:
                if _ASSEMBLY.fullmatch(q.upper()):
                    if assembly_accession and assembly_accession != q.upper():
                        raise CatalogQueryError("the accession query and genome filter differ")
                    assembly_accession, q = q.upper(), None
                elif not assembly_accession and array_id is None:
                    try:
                        candidate_ids = self.search_index.candidates(
                            name, q, max(1, last // _SHARD_ID_BASE),
                            snapshot_id=summary.snapshot_id,
                            shards_sha256=summary.provenance.get("shards_sha256", ""),
                            similar=similar,
                        )
                    except SearchIndexUnavailable as error:
                        raise CatalogUnavailable("global exact search is temporarily unavailable") from error
                    # Initial snapshot digest verification is outside the SQL
                    # instruction/deadline budget; all subsequent SQL is bounded.
                    _configure(connection)
            select, source = _query_parts(entity)
            where = [f"r.{entity.primary_key} > ?"]
            parameters: list[object] = [last]
            order = f"r.{entity.primary_key}"
            if name == "assemblies" and q:
                # Prefix matches are paginated in accession order so even a
                # broad prefix uses its index without sorting the whole corpus.
                where = [
                    "r.accession > COALESCE((SELECT accession FROM assemblies "
                    "WHERE assembly_pk = ?), '')"
                ]
                order = "r.accession"
            assembly_alias = "r" if name == "assemblies" else "a"
            if assembly_accession:
                where.append(f"{assembly_alias}.accession = ?")
                parameters.append(assembly_accession)
            if array_id is not None:
                if name not in {"repeats", "spacers"}:
                    raise CatalogQueryError("array_id is supported for repeats and spacers only")
                where.append("r.array_pk = ?")
                parameters.append(array_id)
            if q:
                if similar:
                    where.append("catalog_sequence_similar(r.sequence)")
                    if name == "spacers":
                        where.append("r.kind = 'observed'")
                elif name == "assemblies" or q.upper().startswith(("GCA_", "GCF_")):
                    # Range predicates provide literal prefix search, including '_'.
                    prefix = q.upper()
                    where.extend([
                        f"{assembly_alias}.accession >= ?",
                        f"{assembly_alias}.accession < ?",
                    ])
                    parameters.extend([prefix, prefix + "\uffff"])
                elif name in {"arrays", "repeats", "spacers"} and _DNA.fullmatch(q):
                    column = "repeat_consensus" if name == "arrays" else "sequence"
                    where.append(f"upper(r.{column}) = ?")
                    parameters.append(q.upper())
                elif name == "cas-genes":
                    where.append("catalog_gene_match(r.annotation, r.profile, ?)")
                    parameters.append(q)
                elif name == "cas-systems":
                    where.append("r.subtype = ?")
                    parameters.append(q)
                elif name == "arrays":
                    where.append("r.call_id = ?")
                    parameters.append(q)
                else:
                    # Non-DNA occurrence queries cannot match a biological sequence.
                    where.append("0")
            parameters.append(limit + 1)
            query = (
                f"SELECT {select} FROM {source} WHERE {' AND '.join(where)} "
                f"ORDER BY {order} LIMIT ?"
            )
            continuation = None
            if not sharded:
                if matcher:
                    connection.create_function("catalog_sequence_similar", 1, matcher, deterministic=True)
                rows, _ = _bounded_rows(connection.execute(query, parameters))
            else:
                shard_where = ["id >= ?", "json_extract(counts_json, ?) > 0"]
                shard_parameters: list[object] = [
                    max(1, last // _SHARD_ID_BASE), "$." + entity.table
                ]
                if candidate_ids is not None:
                    if not candidate_ids:
                        return {"items": [], "next_cursor": None, "limit": limit}
                    shard_where.append("id IN (" + ",".join("?" for _ in candidate_ids) + ")")
                    shard_parameters.extend(candidate_ids)
                if assembly_accession:
                    shard_where.append(
                        "id IN (SELECT ash.shard_id FROM assembly_shards AS ash "
                        "JOIN assemblies AS a ON a.assembly_pk = ash.assembly_pk "
                        "WHERE a.accession = ?)"
                    )
                    shard_parameters.append(assembly_accession)
                if array_id is not None:
                    shard_where.append("id = ?")
                    shard_parameters.append(array_id // _SHARD_ID_BASE)
                shards = connection.execute(
                    "SELECT id, relative_path, sha256, uncompressed_bytes FROM shards WHERE "
                    + " AND ".join(shard_where) + " ORDER BY id LIMIT 5",
                    shard_parameters,
                ).fetchall()
                if candidate_ids is not None and {row["id"] for row in shards} != set(candidate_ids):
                    raise CatalogUnavailable("exact-search index references unavailable shards")
                rows = []
                used = 0
                for shard in shards[:4]:
                    with self._shard(shard) as shard_connection:
                        if matcher:
                            shard_connection.create_function("catalog_sequence_similar", 1, matcher, deterministic=True)
                        batch, batch_bytes = _bounded_rows(shard_connection.execute(
                            query, [*parameters[:-1], limit + 1 - len(rows)]
                        ), _MAX_RESPONSE_BYTES - used)
                        rows.extend(batch)
                        used += batch_bytes
                    if len(rows) > limit:
                        break
                if len(rows) <= limit and len(shards) > 4:
                    # An empty page can still advance over shards without matches.
                    continuation = (shards[3]["id"] + 1) * _SHARD_ID_BASE
            items = [entity.model.model_validate(dict(row)).model_dump() for row in rows[:limit]]
            if len(json.dumps(items).encode()) > _MAX_RESPONSE_BYTES:
                raise CatalogUnavailable("catalogue page exceeds the response limit")
            next_cursor = None
            if len(rows) > limit or continuation is not None:
                next_cursor = base64.urlsafe_b64encode(json.dumps([
                    summary.snapshot_id, name, filters,
                    continuation if continuation is not None else items[-1]["id"]
                ]).encode()).decode()
            return {"items": items, "next_cursor": next_cursor, "limit": limit}

    def detail(self, name: str, identifier: int | str) -> CatalogModel | None:
        entity = ENTITIES[name]
        select, source = _query_parts(entity)
        column = "accession" if name == "assemblies" else entity.primary_key
        with self._connect() as connection:
            self._summary(connection)
            query = f"SELECT {select} FROM {source} WHERE r.{column} = ? LIMIT 1"
            sharded = name != "assemblies" and connection.execute(
                "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'shards'"
            ).fetchone() is not None
            if sharded:
                shard = connection.execute(
                    "SELECT id, relative_path, sha256, uncompressed_bytes FROM shards WHERE id = ?",
                    (int(identifier) // _SHARD_ID_BASE,),
                ).fetchone()
                if shard is None:
                    return None
                with self._shard(shard) as shard_connection:
                    row = shard_connection.execute(query, (identifier,)).fetchone()
            else:
                row = connection.execute(query, (identifier,)).fetchone()
            return None if row is None else entity.model.model_validate(dict(row))
