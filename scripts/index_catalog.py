#!/usr/bin/env python3
"""Build an immutable exact-search sidecar without changing a published catalogue.

The index contains candidate shard IDs, not scientific results. Query hashes
are verified against exact original terms by the API, including source-derived
Cas label aliases. Four bounded workers read and verify shards; one writer
commits search keys and completion receipts together for safe resumption.
"""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import io
from itertools import chain
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import sqlite3
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "backend" / "src"))
from crispr_web.catalog_terms import TERM_VERSION, query_hash, terms_for

from catalog_snapshot import MAX_SHARD, TABLES


SCHEMA = "crisprloci.exact-search.v1"
FILENAME = "exact-search.sqlite3"
BUILDING = "exact-search.building.sqlite3"
MAX_COMPRESSED = 64 * 1024**2
MIN_FREE = 2 * 1024**3
QUERIES = {
    1: (
        "SELECT DISTINCT repeat_consensus FROM arrays WHERE repeat_consensus IS NOT NULL",
        "SELECT DISTINCT call_id FROM arrays WHERE call_id IS NOT NULL",
    ),
    2: ("SELECT DISTINCT sequence FROM repeats WHERE sequence IS NOT NULL",),
    3: ("SELECT DISTINCT sequence FROM spacers WHERE kind='observed' AND sequence IS NOT NULL",),
    4: ("SELECT DISTINCT subtype FROM cas_systems WHERE subtype IS NOT NULL",),
    5: (
        "SELECT DISTINCT annotation FROM cas_genes WHERE annotation IS NOT NULL",
        "SELECT DISTINCT profile FROM cas_genes WHERE profile IS NOT NULL",
    ),
}


def file_digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def safe_file(path, maximum=None):
    for candidate in (*path.parents, path):
        if candidate.is_symlink():
            raise ValueError(f"Catalogue source cannot contain symlinks: {path}")
    if not path.is_file() or (maximum is not None and path.stat().st_size > maximum):
        raise ValueError(f"Missing or oversized catalogue source: {path}")


def read_registry(directory):
    path = directory / "catalog.sqlite3"
    safe_file(path)
    connection = sqlite3.connect(path.as_uri() + "?mode=ro&immutable=1", uri=True)
    try:
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        row = connection.execute("SELECT value FROM export_metadata WHERE key='summary'").fetchone()
        if row is None:
            raise ValueError("Catalogue has no published summary")
        summary = json.loads(row[0])
        shards = connection.execute("SELECT id,relative_path,sha256,uncompressed_bytes,counts_json FROM shards ORDER BY id").fetchall()
        combined = hashlib.sha256()
        for shard_id, relative, digest, _, _ in shards:
            if relative != f"shards/{shard_id:08d}.sqlite3.gz":
                raise ValueError("Unexpected catalogue shard path")
            combined.update(f"{shard_id}:{digest}\n".encode())
        if combined.hexdigest() != summary["provenance"]["shards_sha256"]:
            raise ValueError("Catalogue shard inventory does not match its publication")
        return summary, shards, file_digest(path)
    finally:
        connection.close()


def scan_shard(directory, row):
    shard_id, relative, expected, expanded_size, raw_counts = row
    path = directory / relative
    safe_file(path, MAX_COMPRESSED)
    with path.open("rb") as handle:
        compressed = handle.read(MAX_COMPRESSED + 1)
    if len(compressed) > MAX_COMPRESSED or hashlib.sha256(compressed).hexdigest() != expected:
        raise ValueError(f"Catalogue shard checksum mismatch: {shard_id}")
    if not 0 < expanded_size <= MAX_SHARD:
        raise ValueError("Catalogue shard expansion is outside its bound")
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as handle:
        data = handle.read(MAX_SHARD + 1)
    if len(data) != expanded_size:
        raise ValueError(f"Catalogue shard expansion disagrees: {shard_id}")
    connection = sqlite3.connect(":memory:")
    try:
        connection.deserialize(data)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        connection.execute("PRAGMA temp_store=MEMORY")
        connection.setlimit(sqlite3.SQLITE_LIMIT_LENGTH, 1024**2)
        counts = {table: connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0] for table in TABLES}
        if counts != json.loads(raw_counts):
            raise ValueError(f"Catalogue source row counts disagree: {shard_id}")
        hits = []
        hit_counts = {}
        for kind, queries in QUERIES.items():
            hashes = set()
            for query in queries:
                for value, in connection.execute(query):
                    if not isinstance(value, str) or len(value) > 65536 or "\x00" in value:
                        raise ValueError("Invalid scientific search term")
                    hashes.update(query_hash(kind, term) for term in terms_for(kind, value))
            hit_counts[kind] = len(hashes)
            hits.extend((kind, digest, shard_id) for digest in sorted(hashes))
        return shard_id, expected, counts, hit_counts, hits
    finally:
        connection.close()


def build_index(directory, workers=4, min_free=MIN_FREE, commit_batch=50, cache_mib=512):
    if not 1 <= workers <= 4:
        raise ValueError("Use one to four bounded index workers")
    if not 1 <= commit_batch <= 1000:
        raise ValueError("Commit batches must contain one to 1,000 shards")
    if not 16 <= cache_mib <= 4096:
        raise ValueError("Offline SQLite cache must be between 16 and 4,096 MiB")
    directory = Path(directory).absolute()
    summary, shards, registry_digest = read_registry(directory)
    published = directory / FILENAME
    building = directory / BUILDING
    if published.exists():
        raise ValueError("Search index already published; choose a separate snapshot copy to rebuild")
    for path in (building, directory / ".exact-search.lock"):
        if path.is_symlink():
            raise ValueError("Search index paths cannot be symlinks")
    identity = {
        "schema": SCHEMA, "term_version": TERM_VERSION,
        "snapshot_id": summary["snapshot_id"], "registry_sha256": registry_digest,
        "shards_sha256": summary["provenance"]["shards_sha256"],
    }
    with (directory / ".exact-search.lock").open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        db = sqlite3.connect(building)
        try:
            db.execute("PRAGMA journal_mode=WAL")
            db.execute("PRAGMA synchronous=FULL")
            db.execute(f"PRAGMA cache_size=-{cache_mib * 1024}")
            db.execute("PRAGMA wal_autocheckpoint=32768")
            db.executescript("""
                CREATE TABLE IF NOT EXISTS metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL);
                CREATE TABLE IF NOT EXISTS hits(
                    kind INTEGER NOT NULL,key_hash BLOB NOT NULL,shard_id INTEGER NOT NULL,
                    PRIMARY KEY(kind,key_hash,shard_id)
                ) WITHOUT ROWID;
                CREATE TABLE IF NOT EXISTS completed_shards(
                    shard_id INTEGER PRIMARY KEY,sha256 TEXT NOT NULL,
                    row_counts_json TEXT NOT NULL,hit_counts_json TEXT NOT NULL
                );
            """)
            previous = dict(db.execute("SELECT key,value FROM metadata"))
            if previous and any(previous.get(key) != value for key, value in identity.items()):
                raise ValueError("Resume requires the same catalogue and search term semantics")
            with db:
                db.executemany("INSERT OR IGNORE INTO metadata VALUES (?,?)", identity.items())
                db.execute("INSERT OR REPLACE INTO metadata VALUES ('status','building')")
            completed = dict(db.execute("SELECT shard_id,sha256 FROM completed_shards"))
            source_hashes = {row[0]: row[2] for row in shards}
            if any(source_hashes.get(shard_id) != digest for shard_id, digest in completed.items()):
                raise ValueError("Completed search shards differ from the source inventory")
            pending_rows = iter(row for row in shards if row[0] not in completed)
            processed = len(completed)
            inserted = 0
            started = time.monotonic()
            last_reported = processed
            batch = []

            def commit_results():
                nonlocal processed, inserted, last_reported
                # A large sorted transaction may touch every existing page.
                # Reserve room for that WAL and conservatively estimated growth,
                # then truncate its checkpoint after each large transaction.
                working_headroom = max(
                    256 * 1024**2,
                    int(building.stat().st_size * 1.1)
                    + sum(len(item[4]) for item in batch) * 128 + 64 * 1024**2,
                )
                if shutil.disk_usage(directory).free < min_free + working_headroom:
                    raise ValueError("Search index paused at the free-space reserve; committed shards can resume")
                with db:
                    before = db.total_changes
                    # Sorting across the entire batch avoids repeatedly walking
                    # distant B-tree pages for each individual source shard.
                    ordered = sorted(chain.from_iterable(item[4] for item in batch))
                    db.executemany("INSERT INTO hits VALUES (?,?,?)", ordered)
                    if db.total_changes - before != sum(sum(item[3].values()) for item in batch):
                        raise ValueError("Search candidate counts disagree with inserted rows")
                    for shard_id, digest, counts, hit_counts, hits in batch:
                        db.execute("INSERT INTO completed_shards VALUES (?,?,?,?)", (
                            shard_id, digest, json.dumps(counts, sort_keys=True), json.dumps(hit_counts, sort_keys=True),
                        ))
                processed += len(batch)
                inserted += sum(len(item[4]) for item in batch)
                batch.clear()
                if processed - last_reported >= 250:
                    db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
                    print(json.dumps({"shards": processed, "total": len(shards), "new_hits": inserted,
                                      "bytes": building.stat().st_size, "seconds": round(time.monotonic() - started, 1)}), flush=True)
                    last_reported = processed

            with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as pool:
                pending = set()
                for _ in range(workers):
                    row = next(pending_rows, None)
                    if row is not None:
                        pending.add(pool.submit(scan_shard, directory, row))
                while pending:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        batch.append(future.result())
                        if len(batch) >= commit_batch:
                            commit_results()
                        row = next(pending_rows, None)
                        if row is not None:
                            pending.add(pool.submit(scan_shard, directory, row))
            if batch:
                commit_results()
            totals = Counter()
            hit_totals = Counter()
            completed_count = 0
            for counts, hits in db.execute("SELECT row_counts_json,hit_counts_json FROM completed_shards"):
                totals.update(json.loads(counts))
                hit_totals.update({int(kind): count for kind, count in json.loads(hits).items()})
                completed_count += 1
            if completed_count != len(shards):
                raise ValueError("Search index does not cover every source shard")
            if any(totals[table] != count for table, count in summary["counts"].items()):
                raise ValueError("Search index source totals disagree with the published catalogue")
            actual_hits = dict(db.execute("SELECT kind,count(*) FROM hits GROUP BY kind"))
            if {kind: actual_hits.get(kind, 0) for kind in QUERIES} != dict(hit_totals):
                raise ValueError("Search index candidate counts disagree with its completion receipts")
            if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
                raise ValueError("Search index integrity check failed")
            if file_digest(directory / "catalog.sqlite3") != registry_digest:
                raise ValueError("Published catalogue changed while building its search index")
            result = {**identity, "source_counts": dict(totals), "candidate_counts": dict(hit_totals),
                      "shards": completed_count, "created_at": datetime.now(timezone.utc).isoformat(), "status": "complete"}
            with db:
                db.execute("INSERT OR REPLACE INTO metadata VALUES ('status','complete')")
                for key in ("source_counts", "candidate_counts", "shards", "created_at"):
                    value = result[key]
                    db.execute("INSERT OR REPLACE INTO metadata VALUES (?,?)", (key, value if isinstance(value, str) else json.dumps(value, sort_keys=True)))
            db.execute("PRAGMA wal_checkpoint(TRUNCATE)")
            db.execute("PRAGMA journal_mode=DELETE")
        finally:
            db.close()
        with building.open("rb") as handle:
            os.fsync(handle.fileno())
        os.replace(building, published)
        descriptor = os.open(directory, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        result["bytes"] = published.stat().st_size
        return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--workers", type=int, choices=range(1, 5), default=4)
    parser.add_argument("--cache-mib", type=int, default=512,
                        help="Bounded offline SQLite cache (16–4096 MiB; does not change API limits)")
    parser.add_argument("--commit-batch", type=int, default=50,
                        help="Source shards per atomic index transaction (1–1000)")
    arguments = parser.parse_args()
    print(json.dumps(build_index(arguments.snapshot, arguments.workers,
                                 cache_mib=arguments.cache_mib, commit_batch=arguments.commit_batch), indent=2))
