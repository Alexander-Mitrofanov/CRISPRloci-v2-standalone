#!/usr/bin/env python3
"""Compact auxiliary indexes in a stopped, unpublished catalogue import.

Every table row, scientific value, ID, and relationship is retained. At most
four shard images are processed concurrently; durable file, receipt and registry
updates are serialized through the existing recoverable compaction journal.
Never run this concurrently with import_catalog.py or the sequence pruner.
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import fcntl
import gzip
import io
import json
from pathlib import Path
import sqlite3

from catalog_prune_sequences import JOURNAL, atomic_bytes, finish_pending, sha
from catalog_snapshot import MAX_SHARD, TABLES


INDEX_POLICY = "genome-scoped-relational-v1"
REMOVABLE_INDEXES = frozenset({
    "idx_arrays_call_id", "idx_arrays_repeat_consensus", "idx_repeats_sequence",
    "idx_spacers_sequence", "idx_cas_systems_call_id", "idx_cas_systems_subtype",
    "idx_cas_genes_annotation", "idx_cas_genes_profile",
})
RETAINED_INDEXES = frozenset({
    "idx_sequences_assembly_pk", "idx_arrays_sequence_pk", "idx_repeats_array_pk_ordinal",
    "idx_spacers_array_pk_ordinal", "idx_cas_systems_sequence_pk",
    "idx_cas_genes_sequence_pk", "idx_cas_genes_cas_system_pk",
})


def compact_shard(compressed):
    """Return compressed bytes, unchanged table counts, byte size, removed names."""
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:
        content = stream.read(MAX_SHARD + 1)
    if len(content) > MAX_SHARD:
        raise ValueError("Catalogue shard exceeds the decompression bound")
    db = sqlite3.connect(":memory:")
    try:
        db.deserialize(content)
        db.execute("PRAGMA trusted_schema=OFF")
        db.execute("PRAGMA foreign_keys=ON")
        if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
            raise ValueError("Invalid source shard")
        indexes = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='index'")}
        if not RETAINED_INDEXES <= indexes:
            raise ValueError("Source shard is missing an essential relational index")
        counts = {table: db.execute(f"SELECT count(*) FROM {table}").fetchone()[0] for table in TABLES}
        removed = sorted(indexes & REMOVABLE_INDEXES)
        for name in removed:
            db.execute(f"DROP INDEX {name}")
        if db.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise ValueError("Source shard has broken scientific relationships")
        db.commit()
        db.execute("VACUUM")
        after = {table: db.execute(f"SELECT count(*) FROM {table}").fetchone()[0] for table in TABLES}
        if counts != after:
            raise ValueError("Index compaction changed a scientific table count")
        remaining = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='index'")}
        if remaining != indexes - set(removed):
            raise ValueError("Index compaction changed an unselected index")
        normalized = db.serialize()
        return gzip.compress(normalized, compresslevel=6, mtime=0), counts, len(normalized), removed
    finally:
        db.close()


def prepare_shard(directory, row):
    """Read and compute one candidate; workers never write files or the registry."""
    shard_id, expected, raw_counts = row
    receipt_path = directory / "receipts" / f"{shard_id:08d}.json"
    if receipt_path.is_symlink() or receipt_path.stat().st_size > 1024**2:
        raise ValueError("Invalid shard receipt")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("index_policy") == INDEX_POLICY:
        return None
    path = directory / "shards" / f"{shard_id:08d}.sqlite3.gz"
    if path.is_symlink() or not path.is_file() or path.stat().st_size > MAX_SHARD:
        raise ValueError("Invalid stored shard")
    data = path.read_bytes()
    if sha(data) != expected or receipt["sha256"] != expected or receipt.get("id") != shard_id:
        raise ValueError("Catalogue shard or receipt checksum mismatch")
    compressed, counts, size, removed = compact_shard(data)
    if counts != json.loads(raw_counts) or counts != receipt["counts"]:
        raise ValueError("Catalogue scientific counts disagree with the stored shard")
    result = {**receipt, "sha256": sha(compressed), "bytes": len(compressed),
              "uncompressed_bytes": size, "index_policy": INDEX_POLICY,
              "removed_auxiliary_indexes": removed}
    return path, expected, result, compressed, len(data) - len(compressed)


def optimize(directory, workers=4, limit=None):
    if not 1 <= workers <= 4:
        raise ValueError("Index compaction supports one to four workers")
    if limit is not None and limit < 1:
        raise ValueError("Compaction limit must be positive")
    directory = Path(directory).absolute()
    for path in (*directory.parents, directory):
        if path.is_symlink():
            raise ValueError("Catalogue directory cannot contain a symlink")
    if (directory / "catalog.sqlite3").exists():
        raise ValueError("Only unpublished catalogue imports can be compacted")
    database = directory / "catalog.building.sqlite3"
    if not database.is_file() or database.is_symlink():
        raise ValueError("No unpublished catalogue registry found")
    for name in ("shards", "receipts", ".sequence-pruning.lock"):
        if (directory / name).is_symlink():
            raise ValueError("Catalogue members cannot be symlinks")
    with (directory / ".sequence-pruning.lock").open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        registry = sqlite3.connect(database)
        try:
            finish_pending(directory, registry)
            rows = registry.execute("SELECT id,sha256,counts_json FROM shards ORDER BY id").fetchall()
            if limit is not None:
                rows = rows[:limit]
            iterator = iter(rows)
            processed = skipped = saved = 0
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending = set()
                for _ in range(workers):
                    row = next(iterator, None)
                    if row is not None:
                        pending.add(pool.submit(prepare_shard, directory, row))
                while pending:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for future in done:
                        result = future.result()
                        if result is None:
                            skipped += 1
                        else:
                            path, previous, receipt, compressed, reclaimed = result
                            atomic_bytes(path.with_name(path.name + ".pruned"), compressed)
                            journal = {"previous_sha256": previous, "receipt": receipt}
                            atomic_bytes(directory / JOURNAL, (json.dumps(journal, sort_keys=True) + "\n").encode())
                            finish_pending(directory, registry)
                            processed += 1
                            saved += reclaimed
                            if processed % 100 == 0:
                                print(json.dumps({"processed": processed, "skipped": skipped, "saved_bytes": saved}), flush=True)
                        row = next(iterator, None)
                        if row is not None:
                            pending.add(pool.submit(prepare_shard, directory, row))
            if limit is None:
                with registry:
                    registry.execute("INSERT OR REPLACE INTO export_metadata VALUES ('index_policy',?)", (INDEX_POLICY,))
            return {"processed": processed, "skipped": skipped, "saved_bytes": saved,
                    "index_policy": INDEX_POLICY}
        finally:
            registry.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--workers", type=int, choices=range(1, 5), default=4)
    parser.add_argument("--limit", type=int)
    arguments = parser.parse_args()
    print(json.dumps(optimize(arguments.snapshot, arguments.workers, arguments.limit)))
