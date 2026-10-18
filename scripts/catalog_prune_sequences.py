#!/usr/bin/env python3
"""Shrink an unpublished, stopped catalogue import without losing feature rows.

Only contigs unreferenced by arrays, Cas systems, and Cas genes are omitted.
Assembly sequence_count remains the original source count. Per-shard journals
make the checksum/receipt/registry replacement recoverable after interruption.
Never run this concurrently with import_catalog.py.
"""

from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import sqlite3

from catalog_snapshot import MAX_SHARD


RETENTION = "contigs-with-arrays-or-cas-features"
NORMALIZER_VERSION = "2-feature-context"
JOURNAL = ".sequence-pruning-journal.json"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fsync_directory(path):
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_bytes(path, data):
    temporary = path.with_name(path.name + ".pruning-tmp")
    if path.is_symlink() or temporary.is_symlink():
        raise ValueError("Refusing symlinked catalogue metadata")
    with temporary.open("wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)


def prune_shard(compressed):
    """Return new compressed bytes, retained and omitted contig counts, size."""
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
        source_count = db.execute("SELECT count(*) FROM sequences").fetchone()[0]
        db.execute("""DELETE FROM sequences
            WHERE NOT EXISTS (SELECT 1 FROM arrays WHERE arrays.sequence_pk=sequences.sequence_pk)
              AND NOT EXISTS (SELECT 1 FROM cas_systems WHERE cas_systems.sequence_pk=sequences.sequence_pk)
              AND NOT EXISTS (SELECT 1 FROM cas_genes WHERE cas_genes.sequence_pk=sequences.sequence_pk)
        """)
        kept = db.execute("SELECT count(*) FROM sequences").fetchone()[0]
        if db.execute("PRAGMA foreign_key_check").fetchone() is not None:
            raise ValueError("Sequence pruning broke a scientific relation")
        db.commit()
        db.execute("VACUUM")
        normalized = db.serialize()
        return gzip.compress(normalized, compresslevel=6, mtime=0), kept, source_count - kept, len(normalized)
    finally:
        db.close()


def finish_pending(directory, registry):
    """Finish a previously staged replacement, including after a process crash."""
    journal_path = directory / JOURNAL
    if not journal_path.exists():
        return
    if journal_path.is_symlink():
        raise ValueError("Invalid pruning journal")
    journal = json.loads(journal_path.read_text())
    shard_id = journal["receipt"]["id"]
    path = directory / "shards" / f"{shard_id:08d}.sqlite3.gz"
    temporary = path.with_name(path.name + ".pruned")
    if path.is_symlink() or temporary.is_symlink():
        raise ValueError("Invalid catalogue shard path")
    observed = sha(path.read_bytes())
    if observed == journal["previous_sha256"]:
        if not temporary.is_file() or sha(temporary.read_bytes()) != journal["receipt"]["sha256"]:
            raise ValueError("Pending compacted shard is missing or changed")
        os.replace(temporary, path)
        fsync_directory(path.parent)
    elif observed != journal["receipt"]["sha256"]:
        raise ValueError("Interrupted pruning found an unexpected shard checksum")
    receipt = journal["receipt"]
    receipt_path = directory / "receipts" / f"{shard_id:08d}.json"
    atomic_bytes(receipt_path, (json.dumps(receipt, sort_keys=True) + "\n").encode())
    with registry:
        cursor = registry.execute(
            "UPDATE shards SET sha256=?,uncompressed_bytes=?,counts_json=? WHERE id=? AND sha256 IN (?,?)",
            (receipt["sha256"], receipt["uncompressed_bytes"], json.dumps(receipt["counts"]),
             shard_id, journal["previous_sha256"], receipt["sha256"]),
        )
        if cursor.rowcount != 1:
            raise ValueError("Catalogue registry changed during sequence pruning")
    journal_path.unlink()
    fsync_directory(directory)


def optimize(directory, limit=None):
    directory = Path(directory).absolute()
    if directory.is_symlink() or (directory / "catalog.sqlite3").exists():
        raise ValueError("Only unpublished catalogue imports can be compacted")
    database = directory / "catalog.building.sqlite3"
    if not database.is_file() or database.is_symlink():
        raise ValueError("No unpublished catalogue registry found")
    lock_path = directory / ".sequence-pruning.lock"
    if lock_path.is_symlink():
        raise ValueError("Invalid pruning lock")
    with lock_path.open("a+b") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        registry = sqlite3.connect(database)
        try:
            finish_pending(directory, registry)
            rows = registry.execute("SELECT id,sha256,counts_json FROM shards ORDER BY id").fetchall()
            processed = removed = saved = 0
            for shard_id, expected, raw_counts in rows:
                receipt_path = directory / "receipts" / f"{shard_id:08d}.json"
                if receipt_path.is_symlink():
                    raise ValueError("Invalid shard receipt")
                receipt = json.loads(receipt_path.read_text())
                if receipt.get("sequence_retention") == RETENTION:
                    continue
                path = directory / "shards" / f"{shard_id:08d}.sqlite3.gz"
                if path.is_symlink() or not path.is_file() or path.stat().st_size > MAX_SHARD:
                    raise ValueError("Invalid stored shard")
                data = path.read_bytes()
                if sha(data) != expected or receipt["sha256"] != expected:
                    raise ValueError("Catalogue shard or receipt checksum mismatch")
                compressed, kept, omitted, size = prune_shard(data)
                counts = json.loads(raw_counts)
                if counts["sequences"] != kept + omitted:
                    raise ValueError("Catalogue sequence counts disagree")
                counts["sequences"] = kept
                receipt = {**receipt, "sha256": sha(compressed), "bytes": len(compressed),
                           "uncompressed_bytes": size, "counts": counts,
                           "sequence_retention": RETENTION, "omitted_feature_free_contigs": omitted}
                temporary = path.with_name(path.name + ".pruned")
                atomic_bytes(temporary, compressed)
                journal = {"previous_sha256": expected, "receipt": receipt}
                atomic_bytes(directory / JOURNAL, (json.dumps(journal, sort_keys=True) + "\n").encode())
                finish_pending(directory, registry)
                processed += 1
                removed += omitted
                saved += len(data) - len(compressed)
                if processed % 50 == 0:
                    print(json.dumps({"processed": processed, "omitted_contigs": removed, "saved_bytes": saved}), flush=True)
                if limit is not None and processed >= limit:
                    break
            if limit is None:
                with registry:
                    registry.execute("INSERT OR REPLACE INTO export_metadata VALUES ('normalizer_version',?)", (NORMALIZER_VERSION,))
            return {"processed": processed, "omitted_contigs": removed, "saved_bytes": saved}
        finally:
            registry.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--limit", type=int)
    arguments = parser.parse_args()
    if arguments.limit is not None and arguments.limit < 1:
        parser.error("--limit must be positive")
    print(json.dumps(optimize(arguments.snapshot, arguments.limit)))
