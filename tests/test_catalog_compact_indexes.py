from __future__ import annotations

import gzip
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from catalog_compact_indexes import INDEX_POLICY, REMOVABLE_INDEXES, RETAINED_INDEXES, compact_shard, optimize
from catalog_prune_sequences import sha
from catalog_snapshot import TABLES, create_tables


def legacy_shard():
    db = sqlite3.connect(":memory:")
    create_tables(db)
    db.execute("INSERT INTO assemblies VALUES (1,'GCA_000000001.1','complete','valid',2,1,1,1)")
    db.executemany("INSERT INTO sequences VALUES (?,1,?,NULL,300,NULL,?)", [(1, "contig-1", 1), (2, "feature-free-retained", 2)])
    db.execute("INSERT INTO arrays(array_pk,sequence_pk,call_id,repeat_consensus,strand) VALUES (1,1,'array-1','ACTG',NULL)")
    db.execute("INSERT INTO repeats(repeat_pk,array_pk,ordinal,start,end,sequence) VALUES (1,1,1,5,8,'ACTG')")
    db.execute("INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,boundary_position) VALUES (1,1,1,'deletion',8)")
    db.execute("INSERT INTO cas_systems(cas_system_pk,sequence_pk,call_id,subtype) VALUES (1,1,'system-1','I-A')")
    db.execute("INSERT INTO cas_genes(cas_gene_pk,sequence_pk,cas_system_pk,call_id,annotation,profile,score) VALUES (1,1,1,'gene-1','Cas1','Cas1-profile',0)")
    # Deliberately independent of the current normalizer's INDEXES constant.
    legacy_indexes = {
        "sequences": [("assembly_pk",)],
        "arrays": [("sequence_pk",), ("call_id",), ("repeat_consensus",)],
        "repeats": [("array_pk", "ordinal"), ("sequence",)],
        "spacers": [("array_pk", "ordinal"), ("sequence",)],
        "cas_systems": [("sequence_pk",), ("call_id",), ("subtype",)],
        "cas_genes": [("sequence_pk",), ("cas_system_pk",), ("annotation",), ("profile",)],
    }
    for table, groups in legacy_indexes.items():
        for columns in groups:
            db.execute(f"CREATE INDEX idx_{table}_{'_'.join(columns)} ON {table} ({','.join(columns)})")
    db.commit()
    values = {table: db.execute(f"SELECT * FROM {table} ORDER BY 1").fetchall() for table in TABLES}
    data = db.serialize()
    db.close()
    return gzip.compress(data, mtime=0), len(data), values


def test_compaction_changes_indexes_only_and_preserves_all_rows():
    compressed, _, values = legacy_shard()
    result, counts, size, removed = compact_shard(compressed)
    assert set(removed) == REMOVABLE_INDEXES
    db = sqlite3.connect(":memory:")
    try:
        db.deserialize(gzip.decompress(result))
        assert len(gzip.decompress(result)) == size
        assert {table: db.execute(f"SELECT * FROM {table} ORDER BY 1").fetchall() for table in TABLES} == values
        assert counts == {table: len(rows) for table, rows in values.items()}
        indexes = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='index'")}
        assert RETAINED_INDEXES <= indexes
        assert not indexes & REMOVABLE_INDEXES
        assert "sqlite_autoindex_assemblies_1" in indexes
        assert db.execute("PRAGMA foreign_key_check").fetchone() is None
    finally:
        db.close()


def snapshot(directory, total=8):
    (directory / "shards").mkdir()
    (directory / "receipts").mkdir()
    data, size, values = legacy_shard()
    counts = {table: len(rows) for table, rows in values.items()}
    db = sqlite3.connect(directory / "catalog.building.sqlite3")
    db.execute("CREATE TABLE shards(id INTEGER PRIMARY KEY,sha256 TEXT,uncompressed_bytes INTEGER,counts_json TEXT)")
    db.execute("CREATE TABLE export_metadata(key TEXT PRIMARY KEY,value TEXT)")
    db.execute("INSERT INTO export_metadata VALUES ('normalizer_version','2-feature-context')")
    for shard_id in range(1, total + 1):
        (directory / f"shards/{shard_id:08d}.sqlite3.gz").write_bytes(data)
        receipt = {"id": shard_id, "sha256": sha(data), "bytes": len(data),
                   "uncompressed_bytes": size, "counts": counts, "source_sha256": "a" * 64}
        (directory / f"receipts/{shard_id:08d}.json").write_text(json.dumps(receipt))
        db.execute("INSERT INTO shards VALUES (?,?,?,?)", (shard_id, sha(data), size, json.dumps(counts)))
    db.commit()
    db.close()


def test_four_worker_compaction_updates_durable_receipts_and_resumes(tmp_path):
    snapshot(tmp_path)
    result = optimize(tmp_path, workers=4)
    assert result["processed"] == 8
    db = sqlite3.connect(tmp_path / "catalog.building.sqlite3")
    try:
        assert dict(db.execute("SELECT key,value FROM export_metadata")) == {
            "normalizer_version": "2-feature-context", "index_policy": INDEX_POLICY}
        for shard_id, expected, size, counts in db.execute("SELECT * FROM shards"):
            path = tmp_path / f"shards/{shard_id:08d}.sqlite3.gz"
            receipt = json.loads((tmp_path / f"receipts/{shard_id:08d}.json").read_text())
            assert expected == sha(path.read_bytes()) == receipt["sha256"]
            assert receipt["counts"] == json.loads(counts)
            assert receipt["uncompressed_bytes"] == size
            assert receipt["source_sha256"] == "a" * 64
    finally:
        db.close()
    rerun = optimize(tmp_path, workers=4)
    assert (rerun["processed"], rerun["skipped"]) == (0, 8)
    assert not (tmp_path / ".sequence-pruning-journal.json").exists()


def test_compactor_rejects_corrupt_and_published_snapshots(tmp_path):
    snapshot(tmp_path, total=1)
    (tmp_path / "shards/00000001.sqlite3.gz").write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        optimize(tmp_path)
    (tmp_path / "catalog.sqlite3").write_bytes(b"published")
    with pytest.raises(ValueError, match="Only unpublished"):
        optimize(tmp_path)
