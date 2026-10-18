from __future__ import annotations

import gzip
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from catalog_prune_sequences import JOURNAL, RETENTION, finish_pending, optimize, prune_shard, sha
from catalog_snapshot import create_tables


def unpruned_shard():
    db = sqlite3.connect(":memory:")
    create_tables(db)
    db.execute("INSERT INTO assemblies VALUES (1,'GCA_000000001.1','complete','valid',4,1,1,1)")
    for key in range(1, 5):
        db.execute("INSERT INTO sequences VALUES (?,?,?,NULL,300,?,?)", (key, 1, f"contig-{key}", "a" * 64, key))
    db.execute("INSERT INTO arrays(array_pk,sequence_pk,call_id) VALUES (1,1,'array-1')")
    db.execute("INSERT INTO cas_systems(cas_system_pk,sequence_pk,call_id) VALUES (1,2,'system-1')")
    db.execute("INSERT INTO cas_genes(cas_gene_pk,sequence_pk,call_id) VALUES (1,3,'gene-1')")
    db.commit()
    data = db.serialize()
    db.close()
    return gzip.compress(data, mtime=0), len(data)


def snapshot(tmp_path):
    data, size = unpruned_shard()
    (tmp_path / "shards").mkdir()
    (tmp_path / "receipts").mkdir()
    path = tmp_path / "shards/00000001.sqlite3.gz"
    path.write_bytes(data)
    receipt = {"id": 1, "sha256": sha(data), "bytes": len(data), "uncompressed_bytes": size,
               "counts": {"assemblies": 1, "sequences": 4, "arrays": 1, "cas_systems": 1, "cas_genes": 1}}
    (tmp_path / "receipts/00000001.json").write_text(json.dumps(receipt))
    db = sqlite3.connect(tmp_path / "catalog.building.sqlite3")
    db.execute("CREATE TABLE export_metadata(key TEXT PRIMARY KEY,value TEXT)")
    db.execute("CREATE TABLE shards(id INTEGER PRIMARY KEY,sha256 TEXT,uncompressed_bytes INTEGER,counts_json TEXT)")
    db.execute("INSERT INTO shards VALUES (?,?,?,?)", (1, receipt["sha256"], size, json.dumps(receipt["counts"])))
    db.commit()
    db.close()
    return path, receipt


def test_pruning_preserves_every_feature_and_original_sequence_count():
    compressed, _ = unpruned_shard()
    result, kept, omitted, size = prune_shard(compressed)
    assert (kept, omitted) == (3, 1)
    assert size == len(gzip.decompress(result))
    db = sqlite3.connect(":memory:")
    db.deserialize(gzip.decompress(result))
    assert db.execute("SELECT sequence_count FROM assemblies").fetchone()[0] == 4
    assert db.execute("SELECT sequence_pk FROM sequences ORDER BY sequence_pk").fetchall() == [(1,), (2,), (3,)]
    assert db.execute("PRAGMA foreign_key_check").fetchone() is None
    db.close()


def test_snapshot_compaction_updates_receipt_and_registry_and_is_resumable(tmp_path):
    path, _ = snapshot(tmp_path)
    result = optimize(tmp_path)
    assert result["omitted_contigs"] == 1
    receipt = json.loads((tmp_path / "receipts/00000001.json").read_text())
    assert receipt["sequence_retention"] == RETENTION
    assert receipt["sha256"] == sha(path.read_bytes())
    db = sqlite3.connect(tmp_path / "catalog.building.sqlite3")
    digest, counts = db.execute("SELECT sha256,counts_json FROM shards").fetchone()
    assert digest == receipt["sha256"]
    assert json.loads(counts)["sequences"] == 3
    assert db.execute("SELECT value FROM export_metadata WHERE key='normalizer_version'").fetchone()[0] == "2-feature-context"
    db.close()
    assert optimize(tmp_path)["processed"] == 0
    (tmp_path / "catalog.sqlite3").write_bytes(b"published")
    with pytest.raises(ValueError, match="Only unpublished"):
        optimize(tmp_path)


def test_interrupted_shard_replacement_recovers_before_registry_commit(tmp_path):
    path, receipt = snapshot(tmp_path)
    original = receipt["sha256"]
    data, kept, _, size = prune_shard(path.read_bytes())
    receipt.update(sha256=sha(data), bytes=len(data), uncompressed_bytes=size, sequence_retention=RETENTION)
    receipt["counts"]["sequences"] = kept
    (tmp_path / JOURNAL).write_text(json.dumps({"previous_sha256": original, "receipt": receipt}))
    # Simulate interruption after shard replacement but before receipt/SQL writes.
    path.write_bytes(data)
    db = sqlite3.connect(tmp_path / "catalog.building.sqlite3")
    finish_pending(tmp_path, db)
    assert db.execute("SELECT sha256 FROM shards").fetchone()[0] == receipt["sha256"]
    assert not (tmp_path / JOURNAL).exists()
    db.close()
