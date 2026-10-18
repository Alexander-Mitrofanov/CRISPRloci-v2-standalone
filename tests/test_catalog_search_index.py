from __future__ import annotations

import hashlib
from pathlib import Path
import sqlite3
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from catalog_snapshot import create_tables, normalize_database
from import_catalog import SnapshotWriter
from index_catalog import BUILDING, FILENAME, build_index, query_hash


@pytest.fixture(autouse=True)
def fixture_disk_capacity(monkeypatch):
    monkeypatch.setattr("import_catalog.shutil.disk_usage", lambda _path: SimpleNamespace(free=10 * 1024**3))


def publish_fixture(directory, total=1):
    writer = SnapshotWriter(directory, "search-test")
    for shard_id in range(1, total + 1):
        db = sqlite3.connect(":memory:")
        create_tables(db)
        db.execute("INSERT INTO assemblies VALUES (1,?,'complete','valid',1,1,1,2)", (f"GCA_{shard_id:09d}.1",))
        db.execute("INSERT INTO sequences(sequence_pk,assembly_pk,sequence_id) VALUES (1,1,'contig-1')")
        db.execute("INSERT INTO arrays(array_pk,sequence_pk,call_id,repeat_consensus) VALUES (1,1,'array-1','ACGT')")
        db.execute("INSERT INTO repeats(repeat_pk,array_pk,ordinal,sequence) VALUES (1,1,1,'ACGT')")
        db.execute("INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,sequence) VALUES (1,1,1,'observed','GATTACA')")
        db.execute("INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,boundary_position) VALUES (2,1,2,'deletion',10)")
        db.execute("INSERT INTO cas_systems(cas_system_pk,sequence_pk,subtype) VALUES (1,1,'II-A')")
        db.execute("INSERT INTO cas_genes(cas_gene_pk,sequence_pk,annotation,profile) VALUES (1,1,'C25_Cas9_1','C25_Cas9_1')")
        db.execute("INSERT INTO cas_genes(cas_gene_pk,sequence_pk,annotation,profile) VALUES (2,1,'C25_Cas9b_1','C25_Cas9b_1')")
        db.commit()
        compressed, registry, counts, size = normalize_database(db.serialize(), shard_id)
        db.close()
        writer.accept({"id": shard_id, "sha256": hashlib.sha256(compressed).hexdigest(),
                       "uncompressed_bytes": size, "counts": counts, "assemblies": registry}, compressed)
    writer.publish(expected=total)


def test_global_index_preserves_source_files_and_exact_source_aliases(tmp_path):
    publish_fixture(tmp_path)
    immutable = [tmp_path / "catalog.sqlite3", tmp_path / "manifest.json", tmp_path / "shards/00000001.sqlite3.gz"]
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in immutable}
    result = build_index(tmp_path, workers=2, min_free=0)
    assert result["source_counts"]["spacers"] == 2
    assert result["candidate_counts"][3] == 1  # No fictitious sequence for a deletion.
    db = sqlite3.connect(tmp_path / FILENAME)
    try:
        metadata = dict(db.execute("SELECT key,value FROM metadata"))
        assert metadata["status"] == "complete"
        assert metadata["registry_sha256"] == before[tmp_path / "catalog.sqlite3"]
        for kind, value in [(1, "array-1"), (2, "acgt"), (3, "GATTACA"), (4, "II-A"),
                            (5, "Cas9"), (5, "c25_cas9_1"), (5, "Cas9b")]:
            assert db.execute("SELECT shard_id FROM hits WHERE kind=? AND key_hash=?", (kind, query_hash(kind, value))).fetchall() == [(1,)]
        assert db.execute("SELECT shard_id FROM hits WHERE kind=5 AND key_hash=?", (query_hash(5, "Cas"),)).fetchall() == []
    finally:
        db.close()
    assert {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in immutable} == before


def test_index_resumes_committed_shards_after_source_checksum_failure(tmp_path):
    publish_fixture(tmp_path, total=2)
    second = tmp_path / "shards/00000002.sqlite3.gz"
    original = second.read_bytes()
    second.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        build_index(tmp_path, workers=1, min_free=0, commit_batch=1)
    assert not (tmp_path / FILENAME).exists()
    db = sqlite3.connect(tmp_path / BUILDING)
    assert db.execute("SELECT shard_id FROM completed_shards").fetchall() == [(1,)]
    db.close()
    second.write_bytes(original)
    assert build_index(tmp_path, workers=1, min_free=0)["shards"] == 2


def test_index_cannot_publish_source_totals_mismatch(tmp_path):
    publish_fixture(tmp_path)
    db = sqlite3.connect(tmp_path / "catalog.sqlite3")
    db.execute("UPDATE export_metadata SET value=json_set(value,'$.counts.cas_genes',999) WHERE key='summary'")
    db.commit()
    db.close()
    with pytest.raises(ValueError, match="source totals disagree"):
        build_index(tmp_path, workers=1, min_free=0)
    assert not (tmp_path / FILENAME).exists()


def test_uncommitted_batch_has_neither_keys_nor_completion_receipts(tmp_path):
    publish_fixture(tmp_path, total=2)
    second = tmp_path / "shards/00000002.sqlite3.gz"
    original = second.read_bytes()
    second.write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        build_index(tmp_path, workers=1, min_free=0, commit_batch=2)
    db = sqlite3.connect(tmp_path / BUILDING)
    assert db.execute("SELECT count(*) FROM completed_shards").fetchone()[0] == 0
    assert db.execute("SELECT count(*) FROM hits").fetchone()[0] == 0
    db.close()
    second.write_bytes(original)
    assert build_index(tmp_path, workers=2, min_free=0, commit_batch=2)["shards"] == 2
    db = sqlite3.connect(tmp_path / FILENAME)
    assert db.execute("SELECT shard_id FROM hits WHERE kind=5 AND key_hash=? ORDER BY shard_id",
                      (query_hash(5, "Cas9"),)).fetchall() == [(1,), (2,)]
    db.close()
