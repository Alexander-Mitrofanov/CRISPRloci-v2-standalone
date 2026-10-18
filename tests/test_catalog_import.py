from __future__ import annotations

import gzip
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import sys
import tarfile
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from catalog_snapshot import create_tables, normalize_database, read_frame, transform_bundle
from import_catalog import SnapshotWriter


@pytest.fixture(autouse=True)
def fixture_disk_capacity(monkeypatch):
    # Tiny unit fixtures must not depend on the developer volume's free space.
    monkeypatch.setattr("import_catalog.shutil.disk_usage", lambda _path: SimpleNamespace(free=10 * 1024**3))


def fixture_database():
    db = sqlite3.connect(":memory:")
    create_tables(db)
    db.execute("ALTER TABLE arrays ADD COLUMN native_json TEXT")
    db.execute("INSERT INTO assemblies VALUES (1,'GCA_000001.1','complete','valid',1,1,0,0)")
    db.execute("INSERT INTO assemblies VALUES (2,'GCA_000002.1','failed',NULL,0,0,0,0)")
    db.execute("INSERT INTO sequences(sequence_pk,assembly_pk,sequence_id,length) VALUES (1,1,'contig_1',1000)")
    db.execute("INSERT INTO arrays(array_pk,sequence_pk,tool,tool_version,call_id,start,end,repeat_consensus,category,certainty_score,native_json) VALUES (1,1,'CRISPRidentify','2.3.2','array-1',10,50,'ACGT','Bona-fide',0,'private payload')")
    db.execute("INSERT INTO repeats(repeat_pk,array_pk,ordinal,start,end,length,sequence) VALUES (1,1,1,10,13,4,'ACGT')")
    db.execute("INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,start,end,length,boundary_position,sequence) VALUES (1,1,1,'deletion',NULL,NULL,0,14,NULL)")
    db.commit()
    data = db.serialize()
    db.close()
    return data


def test_normalization_preserves_uncertainty_occurrences_and_relations():
    compressed, registry, counts, size = normalize_database(fixture_database(), 8)
    db = sqlite3.connect(":memory:")
    db.deserialize(gzip.decompress(compressed))
    assert size == len(gzip.decompress(compressed))
    assert counts["assemblies"] == 2
    assert registry[1][4:] == [None, None, None, None]
    assert db.execute("SELECT certainty_score,strand FROM arrays").fetchone() == (0, None)
    assert db.execute("SELECT array_pk,kind,start,end,length,boundary_position,sequence FROM spacers").fetchone() == (8_000_000_001, "deletion", None, None, 0, 14, None)
    assert "native_json" not in {row[1] for row in db.execute("PRAGMA table_info(arrays)")}
    assert db.execute("PRAGMA foreign_key_check").fetchall() == []
    db.close()


def test_selection_drops_superseded_outcome_with_children():
    compressed, registry, counts, _ = normalize_database(fixture_database(), 2, selected_accessions={"GCA_000002.1"})
    assert [row[1] for row in registry] == ["GCA_000002.1"]
    assert counts["arrays"] == counts["repeats"] == counts["spacers"] == 0
    assert gzip.decompress(compressed).startswith(b"SQLite format 3")


def test_archive_checksum_and_manifest_are_independently_validated(tmp_path):
    data = fixture_database()
    manifest = {"database_schema_version": 1, "artifacts": {"results.sqlite3": {"bytes": len(data), "sha256": "0" * 64}}, "counts": {"export_issues": 0}}
    path = tmp_path / "bundle.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, content in {"results.sqlite3": data, "export_manifest.json": json.dumps(manifest).encode(), "collection_index.json": b"{}", "analysis-profile.json": b"{}"}.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="bundle checksum"):
        transform_bundle(path, "0" * 64, 1)
    with pytest.raises(ValueError, match="database checksum"):
        transform_bundle(path, digest, 1)


def test_transfer_rejects_truncation_and_corruption():
    receipt = {"bytes": 4, "sha256": hashlib.sha256(b"abcd").hexdigest()}
    header = json.dumps(receipt).encode() + b"\n"
    with pytest.raises(ValueError, match="truncated"):
        read_frame(io.BytesIO(header + b"ab"))
    with pytest.raises(ValueError, match="checksum"):
        read_frame(io.BytesIO(header + b"xxxx"))


def test_incomplete_snapshot_cannot_be_published_as_complete(tmp_path):
    compressed, registry, counts, size = normalize_database(fixture_database(), 1)
    writer = SnapshotWriter(tmp_path, "test-snapshot")
    writer.accept({"id": 1, "sha256": hashlib.sha256(compressed).hexdigest(), "uncompressed_bytes": size, "counts": counts, "assemblies": registry}, compressed)
    with pytest.raises(ValueError, match="coverage mismatch"):
        writer.publish(expected=3)
    assert not (tmp_path / "catalog.sqlite3").exists()
    summary = writer.publish(expected=2)
    assert summary["counts"]["assemblies"] == 2
    assert "partial" in summary["limitations"][-1]
    assert (tmp_path / "catalog.sqlite3").is_file()
    with pytest.raises(ValueError, match="already published"):
        SnapshotWriter(tmp_path, "another-snapshot")


def test_committed_shard_cannot_be_overwritten_or_published_after_corruption(tmp_path):
    compressed, registry, counts, size = normalize_database(fixture_database(), 1)
    writer = SnapshotWriter(tmp_path, "test-snapshot")
    receipt = {"id": 1, "sha256": hashlib.sha256(compressed).hexdigest(), "uncompressed_bytes": size, "counts": counts, "assemblies": registry}
    writer.accept(receipt, compressed)
    path = tmp_path / "shards/00000001.sqlite3.gz"
    with pytest.raises(ValueError, match="committed shard"):
        writer.accept(receipt, b"bad")
    assert path.read_bytes() == compressed
    path.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum"):
        writer.publish(expected=2)
    assert not (tmp_path / "catalog.sqlite3").exists()
    writer.db.close()


def test_transfer_preserves_production_disk_reserve(tmp_path, monkeypatch):
    compressed, registry, counts, size = normalize_database(fixture_database(), 1)
    writer = SnapshotWriter(tmp_path, "low-disk")
    monkeypatch.setattr("import_catalog.shutil.disk_usage", lambda _path: SimpleNamespace(free=5 * 1024**3))
    with pytest.raises(ValueError, match="5 GiB free-space reserve"):
        writer.accept({"id": 1, "sha256": hashlib.sha256(compressed).hexdigest(), "uncompressed_bytes": size,
                       "counts": counts, "assemblies": registry}, compressed)
    assert writer.db.execute("SELECT count(*) FROM shards").fetchone()[0] == 0
    assert not (tmp_path / "shards/00000001.sqlite3.gz").exists()
    writer.db.close()
