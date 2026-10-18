"""Catalogue provenance, occurrence reconstruction, merge and recovery fixtures."""
import gzip
import hashlib
import json
import sqlite3

import pytest

from crisprmap.catalog import (build_partition, load_catalog, merge_catalog, normalized_sequence, reconstruct_sequence)
from crisprmap.io import InputError


def raw_plan(tmp_path, accession="GCF_000001.1", category="Bona-fide", source_name="source"):
    profile = "1" * 64
    payload = {"assembly": {"accession": accession, "ncbi_metadata": {"organism": {"taxId": 123, "organismName": "Example species"}}},
               "status": "complete", "validation": {"status": "valid"},
               "provenance": {"analysis_profile_sha256": profile, "private_root": "/never/published"},
               "summary": {"array_count": 1}, "sequences": [{"id": "NC_000001.1", "length": 100, "sha256": "a" * 64,
               "tools": [{"name": "CRISPRidentify", "version": "2.3.2", "arrays": [{"call_id": "a1", "interval": {"start": 10, "end": 23},
                   "strand": None, "repeat_consensus": "AACG", "repeat_count": 2, "spacer_count": 1,
                   "confidence": {"category": category, "score": 0.9},
                   "repeats": [{"index": 1, "interval": {"start": 10, "end": 13}, "sequence": "AACG"},
                               {"index": 2, "interval": {"start": 20, "end": 23}, "sequence": "CGTT"}]}]}]}]}
    path = tmp_path / f"{source_name}.json.gz"
    path.write_bytes(gzip.compress(json.dumps(payload).encode(), mtime=0))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    source = {"kind": "raw", "path": str(path), "sha256": sha,
              "entries": [{"accession": accession, "status": "complete", "sha256": sha}]}
    plan = tmp_path / f"{source_name}-plan.json"
    plan.write_text(json.dumps({"profile_sha256": profile, "sources": [source]}))
    return plan, profile, source


def test_all_occurrences_share_canonical_node_without_losing_orientation(tmp_path):
    plan, profile, _ = raw_plan(tmp_path)
    output = tmp_path / "catalog"
    result = build_partition(plan, output, profile, min_free_bytes=0)
    assert result["counts"]["arrays"] == 1
    assert result["counts"]["repeat_occurrences"] == 2
    assert result["counts"]["nodes"] == 1
    db, _ = load_catalog(output)
    try:
        assert db.execute("SELECT sequence,bona_fide_consensus_count,bona_fide_repeat_count FROM nodes").fetchone() == ("AACG", 1, 2)
        rows = db.execute("SELECT n.sequence,r.orientation,r.raw_sequence FROM repeat_occurrences r JOIN nodes n USING(node_id) ORDER BY r.ordinal").fetchall()
        assert [reconstruct_sequence(*row) for row in rows] == ["AACG", "CGTT"]
        assert db.execute("SELECT tax_id,organism_name FROM assemblies").fetchone() == (123, "Example species")
        with pytest.raises(sqlite3.OperationalError):
            db.execute("DELETE FROM nodes")
    finally:
        db.close()
    assert "/never/published" not in (output / "catalog.sqlite3").read_bytes().decode("latin1")
    assert build_partition(plan, output, profile, min_free_bytes=0) == result


def test_merge_rebinds_ids_and_keeps_tiers_separate(tmp_path):
    first_plan, profile, _ = raw_plan(tmp_path, source_name="first")
    second_plan, _, _ = raw_plan(tmp_path, "GCF_000002.1", "Possible", "second")
    first, second = tmp_path / "first-part", tmp_path / "second-part"
    build_partition(first_plan, first, profile, min_free_bytes=0)
    build_partition(second_plan, second, profile, min_free_bytes=0)
    merged = tmp_path / "merged"
    result = merge_catalog([first, second], merged, min_free_bytes=0)
    assert result["counts"]["assemblies"] == 2
    assert result["counts"]["repeat_occurrences"] == 4
    assert result["counts"]["nodes"] == 1
    db, _ = load_catalog(merged)
    try:
        assert db.execute("SELECT bona_fide_consensus_count,possible_consensus_count,bona_fide_repeat_count,possible_repeat_count FROM nodes").fetchone() == (1, 1, 2, 2)
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
    finally:
        db.close()


def test_writable_catalogue_attaches_escaped_partition_uri_readonly(tmp_path):
    from crisprmap.catalog import _connect
    partition = tmp_path / "partition #? with spaces.sqlite3"
    with sqlite3.connect(partition) as source:
        source.execute("CREATE TABLE sample(value INTEGER)")
        source.execute("INSERT INTO sample VALUES (7)")
    connection = _connect(tmp_path / "building.sqlite3")
    try:
        connection.execute("ATTACH DATABASE ? AS incoming", (partition.as_uri() + "?mode=ro&immutable=1",))
        assert connection.execute("SELECT value FROM incoming.sample").fetchone() == (7,)
        with pytest.raises(sqlite3.OperationalError, match="readonly"):
            connection.execute("INSERT INTO incoming.sample VALUES (8)")
    finally:
        connection.close()
    with sqlite3.connect(partition) as source:
        assert source.execute("SELECT value FROM sample").fetchall() == [(7,)]


def test_changed_source_fails_before_commit_then_resumes(tmp_path):
    plan, profile, source = raw_plan(tmp_path)
    from pathlib import Path
    path = Path(source["path"])
    original = path.read_bytes()
    path.write_bytes(original + b"corrupt")
    output = tmp_path / "catalog"
    with pytest.raises(InputError, match="checksum"):
        build_partition(plan, output, profile, min_free_bytes=0)
    assert not (output / "manifest.json").exists()
    db = sqlite3.connect(output / "catalog.building.sqlite3")
    assert db.execute("SELECT count(*) FROM sources").fetchone()[0] == 0
    db.close()
    path.write_bytes(original)
    assert build_partition(plan, output, profile, min_free_bytes=0)["counts"]["repeat_occurrences"] == 2


def test_sequence_normalization_retains_original_ambiguity_and_rna():
    for original in ["aacg", "ACGU", "ACNG", "ACGT"]:
        canonical, orientation, raw = normalized_sequence(original)
        assert reconstruct_sequence(canonical, orientation, raw) == original
    with pytest.raises(InputError):
        normalized_sequence("AC-G")


def test_duplicate_assembly_across_parts_is_rejected(tmp_path):
    plan, profile, _ = raw_plan(tmp_path)
    first = tmp_path / "part"
    build_partition(plan, first, profile, min_free_bytes=0)
    with pytest.raises(InputError, match="Duplicate partition"):
        merge_catalog([first, first], tmp_path / "merged", min_free_bytes=0)


def test_resume_after_database_rename_before_final_manifest(tmp_path, monkeypatch):
    import crisprmap.catalog as catalog
    plan, profile, _ = raw_plan(tmp_path)
    output = tmp_path / "catalog"
    real_write = catalog.atomic_json
    def interrupted(*_args, **_kwargs):
        raise OSError("simulated publication interruption")
    monkeypatch.setattr(catalog, "atomic_json", interrupted)
    with pytest.raises(OSError, match="interruption"):
        build_partition(plan, output, profile, min_free_bytes=0)
    assert (output / "catalog.sqlite3").exists()
    assert not (output / "manifest.json").exists()
    monkeypatch.setattr(catalog, "atomic_json", real_write)
    assert build_partition(plan, output, profile, min_free_bytes=0)["counts"]["repeat_occurrences"] == 2


@pytest.mark.parametrize("reserve", [-1, True, float("nan"), float("inf"), 1.5, 2**63])
def test_invalid_free_reserve_rejected_before_build_or_merge(tmp_path, reserve):
    from crisprmap.source_import import check_free
    with pytest.raises(InputError, match="budgets"):
        build_partition(tmp_path / "absent-plan", tmp_path / "absent-output", "1" * 64, min_free_bytes=reserve)
    with pytest.raises(InputError, match="budgets"):
        merge_catalog([], tmp_path / "absent-output", min_free_bytes=reserve)
    with pytest.raises(InputError, match="budgets"):
        check_free(tmp_path, reserve)
    assert not (tmp_path / "absent-output").exists()


def test_index_publication_reserve_failure_remains_resumable(tmp_path, monkeypatch):
    import crisprmap.catalog as catalog
    plan, profile, _ = raw_plan(tmp_path)
    output = tmp_path / "catalog"
    original_connect = catalog._connect
    original_check = catalog.check_free
    state = {"index_started": False, "raised": False}

    def traced_connect(path, readonly=False):
        connection = original_connect(path, readonly=readonly)
        if not readonly:
            connection.set_trace_callback(lambda statement: state.update(index_started=True)
                                          if statement.lstrip().upper().startswith("CREATE INDEX") else None)
        return connection

    def checked(directory, minimum, extra=0):
        if state["index_started"] and not state["raised"]:
            state["raised"] = True
            raise InputError("simulated publication reserve exhaustion")
        return original_check(directory, minimum, extra)

    monkeypatch.setattr(catalog, "_connect", traced_connect)
    monkeypatch.setattr(catalog, "check_free", checked)
    with pytest.raises(InputError, match="reserve exhaustion"):
        build_partition(plan, output, profile, min_free_bytes=0)
    assert not (output / "manifest.json").exists()
    assert (output / "catalog.building.sqlite3").exists()
    monkeypatch.setattr(catalog, "check_free", original_check)
    result = build_partition(plan, output, profile, min_free_bytes=0)
    assert result["counts"]["assemblies"] == 1
    assert result["counts"]["arrays"] == 1
    assert result["counts"]["repeat_occurrences"] == 2


@pytest.mark.parametrize("status,validation", [("failed", None), ("failed", "valid"), ("complete", None), ("complete", "invalid")])
def test_unvalidated_assembly_features_roll_back_source_and_node_counts(tmp_path, status, validation):
    from crisprmap.catalog import CatalogWriter
    from test_source_import import make_bundle
    first, second = tmp_path / "first", tmp_path / "second"
    first.mkdir()
    second.mkdir()
    good, profile = make_bundle(first)
    bad, _ = make_bundle(second, accession="GCF_000002.1", assembly_status=status, validation_status=validation)
    writer = CatalogWriter(tmp_path / "catalog", profile, "3" * 64, "partition", 0)
    try:
        writer.ingest(good)
        before = writer.db.execute("SELECT sequence,bona_fide_consensus_count,bona_fide_repeat_count FROM nodes ORDER BY node_id").fetchall()
        with pytest.raises(InputError, match="validat"):
            writer.ingest(bad)
        assert writer.db.execute("SELECT count(*) FROM sources").fetchone()[0] == 1
        assert writer.db.execute("SELECT count(*) FROM assemblies").fetchone()[0] == 1
        assert writer.db.execute("SELECT count(*) FROM arrays").fetchone()[0] == 1
        assert writer.db.execute("SELECT count(*) FROM repeat_occurrences").fetchone()[0] == 2
        assert writer.db.execute("SELECT sequence,bona_fide_consensus_count,bona_fide_repeat_count FROM nodes ORDER BY node_id").fetchall() == before
        assert writer.db.execute("PRAGMA foreign_key_check").fetchall() == []
    finally:
        writer.close()


def test_failed_receipt_without_features_remains_an_accepted_outcome(tmp_path):
    from crisprmap.catalog import CatalogWriter
    from test_source_import import make_bundle
    source, profile = make_bundle(tmp_path, assembly_status="failed", validation_status=None, include_features=False)
    writer = CatalogWriter(tmp_path / "catalog", profile, "3" * 64, "partition", 0)
    try:
        writer.ingest(source)
        assert writer.db.execute("SELECT status FROM assemblies").fetchone()[0] == "failed"
        assert writer.db.execute("SELECT count(*) FROM nodes").fetchone()[0] == 0
        assert writer.db.execute("SELECT count(*) FROM arrays").fetchone()[0] == 0
    finally:
        writer.close()



def test_bundle_import_accepts_native_keys_out_of_ordinal_order(tmp_path):
    from crisprmap.catalog import CatalogWriter
    from test_source_import import make_bundle

    source, profile = make_bundle(tmp_path, repeat_pks=(2, 1))
    writer = CatalogWriter(tmp_path / "catalog", profile, "3" * 64, "partition", 0)
    try:
        assert writer.ingest(source)
        rows = writer.db.execute("SELECT source_native_pk,ordinal FROM repeat_occurrences ORDER BY ordinal").fetchall()
        assert rows == [(2, 1), (1, 2)]
        assert writer.db.execute("SELECT bona_fide_repeat_count FROM nodes").fetchone()[0] == 2
    finally:
        writer.close()
