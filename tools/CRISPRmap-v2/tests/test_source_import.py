"""Synthetic source receipts; tests never require the source project or network."""
import gzip
import hashlib
import io
import json
import sqlite3
import tarfile

import pytest

from crisprmap.io import InputError
from crisprmap.source_import import (bounded_json, canonical_digest, open_source, plan_sources,
                                    source_identity, taxonomy)


def make_bundle(tmp_path, first_repeat="AACG", *, accession="GCF_000001.1", assembly_status="complete",
                validation_status="valid", include_features=True, repeat_pks=(1, 2)):
    database = tmp_path / "source.sqlite3"
    connection = sqlite3.connect(database)
    connection.executescript("""
      PRAGMA user_version=1;
      CREATE TABLE assemblies(assembly_pk INTEGER PRIMARY KEY,accession TEXT,status TEXT,validation_status TEXT,sequence_count INTEGER,array_count INTEGER,metadata_json TEXT);
      CREATE TABLE sequences(sequence_pk INTEGER PRIMARY KEY,assembly_pk INTEGER,sequence_id TEXT,length INTEGER,sha256 TEXT);
      CREATE TABLE arrays(array_pk INTEGER PRIMARY KEY,sequence_pk INTEGER,tool TEXT,tool_version TEXT,call_id TEXT,start INTEGER,end INTEGER,strand TEXT,repeat_consensus TEXT,repeat_count INTEGER,spacer_count INTEGER,category TEXT,evidence_level TEXT,certainty_score REAL);
      CREATE TABLE repeats(repeat_pk INTEGER PRIMARY KEY,array_pk INTEGER,ordinal INTEGER,start INTEGER,end INTEGER,sequence TEXT);
    """)
    connection.execute("INSERT INTO assemblies VALUES (1,?,?,?,1,1,?)", (accession, assembly_status, validation_status, json.dumps({"organism": {"taxId": 123, "organismName": "Example species"}, "private": "/not/exported"})))
    connection.execute("INSERT INTO sequences VALUES (1,1,'NC_000001.1',100,?)", ("a" * 64,))
    connection.execute("INSERT INTO arrays VALUES (1,1,'CRISPRidentify','2.3.2','a1',10,23,NULL,'AACG',2,1,'Bona-fide',NULL,0.9)")
    connection.executemany("INSERT INTO repeats VALUES (?,?,?,?,?,?)", [(repeat_pks[0], 1, 1, 10, 13, first_repeat), (repeat_pks[1], 1, 2, 20, 23, "CGTT")])
    if not include_features:
        connection.executescript("DELETE FROM repeats; DELETE FROM arrays; DELETE FROM sequences; UPDATE assemblies SET sequence_count=0,array_count=0;")
    connection.commit()
    connection.close()
    profile = {"schema": "synthetic-profile", "orientation": False}
    profile_sha = canonical_digest(profile)
    data = database.read_bytes()
    collection = {"analysis_profile_sha256": profile_sha, "entries": [
        {"accession": accession, "status": assembly_status, "sha256": "b" * 64},
    ]}
    files = {
        "results.sqlite3": data,
        "export_manifest.json": json.dumps({"database_schema_version": 1, "status": "complete", "commit_marker": "export_manifest.json", "counts": {"export_issues": 0},
                                             "artifacts": {"results.sqlite3": {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}}}).encode(),
        "collection_index.json": json.dumps(collection).encode(),
        "analysis-profile.json": json.dumps({"analysis_profile": profile, "analysis_profile_sha256": profile_sha}).encode(),
    }
    path = tmp_path / "bundle.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, content in files.items():
            member = tarfile.TarInfo(name)
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))
    source = {"kind": "bundle", "path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "selected": [accession]}
    return source, profile_sha


def test_bundle_preserves_selected_occurrences_and_allowlists_metadata(tmp_path):
    source, profile = make_bundle(tmp_path)
    with open_source(source, profile, tmp_path, 0) as opened:
        assemblies = list(opened.assemblies())
        assert assemblies[0]["taxonomy"] == (123, "Example species")
        assert assemblies[0]["source_document_sha256"] == "b" * 64
        assert [row["sequence"] for row in opened.repeats()] == ["AACG", "CGTT"]
        assert list(opened.arrays())[0]["strand"] is None
        assert set(opened.metadata) == {"source_key", "kind", "source_sha256", "selection_sha256", "database_sha256", "collection_index_sha256", "profile_sha256", "profile_file_sha256"}
        assert "/" not in json.dumps(opened.metadata)
    assert not list(tmp_path.glob(".crisprmap-source-*"))


def test_source_digest_and_profile_fail_closed(tmp_path):
    source, profile = make_bundle(tmp_path)
    with pytest.raises(InputError, match="profile"):
        with open_source(source, "0" * 64, tmp_path, 0):
            pass
    source["sha256"] = "0" * 64
    with pytest.raises(InputError, match="checksum"):
        with open_source(source, profile, tmp_path, 0):
            pass


def test_source_selection_identity_ignores_private_location_and_host():
    source = {"kind": "bundle", "path": "/operator/a", "sha256": "1" * 64, "selected": ["GCF_000001.1"], "node": "private-host"}
    assert source_identity(source) == source_identity({**source, "path": "/operator/b", "node": "other-host"})
    assert source_identity(source) != source_identity({**source, "selected": ["GCF_000002.1"]})


def test_plan_flattens_transfer_lanes_and_raw_groups():
    source = {"kind": "raw", "path": "/operator/raw", "sha256": "1" * 64,
              "entries": [{"accession": "GCF_000001.1", "status": "complete", "sha256": "1" * 64}]}
    assert list(plan_sources({"nodes": {"private-name": [{"kind": "raw_group", "items": [source]}]}})) == [source]


def test_unavailable_receipts_are_outcomes_not_empty_successes(tmp_path):
    payload = {"unavailable": {"GCF_000003.1": {"reason": "missing"}}}
    path = tmp_path / "receipt.json.gz"
    path.write_bytes(gzip.compress(json.dumps(payload).encode(), mtime=0))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    source = {"kind": "unavailable", "path": str(path), "sha256": sha,
              "entries": [{"accession": "GCF_000003.1", "status": "unavailable", "sha256": sha}]}
    with open_source(source, "1" * 64, tmp_path, 0) as opened:
        assert list(opened.assemblies())[0]["status"] == "unavailable"
        assert list(opened.arrays()) == []
        assert list(opened.repeats()) == []


def test_taxonomy_does_not_infer_missing_taxon_or_copy_native_fields():
    assert taxonomy({"organism": {"organismName": "Unknown", "taxId": None}, "path": "/secret"}) == (None, "Unknown")
    assert taxonomy({"organism": {"taxId": True}}) == (None, None)


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("payload", [b'{"status":"failed","status":"complete"}',
                                    b'{"certainty_score":NaN}', b'{"certainty_score":1e999}',
                                    b'[]'])
def test_source_json_rejects_ambiguous_or_nonfinite_content(tmp_path, payload, compressed):
    path = tmp_path / ("source.json.gz" if compressed else "source.json")
    path.write_bytes(gzip.compress(payload, mtime=0) if compressed else payload)
    with pytest.raises(InputError):
        bounded_json(path)


@pytest.mark.parametrize("artifact", ["export_manifest.json", "collection_index.json", "analysis-profile.json"])
@pytest.mark.parametrize("mutation", ["duplicate", "nonfinite", "not_object"])
def test_bundle_metadata_uses_strict_json_and_cleans_scratch(tmp_path, artifact, mutation):
    source, profile = make_bundle(tmp_path)
    path = tmp_path / "rewritten.tar.gz"
    with tarfile.open(source["path"], "r:gz") as original, tarfile.open(path, "w:gz") as rewritten:
        for member in original:
            payload = original.extractfile(member).read()
            if member.name == artifact:
                if mutation == "duplicate":
                    document = json.loads(payload)
                    key = next(iter(document))
                    prefix = json.dumps({key: document[key]}).encode()[:-1] + b","
                    payload = prefix + payload[1:]
                elif mutation == "nonfinite":
                    payload = b'{"unused":NaN,' + payload[1:]
                else:
                    payload = b"[]"
                member.size = len(payload)
            rewritten.addfile(member, io.BytesIO(payload))
    source.update(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    with pytest.raises(InputError):
        with open_source(source, profile, tmp_path, 0):
            pytest.fail("Malformed source evidence reached the catalogue importer")
    assert not list(tmp_path.glob(".crisprmap-source-*"))


@pytest.mark.parametrize("payload", ['{"taxId":1,"taxId":2}', '{"taxId":NaN}'])
def test_serialized_taxonomy_rejects_conflicting_or_nonfinite_metadata(payload):
    with pytest.raises(InputError, match="source assembly metadata"):
        taxonomy(payload)


@pytest.mark.parametrize("selected", [[[]], [{}], [True], ["GCF_000001.1", "GCF_000001.1"]])
def test_malformed_selection_raises_input_error(selected):
    with pytest.raises(InputError):
        source_identity({"kind": "bundle", "sha256": "1" * 64, "selected": selected})


class ConnectionWithoutSetlimit:
    """Exercise the Python3.10 fallback on newer Python runtimes."""
    def __init__(self, connection):
        self.connection = connection

    def __getattr__(self, name):
        if name == "setlimit":
            raise AttributeError(name)
        return getattr(self.connection, name)


def hide_connection_setlimit(monkeypatch):
    original = sqlite3.connect
    monkeypatch.setattr(sqlite3, "connect", lambda *args, **kwargs: ConnectionWithoutSetlimit(original(*args, **kwargs)))


def test_python310_sqlite_fallback_preserves_source_rows(tmp_path, monkeypatch):
    source, profile = make_bundle(tmp_path)
    hide_connection_setlimit(monkeypatch)
    with open_source(source, profile, tmp_path, 0) as opened:
        assert list(opened.assemblies())[0]["taxonomy"] == (123, "Example species")
        assert [row["sequence"] for row in opened.repeats()] == ["AACG", "CGTT"]
        assert list(opened.contigs())[0]["sequence_sha256"] == "a" * 64
        assert list(opened.arrays())[0]["repeat_consensus"] == "AACG"


@pytest.mark.parametrize("value", ["é" * (4 * 1024**2 + 1), b"x" * (8 * 1024**2 + 1)])
def test_python310_fallback_rejects_oversized_text_and_blob_before_row_fetch(tmp_path, monkeypatch, value):
    source, profile = make_bundle(tmp_path, first_repeat=value)
    hide_connection_setlimit(monkeypatch)
    with pytest.raises(InputError, match="projected row exceeds"):
        with open_source(source, profile, tmp_path, 0):
            pytest.fail("Oversized scientific values reached a Python source reader")
    assert not list(tmp_path.glob(".crisprmap-source-*"))


@pytest.mark.parametrize("taxon,expected", [("9" * 5000, None), (str(2**63), None),
                                           (str(2**63 - 1), 2**63 - 1), ("0" * 5000 + "123", 123)])
def test_decimal_taxonomy_conversion_is_bounded_and_preserves_valid_ids(taxon, expected):
    assert taxonomy({"taxId": taxon, "organism_name": "Example"}) == (expected, "Example")


@pytest.mark.parametrize("plan", [[], {"nodes": {"lane": None}}, {"nodes": {"lane": {}}},
                                   {"sources": [{"kind": "raw_group", "items": [None]}]}])
def test_malformed_nested_plan_is_an_input_error(plan):
    with pytest.raises(InputError):
        list(plan_sources(plan))


def test_bundle_repeat_order_uses_array_ordinals_and_preserves_native_keys(tmp_path):
    source, profile = make_bundle(tmp_path, repeat_pks=(2, 1))
    with open_source(source, profile, tmp_path, 0) as opened:
        assert [(row["native_pk"], row["ordinal"], row["sequence"]) for row in opened.repeats()] == [
            (2, 1, "AACG"), (1, 2, "CGTT")]


@pytest.mark.parametrize("location", ["assembly", "provenance", "validation", "summary", "sequences",
    "contig", "tools", "tool", "arrays", "array", "array_interval", "confidence", "repeats", "repeat", "repeat_interval"])
def test_malformed_raw_shape_is_an_input_error_and_rolls_back(tmp_path, location):
    from pathlib import Path
    from crisprmap.catalog import CatalogWriter
    from test_catalog import raw_plan

    _, profile, source = raw_plan(tmp_path)
    path = Path(source["path"])
    payload = json.loads(gzip.decompress(path.read_bytes()))
    contig = payload["sequences"][0]
    tool = contig["tools"][0]
    array = tool["arrays"][0]
    repeat = array["repeats"][0]
    parent, key = {
        "assembly": (payload, "assembly"), "provenance": (payload, "provenance"),
        "validation": (payload, "validation"), "summary": (payload, "summary"),
        "sequences": (payload, "sequences"), "contig": (payload["sequences"], 0),
        "tools": (contig, "tools"), "tool": (contig["tools"], 0), "arrays": (tool, "arrays"),
        "array": (tool["arrays"], 0), "array_interval": (array, "interval"),
        "confidence": (array, "confidence"), "repeats": (array, "repeats"),
        "repeat": (array["repeats"], 0), "repeat_interval": (repeat, "interval"),
    }[location]
    parent[key] = None
    path.write_bytes(gzip.compress(json.dumps(payload).encode(), mtime=0))
    source["sha256"] = source["entries"][0]["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    writer = CatalogWriter(tmp_path / "catalog", profile, "3" * 64, "partition", 0)
    try:
        with pytest.raises(InputError, match="Raw source"):
            writer.ingest(source)
        for table in ("sources", "assemblies", "contigs", "nodes", "arrays", "repeat_occurrences"):
            assert writer.db.execute(f"SELECT count(*) FROM {table}").fetchone()[0] == 0
    finally:
        writer.close()
