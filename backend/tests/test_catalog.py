from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import sqlite3
import threading

from fastapi.testclient import TestClient
import pytest

from crispr_web.api import create_app
from crispr_web.catalog import Catalog, CatalogCursorError, CatalogUnavailable, _MAX_SHARD_BYTES, _configure
from crispr_web.catalog_terms import TERM_VERSION, gene_matches, query_hash, terms_for


def _summary(snapshot="test-snapshot"):
    return {
        "available": True,
        "snapshot_id": snapshot,
        "title": "Independent test reference",
        "counts": {"assemblies": 2, "arrays": 2, "repeats": 1, "spacers": 2,
                   "cas_systems": 0, "cas_genes": 1},
        "coverage": {"catalog_assemblies": 10, "imported_assemblies": 2,
                     "successful_assemblies": 9, "failed_assemblies": 1},
        "limitations": ["Imported coverage is partial."],
        "provenance": {"release": "test", "shards_sha256": "c" * 64},
    }


def _database(path, offset=0):
    # The in-repository importer owns the normalized schema; no predecessor data
    # or scientific runtime is needed to create this fully synthetic fixture.
    schema = runpy.run_path(str(Path(__file__).resolve().parents[2] / "scripts/catalog_snapshot.py"))
    with sqlite3.connect(path) as connection:
        schema["create_tables"](connection)
        connection.execute("INSERT INTO export_metadata VALUES ('summary',?)", (json.dumps(_summary()),))
        connection.executemany(
            "INSERT INTO assemblies VALUES (?,?,?,?,?,?,?,?)",
            [(offset + 1, "GCA_000000001.1", "complete", "valid", 1, 2, 0, 1),
             (offset + 2, "GCA_000000002.1", "failed", None, 0, 0, 0, 0)],
        )
        connection.execute(
            "INSERT INTO sequences(sequence_pk,assembly_pk,sequence_id) VALUES (?,?,?)",
            (offset + 1, offset + 1, "contig_1"),
        )
        connection.executemany(
            "INSERT INTO arrays(array_pk,sequence_pk,tool,tool_version,call_id,start,end,"
            "repeat_consensus,repeat_count,spacer_count,category,certainty_score) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            [(offset + 1, offset + 1, "CRISPRidentify", "2.3.2", "call-1", 1, 15,
              "ACGT", 2, 1, "Bona-fide", 0.84),
             (offset + 2, offset + 1, "CRISPRidentify", "2.3.2", "call-2", 30, 50,
              None, None, None, "Possible", None)],
        )
        connection.execute(
            "INSERT INTO repeats(repeat_pk,array_pk,ordinal,start,end,length,sequence) "
            "VALUES (?,?,0,1,4,4,'ACGT')", (offset + 1, offset + 1),
        )
        connection.executemany(
            "INSERT INTO spacers(spacer_pk,array_pk,ordinal,kind,start,end,length,boundary_position,sequence) "
            "VALUES (?,?,?,?,?,?,?,?,?)",
            [(offset + 1, offset + 1, 0, "observed", 5, 7, 3, None, "TTT"),
             (offset + 2, offset + 1, 1, "deletion", None, None, None, 11, None)],
        )
        connection.execute(
            "INSERT INTO cas_genes(cas_gene_pk,sequence_pk,tool,tool_version,call_id,annotation,"
            "profile,is_cas,score_kind) VALUES (?,?, 'CasAndra','1','gene-1','cas1','Cas1',0,'bitscore')",
            (offset + 1, offset + 1),
        )
        schema["add_indexes"](connection)
    return path


@pytest.fixture
def catalog_path(tmp_path):
    return _database(tmp_path / "catalog.sqlite3")


def test_catalog_is_public_bounded_and_preserves_scientific_missingness(settings_factory, catalog_path):
    with TestClient(create_app(settings_factory(catalog_database_path=catalog_path))) as client:
        summary = client.get("/api/v1/catalog")
        assert summary.status_code == 200
        assert summary.json()["coverage"]["catalog_assemblies"] == 10
        first = client.get("/api/v1/catalog/arrays", params={"limit": 1}).json()
        assert first["items"][0]["id"] == 1
        assert first["items"][0]["certainty_score"] == 0.84
        assert "native_json" not in first["items"][0]
        second = client.get("/api/v1/catalog/arrays", params={"cursor": first["next_cursor"], "limit": 1}).json()
        assert second["items"][0]["repeat_consensus"] is None
        assert second["items"][0]["strand"] is None
        assert second["next_cursor"] is None
        spacers = client.get("/api/v1/catalog/spacers", params={"array_id": 1}).json()["items"]
        assert spacers[1]["kind"] == "deletion"
        assert spacers[1]["sequence"] is None
        assert spacers[1]["boundary_position"] == 11
        gene = client.get("/api/v1/catalog/cas-genes").json()["items"][0]
        assert gene["is_cas"] is False
        assert gene["cas_system_id"] is None
        assert gene["score"] is None
        assert client.get("/api/v1/catalog/cas-systems").json()["items"] == []
        failed = client.get("/api/v1/catalog/assemblies/GCA_000000002.1").json()
        assert failed["status"] == "failed"
        assert failed["validation_status"] is None
        assert client.get("/api/v1/catalog/arrays/1").status_code == 200
        assert client.get("/api/v1/catalog/arrays/999").status_code == 404
        assert client.get("/api/v1/catalog/repeats", params={"q": "ACGT"}).json()["items"][0]["array_id"] == 1
        assert client.get("/api/v1/catalog/arrays", params={"q": "' OR 1=1 --"}).json()["items"] == []
        assert client.get("/api/v1/catalog/arrays", params={"limit": 101}).status_code == 422
        assert client.get("/api/v1/catalog/arrays", params={"q": "a" * 1001}).status_code == 422
        assert client.get("/api/v1/catalog/arrays", params={"cursor": "malformed"}).status_code == 400
        assert client.get("/api/v1/catalog/spacers", params={"cursor": first["next_cursor"]}).status_code == 400
        assert client.post("/api/v1/catalog/arrays", json={}).status_code == 405


def test_catalog_absent_and_invalid_never_exposes_local_paths(settings_factory, tmp_path):
    with TestClient(create_app(settings_factory())) as client:
        assert client.get("/api/v1/catalog").json()["available"] is False
        assert client.get("/api/v1/catalog/arrays").status_code == 503
    private = tmp_path / "private-name.sqlite3"
    private.write_text("not sqlite")
    with TestClient(create_app(settings_factory(catalog_database_path=private))) as client:
        response = client.get("/api/v1/catalog")
        assert response.status_code == 503
        assert "private-name" not in response.text
        assert response.headers["retry-after"] == "30"


def test_catalog_cursor_binds_snapshot_filters_and_symlink_rejection(catalog_path, tmp_path):
    catalog = Catalog(catalog_path)
    cursor = catalog.page("arrays", limit=1)["next_cursor"]
    with pytest.raises(CatalogCursorError):
        catalog.page("arrays", cursor=cursor, assembly_accession="GCA_000000001.1")
    with sqlite3.connect(catalog_path) as connection:
        connection.execute("UPDATE export_metadata SET value=?", (json.dumps(_summary("new-snapshot")),))
    with pytest.raises(CatalogCursorError):
        catalog.page("arrays", cursor=cursor)
    symlink = tmp_path / "linked.sqlite3"
    symlink.symlink_to(catalog_path)
    with pytest.raises(CatalogUnavailable):
        Catalog(symlink).summary()


def test_accession_prefix_pagination_uses_accession_order(catalog_path):
    with sqlite3.connect(catalog_path) as connection:
        connection.execute("INSERT INTO assemblies VALUES (3,'GCA_000000000.1','unavailable',NULL,NULL,NULL,NULL,NULL)")
    catalog = Catalog(catalog_path)
    first = catalog.page("assemblies", q="GCA_", limit=1)
    assert first["items"][0]["id"] == 3
    assert first["items"][0]["array_count"] is None
    second = catalog.page("assemblies", q="GCA_", limit=1, cursor=first["next_cursor"])
    assert second["items"][0]["id"] == 1
    third = catalog.page("assemblies", q="GCA_", limit=1, cursor=second["next_cursor"])
    assert third["items"][0]["id"] == 2
    assert third["next_cursor"] is None


def _sharded(tmp_path, gene_annotation=None):
    root = _database(tmp_path / "root.sqlite3", 1_000_000_000)
    shard = _database(tmp_path / "source.sqlite3", 1_000_000_000)
    if gene_annotation is not None:
        with sqlite3.connect(shard) as connection:
            connection.execute("UPDATE cas_genes SET annotation=?,profile=?", (gene_annotation, gene_annotation))
    raw = shard.read_bytes()
    compressed = gzip.compress(raw)
    (tmp_path / "shard.sqlite3.gz").write_bytes(compressed)
    with sqlite3.connect(root) as connection:
        connection.execute("CREATE TABLE shards(id INTEGER PRIMARY KEY,relative_path TEXT,sha256 TEXT,uncompressed_bytes INTEGER,counts_json TEXT)")
        connection.execute("CREATE TABLE assembly_shards(assembly_pk INTEGER PRIMARY KEY,shard_id INTEGER)")
        connection.execute("INSERT INTO shards VALUES (1,?,?,?,?)", (
            "shard.sqlite3.gz", hashlib.sha256(compressed).hexdigest(), len(raw), json.dumps(_summary()["counts"]),
        ))
        connection.execute("INSERT INTO assembly_shards VALUES (1000000001,1)")
    return root


def test_compressed_shards_work_with_scoped_search_and_verified_bytes(settings_factory, tmp_path):
    root = _sharded(tmp_path)
    with TestClient(create_app(settings_factory(catalog_database_path=root))) as client:
        assert client.get("/api/v1/catalog/arrays").json()["items"][0]["id"] == 1_000_000_001
        assert client.get("/api/v1/catalog/arrays/1000000001").status_code == 200
        assert client.get("/api/v1/catalog/arrays", params={"q": "ACGT"}).status_code == 503
        response = client.get("/api/v1/catalog/arrays", params={"q": "ACGT", "assembly_accession": "GCA_000000001.1"})
        assert response.status_code == 200
        assert len(response.json()["items"]) == 1
        assert len(client.get("/api/v1/catalog/arrays", params={"q": "GCA_000000001.1"}).json()["items"]) == 2
        assert client.get("/api/v1/catalog/spacers", params={"array_id": 1_000_000_001}).json()["items"][1]["kind"] == "deletion"
    (tmp_path / "shard.sqlite3.gz").write_bytes(b"corrupted")
    with pytest.raises(CatalogUnavailable):
        Catalog(root).page("arrays")


def test_shard_traversal_and_expansion_are_rejected(tmp_path):
    root = _sharded(tmp_path)
    with sqlite3.connect(root) as connection:
        connection.execute("UPDATE shards SET relative_path='../outside.sqlite3.gz'")
    with pytest.raises(CatalogUnavailable):
        Catalog(root).page("arrays")
    with sqlite3.connect(root) as connection:
        connection.execute("UPDATE shards SET relative_path='shard.sqlite3.gz',uncompressed_bytes=999999999")
    with pytest.raises(CatalogUnavailable):
        Catalog(root).page("arrays")


@pytest.mark.parametrize("difference", [-1, 1])
def test_shard_streaming_requires_exact_declared_expansion(tmp_path, difference):
    root = _sharded(tmp_path)
    with sqlite3.connect(root) as connection:
        connection.execute("UPDATE shards SET uncompressed_bytes=uncompressed_bytes+?", (difference,))
    with pytest.raises(CatalogUnavailable, match="expansion"):
        Catalog(root).page("arrays")


def test_largest_shards_are_supported_but_cannot_allocate_concurrently(tmp_path):
    catalog = Catalog(tmp_path / "unused.sqlite3")
    started, admitted = threading.Event(), threading.Event()

    def another_large_reader():
        started.set()
        with catalog._admit_shard(_MAX_SHARD_BYTES):
            admitted.set()

    with ThreadPoolExecutor(max_workers=1) as executor:
        with catalog._admit_shard(_MAX_SHARD_BYTES):
            waiting = executor.submit(another_large_reader)
            assert started.wait(1)
            assert not admitted.wait(0.05)
        waiting.result(timeout=2)
    assert admitted.is_set()
    assert catalog._active_shard_bytes == 0


def test_shard_pagination_continues_after_bounded_scan(tmp_path):
    root = _sharded(tmp_path)
    for shard_id in range(2, 6):
        source = _database(tmp_path / f"source-{shard_id}.sqlite3", shard_id * 1_000_000_000)
        raw = source.read_bytes()
        compressed = gzip.compress(raw)
        name = f"shard-{shard_id}.sqlite3.gz"
        (tmp_path / name).write_bytes(compressed)
        with sqlite3.connect(root) as connection:
            connection.execute("INSERT INTO shards VALUES (?,?,?,?,?)", (
                shard_id, name, hashlib.sha256(compressed).hexdigest(), len(raw), json.dumps(_summary()["counts"]),
            ))
    catalog = Catalog(root)
    first = catalog.page("arrays", limit=100)
    assert len(first["items"]) == 8
    assert first["next_cursor"] is not None
    second = catalog.page("arrays", limit=100, cursor=first["next_cursor"])
    assert [item["id"] for item in second["items"]] == [5_000_000_001, 5_000_000_002]
    assert second["next_cursor"] is None
    assert len(catalog._shard_cache) <= 2


def test_catalog_query_work_and_concurrent_shard_memory_are_bounded(tmp_path):
    with sqlite3.connect(":memory:") as connection:
        _configure(connection)
        with pytest.raises(sqlite3.OperationalError, match="interrupted"):
            connection.execute(
                "WITH RECURSIVE long_query(n) AS (VALUES(0) UNION ALL "
                "SELECT n + 1 FROM long_query WHERE n < 10000000) SELECT sum(n) FROM long_query"
            ).fetchone()
    catalog = Catalog(_sharded(tmp_path))
    # Simulate two readers occupying the maximum shard deserialization slots.
    assert catalog._shard_slots.acquire(blocking=False)
    assert catalog._shard_slots.acquire(blocking=False)
    try:
        with pytest.raises(CatalogUnavailable, match="capacity"):
            catalog.page("arrays")
    finally:
        catalog._shard_slots.release()
        catalog._shard_slots.release()
    assert len(catalog.page("arrays")["items"]) == 2


def _search_index(root, hits):
    path = root.with_name("exact-search.sqlite3")
    metadata = {"schema": "crisprloci.exact-search.v1", "term_version": TERM_VERSION,
                "status": "complete", "snapshot_id": "test-snapshot",
                "registry_sha256": hashlib.sha256(root.read_bytes()).hexdigest(), "shards_sha256": "c" * 64}
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT NOT NULL)")
        connection.execute("CREATE TABLE hits(kind INTEGER,key_hash BLOB,shard_id INTEGER,PRIMARY KEY(kind,key_hash,shard_id)) WITHOUT ROWID")
        connection.executemany("INSERT INTO metadata VALUES (?,?)", metadata.items())
        connection.executemany("INSERT INTO hits VALUES (?,?,?)", [(kind, query_hash(kind, term), shard_id) for kind, term, shard_id in hits])
    return path


def test_global_exact_search_verifies_sequences_and_source_gene_aliases(settings_factory, tmp_path):
    root = _sharded(tmp_path, "C25_Cas9_1")
    _search_index(root, [(1, "call-1", 1), (2, "ACGT", 1), (3, "TTT", 1),
                         (5, "cas9", 1), (5, "C25_Cas9_1", 1), (2, "CCCC", 1)])
    with TestClient(create_app(settings_factory(catalog_database_path=root))) as client:
        repeat = client.get("/api/v1/catalog/repeats", params={"q": "acgt"})
        assert repeat.status_code == 200
        assert repeat.json()["items"][0]["sequence"] == "ACGT"
        spacer = client.get("/api/v1/catalog/spacers", params={"q": "ttt"}).json()["items"][0]
        assert spacer["kind"] == "observed" and spacer["sequence"] == "TTT"
        for query in ("Cas9", "c25_CAS9_1"):
            genes = client.get("/api/v1/catalog/cas-genes", params={"q": query})
            assert genes.status_code == 200
            assert genes.json()["items"][0]["annotation"] == "C25_Cas9_1"
        assert client.get("/api/v1/catalog/cas-genes", params={"q": "cas"}).json()["items"] == []
        assert client.get("/api/v1/catalog/cas-genes", params={"q": "cas90"}).json()["items"] == []
        # Even a false candidate or a truncated-hash collision cannot produce
        # a scientific match without exact verification in the original shard.
        assert client.get("/api/v1/catalog/repeats", params={"q": "CCCC"}).json()["items"] == []
        assert client.get("/api/v1/catalog/arrays", params={"q": "call-1"}).json()["items"][0]["call_id"] == "call-1"


@pytest.mark.parametrize("match", ["exact", "similar"])
def test_global_search_advances_over_bounded_false_candidates(tmp_path, match):
    root = _sharded(tmp_path)
    for shard_id in range(2, 6):
        source = _database(tmp_path / f"search-source-{shard_id}.sqlite3", shard_id * 1_000_000_000)
        if shard_id == 5:
            with sqlite3.connect(source) as connection:
                connection.execute("UPDATE repeats SET sequence='CCCC'")
        raw = source.read_bytes()
        compressed = gzip.compress(raw)
        name = f"search-shard-{shard_id}.sqlite3.gz"
        (tmp_path / name).write_bytes(compressed)
        with sqlite3.connect(root) as connection:
            connection.execute("INSERT INTO shards VALUES (?,?,?,?,?)", (
                shard_id, name, hashlib.sha256(compressed).hexdigest(), len(raw), json.dumps(_summary()["counts"]),
            ))
    _search_index(root, [(2, "CCCC", shard_id) for shard_id in range(1, 6)])
    catalog = Catalog(root)
    first = catalog.page("repeats", q="CCCC", match=match)
    assert first["items"] == [] and first["next_cursor"] is not None
    second = catalog.page("repeats", q="CCCC", cursor=first["next_cursor"], match=match)
    assert second["items"][0]["id"] == 5_000_000_001
    assert second["next_cursor"] is None
    with pytest.raises(CatalogCursorError):
        catalog.page("repeats", q="ACGT", cursor=first["next_cursor"])


def test_similar_search_uses_index_and_verifies_each_candidate(settings_factory, tmp_path):
    root = _sharded(tmp_path)
    _search_index(root, [(2, "ACGT", 1), (2, "CCCC", 1), (3, "TTT", 1)])
    with TestClient(create_app(settings_factory(catalog_database_path=root))) as client:
        def search(kind, query, **kwargs):
            return client.get(f"/api/v1/catalog/{kind}", params={"q": query, "match": "similar", **kwargs})

        for query in ("ACGT", "acga", "ACGN"):
            response = search("repeats", query)
            assert response.status_code == 200
            assert [item["sequence"] for item in response.json()["items"]] == ["ACGT"]
        # Only one literal symbol may differ; N is not a wildcard.
        for query in ("ACNN", "CCCC", "ACG", "ACGTT"):
            assert search("repeats", query).json()["items"] == []
        for query in ("AAA", "AAT"):
            rows = search("spacers", query).json()["items"]
            assert [row["sequence"] for row in rows] == ["TTT"]  # Reverse complement.
            assert all(row["kind"] == "observed" for row in rows)
        for kind, query in (("arrays", "ACGT"), ("repeats", "A-CG"), ("spacers", "A" * 1001), ("spacers", "")):
            assert search(kind, query).status_code == 422
        assert client.get("/api/v1/catalog/repeats", params={"match": "similar"}).status_code == 422
        assert search("repeats", "ACGA", assembly_accession="GCA_000000002.1").json()["items"] == []


def test_similar_pagination_binds_match_policy_and_preserves_occurrences(catalog_path):
    with sqlite3.connect(catalog_path) as connection:
        connection.executemany(
            "INSERT INTO repeats(repeat_pk,array_pk,ordinal,sequence) VALUES (?,1,?,'ACGA')",
            [(2, 1), (3, 2)],
        )
    catalog = Catalog(catalog_path)
    first = catalog.page("repeats", q="ACGT", match="similar", limit=1)
    second = catalog.page("repeats", q="ACGT", match="similar", limit=1, cursor=first["next_cursor"])
    third = catalog.page("repeats", q="ACGT", match="similar", limit=1, cursor=second["next_cursor"])
    assert [page["items"][0]["id"] for page in (first, second, third)] == [1, 2, 3]
    assert third["next_cursor"] is None
    with pytest.raises(CatalogCursorError):
        catalog.page("repeats", q="ACGT", cursor=first["next_cursor"])


def test_similar_search_fails_closed_without_snapshot_index(tmp_path):
    catalog = Catalog(_sharded(tmp_path))
    with pytest.raises(CatalogUnavailable):
        catalog.page("repeats", q="ACGA", match="similar")


def test_global_search_accepts_observed_spacers_longer_than_120_bases(settings_factory, tmp_path):
    root = _sharded(tmp_path)
    sequence = "ACGT" * 35
    source = tmp_path / "source.sqlite3"
    with sqlite3.connect(source) as connection:
        connection.execute("UPDATE spacers SET sequence=?,length=140,end=144 WHERE kind='observed'", (sequence,))
        connection.execute("UPDATE arrays SET end=152 WHERE array_pk=1000000001")
    raw = source.read_bytes()
    compressed = gzip.compress(raw)
    (tmp_path / "shard.sqlite3.gz").write_bytes(compressed)
    with sqlite3.connect(root) as connection:
        connection.execute("UPDATE shards SET sha256=?,uncompressed_bytes=?", (hashlib.sha256(compressed).hexdigest(), len(raw)))
    _search_index(root, [(3, sequence, 1)])
    with TestClient(create_app(settings_factory(catalog_database_path=root))) as client:
        response = client.get("/api/v1/catalog/spacers", params={"q": sequence.lower()})
        assert response.status_code == 200
        assert response.json()["items"][0]["sequence"] == sequence
        assert client.get("/api/v1/catalog/spacers", params={"q": sequence[:120]}).json()["items"] == []


def test_global_dna_search_normalizes_case_without_ambiguity_wildcards(settings_factory, tmp_path):
    root = _sharded(tmp_path)
    source = tmp_path / "source.sqlite3"
    with sqlite3.connect(source) as connection:
        connection.execute("UPDATE repeats SET sequence='aCgT'")
        connection.execute("UPDATE arrays SET repeat_consensus='aCgT' WHERE array_pk=1000000001")
        connection.execute("UPDATE spacers SET sequence='tTt' WHERE kind='observed'")
    raw = source.read_bytes()
    compressed = gzip.compress(raw)
    (tmp_path / "shard.sqlite3.gz").write_bytes(compressed)
    with sqlite3.connect(root) as connection:
        connection.execute("UPDATE shards SET sha256=?,uncompressed_bytes=?", (hashlib.sha256(compressed).hexdigest(), len(raw)))
    _search_index(root, [(1, "ACGT", 1), (2, "ACGT", 1), (3, "TTT", 1), (2, "ACGN", 1)])
    with TestClient(create_app(settings_factory(catalog_database_path=root))) as client:
        for entity, query, field, expected in (
            ("arrays", "ACGT", "repeat_consensus", "aCgT"),
            ("repeats", "ACGT", "sequence", "aCgT"),
            ("spacers", "TTT", "sequence", "tTt"),
        ):
            response = client.get("/api/v1/catalog/" + entity, params={"q": query})
            assert response.status_code == 200
            assert response.json()["items"][0][field] == expected
        assert client.get("/api/v1/catalog/repeats", params={"q": "ACGN"}).json()["items"] == []


@pytest.mark.parametrize("key,value", [
    ("status", "building"), ("schema", "unknown"), ("term_version", "old"),
    ("snapshot_id", "other"), ("registry_sha256", "0" * 64), ("shards_sha256", "0" * 64),
])
def test_global_search_rejects_unpinned_or_incomplete_sidecar(tmp_path, key, value):
    root = _sharded(tmp_path)
    index = _search_index(root, [(2, "ACGT", 1)])
    with sqlite3.connect(index) as connection:
        connection.execute("UPDATE metadata SET value=? WHERE key=?", (value, key))
    catalog = Catalog(root)
    with pytest.raises(CatalogUnavailable):
        catalog.page("repeats", q="ACGT")
    assert catalog.page("repeats", q="ACGT", assembly_accession="GCA_000000001.1")["items"]


def test_global_search_rejects_missing_referenced_shards_and_symlink_sidecars(tmp_path):
    root = _sharded(tmp_path)
    index = _search_index(root, [(2, "ACGT", 99)])
    with pytest.raises(CatalogUnavailable, match="unavailable shards"):
        Catalog(root).page("repeats", q="ACGT")
    alternate = index.with_name("alternate.sqlite3")
    index.rename(alternate)
    index.symlink_to(alternate)
    with pytest.raises(CatalogUnavailable):
        Catalog(root).page("repeats", q="ACGT")


def test_search_summary_updates_runtime_limitation_without_modifying_snapshot(tmp_path):
    root = _sharded(tmp_path)
    obsolete = "Exact repeat, spacer and Cas profile searches require a selected genome. No sequence alignment or similarity search is performed."
    summary = _summary()
    summary["limitations"] = [obsolete]
    with sqlite3.connect(root) as connection:
        connection.execute("UPDATE export_metadata SET value=? WHERE key='summary'", (json.dumps(summary),))
    assert Catalog(root).summary().limitations == [obsolete]
    before = root.read_bytes()
    _search_index(root, [(2, "ACGT", 1)])
    assert "across the snapshot" in Catalog(root).summary().limitations[0]
    assert root.read_bytes() == before


def test_gene_aliases_preserve_source_tokens_without_inferred_family_expansion():
    assert terms_for(5, "C25_Cas5_cmr3_1") == {"c25_cas5_cmr3_1", "cas5", "cmr3"}
    assert gene_matches("C25_Cas9_1", None, "cAs9")
    assert not gene_matches("C25_Cas5f1_1", None, "Cas5")
    assert gene_matches("C25_Cas5f1_1", None, "Cas5f1")
    assert not gene_matches("C25_Cas9_1", None, "Cas")
