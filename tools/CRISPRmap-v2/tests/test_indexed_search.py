import json
from itertools import product
import sqlite3

import pytest

from crisprmap import catalog
from crisprmap.catalog_index import build_index, file_sha256
from crisprmap.core import canonical_sequence
from crisprmap.indexed_search import _distance_backend, search_catalog
from crisprmap.io import InputError, LimitError, digest


def fixture(tmp_path, monkeypatch, sequences=None, core_ids=None):
    source = tmp_path / "catalog"
    source.mkdir()
    path = source / "catalog.sqlite3"
    with sqlite3.connect(path) as connection:
        connection.execute("""CREATE TABLE nodes(node_id INTEGER PRIMARY KEY,sequence TEXT UNIQUE,
            sequence_sha256 TEXT UNIQUE,length INTEGER,eligible INTEGER,bona_fide_consensus_count INTEGER,
            possible_consensus_count INTEGER,bona_fide_repeat_count INTEGER,possible_repeat_count INTEGER)""")
        for node_id, raw in enumerate(sequences or ["AACCG", "AACCT", "AAGCT", "CGCGC"], 1):
            sequence = canonical_sequence(raw)
            connection.execute("INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?)",
                               (node_id, sequence, digest(sequence), len(sequence), 1,
                                int(core_ids is None or node_id in core_ids), 2, 3, 4))
    manifest = {"database": {"sha256": file_sha256(path)}}
    (source / "manifest.json").write_text(json.dumps(manifest))

    def load(directory, verify_hash=True):
        return sqlite3.connect(directory / "catalog.sqlite3"), json.loads((directory / "manifest.json").read_text())

    monkeypatch.setattr(catalog, "load_catalog", load)
    index = tmp_path / "index"
    build_index(source, index, min_free_bytes=0)
    return source, index


def test_selected_distance_backend_matches_independent_dp_and_cutoffs():
    def oracle(a, b):
        previous = list(range(len(b) + 1))
        for i, x in enumerate(a, 1):
            current = [i]
            for j, y in enumerate(b, 1):
                current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (x != y)))
            previous = current
        return previous[-1]

    distance, identity = _distance_backend()
    assert identity["name"] in {"crisprmap.myers", "rapidfuzz.Levenshtein"}
    assert identity["version"]
    words = ["".join(letters) for length in range(4) for letters in product("AC", repeat=length)]
    for a in words:
        for b in words:
            expected = oracle(a, b)
            assert distance(a, b) == expected
            for cutoff in (0, 1, 2):
                observed = distance(a, b, cutoff)
                assert (observed <= cutoff) == (expected <= cutoff)
                if expected <= cutoff:
                    assert observed == expected


def test_streamed_exact_reverse_near_and_unsupported_results(tmp_path, monkeypatch):
    source, index = fixture(tmp_path, monkeypatch)
    destination = tmp_path / "results"
    queries = [{"id": "exact", "sequence": "AACCG"}, {"id": "reverse", "sequence": "CGGTT"},
               {"id": "near", "sequence": "AACCA"}, {"id": "n", "sequence": "NNNNN"}]
    manifest = search_catalog(iter(queries), source, index, destination, top_k=1, max_distance=1)
    results = [json.loads(line) for line in (destination / "results.jsonl").read_text().splitlines()]
    assert [row["status"] for row in results] == ["exact", "exact", "near", "unsupported"]
    assert results[1]["hits"][0]["orientation"] == "reverse"
    assert results[2]["best_node_count"] == 2
    assert results[2]["hits_truncated"] is True
    assert results[2]["best_occurrence_counts"]["bona_fide_repeat_count"] == 6
    assert results[2]["biological_family"] is None
    assert results[2]["structural_motif"] is None
    with sqlite3.connect(destination / "results.sqlite3") as connection:
        assert connection.execute("SELECT count(*) FROM matches WHERE query_id='near'").fetchone()[0] == 2
    assert manifest["counts"]["queries"] == 4
    assert manifest["status"] == "complete"
    for name, descriptor in manifest["files"].items():
        assert file_sha256(destination / name) == descriptor["sha256"]


def test_duplicate_identifiers_and_budget_fail_without_completed_artifact(tmp_path, monkeypatch):
    source, index = fixture(tmp_path, monkeypatch)
    destination = tmp_path / "results"
    with pytest.raises(InputError, match="Duplicate"):
        search_catalog([{"id": "q", "sequence": "AACCG"}] * 2, source, index, destination)
    assert not destination.exists()
    with pytest.raises(LimitError, match="query budget"):
        search_catalog([{"id": "a", "sequence": "AACCG"}, {"id": "b", "sequence": "AACCG"}],
                       source, index, destination, max_queries=1)
    assert not destination.exists()
    with pytest.raises(LimitError, match="candidate budget"):
        search_catalog([{"id": "q", "sequence": "AACCA"}], source, index, destination,
                       max_candidates=1, max_distance=6)
    assert not destination.exists()


def test_output_is_exclusive_and_no_hit_is_not_a_failed_search(tmp_path, monkeypatch):
    source, index = fixture(tmp_path, monkeypatch)
    destination = tmp_path / "results"
    query = [{"id": "q", "sequence": "AAAAA"}]
    report = search_catalog(query, source, index, destination, max_distance=0)
    assert report["counts"]["no_hit"] == 1
    before = (destination / "manifest.json").read_bytes()
    with pytest.raises(InputError):
        search_catalog(query, source, index, destination)
    assert (destination / "manifest.json").read_bytes() == before


def test_all_tied_cluster_groups_survive_hit_display_cap(tmp_path, monkeypatch):
    source, index = fixture(tmp_path, monkeypatch)
    clusters = tmp_path / "clusters"
    clusters.mkdir()
    with sqlite3.connect(clusters / "clusters.sqlite3") as connection:
        connection.execute("CREATE TABLE members(node_id INTEGER,threshold TEXT,cluster_id TEXT)")
        connection.executemany("INSERT INTO members VALUES (?,?,?)", [(1, "0.9", "cluster-a"), (2, "0.9", "cluster-b")])
    data = {"schema": {"name": "crisprmap.sequence-clusters", "version": "1.0.0"}, "status": "complete",
            "thresholds": {"0.9": {}},
            "catalog": {"manifest_sha256": file_sha256(source / "manifest.json"),
                        "database_sha256": file_sha256(source / "catalog.sqlite3")},
            "policy": {"index_manifest_sha256": file_sha256(index / "manifest.json")},
            "database": {"file": "clusters.sqlite3", "sha256": file_sha256(clusters / "clusters.sqlite3"),
                         "bytes": (clusters / "clusters.sqlite3").stat().st_size}}
    (clusters / "manifest.json").write_text(json.dumps(data))
    destination = tmp_path / "result"
    search_catalog([{"id": "q", "sequence": "AACCA"}], source, index, destination,
                   clusters_dir=clusters, top_k=1)
    result = json.loads((destination / "results.jsonl").read_text())
    assert len(result["hits"]) == 1
    assert result["operational_clusters"] == [
        {"threshold": "0.9", "cluster_id": "cluster-a", "best_node_count": 1},
        {"threshold": "0.9", "cluster_id": "cluster-b", "best_node_count": 1}]
    assert result["cluster_relation"] == "nearest_reference_cluster_evidence"


def test_exact_noncore_repeat_maps_all_qualifying_clusters_without_merging(tmp_path, monkeypatch):
    source, index = fixture(tmp_path, monkeypatch,
                            sequences=["CAAAAAAAAA", "AACAAAAAAA", "AAAAAAAAAA"], core_ids={1, 2})
    clusters = tmp_path / "clusters"
    clusters.mkdir()
    with sqlite3.connect(clusters / "clusters.sqlite3") as connection:
        connection.execute("CREATE TABLE members(node_id INTEGER,threshold TEXT,cluster_id TEXT)")
        connection.executemany("INSERT INTO members VALUES (?,?,?)", [(1, "0.9", "left"), (2, "0.9", "right")])
    data = {"schema": {"name": "crisprmap.sequence-clusters", "version": "1.0.0"}, "status": "complete",
            "thresholds": {"0.9": {}},
            "catalog": {"manifest_sha256": file_sha256(source / "manifest.json"),
                        "database_sha256": file_sha256(source / "catalog.sqlite3")},
            "policy": {"index_manifest_sha256": file_sha256(index / "manifest.json")},
            "database": {"file": "clusters.sqlite3", "sha256": file_sha256(clusters / "clusters.sqlite3"),
                         "bytes": (clusters / "clusters.sqlite3").stat().st_size}}
    (clusters / "manifest.json").write_text(json.dumps(data))
    original_sha = file_sha256(clusters / "clusters.sqlite3")
    destination = tmp_path / "result"
    search_catalog([{"id": "variant", "sequence": "AAAAAAAAAA"}], source, index, destination,
                   clusters_dir=clusters, top_k=1, max_distance=0)
    result = json.loads((destination / "results.jsonl").read_text())
    assert result["status"] == "exact"
    assert result["hits"][0]["node_id"] == 3
    assert result["operational_clusters"] == []
    assert result["cluster_mapping_status"] == "completed_with_candidates"
    assert result["operational_cluster_candidates"] == [
        {"threshold": "0.9", "cluster_id": "left", "best_distance": 1, "qualifying_core_node_count": 1},
        {"threshold": "0.9", "cluster_id": "right", "best_distance": 1, "qualifying_core_node_count": 1}]
    assert file_sha256(clusters / "clusters.sqlite3") == original_sha
    with sqlite3.connect(destination / "results.sqlite3") as connection:
        assert connection.execute("SELECT count(*) FROM cluster_candidates").fetchone()[0] == 2


def test_real_catalog_validation_integrates_index_clusters_and_search(tmp_path):
    from crisprmap.catalog import build_partition, load_catalog
    from crisprmap.catalog_index import validate_index
    from crisprmap.clusters import assemble_clusters, build_forest
    from test_catalog import raw_plan

    # Exercise the real schema, metadata/provenance checks, query-only catalogue
    # connection and published byte hashes; no load_catalog monkeypatch here.
    plan, profile, _ = raw_plan(tmp_path)
    source = tmp_path / "real-catalog"
    build_partition(plan, source, profile, min_free_bytes=0)
    connection, catalog_manifest = load_catalog(source)
    try:
        assert connection.execute("PRAGMA query_only").fetchone()[0] == 1
        assert catalog_manifest["counts"]["nodes"] == 1
    finally:
        connection.close()
    index = tmp_path / "real-index"
    build_index(source, index, min_free_bytes=0)
    validate_index(source, index, deep=True)
    forest = tmp_path / "real-forest"
    build_forest(source, index, forest, thresholds=["0.9"])
    clusters = tmp_path / "real-clusters"
    assemble_clusters(source, [forest], clusters)
    output = tmp_path / "real-results"
    search_catalog([{"id": "original", "sequence": "AACG"},
                    {"id": "reverse", "sequence": "CGTT"}],
                   source, index, output, clusters_dir=clusters)
    results = [json.loads(line) for line in (output / "results.jsonl").read_text().splitlines()]
    assert [row["status"] for row in results] == ["exact", "exact"]
    assert results[1]["hits"][0]["orientation"] == "reverse"
    assert results[0]["best_occurrence_counts"]["bona_fide_repeat_count"] == 2
    assert results[0]["cluster_mapping_status"] == "not_needed_exact_core_member"
    assert len(results[0]["operational_clusters"]) == 1
