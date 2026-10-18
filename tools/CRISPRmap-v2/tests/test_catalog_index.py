import json
from itertools import product
import random
import sqlite3
from types import SimpleNamespace

import pytest

from crisprmap import catalog
from crisprmap import catalog_index
from crisprmap.catalog_index import CandidateIndex, build_index, file_sha256, validate_index
from crisprmap.core import canonical_sequence, reverse_complement
from crisprmap.io import InputError, LimitError, digest


def make_catalog(tmp_path, monkeypatch, sequences):
    path = tmp_path / "catalog"
    path.mkdir()
    database = path / "catalog.sqlite3"
    with sqlite3.connect(database) as connection:
        connection.execute("""CREATE TABLE nodes(node_id INTEGER PRIMARY KEY,sequence TEXT UNIQUE,
            sequence_sha256 TEXT UNIQUE,length INTEGER,eligible INTEGER,bona_fide_consensus_count INTEGER,
            possible_consensus_count INTEGER,bona_fide_repeat_count INTEGER,possible_repeat_count INTEGER)""")
        for index, sequence in enumerate(sorted({canonical_sequence(s) for s in sequences}), 1):
            connection.execute("INSERT INTO nodes VALUES (?,?,?,?,?,?,?,?,?)",
                               (index, sequence, digest(sequence), len(sequence), 1, index % 2, 1, 4, 2))
    manifest = {"schema": {"name": "crisprmap.catalog", "version": "1.0.0"}, "status": "complete",
                "database": {"file": "catalog.sqlite3", "sha256": file_sha256(database),
                             "bytes": database.stat().st_size}}
    (path / "manifest.json").write_text(json.dumps(manifest))

    def load(directory, verify_hash=True):
        data = json.loads((directory / "manifest.json").read_text())
        if verify_hash and file_sha256(directory / "catalog.sqlite3") != data["database"]["sha256"]:
            raise InputError("test catalog checksum failure")
        connection = sqlite3.connect((directory / "catalog.sqlite3").as_uri() + "?mode=ro", uri=True)
        return connection, data

    monkeypatch.setattr(catalog, "load_catalog", load)
    return path


def oracle(a, b):
    previous = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (x != y)))
        previous = current
    return previous[-1]


def test_partition_candidates_are_complete_against_independent_oracle(tmp_path, monkeypatch):
    rng = random.Random(20261001)
    sequences = ["".join(rng.choices("ACGT", k=rng.randrange(1, 61))) for _ in range(70)]
    source = make_catalog(tmp_path, monkeypatch, sequences)
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    validate_index(source, destination, deep=True)
    queries = sequences[:20] + [reverse_complement(s) for s in sequences[20:30]]
    queries += [s[:len(s) // 2] + "A" + s[len(s) // 2:] for s in sequences[:20]]
    queries += [s[:-1] if len(s) > 1 else "T" for s in sequences[:10]]
    with CandidateIndex(source, destination) as reader:
        nodes = list(reader.catalog.execute("SELECT node_id,sequence FROM nodes"))
        for distance in (0, 1, 3, 6):
            for query in queries:
                expected = {node_id for node_id, sequence in nodes
                            if min(oracle(query, sequence), oracle(query, reverse_complement(sequence))) <= distance}
                candidates = {row["node_id"] for row in reader.iter_candidates(query, max_distance=distance)}
                assert expected <= candidates


def test_tight_offsets_exhaustively_preserve_short_edit_neighbors(tmp_path, monkeypatch):
    sequences = ["".join(bases) for length in range(1, 7) for bases in product("AC", repeat=length)]
    source = make_catalog(tmp_path, monkeypatch, sequences)
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with CandidateIndex(source, destination) as reader:
        nodes = list(reader.catalog.execute("SELECT node_id,sequence FROM nodes"))
        for query in sequences:
            distances = {node: min(oracle(query, target), oracle(query, reverse_complement(target)))
                         for node, target in nodes}
            for distance in (0, 1, 2, 3, 4):
                expected = {node for node, observed in distances.items() if observed <= distance}
                # Alternating input orientation exercises both sides of the
                # reverse-complement equivalence without changing the oracle.
                oriented = reverse_complement(query) if distance % 2 else query
                rows = list(reader.iter_candidates(oriented, max_distance=distance))
                assert expected <= {row["node_id"] for row in rows}
                assert [row["node_id"] for row in rows] == sorted(row["node_id"] for row in rows)


def test_tight_offsets_preserve_boundary_indels_under_fixed_and_relative_policies(tmp_path, monkeypatch):
    rng = random.Random(991307)
    references = [canonical_sequence("".join(rng.choices("ACGT", k=length))) for length in (30, 59)]
    source = make_catalog(tmp_path, monkeypatch, references)
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    queries = set()
    for reference in references:
        # Ends and positions on either side of balanced segment boundaries.
        # Odd indel widths exercise rounding with both signs of length change.
        positions = {0, len(reference), len(reference) // 2}
        for count in (4, 7, 11):
            boundary = len(reference) // count
            positions.update((boundary - 1, boundary, boundary + 1))
        for position in positions:
            for width in (1, 3, 5):
                queries.add(reference[:position] + "GATAC"[:width] + reference[position:])
                if position + width <= len(reference):
                    queries.add(reference[:position] + reference[position + width:])
            if position < len(reference):
                # Compensating insertion/deletion has equal total length but
                # can shift an intact middle segment in either direction.
                queries.add(reference[:position] + "G" + reference[position:-1])
                queries.add(reference[:position] + reference[position + 1:] + "T")
    with CandidateIndex(source, destination) as reader:
        nodes = list(reader.catalog.execute("SELECT node_id,sequence FROM nodes"))
        for number, query in enumerate(sorted(queries)):
            observed = {node: min(oracle(query, target), oracle(query, reverse_complement(target)))
                        for node, target in nodes}
            oriented = reverse_complement(query) if number % 2 else query
            for distance in (1, 3, 5):
                candidates = {row["node_id"] for row in reader.iter_candidates(oriented, max_distance=distance)}
                assert {node for node, value in observed.items() if value <= distance} <= candidates
            for percent in (85, 90, 95):
                candidates = {row["node_id"] for row in reader.iter_candidates(
                    oriented, max_distance=None, min_similarity=percent / 100)}
                expected = {node for node, target in nodes
                            if 100 * observed[node] <= (100 - percent) * max(len(query), len(target))}
                assert expected <= candidates


@pytest.mark.parametrize("difference", [-3, -2, -1, 0, 1, 2, 3])
def test_tight_offset_keys_match_independent_prefix_suffix_bound(tmp_path, monkeypatch, difference):
    rng = random.Random(7283)
    reference = "".join(rng.choices("ACGT", k=30))
    query = "".join(rng.choices("ACGT", k=30 + difference))
    source = make_catalog(tmp_path, monkeypatch, [reference])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with CandidateIndex(source, destination) as reader:
        list(reader.iter_candidates(query, max_distance=3))
        expected, previous = set(), set()
        k = min(value for value in reader.profiles[30] if value >= 3)
        for orientation in (query, reverse_complement(query)):
            for ordinal, start, size in reader.profiles[30][k]:
                for offset in range(len(query) - size + 1):
                    item = (30, k, ordinal, orientation[offset:offset + size])
                    prefix_edits = abs(offset - start)
                    suffix_edits = abs((len(query) - offset - size) - (30 - start - size))
                    if prefix_edits + suffix_edits <= 3:
                        expected.add(item)
                    if prefix_edits <= 3:
                        previous.add(item)
        actual = {tuple(row) for row in reader.connection.execute("SELECT * FROM temp.query_keys")}
        assert actual == expected
        assert actual <= previous
        if difference == 0:
            assert len(actual) < len(previous)


@pytest.mark.parametrize("similarity", [0.85, 0.90, 0.95, 1.0])
def test_relative_graph_candidates_cover_length_boundaries(tmp_path, monkeypatch, similarity):
    sequences = ["A" * length for length in (1, 2, 6, 9, 17, 18, 19, 29, 30, 59, 60, 61, 100, 199, 200)]
    sequences += ["ACGT" * 10, "ACGT" * 12 + "AC", "A" * 29 + "C", "A" * 59 + "C"]
    source = make_catalog(tmp_path, monkeypatch, sequences)
    destination = tmp_path / "index"
    build_index(source, destination, min_similarity=0.85, min_free_bytes=0)
    numerator = round((1 - similarity) * 100)
    with CandidateIndex(source, destination) as reader:
        nodes = list(reader.catalog.execute("SELECT node_id,sequence FROM nodes"))
        for query in sequences:
            expected = {node_id for node_id, sequence in nodes
                        if 100 * min(oracle(query, sequence), oracle(query, reverse_complement(sequence)))
                        <= numerator * max(len(query), len(sequence))}
            candidates = {row["node_id"] for row in reader.iter_candidates(
                query, max_distance=None, min_similarity=similarity)}
            assert expected <= candidates


def test_core_and_node_order_filters_are_complete(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["AAAA", "AAAC", "AACC", "ACCC", "CCCC"])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with CandidateIndex(source, destination) as reader:
        rows = list(reader.iter_candidates("AAAA", max_distance=6, core_only=True, min_node_id=1))
        assert [row["node_id"] for row in rows] == [3, 5]
        assert all(row["bona_fide_consensus_count"] > 0 for row in rows)


def test_sparse_candidates_lookup_nodes_without_scanning_catalogue(tmp_path, monkeypatch):
    # Every ten-base partition is unique. One exact query therefore nominates
    # one node in a much larger catalogue, exposing accidental catalogue scans.
    codes = ["".join("ACGT"[(value >> shift) & 3] for shift in range(0, 12, 2))
             for value in range(512)]
    sequences = ["A" + code + "CCC" + code + "GCTA" for code in codes]
    source = make_catalog(tmp_path, monkeypatch, sequences)
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with CandidateIndex(source, destination) as reader:
        statements = []
        reader.connection.set_trace_callback(statements.append)
        rows = list(reader.iter_candidates(sequences[0], max_distance=0))
        reader.connection.set_trace_callback(None)
        assert [row["sequence"] for row in rows] == [sequences[0]]
        assert reader.last_stats["candidates"] == 1
        # Inspect the actual SELECT executed by the production iterator, not a
        # test-side substitute. Candidate enumeration should drive PK lookups;
        # table-size growth must not turn every sparse query into a full scan.
        selects = [statement for statement in statements if "FROM temp.candidates c" in statement]
        assert len(selects) == 1
        details = [row[3] for row in reader.connection.execute("EXPLAIN QUERY PLAN " + selects[0])]
        assert any(detail.startswith("SCAN c") for detail in details)
        assert any(detail.startswith("SEARCH n USING INTEGER PRIMARY KEY") for detail in details)
        assert not any(detail.startswith("SCAN n") for detail in details)


def test_core_posting_lookup_preserves_all_matches_and_exact_budgets(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["A" * 19 + base for base in "ACG"])
    destination = tmp_path / "index"
    manifest = build_index(source, destination, min_free_bytes=0)
    assert manifest["counts"]["core_nodes"] == 2
    validate_index(source, destination, deep=True)
    with CandidateIndex(source, destination) as reader:
        # k=1 partitions have two ten-base segments. The exact core node
        # contributes two postings and each alternative contributes one.
        rows = list(reader.iter_candidates("A" * 20, max_distance=0, core_only=True, max_postings=3))
        assert [row["node_id"] for row in rows] == [1, 3]
        assert reader.last_stats["postings"] == 3
        with pytest.raises(LimitError, match="posting budget"):
            list(reader.iter_candidates("A" * 20, max_distance=0, max_postings=3))
        rows = list(reader.iter_candidates("A" * 20, max_distance=0, core_only=True,
                                           min_node_id=1, max_postings=1))
        assert [row["node_id"] for row in rows] == [3]
        assert reader.last_stats["postings"] == 1


def test_deep_validation_rejects_rehashed_incorrect_core_membership(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["A" * 20, "A" * 19 + "C"])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with sqlite3.connect(destination / "index.sqlite3") as connection:
        connection.execute("UPDATE postings SET is_core=1 WHERE node_id=2")
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["database"]["sha256"] = file_sha256(destination / "index.sqlite3")
    manifest["database"]["bytes"] = (destination / "index.sqlite3").stat().st_size
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(InputError, match="inventory"):
        validate_index(source, destination, deep=True)


def test_previous_index_schema_requires_explicit_rebuild(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["A" * 20])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    manifest_path = destination / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["schema"]["version"] = "1.0.0"
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(InputError, match="Unsupported candidate-index schema"):
        CandidateIndex(source, destination)


def test_candidate_caps_and_unsupported_policy_fail_explicitly(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["A", "C", "AA", "AC", "CC"])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with CandidateIndex(source, destination) as reader:
        with pytest.raises(LimitError, match="candidate budget"):
            list(reader.iter_candidates("AAAA", max_distance=6, max_candidates=1))
        with pytest.raises(LimitError, match="posting budget"):
            list(reader.iter_candidates("AAAA", max_distance=6, max_postings=1))
        with pytest.raises(InputError, match="below"):
            list(reader.iter_candidates("AAAA", max_distance=None, min_similarity=0.8))
        with pytest.raises(InputError):
            list(reader.iter_candidates("ANAA"))


def test_bulk_candidates_count_duplicate_postings_and_fallback_at_exact_caps(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["A" * length for length in (1, 2, 6, 7, 8, 9)])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with CandidateIndex(source, destination) as reader:
        # Lengths1,2,6 each nominate one fallback node; lengths7,8,9
        # each supply seven matching segments under the k=6 profile.
        rows = list(reader.iter_candidates("AAAAA", max_distance=6, max_postings=24, max_candidates=6))
        assert [row["node_id"] for row in rows] == list(range(1, 7))
        assert reader.last_stats["postings"] == 24
        assert reader.last_stats["candidates"] == 6
        assert reader.last_stats["fallback_length_count"] == 3
        with pytest.raises(LimitError, match="posting budget"):
            list(reader.iter_candidates("AAAAA", max_distance=6, max_postings=23))
        with pytest.raises(LimitError, match="candidate budget"):
            list(reader.iter_candidates("AAAAA", max_distance=6, max_candidates=5))
        rows = list(reader.iter_candidates("AAAAA", max_distance=6, core_only=True, min_node_id=1,
                                           max_postings=8, max_candidates=2))
        assert [row["node_id"] for row in rows] == [3, 5]
        assert reader.last_stats["postings"] == 8


def test_index_hash_and_rehashed_missing_postings_are_detected(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["AACCGG", "AAAAAA"])
    destination = tmp_path / "index"
    build_index(source, destination, min_free_bytes=0)
    with sqlite3.connect(destination / "index.sqlite3") as connection:
        connection.execute("DELETE FROM postings WHERE node_id=1")
    with pytest.raises(InputError, match="checksum"):
        CandidateIndex(source, destination)
    manifest = json.loads((destination / "manifest.json").read_text())
    manifest["database"]["sha256"] = file_sha256(destination / "index.sqlite3")
    manifest["database"]["bytes"] = (destination / "index.sqlite3").stat().st_size
    with sqlite3.connect(destination / "index.sqlite3") as connection:
        manifest["counts"]["postings"] = connection.execute("SELECT count(*) FROM postings").fetchone()[0]
    (destination / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(InputError, match="inventory"):
        validate_index(source, destination, deep=True)


def test_index_build_failure_is_atomic_and_preserves_existing_output(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["AAAA", "AAAC"])
    destination = tmp_path / "index"
    with pytest.raises(LimitError, match="node budget"):
        build_index(source, destination, max_nodes=1, min_free_bytes=0)
    assert not destination.exists()
    build_index(source, destination, min_free_bytes=0)
    before = (destination / "manifest.json").read_bytes()
    with pytest.raises(InputError):
        build_index(source, destination, min_free_bytes=0)
    assert (destination / "manifest.json").read_bytes() == before


def test_index_free_space_reserve_stops_before_start_and_cleans_partial_output(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["AAAA", "AAAC"])
    destination = tmp_path / "index"
    monkeypatch.setattr(catalog_index.shutil, "disk_usage", lambda _: SimpleNamespace(free=0))
    with pytest.raises(LimitError, match="free-space reserve"):
        build_index(source, destination, min_free_bytes=1024)
    assert not destination.exists()
    calls = 0

    def exhaust_after_start(_):
        nonlocal calls
        calls += 1
        return SimpleNamespace(free=2048 if calls == 1 else 0)

    monkeypatch.setattr(catalog_index.shutil, "disk_usage", exhaust_after_start)
    with pytest.raises(LimitError, match="free-space reserve"):
        build_index(source, destination, min_free_bytes=1024)
    assert not destination.exists()



def test_default_index_reserve_remains_two_gib_and_rejects_before_output(tmp_path, monkeypatch):
    source = make_catalog(tmp_path, monkeypatch, ["AAAA"])
    destination = tmp_path / "index"
    monkeypatch.setattr(catalog_index.shutil, "disk_usage", lambda _: SimpleNamespace(free=2 * 1024**3 - 1))
    with pytest.raises(LimitError, match="free-space reserve"):
        build_index(source, destination)
    assert not destination.exists()
