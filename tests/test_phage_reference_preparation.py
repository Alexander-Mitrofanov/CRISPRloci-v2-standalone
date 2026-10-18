"""Offline operator acquisition tests with synthetic records and no network access."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sqlite3
import sys
import threading
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import prepare_phage_reference as builder


@pytest.fixture(autouse=True)
def enough_operator_capacity(monkeypatch):
    monkeypatch.setattr(builder.shutil, "disk_usage", lambda _: SimpleNamespace(free=100 * 1024**3))


def record(uid, accession, sequence="A" * 1000, **changes):
    values = {
        "accession-version": accession,
        "definition": "Synthetic fixture phage, complete genome",
        "organism": "Fixture phage",
        "taxonomy": "Viruses; fixture",
        "moltype": "DNA",
        "length": str(len(sequence)),
        "sequence": sequence.lower(),
        "update-date": "01-JAN-2026",
    }
    values.update(changes)
    node = ET.Element("GBSeq")
    for name, value in values.items():
        if value is not None:
            ET.SubElement(node, "GBSeq_" + name).text = value
    ids = ET.SubElement(node, "GBSeq_other-seqids")
    ET.SubElement(ids, "GBSeqid").text = "gi|" + uid
    table = ET.SubElement(node, "GBSeq_feature-table")
    feature = ET.SubElement(table, "GBFeature")
    ET.SubElement(feature, "GBFeature_key").text = "source"
    qualifiers = ET.SubElement(feature, "GBFeature_quals")
    for name, value in (("mol_type", "genomic DNA"), ("db_xref", "taxon:10239")):
        qualifier = ET.SubElement(qualifiers, "GBQualifier")
        ET.SubElement(qualifier, "GBQualifier_name").text = name
        ET.SubElement(qualifier, "GBQualifier_value").text = value
    return node


def document(records):
    root = ET.Element("GBSet")
    root.extend(records)
    return ET.tostring(root, encoding="utf-8", xml_declaration=True)


def snapshot(tmp_path, records, *, ids=None, host_ids=None):
    root = tmp_path / "snapshot"
    sources = root / "sources"
    sources.mkdir(parents=True)
    ids = ids or [builder.parse_record(row)["uid"] for row in records]
    host_ids = host_ids or ids[:1]
    for name, values in (("candidates.json", ids), ("host-qualified.json", host_ids)):
        builder.save_json(sources / name, {"count": str(len(values)), "idlist": values})
    metadata = {
        "schema_version": "crisprloci.ncbi-phage-snapshot/1.0.0",
        "query": builder.QUERY,
        "host_query": builder.HOST_QUERY,
        "retrieved_at": "2026-10-10T12:00:00+00:00",
        "acquisition_completed_at": "2026-10-10T13:00:00+00:00",
        "batch_size": builder.BATCH_SIZE,
        "candidate_count": len(ids),
        "search_sha256": {
            name: builder.sha256(sources / name)
            for name in ("candidates.json", "host-qualified.json")
        },
    }
    builder.save_json(root / "snapshot.json", metadata)
    for start in range(0, len(ids), builder.BATCH_SIZE):
        path = sources / f"batch-{start:06d}.xml.gz"
        payload = document(records[start:start + builder.BATCH_SIZE])
        path.write_bytes(gzip.compress(payload, mtime=0))
        builder.save_json(path.with_suffix(".json"), {
            "ids": ids[start:start + builder.BATCH_SIZE],
            "sha256": builder.sha256(path),
            "content_sha256": hashlib.sha256(payload).hexdigest(),
            "content_bytes": len(payload),
        })
    return root


def test_selection_accounts_for_all_candidates_and_keeps_preferred_aliases(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "BATCH_SIZE", 2)
    records = [
        record("1", "AB000001.1", definition="Fixture phage, complete genome"),
        record("2", "NC_000002.1", definition="Fixture phage, complete genome"),
        record("3", "NC_000001.1", definition="Fixture phage, complete genome"),
        record("4", "AB000004.1", "C" * 1000, definition="Fixture phage, complete genome"),
        record("5", "AB000005.1", "N" * 1000, definition="Fixture phage, complete genome"),
    ]
    root = snapshot(tmp_path, records)
    receipt = builder.select(root)
    assert receipt["records"] == 5
    assert receipt["statuses"] == {"duplicate": 2, "excluded": 1, "selected": 2}
    assert receipt["exclusions"] == {"no_called_bases": 1}
    assert receipt["genome_count"] == 2 and receipt["total_bases"] == 2000
    output = root / "selected"
    with sqlite3.connect(output / "selection.sqlite3") as db:
        aliases = db.execute("SELECT accession,representative FROM records WHERE status='duplicate' ORDER BY accession").fetchall()
        assert aliases == [("AB000001.1", "NC_000001.1"), ("NC_000002.1", "NC_000001.1")]
        assert db.execute("SELECT count(DISTINCT uid) FROM records").fetchone()[0] == 5
        metadata = json.loads(db.execute("SELECT metadata_json FROM records WHERE uid='1'").fetchone()[0])
        assert metadata["inclusion_evidence"] == ["ncbi_prokaryotic_vhost_filter", "explicit_phage_title"]
    files = (output / "genomes.txt").read_text().splitlines()
    assert [Path(name).name for name in files] == ["NC_000001.1.fasta", "AB000004.1.fasta"]
    assert all(Path(name).is_file() for name in files)
    assert len(list((output / "genomes").glob("*.fasta"))) == 2
    assert receipt["selection_database_sha256"] == builder.sha256(output / "selection.sqlite3")
    source = builder.read_json(output / "source.json")
    assert set(source) == {"collection", "release", "retrieved_at", "selection_policy"}
    assert all(isinstance(value, str) and 1 <= len(value) <= 1000 and all(ord(c) >= 32 for c in value) for value in source.values())
    assert builder.sha256(output / "selection.json") in source["selection_policy"]
    previous = {path.name: builder.sha256(path) for path in output.iterdir() if path.is_file()}
    with pytest.raises(ValueError, match="never overwrite"):
        builder.select(root)
    assert previous == {path.name: builder.sha256(path) for path in output.iterdir() if path.is_file()}


@pytest.mark.parametrize("change,reason", [
    ({"taxonomy": "Bacteria; Bacillota"}, "not_viral_taxonomy"),
    ({"taxonomy": None}, "not_viral_taxonomy"),
    ({"definition": "Fixture phage, partial genome"}, "not_complete_genome_title"),
    ({"definition": "Fixture virus, complete genome"}, "no_phage_or_prokaryotic_vhost_evidence"),
    ({"moltype": "RNA"}, "not_single_genomic_DNA_record"),
    ({"moltype": None}, "not_single_genomic_DNA_record"),
    ({"definition": "Fixture phage segment, complete genome"}, "segmented_or_multiple_source_features"),
    ({"definition": "Fixture phage satellite, complete genome"}, "satellite_vector_partial_or_synthetic"),
])
def test_missing_or_ineligible_metadata_is_explicitly_excluded(tmp_path, change, reason):
    root = snapshot(tmp_path, [
        record("1", "NC_000001.1", definition="Fixture phage, complete genome"),
        record("2", "AB000002.1", **change),
    ])
    receipt = builder.select(root)
    assert receipt["records"] == 2 and receipt["statuses"] == {"excluded": 1, "selected": 1}
    assert receipt["exclusions"] == {reason: 1}


@pytest.mark.parametrize("records,ids", [
    ([record("1", "AB000001.1")], ["1", "2"]),
    ([record("3", "AB000001.1")], ["1"]),
    ([record("1", "AB000001.1"), record("1", "AB000002.1")], ["1", "2"]),
])
def test_selection_rejects_missing_wrong_or_duplicate_uid_without_publish(tmp_path, records, ids):
    root = snapshot(tmp_path, records, ids=ids)
    with pytest.raises(ValueError, match="coverage mismatch|identity mismatch"):
        builder.select(root)
    assert not (root / "selected").exists()
    assert not list(root.glob(".selection-*"))


def test_truncated_xml_cannot_publish_even_with_matching_file_receipt(tmp_path):
    root = snapshot(tmp_path, [record("1", "AB000001.1")])
    path = root / "sources/batch-000000.xml.gz"
    path.write_bytes(gzip.compress(document([record("1", "AB000001.1")])[:-10]))
    receipt_path = path.with_suffix(".json")
    receipt = builder.read_json(receipt_path)
    receipt["sha256"] = builder.sha256(path)
    builder.save_json(receipt_path, receipt)
    with pytest.raises(ET.ParseError):
        builder.select(root)
    assert not (root / "selected").exists()
    assert not list(root.glob(".selection-*"))


def test_frozen_metadata_and_compressed_batches_must_match_receipts(tmp_path):
    root = snapshot(tmp_path, [record("1", "AB000001.1")])
    path = root / "sources/batch-000000.xml.gz"
    path.write_bytes(path.read_bytes() + b"changed")
    with pytest.raises(ValueError, match="batch changed"):
        builder.select(root)
    (root / "sources/candidates.json").write_text("{}")
    with pytest.raises(ValueError, match="search receipt changed"):
        builder.select(root)
    assert not (root / "selected").exists()


def test_acquisition_resumes_receipted_batches_without_network(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "BATCH_SIZE", 2)
    monkeypatch.setattr(builder.shutil, "disk_usage", lambda _: SimpleNamespace(free=100 * 1024**3))
    monkeypatch.setattr(builder, "search", lambda query: {"count": "3" if query == builder.QUERY else "1", "idlist": ["1", "2", "3"] if query == builder.QUERY else ["1"]})
    calls = []

    def fetch(endpoint, parameters):
        assert parameters["rettype"] == "gbwithparts" and parameters["retmode"] == "xml"
        calls.append(parameters["id"])
        return io.BytesIO(document([record(uid, f"AB{int(uid):06d}.1") for uid in parameters["id"].split(",")]))

    monkeypatch.setattr(builder, "request", fetch)
    root = tmp_path / "acquire"
    builder.acquire(root)
    assert calls == ["1,2", "3"]
    assert builder.read_json(root / "sources/batch-000000.xml.json")["rettype"] == "gbwithparts"
    monkeypatch.setattr(builder, "search", lambda _: pytest.fail("Resume repeated frozen search"))
    monkeypatch.setattr(builder, "request", lambda *_: pytest.fail("Resume downloaded a completed batch"))
    builder.acquire(root)
    path = root / "sources/batch-000000.xml.json"
    receipt = builder.read_json(path)
    receipt["sha256"] = "0" * 64
    builder.save_json(path, receipt)
    with pytest.raises(ValueError, match="batch changed"):
        builder.acquire(root)


@pytest.mark.parametrize("payload", [
    document([record("2", "AB000002.1")]),
    b"<html>Service unavailable</html>",
    document([record("1", "AB000001.1")])[:-10],
])
def test_invalid_download_is_never_receipted_or_marked_complete(tmp_path, monkeypatch, payload):
    monkeypatch.setattr(builder.shutil, "disk_usage", lambda _: SimpleNamespace(free=100 * 1024**3))
    monkeypatch.setattr(builder, "search", lambda _: {"count": "1", "idlist": ["1"]})
    monkeypatch.setattr(builder, "request", lambda *_: io.BytesIO(payload))
    root = tmp_path / "invalid"
    with pytest.raises((ValueError, ET.ParseError)):
        builder.acquire(root)
    assert not (root / "sources/batch-000000.xml.json").exists()
    assert "acquisition_completed_at" not in builder.read_json(root / "snapshot.json")


def test_cached_xml_and_json_reads_remain_bounded(tmp_path, monkeypatch):
    path = tmp_path / "large.xml.gz"
    path.write_bytes(gzip.compress(document([record("1", "AB000001.1")])) )
    monkeypatch.setattr(builder, "MAX_BATCH_BYTES", 100)
    with pytest.raises(ValueError, match="bounded decompressed size"):
        list(builder.iter_records(path))
    metadata = tmp_path / "oversize.json"
    metadata.write_text('{"unused":"' + "x" * 100 + '"}')
    monkeypatch.setattr(builder, "MAX_METADATA_BYTES", 50)
    with pytest.raises(ValueError, match="metadata exceeds"):
        builder.read_json(metadata)


def test_selection_stops_before_writing_into_disk_reserve_and_cleans_staging(tmp_path, monkeypatch):
    root = snapshot(tmp_path, [record("1", "AB000001.1", definition="Fixture phage, complete genome")])
    monkeypatch.setattr(builder.shutil, "disk_usage", lambda _: SimpleNamespace(free=builder.MIN_FREE_BYTES + 100))
    with pytest.raises(ValueError, match="deployment reserve"):
        builder.select(root)
    assert not (root / "selected").exists()
    assert not list(root.glob(".selection-*"))


@pytest.mark.parametrize("result", [
    {"count": "2", "idlist": ["1"]},
    {"count": "2", "idlist": ["1", "1"]},
    {"count": "1", "idlist": [1]},
    {"count": "1", "idlist": ["wrong"]},
])
def test_incomplete_or_malformed_search_cannot_define_snapshot(result):
    with pytest.raises(ValueError, match="every unique candidate"):
        builder.search_ids(result)


def test_three_download_workers_cover_every_frozen_uid_and_resume(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "BATCH_SIZE", 2)
    ids = [str(value) for value in range(1, 7)]
    monkeypatch.setattr(builder, "search", lambda query: {
        "count": "6" if query == builder.QUERY else "1",
        "idlist": ids if query == builder.QUERY else ids[:1],
    })
    barrier = threading.Barrier(3)
    calls = []

    def fetch(endpoint, parameters):
        batch = parameters["id"].split(",")
        calls.append(batch)
        barrier.wait(timeout=5)
        return io.BytesIO(document([record(uid, f"AB{int(uid):06d}.1") for uid in batch]))

    monkeypatch.setattr(builder, "request", fetch)
    root = tmp_path / "parallel"
    result = builder.acquire(root, download_workers=3)
    assert result["acquisition_completed_at"]
    assert sorted(uid for batch in calls for uid in batch) == ids
    receipts = [builder.read_json(path) for path in sorted((root / "sources").glob("batch-*.json"))]
    assert len(receipts) == 3
    assert [uid for receipt in receipts for uid in receipt["ids"]] == ids
    for path in (root / "sources").glob("batch-*.xml.gz"):
        assert builder.read_json(path.with_suffix(".json"))["sha256"] == builder.sha256(path)
    monkeypatch.setattr(builder, "search", lambda _: pytest.fail("Frozen queries must not repeat"))
    monkeypatch.setattr(builder, "request", lambda *_: pytest.fail("Completed batches must not download"))
    builder.acquire(root, download_workers=2)


def test_parallel_failure_keeps_completed_checkpoints_and_resumes_missing_batch(tmp_path, monkeypatch):
    monkeypatch.setattr(builder, "BATCH_SIZE", 2)
    ids = [str(value) for value in range(1, 7)]
    monkeypatch.setattr(builder, "search", lambda query: {
        "count": "6" if query == builder.QUERY else "1",
        "idlist": ids if query == builder.QUERY else ids[:1],
    })
    checkpoint = threading.Event()
    save = builder.save_json

    def save_and_signal(path, value):
        save(path, value)
        if path.name == "batch-000000.xml.json":
            checkpoint.set()

    monkeypatch.setattr(builder, "save_json", save_and_signal)

    def fetch(endpoint, parameters):
        batch = parameters["id"].split(",")
        if batch[0] == "3":
            assert checkpoint.wait(timeout=5), "A completed batch must precede the injected failure"
            raise TimeoutError("Synthetic interrupted download")
        return io.BytesIO(document([record(uid, f"AB{int(uid):06d}.1") for uid in batch]))

    monkeypatch.setattr(builder, "request", fetch)
    root = tmp_path / "partial"
    with pytest.raises(TimeoutError, match="Synthetic interrupted"):
        builder.acquire(root, download_workers=3)
    assert "acquisition_completed_at" not in builder.read_json(root / "snapshot.json")
    checkpoint_path = root / "sources/batch-000000.xml.json"
    before = builder.sha256(checkpoint_path)
    assert not (root / "sources/batch-000002.xml.json").exists()
    calls = []

    def retry(endpoint, parameters):
        batch = parameters["id"].split(",")
        assert batch[0] != "1", "Verified checkpoint must not be fetched again"
        calls.append(batch)
        return io.BytesIO(document([record(uid, f"AB{int(uid):06d}.1") for uid in batch]))

    monkeypatch.setattr(builder, "request", retry)
    monkeypatch.setattr(builder, "search", lambda _: pytest.fail("Resume must preserve the search snapshot"))
    assert builder.acquire(root, download_workers=3)["acquisition_completed_at"]
    assert ["3", "4"] in calls
    assert builder.sha256(checkpoint_path) == before


def test_request_rate_budget_is_shared_across_download_threads(monkeypatch):
    clock = {"now": 100.0}
    sleeps = []

    def sleep(delay):
        sleeps.append(delay)
        clock["now"] += delay

    monkeypatch.setattr(builder, "time", SimpleNamespace(monotonic=lambda: clock["now"], sleep=sleep))
    monkeypatch.setattr(builder, "_next_request_at", 0.0)
    starts = []
    barrier = threading.Barrier(3)

    def network(request, timeout):
        starts.append(clock["now"])
        response = io.BytesIO(b"response")
        response.url = builder.BASE + "efetch.fcgi"
        return response

    monkeypatch.setattr(builder.urllib.request, "urlopen", network)

    def worker(_):
        barrier.wait(timeout=5)
        with builder.request("efetch.fcgi", {"id": "1"}) as response:
            assert response.read() == b"response"

    with ThreadPoolExecutor(max_workers=3) as executor:
        list(executor.map(worker, range(9)))
    assert len(starts) == 9 and len(sleeps) == 8
    assert all(later - earlier >= builder.REQUEST_INTERVAL - 1e-8 for earlier, later in zip(sorted(starts), sorted(starts)[1:]))


@pytest.mark.parametrize("workers", [0, 4, -1, True, 1.5])
def test_download_worker_count_is_bounded_before_any_acquisition(tmp_path, workers):
    with pytest.raises(ValueError, match="between one and three"):
        builder.acquire(tmp_path / "absent", download_workers=workers)
    assert not (tmp_path / "absent").exists()


def test_nonviral_contig_without_inline_sequence_is_accounted_as_excluded(tmp_path):
    contig = record(
        "2", "NZ_CP011102.1", sequence="",
        definition="Fixture bacterial chromosome, complete genome",
        taxonomy="Bacteria; Bacillota", length="3406292", contig="join(CP011102.1:1..3406292)",
    )
    parsed = builder.parse_record(contig)
    assert parsed["length"] == 3406292 and parsed["sequence_available"] is False
    assert parsed["sequence_sha256"] is None and parsed["sequence"] == ""
    root = snapshot(tmp_path, [
        record("1", "NC_000001.1", definition="Fixture phage, complete genome"), contig,
    ], host_ids=["1", "2"])
    receipt = builder.select(root)
    assert receipt["records"] == 2
    assert receipt["statuses"] == {"excluded": 1, "selected": 1}
    assert receipt["exclusions"] == {"not_viral_taxonomy": 1}
    with sqlite3.connect(root / "selected/selection.sqlite3") as db:
        missing = db.execute("SELECT sequence_sha256,metadata_json FROM records WHERE uid='2'").fetchone()
    assert missing[0] is None and json.loads(missing[1])["length"] == 3406292


@pytest.mark.parametrize("sequence,changes", [
    ("", {"length": "1000"}),
    ("A" * 999, {"length": "1000", "contig": "join(AB000001.1:1..1000)"}),
])
def test_missing_without_contig_and_nonempty_length_mismatches_still_fail(sequence, changes):
    with pytest.raises(ValueError, match="sequence length"):
        builder.parse_record(record("1", "AB000001.1", sequence=sequence, **changes))


def test_viral_contig_missing_sequence_is_not_silently_excluded_or_receipted(tmp_path, monkeypatch):
    contig = record("1", "AB000001.1", sequence="", length="1000", contig="join(AB000002.1:1..1000)")
    with pytest.raises(ValueError, match="requires expanded NCBI sequence"):
        builder.exclusion(builder.parse_record(contig), {"1"})
    monkeypatch.setattr(builder, "search", lambda _: {"count": "1", "idlist": ["1"]})
    monkeypatch.setattr(builder, "request", lambda *_: io.BytesIO(document([contig])))
    root = tmp_path / "missing-viral"
    with pytest.raises(ValueError, match="requires expanded NCBI sequence"):
        builder.acquire(root)
    assert not (root / "sources/batch-000000.xml.json").exists()
    assert "acquisition_completed_at" not in builder.read_json(root / "snapshot.json")
