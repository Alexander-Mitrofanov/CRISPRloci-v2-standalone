from __future__ import annotations

import gzip
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import catalog_source_inventory as catalog


def write_json(path, value, compressed=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value).encode()
    path.write_bytes(gzip.compress(data, mtime=0) if compressed else data)
    return path.read_bytes()


def test_inventory_includes_unpublished_and_unavailable_outcomes(tmp_path):
    public = tmp_path / "part-01/public"
    bundle = public / "shards/shard-000001/publication/partial-results.tar.gz"
    bundle.parent.mkdir(parents=True)
    bundle.write_bytes(b"verified later by fetcher")
    first = {"accession": "GCA_000000001.1", "status": "complete", "sha256": "a" * 64,
             "file": "assemblies/GCA_000000001.1.json.gz"}
    collection = {"analysis_profile_sha256": catalog.PROFILE, "entries": [first]}
    collection_path = tmp_path / "part-01/shards/shard-000001/collection_index.json"
    raw = write_json(collection_path, collection)
    write_json(public / "index.json", {"publications": [{
        "analysis_profile_sha256": catalog.PROFILE, "shard_name": "shard-000001",
        "collection_index_sha256": catalog.digest(raw),
        "bundle": {"download_path": str(bundle.relative_to(public)),
                   "bytes": bundle.stat().st_size, "sha256": catalog.digest(bundle.read_bytes())},
    }]})
    # Published raw files must not be inventoried again, even when retained.
    write_json(collection_path.parent / first["file"], {"ignored": True}, True)
    write_json(tmp_path / "unfinished/assemblies/GCA_000000002.1.json.gz", {
        "assembly": {"accession": "GCA_000000002.1"}, "status": "complete",
        "provenance": {"analysis_profile_sha256": catalog.PROFILE},
    }, True)
    write_json(tmp_path / "unfinished/failures/GCA_000000003.1.json.gz", {
        "assembly": {"accession": "GCA_000000003.1"}, "status": "failed",
        "error": {"type": "TimeoutError", "message": "time limit exceeded"},
        "provenance": {"analysis_profile_sha256": catalog.PROFILE},
    }, True)
    write_json(tmp_path / "tail/job-receipt.json", {
        "unavailable": {"GCA_000000004.1": {"reason": "NCBI input unavailable"}},
    })
    result = catalog.inventory([("original", str(tmp_path), 0)], "vm01")
    assert result["candidate_counts"] == {"complete": 2, "failed": 1, "unavailable": 1}
    assert sorted(source["kind"] for source in result["sources"]) == ["bundle", "raw", "raw", "unavailable"]
    assert result["bundle_hash_verification"] == "required_at_fetch"
    # Corruption must fail closed before any source is imported.
    collection_path.write_bytes(raw + b" ")
    with pytest.raises(ValueError, match="checksum mismatch"):
        catalog.inventory([("original", str(tmp_path), 0)])


def test_merge_keeps_success_and_latest_failed_attempt_without_inflating_scope():
    def source(status, accession, sha, precedence):
        return {"kind": "raw", "path": "/example", "precedence": precedence,
                "entries": [{"accession": accession, "status": status, "sha256": sha * 64}]}

    complete = source("complete", "GCA_000000001.1", "a", 0)
    sources = [complete, source("failed", "GCA_000000001.1", "b", 40),
               source("failed", "GCA_000000002.1", "c", 0),
               source("failed", "GCA_000000002.1", "d", 40),
               source("complete", "GCA_000000001.1", "a", 5)]
    result = catalog.merge_inventories([{"profile_sha256": catalog.PROFILE, "sources": sources}])
    assert result["total"] == 2
    assert result["counts"] == {"complete": 1, "failed": 1}
    assert {entry["sha256"] for item in result["sources"] for entry in item["entries"]} == {"a" * 64, "d" * 64}
    sources.append(source("complete", "GCA_000000001.1", "e", 50))
    with pytest.raises(ValueError, match="Conflicting successful"):
        catalog.merge_inventories([{"profile_sha256": catalog.PROFILE, "sources": sources}])


def test_rejects_unexpected_profile_and_escaping_source_path(tmp_path):
    write_json(tmp_path / "assemblies/GCA_000000001.1.json.gz", {
        "assembly": {"accession": "GCA_000000001.1"}, "status": "complete",
        "provenance": {"analysis_profile_sha256": "b" * 64},
    }, True)
    with pytest.raises(ValueError, match="Unexpected scientific profile"):
        catalog.inventory([("original", str(tmp_path), 0)])
    with pytest.raises(ValueError, match="Unsafe source member path"):
        catalog.safe_child(tmp_path, "../other-project/private.json")
