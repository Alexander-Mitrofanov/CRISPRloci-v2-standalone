from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import catalog_raw
from catalog_snapshot import STRIDE


def source_item(tmp_path, payload, kind="raw", entries=None):
    accession = payload.get("assembly", {}).get("accession", "unavailable")
    path = tmp_path / f"{accession}.json.gz"
    data = gzip.compress(json.dumps(payload).encode(), mtime=0)
    path.write_bytes(data)
    digest = hashlib.sha256(data).hexdigest()
    if entries is None:
        entries = [{"accession": accession, "status": payload["status"], "sha256": digest}]
    else:
        entries = [{**entry, "sha256": digest} for entry in entries]
    return {"path": str(path), "sha256": digest, "kind": kind, "entries": entries}


def specimen():
    gene = {"call_id": "protein-1", "sequence_id": "contig-1", "interval": {"start": 100, "end": 150},
            "native": {"annotation": "Cas1", "profile": "Cas1-model", "score": 0.0,
                       "score_kind": "profile score", "strand": "+", "translation_table": 11,
                       "protein_length": 16, "is_cas": True}}
    system = {"call_id": "cassette-1", "interval": {"start": 100, "end": 150},
              "native": {"classification": {"type": "I", "subtype": None, "confidence": {"vote_fraction": 0.8}},
                         "genes": [{"protein_id": "protein-1"}]}}
    array = {"call_id": "array-1", "interval": {"start": 1, "end": 20}, "repeat_consensus": "ACGT",
             "confidence": {"category": "Possible", "score": 0.0},
             "repeats": [{"index": 1, "interval": {"start": 1, "end": 4}, "sequence": "ACGT"}],
             "spacers": [{"index": 1, "interval": {"start": 5, "end": 8}, "sequence": "GGTT"}],
             "spacer_deletions": [{"index": 2, "source_boundary": {"position": 8}, "length": 0}]}
    return {"assembly": {"accession": "GCA_000000001.1"}, "status": "complete",
            "validation": {"status": "valid"}, "provenance": {"analysis_profile_sha256": catalog_raw.RAW_PROFILE},
            "sequences": [{"id": "contig-1", "length": 300, "tools": [
                {"name": "CRISPRidentify", "version": "2.3.2", "arrays": [array]},
                {"name": "CasAndra", "version": "campaign", "cas_systems": [system], "cas_genes": [gene]},
            ]}]}


def test_raw_retained_success_preserves_relations_unknowns_and_deletions(tmp_path):
    item = source_item(tmp_path, specimen())
    receipt, compressed = catalog_raw.transform_raw([item], 8)
    assert receipt["counts"] == {"assemblies": 1, "sequences": 1, "arrays": 1, "repeats": 1,
                                 "spacers": 2, "cas_systems": 1, "cas_genes": 1}
    db = sqlite3.connect(":memory:")
    try:
        db.deserialize(gzip.decompress(compressed))
        assert db.execute("SELECT strand,certainty_score FROM arrays").fetchone() == (None, 0.0)
        assert db.execute("SELECT start,end,length FROM arrays").fetchone() == (1, 20, 20)
        assert db.execute("SELECT kind,start,end,length,boundary_position,sequence FROM spacers WHERE kind='deletion'").fetchone() == ("deletion", None, None, 0, 8, None)
        assert db.execute("SELECT confidence,experimental,subtype FROM cas_systems").fetchone() == (None, None, None)
        assert db.execute("SELECT cas_system_pk,score FROM cas_genes").fetchone() == (8 * STRIDE + 1, 0.0)
        assert db.execute("PRAGMA foreign_key_check").fetchone() is None
    finally:
        db.close()


def test_raw_failure_and_unavailability_stay_explicit(tmp_path):
    failure = source_item(tmp_path, {"assembly": {"accession": "GCA_000000002.1"}, "status": "failed",
                                    "provenance": {"analysis_profile_sha256": catalog_raw.RAW_PROFILE}})
    unavailable = source_item(tmp_path, {"unavailable": {"GCA_000000003.1": {"reason": "no NCBI input"}}},
                              "unavailable", [{"accession": "GCA_000000003.1", "status": "unavailable"}])
    receipt, _ = catalog_raw.transform_raw([failure, unavailable], 9)
    assert [(row[1], row[2], row[3]) for row in receipt["assemblies"]] == [
        ("GCA_000000002.1", "failed", None), ("GCA_000000003.1", "unavailable", None)]
    assert receipt["counts"]["arrays"] == 0
    assert receipt["source_documents"][1]["kind"] == "unavailable"


def test_raw_refuses_changed_or_unvalidated_source(tmp_path):
    item = source_item(tmp_path, specimen())
    Path(item["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError, match="checksum mismatch"):
        catalog_raw.transform_raw([item], 1)
    payload = specimen()
    payload["provenance"]["analysis_profile_sha256"] = "incorrect"
    item = source_item(tmp_path, payload)
    with pytest.raises(ValueError, match="Unexpected raw scientific profile"):
        catalog_raw.transform_raw([item], 1)
    payload = specimen()
    payload["validation"]["status"] = "unknown"
    item = source_item(tmp_path, payload)
    with pytest.raises(ValueError, match="not passed source validation"):
        catalog_raw.transform_raw([item], 1)
