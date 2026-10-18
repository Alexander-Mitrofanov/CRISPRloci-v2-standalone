"""Large references remain immutable without loading or publishing every genome."""

import hashlib
import json
import os
import sqlite3

import pytest

from crispr_integration import association_reference as adapter
from crispr_integration.association import validate_report as adapter_report
from crispr_web.execution import validation_association_reference as backend
from crispr_web.execution.validation_association import (
    validate_report as backend_report,
)
from crispr_web.summaries.association import public_report


def sha(data):
    return hashlib.sha256(data).hexdigest()


def identity(path):
    data = path.read_bytes()
    return {"size_bytes": len(data), "sha256": sha(data)}


def save_manifest(root, doc):
    data = json.dumps(doc, sort_keys=True).encode()
    (root / "reference.json").write_bytes(data)
    return sha(data)


def catalogue(root, *, count=130, length=400_001):
    (root / "genomes").mkdir(parents=True)
    (root / "db").mkdir()
    # The catalogue is intentionally larger than the former 100-genome/50-Mb
    # panel. Only selected genomes contain full source bytes in this validator
    # fixture: reference checks must stat but not read FASTAs without a hit.
    with sqlite3.connect(root / "catalogue.sqlite3") as db:
        db.execute(
            "CREATE TABLE genomes(id TEXT PRIMARY KEY, length INTEGER, sha256 TEXT, file TEXT)"
        )
        for i in range(count):
            sequence = "N" + "A" * (length - 1)
            db.execute(
                "INSERT INTO genomes VALUES (?,?,?,?)",
                (f"NC_{i:06}.1", length, sha(sequence.encode()), f"genomes/{i}.fasta"),
            )
            source = sequence if i < 10 else "placeholder"
            (root / "genomes" / f"{i}.fasta").write_text(f">NC_{i:06}.1\n{source}\n")
    members = {}
    for prefix in ("target", "control"):
        for suffix in ("", ".index", ".dbtype"):
            path = root / "db" / (prefix + suffix)
            path.write_bytes(b"native fixture bytes\n")
            members[f"db/{path.name}"] = identity(path)
    doc = {
        "schema_version": "spacepharer.reference/2.0.0",
        "name": "ncbi-complete-phage-snapshot",
        "engine_version": "2.1.0",
        "control_method": "reversed_translated_orfs",
        "genome_count": count,
        "total_bases": count * length,
        "scope": "broad_phage_reference",
        "source": {
            "collection": "NCBI RefSeq",
            "release": "test",
            "retrieved_at": "2026-10-10",
            "selection_policy": "Fixture for bounded complete-reference validation.",
        },
        "catalogue": {
            "file": "catalogue.sqlite3",
            **identity(root / "catalogue.sqlite3"),
        },
        "files": members,
    }
    return doc, save_manifest(root, doc)


@pytest.fixture(params=[adapter, backend], ids=["Integration", "backend"])
def layer(request):
    return request.param


def test_broad_reference_lazily_validates_only_used_targets(
    tmp_path, layer, monkeypatch
):
    doc, digest = catalogue(tmp_path)
    calls = []
    original = layer._read_genome
    monkeypatch.setattr(
        layer, "_read_genome", lambda path: (calls.append(path), original(path))[1]
    )
    public, targets = layer.reference_data(tmp_path, digest)
    try:
        assert public["genome_count"] == 130 and public["total_bases"] > 50_000_000
        assert (
            "genomes" not in public
            and "files" not in public
            and "catalogue" not in public
        )
        assert len(json.dumps(public)) < 1000 and calls == []
        assert "NC_000129.1" in targets  # This metadata lookup never loads a genome.
        assert "missing" not in targets
        assert targets["NC_000000.1"].startswith("N")
        assert targets["NC_000000.1"].startswith("N") and len(calls) == 1
        for i in range(1, 10):
            targets[f"NC_{i:06}.1"]
        assert len(targets.cache) == 8
        targets["NC_000000.1"]
        assert len(calls) == 11
        with pytest.raises(KeyError):
            targets["missing"]
    finally:
        targets.close()


@pytest.mark.parametrize("member", ["db/target", "catalogue.sqlite3"])
def test_reference_hashes_reject_same_size_mtime_restored_changes(
    tmp_path, layer, member
):
    _, digest = catalogue(tmp_path, count=2, length=30)
    path = tmp_path / member
    before = path.stat()
    data = bytearray(path.read_bytes())
    data[-1] ^= 1
    path.write_bytes(data)
    os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
    with pytest.raises(ValueError, match="checksum"):
        layer.reference_data(tmp_path, digest)


@pytest.mark.parametrize(
    "change",
    [
        "count",
        "bases",
        "boolean",
        "path",
        "scope",
        "metadata",
        "inventory",
        "catalogue-size",
        "index-size",
        "missing-genome",
        "duplicate-source",
    ],
)
def test_reference_manifest_and_catalogue_fail_closed(tmp_path, layer, change):
    doc, _ = catalogue(tmp_path, count=2, length=30)
    if change == "count":
        doc["genome_count"] += 1
    elif change == "bases":
        doc["total_bases"] += 1
    elif change == "boolean":
        doc["genome_count"] = True
    elif change == "path":
        with sqlite3.connect(tmp_path / "catalogue.sqlite3") as db:
            db.execute("UPDATE genomes SET file='../outside.fasta'")
        doc["catalogue"].update(identity(tmp_path / "catalogue.sqlite3"))
    elif change == "scope":
        doc["scope"] = "all_phages_ever"
    elif change == "metadata":
        doc["source"]["collection"] = "NCBI\nprivate implementation detail"
    elif change == "inventory":
        (tmp_path / "db/unbound").write_bytes(b"extra")
    elif change == "catalogue-size":
        doc["catalogue"]["size_bytes"] = 129 * 1024 * 1024
    elif change == "index-size":
        doc["files"]["db/target"]["size_bytes"] = 17 * 1024**3
    elif change == "missing-genome":
        (tmp_path / "genomes/1.fasta").unlink()
    elif change == "duplicate-source":
        with sqlite3.connect(tmp_path / "catalogue.sqlite3") as db:
            db.execute("UPDATE genomes SET file='genomes/0.fasta'")
        doc["catalogue"].update(identity(tmp_path / "catalogue.sqlite3"))
    with pytest.raises(ValueError):
        layer.reference_data(tmp_path, save_manifest(tmp_path, doc))


@pytest.mark.parametrize("change", ["sequence", "identifier", "symlink"])
def test_lazy_target_source_binding_and_symlink_rejection(tmp_path, layer, change):
    _, digest = catalogue(tmp_path, count=2, length=30)
    _, targets = layer.reference_data(tmp_path, digest)
    try:
        path = tmp_path / "genomes/0.fasta"
        if change == "sequence":
            path.write_text(">NC_000000.1\n" + "T" * 30 + "\n")
        elif change == "identifier":
            path.write_text(">NC_000001.1\nN" + "A" * 29 + "\n")
        else:
            path.unlink()
            path.symlink_to(tmp_path / "genomes/1.fasta")
        with pytest.raises(ValueError):
            targets["NC_000000.1"]
    finally:
        targets.close()


def test_reference_verification_obeys_cancellation_checkpoint(tmp_path, layer):
    _, digest = catalogue(tmp_path, count=2, length=30)

    def stop():
        raise InterruptedError("cancelled")

    with pytest.raises(InterruptedError, match="cancelled"):
        layer.reference_data(tmp_path, digest, checkpoint=stop)


def test_public_reference_is_bounded_and_legacy_evidence_still_exports(tmp_path):
    _, digest = catalogue(tmp_path, count=2, length=30)
    reference, targets = backend.reference_data(tmp_path, digest)
    targets.close()
    report = {
        "schema_version": "spacepharer.association/2.0.0",
        "reference": reference,
        "requested_fdr": 0.05,
        "calibration_status": "not_run",
        "calibration": None,
        "records": [],
        "groups": [],
        "associations": [],
    }
    public = public_report(report)
    assert public["schema"]["version"] == "2.0.0"
    assert public["reference"] == reference
    assert "catalogue.sqlite3" not in json.dumps(public)
    report["schema_version"] = "spacepharer.association/1.0.0"
    report["reference"] = {
        key: reference[key]
        for key in ("name", "engine_version", "control_method", "manifest_sha256")
    } | {"genomes": [{"id": "legacy", "length": 30, "sha256": "a" * 64}]}
    assert public_report(report)["schema"]["version"] == "1.0.0"
    assert public_report(report)["reference"]["genomes"][0]["id"] == "legacy"


@pytest.mark.parametrize("validate", [adapter_report, backend_report])
@pytest.mark.parametrize("population", [10_000_000, 10_000_001])
def test_calibration_population_bound_tracks_reference_and_group_limits(
    validate, population
):
    records = {"spacer": "ACGT" * 8}
    groups = {
        "schema_version": "spacepharer.groups/1.0.0",
        "groups": [{"id": "group", "members": ["spacer"]}],
    }
    reference = {"genome_count": 100_000}
    doc = {
        "schema_version": "spacepharer.association/2.0.0",
        "tool_version": "2.1.0",
        "input_sha256": "a" * 64,
        "groups_sha256": "b" * 64,
        "requested_fdr": 0.05,
        "calibration_status": "estimated",
        "calibration": {
            "schema_version": "spacepharer.calibration/1.0.0",
            "status": "estimated",
            "target_score_count": population,
            "control_score_count": population,
        },
        "warnings": [],
        "reference": reference,
        "records": [{"id": "spacer", "sequence": records["spacer"]}],
        "groups": [
            {
                "id": "group",
                "members": ["spacer"],
                "supported_occurrence_count": 1,
                "unique_sequence_count": 1,
                "unsupported_occurrence_ids": [],
            }
        ],
        "associations": [],
        "interpretation": "Spacer similarity supports a group–target association hypothesis, not active infection, cleavage or a PAM prediction.",
    }
    args = (doc, records, groups, reference, {}, "a" * 64, "b" * 64, 0.05)
    if population == 10_000_000:
        assert validate(*args)["groups"] == 1
    else:
        with pytest.raises(ValueError, match="calibration"):
            validate(*args)
