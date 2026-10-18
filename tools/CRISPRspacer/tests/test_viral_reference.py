"""Tiny real-BLAST references plus independent integrity/budget regressions."""
import gzip
import hashlib
import json
import shutil
import sqlite3

import pytest

from crisprspacer.io import InputError, LimitError, write_json
from crisprspacer import viral_reference as vr


REAL_BLAST = shutil.which("makeblastdb") and shutil.which("blastdbcmd")


def source(tmp_path, body=">NC_000001.1 First genome\naCgTNNACGTTGCAACGT\n>NC_000002.2 Second genome\nGATTACAACCTG\n", name="viral.fna"):
    path = tmp_path / name
    path.write_text(body)
    return path


def build(paths, destination, **kw):
    return vr.build_viral_reference(paths, destination, name="synthetic-viral", source_release="test-1",
                                    source_url="https://example.org/viral/", **kw)


@pytest.fixture
def reference(tmp_path):
    if not REAL_BLAST:
        pytest.skip("Real local BLAST+ executables are not installed")
    directory = tmp_path / "reference"
    ref = build([source(tmp_path)], directory)
    return directory, ref


def rehash(directory):
    ref = json.loads((directory / "reference.json").read_text())
    ref["files"] = {p.name: vr._file_identity(p) for p in directory.iterdir() if p.name != "reference.json"}
    write_json(directory / "reference.json", ref)
    return ref


def test_real_reference_preserves_accession_order_content_and_immutable_reads(reference):
    directory, ref = reference
    loaded, identity = vr.load_viral_reference(directory)
    assert loaded == ref
    assert identity == hashlib.sha256((directory / "reference.json").read_bytes()).hexdigest()
    assert ref["record_count"] == 2 and ref["total_bases"] == 30
    assert ref["build"]["database_format"] == 5
    connection = vr.open_subjects(directory)
    try:
        rows = [dict(r) for r in connection.execute("SELECT * FROM subjects ORDER BY id")]
        assert [(r["id"], r["accession"], r["title"]) for r in rows] == [
            ("v000000001", "NC_000001.1", "First genome"), ("v000000002", "NC_000002.2", "Second genome")]
        assert rows[0]["sequence"] == "ACGTNNACGTTGCAACGT"
        with pytest.raises(sqlite3.OperationalError):
            connection.execute("DELETE FROM subjects")
    finally:
        connection.close()
    assert vr.validate_viral_reference(directory, identity)["status"] == "valid"


def test_gzip_multisource_and_receipt_bound_to_exact_source_bytes(tmp_path):
    if not REAL_BLAST:
        pytest.skip("BLAST+ unavailable")
    plain = source(tmp_path)
    zipped = tmp_path / "more.fna.gz"
    zipped.write_bytes(gzip.compress(b">NC_000003.3 Third genome\nAACCGGTTAACCGGTT\n"))
    rows = []
    for path in (plain, zipped):
        raw = path.read_bytes()
        rows.append({"name": path.name, "url": "https://example.org/viral/" + path.name,
                     "bytes": len(raw), "md5": hashlib.md5(raw).hexdigest(), "sha256": hashlib.sha256(raw).hexdigest()})
    receipt = tmp_path / "download.json"
    write_json(receipt, {"schema": {"name": "crisprspacer.ncbi-viral-download", "version": "1.0.0"},
                        "source_release": "test-1", "source_url": "https://example.org/viral/", "files": rows})
    ref = build([plain, zipped], tmp_path / "ref", source_manifest=receipt)
    assert ref["record_count"] == 3
    assert ref["source_files"] == rows
    assert (tmp_path / "ref/source-manifest.json").read_bytes() == receipt.read_bytes()
    assert vr.load_viral_reference(tmp_path / "ref")[0] == ref


@pytest.mark.parametrize("body,pattern", [
    (">dup.1\nACGT\n>dup.1\nTGCA\n", "Duplicate"),
    (">empty.1\n>next.1\nACGT\n", "Empty"),
    ("ACGT\n", "precedes"),
    (">rna.1\nACGU\n", "IUPAC"),
    (">" + "a" * 2000 + "\nACGT\n", "header"),
])
def test_import_errors_are_atomic(tmp_path, monkeypatch, body, pattern):
    monkeypatch.setattr(vr, "executable_identity", lambda p: {"path": "/unused/" + p, "version": "test", "sha256": "0" * 64})
    with pytest.raises(InputError, match=pattern):
        build([source(tmp_path, body)], tmp_path / "ref")
    assert not (tmp_path / "ref").exists()
    assert not (tmp_path / "ref.partial").exists()
    assert not list(tmp_path.glob(".crisprspacer-*"))


@pytest.mark.parametrize("budget,value,pattern", [("MAX_RECORDS", 1, "record budget"),
                                                  ("MAX_SEQUENCE_BASES", 10, "length budget"),
                                                  ("MAX_CORPUS_BYTES", 20, "compressed byte")])
def test_independent_operator_import_limits(tmp_path, monkeypatch, budget, value, pattern):
    monkeypatch.setattr(vr, "executable_identity", lambda p: {"path": "/unused/" + p, "version": "test", "sha256": "0" * 64})
    monkeypatch.setattr(vr, budget, value)
    with pytest.raises(LimitError, match=pattern):
        build([source(tmp_path)], tmp_path / "ref")
    assert not (tmp_path / "ref").exists()


def test_decompressed_limit_independent_of_small_gzip_file(tmp_path, monkeypatch):
    monkeypatch.setattr(vr, "executable_identity", lambda p: {"path": "/unused/" + p, "version": "test", "sha256": "0" * 64})
    monkeypatch.setattr(vr, "MAX_CORPUS_BYTES", 150)
    path = tmp_path / "s.gz"
    path.write_bytes(gzip.compress(b">accession.1\n" + b"A" * 200 + b"\n"))
    with pytest.raises(LimitError, match="decompressed"):
        build([path], tmp_path / "ref")


def test_disk_reserve_aborts_before_publication(tmp_path, monkeypatch):
    monkeypatch.setattr(vr, "executable_identity", lambda p: {"path": "/unused/" + p, "version": "test", "sha256": "0" * 64})
    monkeypatch.setattr(vr.shutil, "disk_usage", lambda p: type("Capacity", (), {"free": 1})())
    with pytest.raises(LimitError, match="free disk"):
        build([source(tmp_path)], tmp_path / "ref")
    assert not (tmp_path / "ref").exists()


def test_pinned_digest_extra_file_and_symlink_rejected(reference, tmp_path):
    directory, _ = reference
    with pytest.raises(InputError, match="pinned"):
        vr.load_viral_reference(directory, "0" * 64)
    extra = directory / "unexpected"
    extra.write_text("bad")
    with pytest.raises(InputError, match="Unexpected"):
        vr.load_viral_reference(directory)
    extra.unlink()
    link = tmp_path / "link"
    link.symlink_to(directory, target_is_directory=True)
    with pytest.raises(InputError, match="Symlinks"):
        vr.load_viral_reference(link)


def test_tampered_native_volume_rejected(reference):
    directory, _ = reference
    with (directory / "db.nsq").open("r+b") as f:
        first = f.read(1)
        f.seek(0)
        f.write(bytes([first[0] ^ 1]))
    with pytest.raises(InputError, match="integrity"):
        vr.load_viral_reference(directory)


def test_regenerated_manifest_cannot_hide_catalog_source_disagreement(reference):
    directory, _ = reference
    connection = sqlite3.connect(directory / "subjects.sqlite3")
    connection.execute("UPDATE subjects SET sequence='T'||substr(sequence,2) WHERE id='v000000001'")
    connection.commit()
    connection.close()
    rehash(directory)
    with pytest.raises(InputError, match="sequence/hash"):
        vr.validate_viral_reference(directory)


def test_regenerated_source_hash_cannot_hide_wrong_native_subjects(reference):
    directory, _ = reference
    connection = sqlite3.connect(directory / "subjects.sqlite3")
    seq = "TTTTNNACGTTGCAACGT"
    connection.execute("UPDATE subjects SET sequence=?,sha256=? WHERE id='v000000001'", (seq, hashlib.sha256(seq.encode()).hexdigest()))
    connection.commit()
    connection.row_factory = sqlite3.Row
    content = hashlib.sha256()
    with (directory / "subjects.fasta").open("wb") as f:
        for row in connection.execute("SELECT * FROM subjects ORDER BY id"):
            content.update(vr.canonical(vr._content_row(row)))
            f.writelines(vr._normalized_chunks(row["id"], row["sequence"]))
    connection.close()
    ref = rehash(directory)
    ref["content_sha256"] = content.hexdigest()
    write_json(directory / "reference.json", ref)
    with pytest.raises(InputError, match="roundtrip differs"):
        vr.validate_viral_reference(directory)


def test_native_alias_cannot_escape_even_after_rehash(reference):
    directory, _ = reference
    (directory / "db.nal").write_text("TITLE hostile\nDBLIST /tmp/other\n")
    rehash(directory)
    with pytest.raises(InputError, match="escapes"):
        vr.load_viral_reference(directory)


def test_receipt_mismatch_and_changed_binary_are_rejected(reference, tmp_path):
    directory, _ = reference
    receipt = tmp_path / "download.json"
    write_json(receipt, {"schema": {"name": "crisprspacer.ncbi-viral-download", "version": "1.0.0"},
                        "source_release": "test-1", "source_url": "https://example.org/viral/", "files": []})
    with pytest.raises(InputError, match="receipt"):
        build([source(tmp_path, name="other.fna")], tmp_path / "other-ref", source_manifest=receipt)
    ref = json.loads((directory / "reference.json").read_text())
    ref["executables"]["blastdbcmd"]["sha256"] = "0" * 64
    write_json(directory / "reference.json", ref)
    with pytest.raises(InputError, match="receipt"):
        vr.validate_viral_reference(directory)


def test_normalized_fasta_wraps_long_sequences(tmp_path):
    if not REAL_BLAST:
        pytest.skip("BLAST+ unavailable")
    sequence = "ACGT" * 26
    ref = build([source(tmp_path, ">NC_1.1 title\n" + sequence + "\n")], tmp_path / "ref")
    lines = (tmp_path / "ref/subjects.fasta").read_text().splitlines()
    assert lines == [">v000000001", sequence[:80], sequence[80:]]
    assert ref["total_bases"] == 104


def test_native_failure_cleans_build_and_reservation(tmp_path, monkeypatch):
    monkeypatch.setattr(vr, "executable_identity", lambda p: {"path": "/unused/" + p, "version": "test", "sha256": "0" * 64})

    def fail(*args, **kwargs):
        (kwargs["cwd"] / "db.nsq").write_text("unfinished native output")
        raise LimitError("native timeout")

    monkeypatch.setattr(vr, "run_blast", fail)
    with pytest.raises(LimitError, match="native timeout"):
        build([source(tmp_path)], tmp_path / "ref")
    assert not (tmp_path / "ref").exists()
    assert not (tmp_path / "ref.partial").exists()
    assert not list(tmp_path.glob(".crisprspacer-*"))


def test_no_overwrite_or_symlink_source(tmp_path, reference):
    directory, _ = reference
    before = (directory / "reference.json").read_bytes()
    with pytest.raises(InputError, match="already exists"):
        build([source(tmp_path, name="another.fna")], directory)
    assert (directory / "reference.json").read_bytes() == before
    link = tmp_path / "linked.fna"
    link.symlink_to(tmp_path / "another.fna")
    with pytest.raises(InputError, match="Symlinks"):
        build([link], tmp_path / "new")
