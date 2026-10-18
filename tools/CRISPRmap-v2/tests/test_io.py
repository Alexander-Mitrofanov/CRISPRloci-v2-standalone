import json

import pytest

from crisprmap import io
from crisprmap.io import InputError, LimitError, atomic_json, decode_json, digest, read_fasta


def test_fasta_snapshot_normalization_description_and_hash(tmp_path):
    path = tmp_path / "repeats.fa"
    raw = b">one annotation words\naCuT\n>two\nNNRY\n"
    path.write_bytes(raw)
    records, sha = read_fasta(path)
    assert records == [{"id": "one", "description": "annotation words", "sequence": "ACTT"},
                       {"id": "two", "description": "", "sequence": "NNRY"}]
    assert sha == digest(raw)


@pytest.mark.parametrize("payload", [b"", b"ACGT\n", b">\nACGT\n", b">a\n>b\nACGT\n",
                                       b">a\nACGT\n>a\nAAAA\n", b">a\nA-C\n", b">a\nA C\n",
                                       b">a\nA\xffC\n", b">a\n" + b"A" * 201])
def test_invalid_fasta_rejected(tmp_path, payload):
    path = tmp_path / "bad.fa"
    path.write_bytes(payload)
    with pytest.raises(InputError):
        read_fasta(path)


@pytest.mark.parametrize("payload", [b'{"a":1,"a":2}', b'{"a":NaN}', b'{"a":Infinity}',
                                      b'{"a":1e999}', b'[' * 2000, b'{"a":' + b'9' * 200 + b'}'])
def test_strict_json_rejected(payload):
    with pytest.raises(InputError):
        decode_json(payload)


def test_symlink_and_byte_budget_rejection(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">a\nACGT\n")
    link = tmp_path / "link.fa"
    link.symlink_to(source)
    with pytest.raises(InputError, match="symlink"):
        read_fasta(link)
    with pytest.raises(LimitError):
        io.read_bytes(source, limit=2)


def test_atomic_file_refuses_clobber_and_cleans_failed_staging(tmp_path):
    output = tmp_path / "reference.json"
    atomic_json(output, {"first": True})
    with pytest.raises(InputError, match="refusing replacement"):
        atomic_json(output, {"second": True})
    assert json.loads(output.read_text()) == {"first": True}
    assert not list(tmp_path.glob(".crisprmap-*"))


def test_atomic_file_cleans_failed_staging(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "MAX_BYTES", 2)
    with pytest.raises(LimitError):
        atomic_json(tmp_path / "reference.json", {"cannot": "fit"})
    assert list(tmp_path.iterdir()) == []



def test_multiline_fasta_length_accounting_resets_at_record_boundaries(tmp_path):
    path = tmp_path / "wrapped.fa"
    raw = ">one\n" + "a\n" * 200 + ">two\n" + "u\n" * 200
    path.write_text(raw)
    records, sha = read_fasta(path)
    assert [row["sequence"] for row in records] == ["A" * 200, "T" * 200]
    assert sha == digest(raw)
    path.write_text(raw + "c\n")
    with pytest.raises(LimitError, match="length budget"):
        read_fasta(path)
