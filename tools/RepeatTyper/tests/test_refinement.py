"""Independent controls for exact features, bounded parsing and projection."""
import hashlib
from collections import Counter
from itertools import product

import numpy as np
import pytest

from repeattyper import Record, classify, inference, model, output
from repeattyper import input as fasta
from repeattyper.errors import InputError
from repeattyper.features import FEATURE_NAMES, matrix

DIGEST = "0" * 64


def test_all_4mer_features_match_independent_counting():
    sequences = ["".join(bases) for bases in product("ACGT", repeat=4)]
    sequences += ["A", "AC", "ACG", "ACGT" * 50, "A" * 200, "C" * 200]
    complement = dict(zip("ACGT", "TGCA", strict=True))
    expected = np.zeros((len(sequences) + 1, len(FEATURE_NAMES)), dtype=np.float32)
    for row, sequence in enumerate(sequences, start=1):
        windows = [sequence[offset:offset + 4] for offset in range(len(sequence) - 3)]
        canonical = [min(window, "".join(complement[base] for base in window[::-1]))
                     for window in windows]
        counts = Counter(canonical)
        counts["GC"] = (sequence.count("G") + sequence.count("C")) / len(sequence)
        counts["Length"] = len(sequence)
        expected[row] = [counts[name] for name in FEATURE_NAMES]
    np.testing.assert_array_equal(matrix(sequences), expected)


@pytest.mark.parametrize("threshold,status", [(0.5, "predicted"),
    (np.nextafter(0.5, 1.0), "uncertain")])
def test_ties_and_threshold_boundary_preserve_first_label(monkeypatch, threshold, status):
    scores = np.zeros((2, 37), dtype=np.float32)
    scores[:, :2] = 0.5
    booster = type("Booster", (), {"predict": lambda *args, **kwargs: scores})()
    monkeypatch.setattr(inference, "load_model", lambda: (booster, {"labels": ["I-A", "I-B"] + ["I-C"] * 35}))
    row = classify([Record("a", "ACGT")], input_sha256=DIGEST, threshold=threshold)["results"][0]
    assert row["winner"] == "I-A"
    assert row["native_score"] == 0.5
    assert row["status"] == status
    assert row["prediction"] == ("I-A" if status == "predicted" else None)


def test_duplicate_occurrences_hash_each_unique_sequence_once(monkeypatch):
    booster, provenance = model.load_model()
    monkeypatch.setattr(inference, "load_model", lambda: (booster, provenance))
    real_sha256 = hashlib.sha256
    hashed = []

    def sha256(payload):
        hashed.append(payload)
        return real_sha256(payload)

    monkeypatch.setattr(inference.hashlib, "sha256", sha256)
    report = classify([Record("a", "ACGT"), Record("b", "ACGT"),
                       Record("c", "ACGTN"), Record("d", "ACGTN")], input_sha256=DIGEST)
    assert hashed == [b"ACGT", b"ACGTN"]
    assert [row["id"] for row in report["results"]] == ["a", "b", "c", "d"]
    assert report["counts"]["distinct_classified_sequences"] == 1
    assert report["counts"]["unsupported"] == 2
    for row in report["results"]:
        assert row["sequence_sha256"] == real_sha256(row["sequence"].encode("ascii")).hexdigest()


@pytest.mark.parametrize("text,limit,created_limit", [
    (">a\n" + "A" * 201 + "\n", "1–200", 0),
    ("".join(f">r{i}\nA\n" for i in range(10001)), "10000", 10000),
])
def test_fasta_stops_before_materializing_over_limit_records(tmp_path, monkeypatch, text, limit, created_limit):
    path = tmp_path / "repeats.fa"
    path.write_text(text)
    real_record = fasta.Record
    created = []

    class TrackedRecord(real_record):
        def __init__(self, identifier, sequence):
            created.append(identifier)
            super().__init__(identifier, sequence)

    monkeypatch.setattr(fasta, "Record", TrackedRecord)
    with pytest.raises(InputError, match=limit):
        fasta.read_fasta(path)
    assert len(created) == created_limit


def test_fragmented_fasta_limit_and_empty_records(tmp_path):
    path = tmp_path / "repeats.fa"
    path.write_text(">a description\n" + "acgt\n" * 50 + ">b\nN\n")
    records, digest = fasta.read_fasta(path)
    assert records == [Record("a", "ACGT" * 50), Record("b", "N")]
    assert digest == hashlib.sha256(path.read_bytes()).hexdigest()
    path.write_text(">a\n" + "a\n" * 201)
    with pytest.raises(InputError, match="1–200"):
        fasta.read_fasta(path)
    path.write_text(">a\n>b\nACGT\n")
    with pytest.raises(InputError, match="1–200"):
        fasta.read_fasta(path)


def test_fasta_record_and_byte_limits_accept_boundaries(tmp_path):
    path = tmp_path / "repeats.fa"
    path.write_text("".join(f">r{i}\nA\n" for i in range(fasta.MAX_RECORDS)))
    records, _ = fasta.read_fasta(path)
    assert len(records) == fasta.MAX_RECORDS
    assert records[-1] == Record("r9999", "A")
    payload = b">a\nACGT\n"
    path.write_bytes(payload + b"\n" * (fasta.MAX_BYTES - len(payload)))
    records, digest = fasta.read_fasta(path)
    assert records == [Record("a", "ACGT")]
    assert digest == hashlib.sha256(path.read_bytes()).hexdigest()


def test_publication_write_failure_cleans_temporary_output(tmp_path, monkeypatch):
    report = classify([], input_sha256=DIGEST)
    destination = tmp_path / "results"

    def failed_write(*args, **kwargs):
        raise OSError("simulated TSV write failure")

    monkeypatch.setattr(type(destination), "write_text", failed_write)
    with pytest.raises(OSError, match="simulated TSV"):
        output.publish(report, destination)
    assert not destination.exists()
    assert not list(tmp_path.glob(".repeattyper-*"))
