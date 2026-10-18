import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from crisprrepeat.analysis import analyze
from crisprrepeat import inputs
from crisprrepeat.inputs import load_arrays, sha
from crisprrepeat.references import References


def small_report():
    sequence = "GCGCGAAAACGCGC"
    return sequence, {
        "schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
        "source": {"id": "source", "length": len(sequence), "sha256": sha(sequence),
                   "sha256_scope": "uppercase_sequence"},
        "arrays": [{"id": "array", "category": "Bona-fide", "repeat_count": 1,
                    "source_interval": {"start": 1, "end": len(sequence)},
                    "source_sequence": sequence, "consensus_repeat": sequence,
                    "repeats": [{"index": 1, "sequence": sequence,
                                 "source_interval": {"start": 1, "end": len(sequence)}}]}],
    }


@pytest.mark.parametrize("field", ["schema", "source", "array", "repeat", "orientation", "index"])
def test_malformed_detector_objects_fail_as_validation_errors(tmp_path, field):
    sequence, report = small_report()
    if field in {"schema", "source"}:
        report[field] = []
    elif field == "array":
        report["arrays"][0] = None
    elif field == "repeat":
        report["arrays"][0]["repeats"][0] = 42
    elif field == "index":
        report["arrays"][0]["repeats"][0]["index"] = True
    else:
        report["arrays"][0]["orientation"] = []
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        load_arrays(path, {"source": sequence}, ["Bona-fide"])


def test_ambiguous_consensus_preserves_observed_repeat_analysis(tmp_path):
    sequence, report = small_report()
    report["arrays"][0]["consensus_repeat"] = "N" * len(sequence)
    path, fasta = tmp_path / "report.json", tmp_path / "sequence.fa"
    path.write_text(json.dumps(report))
    fasta.write_text(">source\n" + sequence + "\n")
    result = analyze(fasta, reports=path, orientation="forward")
    assert result["consensus_folds"][0]["isolated"] == {
        "status": "unsupported", "reason": "ambiguous_bases"}
    assert result["instances"][0]["isolated"]["status"] == "completed"
    assert result["contexts"][0]["status"] == "completed"


def reference_document():
    return {
        "schema": {"name": "CRISPRrepeat-reference", "version": "1.0.0"},
        "release": "test", "license": "CC0", "source_url": "https://example.org/fixture",
        "build_date": "2026-09-27", "entries": [{"id": "α", "sequence": "ACGT", "sha256": sha("ACGT"),
        "occurrences": [{"accession": "測試", "source_url": "https://example.org/fixture",
                         "evidence": "δοκιμή"}]}],
    }


@pytest.mark.parametrize("field", ["entry", "occurrence", "id"])
def test_reference_nested_types_are_validated(tmp_path, field):
    reference = reference_document()
    if field == "entry":
        reference["entries"][0] = 42
    elif field == "occurrence":
        reference["entries"][0]["occurrences"][0] = []
    else:
        reference["entries"][0]["id"] = []
    path = tmp_path / "reference.json"
    path.write_text(json.dumps(reference))
    with pytest.raises(ValueError):
        References(path)


def test_unicode_reference_hash_binds_the_single_parsed_byte_snapshot(tmp_path, monkeypatch):
    payload = json.dumps(reference_document(), ensure_ascii=False, indent=2).replace("\n", "\r\n").encode()
    path = tmp_path / "reference.json"
    path.write_bytes(payload)
    # No second unbounded read may provide a hash for different file contents.
    monkeypatch.setattr(Path, "read_bytes", lambda path: pytest.fail("Reference was reread"))
    result = References(path).search("ACGU")
    assert result["reference_sha256"] == hashlib.sha256(payload).hexdigest()
    assert result["hits"][0]["occurrences"][0]["accession"] == "測試"


def test_directory_inventory_stops_at_the_bound(monkeypatch):
    monkeypatch.setattr(inputs, "MAX_RECORDS", 2)

    def inventory():
        for i in range(3):
            yield Path(str(i))
        pytest.fail("Unbounded inventory was consumed")

    with pytest.raises(ValueError, match="Too many"):
        inputs.bounded_paths(inventory(), "files")


def test_read_limit_counts_utf8_bytes(tmp_path):
    path = tmp_path / "unicode.txt"
    path.write_bytes("éé".encode())
    assert inputs.bounded_text(path, 4) == "éé"
    with pytest.raises(ValueError, match="oversized"):
        inputs.bounded_text(path, 3)


def test_repeat_length_limit_stops_consuming_an_oversized_record(tmp_path, monkeypatch):
    class Input:
        def splitlines(self):
            yield ">repeat"
            yield "A" * 201
            pytest.fail("Oversized repeat input was still consumed")

    monkeypatch.setattr(inputs, "bounded_text", lambda path: Input())
    with pytest.raises(ValueError, match="exceeds 200"):
        inputs.read_fasta(tmp_path / "repeat.fa", repeats=True)


def test_cli_rna_defaults_to_transcribed_forward(tmp_path):
    fasta, output = tmp_path / "rna.fa", tmp_path / "results"
    fasta.write_text(">rna\nGUUCAACACCCUCUUUUCCCCGUCAGGGGACUGAAAC\n")
    subprocess.run([sys.executable, "-m", "crisprrepeat.cli", "--input", str(fasta),
                    "--output", str(output), "--molecule", "RNA"], check=True, capture_output=True)
    report = json.loads((output / "repeat-analysis.json").read_text())
    assert report["orientation_policy"] == "forward"
    assert report["counts"]["instances"] == 1
