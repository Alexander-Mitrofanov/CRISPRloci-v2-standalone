import json
import subprocess
import sys

import pytest

from crisprmap.cli import main
from crisprmap.reference import build_reference
from crisprmap.io import atomic_json


@pytest.fixture
def inputs(tmp_path):
    query = tmp_path / "query.fa"
    query.write_text(">long_identifier_preserved_in_full\nACGTCGATGCTAGCTAGTCGATCG\n")
    reference = tmp_path / "reference.json"
    atomic_json(reference, build_reference(
        [{"id": "reference1", "sequence": "ACGTCGATGCTAGCTAGTCGATCG"}],
        {"reference1": {"family": "F1", "motif": None}},
        {"namespace": "synthetic-test", "source_url": "https://example.org/test", "source_sha256": "a" * 64},
    ))
    return query, reference, tmp_path / "result"


def test_cli_offline_complete_artifacts_preserve_missing_and_ids(inputs):
    query, reference, output = inputs
    assert main(["classify", "--input", str(query), "--reference", str(reference), "--output", str(output)]) == 0
    report = json.loads((output / "report.json").read_text())
    row = report["results"][0]
    assert row["id"] == "long_identifier_preserved_in_full"
    assert row["annotations"]["family"]["value"] == "F1"
    assert row["annotations"]["motif"]["value"] is None
    assert len(report["input"]["sha256"]) == 64
    assert len(report["reference"]["file_sha256"]) == 64
    assert (output / "report.tsv").is_file()
    html = (output / "report.html").read_text()
    assert "No structural motif is assigned".lower() in html.lower()
    assert "<script" not in html
    assert "Content-Security-Policy" in html


def test_cli_refuses_overwrite(inputs):
    query, reference, output = inputs
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("unchanged")
    assert main(["classify", "--input", str(query), "--reference", str(reference), "--output", str(output)]) == 2
    assert list(output.iterdir()) == [sentinel]
    assert sentinel.read_text() == "unchanged"


@pytest.mark.parametrize("option,value", [("--max-distance", "-1"), ("--top-k", "0"),
                                            ("--deadline-seconds", "nan"), ("--max-comparisons", "0")])
def test_failed_analysis_publishes_no_bundle(inputs, option, value):
    query, reference, output = inputs
    assert main(["classify", "--input", str(query), "--reference", str(reference),
                 "--output", str(output), option, value]) == 2
    assert not output.exists()


def test_build_reference_records_metadata_hash_and_rejects_unknown_ids(tmp_path):
    fasta, metadata, reference = (tmp_path / name for name in ("ref.fa", "labels.tsv", "ref.json"))
    fasta.write_text(">known\nACGT\n")
    metadata.write_text("id\tfamily\tmotif\nknown\tF1\t\n")
    argv = ["build-reference", "--fasta", str(fasta), "--metadata", str(metadata),
            "--namespace", "test", "--source-url", "https://example.org/test", "--output", str(reference)]
    assert main(argv) == 0
    assert len(json.loads(reference.read_text())["provenance"]["metadata_sha256"]) == 64
    assert main(argv) == 2
    reference.unlink()
    metadata.write_text("id\tfamily\tmotif\nunknown\tF1\t\n")
    assert main(argv) == 2
    assert not reference.exists()


@pytest.mark.parametrize("payload", [
    "", "id\tfamily\n", "id\tfamily\tmotif\tfamily\nknown\tF1\tM1\tF2\n",
    "id\tfamily\tmotif\nknown\tF1\n", "id\tfamily\tmotif\nknown\tF1\tM1\textra\n",
    "id\tfamily\tmotif\nknown\t\"F1\tM1\n",
    "id\tfamily\tmotif\tunknown\n",
])
def test_malformed_metadata_is_rejected_before_reference_publication(tmp_path, payload):
    fasta, metadata, reference = (tmp_path / name for name in ("ref.fa", "labels.tsv", "ref.json"))
    fasta.write_text(">known\nACGT\n")
    metadata.write_text(payload)
    assert main(["build-reference", "--fasta", str(fasta), "--metadata", str(metadata),
                 "--namespace", "test", "--source-url", "https://example.org/test",
                 "--output", str(reference)]) == 2
    assert not reference.exists()


def test_metadata_quoted_labels_explicit_missing_values_and_reordered_header(tmp_path):
    fasta, metadata, reference = (tmp_path / name for name in ("ref.fa", "labels.tsv", "ref.json"))
    fasta.write_text(">known\nACGT\n>unannotated\nACGC\n")
    metadata.write_text('motif\tid\taccession\tfamily\n\n\tknown\tNC_000001\t"F""1"\n')
    assert main(["build-reference", "--fasta", str(fasta), "--metadata", str(metadata),
                 "--namespace", "test", "--source-url", "https://example.org/test",
                 "--output", str(reference)]) == 0
    entries = json.loads(reference.read_text())["entries"]
    assert entries[0]["family"] == 'F"1'
    assert entries[0]["motif"] is None
    assert entries[0]["accession"] == "NC_000001"
    assert all(entries[1][field] is None for field in ("family", "motif", "accession"))


def test_installed_cli_module_from_unrelated_directory(tmp_path):
    result = subprocess.run([sys.executable, "-m", "crisprmap", "--version"],
                            cwd=tmp_path, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert "CRISPRmap-v2 2.1.0a2" in result.stdout


def test_escaping_and_spreadsheet_formula_neutralization(inputs):
    query, reference, output = inputs
    query.write_text(">=danger<script>\nACGTCGATGCTAGCTAGTCGATCG\n")
    assert main(["classify", "--input", str(query), "--reference", str(reference), "--output", str(output)]) == 0
    text = (output / "report.html").read_text()
    assert "=danger&lt;script&gt;" in text
    assert "<script>" not in text
    assert "'=danger<script>" in (output / "report.tsv").read_text()


def test_output_does_not_follow_symlink_parent(inputs, tmp_path):
    query, reference, _ = inputs
    destination = tmp_path / "real"
    destination.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(destination, target_is_directory=True)
    assert main(["classify", "--input", str(query), "--reference", str(reference),
                 "--output", str(alias / "result")]) == 2
    assert not list(destination.iterdir())
