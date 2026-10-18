"""Opt-in real-tool tests, configured through CRISPRTRACRRNA_TEST_CONFIG."""

import csv
import os
from pathlib import Path

import pytest

from crisprtracrrna.config import load_settings
from crisprtracrrna.models import Interval
from crisprtracrrna.pipeline import run
from crisprtracrrna.sequence import read_fasta, reverse_complement

ROOT = Path(__file__).resolve().parents[1]
CONFIG = os.environ.get("CRISPRTRACRRNA_TEST_CONFIG")
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not CONFIG, reason="Real-tool config not set"),
]


def test_real_models_multifasta_and_reverse_complement(tmp_path):
    settings = load_settings(Path(CONFIG), {"mode": "model"})
    original = next(read_fasta(ROOT / "examples/model-II.fasta"))
    # A real CM-positive region from the complete supplied genome is found first.
    full = run(ROOT / "examples/model-II.fasta", tmp_path / "full", settings)
    assert full["summary"]["prediction_count"] > 0
    first = full["records"][0]["predictions"][0]["interval"]
    interval = Interval(max(0, first["start"] - 101), min(len(original.bases), first["end"] + 100))
    sequence = original.bases[interval.start : interval.end]
    fasta = tmp_path / "two records.fasta"
    fasta.write_text(f">forward.1\n{sequence}\n>reverse.1\n{reverse_complement(sequence)}\n")
    report = run(fasta, tmp_path / "multi", settings)
    assert report["summary"]["sequence_count"] == 2
    forward, reverse = report["records"]
    assert forward["predictions"] and reverse["predictions"]
    assert forward["predictions"][0]["sequence"] == reverse["predictions"][0]["sequence"]
    assert forward["predictions"][0]["strand"] != reverse["predictions"][0]["strand"]


def test_real_complete_pipeline(tmp_path):
    settings = load_settings(Path(CONFIG), {"mode": "complete"})
    report = run(ROOT / "examples/complete.fasta", tmp_path / "complete", settings)
    record = report["records"][0]
    assert record["arrays"]
    assert record["cas"]["status"] == "computed"
    assert record["predictions"]
    assert any(p["interaction"] for p in record["predictions"])
    assert any(p["model_hit"] for p in record["predictions"])
    assert any(p["terminator"] for p in record["predictions"])
    assert all(p["interval"]["length"] == len(p["sequence"]) for p in record["predictions"])
    with (tmp_path / "complete" / "summary.csv").open() as handle:
        exported = list(csv.DictReader(handle))
    assert len(exported) == len(record["predictions"])
    for row, prediction in zip(exported, record["predictions"], strict=True):
        assert row["candidate_origin"] == prediction["candidate_origin"]
        assert row["boundary_kind"] == prediction["boundary_kind"]
        assert row["prediction_id"] == prediction["export_id"]
    original = next(read_fasta(ROOT / "examples/complete.fasta")).bases
    for prediction in record["predictions"]:
        interval = prediction["interval"]
        claimed = original[interval["start"] - 1 : interval["end"]]
        if prediction["strand"] == "-":
            claimed = reverse_complement(claimed)
        assert claimed == prediction["sequence"]
        for hit in prediction["model_hits"] + prediction["terminator_hits"]:
            local, genomic = hit["interval"], hit["genomic_interval"]
            local_sequence = claimed[local["start"] - 1 : local["end"]]
            genomic_sequence = original[genomic["start"] - 1 : genomic["end"]]
            if hit["strand"] == "-":
                local_sequence = reverse_complement(local_sequence)
            if hit["genomic_strand"] == "-":
                genomic_sequence = reverse_complement(genomic_sequence)
            assert local_sequence == genomic_sequence


def test_real_type_v_models_and_motif_flag(tmp_path):
    settings = load_settings(
        Path(CONFIG), {"mode": "model", "model_type": "V", "type_v_motifs": True}
    )
    report = run(ROOT / "examples/model-V.fasta", tmp_path / "type-v", settings)
    assert report["summary"]["prediction_count"] > 0
    assert len(report["model_sha256"]) == 4
    assert all(p["type_v_motifs"] is not None for p in report["records"][0]["predictions"])


def test_real_zero_hit_report(tmp_path):
    settings = load_settings(Path(CONFIG), {"mode": "model"})
    fasta = tmp_path / "negative.fasta"
    fasta.write_text(">no_hits\n" + "A" * 200 + "\n")
    output = tmp_path / "negative"
    report = run(fasta, output, settings)
    assert report["summary"]["prediction_count"] == 0
    assert len((output / "summary.csv").read_text().splitlines()) == 1
