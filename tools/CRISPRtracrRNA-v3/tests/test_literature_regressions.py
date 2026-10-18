"""Failures observed in the literature benchmark, with deterministic contracts."""

import csv
import hashlib
import json
from pathlib import Path

import pytest

from crisprtracrrna.config import Settings
from crisprtracrrna.errors import TracrError
from crisprtracrrna.evidence import model_rescue_predictions
from crisprtracrrna.models import Array, Hit, Interval, Sequence
from crisprtracrrna.pipeline import write_exports
from crisprtracrrna.search import search_anti_repeats
from crisprtracrrna.sequence import read_fasta, reverse_complement

DATA = Path(__file__).parent / "data" / "literature"


def test_literature_fixtures_reconstruct_exact_reference_sequences():
    for row in json.loads((DATA / "fixtures.json").read_text())["records"]:
        fasta = DATA / row["fasta"]
        assert hashlib.sha256(fasta.read_bytes()).hexdigest() == row["file_sha256"]
        source = next(read_fasta(fasta))
        assert source.sha256 == row["sequence_sha256"]
        truth = row["truth"]
        assert source.id == truth["contig"]
        found = source.bases[truth["start"] - 1 : truth["end"]]
        if truth["strand"] == "-":
            found = reverse_complement(found)
        assert found == row["tracr_sequence"]


def test_fasta_uses_short_names_in_long_directory(tmp_path):
    work = tmp_path / ("long output with spaces " * 5)
    work.mkdir()
    fasta = work / "input.fasta"
    fasta.write_text(">s\n" + "ACGT" * 100 + "\n")
    source = Sequence("s", "s", "ACGT" * 100)
    arrays = [Array("a", Interval(200, 300), "ACGT" * 9, "+", "Bona-fide", 1.0)]

    class LegacyPathLimitRunner:
        calls = 0

        def run(self, command, *, cwd, stdout):
            self.calls += 1
            for name in command[-2:]:
                assert len(name) < 120
                assert (cwd / name).is_file()
            assert command[-2] == "repeats.fasta"
            assert cwd == work
            stdout.write_text("")

    runner = LegacyPathLimitRunner()
    assert search_anti_repeats(source, fasta, arrays, work, Settings(), runner) == []
    assert runner.calls == 1
    masked = next(read_fasta(work / "anti-search.fasta")).bases
    assert masked[200:300] == "N" * 100
    assert masked[:200] == source.bases[:200]
    assert masked[300:] == source.bases[300:]
    assert fasta.read_text().splitlines()[1] == source.bases


@pytest.mark.parametrize("strand", ["+", "-"])
def test_model_rescue_without_arrays_preserves_local_and_genomic_coordinates(strand, tmp_path):
    source = Sequence("s", "s", "ACGTTGCA" * 50)
    hits = [
        Hit("s", Interval(20, 75), strand, 40, 1e-20, "a"),
        Hit("s", Interval(35, 80), strand, 30, 1e-12, "b"),
    ]
    rows = model_rescue_predictions(source, hits, [], [], {"cassettes": []}, Settings())
    assert len(rows) == 1
    row = rows[0]
    assert row["candidate_origin"] == "covariance_model"
    assert row["boundary_kind"] == "model_supported_region"
    assert row["array"] is None and row["anti_repeat"] is None
    for hit in row["model_hits"]:
        g = hit["genomic_interval"]
        local = hit["interval"]
        genome = source.bases[g["start"] - 1 : g["end"]]
        if strand == "-":
            genome = reverse_complement(genome)
        assert row["sequence"][local["start"] - 1 : local["end"]] == genome
    assert all(v == 0 for k, v in row["ranking"]["features"].items() if k != "model")
    write_exports(
        tmp_path,
        [
            {
                "key": "s",
                "source": {"file": "input.fasta", "id": "s", "length": 400},
                "predictions": rows,
            }
        ],
    )
    with (tmp_path / "summary.csv").open() as handle:
        exported = next(csv.DictReader(handle))
    assert exported["candidate_origin"] == row["candidate_origin"]
    assert exported["boundary_kind"] == row["boundary_kind"]
    assert exported["strand"] == strand
    assert exported["sequence"] == row["sequence"]
    for name in (
        "array_id",
        "array_category",
        "anti_repeat_identity",
        "anti_repeat_coverage",
        "interaction_energy_kcal_mol",
        "terminator_evalue",
    ):
        assert exported[name] == ""


def test_rescue_excludes_weak_hits_array_hits_and_already_covered_strand():
    source = Sequence("s", "s", "ACGT" * 100)
    array = Array("a", Interval(200, 300), "ACGT", "+", "Bona-fide", 1.0)
    hits = [
        Hit("s", Interval(20, 60), "+", 40, 1e-20),
        Hit("s", Interval(20, 60), "-", 40, 1e-20),
        Hit("s", Interval(100, 150), "+", 10, 0.1),
        Hit("s", Interval(190, 220), "+", 40, 1e-20),
    ]
    existing = [{"interval": Interval(10, 80).public(), "strand": "+"}]
    rows = model_rescue_predictions(source, hits, [array], existing, {"cassettes": []}, Settings())
    assert len(rows) == 1 and rows[0]["strand"] == "-"


@pytest.mark.parametrize("value", [True, 0, -1, float("nan"), float("inf")])
def test_invalid_rescue_threshold_rejected(value):
    with pytest.raises(TracrError):
        Settings(model_rescue_evalue=value)
