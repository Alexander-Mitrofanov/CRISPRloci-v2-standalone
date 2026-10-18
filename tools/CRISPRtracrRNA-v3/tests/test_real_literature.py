"""Native FASTA36 regressions using independently sourced natural loci."""

import json
import os
import shutil
from pathlib import Path

import pytest

from crisprtracrrna.config import load_settings
from crisprtracrrna.models import Array, Interval
from crisprtracrrna.process import Runner
from crisprtracrrna.search import make_candidates, search_anti_repeats
from crisprtracrrna.sequence import read_fasta, reverse_complement

DATA = Path(__file__).parent / "data" / "literature"
CONFIG = os.environ.get("CRISPRTRACRRNA_TEST_CONFIG")
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not CONFIG, reason="Real-tool config not set"),
]


@pytest.mark.parametrize("label", ["Ain", "Kki", "Lsp1"])
def test_real_fasta_recovers_reference_in_long_directory(label, tmp_path):
    row = next(
        r
        for r in json.loads((DATA / "fixtures.json").read_text())["records"]
        if r["label"] == label
    )
    # Lsp1 also exercises the unmasked short-filename path.
    settings = load_settings(
        Path(CONFIG),
        {"threads": 1, "search_backend": "fasta36", "mask_arrays": label != "Lsp1"},
    )
    work = tmp_path / ("long output path with spaces " * 5)
    work.mkdir()
    fasta = work / "input.fasta"
    shutil.copyfile(DATA / row["fasta"], fasta)
    original_bytes = fasta.read_bytes()
    assert len(str(fasta)) > 120
    source = next(read_fasta(fasta))
    arrays = [
        Array(
            a["id"],
            Interval.from_public(a["interval"]["start"], a["interval"]["end"]),
            a["repeat"],
            a["strand"],
            a["category"],
            a["score"],
        )
        for a in row["arrays"]
    ]
    runner = Runner(work / "commands", settings.timeout, settings.environment)
    hits = search_anti_repeats(source, fasta, arrays, work, settings, runner)
    candidates = make_candidates(source, hits, settings)
    truth = row["truth"]
    interval = Interval.from_public(truth["start"], truth["end"])
    assert any(
        candidate.anti.strand == truth["strand"]
        and candidate.interval.overlap(interval) >= 0.9 * interval.length
        for candidate in candidates
    )
    assert fasta.read_bytes() == original_bytes
    assert len(runner.records) == 1
    assert runner.records[0]["argv"][-2] == "repeats.fasta"
    assert all(hit.interval.overlap(a.interval) == 0 for hit in hits for a in arrays)
    for candidate in candidates:
        sequence = source.bases[candidate.interval.start : candidate.interval.end]
        if candidate.anti.strand == "-":
            sequence = reverse_complement(sequence)
        assert candidate.sequence == sequence
    if settings.mask_arrays:
        masked = next(read_fasta(work / "anti-search.fasta"))
        for array in arrays:
            assert masked.bases[array.interval.start : array.interval.end] == (
                "N" * array.interval.length
            )
