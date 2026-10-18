"""Protect scientific benchmark denominators and coordinate interpretation."""

import runpy
from pathlib import Path

import pytest

# These are checkout utilities, not installed runtime modules. Loading by their
# absolute paths avoids sys.path changes and requires only the standard library.
SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
family = runpy.run_path(str(SCRIPTS / "summarize_sequence_benchmark.py"))["family"]
localization_match = runpy.run_path(str(SCRIPTS / "summarize_literature_benchmark.py"))[
    "localization_match"
]


def test_model_family_grouping_uses_exact_major_type():
    expected = {
        "II": "II",
        "II-A": "II",
        "II-B": "II",
        "II-C": "II",
        "unknown": "II",  # Explicit exploratory screen, not a subtype assignment.
        "V": "V",
        "V-B": "V",
        "V-F1": "V",
        "V-K_V1": "V",
        "III": None,
        "III-A": None,
        "III-D": None,
        "VI": None,
        "VI-A": None,
        "I-C": None,
        "IV-A1": None,
    }
    assert {subtype: family(subtype) for subtype in expected} == expected


@pytest.mark.parametrize("contig,strand", [("other_contig", "+"), ("reference", "-")])
def test_localization_rejects_wrong_contig_or_strand(contig, strand):
    truth = {"contig": "reference", "start": 101, "end": 200, "strand": "+"}
    prediction = {"interval": {"start": 101, "end": 200}, "strand": strand}
    assert localization_match(prediction, truth, contig) == 0


@pytest.mark.parametrize(
    "start,end,expected_overlap",
    [
        (141, 180, 40),  # A conserved region can localize without full transcript boundaries.
        (80, 230, 100),  # A covering search window localizes but is not an exact boundary call.
        (151, 250, 50),  # Exactly half the shorter interval, with inclusive public coordinates.
        (181, 280, 0),  # Twenty incidental bases are less than half of either interval.
        (122, 140, 0),  # A fully contained 19 nt fragment fails the minimum overlap requirement.
        (201, 260, 0),  # Adjacent intervals do not overlap.
    ],
)
def test_localization_requires_meaningful_reference_overlap(start, end, expected_overlap):
    truth = {"contig": "reference", "start": 101, "end": 200, "strand": "+"}
    prediction = {"interval": {"start": start, "end": end}, "strand": "+"}
    assert localization_match(prediction, truth, "reference") == expected_overlap
