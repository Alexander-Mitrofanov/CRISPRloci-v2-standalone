from __future__ import annotations

import json
from pathlib import Path

from crispr_web.summaries.detection import _array_view
from crispr_web.summaries.detection_sequences import observed_units


def test_native_array_summary_preserves_every_observed_sequence_and_coordinates():
    report = json.loads((Path(__file__).parent / "fixtures/repeats/identify.json").read_text())
    for array in report["arrays"]:
        view = _array_view(array, report["source"]["id"])
        for kind in ("repeats", "spacers"):
            assert len(view[kind]) == len(array[kind])
            for actual, expected in zip(view[kind], array[kind]):
                assert actual["sequence"] == expected["sequence"]
                assert actual["ordinal"] == expected["index"]
                assert actual["start"] == expected["source_interval"]["start"]
                assert actual["end"] == expected["source_interval"]["end"]
        assert view["sequences_truncated"] is False


def test_missing_and_unsearchable_sequences_are_never_fabricated_or_shortened():
    assert _array_view({"consensus_repeat": "ACGT"}, "record")["repeats"] == []
    rows, truncated = observed_units([
        "acgt", {}, {"sequence": "A" * 1001}, {"kind": "deletion", "sequence": "ACGT"},
        {"sequence": "/tmp/private"},
    ])
    assert [row["sequence"] for row in rows] == ["ACGT", None, None, None, None]
    assert [row["sequence_status"] for row in rows] == ["reported", "not_reported", "too_long", "deletion", "invalid"]
    assert truncated is False
    assert observed_units(["A"] * 1001)[1] is True
    rows, truncated = observed_units(["A" * 1000] * 100)
    assert len(rows) == 50 and truncated
