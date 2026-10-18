"""Threshold filtering retains the full edit-alignment oracle and reference ties."""

from itertools import product
import json

import pytest

from crisprrepeat import references
from crisprrepeat.inputs import sha
from crisprrepeat.references import References, align, within_distance


def test_banded_threshold_matches_full_alignment_exhaustively():
    sequences = ["".join(items) for n in range(1, 5) for items in product("AC", repeat=n)]
    for query in sequences:
        for target in sequences:
            distance = align(query, target)["edit_distance"]
            for maximum in range(4):
                assert within_distance(query, target, maximum) == (distance <= maximum)


def reference_file(path, sequences):
    path.write_text(json.dumps({
        "schema": {"name": "CRISPRrepeat-reference", "version": "1.0.0"},
        "release": "test", "license": "CC0", "source_url": "https://example.org/fixture",
        "build_date": "2026-10-11",
        "entries": [{"id": f"ref{i}", "sequence": sequence, "sha256": sha(sequence),
                     "occurrences": [{"accession": f"occurrence{i}", "source_url": "https://example.org/fixture",
                                      "evidence": "synthetic control"}]}
                    for i, sequence in enumerate(sequences)],
    }))
    return path


def test_far_reference_matches_do_not_allocate_full_traceback(tmp_path, monkeypatch):
    ref = References(reference_file(tmp_path / "references.json", ["C" * 200]))
    monkeypatch.setattr(references, "align", lambda *a: pytest.fail("Far match allocated full traceback"))
    assert ref.search("A" * 200)["status"] == "no_hit"


def test_duplicate_reference_sequences_share_alignments_and_preserve_occurrences(tmp_path, monkeypatch):
    sequence = "GCGCGAAAACGCGC"
    ref = References(reference_file(tmp_path / "references.json", [sequence, sequence]))
    calls = []
    original = references.align

    def tracked(query, target):
        calls.append(target)
        return original(query, target)

    monkeypatch.setattr(references, "align", tracked)
    result = ref.search(sequence)
    assert len(calls) == len(set(calls))
    assert {hit["reference_id"] for hit in result["hits"]} == {"ref0", "ref1"}
    assert {hit["occurrences"][0]["accession"] for hit in result["hits"]} == {"occurrence0", "occurrence1"}


@pytest.mark.parametrize("maximum", [-1, 4, True, 1.5, float("inf")])
def test_reference_distance_bounds_fail_before_reading(maximum):
    with pytest.raises(ValueError, match="maximum distance"):
        References(None, max_distance=maximum)


def test_indexed_reference_search_matches_full_alignment_oracle_and_report_cap(tmp_path):
    sequences = ["ACGT", "AGT", "CGTA", "AAA", "GCGCGAAAACGCGC"] * 22
    reference = References(reference_file(tmp_path / "references.json", sequences))
    queries = ["".join(items) for n in range(1, 5) for items in product("AC", repeat=n)]
    queries.extend(["ACGU", "GCGCGAAAUCGCGC", "GCGCGAAAACGCGC", "U" * 40])
    for query in queries:
        expected = []
        for entry in reference.entries:
            for orientation, target in (("forward", entry["sequence"]),
                                        ("reverse", references.reverse_complement(entry["sequence"]))):
                aligned = align(query.replace("U", "T"), target)
                if aligned["edit_distance"] <= reference.max_distance:
                    expected.append({"reference_id": entry["id"], "orientation": orientation,
                                     "match": "exact" if aligned["edit_distance"] == 0 else "near",
                                     "occurrences": entry["occurrences"], **aligned})
        expected.sort(key=lambda hit: (hit["edit_distance"], hit["reference_id"], hit["orientation"]))
        result = reference.search(query)
        assert result["total_hit_count"] == len(expected)
        assert result["reported_hit_count"] == min(len(expected), 50)
        assert result["status"] == ("truncated" if len(expected) > 50 else "completed" if expected else "no_hit")
        assert result["hits"] == expected[:50]
