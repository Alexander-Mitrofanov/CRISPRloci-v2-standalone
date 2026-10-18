from itertools import product
import random

import pytest

from crisprmap import core
from crisprmap.core import canonical_sequence, global_alignment, levenshtein, reverse_complement, search
from crisprmap.io import InputError, LimitError
from crisprmap.reference import build_reference


def reference(rows):
    return build_reference(
        [{"id": ident, "sequence": sequence} for ident, sequence, _, _ in rows],
        {ident: {"family": family, "motif": motif} for ident, _, family, motif in rows},
        {"namespace": "synthetic", "source_url": "https://example.org/synthetic", "source_sha256": "0" * 64})


def oracle(a, b):
    row = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        current = [i]
        for j, y in enumerate(b, 1):
            current.append(min(row[j] + 1, current[j - 1] + 1, row[j - 1] + (x != y)))
        row = current
    return row[-1]


def test_myers_exhaustive_and_long_integer_words():
    words = ["".join(bases) for length in range(6) for bases in product("AC", repeat=length)]
    rng = random.Random(20260930)
    pairs = [(a, b) for a in words for b in words]
    pairs += [("".join(rng.choices("ACGT", k=rng.randrange(1, 201))),
               "".join(rng.choices("ACGT", k=rng.randrange(1, 201)))) for _ in range(200)]
    for a, b in pairs:
        expected = oracle(a, b)
        assert levenshtein(a, b) == expected
        aligned = global_alignment(a, b)
        aq, ar = aligned["query_aligned"], aligned["reference_aligned"]
        assert aq.replace("-", "") == a
        assert ar.replace("-", "") == b
        assert sum(x != y for x, y in zip(aq, ar)) == expected


def test_retrieval_against_independent_both_orientation_oracle():
    rng = random.Random(63026)
    rows = [(f"r{i}", "".join(rng.choices("ACGT", k=rng.randrange(4, 13))), f"F{i % 3}", None)
            for i in range(35)]
    ref = reference(rows)
    queries = [{"id": f"q{i}", "sequence": "".join(rng.choices("ACGT", k=rng.randrange(4, 13)))}
               for i in range(90)]
    result = search(queries, ref, max_distance=3, top_k=1000)
    for query, observed in zip(queries, result["results"]):
        scores = {}
        for ident, sequence, _, _ in rows:
            fwd = oracle(query["sequence"], sequence)
            rev = oracle(query["sequence"], reverse_complement(sequence))
            scores[ident] = min(fwd, rev), "both" if fwd == rev else ("forward" if fwd < rev else "reverse")
        best = min(value[0] for value in scores.values())
        expected = {key: value for key, value in scores.items() if value[0] == best and best <= 3}
        assert {hit["reference_id"]: (hit["distance"], hit["orientation"])
                for hit in observed["hits"]} == expected
        assert observed["best_hit_count"] == len(expected)


def test_reverse_complement_duplicates_and_palindromes():
    ref = reference([("a", "AAACC", "F1", "M1"), ("b", "GGTTT", "F1", "M1"),
                     ("p", "ACGT", None, None)])
    out = search([{"id": "q", "sequence": "AAACC"}, {"id": "p", "sequence": "ACGT"}], ref)
    assert [(h["reference_id"], h["orientation"]) for h in out["results"][0]["hits"]] == [
        ("a", "forward"), ("b", "reverse")]
    assert out["results"][1]["hits"][0]["orientation"] == "both"
    assert out["counts"]["distance_comparisons"] == 0
    assert canonical_sequence("AAACC") == canonical_sequence("GGTTT")


def test_ties_are_pooled_before_display_cap():
    ref = reference([("a", "AACCG", "F1", "M1"), ("b", "AACCG", "F2", "M1"),
                     ("c", "AACCG", None, None)])
    result = search([{"id": "q", "sequence": "AACCG"}], ref, top_k=1)["results"][0]
    assert result["best_hit_count"] == 3
    assert result["reported_hit_count"] == 1
    assert result["hits_truncated"] is True
    assert result["annotations"]["family"] == {
        "status": "conflicting", "values": ["F1", "F2"], "value": None, "missing_reference_count": 1}
    assert result["annotations"]["motif"]["status"] == "incomplete"
    assert result["annotations"]["motif"]["value"] is None


def test_near_candidates_do_not_assign_family_or_motif():
    ref = reference([("a", "AACCG", "F1", "M1")])
    result = search([{"id": "q", "sequence": "AACCA"}], ref)["results"][0]
    assert result["status"] == "near"
    assert result["annotations"]["family"]["status"] == "candidate"
    assert result["annotations"]["motif"]["status"] == "reference_evidence"
    assert all(value["value"] is None for value in result["annotations"].values())
    assert result["hits"][0]["similarity"] == 0.8


def test_exact_unanimous_historical_annotations_only():
    ref = reference([("a", "AACCG", "F1", None)])
    result = search([{"id": "q", "sequence": "AACCG"}], ref)["results"][0]
    assert result["annotations"]["family"]["value"] == "F1"
    assert result["annotations"]["family"]["status"] == "historical_association"
    assert result["annotations"]["motif"]["status"] == "unannotated"


def test_missing_and_ambiguous_remain_explicit():
    ref = reference([("a", "AAAA", None, None), ("n", "NNNN", "F1", "M1")])
    out = search([{"id": "unknown", "sequence": "ANAA"}, {"id": "far", "sequence": "CGCG"}], ref,
                 max_distance=0)
    assert [r["status"] for r in out["results"]] == ["unsupported", "no_hit"]
    assert out["reference"]["eligible_entry_count"] == 1
    assert out["results"][1]["best_distance"] is None
    ref = reference([("n", "NNNN", "F1", "M1")])
    assert search([{"id": "q", "sequence": "AAAA"}], ref)["results"][0]["reason"] == (
        "no_unambiguous_reference_sequences")


def test_duplicate_query_occurrences_preserved_and_cached():
    ref = reference([("a", "AACCG", "F1", "M1")])
    rows = [{"id": "one", "sequence": "AACCA"}, {"id": "two", "sequence": "AACCA"}]
    out = search(rows, ref)
    assert [row["id"] for row in out["results"]] == ["one", "two"]
    assert out["counts"]["unique_sequence_count"] == 1
    assert out["counts"]["distance_comparisons"] == 2
    out["results"][0]["hits"][0]["family"] = "changed"
    assert out["results"][1]["hits"][0]["family"] == "F1"


def test_comparison_preflight_and_deadline_fail_instead_of_completing():
    ref = reference([("a", "AACCG", "F1", None)])
    with pytest.raises(LimitError, match="before execution"):
        search([{"id": "q", "sequence": "AAACC"}], ref, max_comparisons=1)
    with pytest.raises(LimitError, match="deadline"):
        search([{"id": "q", "sequence": "AAACC"}], ref, deadline_seconds=1e-20)


def test_report_budget_checked_before_duplicate_result_assembly(monkeypatch):
    ref = reference([("a", "AACCG", "F1", None)])
    monkeypatch.setattr(core, "MAX_REPORT_BYTES", 6000)
    with pytest.raises(LimitError, match="report byte budget"):
        search([{"id": f"q{i}", "sequence": "AACCG"} for i in range(100)], ref)


@pytest.mark.parametrize("options", [{"top_k": 0}, {"top_k": True}, {"max_distance": -1},
                                     {"max_comparisons": 0}, {"deadline_seconds": float("nan")},
                                     {"deadline_seconds": True}])
def test_invalid_policy(options):
    with pytest.raises(InputError):
        search([{"id": "q", "sequence": "AAAA"}], reference([("a", "AAAA", None, None)]), **options)



def test_tied_occurrences_reuse_alignment_without_losing_labels(monkeypatch):
    ref = reference([("a", "AACCG", "F1", "M1"), ("b", "CGGTT", "F2", None),
                     ("c", "AACCG", None, "M2")])
    calls = []
    original = core.global_alignment

    def counted(query, target):
        calls.append((query, target))
        return original(query, target)

    monkeypatch.setattr(core, "global_alignment", counted)
    row = search([{"id": "q", "sequence": "AACCA"}], ref)["results"][0]
    assert calls == [("AACCA", "AACCG")]
    assert [(hit["reference_id"], hit["orientation"]) for hit in row["hits"]] == [
        ("a", "forward"), ("b", "reverse"), ("c", "forward")]
    assert row["annotations"]["family"]["values"] == ["F1", "F2"]
    assert row["annotations"]["motif"]["values"] == ["M1", "M2"]
    assert row["annotations"]["family"]["missing_reference_count"] == 1
    assert all(sum(a != b for a, b in zip(hit["query_aligned"], hit["reference_aligned"])) == 1
               for hit in row["hits"])


@pytest.mark.parametrize("sequence", ["", "ACGT", "A" * 200])
def test_exact_alignment_is_the_unique_zero_edit_traceback(sequence):
    assert global_alignment(sequence, sequence) == {
        "query_aligned": sequence, "reference_aligned": sequence,
        "alignment_tie_policy": "diagonal_deletion_insertion"}
