"""Generic PWM engineering checks against independent exhaustive enumeration."""

from __future__ import annotations

from dataclasses import replace
from itertools import product
import math

import pytest

from crisprleader.errors import InputError, ReferenceError
from crisprleader.inference.motifs import benjamini_hochberg, scan_motifs
from crisprleader.models import FlankCandidate, RepeatAnchor, Segment

BASES = "ACGT"


def candidate(sequence: str = "AACG", repeat: str = "AT") -> FlankCandidate:
    return FlankCandidate(
        "synthetic-left", "synthetic-source", "synthetic-array", "left", "+",
        (Segment(0, len(sequence)),) if sequence else (), sequence, 500, None,
        RepeatAnchor((Segment(100, 100 + len(repeat)),), repeat), repeat, "CC",
    )


def motif(word: str = "AC", domain: str = "leader") -> dict:
    return {
        "id": "synthetic-motif", "domain": domain,
        "pwm": [{base: (0.7 if base == preferred else 0.1) for base in BASES} for preferred in word],
        "background": {base: 0.25 for base in BASES}, "strands": ["+"],
        "score_scale": 100, "report_pvalue": 1.0,
        "citations": ["synthetic:software-test"], "metadata": {"synthetic": True},
    }


def brute_force_tail(model: dict, target: str, strand: str) -> float:
    """Enumerate actual forward-base strings with an independent implementation."""
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    width = len(model["pwm"])

    def score(word: str | tuple[str, ...]) -> int:
        values = []
        for index, base in enumerate(word):
            probability = model["pwm"][index][base] if strand == "+" else (
                model["pwm"][width - index - 1][complement[base]])
            values.append(round(math.log2(probability / model["background"][base]) * model["score_scale"]))
        return sum(values)

    threshold = score(target)
    masses = []
    for word in product(BASES, repeat=width):
        if score(word) >= threshold:
            masses.append(math.prod(model["background"][base] for base in word))
    return math.fsum(masses)


def test_exact_null_tails_match_exhaustive_enumeration_both_strands() -> None:
    model = motif("ACGT")
    model["background"] = {"A": 0.1, "C": 0.2, "G": 0.3, "T": 0.4}
    model["strands"] = ["+", "-"]
    sequence = "ACTGACGTATTA"
    hits = scan_motifs([candidate(sequence)], [model])
    assert len(hits) == 2 * (len(sequence) - 4 + 1)
    for hit in hits:
        expected = brute_force_tail(model, sequence[hit.start:hit.end], hit.strand)
        assert hit.pvalue == pytest.approx(expected, abs=1e-13)


def test_bh_all_positions_sides_models_even_when_reporting_filtered() -> None:
    first = motif("AC")
    second = motif("TA")
    second["id"] = "other-synthetic-motif"
    candidates = [candidate("AACG"), replace(candidate("TAAA"), id="synthetic-right", side="right", strand="-")]
    all_hits = scan_motifs(candidates, [first, second])
    expected = benjamini_hochberg([hit.pvalue for hit in all_hits])
    assert len(all_hits) == 12
    for hit, qvalue in zip(all_hits, expected, strict=True):
        assert hit.qvalue == pytest.approx(qvalue)
        assert hit.details["testing_scope"]["number_tests"] == 12
    first["report_pvalue"] = second["report_pvalue"] = 0.1
    filtered = scan_motifs(candidates, [first, second])
    lookup = {(h.candidate_id, h.model_id, h.start): h.qvalue for h in all_hits}
    assert filtered
    assert all(hit.pvalue <= 0.1 for hit in filtered)
    assert all(hit.qvalue == lookup[(hit.candidate_id, hit.model_id, hit.start)] for hit in filtered)


def test_ambiguous_positions_are_counted_but_never_reported() -> None:
    hits = scan_motifs([candidate("ACNAC")], [motif("AC")])
    assert len(hits) == 2
    assert all(hit.details["testing_scope"]["number_tests"] == 4 for hit in hits)
    assert all(hit.details["testing_scope"]["number_ambiguous_tests"] == 2 for hit in hits)
    assert all(hit.qvalue == pytest.approx(0.125) for hit in hits)


def test_domains_use_observed_repeat_and_true_spanning_junction() -> None:
    one = candidate("AACG", "TA")
    model = motif("GT", "junction")
    hits = scan_motifs([one], [model])
    assert [(hit.start, hit.end) for hit in hits] == [(3, 5)]
    assert hits[0].details["matched_sequence"] == "GT"
    assert hits[0].pvalue == pytest.approx(1 / 16)
    repeat_model = motif("TA", "repeat")
    repeat_hits = scan_motifs([one], [repeat_model])
    assert [(hit.start, hit.end) for hit in repeat_hits] == [(4, 6)]
    assert repeat_hits[0].details["matched_sequence"] == "TA"
    assert repeat_hits[0].details["matched_sequence"] != one.consensus_repeat
    assert not scan_motifs([one], [motif("A", "junction")])


def test_predeclared_position_range_changes_actual_testing_universe() -> None:
    model = motif("AC")
    model["junction_end_range"] = [-2, -1]
    hits = scan_motifs([candidate("AACG")], [model])
    assert [(hit.start, hit.end) for hit in hits] == [(0, 2), (1, 3)]
    assert hits[0].details["testing_scope"]["number_tests"] == 2


def test_scope_bound_fails_before_returning_partial_significance() -> None:
    with pytest.raises(InputError, match="testing scope has 3 tests"):
        scan_motifs([candidate()], [motif()], max_tests=2)


def test_missing_observed_repeat_cannot_fall_back_to_consensus() -> None:
    with pytest.raises(InputError, match="actual observed"):
        scan_motifs([replace(candidate(), observed_repeat="")], [motif()])


@pytest.mark.parametrize("value", [0.0, -0.1, float("nan"), float("inf"), True])
def test_invalid_matrix_values_rejected(value: float) -> None:
    model = motif()
    model["pwm"][0]["A"] = value
    with pytest.raises(ReferenceError):
        scan_motifs([candidate()], [model])


def test_null_distribution_state_cap_rejects_not_approximates(monkeypatch: pytest.MonkeyPatch) -> None:
    import crisprleader.inference.motifs as module

    module._score_table.cache_clear()
    monkeypatch.setattr(module, "MAX_NULL_STATES", 1)
    with pytest.raises(ReferenceError, match="state bound"):
        scan_motifs([candidate()], [motif()])
    module._score_table.cache_clear()


def test_tied_pvalues_have_equal_adjustments() -> None:
    assert benjamini_hochberg([0.01, 1, 0.01, 0.2]) == pytest.approx([0.02, 1, 0.02, 4 * 0.2 / 3])
    assert benjamini_hochberg([]) == []
    with pytest.raises(ValueError):
        benjamini_hochberg([float("nan")])


def test_duplicate_ids_rejected_instead_of_double_counting() -> None:
    with pytest.raises(InputError, match="duplicate candidate"):
        scan_motifs([candidate(), candidate()], [motif()])
    with pytest.raises(ReferenceError, match="duplicate model"):
        scan_motifs([candidate()], [motif(), motif()])


def test_empty_scopes_are_valid() -> None:
    assert scan_motifs([], [motif()]) == []
    assert scan_motifs([candidate()], []) == []


def test_extreme_null_cannot_silently_underflow_to_zero() -> None:
    model = motif("AC")
    model["background"] = {"A": 1e-200, "C": 0.25, "G": 0.25, "T": 0.5}
    with pytest.raises(ReferenceError, match="underflowed zero"):
        scan_motifs([candidate()], [model])


def test_nonstring_strand_has_actionable_reference_error() -> None:
    model = motif()
    model["strands"] = [["+"]]
    with pytest.raises(ReferenceError, match="strands"):
        scan_motifs([candidate()], [model])


@pytest.mark.parametrize("domain", ["leader", "junction", "repeat"])
def test_mixed_contexts_match_complete_exhaustive_testing_family(domain: str) -> None:
    """Overlapping ambiguity, clipping and both strands retain the entire null family."""
    models = [motif("AC", domain), motif("TAG", domain)]
    for index, model in enumerate(models):
        model["id"] = f"mixed-{index}"
        model["background"] = dict(zip(BASES, (0.1, 0.2, 0.3, 0.4), strict=True))
        model["strands"] = ["-", "+"]
        model["junction_end_range"] = [-4, 4]
    candidates = [candidate("NACNAN", "GTNAC"),
                  replace(candidate("", "ACGT"), id="clipped"),
                  replace(candidate("TACGTA", "TAC"), id="right", side="right", strand="-")]
    pvalues = []
    concrete_tests = []
    ambiguous = 0
    for context in candidates:
        junction = context.sequence + context.observed_repeat
        boundary = len(context.sequence)
        for model in models:
            width = len(model["pwm"])
            for strand in model["strands"]:
                for start in range(len(junction) - width + 1):
                    end = start + width
                    eligible = ((domain == "leader" and end <= boundary)
                                or (domain == "repeat" and start >= boundary)
                                or (domain == "junction" and start < boundary < end))
                    if not eligible or not -4 <= end - boundary <= 4:
                        continue
                    word = junction[start:end]
                    if set(word) <= set(BASES):
                        pvalue = brute_force_tail(model, word, strand)
                        concrete_tests.append((len(pvalues), context.id, model["id"], start, strand, pvalue))
                        pvalues.append(pvalue)
                    else:
                        pvalues.append(1.0)
                        ambiguous += 1
    adjusted = benjamini_hochberg(pvalues)
    hits = scan_motifs(candidates, models)
    assert hits and ambiguous
    assert len(hits) == len(concrete_tests)
    for hit, (index, context_id, model_id, start, strand, pvalue) in zip(hits, concrete_tests, strict=True):
        assert (hit.candidate_id, hit.model_id, hit.start, hit.strand) == (context_id, model_id, start, strand)
        assert hit.pvalue == pytest.approx(pvalue, abs=1e-13)
        assert hit.qvalue == pytest.approx(adjusted[index], abs=1e-13)
        scope = hit.details["testing_scope"]
        assert scope["number_tests"] == len(pvalues)
        assert scope["number_ambiguous_tests"] == ambiguous
