"""Hand-computed and adversarial fixtures; these contain no biological labels."""
from __future__ import annotations

from copy import deepcopy
import json
import math
import random

import pytest

import crisprleader.evaluation as evaluation
from crisprleader.errors import InputError
from crisprleader.evaluation import evaluate


def case(identifier: str, group: str = "g1", truth: str | None = "+",
         kind: str = "experimental", stratum: str = "synthetic") -> dict:
    return {"id": identifier, "group_id": group, "stratum": stratum, "truth": truth,
            "label_kind": kind, "label_source": "synthetic fixture; no biological claim"}


def prediction(identifier: str, status: str = "called", orientation: str | None = "+") -> dict:
    return {"case_id": identifier, "status": status, "orientation": orientation, "reason": "fixture"}


def documents(cases: list[dict], rows: dict[str, list[dict]]) -> tuple[dict, dict]:
    dataset = {"schema_version": "1.0.0", "dataset_id": "synthetic-comparison",
               "data_kind": "synthetic", "orientation_semantics": "transcription_orientation", "cases": cases}
    predictions = {"schema_version": "1.0.0", "dataset_id": dataset["dataset_id"],
                   "tools": [{"id": identifier, "version": "fixture-v1", "reference_id": None,
                              "predictions": values} for identifier, values in rows.items()]}
    return dataset, predictions


def mixed_documents() -> tuple[dict, dict]:
    cases = [case(str(index), f"g{index // 2}", "+" if index % 2 == 0 else "-",
                  "experimental" if index < 4 else "reviewed_comparative", "A" if index < 3 else "B")
             for index in range(6)]
    cases.extend([case("computational", truth="+", kind="computational"),
                  case("unknown", truth=None, kind="unknown"),
                  case("missing-truth", truth=None, kind="experimental")])
    rows = [prediction("0"), prediction("1"), prediction("2", "abstained", None),
            prediction("3", "unsupported", None), prediction("4", "failed", None),
            prediction("computational"), prediction("unknown"), prediction("missing-truth")]
    return documents(cases, {"tool-a": rows, "tool-b": [prediction(str(index)) for index in range(6)]})


def test_full_denominator_hand_computed_and_unscored_cases_remain_separate() -> None:
    report = evaluate(*mixed_documents(), bootstrap_replicates=2000)
    result = report["tools"][0]
    overall = result["eligible"]["overall"]
    assert overall["denominator"] == 6
    assert overall["counts"] == dict.fromkeys(("correct", "incorrect", "abstained", "unsupported", "failed", "missing"), 1)
    assert overall["coverage"] == pytest.approx(2 / 6)
    assert overall["conditional_accuracy"] == overall["conditional_error"] == 0.5
    assert overall["correct_yield"] == overall["incorrect_yield"] == pytest.approx(1 / 6)
    assert overall["completion"] == pytest.approx(4 / 6)
    assert result["eligible"]["by_stratum"]["A"]["denominator"] == 3
    assert result["eligible"]["by_label_kind"]["experimental"]["denominator"] == 4
    assert result["eligible"]["by_label_kind"]["reviewed_comparative"]["denominator"] == 2
    assert result["eligible"]["by_label_kind"]["computational"]["denominator"] == 0
    assert result["unscored"]["overall"]["case_count"] == 3
    assert result["unscored"]["overall"]["status_counts"]["called_unscored"] == 3
    assert result["unscored"]["overall"]["accuracy"] is None
    unscored = [row for row in report["case_audit"] if not row["eligible"]]
    assert {row["id"] for row in unscored} == {"computational", "unknown", "missing-truth"}
    assert all(row["exclusion_reasons"] for row in unscored)
    assert all(row["outcome"] == "called_unscored" for row in result["case_results"] if not row["eligible"])
    comparison = report["pairwise_comparisons"][0]
    assert comparison["difference_a_minus_b"] == pytest.approx(-1 / 3)
    assert comparison["confidence_interval"]["lower"] == -0.5
    assert comparison["confidence_interval"]["upper"] == 0.0
    assert comparison["all_eligible_predictions_completed"] is False


def test_equal_group_weighting_is_distinct_from_case_weighting() -> None:
    cases = [case("small", "small")] + [case(str(index), "large") for index in range(9)]
    tool_a = [prediction("small")] + [prediction(str(index), orientation="-") for index in range(9)]
    tool_b = [prediction("small", orientation="-")] + [prediction(str(index)) for index in range(9)]
    report = evaluate(*documents(cases, {"a": tool_a, "b": tool_b}), bootstrap_replicates=500)
    assert [tool["eligible"]["overall"]["correct_yield"] for tool in report["tools"]] == [0.1, 0.9]
    assert [tool["primary_group_correct_yield"] for tool in report["tools"]] == [0.5, 0.5]
    comparison = report["pairwise_comparisons"][0]
    assert comparison["difference_a_minus_b"] == 0
    assert comparison["confidence_interval"]["lower"] == -1
    assert comparison["confidence_interval"]["upper"] == 1
    assert comparison["all_eligible_predictions_completed"] is True


def test_paired_resampling_preserves_identical_tool_predictions_exactly() -> None:
    cases = [case(str(index), f"group{index // 3}") for index in range(12)]
    rows = [prediction(str(index), orientation="+" if index % 2 else "-") for index in range(12)]
    report = evaluate(*documents(cases, {"b": rows, "a": deepcopy(rows)}), bootstrap_replicates=2000, seed=7)
    comparison = report["pairwise_comparisons"][0]
    assert comparison["difference_a_minus_b"] == 0
    assert comparison["confidence_interval"]["lower"] == comparison["confidence_interval"]["upper"] == 0


def test_bootstrap_matches_independent_whole_group_oracle() -> None:
    cases = [case("a1", "a"), case("b1", "b"), case("b2", "b"),
             case("c1", "c"), case("c2", "c"), case("c3", "c")]
    a = [prediction("a1"), prediction("b1"), prediction("b2", orientation="-"),
         prediction("c1", "abstained", None), prediction("c2", orientation="-"), prediction("c3", "failed", None)]
    b = [prediction("a1", orientation="-"), prediction("b1"), prediction("b2"), prediction("c1")]
    report = evaluate(*documents(cases, {"a": a, "b": b}), bootstrap_replicates=73, seed=17)
    # Whole-group correct-yield differences: 1, -1/2, -1/3. Each group keeps all cases.
    differences = [1.0, -0.5, -1 / 3]
    rng = random.Random(17)
    draws = sorted(math.fsum(differences[rng.randrange(3)] for _ in range(3)) / 3 for _ in range(73))
    def percentile(probability: float) -> float:
        location = (len(draws) - 1) * probability
        low, high = math.floor(location), math.ceil(location)
        return draws[low] * (high - location) + draws[high] * (location - low) if low != high else draws[low]
    interval = report["pairwise_comparisons"][0]["confidence_interval"]
    assert interval["lower"] == pytest.approx(percentile(0.025))
    assert interval["upper"] == pytest.approx(percentile(0.975))
    assert report["primary_comparison"]["group_ids"] == ["a", "b", "c"]


def test_input_order_does_not_change_output_receipts_or_resamples() -> None:
    dataset, predictions = mixed_documents()
    before = evaluate(dataset, predictions, bootstrap_replicates=47, seed=8)
    dataset["cases"].reverse()
    predictions["tools"].reverse()
    for tool in predictions["tools"]:
        tool["predictions"].reverse()
    after = evaluate(dataset, predictions, bootstrap_replicates=47, seed=8)
    assert after == before
    json.dumps(after, allow_nan=False)


@pytest.mark.parametrize("cases", [[], [case("weak", kind="computational")], [case("none", truth=None)]])
def test_no_eligible_cases_have_nullable_denominators(cases: list[dict]) -> None:
    report = evaluate(*documents(cases, {"a": [], "b": []}))
    for tool in report["tools"]:
        summary = tool["eligible"]["overall"]
        assert summary["denominator"] == 0
        for name in ("coverage", "conditional_accuracy", "conditional_error", "correct_yield", "incorrect_yield", "completion"):
            assert summary[name] is None
        assert tool["primary_group_correct_yield"] is None
    pair = report["pairwise_comparisons"][0]
    assert pair["difference_a_minus_b"] is None
    assert pair["confidence_interval"] is None
    assert pair["interval_unavailable_reason"] == "no_eligible_cases"


def test_all_abstained_or_missing_retains_denominator_and_null_conditional_accuracy() -> None:
    report = evaluate(*documents([case("one"), case("two")], {"a": [prediction("one", "abstained", None)]}))
    summary = report["tools"][0]["eligible"]["overall"]
    assert summary["denominator"] == 2
    assert summary["coverage"] == summary["correct_yield"] == summary["incorrect_yield"] == 0
    assert summary["conditional_accuracy"] is summary["conditional_error"] is None
    assert summary["completion"] == 0.5


def test_one_group_has_point_comparison_but_no_interval() -> None:
    report = evaluate(*documents([case("one"), case("two")], {"a": [prediction("one"), prediction("two")],
             "b": [prediction("one", orientation="-"), prediction("two", orientation="-")]}))
    comparison = report["pairwise_comparisons"][0]
    assert comparison["difference_a_minus_b"] == 1
    assert comparison["confidence_interval"] is None
    assert comparison["eligible_groups"] == 1


def test_zero_replicates_explicitly_disables_intervals() -> None:
    report = evaluate(*mixed_documents(), bootstrap_replicates=0)
    comparison = report["pairwise_comparisons"][0]
    assert comparison["confidence_interval"] is None
    assert comparison["interval_unavailable_reason"] == "bootstrap_disabled"
    assert report["primary_comparison"]["resampling_work_units"] == 0


@pytest.mark.parametrize("data_kind", ["synthetic", "biological"])
def test_declared_data_kind_never_establishes_scientific_superiority(data_kind: str) -> None:
    dataset, predictions = documents([case("a", "a"), case("b", "b")], {"perfect": [prediction("a"), prediction("b")],
         "wrong": [prediction("a", orientation="-"), prediction("b", orientation="-")]})
    dataset["data_kind"] = data_kind
    report = evaluate(dataset, predictions)
    assert report["data_kind"] == data_kind
    assert report["synthetic"] is (data_kind == "synthetic")
    assert report["scientific_superiority_established"] is False
    assert len(report["prerequisites_for_confirmatory_claims"]) == 4
    assert report["pairwise_comparisons"][0]["confidence_interval"]["multiplicity_adjusted"] is False
    assert any("not independently verified" in warning for warning in report["warnings"])


@pytest.mark.parametrize("mutate", [
    lambda d, p: d.update(unknown="extra"),
    lambda d, p: d.update(schema_version="2.0.0"),
    lambda d, p: d.update(data_kind="assumed"),
    lambda d, p: d.update(dataset_id=False),
    lambda d, p: d.update(orientation_semantics="acquisition_orientation"),
    lambda d, p: d["cases"].append(deepcopy(d["cases"][0])),
    lambda d, p: d["cases"][0].update(id=[]),
    lambda d, p: d["cases"][0].update(group_id=""),
    lambda d, p: d["cases"][0].update(stratum=None),
    lambda d, p: d["cases"][0].update(truth=True),
    lambda d, p: d["cases"][0].update(truth=float("nan")),
    lambda d, p: d["cases"][0].update(label_kind="automatic_verified"),
    lambda d, p: d["cases"][0].update(label_source=""),
    lambda d, p: d["cases"][0].update(label_source=float("inf")),
    lambda d, p: p.update(dataset_id="wrong-dataset"),
    lambda d, p: p.update(unknown="extra"),
    lambda d, p: p["tools"].append(deepcopy(p["tools"][0])),
    lambda d, p: p["tools"][0].update(reference_id=False),
    lambda d, p: p["tools"][0].update(version=""),
    lambda d, p: p["tools"][0]["predictions"].append(deepcopy(p["tools"][0]["predictions"][0])),
    lambda d, p: p["tools"][0]["predictions"][0].update(case_id="unknown-case"),
    lambda d, p: p["tools"][0]["predictions"][0].update(status="missing"),
    lambda d, p: p["tools"][0]["predictions"][0].update(orientation=None),
    lambda d, p: p["tools"][0]["predictions"][0].update(status="failed"),
    lambda d, p: p["tools"][0]["predictions"][0].update(reason={}),
    lambda d, p: p["tools"][0]["predictions"][0].update(reason="\ud800"),
])
def test_malformed_or_unbound_inputs_reject_with_input_error(mutate) -> None:
    dataset, predictions = mixed_documents()
    mutate(dataset, predictions)
    with pytest.raises(InputError):
        evaluate(dataset, predictions, bootstrap_replicates=0)


@pytest.mark.parametrize("settings", [
    {"bootstrap_replicates": True}, {"bootstrap_replicates": -1}, {"bootstrap_replicates": 10001},
    {"bootstrap_replicates": float("nan")}, {"seed": False}, {"seed": -1}, {"seed": 2**64},
])
def test_invalid_resampling_settings_reject(settings: dict) -> None:
    with pytest.raises(InputError):
        evaluate(*mixed_documents(), **settings)


@pytest.mark.parametrize("constant,value", [
    ("MAX_CASES", 1), ("MAX_TOOLS", 1), ("MAX_TOOL_CASE_CELLS", 1),
    ("MAX_TEXT_BYTES", 1), ("MAX_RESAMPLING_WORK", 1),
])
def test_declared_resource_bounds_fail_without_silent_subsampling(monkeypatch, constant: str, value: int) -> None:
    monkeypatch.setattr(evaluation, constant, value)
    dataset, predictions = mixed_documents()
    original = deepcopy((dataset, predictions))
    with pytest.raises(InputError):
        evaluate(dataset, predictions)
    assert (dataset, predictions) == original


def test_inputs_are_not_mutated() -> None:
    dataset, predictions = mixed_documents()
    original = deepcopy((dataset, predictions))
    evaluate(dataset, predictions, bootstrap_replicates=31)
    assert (dataset, predictions) == original


@pytest.mark.parametrize("rows", [[], [prediction("one", "failed", None), prediction("two", "failed", None)]])
def test_entirely_missing_or_failed_run_cannot_be_used_as_zero_accuracy_comparator(rows: list[dict]) -> None:
    dataset, predictions = documents(
        [case("one", "g1"), case("two", "g2")],
        {"complete": [prediction("one"), prediction("two")], "unavailable": rows},
    )
    report = evaluate(dataset, predictions)
    unavailable = report["tools"][1]
    assert unavailable["eligible"]["overall"]["denominator"] == 2
    assert unavailable["eligible"]["overall"]["counts"]["missing" if not rows else "failed"] == 2
    assert unavailable["eligible"]["overall"]["conditional_accuracy"] is None
    assert unavailable["primary_group_correct_yield"] is None
    assert unavailable["comparison_eligible"] is False
    comparison = report["pairwise_comparisons"][0]
    assert comparison["difference_a_minus_b"] is None
    assert comparison["confidence_interval"] is None
    assert comparison["comparison_status"] == "not_comparable"
    assert comparison["comparison_unavailable_reason"] == "no_completed_eligible_outcomes_for_tool"
    assert comparison["tools_without_completed_eligible_outcomes"] == ["unavailable"]
    assert report["primary_comparison"]["bootstrap_replicates_computed"] == 0


def test_partial_missing_run_keeps_common_denominator_and_remains_descriptive() -> None:
    report = evaluate(*documents(
        [case("one", "g1"), case("two", "g2")],
        {"complete": [prediction("one"), prediction("two")], "partial": [prediction("one")]},
    ))
    partial = report["tools"][1]
    assert partial["eligible"]["overall"]["denominator"] == 2
    assert partial["eligible"]["overall"]["correct_yield"] == 0.5
    assert partial["comparison_eligible"] is True
    comparison = report["pairwise_comparisons"][0]
    assert comparison["difference_a_minus_b"] == 0.5
    assert comparison["comparison_status"] == "descriptive"
    assert comparison["all_eligible_predictions_completed"] is False


def test_group_weighted_metrics_use_ratios_of_means_with_imbalanced_and_zero_call_groups() -> None:
    cases = ([case("small", "small")] + [case(f"large-{i}", "large") for i in range(4)]
             + [case(f"zero-{i}", "zero-calls") for i in range(2)])
    rows = [prediction("small"), prediction("large-0", orientation="-"),
            prediction("large-1", "abstained", None), prediction("large-2", "failed", None),
            prediction("zero-0", "unsupported", None)]
    report = evaluate(*documents(cases, {"tool": rows}), bootstrap_replicates=0)
    tool = report["tools"][0]
    weighted = tool["group_weighted"]
    assert weighted["group_count"] == 3
    assert weighted["eligible_cases"] == 7
    # Groups: (correct, wrong, call, completion) fractions:
    # small (1,0,1,1); large (0,1/4,1/4,1/2); zero-calls (0,0,0,1/2).
    assert weighted["correct_yield"] == pytest.approx(1 / 3)
    assert weighted["incorrect_yield"] == pytest.approx(1 / 12)
    assert weighted["coverage"] == pytest.approx(5 / 12)
    assert weighted["completion"] == pytest.approx(2 / 3)
    assert weighted["conditional_accuracy"] == pytest.approx(4 / 5)
    assert weighted["conditional_error"] == pytest.approx(1 / 5)
    assert weighted["correct_yield"] == tool["primary_group_correct_yield"]
    # Pooled calls and averaging group accuracies would both give 1/2, a different estimand.
    assert tool["eligible"]["overall"]["conditional_accuracy"] == 0.5
    assert tool["eligible"]["overall"]["correct_yield"] == pytest.approx(1 / 7)
    zero = next(row for row in tool["group_yields"] if row["group_id"] == "zero-calls")
    assert zero["eligible_cases"] == 2
    assert zero["called"] == zero["correct"] == zero["incorrect"] == 0
    assert zero["completed"] == 1
    assert zero["coverage"] == 0
    assert zero["completion"] == 0.5
    assert report["scientific_superiority_established"] is False


def test_all_zero_call_groups_remain_in_group_denominator_without_conditional_accuracy() -> None:
    cases = [case("a", "g1"), case("b", "g2"), case("c", "g2")]
    rows = [prediction("a", "abstained", None), prediction("b", "unsupported", None)]
    report = evaluate(*documents(cases, {"tool": rows}))
    weighted = report["tools"][0]["group_weighted"]
    assert weighted["group_count"] == 2
    assert weighted["eligible_cases"] == 3
    assert weighted["coverage"] == weighted["correct_yield"] == weighted["incorrect_yield"] == 0
    assert weighted["completion"] == 0.75
    assert weighted["conditional_accuracy"] is weighted["conditional_error"] is None


def test_group_weighted_summary_without_eligible_groups_is_nullable() -> None:
    report = evaluate(*documents([case("weak", kind="computational")], {"tool": []}))
    weighted = report["tools"][0]["group_weighted"]
    assert weighted["group_count"] == weighted["eligible_cases"] == 0
    assert all(value is None for name, value in weighted.items() if name not in ("group_count", "eligible_cases"))


def test_group_bookkeeping_does_not_override_unavailable_comparator_gate() -> None:
    report = evaluate(*documents([case("a", "g1"), case("b", "g2")],
                                {"complete": [prediction("a"), prediction("b")], "missing": []}))
    missing = report["tools"][1]
    assert missing["group_weighted"]["correct_yield"] == 0
    assert missing["group_weighted"]["coverage"] == 0
    assert missing["group_weighted"]["completion"] == 0
    assert missing["group_weighted"]["conditional_accuracy"] is None
    assert missing["primary_group_correct_yield"] is None
    assert missing["comparison_eligible"] is False
    assert report["pairwise_comparisons"][0]["difference_a_minus_b"] is None
    assert report["scientific_superiority_established"] is False
