"""Bounded descriptive orientation evaluation; no predictor or biological labels.

Truth labels and group membership are assertions supplied by the caller. This
module checks their structure, keeps weak/unlabeled cases unscored, and never
interprets a bootstrap interval as proof of scientific superiority.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import itertools
import json
import math
import random
from typing import Any

from .errors import InputError

MAX_CASES = 10_000
MAX_TOOLS = 16
MAX_TOOL_CASE_CELLS = 100_000
MAX_TEXT_BYTES = 16 * 1024 * 1024
MAX_BOOTSTRAP_REPLICATES = 10_000
MAX_RESAMPLING_WORK = 20_000_000
_ELIGIBLE_KINDS = ("experimental", "reviewed_comparative")
_LABEL_KINDS = (*_ELIGIBLE_KINDS, "computational", "unknown")
_STATUSES = ("called", "abstained", "unsupported", "failed")
_OUTCOMES = ("correct", "incorrect", "abstained", "unsupported", "failed", "missing")


@dataclass(frozen=True)
class _Case:
    id: str
    group_id: str
    stratum: str
    truth: str | None
    label_kind: str
    label_source: str

    @property
    def eligible(self) -> bool:
        return self.truth is not None and self.label_kind in _ELIGIBLE_KINDS

    def record(self) -> dict[str, Any]:
        return {"id": self.id, "group_id": self.group_id, "stratum": self.stratum,
                "truth": self.truth, "label_kind": self.label_kind, "label_source": self.label_source}


@dataclass(frozen=True)
class _Prediction:
    case_id: str
    status: str
    orientation: str | None
    reason: str

    def record(self) -> dict[str, Any]:
        return {"case_id": self.case_id, "status": self.status,
                "orientation": self.orientation, "reason": self.reason}


@dataclass(frozen=True)
class _Tool:
    id: str
    version: str
    reference_id: str | None
    predictions: dict[str, _Prediction]

    def record(self) -> dict[str, Any]:
        return {"id": self.id, "version": self.version, "reference_id": self.reference_id,
                "predictions": [self.predictions[key].record() for key in sorted(self.predictions)]}


class _TextBudget:
    def __init__(self) -> None:
        self.bytes = 0

    def text(self, value: Any, name: str, limit: int = 256, *, empty: bool = False) -> str:
        if not isinstance(value, str) or len(value) > limit or (not empty and not value.strip()):
            raise InputError(f"{name} must be a {'possibly empty' if empty else 'nonempty'} string of at most {limit} characters")
        try:
            self.bytes += len(value.encode("utf-8"))
        except UnicodeError as error:
            raise InputError(f"{name} contains invalid Unicode") from error
        if self.bytes > MAX_TEXT_BYTES:
            raise InputError(f"Evaluation inputs exceed the {MAX_TEXT_BYTES}-byte text budget")
        return value


def _object(value: Any, fields: set[str], name: str) -> dict[str, Any]:
    if not isinstance(value, dict) or any(not isinstance(key, str) for key in value):
        raise InputError(f"{name} must be an object with string keys")
    missing, unknown = fields - value.keys(), value.keys() - fields
    if missing or unknown:
        raise InputError(f"{name}: missing fields {sorted(missing)}, unknown fields {sorted(unknown)}")
    return value


def _list(value: Any, limit: int, name: str) -> list[Any]:
    if not isinstance(value, list) or len(value) > limit:
        raise InputError(f"{name} must be a list with at most {limit} entries")
    return value


def _enum(value: Any, values: tuple[str, ...], name: str) -> str:
    if not isinstance(value, str) or value not in values:
        raise InputError(f"{name} must be one of {', '.join(values)}")
    return value


def _strand(value: Any, name: str) -> str | None:
    if value is None:
        return None
    return _enum(value, ("+", "-"), name)


def _parse(dataset: Any, predictions: Any) -> tuple[str, str, list[_Case], list[_Tool]]:
    budget = _TextBudget()
    raw = _object(dataset, {"schema_version", "dataset_id", "data_kind", "orientation_semantics", "cases"}, "dataset")
    if raw["schema_version"] != "1.0.0":
        raise InputError("Unsupported evaluation dataset schema_version")
    if raw["orientation_semantics"] != "transcription_orientation":
        raise InputError("orientation_semantics must be transcription_orientation; leader-side labels are not interchangeable")
    dataset_id = budget.text(raw["dataset_id"], "dataset_id")
    data_kind = _enum(raw["data_kind"], ("synthetic", "biological"), "data_kind")
    cases: list[_Case] = []
    case_ids: set[str] = set()
    for index, row in enumerate(_list(raw["cases"], MAX_CASES, "cases")):
        row = _object(row, {"id", "group_id", "stratum", "truth", "label_kind", "label_source"}, f"case {index}")
        identifier = budget.text(row["id"], "case id")
        if identifier in case_ids:
            raise InputError(f"Duplicate dataset case ID: {identifier}")
        case_ids.add(identifier)
        truth = _strand(row["truth"], "truth")
        label_kind = _enum(row["label_kind"], _LABEL_KINDS, "label_kind")
        cases.append(_Case(identifier, budget.text(row["group_id"], "group_id"),
                           budget.text(row["stratum"], "stratum"), truth, label_kind,
                           budget.text(row["label_source"], "label_source", 4096,
                                       empty=truth is None or label_kind not in _ELIGIBLE_KINDS)))
    raw = _object(predictions, {"schema_version", "dataset_id", "tools"}, "predictions")
    if raw["schema_version"] != "1.0.0":
        raise InputError("Unsupported prediction-table schema_version")
    if budget.text(raw["dataset_id"], "prediction dataset_id") != dataset_id:
        raise InputError("Prediction dataset_id does not match dataset")
    raw_tools = _list(raw["tools"], MAX_TOOLS, "tools")
    if not raw_tools:
        raise InputError("At least one named tool is required")
    if len(raw_tools) * len(cases) > MAX_TOOL_CASE_CELLS:
        raise InputError(f"Evaluation exceeds {MAX_TOOL_CASE_CELLS} tool-by-case cells")
    tools: list[_Tool] = []
    tool_ids: set[str] = set()
    for row in raw_tools:
        row = _object(row, {"id", "version", "reference_id", "predictions"}, "tool")
        identifier = budget.text(row["id"], "tool id")
        if identifier in tool_ids:
            raise InputError(f"Duplicate tool ID: {identifier}")
        tool_ids.add(identifier)
        version = budget.text(row["version"], "tool version")
        reference_id = None if row["reference_id"] is None else budget.text(row["reference_id"], "reference_id")
        rows: dict[str, _Prediction] = {}
        for prediction in _list(row["predictions"], len(cases), "tool predictions"):
            prediction = _object(prediction, {"case_id", "status", "orientation", "reason"}, "prediction")
            case_id = budget.text(prediction["case_id"], "prediction case_id")
            if case_id not in case_ids:
                raise InputError(f"Prediction refers to unknown dataset case: {case_id}")
            if case_id in rows:
                raise InputError(f"Duplicate prediction for tool {identifier}, case {case_id}")
            status = _enum(prediction["status"], _STATUSES, "prediction status")
            orientation = _strand(prediction["orientation"], "prediction orientation")
            if (status == "called") != (orientation is not None):
                raise InputError("Called predictions require +/- orientation; all other statuses require null")
            rows[case_id] = _Prediction(case_id, status, orientation,
                                        budget.text(prediction["reason"], "prediction reason", 4096, empty=True))
        tools.append(_Tool(identifier, version, reference_id, rows))
    return dataset_id, data_kind, sorted(cases, key=lambda case: case.id), sorted(tools, key=lambda tool: tool.id)


def _outcome(case: _Case, tool: _Tool) -> str:
    prediction = tool.predictions.get(case.id)
    if prediction is None:
        return "missing"
    if prediction.status != "called":
        return prediction.status
    if not case.eligible:
        return "called_unscored"
    return "correct" if prediction.orientation == case.truth else "incorrect"


def _ratio(numerator: int, denominator: int) -> float | None:
    return numerator / denominator if denominator else None


def _eligible_summary(cases: list[_Case], tool: _Tool) -> dict[str, Any]:
    counts = dict.fromkeys(_OUTCOMES, 0)
    for case in cases:
        counts[_outcome(case, tool)] += 1
    denominator = len(cases)
    called = counts["correct"] + counts["incorrect"]
    completed = denominator - counts["failed"] - counts["missing"]
    return {
        "denominator": denominator, "counts": counts,
        "coverage": _ratio(called, denominator),
        "conditional_accuracy": _ratio(counts["correct"], called),
        "conditional_error": _ratio(counts["incorrect"], called),
        "correct_yield": _ratio(counts["correct"], denominator),
        "incorrect_yield": _ratio(counts["incorrect"], denominator),
        "completion": _ratio(completed, denominator),
    }


def _unscored_summary(cases: list[_Case], tool: _Tool) -> dict[str, Any]:
    counts = dict.fromkeys(("called_unscored", "abstained", "unsupported", "failed", "missing"), 0)
    for case in cases:
        counts[_outcome(case, tool)] += 1
    return {"case_count": len(cases), "status_counts": counts,
            "accuracy": None, "reason": "Not eligible for ground-truth scoring"}


def _group_weighted(groups: list[dict[str, Any]]) -> dict[str, Any]:
    """Equal-group fractions, with conditional metrics formed as ratios of means."""
    fields = ("coverage", "correct_yield", "incorrect_yield", "completion")
    base = {"group_count": len(groups),
            "eligible_cases": sum(row["eligible_cases"] for row in groups)}
    if not groups:
        return {**base, **dict.fromkeys((*fields, "conditional_accuracy", "conditional_error"))}
    means = {name: math.fsum(row[name] for row in groups) / len(groups) for name in fields}
    coverage = means["coverage"]
    return {**base, **means,
            "conditional_accuracy": means["correct_yield"] / coverage if coverage else None,
            "conditional_error": means["incorrect_yield"] / coverage if coverage else None}


def _quantile(sorted_values: list[float], probability: float) -> float:
    position = (len(sorted_values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (position - lower)


def _comparisons(
    cases: list[_Case], tools: list[_Tool], replicates: int, seed: int,
) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]], int]:
    grouped: dict[str, list[_Case]] = {}
    for case in cases:
        grouped.setdefault(case.group_id, []).append(case)
    group_ids = sorted(grouped)
    yields: list[list[float]] = []
    group_records: dict[str, list[dict[str, Any]]] = {}
    for tool in tools:
        rows: list[dict[str, Any]] = []
        for group_id in group_ids:
            group_cases = grouped[group_id]
            summary = _eligible_summary(group_cases, tool)
            counts = summary["counts"]
            rows.append({
                "group_id": group_id, "eligible_cases": len(group_cases),
                "correct": counts["correct"], "incorrect": counts["incorrect"],
                "called": counts["correct"] + counts["incorrect"],
                "completed": len(group_cases) - counts["failed"] - counts["missing"],
                **{name: summary[name] for name in
                   ("coverage", "correct_yield", "incorrect_yield", "completion")},
            })
        group_records[tool.id] = rows
        yields.append([row["correct_yield"] for row in rows])
    pairs = list(itertools.combinations(range(len(tools)), 2))
    size = len(group_ids)
    completed = [sum(_outcome(case, tool) not in ("failed", "missing") for case in cases)
                 for tool in tools]
    comparable = [completed[first] > 0 and completed[second] > 0 for first, second in pairs]
    work = replicates * (size * (len(tools) + 1) + len(pairs)) if size >= 2 and any(comparable) else 0
    if work > MAX_RESAMPLING_WORK:
        raise InputError(
            f"Requested paired bootstrap work ({work}) exceeds {MAX_RESAMPLING_WORK}; "
            "reduce bootstrap_replicates explicitly (zero disables intervals). No cases or groups were dropped."
        )
    samples: list[list[float]] = [[] for _ in pairs]
    if size >= 2 and replicates and any(comparable):
        rng = random.Random(seed)
        for _ in range(replicates):
            weights = Counter(rng.randrange(size) for _ in range(size))
            means = [math.fsum(values[index] * count for index, count in weights.items()) / size
                     for values in yields]
            for sample, (first, second), allowed in zip(samples, pairs, comparable):
                if allowed:
                    sample.append(means[first] - means[second])
    results = []
    for sample, (first, second), allowed in zip(samples, pairs, comparable):
        interval = None
        unavailable = None
        comparison_reason = ("no_eligible_cases" if not cases else
                             "no_completed_eligible_outcomes_for_tool" if not allowed else None)
        if comparison_reason is not None:
            unavailable = comparison_reason
        elif size < 2:
            unavailable = "fewer_than_two_eligible_groups"
        elif replicates == 0:
            unavailable = "bootstrap_disabled"
        else:
            sample.sort()
            interval = {"level": 0.95, "lower": _quantile(sample, 0.025),
                        "upper": _quantile(sample, 0.975),
                        "method": "paired_whole_group_percentile_bootstrap", "multiplicity_adjusted": False}
        difference = ((math.fsum(yields[first]) - math.fsum(yields[second])) / size
                      if size and allowed else None)
        results.append({
            "tool_a": tools[first].id, "tool_b": tools[second].id,
            "difference_a_minus_b": difference,
            "comparison_status": "descriptive" if comparison_reason is None else "not_comparable",
            "comparison_unavailable_reason": comparison_reason,
            "tools_without_completed_eligible_outcomes": [tools[index].id for index in (first, second)
                                                           if completed[index] == 0],
            "eligible_cases": len(cases), "eligible_groups": size,
            "confidence_interval": interval, "interval_unavailable_reason": unavailable,
            "all_eligible_predictions_completed": all(
                _outcome(case, tool) not in ("failed", "missing")
                for tool in (tools[first], tools[second]) for case in cases
            ) if cases else None,
        })
    return results, group_records, work


def evaluate(
    dataset: dict[str, Any], predictions: dict[str, Any], *,
    bootstrap_replicates: int = 2000, seed: int = 0,
) -> dict[str, Any]:
    """Return deterministic descriptive metrics on a common eligible denominator.

    ``bootstrap_replicates=0`` explicitly disables intervals. Both arguments must
    be plain integers; bounded work failures never trigger silent subsampling.
    """
    if type(bootstrap_replicates) is not int or not 0 <= bootstrap_replicates <= MAX_BOOTSTRAP_REPLICATES:
        raise InputError(f"bootstrap_replicates must be an integer between 0 and {MAX_BOOTSTRAP_REPLICATES}")
    if type(seed) is not int or not 0 <= seed < 2**64:
        raise InputError("seed must be an integer in [0, 2**64)")
    dataset_id, data_kind, cases, tools = _parse(dataset, predictions)
    eligible = [case for case in cases if case.eligible]
    unscored = [case for case in cases if not case.eligible]
    strata = sorted({case.stratum for case in cases})
    comparisons, group_records, work = _comparisons(eligible, tools, bootstrap_replicates, seed)
    eligible_strata: dict[str, list[_Case]] = {stratum: [] for stratum in strata}
    unscored_strata: dict[str, list[_Case]] = {stratum: [] for stratum in strata}
    eligible_labels: dict[str, list[_Case]] = {kind: [] for kind in _LABEL_KINDS}
    unscored_labels: dict[str, list[_Case]] = {kind: [] for kind in _LABEL_KINDS}
    for case in cases:
        (eligible_strata if case.eligible else unscored_strata)[case.stratum].append(case)
        (eligible_labels if case.eligible else unscored_labels)[case.label_kind].append(case)
    tool_rows = []
    for tool in tools:
        groups = group_records[tool.id]
        has_completed_eligible = any(_outcome(case, tool) not in ("failed", "missing") for case in eligible)
        tool_rows.append({
            "id": tool.id, "version": tool.version, "reference_id": tool.reference_id,
            "eligible": {
                "overall": _eligible_summary(eligible, tool),
                "by_stratum": {stratum: _eligible_summary(eligible_strata[stratum], tool)
                               for stratum in strata},
                "by_label_kind": {kind: _eligible_summary(eligible_labels[kind], tool)
                                  for kind in _LABEL_KINDS},
            },
            "unscored": {
                "overall": _unscored_summary(unscored, tool),
                "by_stratum": {stratum: _unscored_summary(unscored_strata[stratum], tool)
                               for stratum in strata},
                "by_label_kind": {kind: _unscored_summary(unscored_labels[kind], tool)
                                  for kind in _LABEL_KINDS},
            },
            "comparison_eligible": has_completed_eligible,
            "primary_group_correct_yield": (math.fsum(row["correct_yield"] for row in groups) / len(groups)
                                            if groups and has_completed_eligible else None),
            "group_yields": groups,
            "group_weighted": _group_weighted(groups),
            "case_results": [{"case_id": case.id, "eligible": case.eligible,
                              "outcome": _outcome(case, tool),
                              "status": tool.predictions[case.id].status if case.id in tool.predictions else "missing",
                              "orientation": tool.predictions[case.id].orientation if case.id in tool.predictions else None,
                              "reason": tool.predictions[case.id].reason if case.id in tool.predictions else "No prediction supplied"}
                             for case in cases],
        })
    canonical_dataset = {"schema_version": "1.0.0", "dataset_id": dataset_id,
                         "data_kind": data_kind, "orientation_semantics": "transcription_orientation",
                         "cases": [case.record() for case in cases]}
    canonical_predictions = {"schema_version": "1.0.0", "dataset_id": dataset_id,
                             "tools": [tool.record() for tool in tools]}
    warnings = [
        "Descriptive pairwise intervals are unadjusted for multiple comparisons and cannot establish superiority.",
        "Label sources, data_kind, group membership and absence of training/test leakage are supplied assertions, not independently verified.",
        "Transcription orientation labels do not establish leader boundaries or biological activity.",
        "Submitted tool tables do not verify completeness against an external comparator registry; entirely missing/failed eligible runs are not comparable.",
    ]
    if data_kind == "synthetic":
        warnings.append("Synthetic data measure reporting behavior only, not biological performance.")
    if 0 < bootstrap_replicates < 1000:
        warnings.append("Fewer than 1000 bootstrap replicates: interval endpoints have substantial Monte Carlo uncertainty.")
    return {
        "schema_version": "1.0.0", "report_kind": "orientation_benchmark_descriptive",
        "dataset_id": dataset_id, "data_kind": data_kind, "synthetic": data_kind == "synthetic",
        "orientation_semantics": "transcription_orientation",
        "scientific_superiority_established": False,
        "eligibility": {"total_cases": len(cases), "eligible_cases": len(eligible),
                        "unscored_cases": len(unscored), "eligible_label_kinds": list(_ELIGIBLE_KINDS),
                        "requires_known_truth": True},
        "case_audit": [{**case.record(), "eligible": case.eligible,
                        "exclusion_reasons": ([] if case.eligible else
                            (["unknown_truth"] if case.truth is None else []) +
                            (["label_kind_not_eligible"] if case.label_kind not in _ELIGIBLE_KINDS else []))}
                       for case in cases],
        "tools": tool_rows, "pairwise_comparisons": comparisons,
        "primary_comparison": {
            "estimand": "difference_in_equal_group_mean_correct_call_yield",
            "direction": "tool_a_minus_tool_b", "resampling_unit": "whole_group_id",
            "group_ids": sorted({case.group_id for case in eligible}),
            "bootstrap_replicates": bootstrap_replicates, "seed": seed,
            "bootstrap_replicates_computed": bootstrap_replicates if work else 0,
            "resampling_work_units": work,
            "percentile_interpolation": "linear_at_p_times_n_minus_one",
            "identical_resamples_used_for_all_tools_and_pairs": True,
            "multiplicity_adjusted": False,
        },
        "input_receipts": {
            "hash_scope": "canonical_semantic_JSON_sorted_case_and_tool_IDs",
            "dataset_sha256": hashlib.sha256(json.dumps(canonical_dataset, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest(),
            "predictions_sha256": hashlib.sha256(json.dumps(canonical_predictions, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest(),
        },
        "prerequisites_for_confirmatory_claims": [
            "Curated, independently defensible orientation labels with traceable primary evidence.",
            "A reviewed grouping and training/calibration/test leakage audit across all tools and references.",
            "Preregistered estimands, effect/noninferiority margins, multiplicity policy and power or precision targets.",
            "Actual complete comparator runs with frozen versions, references, configurations and failure accounting.",
        ],
        "warnings": warnings,
    }
