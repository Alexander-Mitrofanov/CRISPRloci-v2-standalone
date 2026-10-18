"""Mode-specific detection, adapter, orientation, and pipeline summaries."""

from __future__ import annotations

from collections import Counter
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping

from .common import (
    _json,
    _safe_file,
)
from .reconstruction import (
    _mapping,
)
from .detection_sequences import observed_units

def _finite_score(value: object) -> int | float | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return value
    if isinstance(value, str):
        try:
            parsed = float(value)
        except ValueError:
            return None
        return parsed if math.isfinite(parsed) else None
    return None


def _first(mapping: Mapping[str, Any], names: Iterable[str]) -> object | None:
    for name in names:
        if name in mapping:
            return mapping[name]
    return None


def _report_arrays(document: object) -> list[Mapping[str, Any]]:
    root = _mapping(document)
    arrays = root.get("arrays")
    if isinstance(arrays, list):
        return [item for item in arrays if isinstance(item, Mapping)]
    results = _mapping(root.get("results"))
    arrays = results.get("crispr_arrays")
    if isinstance(arrays, list):
        return [item for item in arrays if isinstance(item, Mapping)]
    return []


def _source_id(document: object, fallback: str) -> str:
    root = _mapping(document)
    source = _mapping(root.get("source"))
    candidate = source.get("id") or root.get("source_id") or fallback
    return str(candidate)[:200]


def _array_view(array: Mapping[str, Any], source_id: str) -> dict[str, Any]:
    orientation = _mapping(_first(array, ("orientation", "direction", "strand")))
    raw_orientation = _first(array, ("orientation", "direction", "strand"))
    if orientation:
        strand = _first(orientation, ("strand", "label", "direction", "status"))
    else:
        strand = raw_orientation
    spacers = _first(array, ("spacers", "spacer_sequences"))
    spacer_count = array.get("spacer_count")
    if not isinstance(spacer_count, int) and isinstance(spacers, list):
        spacer_count = len(spacers)
    score = _finite_score(_first(array, ("model_score", "certainty_score", "score")))
    interval = _mapping(array.get("source_interval"))
    repeat_units, repeats_truncated = observed_units(_first(array, ("repeats", "repeat_sequences")))
    spacer_units, spacers_truncated = observed_units(spacers)
    return {
        "source_id": source_id,
        "array_id": str(_first(array, ("id", "array_id", "name")) or "unknown")[:240],
        "category": str(
            _first(array, ("category", "classification", "category_name")) or "Unknown"
        )[:120],
        # This is the raw classifier/model value.  It is deliberately never
        # rendered as a probability or percentage.
        "model_score": score,
        "model_score_is_probability": False,
        "strand": None if strand is None else str(strand)[:120],
        "spacer_count": spacer_count if isinstance(spacer_count, int) else None,
        "repeat_count": (
            array.get("repeat_count")
            if isinstance(array.get("repeat_count"), int)
            else None
        ),
        "start": interval.get("start"),
        "end": interval.get("end"),
        "repeats": repeat_units,
        "spacers": spacer_units,
        "sequences_truncated": repeats_truncated or spacers_truncated,
    }


def _detection_summary(
    output_root: Path,
    max_arrays: int,
    max_files: int,
    accepted_categories: frozenset[str],
    warnings: list[str],
    checkpoint: Callable[[], None] | None,
) -> dict[str, Any]:
    report_paths: list[Path] = []
    scanned = 0
    for path in (output_root / "crispridentify").rglob("report.json"):
        if checkpoint is not None:
            checkpoint()
        scanned += 1
        if scanned > max_files:
            warnings.append(
                "CRISPRidentify report discovery reached the configured file limit."
            )
            break
        safe = _safe_file(output_root, path)
        if safe is not None:
            report_paths.append(safe)
    report_paths.sort()
    arrays: list[dict[str, Any]] = []
    categories: Counter[str] = Counter()
    selected_categories: Counter[str] = Counter()
    schemas: Counter[str] = Counter()
    source_count = 0
    total_arrays = 0
    selected_array_count = 0
    for path in report_paths:
        if checkpoint is not None:
            checkpoint()
        try:
            document = _json(output_root, path, checkpoint=checkpoint)
            root = _mapping(document)
            schema = _mapping(root.get("schema"))
            version = schema.get("version")
            if version is not None:
                schemas[str(version)] += 1
            source_id = _source_id(document, path.parent.name)
            source_count += 1
            report_arrays = _report_arrays(document)
            total_arrays += len(report_arrays)
            for raw in report_arrays:
                view = _array_view(raw, source_id)
                category = view["category"]
                categories[category] += 1
                selected = category in accepted_categories
                if selected:
                    selected_categories[category] += 1
                    selected_array_count += 1
                    view["selected_for_analysis"] = True
                    if len(arrays) < max_arrays:
                        arrays.append(view)
        except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
            warnings.append(f"Could not summarize {path.name}: {type(error).__name__}")
    return {
        "source_count": source_count,
        "report_count": len(report_paths),
        "array_count": total_arrays,
        "category_counts": dict(sorted(categories.items())),
        "selected_array_count": selected_array_count,
        "unselected_array_count": total_arrays - selected_array_count,
        "selected_category_counts": dict(sorted(selected_categories.items())),
        "accepted_categories": sorted(accepted_categories),
        "arrays": arrays,
        "arrays_scope": "selected_for_analysis",
        "arrays_truncated": selected_array_count > len(arrays),
        "report_schema_versions": dict(sorted(schemas.items())),
        "score_semantics": {
            "primary_interpretation": "category",
            "model_score_is_probability": False,
            "note": "Raw model scores are not calibrated probabilities.",
        },
    }

