"""Orientation and reconstruction result-table summaries."""

from __future__ import annotations

import csv
import io
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping

from .common import (
    _MAX_SUMMARY_INPUT_BYTES,
    _PUBLIC_RESULT_FIELDS,
    _checkpoint,
    _json,
    _read_bounded,
    _safe_directory,
    _safe_file,
)
from .reconstruction import (
    _mapping,
    _orientation_reconstruction_rows,
    _orientation_tree_rows,
    _reconstruction_tree_rows,
)

def _coerce_csv(value: str | None) -> object:
    if value is None or value == "":
        return None
    text = value.strip()
    try:
        integer = int(text)
        return integer
    except ValueError:
        pass
    try:
        number = float(text)
        return number if math.isfinite(number) else text[:512]
    except ValueError:
        return text[:512]


def _csv_rows(
    root: Path,
    path: Path,
    limit: int,
    checkpoint: Callable[[], None] | None,
    allowed_fields: frozenset[str] = _PUBLIC_RESULT_FIELDS,
    *,
    max_bytes: int = _MAX_SUMMARY_INPUT_BYTES,
) -> tuple[list[dict[str, object]], bool]:
    if _safe_file(root, path) is None:
        return [], False
    payload = _read_bounded(
        root,
        path,
        max_bytes=max_bytes,
        checkpoint=checkpoint,
    ).decode("utf-8")
    rows: list[dict[str, object]] = []
    truncated = False
    with io.StringIO(payload, newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise ValueError("CSV contains missing or duplicate columns")
        selected_fields = [field for field in fields if field in allowed_fields]
        if "name" not in selected_fields:
            raise ValueError("CSV is missing the public group-name column")
        for index, row in enumerate(reader):
            _checkpoint(checkpoint)
            if index >= limit:
                truncated = True
                break
            if None in row:
                raise ValueError("CSV row has fields beyond its header")
            rows.append(
                {
                    str(key)[:240]: _coerce_csv(value)
                    for key, value in row.items()
                    if key is not None and key in selected_fields
                }
            )
    return rows, truncated


def _orientation_summary(
    output_root: Path,
    max_rows: int,
    max_files: int,
    warnings: list[str],
    checkpoint: Callable[[], None] | None,
) -> dict[str, Any] | None:
    root = output_root / "crispr_evor"
    result_path = root / "crispr_evor_results.json"
    csv_path = root / "0_results.csv"
    result_file = _safe_file(output_root, result_path)
    csv_file = _safe_file(output_root, csv_path)
    has_reconstruction_events = any(
        _safe_directory(
            output_root,
            root / orientation / "reconstructed_events",
        )
        is not None
        for orientation in ("0_forward", "0_reversed")
    )
    if result_file is None and csv_file is None and not has_reconstruction_events:
        return None
    document: Mapping[str, Any] = {}
    if result_file is not None:
        try:
            document = _mapping(_json(output_root, result_file, checkpoint=checkpoint))
        except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
            warnings.append(
                f"Could not summarize CRISPR-evOr JSON: {type(error).__name__}"
            )
    comparisons = document.get("results")
    comparisons = comparisons if isinstance(comparisons, list) else []
    safe_comparisons = [
        {
            "group": item.get("group"),
            "prediction": item.get("prediction"),
            "recommended_reverse": item.get("recommended_reverse"),
            "decisive": item.get("decisive"),
            "confidence_threshold": item.get("confidence_threshold"),
            "forward_ln_likelihood_bdm": item.get("forward_ln_likelihood_bdm"),
            "reverse_ln_likelihood_bdm": item.get("reverse_ln_likelihood_bdm"),
            "forward_minus_reverse_ln_likelihood_bdm": item.get(
                "forward_minus_reverse_ln_likelihood_bdm"
            ),
        }
        for item in comparisons[:max_rows]
        if isinstance(item, Mapping)
    ]
    try:
        selected, truncated = _csv_rows(output_root, csv_path, max_rows, checkpoint)
    except (OSError, UnicodeError, ValueError, csv.Error) as error:
        warnings.append(
            f"Could not summarize selected reconstruction CSV: {type(error).__name__}"
        )
        selected, truncated = [], False
    try:
        trees, trees_truncated = _orientation_tree_rows(
            output_root, max_rows, checkpoint
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
        warnings.append(
            f"Could not summarize CRISPR-evOr trees: {type(error).__name__}"
        )
        trees, trees_truncated = [], False
    try:
        reconstructions, reconstructions_truncated = _orientation_reconstruction_rows(
            output_root,
            max_rows,
            max_files,
            warnings,
            checkpoint,
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
        warnings.append(
            "Could not summarize CRISPR-evOr reconstruction events: "
            f"{type(error).__name__}"
        )
        reconstructions, reconstructions_truncated = [], False
    return {
        "status": document.get("status"),
        "tree_policy": document.get("tree_policy"),
        "trees": trees,
        "trees_truncated": trees_truncated,
        "reconstructions": reconstructions,
        "reconstructions_truncated": reconstructions_truncated,
        "comparisons": safe_comparisons,
        "comparisons_truncated": len(comparisons) > len(safe_comparisons),
        # CRISPR-evOr's 0_results.csv contains the reconstruction selected by
        # its forward/reverse likelihood decision and is the definitive table.
        "selected_reconstructions": selected,
        "selected_reconstructions_truncated": truncated,
    }


def _reconstruction_summary(
    output_root: Path,
    max_rows: int,
    warnings: list[str],
    checkpoint: Callable[[], None] | None,
) -> dict[str, Any] | None:
    path = output_root / "spacerplacer" / "0_results.csv"
    if _safe_file(output_root, path) is None:
        return None
    try:
        rows, truncated = _csv_rows(output_root, path, max_rows, checkpoint)
    except (OSError, UnicodeError, ValueError, csv.Error) as error:
        warnings.append(f"Could not summarize SpacerPlacer CSV: {type(error).__name__}")
        return None
    try:
        trees, trees_truncated = _reconstruction_tree_rows(
            output_root, max_rows, checkpoint
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
        warnings.append(
            f"Could not summarize SpacerPlacer trees: {type(error).__name__}"
        )
        trees, trees_truncated = [], False
    return {
        "results": rows,
        "results_truncated": truncated,
        "trees": trees,
        "trees_truncated": trees_truncated,
    }


