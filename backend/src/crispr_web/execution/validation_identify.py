"""CRISPRidentify report-to-input binding validation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path

from .files import _read_json, _safe_nonempty_output
from .validation_common import (
    _SUPPORTED_IDENTIFY_SCHEMA_MAJOR,
    _nonempty_string,
    _supported_schema_major,
)


def _validate_identify_reports(
    output_root: Path,
    expected_source_bindings: Mapping[str, tuple[int, str]],
    checkpoint: Callable[[], None] | None = None,
    *, reports_root: Path | None = None,
) -> bool:
    expected_source_ids = tuple(expected_source_bindings)
    expected_sources = set(expected_source_ids)
    if not expected_sources or len(expected_sources) != len(expected_source_ids):
        return False

    observed_sources: set[str] = set()
    report_count = 0
    for report_path in (reports_root or output_root / "crispridentify").rglob("report.json"):
        report_count += 1
        if report_count > len(expected_source_ids):
            return False
        if not _safe_nonempty_output(output_root, report_path):
            return False
        report = _read_json(output_root, report_path, checkpoint=checkpoint)
        source = report.get("source")
        if (
            not _supported_schema_major(
                report,
                "CRISPRidentify-report",
                _SUPPORTED_IDENTIFY_SCHEMA_MAJOR,
            )
            or not isinstance(source, Mapping)
            or not _nonempty_string(source.get("id"))
        ):
            return False
        source_id = str(source["id"])
        binding = expected_source_bindings.get(source_id)
        if (
            binding is None
            or isinstance(source.get("length"), bool)
            or source.get("length") != binding[0]
            or source.get("sha256") != binding[1]
            or source.get("sha256_scope") != "uppercase_sequence"
        ):
            return False
        arrays = report.get("arrays")
        if not isinstance(arrays, list):
            results = report.get("results")
            arrays = (
                results.get("crispr_arrays") if isinstance(results, Mapping) else None
            )
        if (
            not isinstance(arrays, list)
            or source_id not in expected_sources
            or source_id in observed_sources
        ):
            return False
        observed_sources.add(source_id)

    return observed_sources == expected_sources
