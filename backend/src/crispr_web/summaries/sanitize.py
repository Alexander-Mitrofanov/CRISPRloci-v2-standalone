"""Schema-tolerant, bounded extraction of scientifically useful summaries."""

from __future__ import annotations

import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Mapping

from .contracts import _MAX_PUBLIC_SUMMARY_BYTES


_STATIC_PRIVATE_MARKERS = (
    "/srv",
    "/opt",
    "/etc",
    "/home",
    "/root",
    "/tmp",
    "/var/tmp",
    "/var/lib",
)
_EMBEDDED_ABSOLUTE_PATH = re.compile(
    r"(?<![A-Za-z0-9:/])/(?:[A-Za-z0-9._-]+/)+[A-Za-z0-9._-]+"
)


def _summary_path_markers(private_roots: Iterable[str | Path]) -> tuple[str, ...]:
    markers = set(_STATIC_PRIVATE_MARKERS)
    for root in private_roots:
        rendered = Path(root).expanduser().resolve(strict=False).as_posix().rstrip("/")
        if rendered and rendered != "/":
            markers.add(rendered)
    return tuple(sorted(markers, key=len, reverse=True))


def _path_bearing_string(value: str, markers: tuple[str, ...]) -> bool:
    normalized = value.replace("\\", "/")
    stripped = normalized.strip()
    if (
        stripped.startswith("/")
        or re.match(r"^[A-Za-z]:/", stripped)
        or _EMBEDDED_ABSOLUTE_PATH.search(normalized)
    ):
        return True
    return any(marker in normalized for marker in markers)


def sanitize_public_summary(
    summary: Mapping[str, Any],
    *,
    private_roots: Iterable[str | Path],
    max_depth: int = 20,
    max_items: int = 10_000,
    max_string_characters: int = 10_000,
) -> dict[str, Any]:
    """Recursively bound public scalars and redact path-bearing tool values."""

    markers = _summary_path_markers(private_roots)
    redactions = [0]

    def sanitize(value: object, depth: int) -> object:
        if depth > max_depth:
            redactions[0] += 1
            return None
        if isinstance(value, Mapping):
            cleaned: dict[str, object] = {}
            for index, (raw_key, item) in enumerate(value.items()):
                if index >= max_items:
                    redactions[0] += 1
                    break
                key = str(raw_key)[:240]
                if _path_bearing_string(key, markers) or key in cleaned:
                    redactions[0] += 1
                    continue
                cleaned[key] = sanitize(item, depth + 1)
            return cleaned
        if isinstance(value, (list, tuple)):
            if len(value) > max_items:
                redactions[0] += len(value) - max_items
            return [sanitize(item, depth + 1) for item in value[:max_items]]
        if isinstance(value, str):
            if _path_bearing_string(value, markers):
                redactions[0] += 1
                return None
            if len(value) > max_string_characters:
                redactions[0] += 1
                return value[:max_string_characters]
            return value
        if value is None or isinstance(value, (bool, int)):
            return value
        if isinstance(value, float):
            if math.isfinite(value):
                return value
            redactions[0] += 1
            return None
        redactions[0] += 1
        return None

    sanitized = sanitize(summary, 0)
    if not isinstance(sanitized, dict):
        raise ValueError("summary sanitization did not produce an object")
    if redactions[0]:
        warnings = sanitized.get("warnings")
        if not isinstance(warnings, list):
            warnings = []
            sanitized["warnings"] = warnings
        warnings.append(
            f"{redactions[0]} private, oversized, or unsupported summary value(s) were omitted."
        )
    if (
        len(
            json.dumps(sanitized, separators=(",", ":"), allow_nan=False).encode(
                "utf-8"
            )
        )
        > _MAX_PUBLIC_SUMMARY_BYTES
    ):
        truncation_targets = (
            ("repeats", "instances", "truncated"),
            ("crisprleader", "contexts", "truncated"),
            ("casandra", "genes", "truncated"),
            ("casandra", "cassettes", "truncated"),
            ("tracrrna", "candidates", "truncated"),
            ("detection", "arrays", "arrays_truncated"),
            ("adapter", "groups", "groups_truncated"),
            ("orientation", "comparisons", "comparisons_truncated"),
            (
                "orientation",
                "selected_reconstructions",
                "selected_reconstructions_truncated",
            ),
            ("orientation", "trees", "trees_truncated"),
            (
                "orientation",
                "reconstructions",
                "reconstructions_truncated",
            ),
            ("reconstruction", "results", "results_truncated"),
            ("reconstruction", "trees", "trees_truncated"),
        )
        for section_name, list_name, flag_name in truncation_targets:
            section = sanitized.get(section_name)
            if isinstance(section, dict) and isinstance(section.get(list_name), list):
                section[list_name] = []
                section[flag_name] = True
        warnings = sanitized.get("warnings")
        if not isinstance(warnings, list):
            warnings = []
            sanitized["warnings"] = warnings
        warnings.append(
            "Detailed rows were omitted because the public summary reached its byte limit."
        )
        if (
            len(
                json.dumps(sanitized, separators=(",", ":"), allow_nan=False).encode(
                    "utf-8"
                )
            )
            > _MAX_PUBLIC_SUMMARY_BYTES
        ):
            return {
                "schema_version": sanitized.get("schema_version"),
                "pipeline_status": sanitized.get("pipeline_status"),
                "warnings": [
                    "Public summary details were omitted because they exceeded the byte limit."
                ],
            }
    return sanitized
