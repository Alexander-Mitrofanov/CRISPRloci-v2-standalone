"""Independent validation of Integration and scientific-tool outputs."""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path
import re
from typing import Any, Callable, Mapping

from .files import _read_bounded_output


_SUPPORTED_ADAPTER_SCHEMA_MAJOR = 1
_SUPPORTED_IDENTIFY_SCHEMA_MAJOR = 1
_MAX_RESULT_GROUPS = 10_000
_MAX_RESULT_COLUMNS = 256
_SEMVER = re.compile(
    r"(?P<major>0|[1-9][0-9]*)[.]"
    r"(?P<minor>0|[1-9][0-9]*)[.]"
    r"(?P<patch>0|[1-9][0-9]*)"
    r"(?:-(?:0|[1-9][0-9]*|[0-9]*[A-Za-z-][0-9A-Za-z-]*)"
    r"(?:[.](?:0|[1-9][0-9]*|[0-9]*[A-Za-z-][0-9A-Za-z-]*))*)?"
    r"(?:[+][0-9A-Za-z-]+(?:[.][0-9A-Za-z-]+)*)?"
)


def _schema_version_is(document: Mapping[str, Any], name: str, version: str) -> bool:
    schema = document.get("schema")
    return (
        isinstance(schema, Mapping)
        and schema.get("name") == name
        and schema.get("version") == version
    )


def _supported_schema_major(
    document: Mapping[str, Any], name: str, supported_major: int
) -> bool:
    schema = document.get("schema")
    if not isinstance(schema, Mapping) or schema.get("name") != name:
        return False
    version = schema.get("version")
    if not isinstance(version, str):
        return False
    match = _SEMVER.fullmatch(version)
    return match is not None and int(match.group("major")) == supported_major


def _nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _finite_number(value: object) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(value)
    )


def _read_result_csv(
    root: Path,
    path: Path,
    *,
    required_fields: set[str],
    max_rows: int,
    checkpoint: Callable[[], None] | None = None,
) -> dict[str, dict[str, str]] | None:
    payload = _read_bounded_output(root, path, checkpoint=checkpoint)
    if payload is None:
        return None
    try:
        with io.StringIO(payload.decode("utf-8"), newline="") as stream:
            reader = csv.DictReader(stream)
            fieldnames = reader.fieldnames or []
            if (
                not fieldnames
                or fieldnames[0] != "name"
                or len(fieldnames) != len(set(fieldnames))
                or required_fields.difference(fieldnames)
                or len(fieldnames) > _MAX_RESULT_COLUMNS
            ):
                return None
            rows: dict[str, dict[str, str]] = {}
            for row in reader:
                if checkpoint is not None:
                    checkpoint()
                if len(rows) >= max_rows or None in row:
                    return None
                group = (row.get("name") or "").strip()
                if not group or group in rows:
                    return None
                rows[group] = {key: value or "" for key, value in row.items()}
    except (UnicodeError, csv.Error):
        return None
    return rows


def _csv_finite_float(row: Mapping[str, str], field: str) -> float | None:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _csv_nonnegative_float(row: Mapping[str, str], field: str) -> bool:
    value = _csv_finite_float(row, field)
    return value is not None and value >= 0


def _csv_nonnegative_integer(row: Mapping[str, str], field: str) -> bool:
    value = _csv_finite_float(row, field)
    return value is not None and value >= 0 and value.is_integer()


