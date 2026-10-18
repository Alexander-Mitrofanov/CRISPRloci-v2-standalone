"""Bounded observed array sequences; never substitute consensus for an occurrence."""

from __future__ import annotations

from collections.abc import Mapping
import re

_DNA = re.compile(r"[ACGTRYSWKMBDHVN]+", re.IGNORECASE)


def observed_units(value: object) -> tuple[list[dict], bool]:
    if not isinstance(value, list):
        return [], False
    rows = []
    characters = 0
    for index, raw in enumerate(value[:1000]):
        item = raw if isinstance(raw, Mapping) else {}
        sequence = raw if isinstance(raw, str) else item.get("sequence")
        kind = item.get("kind")
        status = "reported"
        if kind == "deletion":
            sequence, status = None, "deletion"
        elif not isinstance(sequence, str) or not sequence:
            sequence, status = None, "not_reported"
        elif len(sequence) > 1000:
            sequence, status = None, "too_long"
        elif not _DNA.fullmatch(sequence):
            sequence, status = None, "invalid"
        else:
            sequence = sequence.upper()
        characters += len(sequence or "")
        if characters > 50_000:
            break
        interval = item.get("source_interval")
        interval = interval if isinstance(interval, Mapping) else {}
        rows.append({
            "ordinal": item.get("index", index + 1),
            "sequence": sequence,
            "sequence_status": status,
            "start": interval.get("start"),
            "end": interval.get("end"),
        })
    return rows, len(rows) < len(value)
