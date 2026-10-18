"""Parse hierarchical native alignments into sequence-bound evidence."""

from __future__ import annotations

import math
from collections.abc import Mapping

from .io import InputError

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def oriented(sequence: str, start: int, end: int) -> str:
    if not 1 <= min(start, end) <= max(start, end) <= len(sequence):
        raise InputError("Native coordinates outside sequence")
    value = sequence[min(start, end) - 1 : max(start, end)]
    return value if start <= end else value.translate(COMPLEMENT)[::-1]


def number(text: str, low: float = 0, high: float = 1e300) -> float:
    value = float(text)
    if not math.isfinite(value) or not low <= value <= high:
        raise InputError("Invalid native scientific number")
    return value


def parse(
    text: str,
    query: dict[str, dict],
    targets: Mapping[str, str],
    native_groups: dict[str, str],
    calibrated: bool,
    cutoff: float = 0.05,
) -> list[dict]:
    lines = text.splitlines()
    rows: list[dict] = []
    current: dict | None = None
    pairs: set[tuple[str, str]] = set()
    support: set[tuple[str, int, int, int, int]] = set()
    i = 0
    total = 0
    while i < len(lines):
        fields = lines[i].split("\t")
        i += 1
        if fields[0].startswith("#"):
            if len(fields) != 5 or fields[0][1:] not in native_groups or fields[1] not in targets:
                raise InputError("Unexpected native association header")
            fdr = number(fields[3], -1, 1)
            if (calibrated and fdr < 0) or (not calibrated and fdr != -1):
                raise InputError("Native FDR disagrees with calibration receipt")
            pair = native_groups[fields[0][1:]], fields[1]
            if pair in pairs:
                raise InputError("Duplicate native association")
            pairs.add(pair)
            support = set()
            current = {
                "group_id": native_groups[fields[0][1:]],
                "target_id": fields[1],
                "combined_score": number(fields[2]),
                "native_fdr": fdr if calibrated and fdr >= 0 else None,
                "status": ("association" if fdr <= cutoff else "above_cutoff_candidate")
                if calibrated and fdr >= 0
                else "uncalibrated_candidate",
                "native_hit_count": int(fields[4]),
                "hits": [],
            }
            if not 1 <= current["native_hit_count"] <= 10000:
                raise InputError("Invalid native support count")
            rows.append(current)
            if len(rows) > 5000:
                raise InputError("Too many associations")
        elif fields[0].startswith(">") and current is not None:
            if (
                len(fields) != 7
                or i + 1 >= len(lines)
                or fields[0][1:] not in query
                or fields[1] != current["target_id"]
            ):
                raise InputError("Unexpected native support row")
            record = query[fields[0][1:]]
            if record["group_id"] != current["group_id"]:
                raise InputError("Native support crossed query groups")
            qs, qe, ts, te = map(int, fields[3:7])
            qa, ta = lines[i : i + 2]
            i += 2
            if (
                len(qa) != len(ta)
                or not qa
                or len(qa) > 400
                or set(qa + ta) - set("ACGTN-")
                or any(a == b == "-" for a, b in zip(qa, ta))
            ):
                raise InputError("Malformed native alignment")
            if qa.replace("-", "") != oriented(record["sequence"], qs, qe) or ta.replace("-", "") != oriented(
                targets[fields[1]], ts, te
            ):
                raise InputError("Native alignment disagrees with original sequences")
            signature = record["sequence"], qs, qe, ts, te
            if signature in support:
                raise InputError("Duplicate native support alignment")
            support.add(signature)
            current["hits"].append(
                {
                    "occurrence_ids": record["occurrence_ids"],
                    "query_sequence": record["sequence"],
                    "p_best_hit": number(fields[2], 0, 1),
                    "query_start": qs,
                    "query_end": qe,
                    "target_start": ts,
                    "target_end": te,
                    "query_alignment": qa,
                    "target_alignment": ta,
                }
            )
            total += 1
            if total > 10000:
                raise InputError("Too many support alignments")
        else:
            raise InputError("Unexpected native report record")
    for row in rows:
        if not row["hits"] or row["native_hit_count"] != len(row["hits"]):
            raise InputError("Native support count mismatch")
    return rows
