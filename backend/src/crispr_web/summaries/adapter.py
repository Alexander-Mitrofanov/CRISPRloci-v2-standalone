"""Adapter-membership and eligibility summary extraction."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
from typing import Any, Callable, Mapping

from .common import (
    _json,
    _safe_file,
)
from .reconstruction import (
    _mapping,
)

def _adapter_summary(
    output_root: Path,
    warnings: list[str],
    max_rows: int,
    checkpoint: Callable[[], None] | None,
) -> dict[str, Any] | None:
    path = output_root / "adapter" / "manifest.json"
    if _safe_file(output_root, path) is None:
        return None
    try:
        document = _mapping(_json(output_root, path, checkpoint=checkpoint))
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
        warnings.append(f"Could not summarize adapter manifest: {type(error).__name__}")
        return None
    skipped = document.get("skipped_arrays")
    reasons: Counter[str] = Counter()
    skipped_count = len(skipped) if isinstance(skipped, list) else 0
    if isinstance(skipped, list):
        for item in skipped[:max_rows]:
            if isinstance(item, Mapping):
                reasons[str(item.get("reason") or "unspecified")[:240]] += 1
    unknown_orientation = sum(
        count
        for reason, count in reasons.items()
        if reason.lower() == "orientation_not_determined"
        or (
            "unknown" in reason.lower()
            and ("orientation" in reason.lower() or "strand" in reason.lower())
        )
    )
    groups = []
    raw_groups = document.get("groups")
    group_count = len(raw_groups) if isinstance(raw_groups, list) else 0
    if isinstance(raw_groups, list):
        for group in raw_groups[:max_rows]:
            if not isinstance(group, Mapping):
                continue
            raw_arrays = group.get("arrays")
            member_count = len(raw_arrays) if isinstance(raw_arrays, list) else 0
            members = []
            if isinstance(raw_arrays, list):
                for array in raw_arrays[:max_rows]:
                    if not isinstance(array, Mapping):
                        continue
                    members.append(
                        {
                            "source_id": str(array.get("source_id") or "unknown")[:240],
                            "array_id": str(array.get("array_id") or "unknown")[:240],
                            "category": str(array.get("category") or "unknown")[:80],
                            "spacer_count": array.get("spacer_count"),
                            "strand": str(array.get("strand") or "unknown")[:32],
                            "input_sequence_orientation": str(
                                array.get("input_sequence_orientation") or "unknown"
                            )[:32],
                            "ccdb_strand": str(array.get("ccdb_strand") or "unknown")[
                                :32
                            ],
                        }
                    )
            groups.append(
                {
                    "name": str(group.get("name") or "unknown")[:240],
                    "array_count": group.get("array_count"),
                    "repeat_key": (
                        str(group.get("repeat_key"))[:512]
                        if group.get("repeat_key") is not None
                        else None
                    ),
                    "arrays": members,
                    "arrays_truncated": member_count > len(members),
                }
            )
    return {
        "input_array_count": document.get("input_array_count"),
        "selected_array_count": document.get("selected_array_count"),
        "emitted_array_count": document.get("emitted_array_count"),
        "emitted_group_count": document.get("emitted_group_count"),
        "skipped_array_count": skipped_count,
        "skipped_by_reason": dict(sorted(reasons.items())),
        "skipped_by_reason_truncated": skipped_count > sum(reasons.values()),
        "unknown_strand_excluded_count": unknown_orientation,
        "groups": groups,
        "groups_truncated": group_count > len(groups),
    }


