"""Bounded previews of validated exact spacer matches in submitted target DNA."""

import json

from .common import _json, _safe_file


def protospacer_summary(root, *, max_rows, checkpoint=None):
    path = root / "protospacer/analysis.json"
    if _safe_file(root, path) is None:
        return {"available": False, "matches": []}
    report = _json(root, path, max_bytes=20_000_000, checkpoint=checkpoint)
    if report.get("status") != "completed" or report.get("search_complete") is not True:
        return {"available": False, "matches": []}
    matches, size = [], 0
    for row in report["matches"][:min(max_rows, 1000)]:
        if checkpoint:
            checkpoint()
        size += len(json.dumps(row, allow_nan=False).encode())
        if size > 700_000:
            break
        matches.append(dict(row))
    return {
        "available": True, "search_complete": True, "outcome": report["outcome"],
        "reference": report["reference"], "policy": report["policy"], "counts": report["counts"],
        "matches": matches, "matches_truncated": len(matches) < len(report["matches"]),
        "interpretation": "Exact full-length matches to observed reference spacers are sequence evidence. They do not establish functional targeting, infection, immunity or host range. PAMs are not evaluated and array overlap is unknown. Seed-window counts describe 16-base prefix windows; ambiguous bases elsewhere in a spacer-sized target also prevent a match.",
    }
