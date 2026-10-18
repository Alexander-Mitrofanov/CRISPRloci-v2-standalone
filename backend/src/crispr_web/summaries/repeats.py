"""Bounded repeat structure previews from independently validated reports."""

import json

from .common import _json, _safe_file


def repeat_summary(root, *, max_rows, checkpoint=None):
    path = root / "repeat/repeat-analysis.json"
    if _safe_file(root, path) is None:
        return {"status": "not_requested", "instances": []}
    report = _json(root, path, max_bytes=20_000_000, checkpoint=checkpoint)
    contexts = {row["id"]: row for row in report["contexts"]}
    instances = []
    preview_bytes = 0
    for item in report["instances"][: min(max_rows, 20)]:
        if checkpoint:
            checkpoint()
        row = dict(item)
        context = contexts.get(item["context_id"])
        row["context"] = None
        if context:
            offset = item["context_start"]
            end = offset + len(item["rna_sequence"])
            internal, competing = [], []
            for i, j, p in context.get("pairs", []):
                inside_i, inside_j = offset <= i < end, offset <= j < end
                if inside_i and inside_j:
                    internal.append([i - offset, j - offset, p])
                elif inside_i or inside_j:
                    competing.append({"i": i, "j": j, "probability": p})
            competing.sort(
                key=lambda pair: (-pair["probability"], pair["i"], pair["j"])
            )
            row["context"] = {
                "status": context["status"],
                "reason": context.get("reason"),
                "effective_window": context["window"],
                "effective_span": context["span"],
                "intrarepeat_pairs": internal,
                "unpaired_probabilities": context.get("unpaired", [])[offset:end],
                "competitors": competing[:20],
                "competitor_count": len(competing),
                "competitors_truncated": len(competing) > 20,
                "pair_coordinates": "zero-based oriented array indices",
            }
        row_bytes = len(json.dumps(row, separators=(",", ":"), allow_nan=False).encode())
        if preview_bytes + row_bytes > 700_000:
            break
        preview_bytes += row_bytes
        instances.append(row)
    return {
        "status": "completed",
        "model": report["model"],
        "counts": report["counts"],
        "classification": report["classification"],
        "instances": instances,
        "orientation_policy": report["orientation_policy"],
        "molecule": report["molecule"],
        "truncated": len(instances) < len(report["instances"]),
        "interpretation": "Predicted structural support under isolated and local-window ensembles. Differences are descriptive, not cleavage, activity or confidence scores. DNA strand alternatives are not ranked. Context stops at array boundaries; window settings can materially change support.",
    }
