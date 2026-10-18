"""Bounded, public context evidence without native paths or inference claims."""

from __future__ import annotations

import re

from .common import _json, _safe_file

LIMITATION = "Context extraction only. Leader prediction is unavailable; neither side is selected and window length is not a biological boundary."


def leader_summary(root, *, max_rows, checkpoint=None):
    if _safe_file(root, root / "leader/index.json") is None:
        return {"status": "not_requested", "contexts": [], "sources": []}
    index = _json(root, root / "leader/index.json", checkpoint=checkpoint)
    summary = {
        "status": "completed",
        "analysis_mode": "context_extraction",
        "prediction_performed": False,
        "coordinate_system": "zero_based_half_open",
        "flank_length": index["flank_length"],
        "array_count": index["array_count"],
        "context_count": index["context_count"],
        "limitation": LIMITATION,
        "contexts": [],
        "sources": [],
        "truncated": False,
    }
    for record in index["records"]:
        if not re.fullmatch(r"record-[0-9]{6}", record["key"]):
            raise ValueError("Invalid leader record identity")
        report = _json(
            root,
            root / "leader" / record["key"] / "result.json",
            checkpoint=checkpoint,
            max_bytes=20_000_000,
        )
        source = report["sources"][0]
        summary["sources"].append(
            {
                "id": record["source_id"],
                "length": source["length"],
                "sha256": source["sha256"],
                "topology": source["topology"],
            }
        )
        arrays = {row["id"]: row for row in report["arrays"]}
        for row in report["candidates"]:
            if checkpoint:
                checkpoint()
            if len(summary["contexts"]) >= max_rows:
                break
            array = arrays[row["array_id"]]
            summary["contexts"].append(
                {
                    "id": row["id"],
                    "source_id": record["source_id"],
                    "array_id": array["metadata"]["identify_original"]["id"],
                    "array_segments": array["segments"],
                    "category": array["category"],
                    "side": row["side"],
                    "strand": row["strand"],
                    "segments": row["segments"],
                    "sequence": row["sequence"],
                    "available_length": len(row["sequence"]),
                    "requested_length": row["requested_length"],
                    "truncated_reason": row["truncated_reason"],
                    "observed_repeat": row["observed_repeat"],
                    "repeat_segments": row["anchor"]["segments"],
                }
            )
    summary["truncated"] = len(summary["contexts"]) < summary["context_count"]
    return summary
