"""Bounded source-spacer-to-prophage context with complete protected evidence."""

from .common import _json, _safe_file


def public_report(root, checkpoint=None):
    report = _json(
        root, root / "phigaro/report.json", max_bytes=20_000_000, checkpoint=checkpoint
    )
    joined = _json(
        root, root / "phigaro/context.json", max_bytes=20_000_000, checkpoint=checkpoint
    )
    return {
        "schema": {"name": "crisprloci.prophage-context", "version": "1.0.0"},
        "execution_completed": True,
        **{
            k: report[k]
            for k in (
                "tool",
                "coordinate_system",
                "input_sha256",
                "reference",
                "parameters",
                "sources",
                "genes",
                "regions",
                "interpretation",
            )
        },
        **joined,
        "counts": {
            "sources": len(report["sources"]),
            "assessed_sources": sum(
                s["status"] == "assessed" for s in report["sources"]
            ),
            "genes": len(report["genes"]),
            "regions": len(report["regions"]),
            "spacer_matches": len(joined["spacer_matches"]),
            "source_arrays": len(joined["source_arrays"]),
            "source_spacers": len(joined["source_spacers"]),
            "matches_overlapping_regions": sum(
                bool(m["region_ids"]) for m in joined["spacer_matches"]
            ),
        },
    }


def prophage_context_summary(root, *, max_rows, checkpoint=None):
    if _safe_file(root, root / "phigaro/report.json") is None:
        return {"available": False, "regions": [], "spacer_matches": []}
    doc = public_report(root, checkpoint)
    limit = max(0, min(max_rows, 100))
    doc.update(
        available=True,
        regions_truncated=len(doc["regions"]) > limit,
        spacer_matches_truncated=len(doc["spacer_matches"]) > limit,
    )
    doc["regions"] = [
        {k: v for k, v in r.items() if k != "sequence"} for r in doc["regions"][:limit]
    ]
    doc["spacer_matches"] = doc["spacer_matches"][:limit]
    query_ids = {row["query_id"] for row in doc["spacer_matches"]}
    # Keep the crosswalk for displayed hits first, even when earlier queries
    # were ineligible. Each inventory preview remains independently bounded.
    spacers = sorted(doc["source_spacers"], key=lambda row: row["id"] not in query_ids)
    doc["source_spacers"] = spacers[:limit]
    array_ids = {row["array_id"] for row in doc["source_spacers"]}
    arrays = sorted(doc["source_arrays"], key=lambda row: row["id"] not in array_ids)
    doc["source_arrays"] = arrays[:limit]
    doc["source_spacers_truncated"] = len(spacers) > limit
    doc["source_arrays_truncated"] = len(arrays) > limit
    doc["source_inventory_complete"] = len(spacers) <= limit and len(arrays) <= limit
    doc.pop("genes")
    return doc
