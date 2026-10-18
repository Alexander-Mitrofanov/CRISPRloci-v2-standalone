"""Bounded observed-array sharing preview, separate from complete protected exports."""

from .common import _json, _safe_file

ORIGIN_SEMANTICS = "Uploaded source descriptors, origin coordinates, categories, model scores and boundary annotations are client-supplied metadata; source genomes and claimed producers were not independently verified in this comparison."


def public_report(root, checkpoint=None):
    report = _json(
        root,
        root / "array-compare/comparison.json",
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    mapping = _json(
        root,
        root / "array-compare-input/map.json",
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    source = mapping["document"]
    lookup = {(a["source_id"], a["id"]): a for a in source["arrays"]}
    arrays = []
    for row in report["arrays"]:
        original = lookup[row["source_id"], row["array_id"]]
        arrays.append(
            {
                **row,
                **{
                    key: original[key]
                    for key in ("category", "model_score", "interval", "boundary")
                },
                "spacers": [
                    {
                        **s,
                        **{
                            k: original["spacers"][i][k]
                            for k in ("start", "end", "boundary_position")
                            if k in original["spacers"][i]
                        },
                    }
                    for i, s in enumerate(row["spacers"])
                ],
            }
        )
    return {
        "schema": {"name": "crisprloci.array-comparison", "version": "1.0.0"},
        "execution_completed": True,
        "tool": report["tool"],
        "options": report["options"],
        "submitted_sha256": mapping["submitted_sha256"],
        "coordinate_system": source["coordinate_system"],
        "score_semantics": source["score_semantics"],
        "sources": source["sources"],
        "arrays": arrays,
        "edges": report["edges"],
        "components": report["components"],
        "isolates": report["isolates"],
        "counts": {
            "arrays": len(arrays),
            "spacer_occurrences": sum(a["observed_spacer_count"] for a in arrays),
            "deletions": sum(
                len(a["spacers"]) - a["observed_spacer_count"] for a in arrays
            ),
            "edges": len(report["edges"]),
            "components": len(report["components"]),
            "isolates": len(report["isolates"]),
        },
        "interpretation": report["interpretation"],
        "origin_semantics": ORIGIN_SEMANTICS,
    }


def array_compare_summary(root, *, max_rows, checkpoint=None):
    if _safe_file(root, root / "array-compare/comparison.json") is None:
        return {"available": False, "arrays": [], "edges": []}
    doc = public_report(root, checkpoint)
    limit = min(max_rows, 50)
    doc.update(available=True, arrays_truncated=len(doc["arrays"]) > limit)
    doc["arrays"] = doc["arrays"][:limit]
    shown = {a["id"] for a in doc["arrays"]}
    edges = [e for e in doc["edges"] if e["a"] in shown and e["b"] in shown]
    doc["edges_truncated"] = len(doc["edges"]) > min(len(edges), 100)
    doc["edges"] = edges[:100]
    for edge in doc["edges"]:
        edge["shared_symbols_truncated"] = len(edge["shared_symbols"]) > 50
        edge["shared_symbols"] = edge["shared_symbols"][:50]
    for array in doc["arrays"]:
        array["spacers_truncated"] = len(array["spacers"]) > 50
        array["spacers"] = array["spacers"][:50]
    return doc
