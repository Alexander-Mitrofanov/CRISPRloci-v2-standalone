"""Bounded ISEScan preview; complete evidence is retained in protected exports."""

from .common import _json, _safe_file


def public_report(root, checkpoint=None):
    source = _json(
        root, root / "isescan/isescan.json", max_bytes=20_000_000, checkpoint=checkpoint
    )
    elements = source["elements"]
    from .context_annotations import cas_neighbors
    return {
        **source,
        "cas_context": cas_neighbors(root, elements, checkpoint),
        "schema": {"name": "crisprloci.insertion-context", "version": "1.0.0"},
        "execution_completed": True,
        "counts": {
            "sources": len(source["sources"]),
            "elements": len(elements),
            "complete": sum(e["native_type"] == "c" for e in elements),
            "partial": sum(e["native_type"] == "p" for e in elements),
            "with_terminal_inverted_repeats": sum(
                e["terminal_inverted_repeats"] is not None for e in elements
            ),
        },
    }


def insertion_context_summary(root, *, max_rows, checkpoint=None):
    if _safe_file(root, root / "isescan/isescan.json") is None:
        return {"available": False, "elements": []}
    doc = public_report(root, checkpoint)
    limit = min(max_rows, 100)
    doc.update(available=True, elements_truncated=len(doc["elements"]) > limit)
    doc["elements"] = doc["elements"][:limit]
    doc["cas_context"] = doc["cas_context"][:limit]
    return doc
