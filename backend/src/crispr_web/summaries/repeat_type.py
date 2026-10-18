"""Bounded repeat-derived subtype evidence with explicit unknown predictions."""
from .common import _json, _safe_file


def public_report(doc):
    return {"schema": {"name": "crisprloci.repeat-typing", "version": "1.0.0"}, "execution_completed": True,
            "model": {k: doc["model"][k] for k in ("id", "upstream_version", "upstream_commit", "files", "labels", "license")},
            "policy": {k: doc["policy"][k] for k in ("threshold", "score_semantics", "orientation")},
            "counts": {k: doc["counts"][k] for k in ("records", "distinct_classified_sequences", "predicted", "uncertain", "unsupported")},
            "results": [{k: row[k] for k in ("id", "sequence", "sequence_sha256", "status", "prediction", "winner", "native_score", "reason")}
                        for row in doc["results"]]}


def repeat_type_summary(root, *, max_rows, checkpoint=None):
    path = root / "repeat-type/report.json"
    if _safe_file(root, path) is None:
        return {"available": False, "results": []}
    doc = public_report(_json(root, path, max_bytes=4_000_000, checkpoint=checkpoint))
    limit = min(max_rows, 1000)
    doc.update(available=True, results_truncated=len(doc["results"]) > limit)
    doc["results"] = doc["results"][:limit]
    return doc
