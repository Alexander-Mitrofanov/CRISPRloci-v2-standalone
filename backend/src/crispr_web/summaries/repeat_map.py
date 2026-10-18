"""Bounded allowlisted repeat-reference evidence, preserving unknown labels."""

from .common import _json, _safe_file


def public_report(doc):
    result_fields = ("id", "sequence", "sequence_sha256", "status", "reason", "best_distance",
                     "best_hit_count", "reported_hit_count", "hits_truncated")
    hit_fields = ("reference_id", "orientation", "aligned_orientation", "distance", "similarity",
                  "family", "motif", "accession", "query_aligned", "reference_aligned", "alignment_tie_policy")
    annotation_fields = ("status", "values", "value", "missing_reference_count")
    return {
        "schema": {"name": "crisprloci.repeat-mapping", "version": "1.0.0"},
        "execution_completed": True,
        "reference": {k: doc["reference"][k] for k in ("namespace", "content_sha256", "file_sha256", "entry_count", "eligible_entry_count")},
        "policy": {"max_distance": doc["policy"]["max_distance"], "top_k": doc["policy"]["top_k"],
                   "orientation_scope": "relative_to_reference_not_transcription"},
        "counts": {k: doc["counts"][k] for k in ("query_count", "unique_sequence_count", "exact_count", "near_count", "no_hit_count", "unsupported_count")},
        "results": [{**{k: row[k] for k in result_fields},
                     "annotations": {f: {k: row["annotations"][f][k] for k in annotation_fields} for f in ("family", "motif")},
                     "hits": [{k: hit[k] for k in hit_fields} for hit in row["hits"]]} for row in doc["results"]],
    }


def repeat_map_summary(root, *, max_rows, checkpoint=None):
    path = root / "repeat-map/report.json"
    if _safe_file(root, path) is None:
        return {"available": False, "results": []}
    doc = public_report(_json(root, path, max_bytes=20_000_000, checkpoint=checkpoint))
    limit = min(max_rows, 100)
    doc.update(available=True, results_truncated=len(doc["results"]) > limit)
    doc["results"] = doc["results"][:limit]
    return doc
