"""Bounded public previews of independently validated viral matches."""

from .common import _json, _safe_file

PUBLIC_FIELDS = (
    "schema",
    "tool_version",
    "execution_completed",
    "retrieval_completeness",
    "outcome",
    "reference",
    "policy",
    "coordinate_system",
    "counts",
    "queries",
    "matches",
    "candidate_viruses",
)


def viral_summary(root, *, max_rows, checkpoint=None):
    path = root / "viral/result.json"
    if _safe_file(root, path) is None:
        return {
            "available": False,
            "matches": [],
            "candidate_viruses": [],
            "queries": [],
        }
    doc = _json(root, path, max_bytes=20_000_000, checkpoint=checkpoint)
    limit = min(max_rows, 100)
    return {
        **{
            k: doc[k]
            for k in PUBLIC_FIELDS
            if k not in {"queries", "matches", "candidate_viruses"}
        },
        "available": True,
        "queries": doc["queries"][:limit],
        "queries_truncated": len(doc["queries"]) > limit,
        "matches": doc["matches"][:limit],
        "matches_truncated": len(doc["matches"]) > limit,
        "candidate_viruses": doc["candidate_viruses"][:limit],
        "candidates_truncated": len(doc["candidate_viruses"]) > limit,
    }
