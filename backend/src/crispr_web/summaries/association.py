"""Bounded, explicitly calibrated spacer association evidence."""

from .common import _json, _safe_file

INTERPRETATION = "Spacer similarity supports a group–target association hypothesis, not active infection, cleavage or a PAM prediction."


def public_report(doc):
    legacy = doc.get("schema_version") == "spacepharer.association/1.0.0"
    reference_fields = ("name", "engine_version", "control_method", "manifest_sha256")
    reference_fields += (
        ("genomes",) if legacy else ("genome_count", "total_bases", "scope", "source")
    )
    rows = [
        {
            key: row[key]
            for key in (
                "group_id",
                "target_id",
                "combined_score",
                "native_fdr",
                "status",
                "native_hit_count",
            )
        }
        | {
            "hits": [
                {
                    key: hit[key]
                    for key in (
                        "occurrence_ids",
                        "query_sequence",
                        "p_best_hit",
                        "query_start",
                        "query_end",
                        "target_start",
                        "target_end",
                        "query_alignment",
                        "target_alignment",
                    )
                }
                for hit in row["hits"]
            ]
        }
        for row in doc["associations"]
    ]
    return {
        "schema": {
            "name": "crisprloci.spacer-association",
            "version": "1.0.0" if legacy else "2.0.0",
        },
        "execution_completed": True,
        "reference": {key: doc["reference"][key] for key in reference_fields},
        "requested_fdr": doc["requested_fdr"],
        "calibration_status": doc["calibration_status"],
        "calibration": doc["calibration"],
        "groups": [
            {
                k: row[k]
                for k in (
                    "id",
                    "members",
                    "supported_occurrence_count",
                    "unique_sequence_count",
                    "unsupported_occurrence_ids",
                )
            }
            for row in doc["groups"]
        ],
        "associations": rows,
        "counts": {
            "records": len(doc["records"]),
            "groups": len(doc["groups"]),
            "associations": sum(row["status"] == "association" for row in rows),
            "candidates": sum(row["status"] != "association" for row in rows),
            "support_hits": sum(len(row["hits"]) for row in rows),
            "unsupported": sum(
                len(row["unsupported_occurrence_ids"]) for row in doc["groups"]
            ),
        },
        "interpretation": INTERPRETATION,
    }


def association_summary(root, *, max_rows, checkpoint=None):
    path = root / "association/report.json"
    if _safe_file(root, path) is None:
        return {"available": False, "associations": []}
    doc = public_report(_json(root, path, max_bytes=20_000_000, checkpoint=checkpoint))
    limit = min(max_rows, 100)
    doc.update(available=True, associations_truncated=len(doc["associations"]) > limit)
    doc["associations"] = doc["associations"][:limit]
    for row in doc["associations"]:
        row["hits_truncated"] = len(row["hits"]) > 50
        row["hits"] = row["hits"][:50]
    return doc
