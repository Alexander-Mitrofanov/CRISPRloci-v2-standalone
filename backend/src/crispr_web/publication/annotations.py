"""Export all validated annotations, independently of the UI preview limit."""

from __future__ import annotations
import csv
import json

from ..summaries.annotations import annotation_summary
from ..summaries.common import _json, _safe_file


def write_annotation_exports(output_root, summary, checkpoint=None):
    from .leader import write_leader_exports
    write_leader_exports(output_root, summary, checkpoint)
    cas = summary.get("casandra", {})
    tracr = summary.get("tracrrna", {})
    if cas.get("status") != "completed" and tracr.get("status") != "completed":
        return
    warnings = []
    cas, tracr, detection = annotation_summary(
        output_root,
        max_rows=1_000_000,
        accepted_categories=summary.get("detection", {}).get("accepted_categories", []),
        warnings=warnings,
        checkpoint=checkpoint,
    )
    if (
        warnings
        or cas.get("truncated")
        or tracr.get("truncated")
        or (detection or {}).get("arrays_truncated")
    ):
        raise ValueError(
            "Complete annotation export exceeds the supported output contract"
        )
    document = {
        "schema": {"name": "CRISPRloci-annotations", "version": "2.0.0"},
        "coordinate_system": "1-based, end-inclusive, input forward strand",
        "complete": True,
        "casandra": cas,
        "tracrrna": tracr,
        "detection": detection,
    }
    (output_root / "annotations.json").write_text(
        json.dumps(document, indent=2, allow_nan=False) + "\n"
    )
    if tracr.get("mode") == "complete":
        report = _json(
            output_root, output_root / "tracrrna/report.json", checkpoint=checkpoint
        )
        sources = []
        for record in report["records"]:
            arrays = []
            for path in (
                output_root / "tracrrna/work" / record["key"] / "identify"
            ).rglob("report.json"):
                if _safe_file(output_root, path) is not None:
                    arrays.extend(
                        _json(output_root, path, checkpoint=checkpoint)["arrays"]
                    )
            sources.append(
                {
                    "id": record["source"]["id"],
                    "length": record["source"]["length"],
                    "sha256": record["source"]["sha256"],
                    "arrays": arrays,
                }
            )
        (output_root / "array-evidence.json").write_text(
            json.dumps(
                {
                    "schema": {"name": "CRISPRloci-array-evidence", "version": "2.0.0"},
                    "sources": sources,
                    "note": "Native CRISPRidentify array evidence, including categories excluded from the preview.",
                },
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
    for name, rows in [
        ("cas-genes", cas.get("genes", [])),
        ("cas-cassettes", cas.get("cassettes", [])),
        ("tracr-candidates", tracr.get("candidates", [])),
    ]:
        keys = list(rows[0]) if rows else ["id", "source_id", "start", "end", "strand"]
        with (output_root / f"{name}.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=keys, delimiter="\t")
            writer.writeheader()
            for row in rows:
                if checkpoint:
                    checkpoint()
                writer.writerow(row)
    with (output_root / "tracr-predictions.fasta").open("w") as stream:
        for row in tracr.get("candidates", []):
            stream.write(f">{row['id']}\n{row['sequence']}\n")
