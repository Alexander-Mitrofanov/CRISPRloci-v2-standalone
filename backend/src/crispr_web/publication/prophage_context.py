"""Sanitized prophage region and same-record spacer exports."""

import csv
from ..summaries.prophage_context import public_report
from .sequence_ports import write_document
from .context_text import atomic_text
from .tabular import scientific_cell


def write_prophage_context_exports(root, summary, checkpoint=None):
    if summary.get("prophage_context", {}).get("available") is not True:
        return
    doc = public_report(root, checkpoint)
    write_document(
        root, "prophage-context.json", doc, max_bytes=20_000_000, checkpoint=checkpoint
    )
    for filename, key, fields in (
        (
            "prophage-regions.tsv",
            "regions",
            (
                "id",
                "source_id",
                "start",
                "end",
                "gene_ids",
                "supporting_gene_ids",
                "peak_score",
            ),
        ),
        (
            "prophage-spacer-matches.tsv",
            "spacer_matches",
            (
                "query_id",
                "source_id",
                "array_id",
                "spacer_index",
                "start",
                "end",
                "strand",
                "region_ids",
                "match_class",
            ),
        ),
    ):
        with atomic_text(root, filename, checkpoint=checkpoint) as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(fields)
            for row in doc[key]:
                if checkpoint:
                    checkpoint()
                values = [
                    ",".join(row[k]) if isinstance(row[k], list) else row[k]
                    for k in fields
                ]
                writer.writerow([scientific_cell(value, strand=field == "strand") for field, value in zip(fields, values, strict=True)])
