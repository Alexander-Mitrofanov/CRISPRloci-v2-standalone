"""Sanitized, source-bound insertion context; native work directories stay private."""

import csv
from ..summaries.insertion_context import public_report
from .sequence_ports import write_document
from .context_text import atomic_text


def write_insertion_context_exports(root, summary, checkpoint=None):
    if summary.get("insertion_context", {}).get("available") is not True:
        return
    doc = public_report(root, checkpoint)
    write_document(
        root, "insertion-context.json", doc, max_bytes=20_000_000, checkpoint=checkpoint
    )

    def safe(value):
        return (
            "'" + value
            if isinstance(value, str)
            and value.lstrip().startswith(("=", "+", "-", "@"))
            else value
        )

    with atomic_text(root, "insertion-elements.tsv", checkpoint=checkpoint) as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "id",
                "source_id",
                "start",
                "end",
                "family",
                "cluster",
                "native_type",
                "transposase_strand",
                "evalue_this_copy",
                "evalue_best_copy",
                "terminal_inverted_repeats",
                "target_site_duplication_status",
            ]
        )
        for element in doc["elements"]:
            if checkpoint:
                checkpoint()
            writer.writerow(
                [
                    safe(element[key])
                    for key in (
                        "id",
                        "source_id",
                        "start",
                        "end",
                        "family",
                        "cluster",
                        "native_type",
                    )
                ]
                + [
                    (element["transposase_orf"] or {}).get("strand"),
                    element["evalue_this_copy"],
                    element["evalue_best_copy"],
                    element["terminal_inverted_repeats"] is not None,
                    element["target_site_duplication_status"],
                ]
            )
