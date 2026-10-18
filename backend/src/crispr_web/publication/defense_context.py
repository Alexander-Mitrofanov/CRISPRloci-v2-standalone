"""Complete allowlisted defense context downloads; native work stays private."""

import csv

from ..summaries.common import _json
from ..summaries.defense_context import public_report
from .sequence_ports import write_document
from .context_text import atomic_text


def write_defense_context_exports(root, summary, checkpoint=None):
    if summary.get("defense_context", {}).get("available") is not True:
        return
    doc = public_report(
        _json(
            root,
            root / "defense-context/report.json",
            max_bytes=20_000_000,
            checkpoint=checkpoint,
        ),
        root,
        checkpoint,
    )
    write_document(
        root, "defense-context.json", doc, max_bytes=20_000_000, checkpoint=checkpoint
    )
    with atomic_text(root, "defense-systems.tsv", checkpoint=checkpoint) as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        fields = (
            "id",
            "source_id",
            "type",
            "subtype",
            "start",
            "end",
            "wholeness",
            "score",
            "gene_ids",
        )
        writer.writerow(fields)
        for row in doc["systems"]:
            if checkpoint:
                checkpoint()
            values = [
                ",".join(row[k]) if isinstance(row[k], list) else row[k] for k in fields
            ]
            writer.writerow(
                [
                    "'" + v
                    if isinstance(v, str)
                    and v.lstrip().startswith(("=", "+", "-", "@"))
                    else v
                    for v in values
                ]
            )
