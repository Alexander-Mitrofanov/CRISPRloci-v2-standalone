"""Complete capability-protected association evidence exports."""

import csv
import json
from ..summaries.common import _json
from ..summaries.association import public_report


def write_association_exports(root, summary, checkpoint=None):
    if summary.get("spacer_association", {}).get("available") is not True:
        return
    doc = public_report(
        _json(
            root,
            root / "association/report.json",
            max_bytes=20_000_000,
            checkpoint=checkpoint,
        )
    )
    (root / "spacer-associations.json").write_text(
        json.dumps(doc, allow_nan=False) + "\n", encoding="utf8"
    )
    fields = (
        "group_id",
        "target_id",
        "status",
        "combined_score",
        "native_fdr",
        "native_hit_count",
    )
    with (root / "spacer-associations.tsv").open(
        "w", encoding="utf8", newline=""
    ) as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in doc["associations"]:
            if checkpoint:
                checkpoint()
            writer.writerow(
                [
                    "'" + row[k]
                    if isinstance(row[k], str)
                    and row[k].lstrip().startswith(("=", "+", "-", "@"))
                    else row[k]
                    for k in fields
                ]
            )
