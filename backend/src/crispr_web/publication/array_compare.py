"""Complete protected observed-array outputs; native files stay private."""

import csv
from ..summaries.array_compare import public_report
from ..summaries.common import _json
from .sequence_ports import write_document


def write_array_compare_exports(root, summary, checkpoint=None):
    if summary.get("array_compare", {}).get("available") is not True:
        return
    doc = public_report(root, checkpoint)
    write_document(
        root, "array-comparison.json", doc, max_bytes=20_000_000, checkpoint=checkpoint
    )
    mapping = _json(
        root,
        root / "array-compare-input/map.json",
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    write_document(
        root,
        "arrays.json",
        mapping["document"],
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    lookup = {a["id"]: a for a in doc["arrays"]}

    def safe(value):
        return (
            "'" + value
            if isinstance(value, str)
            and value.lstrip().startswith(("=", "+", "-", "@"))
            else value
        )

    with (root / "array-sharing.tsv").open("w", encoding="utf8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["source_a", "array_a", "source_b", "array_b", "shared_spacers", "jaccard"]
        )
        for edge in doc["edges"]:
            if checkpoint:
                checkpoint()
            a, b = lookup[edge["a"]], lookup[edge["b"]]
            writer.writerow(
                [
                    safe(a["source_id"]),
                    safe(a["array_id"]),
                    safe(b["source_id"]),
                    safe(b["array_id"]),
                    edge["nshared"],
                    edge["jaccard"],
                ]
            )
