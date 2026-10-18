"""Explicit public JSON and spreadsheet-safe complete viral match exports."""

import csv
import json
from ..summaries.common import _json
from ..summaries.viral import PUBLIC_FIELDS
from .tabular import scientific_cell


def write_viral_exports(root, summary, checkpoint=None):
    if summary.get("viral_search", {}).get("available") is not True:
        return
    doc = _json(
        root, root / "viral/result.json", max_bytes=20_000_000, checkpoint=checkpoint
    )
    public = {key: doc[key] for key in PUBLIC_FIELDS}
    (root / "viral-matches.json").write_text(
        json.dumps(public, allow_nan=False) + "\n", encoding="utf-8"
    )
    fields = [
        "query_id",
        "query_sequence",
        "subject_accession",
        "subject_title",
        "subject_start",
        "subject_end",
        "relative_strand",
        "query_length",
        "identities",
        "substitutions",
        "query_coverage",
        "evalue",
        "bitscore",
    ]
    with (root / "viral-matches.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for hit in public["matches"]:
            if checkpoint:
                checkpoint()
            writer.writerow([scientific_cell(hit[key], strand=key == "relative_strand") for key in fields])
