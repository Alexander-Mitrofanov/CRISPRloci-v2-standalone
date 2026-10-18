"""Complete public exports regenerated from the independently validated report."""

import csv
import json

from ..summaries.common import _json
from .tabular import scientific_cell


def write_protospacer_exports(root, summary, checkpoint=None):
    if summary.get("protospacer", {}).get("available") is not True:
        return
    report = _json(root, root / "protospacer/analysis.json", max_bytes=20_000_000, checkpoint=checkpoint)
    # Only the public evidence contract is published; retained targets and native
    # command/configuration artifacts stay in the private worker output.
    public = {key: report[key] for key in ("schema", "tool", "status", "search_complete", "outcome", "reference", "policy", "targets", "counts", "matches")}
    (root / "protospacer-analysis.json").write_text(json.dumps(public, allow_nan=False) + "\n", encoding="utf-8")
    fields = ["spacer_id", "spacer_sequence", "spacer_sha256", "target_id", "target_sha256", "target_start", "target_end", "relative_strand", "identities", "query_length", "substitutions", "gaps", "query_coverage", "array_overlap", "pam_status", "confidence", "evalue"]
    with (root / "protospacer-matches.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for row in public["matches"]:
            if checkpoint:
                checkpoint()
            writer.writerow([scientific_cell(row[name], strand=name == "relative_strand") for name in fields])
