"""Complete allowlisted repeat mapping exports; native paths stay private."""

import csv
import json
from ..summaries.common import _json
from ..summaries.repeat_map import public_report


def write_repeat_map_exports(root, summary, checkpoint=None):
    if summary.get("repeat_map", {}).get("available") is not True:
        return
    public = public_report(_json(root, root / "repeat-map/report.json", max_bytes=20_000_000, checkpoint=checkpoint))
    (root / "repeat-mapping.json").write_text(json.dumps(public, allow_nan=False) + "\n", encoding="utf-8")
    with (root / "repeat-mapping.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["query_id", "sequence", "status", "best_distance", "best_hit_count",
                         "family_status", "family_value", "motif_status", "motif_value"])
        for row in public["results"]:
            if checkpoint:
                checkpoint()
            fields = [row[k] for k in ("id", "sequence", "status", "best_distance", "best_hit_count")]
            fields += [row["annotations"][f][k] for f in ("family", "motif") for k in ("status", "value")]
            writer.writerow(["'" + v if isinstance(v, str) and v.lstrip().startswith(("=", "+", "-", "@")) else v for v in fields])
