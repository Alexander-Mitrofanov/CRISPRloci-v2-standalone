"""Complete allowlisted repeat typing exports; raw worker output stays private."""
import csv
import json
from ..summaries.common import _json
from ..summaries.repeat_type import public_report


def write_repeat_type_exports(root, summary, checkpoint=None):
    if summary.get("repeat_type", {}).get("available") is not True:
        return
    public = public_report(_json(root, root / "repeat-type/report.json", max_bytes=4_000_000, checkpoint=checkpoint))
    (root / "repeat-typing.json").write_text(json.dumps(public, allow_nan=False) + "\n", encoding="utf-8")
    with (root / "repeat-typing.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        fields = ("id", "sequence", "status", "prediction", "winner", "native_score", "reason")
        writer.writerow(fields)
        for row in public["results"]:
            if checkpoint:
                checkpoint()
            writer.writerow(["'" + row[k] if isinstance(row[k], str) and row[k].lstrip().startswith(("=", "+", "-", "@")) else row[k] for k in fields])
