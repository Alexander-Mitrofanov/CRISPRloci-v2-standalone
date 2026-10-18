"""Generate complete repeat downloads from independently validated JSON."""

import csv
import json

from ..summaries.common import _json


def write_repeat_exports(root, summary, checkpoint=None):
    if summary.get("repeats", {}).get("status") != "completed":
        return
    report = _json(
        root,
        root / "repeat/repeat-analysis.json",
        max_bytes=20_000_000,
        checkpoint=checkpoint,
    )
    (root / "repeat-analysis.json").write_text(
        json.dumps(report, allow_nan=False) + "\n"
    )
    with (
        (root / "repeat-instances.tsv").open("w", newline="") as handle,
        (root / "repeat-analysis.fasta").open("w") as fasta,
        (root / "repeat-structures.txt").open("w") as structures,
    ):
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "id",
                "source_id",
                "array_id",
                "start_0based",
                "end_exclusive",
                "orientation",
                "rna_sequence",
                "mfe_kcal_mol",
                "isolated_pair_support",
                "context_pair_support",
                "support_change",
                "reason",
            ]
        )
        for row in report["instances"]:
            if checkpoint:
                checkpoint()
            isolated, compared = row["isolated"], row["comparison"]
            cells = [
                row["id"],
                row["source_id"],
                row["array_id"],
                row["start"],
                row["end"],
                row["orientation"],
                row["rna_sequence"],
                isolated.get("mfe_kcal_mol"),
                compared["isolated_pair_support"],
                compared["context_pair_support"],
                compared["support_change"],
                compared["reason"],
            ]
            writer.writerow(
                [
                    "'" + cell
                    if isinstance(cell, str) and cell.startswith(("=", "+", "-", "@"))
                    else cell
                    for cell in cells
                ]
            )
            fasta.write(
                f">{row['id']} source={row['source_id']} orientation={row['orientation']}\n{row['rna_sequence']}\n"
            )
            structures.write(
                f">{row['id']} isolated_MFE\n{row['rna_sequence']}\n{isolated.get('structure', 'unavailable')}\n"
            )
