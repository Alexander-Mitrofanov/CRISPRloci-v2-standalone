"""Complete context downloads generated from independently validated reports."""

from __future__ import annotations

import csv
import json
from urllib.parse import quote

from ..summaries.leader import leader_summary
from .tabular import scientific_cell


def write_leader_exports(root, summary, checkpoint=None):
    if summary.get("crisprleader", {}).get("status") != "completed":
        return
    evidence = leader_summary(root, max_rows=400_000, checkpoint=checkpoint)
    if evidence["truncated"]:
        raise ValueError("Complete leader export exceeds the supported output contract")
    (root / "leader-contexts.json").write_text(
        json.dumps(
            {
                "schema": {"name": "CRISPRloci-leader-contexts", "version": "1.0.0"},
                "complete": True,
                **evidence,
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    with (
        (root / "leader-contexts.tsv").open("w", newline="") as table,
        (root / "leader-contexts.fasta").open("w") as fasta,
        (root / "leader-contexts.gff3").open("w") as gff,
        (root / "leader-contexts.bed").open("w") as bed,
    ):
        writer = csv.writer(table, delimiter="\t")
        writer.writerow(
            [
                "id",
                "source_id",
                "array_id",
                "side",
                "strand",
                "segments_0based_half_open",
                "available_length",
                "requested_length",
                "truncated_reason",
                "sequence",
                "observed_repeat",
            ]
        )
        gff.write("##gff-version 3\n# Context windows only; no leader prediction.\n")
        bed.write(
            "# Context windows only; zero-based half-open; empty flanks have no BED feature.\n"
        )
        for row in evidence["contexts"]:
            if checkpoint:
                checkpoint()
            cells = [
                row["id"],
                row["source_id"],
                row["array_id"],
                row["side"],
                row["strand"],
                json.dumps(row["segments"], separators=(",", ":")),
                row["available_length"],
                row["requested_length"],
                row["truncated_reason"],
                row["sequence"],
                row["observed_repeat"],
            ]
            writer.writerow([scientific_cell(value, strand=index == 4) for index, value in enumerate(cells)])
            if row["sequence"]:
                fasta.write(
                    f">{row['id']} source={row['source_id']} side={row['side']} strand={row['strand']}\n{row['sequence']}\n"
                )
            for number, segment in enumerate(row["segments"], 1):
                identity = row["id"] + f".part{number}"
                gff.write(
                    "\t".join(
                        [
                            quote(row["source_id"], safe=""),
                            "CRISPRleader-v2",
                            "region",
                            str(segment["start"] + 1),
                            str(segment["end"]),
                            ".",
                            row["strand"],
                            ".",
                            f"ID={quote(identity, safe='')};Name=extracted_context;Note=no_leader_prediction",
                        ]
                    )
                    + "\n"
                )
                bed.write(
                    "\t".join(
                        [
                            row["source_id"],
                            str(segment["start"]),
                            str(segment["end"]),
                            identity,
                            "0",
                            row["strand"],
                        ]
                    )
                    + "\n"
                )
