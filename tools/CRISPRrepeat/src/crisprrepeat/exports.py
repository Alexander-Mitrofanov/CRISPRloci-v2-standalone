"""Atomic, bounded, sequence-bound JSON and tabular downloads."""

import csv
import json
import os

MAX_OUTPUT_BYTES = 20_000_000


def write_outputs(root, report):
    root.mkdir(parents=True, exist_ok=False)
    temporary = root / ".report.tmp"
    # Shared cached folds and reference occurrences expand at every occurrence
    # in JSON. Enforce the byte budget while encoding, before allocating the
    # complete serialized report or publishing any derivative artifacts.
    encoder = json.JSONEncoder(sort_keys=True, separators=(",", ":"), allow_nan=False)
    total = 0
    try:
        with temporary.open("wb") as stream:
            for chunk in encoder.iterencode(report):
                data = chunk.encode("utf-8")
                total += len(data)
                if total + 1 > MAX_OUTPUT_BYTES:  # Reserve the terminating newline.
                    raise ValueError(
                        "Complete report exceeds output byte limit; no partial scientific result was published"
                    )
                stream.write(data)
            stream.write(b"\n")
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    with (root / "repeat-instances.tsv").open("w", newline="") as stream:
        fields = [
            "id",
            "source_id",
            "array_id",
            "repeat_index",
            "orientation",
            "start",
            "end",
            "rna_sequence",
            "mfe_kcal_mol",
            "isolated_pair_support",
            "context_pair_support",
            "support_change",
            "reason",
        ]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in report["instances"]:
            merged = {**row, **row["comparison"], "mfe_kcal_mol": row["isolated"].get("mfe_kcal_mol")}
            merged["id"] = row["id"]
            writer.writerow({k: merged.get(k) for k in fields})
    with (
        (root / "repeat-analysis.fasta").open("w") as fasta,
        (root / "repeat-structures.txt").open("w") as structures,
    ):
        for row in report["instances"]:
            fasta.write(f">{row['id']}\n{row['rna_sequence']}\n")
            if row["isolated"]["status"] == "completed":
                structures.write(f">{row['id']}\n{row['rna_sequence']}\n{row['isolated']['structure']}\n")
    os.replace(temporary, root / "repeat-analysis.json")
