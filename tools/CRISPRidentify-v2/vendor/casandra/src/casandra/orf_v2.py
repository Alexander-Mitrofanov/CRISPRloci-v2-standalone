"""Genetic-code-safe revision of Pyrodigal protein extraction."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from Bio.Seq import Seq

from casandra.orf import extract_genome_proteins, sha256_file


def _wrap(sequence: str, width: int = 80):
    for offset in range(0, len(sequence), width):
        yield sequence[offset : offset + width]


def extract_genome_proteins_v2(
    genome_fasta: Path,
    output_dir: Path,
    *,
    mode: str = "auto",
    translation_table: int = 11,
    closed: bool = False,
    mask: bool = True,
):
    """Call genes, then translate each CDS with its caller-selected code."""
    manifest = extract_genome_proteins(
        genome_fasta,
        output_dir,
        mode=mode,
        translation_table=translation_table,
        closed=closed,
        mask=mask,
    )
    rows = []
    with (output_dir / "genes.jsonl").open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            table = int(row["translation_table"])
            protein = str(Seq(row["nucleotide_sequence"]).translate(table=table, cds=False))
            if protein.endswith("*"):
                protein = protein[:-1]
            row["protein_sequence"] = protein
            row["protein_length"] = len(protein)
            row["protein_sha256"] = hashlib.sha256(protein.encode("ascii")).hexdigest()
            row["translation_policy"] = "caller_selected_per_gene"
            rows.append(row)
    with (output_dir / "genes.jsonl").open("w", encoding="utf-8") as handle:
        for row in rows:
            json.dump(row, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
    with (output_dir / "proteins.faa").open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                f">{row['protein_id']} origin={row['genome_source_path']} contig={row['contig_id']} "
                f"start={row['start_1based']} end={row['end_1based_inclusive']} strand={row['strand']} "
                f"table={row['translation_table']} caller=pyrodigal-{row['caller']['version']}\n"
            )
            for line in _wrap(row["protein_sequence"]):
                handle.write(line + "\n")
    audit_path = output_dir / "audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["schema_version"] = 2
    audit["translation_policy"] = "caller_selected_per_gene"
    audit["translation_table_counts"] = {
        str(table): sum(int(row["translation_table"]) == table for row in rows)
        for table in sorted({int(row["translation_table"]) for row in rows})
    }
    audit["internal_stop_genes"] = sum("*" in row["protein_sequence"] for row in rows)
    audit["extraction_ready"] = bool(rows) and audit["internal_stop_genes"] == 0
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = output_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = 2
    manifest["extraction_id"] = "casandra-pyrodigal-proteins-v2"
    manifest["caller"]["translation_policy"] = "caller_selected_per_gene"
    manifest["supersedes"] = "v1 fixed-request translation in metagenomic mode"
    for name in ("genes.jsonl", "proteins.faa", "audit.json"):
        path = output_dir / name
        manifest["artifacts"][name] = {"size": path.stat().st_size, "sha256": sha256_file(path)}
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
