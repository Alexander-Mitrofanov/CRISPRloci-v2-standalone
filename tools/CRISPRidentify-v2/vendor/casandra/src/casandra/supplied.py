"""Reuse supplied GenBank or GFF3 CDS annotations with explicit provenance."""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qsl, quote

from Bio import SeqIO
from Bio.Seq import Seq

from casandra.orf import read_nucleotide_fasta, sha256_file


def _wrap(sequence: str, width: int = 80):
    for offset in range(0, len(sequence), width):
        yield sequence[offset : offset + width]


def _identifier(values: dict[str, list[str]], fallback: str) -> str:
    for key in ("protein_id", "locus_tag", "gene", "ID", "Parent"):
        if values.get(key):
            return str(values[key][0])
    return fallback


def _write_annotation_dataset(records, source: Path, source_format: str, output_dir: Path):
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace annotation extraction: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        with (
            (staging / "genes.jsonl").open("w", encoding="utf-8") as genes,
            (staging / "proteins.faa").open("w", encoding="utf-8") as fasta,
            (staging / "genes.gff3").open("w", encoding="utf-8") as gff,
        ):
            gff.write("##gff-version 3\n")
            for row in records:
                json.dump(row, genes, sort_keys=True, separators=(",", ":"))
                genes.write("\n")
                fasta.write(
                    f">{row['protein_id']} origin={quote(str(source), safe='._-/')} "
                    f"contig={quote(row['contig_id'], safe='._-')} start={row['start_1based']} "
                    f"end={row['end_1based_inclusive']} strand={row['strand']} "
                    f"caller=supplied-{source_format}\n"
                )
                for line in _wrap(row["protein_sequence"]):
                    fasta.write(line + "\n")
                gff.write(
                    f"{row['contig_id']}\tCasAndra-supplied\tCDS\t{row['start_1based']}\t"
                    f"{row['end_1based_inclusive']}\t.\t{row['strand']}\t{row['phase']}\t"
                    f"ID={quote(row['protein_id'], safe='._-|=')}\n"
                )
        audit = {
            "schema_version": 1,
            "source_format": source_format,
            "genes": len(records),
            "internal_stop_genes": sum("*" in row["protein_sequence"] for row in records),
            "coordinates_1based_end_inclusive": True,
            "extraction_ready": bool(records) and not any(
                "*" in row["protein_sequence"] for row in records
            ),
        }
        (staging / "audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        manifest = {
            "schema_version": 1,
            "extraction_id": f"casandra-supplied-{source_format}-v1",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "input": {"path": str(source), "sha256": sha256_file(source)},
            "artifacts": {},
        }
        for path in sorted(staging.iterdir()):
            manifest["artifacts"][path.name] = {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        staging.replace(output_dir)
        return manifest
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def extract_genbank_proteins(genbank: Path, output_dir: Path):
    records = []
    source_hash = sha256_file(genbank)
    for contig_order, sequence_record in enumerate(SeqIO.parse(genbank, "genbank")):
        ordinal = 0
        for feature in sequence_record.features:
            if feature.type != "CDS" or feature.qualifiers.get("pseudo") is not None:
                continue
            ordinal += 1
            nucleotide = str(feature.extract(sequence_record.seq)).upper()
            table = int(feature.qualifiers.get("transl_table", [11])[0])
            codon_start = int(feature.qualifiers.get("codon_start", [1])[0])
            nucleotide = nucleotide[codon_start - 1 :]
            supplied = feature.qualifiers.get("translation", [None])[0]
            protein = str(supplied).replace(" ", "").replace("\n", "") if supplied else str(
                Seq(nucleotide).translate(table=table, cds=False)
            )
            if protein.endswith("*"):
                protein = protein[:-1]
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand != -1 else "-"
            fallback = f"{sequence_record.id}|cds={ordinal:06d}|loc={start}-{end}|strand={strand}"
            protein_id = _identifier(feature.qualifiers, fallback)
            records.append(
                {
                    "protein_id": protein_id,
                    "genome_source_path": str(genbank),
                    "genome_source_sha256": source_hash,
                    "contig_id": sequence_record.id,
                    "contig_description": sequence_record.description,
                    "contig_order": contig_order,
                    "gene_order": ordinal,
                    "start_1based": start,
                    "end_1based_inclusive": end,
                    "strand": strand,
                    "phase": codon_start - 1,
                    "partial_5prime": str(feature.location).startswith("[") and "<" in str(feature.location),
                    "partial_3prime": ">" in str(feature.location),
                    "nucleotide_sequence": nucleotide,
                    "nucleotide_sha256": hashlib.sha256(nucleotide.encode()).hexdigest(),
                    "protein_sequence": protein,
                    "protein_length": len(protein),
                    "protein_sha256": hashlib.sha256(protein.encode()).hexdigest(),
                    "translation_table": table,
                    "caller": {"program": "supplied_annotation", "format": "genbank"},
                }
            )
    return _write_annotation_dataset(records, genbank, "genbank", output_dir)


def _gff_attributes(raw: str) -> dict[str, list[str]]:
    values: dict[str, list[str]] = defaultdict(list)
    for key, value in parse_qsl(raw.replace(";", "&"), keep_blank_values=True):
        values[key].append(value)
    return dict(values)


def extract_gff3_proteins(genome_fasta: Path, gff3: Path, output_dir: Path, translation_table=11):
    contigs = {item.contig_id: item for item in read_nucleotide_fasta(genome_fasta)}
    features = []
    with gff3.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Expected nine GFF3 columns at line {line_number}")
            if fields[2] != "CDS":
                continue
            contig_id, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            if contig_id not in contigs or strand not in {"+", "-"}:
                raise ValueError(f"Invalid GFF3 CDS reference at line {line_number}")
            phase = 0 if fields[7] == "." else int(fields[7])
            attrs = _gff_attributes(fields[8])
            protein_id = _identifier(attrs, f"{contig_id}|cds-line={line_number}")
            sequence = contigs[contig_id].sequence[start - 1 : end]
            if strand == "-":
                sequence = str(Seq(sequence).reverse_complement())
            sequence = sequence[phase:]
            protein = str(Seq(sequence).translate(table=translation_table, cds=False))
            if protein.endswith("*"):
                protein = protein[:-1]
            features.append(
                {
                    "protein_id": protein_id,
                    "genome_source_path": str(genome_fasta),
                    "genome_source_sha256": sha256_file(genome_fasta),
                    "annotation_source_path": str(gff3),
                    "annotation_source_sha256": sha256_file(gff3),
                    "contig_id": contig_id,
                    "contig_description": contigs[contig_id].description,
                    "contig_order": contigs[contig_id].order,
                    "gene_order": 0,
                    "start_1based": start,
                    "end_1based_inclusive": end,
                    "strand": strand,
                    "phase": phase,
                    "partial_5prime": False,
                    "partial_3prime": False,
                    "nucleotide_sequence": sequence,
                    "nucleotide_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                    "protein_sequence": protein,
                    "protein_length": len(protein),
                    "protein_sha256": hashlib.sha256(protein.encode()).hexdigest(),
                    "translation_table": translation_table,
                    "caller": {"program": "supplied_annotation", "format": "gff3"},
                }
            )
    for ordinal, row in enumerate(sorted(features, key=lambda x: (x["contig_order"], x["start_1based"], x["end_1based_inclusive"])), start=1):
        row["gene_order"] = ordinal
    return _write_annotation_dataset(features, gff3, "gff3", output_dir)
