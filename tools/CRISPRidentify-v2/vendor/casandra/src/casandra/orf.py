"""Reproducible CPU-only prokaryotic protein extraction with Pyrodigal."""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, TextIO
from urllib.parse import quote

import pyrodigal


ORF_SCHEMA_VERSION = 1
DNA_ALPHABET = set("ACGTURYSWKMBDHVN")


@dataclass(frozen=True)
class Contig:
    contig_id: str
    description: str
    sequence: str
    order: int


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_nucleotide_fasta(path: Path) -> list[Contig]:
    """Read a nucleotide FASTA, preserving full descriptions and input order."""
    contigs: list[Contig] = []
    header: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        nonlocal header, chunks
        if header is None:
            return
        contig_id = header.split(maxsplit=1)[0]
        sequence = "".join(chunks).upper().replace("U", "T")
        invalid = sorted(set(sequence) - DNA_ALPHABET)
        if invalid:
            raise ValueError(f"Invalid nucleotide symbols in {contig_id}: {''.join(invalid)}")
        if not sequence:
            raise ValueError(f"Empty contig: {contig_id}")
        contigs.append(Contig(contig_id, header, sequence, len(contigs)))

    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                chunks = []
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
            elif header is None:
                raise ValueError(f"Nucleotide sequence before header at {path}:{line_number}")
            else:
                chunks.append(line)
    finish()
    if not contigs:
        raise ValueError(f"No nucleotide records in {path}")
    ids = [contig.contig_id for contig in contigs]
    if len(ids) != len(set(ids)):
        raise ValueError("Nucleotide FASTA contains duplicate primary contig identifiers")
    return contigs


def choose_gene_mode(contigs: list[Contig], requested: str) -> str:
    """Resolve auto mode from the documented Pyrodigal training-size boundary."""
    if requested not in {"auto", "single", "meta"}:
        raise ValueError("Gene mode must be auto, single, or meta")
    if requested != "auto":
        return requested
    total = sum(len(contig.sequence) for contig in contigs)
    return "single" if total >= pyrodigal.IDEAL_SINGLE_GENOME else "meta"


def _write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _write_record(handle: TextIO, value: object) -> None:
    json.dump(value, handle, sort_keys=True, separators=(",", ":"))
    handle.write("\n")


def _protein_id(contig_id: str, ordinal: int, begin: int, end: int, strand: int) -> str:
    encoded = quote(contig_id, safe="._-")
    sign = "+" if strand == 1 else "-"
    return f"casandra|contig={encoded}|cds={ordinal:06d}|loc={begin}-{end}|strand={sign}"


def _wrap(sequence: str, width: int = 80) -> Iterator[str]:
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def extract_genome_proteins(
    genome_fasta: Path,
    output_dir: Path,
    *,
    mode: str = "auto",
    translation_table: int = 11,
    closed: bool = False,
    mask: bool = True,
) -> dict[str, object]:
    """Call prokaryotic CDSs and write explicit proteins, coordinates, and provenance."""
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace protein extraction: {output_dir}")
    contigs = read_nucleotide_fasta(genome_fasta)
    selected_mode = choose_gene_mode(contigs, mode)
    finder = pyrodigal.GeneFinder(meta=selected_mode == "meta", closed=closed, mask=mask)
    if selected_mode == "single":
        finder.train(
            contigs[0].sequence,
            *(contig.sequence for contig in contigs[1:]),
            translation_table=translation_table,
        )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    total_genes = 0
    partial_genes = 0
    internal_stop_genes = 0
    try:
        with (
            (staging / "contigs.jsonl").open("w", encoding="utf-8") as contig_handle,
            (staging / "genes.jsonl").open("w", encoding="utf-8") as gene_handle,
            (staging / "proteins.faa").open("w", encoding="utf-8") as protein_handle,
            (staging / "genes.gff3").open("w", encoding="utf-8") as gff_handle,
        ):
            gff_handle.write("##gff-version 3\n")
            for contig in contigs:
                _write_record(
                    contig_handle,
                    {
                        **asdict(contig),
                        "sequence": contig.sequence,
                        "sequence_length": len(contig.sequence),
                        "sequence_sha256": hashlib.sha256(
                            contig.sequence.encode("ascii")
                        ).hexdigest(),
                    },
                )
                genes = finder.find_genes(contig.sequence)
                for ordinal, gene in enumerate(genes, start=1):
                    protein = str(
                        gene.translate(
                            translation_table=translation_table,
                            include_stop=False,
                            strict=True,
                        )
                    )
                    nucleotide = str(gene.sequence())
                    protein_id = _protein_id(
                        contig.contig_id, ordinal, gene.begin, gene.end, gene.strand
                    )
                    partial = bool(gene.partial_begin or gene.partial_end)
                    internal_stop = "*" in protein
                    total_genes += 1
                    partial_genes += partial
                    internal_stop_genes += internal_stop
                    record = {
                        "protein_id": protein_id,
                        "genome_source_path": str(genome_fasta),
                        "genome_source_sha256": sha256_file(genome_fasta),
                        "contig_id": contig.contig_id,
                        "contig_description": contig.description,
                        "contig_order": contig.order,
                        "gene_order": ordinal,
                        "start_1based": int(gene.begin),
                        "end_1based_inclusive": int(gene.end),
                        "strand": "+" if gene.strand == 1 else "-",
                        "phase": 0,
                        "partial_5prime": bool(
                            gene.partial_begin if gene.strand == 1 else gene.partial_end
                        ),
                        "partial_3prime": bool(
                            gene.partial_end if gene.strand == 1 else gene.partial_begin
                        ),
                        "nucleotide_sequence": nucleotide,
                        "nucleotide_sha256": hashlib.sha256(
                            nucleotide.encode("ascii")
                        ).hexdigest(),
                        "protein_sequence": protein,
                        "protein_length": len(protein),
                        "protein_sha256": hashlib.sha256(protein.encode("ascii")).hexdigest(),
                        "translation_table": int(gene.translation_table),
                        "start_type": gene.start_type,
                        "rbs_motif": gene.rbs_motif,
                        "rbs_spacer": gene.rbs_spacer,
                        "gene_score": float(gene.score),
                        "gene_confidence": float(gene.confidence()),
                        "gc_content": float(gene.gc_cont),
                        "caller": {
                            "program": "Pyrodigal",
                            "version": pyrodigal.__version__,
                            "mode": selected_mode,
                            "requested_mode": mode,
                            "closed": closed,
                            "mask": mask,
                            "translation_table": translation_table,
                        },
                    }
                    _write_record(gene_handle, record)
                    header = (
                        f">{protein_id} origin={quote(str(genome_fasta), safe='._-/')} "
                        f"contig={quote(contig.contig_id, safe='._-')} "
                        f"start={gene.begin} end={gene.end} strand={record['strand']} "
                        f"partial5={int(record['partial_5prime'])} "
                        f"partial3={int(record['partial_3prime'])} caller=pyrodigal-{pyrodigal.__version__}"
                    )
                    protein_handle.write(header + "\n")
                    for line in _wrap(protein):
                        protein_handle.write(line + "\n")
                    attributes = (
                        f"ID={quote(protein_id, safe='._-|=')};"
                        f"translation_table={translation_table};"
                        f"partial5={int(record['partial_5prime'])};"
                        f"partial3={int(record['partial_3prime'])}"
                    )
                    gff_handle.write(
                        f"{contig.contig_id}\tPyrodigal\tCDS\t{gene.begin}\t{gene.end}\t"
                        f"{gene.score:.6g}\t{record['strand']}\t0\t{attributes}\n"
                    )
        audit = {
            "schema_version": ORF_SCHEMA_VERSION,
            "input_contigs": len(contigs),
            "input_nucleotides": sum(len(contig.sequence) for contig in contigs),
            "selected_mode": selected_mode,
            "genes": total_genes,
            "partial_genes": partial_genes,
            "internal_stop_genes": internal_stop_genes,
            "coordinates_1based_end_inclusive": True,
            "protein_headers_include_origin_and_coordinates": True,
            "caller_cpu_only": True,
            "extraction_ready": total_genes > 0 and internal_stop_genes == 0,
        }
        _write_json(staging / "audit.json", audit)
        manifest = {
            "schema_version": ORF_SCHEMA_VERSION,
            "extraction_id": "casandra-pyrodigal-proteins-v1",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "input": {
                "path": str(genome_fasta),
                "sha256": sha256_file(genome_fasta),
            },
            "caller": {
                "program": "Pyrodigal",
                "version": pyrodigal.__version__,
                "upstream_algorithm": "Prodigal v2.6.3+31b300a with Pyrodigal fixes",
                "mode": selected_mode,
                "requested_mode": mode,
                "translation_table": translation_table,
                "closed": closed,
                "mask": mask,
                "cpu_only": True,
            },
            "artifacts": {},
        }
        for path in sorted(staging.iterdir()):
            manifest["artifacts"][path.name] = {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        _write_json(staging / "manifest.json", manifest)
        staging.replace(output_dir)
        return manifest
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
