"""Linear-time Pyrodigal extraction with one input digest per run."""

from __future__ import annotations

import hashlib
import json
import shutil
import tempfile
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import pyrodigal
from Bio.Seq import Seq

from casandra.orf import _protein_id, _wrap, choose_gene_mode, read_nucleotide_fasta, sha256_file


def _write_json(path: Path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def extract_genome_proteins_fast(genome_fasta: Path, output_dir: Path, *, mode="auto", translation_table=11, closed=False, mask=True):
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace extraction: {output_dir}")
    contigs = read_nucleotide_fasta(genome_fasta)
    source_hash = sha256_file(genome_fasta)
    selected_mode = choose_gene_mode(contigs, mode)
    finder = pyrodigal.GeneFinder(meta=selected_mode == "meta", closed=closed, mask=mask)
    if selected_mode == "single":
        finder.train(contigs[0].sequence, *(x.sequence for x in contigs[1:]), translation_table=translation_table)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    total = partial = internal = 0
    table_counts = {}
    try:
        with (
            (staging / "contigs.jsonl").open("w", encoding="utf-8") as ch,
            (staging / "genes.jsonl").open("w", encoding="utf-8") as gh,
            (staging / "proteins.faa").open("w", encoding="utf-8") as ph,
            (staging / "genes.gff3").open("w", encoding="utf-8") as gf,
        ):
            gf.write("##gff-version 3\n")
            for contig in contigs:
                contig_record = {**asdict(contig), "sequence": contig.sequence, "sequence_length": len(contig.sequence), "sequence_sha256": hashlib.sha256(contig.sequence.encode()).hexdigest()}
                ch.write(json.dumps(contig_record, sort_keys=True, separators=(",", ":")) + "\n")
                for ordinal, gene in enumerate(finder.find_genes(contig.sequence), start=1):
                    nucleotide = str(gene.sequence())
                    table = int(gene.translation_table)
                    protein = str(Seq(nucleotide).translate(table=table, cds=False))
                    if protein.endswith("*"):
                        protein = protein[:-1]
                    protein_id = _protein_id(contig.contig_id, ordinal, gene.begin, gene.end, gene.strand)
                    is_partial = bool(gene.partial_begin or gene.partial_end)
                    has_internal = "*" in protein
                    total += 1
                    partial += is_partial
                    internal += has_internal
                    table_counts[str(table)] = table_counts.get(str(table), 0) + 1
                    record = {
                        "protein_id": protein_id, "genome_source_path": str(genome_fasta),
                        "genome_source_sha256": source_hash, "contig_id": contig.contig_id,
                        "contig_description": contig.description, "contig_order": contig.order,
                        "gene_order": ordinal, "start_1based": int(gene.begin),
                        "end_1based_inclusive": int(gene.end), "strand": "+" if gene.strand == 1 else "-",
                        "phase": 0,
                        "partial_5prime": bool(gene.partial_begin if gene.strand == 1 else gene.partial_end),
                        "partial_3prime": bool(gene.partial_end if gene.strand == 1 else gene.partial_begin),
                        "nucleotide_sequence": nucleotide, "nucleotide_sha256": hashlib.sha256(nucleotide.encode()).hexdigest(),
                        "protein_sequence": protein, "protein_length": len(protein),
                        "protein_sha256": hashlib.sha256(protein.encode()).hexdigest(),
                        "translation_table": table, "translation_policy": "caller_selected_per_gene",
                        "start_type": gene.start_type, "rbs_motif": gene.rbs_motif,
                        "rbs_spacer": gene.rbs_spacer, "gene_score": float(gene.score),
                        "gene_confidence": float(gene.confidence()), "gc_content": float(gene.gc_cont),
                        "caller": {"program": "Pyrodigal", "version": pyrodigal.__version__, "mode": selected_mode, "requested_mode": mode, "closed": closed, "mask": mask, "translation_table": table},
                    }
                    gh.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
                    ph.write(f">{protein_id} origin={quote(str(genome_fasta), safe='._-/')} contig={quote(contig.contig_id, safe='._-')} start={gene.begin} end={gene.end} strand={record['strand']} table={table} caller=pyrodigal-{pyrodigal.__version__}\n")
                    for line in _wrap(protein):
                        ph.write(line + "\n")
                    attrs = f"ID={quote(protein_id, safe='._-|=')};translation_table={table};partial5={int(record['partial_5prime'])};partial3={int(record['partial_3prime'])}"
                    gf.write(f"{contig.contig_id}\tPyrodigal\tCDS\t{gene.begin}\t{gene.end}\t{gene.score:.6g}\t{record['strand']}\t0\t{attrs}\n")
        audit = {
            "schema_version": 3, "input_contigs": len(contigs),
            "input_nucleotides": sum(len(x.sequence) for x in contigs), "selected_mode": selected_mode,
            "genes": total, "partial_genes": partial, "internal_stop_genes": internal,
            "translation_policy": "caller_selected_per_gene", "translation_table_counts": table_counts,
            "coordinates_1based_end_inclusive": True, "protein_headers_include_origin_and_coordinates": True,
            "caller_cpu_only": True, "extraction_ready": total > 0 and internal == 0,
            "performance_fix": "genome SHA-256 computed once per run rather than once per gene",
        }
        _write_json(staging / "audit.json", audit)
        manifest = {
            "schema_version": 3, "extraction_id": "casandra-pyrodigal-proteins-v3-linear-time",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "input": {"path": str(genome_fasta), "sha256": source_hash},
            "caller": {"program": "Pyrodigal", "version": pyrodigal.__version__, "upstream_algorithm": "Prodigal v2.6.3+31b300a with Pyrodigal fixes", "mode": selected_mode, "requested_mode": mode, "translation_table": translation_table, "translation_policy": "caller_selected_per_gene", "closed": closed, "mask": mask, "cpu_only": True},
            "artifacts": {},
        }
        for path in sorted(staging.iterdir()):
            manifest["artifacts"][path.name] = {"size": path.stat().st_size, "sha256": sha256_file(path)}
        _write_json(staging / "manifest.json", manifest)
        staging.replace(output_dir)
        return manifest
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
