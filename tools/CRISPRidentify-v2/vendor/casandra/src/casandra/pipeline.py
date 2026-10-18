"""End-to-end offline Cas-only genome prediction pipeline."""

from __future__ import annotations

import csv
import html
import json
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

from casandra.cassette import construct_cassettes
from casandra.cassette_hybrid_v2 import classify_hybrid_cassette_v2
from casandra.orf_v2 import extract_genome_proteins_v2
from casandra.production_v3 import predict_proteins_v3
from casandra.supplied import extract_genbank_proteins, extract_gff3_proteins
from casandra.utils import sha256_file


def _jsonl(path):
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _write_jsonl(path, rows):
    with Path(path).open("w", encoding="utf-8") as handle:
        for row in rows:
            json.dump(row, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")


def _write_outputs(staging, genes, predictions, cassettes):
    prediction_by_id = {row["sequence_id"]: row for row in predictions}
    protein_rows = [{**gene, "prediction": prediction_by_id[gene["protein_id"]]} for gene in genes]
    _write_jsonl(staging / "proteins.jsonl", protein_rows)
    _write_jsonl(staging / "cassettes.jsonl", cassettes)
    with (staging / "cas_proteins.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_id", "contig", "start", "end", "strand", "type", "subtype", "margin", "profile"])
        for gene in genes:
            prediction = prediction_by_id[gene["protein_id"]]
            if not prediction["is_cas"]:
                continue
            writer.writerow([
                gene["protein_id"], gene["contig_id"], gene["start_1based"],
                gene["end_1based_inclusive"], gene["strand"],
                prediction["classification"]["type"], prediction["classification"]["subtype"],
                f"{prediction['score_margin']:.6g}", prediction["best_positive_profile"],
            ])
    with (staging / "cassettes.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["cassette_id", "contig", "start", "end", "cas_genes", "class", "type", "subtype", "method", "confidence"])
        for cassette in cassettes:
            final = cassette["final_classification"]
            writer.writerow([
                cassette["cassette_id"], cassette["contig_id"], cassette["start_1based"],
                cassette["end_1based_inclusive"], cassette["cas_gene_count"], final["class"],
                final["type"], final["subtype"], final["method"], f"{final['confidence']:.6g}",
            ])
    with (staging / "casandra.gff3").open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for cassette in cassettes:
            final = cassette["final_classification"]
            attrs = (
                f"ID={quote(cassette['cassette_id'], safe='._-|=')};"
                f"cas_type={final['type'] or ''};cas_subtype={final['subtype'] or ''};"
                f"method={final['method']};confidence={final['confidence']:.6g}"
            )
            handle.write(
                f"{cassette['contig_id']}\tCasAndra\tCas_cassette\t{cassette['start_1based']}\t"
                f"{cassette['end_1based_inclusive']}\t.\t.\t.\t{attrs}\n"
            )
        for gene in genes:
            prediction = prediction_by_id[gene["protein_id"]]
            if prediction["is_cas"]:
                attrs = (
                    f"ID={quote(gene['protein_id'], safe='._-|=')};"
                    f"cas_type={prediction['classification']['type'] or ''};"
                    f"cas_subtype={prediction['classification']['subtype'] or ''};"
                    f"profile={prediction['best_positive_profile'] or ''};"
                    f"margin={prediction['score_margin']:.6g}"
                )
                handle.write(
                    f"{gene['contig_id']}\tCasAndra\tCas_protein\t{gene['start_1based']}\t"
                    f"{gene['end_1based_inclusive']}\t{prediction['positive_profile_score']:.6g}\t"
                    f"{gene['strand']}\t0\t{attrs}\n"
                )
    rows = "".join(
        "<tr>" + "".join(
            f"<td>{html.escape(str(value if value is not None else ''))}</td>"
            for value in (
                cassette["contig_id"], cassette["start_1based"], cassette["end_1based_inclusive"],
                cassette["cas_gene_count"], cassette["final_classification"]["type"],
                cassette["final_classification"]["subtype"], cassette["final_classification"]["method"],
            )
        ) + "</tr>"
        for cassette in cassettes
    )
    (staging / "report.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>CasAndra report</title>"
        "<style>body{font-family:system-ui;margin:2rem}table{border-collapse:collapse}"
        "th,td{border:1px solid #bbb;padding:.35rem}</style><h1>CasAndra Cas-only report</h1>"
        "<p>CRISPR arrays are outside this tool's scope.</p><table><thead><tr>"
        "<th>Contig</th><th>Start</th><th>End</th><th>Cas genes</th><th>Type</th>"
        "<th>Subtype</th><th>Method</th></tr></thead><tbody>" + rows + "</tbody></table>",
        encoding="utf-8",
    )


def predict_genome(
    model_bundle: Path,
    output_dir: Path,
    *,
    genome_fasta: Path | None = None,
    genbank: Path | None = None,
    gff3: Path | None = None,
    gene_mode: str = "auto",
    translation_table: int = 11,
    threads: int = 6,
):
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace run output: {output_dir}")
    if genbank is not None and (genome_fasta is not None or gff3 is not None):
        raise ValueError("GenBank input is mutually exclusive with FASTA/GFF3 input")
    if genbank is None and genome_fasta is None:
        raise ValueError("A genome FASTA or GenBank input is required")
    if gff3 is not None and genome_fasta is None:
        raise ValueError("GFF3 requires its genome FASTA")
    config = json.loads((model_bundle / "config.json").read_text(encoding="utf-8"))
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    started = time.perf_counter()
    try:
        extraction = staging / "extraction"
        if genbank is not None:
            extract_genbank_proteins(genbank, extraction)
            input_mode = "genbank_supplied_cds"
            source = genbank
        elif gff3 is not None:
            extract_gff3_proteins(genome_fasta, gff3, extraction, translation_table)
            input_mode = "gff3_supplied_cds"
            source = genome_fasta
        else:
            extract_genome_proteins_v2(
                genome_fasta,
                extraction,
                mode=gene_mode,
                translation_table=translation_table,
            )
            input_mode = "pyrodigal_de_novo"
            source = genome_fasta
        predict_proteins_v3(
            model_bundle / config["protein_model_directory"],
            extraction / "proteins.faa",
            staging / "protein_predictions.jsonl",
            threads=threads,
        )
        genes = list(_jsonl(extraction / "genes.jsonl"))
        predictions = list(_jsonl(staging / "protein_predictions.jsonl"))
        parameters = config["cassette_parameters"]
        cassettes = construct_cassettes(genes, predictions, **parameters)
        genes_by_id = {row["protein_id"]: row for row in genes}
        predictions_by_id = {row["sequence_id"]: row for row in predictions}
        for cassette in cassettes:
            cassette["final_classification"] = classify_hybrid_cassette_v2(
                cassette,
                model_bundle / config["cassette_architecture_directory"],
                genes_by_id,
                predictions_by_id,
            )
        _write_outputs(staging, genes, predictions, cassettes)
        run = {
            "schema_version": 1,
            "program": "CasAndra",
            "bundle_id": config["bundle_id"],
            "created_at": datetime.now(timezone.utc).isoformat(),
            "input": {"path": str(source), "sha256": sha256_file(source), "mode": input_mode},
            "cpu_threads": threads,
            "wall_seconds": time.perf_counter() - started,
            "genes": len(genes),
            "cas_proteins": sum(row["is_cas"] for row in predictions),
            "cas_cassettes": len(cassettes),
            "classified_cassettes": sum(
                row["final_classification"]["subtype"] is not None for row in cassettes
            ),
            "crispr_array_prediction": False,
            "offline_inference": True,
        }
        (staging / "run.json").write_text(
            json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        # Keep extraction provenance, but avoid duplicating its large sequence JSON in top-level files.
        manifest = {
            "schema_version": 1,
            "bundle_manifest_sha256": sha256_file(model_bundle / "manifest.json"),
            "files": {
                str(path.relative_to(staging)): {
                    "size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
                for path in sorted(staging.rglob("*")) if path.is_file()
            },
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        staging.replace(output_dir)
        return run
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
