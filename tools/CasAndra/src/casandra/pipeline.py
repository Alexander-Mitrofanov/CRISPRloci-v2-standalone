"""Canonical end-to-end offline Cas-only genome prediction pipeline."""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from casandra import __version__
from casandra.architecture_runtime import ArchitectureModel, load_architecture_model
from casandra.cassette import construct_cassettes
from casandra.cassette_filter import passes_genome_evidence_gate
from casandra.cassette_hybrid import classify_hybrid_cassette
from casandra.model import ProteinModel, validate_threads, verify_bundle
from casandra.orf_fast import extract_genome_proteins_fast
from casandra.prediction import predict_proteins
from casandra.reporting import read_jsonl, write_jsonl, write_outputs
from casandra.supplied import extract_genbank_proteins, extract_gff3_proteins
from casandra.utils import sha256_file

# Historical research code imports these private names. Keep aliases while the
# installed production entry points use the public reporting module.
_jsonl = read_jsonl
_write_jsonl = write_jsonl
_write_outputs = write_outputs


def _input_record(kind: str, path: Path) -> dict[str, object]:
    return {"kind": kind, "name": Path(path).name, "sha256": sha256_file(Path(path))}


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
    protein_model: ProteinModel | None = None,
    architecture_model: ArchitectureModel | None = None,
    verified_bundle: tuple[dict, dict] | None = None,
) -> dict[str, object]:
    """Analyze one genome into structured files suitable for downstream visualization."""
    output_dir = Path(output_dir)
    model_bundle = Path(model_bundle)
    threads = validate_threads(threads)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace run output: {output_dir}")
    if genbank is not None and (genome_fasta is not None or gff3 is not None):
        raise ValueError("GenBank input is mutually exclusive with FASTA/GFF3 input")
    if genbank is None and genome_fasta is None:
        raise ValueError("A genome FASTA or GenBank input is required")
    if gff3 is not None and genome_fasta is None:
        raise ValueError("GFF3 requires its genome FASTA")

    manifest, config = verified_bundle or verify_bundle(model_bundle)
    protein_model = protein_model or ProteinModel.from_directory(
        model_bundle / config["protein_model_directory"], verify=False
    )
    architecture_model = architecture_model or load_architecture_model(
        model_bundle / config["cassette_architecture_directory"]
    )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    started = time.perf_counter()
    try:
        extraction = staging / "extraction"
        if genbank is not None:
            extract_genbank_proteins(Path(genbank), extraction)
            input_mode = "genbank_supplied_cds"
            inputs = [_input_record("genbank", Path(genbank))]
        elif gff3 is not None:
            extract_gff3_proteins(
                Path(genome_fasta), Path(gff3), extraction, translation_table
            )
            input_mode = "gff3_supplied_cds"
            inputs = [
                _input_record("genome_fasta", Path(genome_fasta)),
                _input_record("gff3", Path(gff3)),
            ]
        else:
            extract_genome_proteins_fast(
                Path(genome_fasta),
                extraction,
                mode=gene_mode,
                translation_table=translation_table,
            )
            input_mode = "pyrodigal_de_novo"
            inputs = [_input_record("genome_fasta", Path(genome_fasta))]

        predict_proteins(
            protein_model,
            extraction / "proteins.faa",
            staging / "protein_predictions.jsonl",
            threads=threads,
        )
        genes = list(read_jsonl(extraction / "genes.jsonl"))
        predictions = list(read_jsonl(staging / "protein_predictions.jsonl"))
        candidates = construct_cassettes(genes, predictions, **config["cassette_parameters"])
        genes_by_id = {row["protein_id"]: row for row in genes}
        predictions_by_id = {row["sequence_id"]: row for row in predictions}
        accepted: list[dict] = []
        rejected: list[dict] = []
        gate_policy = config.get("genome_evidence_gate", {})
        for cassette in candidates:
            cassette["final_classification"] = classify_hybrid_cassette(
                cassette,
                architecture_model,
                genes_by_id,
                predictions_by_id,
            )
            keep, reason = passes_genome_evidence_gate(
                cassette, genes_by_id, predictions_by_id, gate_policy
            )
            cassette["genome_evidence_gate"] = {"accepted": keep, "rule": reason}
            (accepted if keep else rejected).append(cassette)

        write_outputs(staging, genes, predictions, accepted)
        write_jsonl(staging / "rejected_cassette_candidates.jsonl", rejected)
        run = {
            "schema_version": 5,
            "program": "CasAndra",
            "program_version": __version__,
            "bundle_id": config["bundle_id"],
            "bundle_role": config.get("bundle_role", manifest.get("role")),
            "created_at": datetime.now(timezone.utc).isoformat(),
            "inputs": inputs,
            "input_mode": input_mode,
            "cpu_threads": threads,
            "wall_seconds": time.perf_counter() - started,
            "genes": len(genes),
            "cas_proteins": sum(bool(row["is_cas"]) for row in predictions),
            "raw_cassette_candidates": len(candidates),
            "rejected_cassette_candidates": len(rejected),
            "cas_cassettes": len(accepted),
            "classified_cassettes": len(accepted),
            "crispr_array_prediction": False,
            "offline_inference": True,
            "genome_evidence_gate": gate_policy,
        }
        (staging / "run.json").write_text(
            json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        output_manifest = {
            "schema_version": 5,
            "bundle_manifest_sha256": sha256_file(model_bundle / "manifest.json"),
            "files": {},
        }
        for path in sorted(staging.rglob("*")):
            if path.is_file():
                output_manifest["files"][str(path.relative_to(staging))] = {
                    "size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
        (staging / "manifest.json").write_text(
            json.dumps(output_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        staging.replace(output_dir)
        return run
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
