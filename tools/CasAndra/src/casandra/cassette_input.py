"""Classification of an explicitly supplied, ordered Cas-protein cassette."""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from casandra import __version__
from casandra.architecture_runtime import ArchitectureModel
from casandra.cassette import classify_cassette
from casandra.cassette_hybrid import classify_hybrid_cassette
from casandra.model import ProteinModel, validate_threads
from casandra.prediction import predict_proteins
from casandra.reporting import read_jsonl
from casandra.runtime import parse_fasta
from casandra.utils import sha256_file


def _ordered_gene_records(
    proteins: list[tuple[str, str]],
) -> list[dict[str, object]]:
    """Build sequence-free ordering records for the architecture model.

    Supplied protein FASTA has a meaningful record order but no nucleotide
    coordinates.  The architecture model consumes order only, so these
    synthetic positions must never be presented as biological coordinates.
    """

    records: list[dict[str, object]] = []
    offset = 1
    for protein_id, sequence in proteins:
        translated_span = max(1, len(sequence.rstrip("*"))) * 3
        records.append(
            {
                "protein_id": protein_id,
                "start_1based": offset,
                "end_1based_inclusive": offset + translated_span - 1,
                "protein_length": len(sequence.rstrip("*")),
            }
        )
        offset += translated_span + 3
    return records


def classify_ordered_predictions(
    proteins: list[tuple[str, str]],
    predictions: list[dict[str, object]],
    architecture_model: ArchitectureModel,
) -> dict[str, object]:
    """Aggregate per-protein calls into one cassette classification."""

    input_ids = [protein_id for protein_id, _sequence in proteins]
    prediction_ids = [str(row.get("sequence_id") or "") for row in predictions]
    if prediction_ids != input_ids:
        raise ValueError("Protein predictions do not preserve the supplied FASTA record order")

    genes = _ordered_gene_records(proteins)
    genes_by_id = {str(row["protein_id"]): row for row in genes}
    predictions_by_id = {str(row["sequence_id"]): row for row in predictions}
    cas_ids = [
        protein_id
        for protein_id in input_ids
        if bool(predictions_by_id[protein_id].get("is_cas"))
    ]
    direct = classify_cassette(predictions)
    cassette = {
        "cassette_id": "supplied-cassette",
        "cas_protein_ids": cas_ids,
        "classification": direct,
    }
    final = classify_hybrid_cassette(
        cassette,
        architecture_model,
        genes_by_id,
        predictions_by_id,
    )
    return {
        "schema_version": 1,
        "cassette_id": "supplied-cassette",
        "input_mode": "ordered_protein_fasta",
        "input_protein_ids": input_ids,
        "cas_protein_ids": cas_ids,
        "non_cas_protein_ids": [value for value in input_ids if value not in set(cas_ids)],
        "protein_count": len(input_ids),
        "cas_protein_count": len(cas_ids),
        "classification": final,
        "direct_profile_classification": direct,
        "order_used_for_architecture": True,
        "coordinates_available": False,
        "crispr_array_evidence_used": False,
    }


def classify_protein_cassette(
    model_bundle: Path,
    input_fasta: Path,
    output_dir: Path,
    *,
    protein_model: ProteinModel,
    architecture_model: ArchitectureModel,
    bundle_manifest: dict[str, object],
    bundle_config: dict[str, object],
    threads: int = 6,
) -> dict[str, object]:
    """Predict an ordered protein FASTA and publish an atomic result directory."""

    model_bundle = Path(model_bundle)
    input_fasta = Path(input_fasta)
    output_dir = Path(output_dir)
    threads = validate_threads(threads)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace run output: {output_dir}")

    proteins = list(parse_fasta(input_fasta))
    if not proteins:
        raise ValueError("Protein FASTA contains no records")

    bundle_manifest_sha256 = sha256_file(model_bundle / "manifest.json")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    started = time.perf_counter()
    try:
        prediction_summary = predict_proteins(
            protein_model,
            input_fasta,
            staging / "proteins.jsonl",
            threads=threads,
        )
        predictions = [dict(row) for row in read_jsonl(staging / "proteins.jsonl")]
        cassette = classify_ordered_predictions(proteins, predictions, architecture_model)
        (staging / "cassette.json").write_text(
            json.dumps(cassette, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        final = cassette["classification"]
        run = {
            "schema_version": 1,
            "program": "CasAndra",
            "program_version": __version__,
            "analysis": "classify_cassette",
            "bundle_id": bundle_config["bundle_id"],
            "bundle_role": bundle_config.get("bundle_role", bundle_manifest.get("role")),
            "bundle_manifest_sha256": bundle_manifest_sha256,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "inputs": [
                {
                    "kind": "protein_fasta",
                    "name": input_fasta.name,
                    "sha256": sha256_file(input_fasta),
                }
            ],
            "model_id": prediction_summary["model_id"],
            "cpu_threads": threads,
            "wall_seconds": time.perf_counter() - started,
            "protein_records": len(proteins),
            "cas_proteins": cassette["cas_protein_count"],
            "classification": {
                "class": final["class"],
                "type": final["type"],
                "subtype": final["subtype"],
                "method": final["method"],
                "confidence": final["confidence"],
            },
            "crispr_array_prediction": False,
            "offline_inference": True,
        }
        (staging / "run.json").write_text(
            json.dumps(run, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        output_manifest: dict[str, object] = {
            "schema_version": 1,
            "bundle_manifest_sha256": bundle_manifest_sha256,
            "files": {},
        }
        manifest_files = output_manifest["files"]
        assert isinstance(manifest_files, dict)
        for path in sorted(staging.iterdir()):
            if path.is_file():
                manifest_files[path.name] = {
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
