"""Atomic, provenance-bound annotation of supplied protein FASTA records."""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

from casandra import __version__
from casandra.model import ProteinModel, validate_threads
from casandra.prediction import predict_proteins, profile_cas_family
from casandra.reporting import read_jsonl
from casandra.runtime import parse_fasta
from casandra.utils import sha256_file


def annotate_protein_fasta(
    model_bundle: Path,
    input_fasta: Path,
    output_dir: Path,
    *,
    protein_model: ProteinModel,
    bundle_manifest: dict[str, object],
    bundle_config: dict[str, object],
    threads: int = 6,
) -> dict[str, object]:
    """Annotate every protein and publish an atomic, checksummed result directory."""

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
            staging / "protein_predictions.jsonl",
            threads=threads,
        )
        predictions = [dict(row) for row in read_jsonl(staging / "protein_predictions.jsonl")]
        input_ids = [protein_id for protein_id, _sequence in proteins]
        prediction_ids = [str(row.get("sequence_id") or "") for row in predictions]
        if prediction_ids != input_ids:
            raise RuntimeError("Protein predictions do not preserve every FASTA record in order")
        model_id = prediction_summary.get("model_id")
        if not isinstance(model_id, str) or not model_id:
            raise RuntimeError("Protein model provenance is unavailable")
        for row in predictions:
            evidence = row.get("evidence")
            if not isinstance(evidence, dict) or evidence.get("model_id") != model_id:
                raise RuntimeError("Protein prediction model provenance is inconsistent")
            is_cas = row.get("is_cas")
            if not isinstance(is_cas, bool):
                raise TypeError("Protein prediction contains an invalid Cas call")
            profile = row.get("best_positive_profile")
            family = row.get("cas_family")
            if is_cas and (
                not isinstance(profile, str) or not profile.strip()
            ):
                raise RuntimeError("A positive Cas call has no protein-family profile identity")
            if is_cas and (
                not isinstance(family, str)
                or family != profile_cas_family(protein_model.metadata, profile)
            ):
                raise RuntimeError("A positive Cas call has no curated Cas-family identity")
            if not is_cas and family is not None:
                raise RuntimeError("A negative protein call contains a Cas-family identity")
            expected_result = family if is_cas else "no cas"
            if row.get("result") != expected_result:
                raise RuntimeError(
                    "Protein prediction result disagrees with its Cas-family call"
                )

        cas_count = sum(row.get("is_cas") is True for row in predictions)
        if (
            prediction_summary.get("input_records") != len(proteins)
            or prediction_summary.get("cas_predictions") != cas_count
        ):
            raise RuntimeError("Protein prediction summary counts are inconsistent")
        run = {
            "schema_version": 1,
            "program": "CasAndra",
            "program_version": __version__,
            "analysis": "annotate_cas_genes",
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
            "model_id": model_id,
            "cpu_threads": threads,
            "wall_seconds": time.perf_counter() - started,
            "protein_records": len(proteins),
            "cas_proteins": cas_count,
            "non_cas_proteins": len(proteins) - cas_count,
            "result_contract": {
                "positive": "cas_family",
                "negative": "no cas",
            },
            "crispr_array_prediction": False,
            "offline_inference": True,
        }
        (staging / "run.json").write_text(
            json.dumps(run, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        manifest_files: dict[str, object] = {}
        for path in sorted(staging.iterdir()):
            if path.is_file():
                manifest_files[path.name] = {
                    "size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
        output_manifest = {
            "schema_version": 1,
            "bundle_manifest_sha256": bundle_manifest_sha256,
            "files": manifest_files,
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
