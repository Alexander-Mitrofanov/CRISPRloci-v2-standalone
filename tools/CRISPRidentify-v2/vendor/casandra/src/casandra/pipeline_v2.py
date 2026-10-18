"""Genome pipeline with the validation-locked high-precision evidence gate."""

from __future__ import annotations

import json
from pathlib import Path

from casandra.cassette_filter import passes_genome_evidence_gate
from casandra.pipeline import _jsonl, _write_jsonl, _write_outputs, predict_genome
from casandra.utils import sha256_file


def predict_genome_v2(model_bundle: Path, output_dir: Path, **kwargs):
    run = predict_genome(model_bundle, output_dir, **kwargs)
    genes = list(_jsonl(output_dir / "extraction" / "genes.jsonl"))
    predictions = list(_jsonl(output_dir / "protein_predictions.jsonl"))
    candidates = list(_jsonl(output_dir / "cassettes.jsonl"))
    genes_by_id = {row["protein_id"]: row for row in genes}
    predictions_by_id = {row["sequence_id"]: row for row in predictions}
    accepted = []
    rejected = []
    for cassette in candidates:
        keep, reason = passes_genome_evidence_gate(cassette, genes_by_id, predictions_by_id)
        cassette["genome_evidence_gate"] = {"accepted": keep, "rule": reason}
        (accepted if keep else rejected).append(cassette)
    _write_outputs(output_dir, genes, predictions, accepted)
    _write_jsonl(output_dir / "rejected_cassette_candidates.jsonl", rejected)
    run.update(
        {
            "pipeline_schema_version": 2,
            "raw_cassette_candidates": len(candidates),
            "rejected_cassette_candidates": len(rejected),
            "cas_cassettes": len(accepted),
            "classified_cassettes": len(accepted),
            "genome_evidence_gate": (
                "supported subtype; Class 1 >=2 Cas genes; Class 2 >=2 genes or "
                "singleton length>=400 aa and positive profile score>=50"
            ),
        }
    )
    (output_dir / "run.json").write_text(
        json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest_path = output_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = 2
    manifest["deployment_gate_applied"] = True
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path != manifest_path:
            manifest["files"][str(path.relative_to(output_dir))] = {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return run
