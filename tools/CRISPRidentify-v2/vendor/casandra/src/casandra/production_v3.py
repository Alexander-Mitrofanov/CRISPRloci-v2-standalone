"""Deployment-safe Type-V-complete production model revision."""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from casandra.decision import is_cas_hit
from casandra.production_v2 import MODEL_ID as V2_MODEL_ID, load_model
from casandra.runtime import parse_fasta, scan_profiles
from casandra.utils import sha256_file


MODEL_ID = "casandra-production-profile-hmm-type-v-v3"


def build_production_model_v3(v2_model_dir: Path, deployment_benchmark: Path, output_dir: Path):
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace production model: {output_dir}")
    output_dir.mkdir(parents=True)
    for name in ("profiles.hmm", "profile_metadata.json", "build_audit.json"):
        shutil.copy2(v2_model_dir / name, output_dir / name)
    prior_card = json.loads((v2_model_dir / "model_card.json").read_text(encoding="utf-8"))
    benchmark = json.loads(deployment_benchmark.read_text(encoding="utf-8"))
    prior_card.update(
        {
            "schema_version": 3,
            "model_id": MODEL_ID,
            "supersedes_model_id": V2_MODEL_ID,
            "decision_rule": "positive_profile_hit_present AND score_margin >= detection_threshold",
            "benchmark_test_detection": benchmark["comparators"]["casandra_reference"]["test_detection"],
        }
    )
    (output_dir / "model_card.json").write_text(
        json.dumps(prior_card, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = {
        "schema_version": 3,
        "model_id": MODEL_ID,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "cpu_only": True,
        "offline_inference": True,
        "inputs": {
            "v2_model_manifest_sha256": sha256_file(v2_model_dir / "manifest.json"),
            "deployment_benchmark_sha256": sha256_file(deployment_benchmark),
        },
        "artifacts": {},
    }
    for path in sorted(output_dir.iterdir()):
        manifest["artifacts"][path.name] = {
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def predict_query_records_v3(model_dir: Path, queries: list[dict[str, str]], threads: int = 6):
    manifest, card, metadata, hmms = load_model(model_dir)
    scans = scan_profiles(queries, hmms, metadata, threads)
    threshold = float(card["detection_threshold"])
    records = []
    for query, scan in zip(queries, scans, strict=True):
        is_cas = is_cas_hit(scan, threshold)
        records.append(
            {
                "sequence_id": query["sequence_id"],
                "sequence_length": len(query["sequence"]),
                "is_cas": is_cas,
                "positive_profile_score": scan["positive_score"],
                "hard_negative_profile_score": scan["negative_score"],
                "score_margin": scan["margin"],
                "best_positive_profile": scan["best_positive_profile"],
                "classification": scan["classification"] if is_cas else {
                    "class": None, "type": None, "subtype": None
                },
                "variant": None,
                "evidence": {
                    "model_id": manifest["model_id"],
                    "profile_hits": scan["hit_count"],
                    "decision_threshold": threshold,
                    "positive_profile_required": True,
                    "type_v_subtype_masked": False,
                },
            }
        )
    return records


def predict_proteins_v3(model_dir: Path, fasta_path: Path, output_path: Path, *, threads=6):
    queries = [{"sequence_id": name, "sequence": seq} for name, seq in parse_fasta(fasta_path)]
    records = predict_query_records_v3(model_dir, queries, threads)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for row in records:
            json.dump(row, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
    return {
        "model_id": MODEL_ID,
        "input_records": len(records),
        "cas_predictions": sum(row["is_cas"] for row in records),
        "output": str(output_path),
        "cpu_threads": threads,
    }
