"""Type-V-complete CPU-only Cas protein model and prediction runtime."""

from __future__ import annotations

import json
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import pyhmmer

from casandra.integrity import verify_artifact_manifest
from casandra.runtime import parse_fasta, scan_profiles
from casandra.utils import sha256_file


MODEL_ID = "casandra-production-profile-hmm-type-v-v2"


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_model(model_dir: Path):
    manifest = verify_artifact_manifest(model_dir)
    card = json.loads((model_dir / "model_card.json").read_text(encoding="utf-8"))
    metadata = json.loads((model_dir / "profile_metadata.json").read_text(encoding="utf-8"))
    hmms = []
    with pyhmmer.plan7.HMMFile(model_dir / "profiles.hmm") as handle:
        while (hmm := handle.read()) is not None:
            hmms.append(hmm)
    return manifest, card, metadata, hmms


def build_production_model_v2(
    dataset_dir: Path,
    profile_catalogue: Path,
    alignments_archive: Path,
    type_v_taxonomy: Path,
    benchmark_metrics: Path,
    comparator_benchmark: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Refit all frozen development folds after the benchmark is locked."""
    from casandra_benchmark.profile_hmm_v2 import build_fold_safe_hmms_v2

    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace production model: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        hmms, metadata, audit = build_fold_safe_hmms_v2(
            dataset_dir,
            profile_catalogue,
            alignments_archive,
            type_v_taxonomy,
            training_splits=frozenset({"train", "validation", "test"}),
        )
        benchmark = json.loads(benchmark_metrics.read_text(encoding="utf-8"))
        comparator = json.loads(comparator_benchmark.read_text(encoding="utf-8"))
        threshold = float(benchmark["test_detection"]["threshold"])
        with (staging / "profiles.hmm").open("wb") as handle:
            for hmm in hmms:
                hmm.write(handle)
        _write_json(staging / "profile_metadata.json", metadata)
        _write_json(staging / "build_audit.json", audit)
        card = {
            "schema_version": 2,
            "model_id": MODEL_ID,
            "purpose": "local CPU-only Cas protein detection and hierarchical classification",
            "detection_threshold": threshold,
            "classification_heads": ["class", "type", "subtype"],
            "operational_type_v_subtypes": [f"V-{letter}" for letter in "ABCDEFGHIJKLMNOPQ"],
            "type_v_subtypes_masked": False,
            "training_splits": ["train", "validation", "test"],
            "training_records_role": "post-evaluation production refit",
            "threshold_source": "frozen validation-selected threshold from benchmark model",
            "benchmark_model_id": benchmark["model_id"],
            "benchmark_test_detection": benchmark["test_detection"],
            "benchmark_test_classification": benchmark["test_classification"],
            "comparator_benchmark_scope": comparator.get("scope", "protein component"),
            "cpu_only": True,
            "offline_inference": True,
            "known_limitations": [
                "no CRISPR-array prediction",
                "variant labels are not emitted without unambiguous supervision",
                "whole-genome SOTA requires an independent end-to-end benchmark",
            ],
        }
        _write_json(staging / "model_card.json", card)
        manifest = {
            "schema_version": 2,
            "model_id": MODEL_ID,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "cpu_only": True,
            "offline_inference": True,
            "inputs": {
                "dataset_manifest_sha256": sha256_file(dataset_dir / "manifest.json"),
                "profile_catalogue_sha256": sha256_file(profile_catalogue),
                "alignments_archive_sha256": sha256_file(alignments_archive),
                "type_v_taxonomy_sha256": sha256_file(type_v_taxonomy),
                "benchmark_metrics_sha256": sha256_file(benchmark_metrics),
                "comparator_benchmark_sha256": sha256_file(comparator_benchmark),
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


def predict_query_records(model_dir: Path, queries: list[dict[str, str]], threads: int = 6):
    manifest, card, metadata, hmms = load_model(model_dir)
    scans = scan_profiles(queries, hmms, metadata, threads)
    threshold = float(card["detection_threshold"])
    records = []
    for query, scan in zip(queries, scans, strict=True):
        is_cas = float(scan["margin"]) >= threshold
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
                    "type_v_subtype_masked": False,
                },
            }
        )
    return records


def predict_proteins_v2(
    model_dir: Path, fasta_path: Path, output_path: Path, *, threads: int = 6
) -> dict[str, object]:
    queries = [
        {"sequence_id": name, "sequence": sequence}
        for name, sequence in parse_fasta(fasta_path)
    ]
    records = predict_query_records(model_dir, queries, threads)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            json.dump(record, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
    return {
        "model_id": MODEL_ID,
        "input_records": len(records),
        "cas_predictions": sum(row["is_cas"] for row in records),
        "type_v_predictions": sum(
            row["classification"]["type"] == "V" for row in records if row["is_cas"]
        ),
        "cpu_threads": threads,
        "output": str(output_path),
    }
