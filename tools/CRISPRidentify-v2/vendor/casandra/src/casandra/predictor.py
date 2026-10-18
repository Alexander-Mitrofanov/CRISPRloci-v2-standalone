"""Production CasAndra CPU-only Cas protein predictor."""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import pyhmmer

from casandra.runtime import parse_fasta, scan_profiles
from casandra.utils import sha256_file


def _write_json(path: Path, value: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def build_production_model(
    dataset_dir: Path,
    profile_catalogue: Path,
    alignments_archive: Path,
    benchmark_metrics: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Refit profiles on all frozen development folds after evaluation is complete."""
    from casandra_benchmark.profile_hmm import build_fold_safe_hmms

    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace existing CasAndra model: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        hmms, metadata, audit = build_fold_safe_hmms(
            dataset_dir,
            profile_catalogue,
            alignments_archive,
            training_splits=frozenset({"train", "validation", "test"}),
        )
        with (temporary / "profiles.hmm").open("wb") as handle:
            for hmm in hmms:
                hmm.write(handle)
        _write_json(temporary / "profile_metadata.json", metadata)
        _write_json(temporary / "build_audit.json", audit)
        benchmark = json.loads(benchmark_metrics.read_text(encoding="utf-8"))
        model_card = {
            "schema_version": 1,
            "model_id": "casandra-production-profile-hmm-v1",
            "detection_threshold": 0.0,
            "classification": ["class", "type", "subtype"],
            "variant_output": "unavailable_no_unambiguous_training_labels",
            "type_v_subtype_output": "masked_pending_authoritative_reconciliation",
            "training_splits": ["train", "validation", "test"],
            "open_set_training": False,
            "benchmark_model_id": benchmark["model_id"],
            "benchmark_test_detection": benchmark["test_detection"],
            "cpu_only": True,
            "offline": True,
        }
        _write_json(temporary / "model_card.json", model_card)
        manifest = {
            "schema_version": 1,
            "model_id": model_card["model_id"],
            "created_at": datetime.now(timezone.utc).isoformat(),
            "dataset_manifest_sha256": sha256_file(dataset_dir / "manifest.json"),
            "benchmark_metrics_sha256": sha256_file(benchmark_metrics),
            "artifacts": {},
        }
        for path in sorted(temporary.iterdir()):
            manifest["artifacts"][path.name] = {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        _write_json(temporary / "manifest.json", manifest)
        temporary.replace(output_dir)
        return manifest
    except Exception:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


def _load_hmms(path: Path) -> list[pyhmmer.plan7.HMM]:
    hmms: list[pyhmmer.plan7.HMM] = []
    with pyhmmer.plan7.HMMFile(path) as handle:
        while True:
            hmm = handle.read()
            if hmm is None:
                break
            hmms.append(hmm)
    return hmms


def predict_proteins(
    model_dir: Path, fasta_path: Path, output_path: Path, *, threads: int = 6
) -> dict[str, object]:
    """Predict Cas membership and current hierarchy labels from an amino-acid FASTA."""
    manifest = json.loads((model_dir / "manifest.json").read_text(encoding="utf-8"))
    metadata = json.loads((model_dir / "profile_metadata.json").read_text(encoding="utf-8"))
    hmms = _load_hmms(model_dir / "profiles.hmm")
    queries = [
        {"sequence_id": name, "sequence": sequence} for name, sequence in parse_fasta(fasta_path)
    ]
    started = time.perf_counter()
    predictions = scan_profiles(queries, hmms, metadata, threads)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cas_count = 0
    with output_path.open("w", encoding="utf-8") as handle:
        for query, prediction in zip(queries, predictions, strict=True):
            is_cas = float(prediction["margin"]) >= 0.0
            cas_count += is_cas
            classification = (
                dict(prediction["classification"])
                if is_cas
                else {"class": None, "type": None, "subtype": None}
            )
            if classification["type"] == "V":
                classification["subtype"] = None
            record = {
                "sequence_id": query["sequence_id"],
                "sequence_length": len(str(query["sequence"])),
                "is_cas": is_cas,
                "positive_profile_score": prediction["positive_score"],
                "hard_negative_profile_score": prediction["negative_score"],
                "score_margin": prediction["margin"],
                "best_positive_profile": prediction["best_positive_profile"],
                "classification": classification,
                "variant": None,
                "evidence": {
                    "model_id": manifest["model_id"],
                    "profile_hits": prediction["hit_count"],
                    "type_v_subtype_masked": classification["type"] == "V",
                },
            }
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")
    elapsed = time.perf_counter() - started
    return {
        "model_id": manifest["model_id"],
        "input_records": len(queries),
        "cas_predictions": cas_count,
        "cpu_threads": threads,
        "wall_seconds": elapsed,
        "proteins_per_second": len(queries) / elapsed if elapsed else None,
        "output": str(output_path),
    }
