"""Paper-era CasAndra CPU bundle with expanded cassette deployment refit."""

from __future__ import annotations

import json
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from casandra.utils import sha256_file


BUNDLE_ID = "casandra-cas-only-cpu-bundle-v3-paper"


def build_bundle_v3(
    protein_model: Path,
    production_architecture: Path,
    evaluation_hybrid: Path,
    evaluation_rule_metrics: Path,
    taxonomy: Path,
    output_dir: Path,
) -> dict[str, object]:
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace bundle: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent))
    try:
        (staging / "protein").mkdir()
        (staging / "cassette_architecture").mkdir()
        for path in protein_model.iterdir():
            if path.is_file():
                shutil.copy2(path, staging / "protein" / path.name)
        for name in ("model.joblib", "model_card.json", "manifest.json"):
            shutil.copy2(production_architecture / name, staging / "cassette_architecture" / name)
        shutil.copy2(evaluation_hybrid / "model_card.json", staging / "cassette_model_card.json")
        shutil.copy2(evaluation_hybrid / "metrics.json", staging / "cassette_evaluation_metrics.json")
        shutil.copy2(taxonomy, staging / "type_v_taxonomy.json")
        rules = json.loads(evaluation_rule_metrics.read_text(encoding="utf-8"))
        config = {
            "schema_version": 3,
            "bundle_id": BUNDLE_ID,
            "protein_model_directory": "protein",
            "cassette_architecture_directory": "cassette_architecture",
            "cassette_parameters": rules["selected_parameters"],
            "hybrid_routing": "Type III architecture ExtraTrees; direct profiles otherwise",
            "genome_evidence_gate": {
                "supported_subtype_required": True,
                "class_1_minimum_cas_genes": 2,
                "class_2_minimum_cas_genes_or_single_effector": 2,
                "class_2_single_effector_minimum_aa": 400,
                "class_2_single_effector_minimum_positive_profile_score": 50.0,
                "selection": "pre-expanded validation and whole-genome acceptance; frozen before paper test",
            },
            "scope": "Cas proteins and Cas cassettes only",
            "crispr_array_prediction": False,
            "cpu_only": True,
            "offline_inference": True,
            "evaluation_model_separate_from_deployment_refit": True,
        }
        (staging / "config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        manifest = {
            "schema_version": 3,
            "bundle_id": BUNDLE_ID,
            "created_at": datetime.now(timezone.utc).isoformat(),
            "cpu_only": True,
            "offline_inference": True,
            "inputs": {
                "protein_manifest_sha256": sha256_file(protein_model / "manifest.json"),
                "production_architecture_manifest_sha256": sha256_file(
                    production_architecture / "manifest.json"
                ),
                "evaluation_hybrid_manifest_sha256": sha256_file(
                    evaluation_hybrid / "manifest.json"
                ),
                "evaluation_rule_metrics_sha256": sha256_file(evaluation_rule_metrics),
                "taxonomy_sha256": sha256_file(taxonomy),
            },
            "files": {},
        }
        for path in sorted(staging.rglob("*")):
            if path.is_file() and path.name != "manifest.json":
                manifest["files"][str(path.relative_to(staging))] = {
                    "size": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
        # Include nested model manifests as files; omit only the top-level one.
        for path in sorted(staging.rglob("manifest.json")):
            manifest["files"][str(path.relative_to(staging))] = {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        staging.replace(output_dir)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return manifest
