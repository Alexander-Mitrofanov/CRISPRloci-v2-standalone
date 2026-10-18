"""Build a single immutable offline CasAndra model bundle."""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from casandra.utils import sha256_file


BUNDLE_ID = "casandra-cas-only-cpu-bundle-v1"


def build_bundle(
    protein_model: Path,
    architecture_model: Path,
    hybrid_model: Path,
    cassette_rule_metrics: Path,
    taxonomy: Path,
    output_dir: Path,
):
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace model bundle: {output_dir}")
    (output_dir / "protein").mkdir(parents=True)
    (output_dir / "cassette_architecture").mkdir()
    for path in protein_model.iterdir():
        if path.is_file():
            shutil.copy2(path, output_dir / "protein" / path.name)
    for path in architecture_model.iterdir():
        if path.is_file():
            shutil.copy2(path, output_dir / "cassette_architecture" / path.name)
    shutil.copy2(hybrid_model / "model_card.json", output_dir / "cassette_model_card.json")
    shutil.copy2(hybrid_model / "metrics.json", output_dir / "cassette_internal_metrics.json")
    shutil.copy2(taxonomy, output_dir / "type_v_taxonomy.json")
    rules = json.loads(cassette_rule_metrics.read_text(encoding="utf-8"))
    config = {
        "schema_version": 1,
        "bundle_id": BUNDLE_ID,
        "protein_model_directory": "protein",
        "cassette_architecture_directory": "cassette_architecture",
        "cassette_parameters": rules["selected_parameters"],
        "hybrid_routing": "Type III architecture ExtraTrees; direct profiles otherwise",
        "scope": "Cas proteins and Cas cassettes only",
        "crispr_array_prediction": False,
        "cpu_only": True,
        "offline_inference": True,
    }
    (output_dir / "config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    files = [path for path in output_dir.rglob("*") if path.is_file()]
    manifest = {
        "schema_version": 1,
        "bundle_id": BUNDLE_ID,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "cpu_only": True,
        "offline_inference": True,
        "files": {
            str(path.relative_to(output_dir)): {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(files)
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest
