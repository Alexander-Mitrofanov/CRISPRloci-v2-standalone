"""Add the locked genome evidence gate to the final offline bundle."""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from casandra.utils import sha256_file


BUNDLE_ID = "casandra-cas-only-cpu-bundle-v2"


def build_bundle_v2(source_bundle: Path, output_dir: Path):
    if output_dir.exists():
        raise FileExistsError(f"Refusing to replace bundle: {output_dir}")
    shutil.copytree(source_bundle, output_dir)
    config_path = output_dir / "config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["bundle_id"] = BUNDLE_ID
    config["genome_evidence_gate"] = {
        "supported_subtype_required": True,
        "class_1_minimum_cas_genes": 2,
        "class_2_minimum_cas_genes_or_single_effector": 2,
        "class_2_single_effector_minimum_aa": 400,
        "class_2_single_effector_minimum_positive_profile_score": 50.0,
        "selection": "validation positive-window recall plus whole-genome acceptance precision",
    }
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = output_dir / "manifest.json"
    manifest = {
        "schema_version": 2,
        "bundle_id": BUNDLE_ID,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "supersedes": "v1 emitted ungated candidate cassettes",
        "cpu_only": True,
        "offline_inference": True,
        "files": {
            str(path.relative_to(output_dir)): {
                "size": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(output_dir.rglob("*"))
            if path.is_file() and path != manifest_path
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
