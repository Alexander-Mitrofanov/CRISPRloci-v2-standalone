"""Reference-backed modes accept reviewed pins at startup, not only at output validation."""

import hashlib
import json
from pathlib import Path

import pytest


@pytest.mark.parametrize(("tool", "version", "supported"), [
    ("crisprspacer", "0.3.0", True), ("crisprspacer", "0.3.1", True),
    ("crisprspacer", "0.2.0", False), ("crisprspacer", "9.9.9", False),
    ("spacepharer", "2.1.0", True), ("spacepharer", "2.1.1", True),
    ("spacepharer", "2.0.0", False), ("spacepharer", "9.9.9", False),
])
def test_reference_startup_accepts_only_supported_tool_patches(settings_factory, tmp_path, tool, version, supported):
    root = Path(__file__).resolve().parents[2]
    lock = json.loads((root / "tools/tools.lock.json").read_text())
    names = {
        "clinker-v2": "clinker", "Phigaro-v2": "phigaro", "ISEScan-v2": "isescan",
        "DefenseFinder-v2": "defensefinder", "RNAalifold-v2": "rnaalifold",
        "CRISPRidentify-v2": "crispridentify", "CasAndra": "casandra",
        "CRISPRtracrRNA-v3": "crispr_tracrrna", "SpacerPlacer": "spacerplacer", "evOr": "crispr_evor",
        "CRISPRleader-v2": "crisprleader", "CRISPRrepeat": "crisprrepeat", "CRISPRspacer": "crisprspacer",
        "CRISPRmap-v2": "crisprmap", "RepeatTyper": "repeattyper", "SpacePHARER": "spacepharer", "CCTK": "cctk",
    }
    artifacts = {
        names[row["name"]]: {
            "name": names[row["name"]], "display_name": row["name"], "version": row["version"],
            "kind": "git_source", "path": "/source/" + row["name"], "sha256": "a" * 64,
            "sha256_scope": "source_inventory_v1", "repository": row["repository"], "commit": row["commit"],
        }
        for row in lock["tools"]
    }
    artifacts[tool]["version"] = version
    artifacts["integration"] = {
        "name": "integration", "display_name": "Integration", "version": "0.14.0", "kind": "wheel",
        "path": "/source/integration.whl", "sha256": "b" * 64, "sha256_scope": "file_bytes",
    }
    manifest = tmp_path / "reviewed-manifest.json"
    payload = json.dumps({"schema": "crispr-web-scientific-release", "schema_version": "10.0.0",
                          "release_id": "reviewed", "scientific_artifacts": artifacts}).encode()
    manifest.write_bytes(payload)
    reference = tmp_path / "reference"
    reference.mkdir()
    field = "viral_reference" if tool == "crisprspacer" else "association_reference"
    options = {"scientific_release_manifest": manifest,
               "scientific_release_manifest_sha256": hashlib.sha256(payload).hexdigest(),
               field + "_path": reference, field + "_sha256": "c" * 64}
    if supported:
        settings = settings_factory(**options)
        assert settings.scientific_release_identities[tool]["version"] == version
    else:
        with pytest.raises(ValueError, match="requires release-pinned"):
            settings_factory(**options)
