"""CCTK source publication and release admission stay exact and offline."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from runtime.common import export_repository, load_lock
import manage_release
from crispr_web.config import _validated_scientific_manifest

ROOT = Path(__file__).resolve().parents[1]

def test_cctk_locked_pin_exports_complete_installed_source(tmp_path):
    lock = load_lock(ROOT / "tools/tools.lock.json")
    assert len(lock["tools"]) >= 12
    entry = next(tool for tool in lock["tools"] if tool["name"] == "CCTK")
    assert entry["version"] == "2.0.2" and entry["entrypoint"] == "cctk"
    export = tmp_path / "CCTK"
    export_repository(entry, export)
    if (ROOT / "tools/bundled-sources.json").is_file():
        from bundled_tools import source_tree, verify_bundled

        assert verify_bundled(ROOT, entry)
        assert source_tree(export) == source_tree(ROOT / entry["path"])
        assert json.loads((export / "UPSTREAM.json").read_text())["commit"] == "fe94fde10a2855cbcc3bf118a253ec097d3a5e46"
        return
    names = subprocess.check_output(
        ["git", "-C", str(ROOT / entry["path"]), "ls-tree", "-r", "--name-only", "-z", entry["commit"]]
    ).decode().rstrip("\0").split("\0")
    assert {str(p.relative_to(export)) for p in export.rglob("*") if p.is_file()} == set(names)
    for name in ("src/cctk/comparison.py", "src/cctk/_exec.py", "requirements.lock", "tests/fixtures/provenance.json"):
        original = subprocess.check_output(
            ["git", "-C", str(ROOT / entry["path"]), "show", entry["commit"] + ":" + name]
        )
        assert (export / name).read_bytes() == original
    assert json.loads((export / "UPSTREAM.json").read_text())["commit"] == "fe94fde10a2855cbcc3bf118a253ec097d3a5e46"


def manifest(*optional):
    artifacts = {}
    names = {"crispridentify", "casandra", "crispr_tracrrna", "spacerplacer", "crispr_evor", "cctk", *optional}
    for name in names:
        repository_name = manage_release.NAMES[name]
        artifacts[name] = {
            "name": name, "display_name": repository_name, "version": "2.0.1",
            "kind": "git_source", "path": "/runtime/sources/" + repository_name,
            "sha256": "a" * 64, "sha256_scope": "source_inventory_v1",
            "repository": f"https://github.com/Alexander-Mitrofanov/{repository_name}.git",
            "commit": "b" * 40,
        }
    artifacts["integration"] = {
        "name": "integration", "display_name": "Integration", "version": "0.12.0",
        "kind": "wheel", "path": "/runtime/integration.whl", "sha256": "c" * 64,
        "sha256_scope": "file_bytes",
    }
    return {"schema": "crispr-web-scientific-release", "schema_version": "9.0.0",
            "release_id": "test-release", "scientific_artifacts": artifacts}


@pytest.mark.parametrize("optional", [(), ("spacepharer",), ("crisprleader", "crisprrepeat", "crisprspacer", "crisprmap", "repeattyper", "spacepharer")])
def test_manifest_nine_requires_cctk_without_requiring_other_optional_tools(optional):
    document = manifest(*optional)
    identities = _validated_scientific_manifest(document)
    assert set(identities) == set(document["scientific_artifacts"])
    assert identities["cctk"]["repository"].endswith("/CCTK.git")
    assert "path" not in identities["cctk"]
    assert manage_release.scientific_manifest_version(identities) == "9.0.0"
    document["scientific_artifacts"].pop("cctk")
    with pytest.raises(ValueError, match="artifact set"):
        _validated_scientific_manifest(document)


@pytest.mark.parametrize("mutation", ["old_schema", "missing_baseline", "bad_commit", "wrong_kind", "extra_identity"])
def test_cctk_release_rejects_incomplete_or_forged_identity(mutation):
    document = manifest("spacepharer")
    if mutation == "old_schema":
        document["schema_version"] = "8.0.0"
    elif mutation == "missing_baseline":
        document["scientific_artifacts"].pop("casandra")
    elif mutation == "bad_commit":
        document["scientific_artifacts"]["cctk"]["commit"] = "main"
    elif mutation == "wrong_kind":
        document["scientific_artifacts"]["cctk"]["kind"] = "wheel"
    else:
        document["scientific_artifacts"]["cctk"]["extra"] = "unexpected"
    with pytest.raises(ValueError):
        _validated_scientific_manifest(document)
