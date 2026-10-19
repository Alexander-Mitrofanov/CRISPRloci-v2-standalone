"""Publication sources must reproduce upstream trees and work without origins."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bundled_tools import export_bundled, source_tree, verify_bundled
from runtime import common
import sync_tools


@pytest.fixture
def bundle(tmp_path):
    root = tmp_path / "publication"
    source = root / "tools/Example"
    source.mkdir(parents=True)
    (source / "build").mkdir()
    (source / "build/native.cmake").write_text("# Required native build recipe\n")
    (source / "main.py").write_text("print('science')\n")
    (source / "run").write_text("#!/bin/sh\nexec python main.py\n")
    (source / "run").chmod(0o755)
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    subprocess.run(["git", "-C", str(source), "add", "."], check=True)
    tree = subprocess.check_output(["git", "-C", str(source), "write-tree"], text=True).strip()
    import shutil

    shutil.rmtree(source / ".git")
    entry = {
        "name": "Example", "path": "tools/Example", "commit": "1" * 40,
        "repository": "https://github.com/example/Example.git",
    }
    manifest = root / "tools/bundled-sources.json"
    manifest.write_text(json.dumps({"schema_version": 1, "tools": [{**entry, "tree": tree}]}))
    return root, source, entry, tree


def test_bundle_matches_git_and_exports_offline_with_native_recipes(bundle, tmp_path, monkeypatch):
    root, source, entry, tree = bundle
    assert source_tree(source).hex() == tree
    assert verify_bundled(root, entry)
    monkeypatch.setattr(sync_tools, "ROOT", root)
    monkeypatch.setattr(common, "ROOT", root)

    def no_git(*args, **kwargs):
        raise AssertionError("Bundled sources must not invoke Git or contact private origins")

    monkeypatch.setattr(subprocess, "check_output", no_git)
    monkeypatch.setattr(subprocess, "run", no_git)
    sync_tools.sync(entry, check=True)
    sync_tools.sync(entry, check=False)
    target = tmp_path / "export"
    common.export_repository(entry, target)
    assert source_tree(target).hex() == tree
    assert (target / "build/native.cmake").is_file()
    assert (target / "run").stat().st_mode & 0o111


@pytest.mark.parametrize("mutation", ["content", "permission", "extra", "missing", "symlink", "pin"])
def test_bundle_rejects_source_and_identity_mutations(bundle, tmp_path, mutation):
    root, source, entry, _ = bundle
    if mutation == "content":
        (source / "main.py").write_text("changed\n")
    elif mutation == "permission":
        (source / "run").chmod(0o644)
    elif mutation == "extra":
        (source / "extra.py").write_text("unrecorded\n")
    elif mutation == "missing":
        (source / "main.py").unlink()
    elif mutation == "symlink":
        (source / "extra").symlink_to(tmp_path)
    else:
        entry = {**entry, "commit": "2" * 40}
    with pytest.raises(ValueError):
        export_bundled(root, entry, tmp_path / "export")
    assert not (tmp_path / "export").exists()


def test_bundle_ignores_only_generated_python_caches(bundle):
    root, source, entry, _ = bundle
    (source / "__pycache__").mkdir()
    (source / "__pycache__/main.pyc").write_bytes(b"cache")
    assert verify_bundled(root, entry)


def test_development_checkout_keeps_original_git_workflow(tmp_path):
    assert not verify_bundled(tmp_path, {})


def test_patched_bundle_cannot_silently_fall_back_to_upstream_git(tmp_path):
    with pytest.raises(ValueError, match="Missing bundled source manifest"):
        verify_bundled(tmp_path, {"bundled_tree": "1" * 40})


def test_bundled_patch_identity_changes_runtime_lock_and_binds_export(bundle, tmp_path):
    root, source, entry, tree = bundle
    pinned = {**entry, "bundled_tree": tree}
    assert verify_bundled(root, pinned)
    before = common.lock_id({"tools": [pinned]})
    (source / "main.py").write_text("print('fixed science')\n")
    changed = source_tree(source).hex()
    manifest = root / "tools/bundled-sources.json"
    document = json.loads(manifest.read_text())
    document["tools"][0]["tree"] = changed
    manifest.write_text(json.dumps(document))
    with pytest.raises(ValueError, match="tree differs"):
        export_bundled(root, pinned, tmp_path / "stale")
    fixed = {**pinned, "bundled_tree": changed}
    assert common.lock_id({"tools": [fixed]}) != before
    assert export_bundled(root, fixed, tmp_path / "fixed")


@pytest.mark.parametrize("tree", [None, 42, "", "0" * 39, "G" * 40])
def test_invalid_bundled_tree_pin_is_rejected_before_sync(bundle, tree):
    _, _, entry, _ = bundle
    with pytest.raises(ValueError, match="bundled tree"):
        sync_tools.validate_entry({**entry, "bundled_tree": tree})
