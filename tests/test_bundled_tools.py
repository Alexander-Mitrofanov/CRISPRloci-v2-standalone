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
