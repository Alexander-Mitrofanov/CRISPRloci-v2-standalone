"""Ensure the distribution check detects missing and damaged scientific assets."""

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("wheel_check", ROOT / "scripts/check_wheel.py")
wheel_check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(wheel_check)


@pytest.fixture
def package(tmp_path):
    package = tmp_path / "package"
    package.mkdir()
    (package / "py.typed").touch()
    names = json.loads((ROOT / "src/crisprtracrrna/data/manifest.json").read_text())
    manifest = {}
    for relative in names:
        path = package / "data" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(relative.encode())
        manifest[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    (package / "data/manifest.json").write_text(json.dumps(manifest))
    (package / "data/erpin/erpin").chmod(0o755)
    return package


def test_resource_check_accepts_complete_intact_distribution(package):
    assert wheel_check.verify_resources(package) == 27


@pytest.mark.parametrize("filename", ["py.typed", "data/models/II.cm.i1f"])
def test_resource_check_rejects_missing_required_files(package, filename):
    (package / filename).unlink()
    with pytest.raises(RuntimeError, match="missing"):
        wheel_check.verify_resources(package)


def test_resource_check_rejects_incomplete_manifest(package):
    path = package / "data/manifest.json"
    manifest = json.loads(path.read_text())
    del manifest["models/II.cm"]
    path.write_text(json.dumps(manifest))
    with pytest.raises(RuntimeError, match="omits required resources"):
        wheel_check.verify_resources(package)


def test_resource_check_rejects_corrupted_asset(package):
    (package / "data/models/II.cm").write_bytes(b"corrupted model")
    with pytest.raises(RuntimeError, match="checksum mismatch: models/II.cm"):
        wheel_check.verify_resources(package)


@pytest.mark.skipif(sys.platform != "linux", reason="Bundled executable targets Linux")
def test_resource_check_rejects_lost_executable_permission(package):
    (package / "data/erpin/erpin").chmod(0o644)
    with pytest.raises(RuntimeError, match="lost executable permission"):
        wheel_check.verify_resources(package)
