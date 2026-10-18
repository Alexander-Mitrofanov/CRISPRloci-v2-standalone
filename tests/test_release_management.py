from __future__ import annotations

import fcntl
import json
from pathlib import Path
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import manage_release
from runtime.common import atomic_json, digest, environment_inventory, source_inventory


def test_runtime_receipt_detects_code_link_and_permission_changes(tmp_path):
    package = tmp_path / "lib"
    package.mkdir()
    code = package / "tool.py"
    code.write_text("print(1)\n")
    link = tmp_path / "entry"
    link.symlink_to(code)
    before = environment_inventory(tmp_path)
    cache = package / "__pycache__"
    cache.mkdir()
    (cache / "tool.pyc").write_bytes(b"cache")
    assert environment_inventory(tmp_path) == before
    code.write_text("print(2)\n")
    assert environment_inventory(tmp_path) != before
    code.write_text("print(1)\n")
    code.chmod(0o755)
    assert environment_inventory(tmp_path) != before


def test_source_snapshots_exclude_local_secrets_and_builds(tmp_path):
    for name in (".env", ".env.production", ".env.example", "code.py"):
        (tmp_path / name).write_text("value")
    (tmp_path / ".venv").mkdir()
    (tmp_path / ".venv/private").write_text("secret")
    assert set(source_inventory(tmp_path)) == {".env.example", "code.py"}


def test_activation_requires_matching_gate_and_stopped_services(tmp_path, monkeypatch):
    root = tmp_path / "project"
    runtime = root / ".runtime"
    release = runtime / "releases" / "candidate"
    release.mkdir(parents=True)
    toolset = runtime / "toolsets" / "example"
    toolset.mkdir(parents=True)
    lock = json.loads(
        (Path(__file__).resolve().parents[1] / "tools/tools.lock.json").read_text()
    )
    atomic_json(toolset / "tools.lock.json", lock)
    atomic_json(release / "release.json", {"toolset": str(toolset)})
    monkeypatch.setattr(manage_release, "ROOT", root)
    monkeypatch.setattr(manage_release, "RUNTIME", runtime)
    monkeypatch.setattr(manage_release, "verify", lambda path: {})
    with pytest.raises(ValueError, match="gate_release"):
        manage_release.activate(release)
    atomic_json(release / "gate.json", {"status": "passed", "release_sha256": "0" * 64})
    with pytest.raises(ValueError, match="does not match"):
        manage_release.activate(release)
    atomic_json(
        release / "gate.json",
        {"status": "passed", "release_sha256": digest(release / "release.json")},
    )
    with (runtime / "active.lock").open("a") as service_lock:
        fcntl.flock(service_lock, fcntl.LOCK_SH)
        with pytest.raises(ValueError, match="running"):
            manage_release.activate(release)
        assert not (runtime / "current").exists()
    manage_release.activate(release)
    assert (runtime / "current").resolve() == release
    assert json.loads((root / "tools/tools.lock.json").read_text()) == lock
    predecessor = runtime / "releases" / "previous-build"
    predecessor.mkdir()
    (runtime / "previous").symlink_to(predecessor)
    manage_release.activate(release)
    assert (runtime / "previous").resolve() == predecessor


def test_candidate_reference_rejects_command_options_before_fetch():
    from update_tools import candidate

    with pytest.raises(ValueError, match="NAME@"):
        candidate(["CasAndra@--upload-pack=bad"])


@pytest.mark.parametrize("version_file,readme,expected", [
    ("2.3.2\n", "# CRISPRidentify v2\n", "2.3.2"),
    (None, "# CRISPRidentify v2.3.1\n", "2.3.1"),
    (None, "# CRISPRidentify v2\n", "2.3.2"),
])
def test_identify_candidate_version(tmp_path, monkeypatch, version_file, readme, expected):
    import update_tools

    lock = json.loads((update_tools.ROOT / "tools/tools.lock.json").read_text())
    entries = {entry["name"]: entry for entry in lock["tools"]}

    def output(command, **kwargs):
        if command[-2:] == ["--get", "remote.origin.url"]:
            return entries[Path(command[2]).name]["repository"]
        if command[-1] == "FETCH_HEAD^{commit}":
            return entries[Path(command[2]).name]["commit"]
        if command[-1].endswith(":upstream-lock.json"):
            return json.dumps(entries)
        raise AssertionError(command)

    def run(command, **kwargs):
        path = command[-1].split(":")[-1]
        if path == "VERSION" and version_file is not None:
            return subprocess.CompletedProcess(command, 0, version_file, "")
        if path == "README.md":
            return subprocess.CompletedProcess(command, 0, readme, "")
        assert path in {"VERSION", "pyproject.toml"}
        return subprocess.CompletedProcess(command, 1, "", "missing")

    monkeypatch.setattr(update_tools, "run", lambda command: None)
    monkeypatch.setattr(update_tools.subprocess, "check_output", output)
    monkeypatch.setattr(update_tools.subprocess, "run", run)
    destination = update_tools.candidate(["CRISPRidentify-v2@main"], tmp_path / "candidate.json")
    saved = json.loads(destination.read_text())
    assert next(e for e in saved["tools"] if e["name"] == "CRISPRidentify-v2")["version"] == expected
