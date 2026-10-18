from __future__ import annotations

import json
from pathlib import Path
import sys
import zipfile

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import manage_release
import provision_tools
from runtime.common import load_lock


def tool_lock(*optional):
    names = ["CRISPRidentify-v2", "CasAndra", "CRISPRtracrRNA-v3", "SpacerPlacer", "evOr", *optional]
    return {
        "schema_version": 1,
        "tools": [
            {
                "name": name,
                "path": f"tools/{name}",
                "repository": f"https://github.com/Alexander-Mitrofanov/{name}.git",
                "commit": "a" * 40,
                "version": "1.0.0",
                "entrypoint": {"CRISPRrepeat": "crisprrepeat", "CRISPRmap-v2": "crisprmap", "RepeatTyper": "repeattyper", "CCTK": "cctk"}.get(name, name),
            }
            for name in names
        ],
    }


@pytest.mark.parametrize("optional", [
    ("CCTK",), ("RepeatTyper",), ("CRISPRleader-v2", "CRISPRrepeat", "CRISPRspacer", "CRISPRmap-v2", "RepeatTyper"),
    ("CRISPRleader-v2", "CRISPRrepeat", "CRISPRspacer", "CRISPRmap-v2", "RepeatTyper", "SpacePHARER", "CCTK"),
    (), ("CRISPRleader-v2",), ("CRISPRrepeat",), ("CRISPRmap-v2",),
    ("CRISPRleader-v2", "CRISPRrepeat"),
    ("CRISPRleader-v2", "CRISPRrepeat", "CRISPRspacer"),
    ("CRISPRleader-v2", "CRISPRrepeat", "CRISPRspacer", "CRISPRmap-v2"),
])
def test_supported_optional_tools_preserve_baseline(tmp_path, optional):
    path = tmp_path / "tools.lock.json"
    lock = tool_lock(*optional)
    path.write_text(json.dumps(lock))
    assert load_lock(path) == lock


@pytest.mark.parametrize("mutation", ["duplicate", "unknown", "missing_baseline", "foreign_repository"])
def test_optional_tool_lock_still_rejects_invalid_identity(tmp_path, mutation):
    lock = tool_lock("CRISPRrepeat")
    if mutation == "duplicate":
        lock["tools"].append(lock["tools"][-1])
    elif mutation == "unknown":
        lock["tools"][-1]["name"] = "UnknownTool"
    elif mutation == "missing_baseline":
        lock["tools"].pop(0)
    else:
        lock["tools"][-1]["repository"] = "https://github.com/another/CRISPRrepeat.git"
    path = tmp_path / "tools.lock.json"
    path.write_text(json.dumps(lock))
    with pytest.raises(ValueError):
        load_lock(path)


@pytest.mark.parametrize("include_leader", [False, True])
@pytest.mark.parametrize("include_map", [False, True])
@pytest.mark.parametrize("include_type", [False, True])
@pytest.mark.parametrize("include_cctk", [False, True])
def test_repeat_toolset_has_pinned_isolated_runtime_and_optional_runner(
    tmp_path, monkeypatch, include_leader, include_map, include_type, include_cctk
):
    lock = tool_lock(*(["CRISPRleader-v2"] if include_leader else []), "CRISPRrepeat", "CRISPRspacer",
                     *(["CRISPRmap-v2"] if include_map else []), *(["RepeatTyper"] if include_type else []),
                     *(["CCTK"] if include_cctk else []))
    path = tmp_path / "tools.lock.json"
    path.write_text(json.dumps(lock))
    runtime = tmp_path / "runtime"
    commands = []

    def export(_tool, destination):
        destination.mkdir()
        if destination.name == "CRISPRtracrRNA-v3":
            (destination / "upstream-lock.json").write_text("{}")

    def environment(prefix, *, conda=False):
        (prefix / "bin").mkdir(parents=True, exist_ok=True)
        for executable in ("python", "IntaRNA"):
            (prefix / "bin" / executable).write_text("test executable")
        (prefix / "lib").mkdir(exist_ok=True)
        (prefix / "lib/_RNA.so").write_text("test native extension")
        if conda:
            (prefix / "conda-meta").mkdir(exist_ok=True)

    def run(command, **_kwargs):
        command = [str(part) for part in command]
        commands.append(command)
        if command[:2] == ["/test/conda", "create"]:
            environment(Path(command[command.index("--prefix") + 1]), conda=True)
        elif "venv" in command:
            environment(Path(command[-1]))
        elif command[0] == "bash":
            environment(Path(command[1]).parent / ".envs/crispridentify-core", conda=True)

    def output(command, **_kwargs):
        if "-I" in command:
            prefix = Path(command[0]).parent.parent
            return json.dumps({"version": "2.7.0", "binding_version": "2.7.0", "path": str(prefix / "lib/_RNA.so")})
        return "# test environment receipt\n"

    monkeypatch.setattr(provision_tools, "RUNTIME", runtime)
    monkeypatch.setattr(provision_tools, "export_repository", export)
    monkeypatch.setattr(provision_tools.shutil, "which", lambda _name: "/test/conda")
    monkeypatch.setattr(provision_tools, "run", run)
    monkeypatch.setattr(provision_tools.subprocess, "check_output", output)
    target = provision_tools.build(path)
    runners = json.loads((target / "runners.json").read_text())
    record = json.loads((target / "toolset.json").read_text())["tools"]["CRISPRrepeat"]
    prefix = target / "envs/CRISPRrepeat"
    assert any(command[-3:] == ["-m", "venv", str(prefix)] for command in commands)
    repeat_install = [command for command in commands if command[0] == str(prefix / "bin/python") and "install" in command]
    assert len(repeat_install) == 1
    assert "ViennaRNA==2.7.0" in repeat_install[0]
    assert not any(str(prefix) in command and "create" in command for command in commands)
    assert record["environment"] == str(prefix)
    spacer_prefix = target / "envs/CRISPRspacer"
    assert any(command[-3:] == ["-m", "venv", str(spacer_prefix)] for command in commands)
    assert not any(str(spacer_prefix) in command and "create" in command for command in commands)
    assert runners["crisprspacer"] == [str(target / "bin/CRISPRspacer")]
    # CRISPRmap gets an isolated web runner only when present in the lock.
    assert ("crisprmap" in runners) is include_map
    assert "CRISPRmap-v2" not in runners
    if include_map:
        assert runners["crisprmap"] == [str(target / "bin/CRISPRmap-v2")]
        map_prefix = target / "envs/CRISPRmap-v2"
        map_record = json.loads((target / "toolset.json").read_text())["tools"]["CRISPRmap-v2"]
        assert any(command[-3:] == ["-m", "venv", str(map_prefix)] for command in commands)
        assert not any(str(map_prefix) in command and "create" in command for command in commands)
        map_install = [command for command in commands
                       if command[0] == str(map_prefix / "bin/python") and "install" in command]
        assert map_install == [[str(map_prefix / "bin/python"), "-m", "pip", "install",
                                str(target / "sources/CRISPRmap-v2")]]
        assert map_record["environment"] == str(map_prefix)
        assert "native_tools" not in map_record
        assert str(map_prefix / "bin/crisprmap") in Path(map_record["command"][0]).read_text()
        launcher = Path(map_record["command"][0]).read_text()
        assert "unset MAFFT_BINARIES\n" in launcher
        assert "mafft" not in launcher.replace("unset MAFFT_BINARIES\n", "").lower()
        assert (target / "CRISPRmap-v2.requirements.txt").is_file()
        assert not (target / "CRISPRmap-v2.conda-explicit.txt").exists()
    assert ("repeattyper" in runners) is include_type
    if include_type:
        typing_prefix = target / "envs/RepeatTyper"
        assert any(command[-3:] == ["-m", "venv", str(typing_prefix)] for command in commands)
        typing_install = [command for command in commands if command[0] == str(typing_prefix / "bin/python") and "install" in command]
        assert typing_install == [[str(typing_prefix / "bin/python"), "-m", "pip", "install", "-r",
                                   str(target / "sources/RepeatTyper/requirements.lock"), str(target / "sources/RepeatTyper")]]
        assert runners["repeattyper"] == [str(target / "bin/RepeatTyper")]
        assert (target / "RepeatTyper.requirements.txt").is_file()
    assert ("cctk" in runners) is include_cctk
    if include_cctk:
        comparison_prefix = target / "envs/CCTK"
        assert any(command[-3:] == ["-m", "venv", str(comparison_prefix)] for command in commands)
        assert not any(str(comparison_prefix) in command and "create" in command for command in commands)
        installs = [command for command in commands if command[0] == str(comparison_prefix / "bin/python") and "install" in command]
        assert installs == [[str(comparison_prefix / "bin/python"), "-m", "pip", "install", "-r",
                             str(target / "sources/CCTK/requirements.lock"), str(target / "sources/CCTK")]]
        assert runners["cctk"] == [str(target / "bin/CCTK")]
        assert (target / "CCTK.requirements.txt").is_file()
        assert not (target / "CCTK.conda-explicit.txt").exists()
        assert str(comparison_prefix / "bin/cctk") in (target / "bin/CCTK").read_text()
    assert runners["crisprrepeat"] == record["command"]
    assert ("crisprleader" in runners) is include_leader
    assert record["native_tools"]["ViennaRNA"]["version"] == "2.7.0"
    assert record["native_tools"]["ViennaRNA"]["sha256"]
    assert str(prefix / "bin/crisprrepeat") in Path(record["command"][0]).read_text()
    assert (target / "CRISPRrepeat.requirements.txt").is_file()
    assert not (runtime / "current").exists()
    monkeypatch.setattr(manage_release, "RUNTIME", runtime)
    document = manage_release.verify_toolset(target)
    # Provisioning an already built toolset must apply the same verification as
    # application release construction, rather than returning a stale runtime.
    installed_code = prefix / "lib/installed.py"
    installed_code.write_text("changed installed code")
    with pytest.raises(ValueError, match="environment changed"):
        provision_tools.build(path)
    installed_code.unlink()
    if include_cctk:
        # Installed code participates in the same immutable receipt as native
        # engines, even though observed comparison needs no external executable.
        code = comparison_prefix / "lib/comparison.py"
        code.write_text("changed code")
        with pytest.raises(ValueError, match="environment changed"):
            manage_release.verify_toolset(target)
        code.unlink()
    document["tools"]["CRISPRrepeat"]["native_tools"]["ViennaRNA"]["version"] = "2.7.2"
    (target / "toolset.json").write_text(json.dumps(document))
    with pytest.raises(ValueError, match="ViennaRNA identity"):
        manage_release.verify_toolset(target)
    (target / "toolset.json").unlink()
    (prefix / "lib/_RNA.so").write_text("changed native extension with the same version label")
    with pytest.raises(ValueError, match="receipt changed before resume"):
        provision_tools.build(path, resume=True)


@pytest.mark.parametrize("mutation", ["commit", "source", "environment", "launcher", "receipt"])
def test_resume_refuses_changed_installation_receipt(tmp_path, monkeypatch, mutation):
    from runtime.common import atomic_json, digest, environment_inventory, lock_id, source_inventory

    lock = tool_lock()
    path = tmp_path / "tools.lock.json"
    path.write_text(json.dumps(lock))
    runtime = tmp_path / "runtime"
    target = runtime / "toolsets" / lock_id(lock)
    tool = lock["tools"][0]
    package = target / "sources" / tool["name"]
    prefix = package / ".envs/crispridentify-core"
    (prefix / "bin").mkdir(parents=True)
    code = prefix / "bin/python"
    code.write_text("original environment")
    source = package / "code.py"
    source.write_text("original source")
    upstream = target / "sources/CRISPRtracrRNA-v3"
    upstream.mkdir(parents=True)
    (upstream / "upstream-lock.json").write_text("{}")
    launcher = target / "bin" / tool["name"]
    launcher.parent.mkdir()
    launcher.write_text("original launcher")
    receipt = target / f".{tool['name']}.environment.json"
    atomic_json(receipt, environment_inventory(prefix))
    record = {**tool, "source": str(package), "environment": str(prefix), "command": [str(launcher)],
              "launcher_sha256": digest(launcher), "source_files": source_inventory(package),
              "environment_sha256": digest(receipt)}
    if mutation == "commit":
        record["commit"] = "b" * 40
    elif mutation == "source":
        source.write_text("changed source")
    elif mutation == "environment":
        code.write_text("changed environment")
    elif mutation == "launcher":
        launcher.write_text("changed launcher")
    else:
        receipt.unlink()
    atomic_json(target / f".{tool['name']}.installed.json", record)
    monkeypatch.setattr(provision_tools, "RUNTIME", runtime)
    monkeypatch.setattr(provision_tools.shutil, "which", lambda _name: "/test/conda")
    monkeypatch.setattr(provision_tools, "export_repository", lambda _tool, destination: destination.mkdir())
    monkeypatch.setattr(provision_tools, "run", lambda *_args, **_kwargs: pytest.fail("Must reject before installation"))
    with pytest.raises(ValueError, match="receipt changed before resume"):
        provision_tools.build(path, resume=True)


@pytest.mark.parametrize("wrong_field", ["version", "binding_version"])
def test_repeat_runtime_refuses_changed_fold_engine(tmp_path, monkeypatch, wrong_field):
    def output(_command, **_kwargs):
        identity = {"version": "2.7.0", "binding_version": "2.7.0", "path": str(tmp_path / "_RNA.so")}
        identity[wrong_field] = "2.7.2"
        return json.dumps(identity)

    monkeypatch.setattr(provision_tools.subprocess, "check_output", output)
    with pytest.raises(ValueError, match="requires the ViennaRNA 2.7.0"):
        provision_tools.verify_repeat_runtime(tmp_path, {})


@pytest.mark.parametrize(
    ("optional", "version"),
    [((), "2.0.0"), (("crisprleader",), "3.0.0"), (("crisprrepeat",), "4.0.0"), (("crisprleader", "crisprrepeat"), "4.0.0"), (("crisprleader", "crisprrepeat", "crisprspacer"), "5.0.0"), (("crisprmap",), "6.0.0"), (("crisprmap", "crisprspacer"), "6.0.0"), (("crisprmap", "repeattyper"), "7.0.0"), (("cctk",), "9.0.0"), (("cctk", "spacepharer", "repeattyper"), "9.0.0")],
)
def test_release_schema_preserves_optional_leader(optional, version):
    assert manage_release.scientific_manifest_version(dict.fromkeys(optional)) == version


def test_release_version_comes_from_exact_integration_wheel(tmp_path):
    path = tmp_path / "integration.whl"
    with zipfile.ZipFile(path, "w") as wheel:
        wheel.writestr("crispr_tool_integration-0.5.0.dist-info/METADATA", "Name: crispr-tool-integration\nVersion: 0.5.0\n")
    assert manage_release.wheel_version(path, "crispr-tool-integration") == "0.5.0"
    with pytest.raises(ValueError, match="identity"):
        manage_release.wheel_version(path, "another-package")


def test_ambiguous_wheel_metadata_is_rejected(tmp_path):
    path = tmp_path / "ambiguous.whl"
    with zipfile.ZipFile(path, "w") as wheel:
        for name in ("one", "two"):
            wheel.writestr(f"{name}.dist-info/METADATA", "Name: crispr-tool-integration\nVersion: 0.5.0\n")
    with pytest.raises(ValueError, match="ambiguous"):
        manage_release.wheel_version(path, "crispr-tool-integration")
