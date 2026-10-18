#!/usr/bin/env python3
"""Build, verify, gate and select CRISPRloci v2 releases. Never changes a live job."""

from __future__ import annotations
import argparse
from email.parser import BytesParser
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time
import zipfile

from runtime.common import (
    ROOT,
    RUNTIME,
    atomic_json,
    digest,
    environment_inventory,
    load_lock,
    lock_id,
    run,
    source_inventory,
)

NAMES = {
    "crispridentify": "CRISPRidentify-v2",
    "casandra": "CasAndra",
    "crispr_tracrrna": "CRISPRtracrRNA-v3",
    "spacerplacer": "SpacerPlacer",
    "crispr_evor": "evOr",
    "crisprleader": "CRISPRleader-v2",
    "crisprrepeat": "CRISPRrepeat",
    "crisprspacer": "CRISPRspacer",
    "crisprmap": "CRISPRmap-v2",
    "repeattyper": "RepeatTyper",
    "spacepharer": "SpacePHARER",
    "cctk": "CCTK",
    "clinker": "clinker-v2",
    "phigaro": "Phigaro-v2",
    "isescan": "ISEScan-v2",
    "defensefinder": "DefenseFinder-v2",
    "rnaalifold": "RNAalifold-v2",
}


def wheel_version(path, expected_name):
    """Use the version of the exact installed wheel, not a release-script constant."""
    with zipfile.ZipFile(path) as wheel:
        metadata_files = [
            item for item in wheel.infolist()
            if item.filename.endswith(".dist-info/METADATA")
        ]
        if len(metadata_files) != 1 or metadata_files[0].file_size > 1_000_000:
            raise ValueError("Application wheel metadata is missing, ambiguous or oversized")
        metadata = BytesParser().parsebytes(wheel.read(metadata_files[0]))
    names, versions = metadata.get_all("Name", []), metadata.get_all("Version", [])
    if (
        len(names) != 1 or names[0].lower().replace("_", "-") != expected_name
        or len(versions) != 1 or not versions[0].strip()
    ):
        raise ValueError("Application wheel identity does not match the release")
    return versions[0].strip()


def scientific_manifest_version(artifacts):
    if any(name in artifacts for name in ("clinker", "phigaro", "isescan", "defensefinder", "rnaalifold")):
        return "10.0.0"
    if "cctk" in artifacts:
        return "9.0.0"
    if "spacepharer" in artifacts:
        return "8.0.0"
    if "repeattyper" in artifacts:
        return "7.0.0"
    if "crisprmap" in artifacts:
        return "6.0.0"
    if "crisprspacer" in artifacts:
        return "5.0.0"
    if "crisprrepeat" in artifacts:
        return "4.0.0"
    return "3.0.0" if "crisprleader" in artifacts else "2.0.0"


def verify_toolset(path):
    path = Path(path).resolve(strict=True)
    path.relative_to(RUNTIME / "toolsets")
    document = json.loads((path / "toolset.json").read_text())
    lock = load_lock(path / "tools.lock.json")
    if document["toolset_id"] != lock_id(lock) or path.name != lock_id(lock):
        raise ValueError("Toolset identity mismatch")
    if digest(path / "runners.json") != document["runners_sha256"]:
        raise ValueError("Runner configuration has changed")
    if digest(path / "environments.json") != document["environments_sha256"]:
        raise ValueError("Installed environment receipt changed")
    for prefix, inventory in json.loads(
        (path / "environments.json").read_text()
    ).items():
        candidate = (path / prefix).resolve(strict=True)
        candidate.relative_to(path)
        if environment_inventory(candidate) != inventory:
            raise ValueError("Installed scientific environment changed: " + prefix)
    for entry in lock["tools"]:
        tool = document["tools"][entry["name"]]
        if (tool["commit"], tool["repository"]) != (
            entry["commit"],
            entry["repository"],
        ):
            raise ValueError("Toolset revision does not match lock")
        source = Path(tool["source"])
        source.relative_to(path)
        if source_inventory(source) != tool["source_files"]:
            raise ValueError("Tool source snapshot changed: " + entry["name"])
        if digest(tool["command"][0]) != tool["launcher_sha256"]:
            raise ValueError("Tool launcher changed: " + entry["name"])
        if entry["name"] in {"clinker-v2", "Phigaro-v2", "ISEScan-v2", "DefenseFinder-v2", "RNAalifold-v2"}:
            from runtime.context_tools import verify_context_native_receipt
            prefix = Path(tool["environment"]).resolve(strict=True)
            prefix.relative_to(path)
            verify_context_native_receipt(tool, prefix)
        if entry["name"] == "SpacePHARER":
            prefix = Path(tool["environment"]).resolve(strict=True)
            prefix.relative_to(path)
            native = tool.get("native_tools", {}).get("spacepharer-native", {})
            expected_binary = prefix / "bin/spacepharer-native"
            if (set(native) != {"version", "path", "sha256"} or native["version"] != entry["version"]
                    or Path(native["path"]) != expected_binary or expected_binary.is_symlink()
                    or digest(expected_binary) != native["sha256"]):
                raise ValueError("SpacePHARER native engine receipt does not match the toolset")
        if entry["name"] == "CRISPRrepeat":
            prefix = Path(tool["environment"]).resolve(strict=True)
            prefix.relative_to(path)
            native_tools = tool.get("native_tools")
            if not isinstance(native_tools, dict) or set(native_tools) != {"ViennaRNA"}:
                raise ValueError("CRISPRrepeat ViennaRNA identity is missing from the toolset")
            native = native_tools["ViennaRNA"]
            if (
                not isinstance(native, dict)
                or set(native) != {"version", "binding_version", "path", "sha256"}
                or native["version"] != "2.7.0"
                or native["binding_version"] != "2.7.0"
            ):
                raise ValueError("CRISPRrepeat ViennaRNA identity does not match the toolset")
            extension = Path(native["path"]).resolve(strict=True)
            extension.relative_to(prefix)
            if digest(extension) != native["sha256"]:
                raise ValueError("CRISPRrepeat fold-engine identity does not match the toolset")
    return document


def build(lock_path, resume=False):
    lock = load_lock(lock_path)
    toolset_path = RUNTIME / "toolsets" / lock_id(lock)
    toolset = verify_toolset(toolset_path)
    sources = {
        name: source_inventory(ROOT / name) for name in ("backend", "Integration")
    }
    identity = hashlib.sha256(
        json.dumps(
            {"sources": sources, "toolset": digest(toolset_path / "toolset.json")},
            sort_keys=True,
        ).encode()
    ).hexdigest()
    release = RUNTIME / "releases" / identity
    if (release / "release.json").exists():
        verify(release)
        print("Existing application release:", release)
        return release
    release.mkdir(parents=True, exist_ok=resume)
    source_root = release / "sources"
    for name, inventory in sources.items():
        for relative in inventory:
            target = source_root / name / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / name / relative, target)
    wheels = release / "wheels"
    wheels.mkdir(exist_ok=resume)
    for name in ("backend", "Integration"):
        run(
            [
                ROOT / "backend/.venv/bin/python",
                "-m",
                "pip",
                "wheel",
                source_root / name,
                "--no-deps",
                "--wheel-dir",
                wheels,
            ]
        )
    backend_wheel = next(wheels.glob("crispr_web_backend-*.whl"))
    integration_wheel = next(wheels.glob("crispr_tool_integration-*.whl"))
    prefix = release / "environment"
    run([sys.executable, "-m", "venv", prefix])
    run(
        [
            prefix / "bin/python",
            "-m",
            "pip",
            "install",
            "--require-hashes",
            "-r",
            source_root / "backend/requirements-production.lock",
        ]
    )
    run(
        [
            prefix / "bin/python",
            "-m",
            "pip",
            "install",
            "--no-deps",
            backend_wheel,
            integration_wheel,
        ]
    )
    run([prefix / "bin/python", "-m", "pip", "check"])
    atomic_json(release / "environment.json", environment_inventory(prefix))
    artifacts = {}
    for key, name in NAMES.items():
        if name not in toolset["tools"]:
            continue
        entry = toolset["tools"][name]
        artifacts[key] = {
            "name": key,
            "display_name": name,
            "version": entry["version"],
            "kind": "git_source",
            "sha256_scope": "source_inventory_v1",
            "path": entry["source"],
            "sha256": hashlib.sha256(
                json.dumps(
                    entry["source_files"], sort_keys=True, separators=(",", ":")
                ).encode()
            ).hexdigest(),
            "repository": entry["repository"],
            "commit": entry["commit"],
        }
    artifacts["integration"] = {
        "name": "integration",
        "display_name": "CRISPRloci integration",
        "version": wheel_version(integration_wheel, "crispr-tool-integration"),
        "kind": "wheel",
        "sha256_scope": "file_bytes",
        "path": str(integration_wheel),
        "sha256": digest(integration_wheel),
    }
    manifest = release / "scientific-manifest.json"
    atomic_json(
        manifest,
        {
            "schema": "crispr-web-scientific-release",
            "schema_version": scientific_manifest_version(artifacts),
            "release_id": identity,
            "scientific_artifacts": artifacts,
            "toolset_manifest_sha256": digest(toolset_path / "toolset.json"),
        },
    )
    record = {
        "schema_version": 2,
        "release_id": identity,
        "toolset": str(toolset_path),
        "toolset_manifest_sha256": digest(toolset_path / "toolset.json"),
        "scientific_manifest_sha256": digest(manifest),
        "backend_wheel_sha256": digest(backend_wheel),
        "integration_wheel_sha256": digest(integration_wheel),
        "sources": sources,
        "environment_sha256": digest(release / "environment.json"),
        "environment": str(prefix),
        "created_at": time.time(),
    }
    atomic_json(release / "release.json", record)
    verify(release)
    print("Application release staged:", release)
    return release


def verify(release):
    release = Path(release).resolve(strict=True)
    release.relative_to(RUNTIME / "releases")
    document = json.loads((release / "release.json").read_text())
    if document["release_id"] != release.name:
        raise ValueError("Application release identity mismatch")
    toolset = Path(document["toolset"])
    if digest(toolset / "toolset.json") != document["toolset_manifest_sha256"]:
        raise ValueError("Scientific toolset manifest changed")
    verify_toolset(toolset)
    if (
        digest(release / "scientific-manifest.json")
        != document["scientific_manifest_sha256"]
    ):
        raise ValueError("Scientific manifest changed")
    for glob, field in [
        ("crispr_web_backend-*.whl", "backend_wheel_sha256"),
        ("crispr_tool_integration-*.whl", "integration_wheel_sha256"),
    ]:
        if digest(next((release / "wheels").glob(glob))) != document[field]:
            raise ValueError("Application wheel changed")
    for name, inventory in document["sources"].items():
        if source_inventory(release / "sources" / name) != inventory:
            raise ValueError("Application source snapshot changed")
    if digest(release / "environment.json") != document[
        "environment_sha256"
    ] or environment_inventory(release / "environment") != json.loads(
        (release / "environment.json").read_text()
    ):
        raise ValueError("Installed application environment changed")
    return document


def environment(
    release, data_root, api_origin="http://127.0.0.1:5173", *, verified=None
):
    release = Path(release).resolve(strict=True)
    document = verify(release) if verified is None else verified
    env = dict(os.environ)
    env.update(
        CRISPR_WEB_DATA_ROOT=str(Path(data_root).resolve()),
        CRISPR_WEB_RUNNER_CONFIG=str(Path(document["toolset"]) / "runners.json"),
        CRISPR_WEB_INTEGRATION_COMMAND=json.dumps(
            [str(release / "environment/bin/crispr-tools")]
        ),
        CRISPR_WEB_CORS_ORIGINS=api_origin,
        CRISPR_WEB_RELEASE_ID=document["release_id"],
        CRISPR_WEB_BACKEND_WHEEL_SHA256=document["backend_wheel_sha256"],
        CRISPR_RELEASE_MANIFEST=str(release / "scientific-manifest.json"),
        CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256=document[
            "scientific_manifest_sha256"
        ],
        CRISPR_WEB_RELEASE=wheel_version(next((release / "wheels").glob("crispr_web_backend-*.whl")), "crispr-web-backend"),
        PYTHONNOUSERSITE="1",
        PYTHONDONTWRITEBYTECODE="1",
    )
    for key in ("PYTHONPATH", "PYTHONHOME"):
        env.pop(key, None)
    return env


def require_gate(release):
    gate_path = release / "gate.json"
    if not gate_path.is_file():
        raise ValueError(
            "Run scripts/gate_release.py for this release before activation"
        )
    gate = json.loads(gate_path.read_text())
    if (
        gate.get("release_sha256") != digest(release / "release.json")
        or gate.get("status") != "passed"
    ):
        raise ValueError("Release gate does not match this build")
    manifest_path = release / "scientific-manifest.json"
    if manifest_path.is_file():
        from runtime.context_gate import CONTROL_CASES, require_coverage
        artifacts = json.loads(manifest_path.read_text())["scientific_artifacts"]
        require_coverage(artifacts, gate)
        if set(artifacts) & set(CONTROL_CASES) and (
            gate.get("gate_script_sha256") != digest(ROOT / "scripts/gate_release.py")
            or gate.get("context_gate_sha256") != digest(ROOT / "scripts/runtime/context_gate.py")
        ):
            raise ValueError("Context release gate scripts changed; rerun the native gate")


def activate(release):
    release = Path(release).resolve(strict=True)
    verify(release)
    require_gate(release)
    RUNTIME.mkdir(exist_ok=True)
    with (RUNTIME / "active.lock").open("a") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError(
                "API/worker are running. Stop them before activation; active jobs must finish or be cancelled."
            ) from error
        current = RUNTIME / "current"
        if current.is_symlink():
            previous = current.resolve(strict=True)
            if previous != release:
                temporary = RUNTIME / ".previous.next"
                temporary.unlink(missing_ok=True)
                temporary.symlink_to(previous)
                os.replace(temporary, RUNTIME / "previous")
        temporary = RUNTIME / ".current.next"
        temporary.unlink(missing_ok=True)
        temporary.symlink_to(release)
        os.replace(temporary, current)
        # Keep the source-reference lock in step with the selected runtime. Git
        # checkouts can then be synchronized explicitly; jobs use the snapshot.
        document = json.loads((release / "release.json").read_text())
        atomic_json(
            ROOT / "tools/tools.lock.json",
            load_lock(Path(document["toolset"]) / "tools.lock.json"),
        )
    print("Selected release:", release.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="action", required=True)
    command = subs.add_parser("build")
    command.add_argument("--lock", type=Path, default=ROOT / "tools/tools.lock.json")
    command.add_argument("--resume", action="store_true", help="Resume an incomplete, unactivated build of these exact sources")
    for action in ("verify", "activate"):
        command = subs.add_parser(action)
        command.add_argument("release", type=Path)
    subs.add_parser("rollback")
    args = parser.parse_args()
    if args.action == "build":
        build(args.lock, args.resume)
    elif args.action == "verify":
        verify(args.release)
        print("Release verified")
    elif args.action == "activate":
        activate(args.release)
    else:
        activate((RUNTIME / "previous").resolve(strict=True))


if __name__ == "__main__":
    main()
