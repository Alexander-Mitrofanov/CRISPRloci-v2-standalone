#!/usr/bin/env python3
"""Execute fixed installer stages under the managed Python 3.12 environment."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess

from installer.engine import atomic_json
from installer.plan import ROOT, choices
from installer.references import build_reference


def selected_release(work):
    release = Path(json.loads((work / "release-path.json").read_text())["release"]).resolve(strict=True)
    release.relative_to(ROOT / ".runtime/releases")
    return release


def execute(stage, work, value):
    lock = ROOT / "tools/tools.lock.json"
    if stage == "application":
        subprocess.run(["bash", str(ROOT / "scripts/setup.sh"), "--development" if value["developer"] else "--runtime-only"], check=True)
    elif stage == "tools":
        from provision_tools import build

        build(lock, resume=True)
    elif stage == "release":
        from manage_release import build

        release = build(lock, resume=True)
        atomic_json(work / "release-path.json", {"release": str(release)})
    elif stage.startswith("reference-"):
        key = stage.removeprefix("reference-")
        build_reference(ROOT, work, key, value["references"][key], selected_release(work))
    elif stage == "gate":
        from gate_release import gate
        from manage_release import require_gate

        release = selected_release(work)
        try:
            require_gate(release)
            print("Existing scientific gate matches this exact release.", flush=True)
        except (ValueError, OSError):
            gate(release)
    elif stage == "tests":
        subprocess.run(["bash", str(ROOT / "scripts/check.sh")], check=True)
    elif stage == "verify-local":
        from installer.verify import verify_local

        verify_local(ROOT, work, selected_release(work), value)
    elif stage == "activate":
        from manage_release import activate, verify

        release = selected_release(work)
        verify(release)
        environment = {}
        for key in value["references"]:
            receipt = json.loads((work / (key + ".reference.json")).read_text())
            environment.update(receipt["environment"])
        activate(release)
        atomic_json(ROOT / ".installer/installation.json", {
            "schema_version": 1, "project": str(ROOT), "release": str(release),
            "environment": environment, "choices": value,
            "data_root": str(ROOT / ".runtime-data/standalone"),
        })
        print("Installation validated and ready. Run bash start.sh.", flush=True)
    else:
        raise ValueError("Unsupported installer stage")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    execute(args.stage, args.work, choices(json.loads(args.config.read_text())))
