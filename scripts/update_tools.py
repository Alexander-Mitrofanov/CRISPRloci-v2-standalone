#!/usr/bin/env python3
"""Resolve GitHub updates into a candidate lock; never modifies a running release."""

from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import re
import subprocess
import tomllib

from runtime.common import ROOT, RUNTIME, atomic_json, load_lock, lock_id, run


def candidate(selections, output=None):
    lock = copy.deepcopy(load_lock(ROOT / "tools/tools.lock.json"))
    entries = {entry["name"]: entry for entry in lock["tools"]}
    for selection in selections:
        name, sep, ref = selection.partition("@")
        if (
            not sep
            or name not in entries
            or not ref
            or ref.startswith("-")
            or not re.fullmatch(r"[A-Za-z0-9_./-]+", ref)
        ):
            raise ValueError("Use --tool NAME@branch-or-tag-or-commit")
        entry = entries[name]
        repository = ROOT / entry["path"]
        configured = subprocess.check_output(
            ["git", "-C", str(repository), "config", "--get", "remote.origin.url"],
            text=True,
        ).strip()
        if configured != entry["repository"]:
            raise ValueError("Repository origin differs from lock")
        run(["git", "-C", repository, "fetch", "--no-tags", "origin", ref])
        commit = subprocess.check_output(
            ["git", "-C", str(repository), "rev-parse", "FETCH_HEAD^{commit}"],
            text=True,
        ).strip()
        entry["commit"] = commit
        result = subprocess.run(
            ["git", "-C", str(repository), "show", commit + ":pyproject.toml"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            metadata = tomllib.loads(result.stdout)
            entry["version"] = metadata.get("project", {}).get(
                "version", entry["version"]
            )
        elif name == "CRISPRidentify-v2":
            result = subprocess.run(
                ["git", "-C", str(repository), "show", commit + ":VERSION"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                version = result.stdout.strip()
                if not re.fullmatch(r"\d+\.\d+\.\d+(?:[A-Za-z0-9.+-]*)", version):
                    raise ValueError("CRISPRidentify VERSION is invalid")
                entry["version"] = version
                continue
            result = subprocess.run(
                ["git", "-C", str(repository), "show", commit + ":README.md"],
                capture_output=True,
                text=True,
                check=True,
            )
            match = re.match(r"# CRISPRidentify v(\d+\.\d+\.\d+\S*)", result.stdout)
            if match:
                entry["version"] = match.group(1)
    tracr = entries["CRISPRtracrRNA-v3"]
    expected = json.loads(
        subprocess.check_output(
            [
                "git",
                "-C",
                str(ROOT / tracr["path"]),
                "show",
                tracr["commit"] + ":upstream-lock.json",
            ],
            text=True,
        )
    )
    for name in ("CRISPRidentify-v2", "CasAndra"):
        if entries[name]["commit"] != expected[name]["commit"]:
            raise ValueError(
                f"{name} and CRISPRtracrRNA dependency pins disagree; propose their compatible revisions together"
            )
    destination = output or RUNTIME / "candidates" / f"{lock_id(lock)}.json"
    if destination.resolve() == (ROOT / "tools/tools.lock.json").resolve():
        raise ValueError("Candidate output must not overwrite the active source lock")
    atomic_json(destination, lock)
    print("Candidate lock:", destination)
    print(
        "Next: provision_tools.py --lock PATH; manage_release.py build --lock PATH; gate_release.py RELEASE; manage_release.py activate RELEASE"
    )
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tool",
        action="append",
        required=True,
        help="Tool@ref; repeat for a compatible update set",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    candidate(args.tool, args.output)
