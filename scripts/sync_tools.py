#!/usr/bin/env python3
"""Restore/check independent tool repositories at tools/tools.lock.json revisions."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

from bundled_tools import verify_bundled

ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "tools" / "tools.lock.json"


def git(path, *args):
    return subprocess.check_output(
        ["git", "-C", str(path), *args], text=True, stderr=subprocess.PIPE
    ).strip()


def validate_entry(entry):
    name = entry["name"]
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name):
        raise ValueError(f"Invalid tool name: {name!r}")
    if entry["path"] != f"tools/{name}":
        raise ValueError(f"Invalid checkout path for {name}")
    if not re.fullmatch(r"[0-9a-f]{40}", entry["commit"]):
        raise ValueError(f"Missing or invalid commit for {name}")
    if not re.fullmatch(
        r"https://github\.com/[A-Za-z0-9_-]+/[A-Za-z0-9_.-]+\.git", entry["repository"]
    ):
        raise ValueError(f"Invalid repository URL for {name}")


def sync(entry, check):
    name, url, revision = entry["name"], entry["repository"], entry["commit"]
    if verify_bundled(ROOT, entry):
        print(f"{name}: OK {revision} (bundled source)")
        return
    path = ROOT / entry["path"]
    if path.is_symlink():
        raise ValueError(f"{name}: checkout path must not be a symlink")
    if not path.exists():
        if check:
            raise ValueError(f"{name}: missing; run python3 scripts/sync_tools.py")
        path.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["git", "clone", "--no-checkout", "--", url, str(path)], check=True
        )
        # Clone intentionally leaves an empty working tree; populate the pinned
        # commit directly before enforcing cleanliness on subsequent invocations.
        try:
            git(path, "cat-file", "-e", revision + "^{commit}")
        except subprocess.CalledProcessError:
            subprocess.run(
                ["git", "-C", str(path), "fetch", "origin", revision], check=True
            )
        subprocess.run(
            ["git", "-C", str(path), "checkout", "--detach", revision], check=True
        )
    if git(path, "rev-parse", "--show-toplevel") != str(path.resolve()):
        raise ValueError(f"{name}: path is not an independent repository")
    # Compare the configured repository identity, before any user transport
    # rewrite (for example, HTTPS to SSH via url.*.insteadOf).
    actual_url = git(path, "config", "--get", "remote.origin.url")
    if actual_url != url:
        raise ValueError(f"{name}: origin differs from lock: {actual_url}")
    if git(path, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError(
            f"{name}: has local changes; preserve/commit them before syncing"
        )
    actual_revision = git(path, "rev-parse", "HEAD")
    if actual_revision != revision:
        if check:
            raise ValueError(f"{name}: expected {revision}, found {actual_revision}")
        try:
            git(path, "cat-file", "-e", revision + "^{commit}")
        except subprocess.CalledProcessError:
            subprocess.run(
                ["git", "-C", str(path), "fetch", "origin", revision], check=True
            )
        subprocess.run(
            ["git", "-C", str(path), "checkout", "--detach", revision], check=True
        )
    if git(path, "rev-parse", "HEAD") != revision:
        raise ValueError(f"{name}: checkout did not resolve to locked revision")
    print(f"{name}: OK {revision}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Verify local checkouts without network access or changes",
    )
    parser.add_argument(
        "--tool",
        action="append",
        help="Only process this tool; repeat to select several",
    )
    args = parser.parse_args()
    lock = json.loads(LOCK.read_text())
    if lock["schema_version"] != 1:
        raise ValueError("Unsupported tools lock schema")
    entries = lock["tools"]
    names = [entry["name"] for entry in entries]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate tool names in lock")
    for entry in entries:
        validate_entry(entry)
    if args.tool and set(args.tool) - set(names):
        parser.error("Unknown tool: " + ", ".join(sorted(set(args.tool) - set(names))))
    failed = False
    for entry in entries:
        if args.tool and entry["name"] not in args.tool:
            continue
        try:
            sync(entry, args.check)
        except (ValueError, subprocess.CalledProcessError) as error:
            failed = True
            print(f"ERROR: {error}", file=sys.stderr)
            if isinstance(error, subprocess.CalledProcessError) and error.stderr:
                print(error.stderr.strip(), file=sys.stderr)
    return int(failed)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"ERROR: {error}")
