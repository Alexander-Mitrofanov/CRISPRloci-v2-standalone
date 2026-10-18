#!/usr/bin/env python3
"""Fetch pinned upstream repositories or advance one clean checkout explicitly.

No upstream source is patched. Updating a dependency requires re-running the
integration tests and its documented installation/build steps before production.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("CRISPRidentify-v2", "CasAndra")


def git(*args: str, cwd: Path | None = None) -> str:
    return subprocess.check_output(["git", *args], cwd=cwd, text=True).strip()


def sync(name: str, entry: dict[str, str]) -> None:
    destination = ROOT / "external" / name
    if destination.exists():
        if git("status", "--porcelain", "--untracked-files=no", cwd=destination):
            raise RuntimeError(f"Refusing to alter modified upstream checkout: {destination}")
        current = git("rev-parse", "HEAD", cwd=destination)
        if current == entry["commit"]:
            print(f"{name}: already at {current}")
            return
        raise RuntimeError(
            f"{name} is at {current}, not locked {entry['commit']}; "
            "preserve that checkout elsewhere before syncing"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{name}-", dir=destination.parent) as staging:
        checkout = Path(staging) / "checkout"
        git("clone", "--no-checkout", "--filter=blob:none", entry["url"], str(checkout))
        git("checkout", "--detach", entry["commit"], cwd=checkout)
        checkout.rename(destination)
    print(f"{name}: {git('rev-parse', 'HEAD', cwd=destination)}")


def update(name: str, ref: str, lock: dict[str, dict[str, str]]) -> None:
    if not ref.strip() or ref.startswith("-") or ":" in ref:
        raise ValueError("Use a nonempty Git tag, branch, or commit, not an option or refspec")
    destination = ROOT / "external" / name
    if git("status", "--porcelain", cwd=destination):
        raise RuntimeError(f"Refusing update: {destination} has local changes")
    if git("remote", "get-url", "origin", cwd=destination) != lock[name]["url"]:
        raise RuntimeError("Repository remote does not match the dependency lock")
    old = git("rev-parse", "HEAD", cwd=destination)
    if old != lock[name]["commit"]:
        raise RuntimeError(
            "Current checkout differs from lock; resolve the mismatch before updating"
        )
    branch = git("rev-parse", "--abbrev-ref", "HEAD", cwd=destination)
    git("fetch", "origin", ref, cwd=destination)
    new = git("rev-parse", "--verify", "FETCH_HEAD^{commit}", cwd=destination)
    proposed = {key: dict(value) for key, value in lock.items()}
    proposed[name].update(commit=new, previous_commit=old, requested_ref=ref)
    temporary = None
    checkout_attempted = False
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".upstream-lock-",
            suffix=".json",
            dir=ROOT,
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(proposed, handle, indent=2)
            handle.write("\n")
        checkout_attempted = True
        git("checkout", "--detach", new, cwd=destination)
        temporary.replace(ROOT / "upstream-lock.json")
        temporary = None  # Publication consumed the staging path; no cleanup remains.
    except Exception as update_error:
        if checkout_attempted:
            try:
                if branch == "HEAD":
                    git("checkout", "--detach", old, cwd=destination)
                else:
                    git("checkout", "--no-guess", branch, cwd=destination)
            except Exception as restoration_error:
                raise RuntimeError(
                    f"Dependency update failed: {update_error}; restoration also failed: "
                    f"{restoration_error}. The lock was not published. Restore {destination} "
                    f"to {old} ({'detached HEAD' if branch == 'HEAD' else f'branch {branch}'}) "
                    "before retrying."
                ) from update_error
        raise
    finally:
        if temporary is not None:
            pending_error = sys.exception()
            try:
                temporary.unlink(missing_ok=True)
            except OSError as cleanup_error:
                detail = f"Could not remove staged lock {temporary}: {cleanup_error}"
                if pending_error is not None:
                    pending_error.add_note(detail)
                else:
                    print(detail, file=sys.stderr)
    lock[name] = proposed[name]
    print(f"{name}: {old} -> {new}")
    print(
        "Reinstall/rebuild this upstream using its instructions, "
        "then run the real integration tests."
    )
    print("This command records a candidate dependency update; it does not certify compatibility.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("sync", help="Download the exact versions in upstream-lock.json")
    change = commands.add_parser("update", help="Explicitly update one unmodified upstream")
    change.add_argument("name", choices=NAMES)
    change.add_argument("--ref", required=True, help="Git tag, branch, or commit")
    args = parser.parse_args()
    lock = json.loads((ROOT / "upstream-lock.json").read_text())
    if args.command == "sync":
        for name in NAMES:
            sync(name, lock[name])
    else:
        update(args.name, args.ref, lock)


if __name__ == "__main__":
    main()
