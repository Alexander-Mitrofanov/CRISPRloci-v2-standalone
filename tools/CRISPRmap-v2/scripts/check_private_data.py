#!/usr/bin/env python3
"""Reject private operations material before commits, pushes and frontend export.

Never print matched content. Git modes inspect stored blobs, not working copies.
Workspace mode works even in source archives without Git metadata.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MAX_BYTES = 8 * 1024 * 1024
PRIVATE_PARTS = {".ops-private", ".ssh", ".runtime", ".runtime-data", ".catalog", ".venv", ".envs", "runs"}
GENERATED = {".git", "node_modules", "dist", "build", "coverage", "__pycache__"}
PATTERNS = {
    "private key": re.compile(rb"-----BEGIN (?:OPENSSH |RSA |EC |DSA |ENCRYPTED )?PRIVATE KEY-----"),
    "SSH connection configuration": re.compile(rb"(?m)^\s*(?:HostName|ProxyJump|IdentityFile)\s+\S+"),
    "SSH connection command": re.compile(rb"(?im)^\s*(?:ssh|scp|sftp)\s+[^\n]*\b[a-z_][\w-]*@[a-z0-9][\w.-]*"),
    "GitHub credential": re.compile(rb"\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,})\b"),
}


def git(root, *args):
    result = subprocess.run(["git", "-C", str(root), *args], capture_output=True)
    if result.returncode:
        raise ValueError("Git inspection failed; run from a real checkout with accessible metadata")
    return result.stdout


def private_path(name):
    path = Path(name)
    return (
        bool(set(path.parts) & PRIVATE_PARTS)
        or path.parts[0] in {"output", "spacer analisis"}
        or (path.name.startswith(".env") and not path.name.endswith(".example"))
        or path.suffix.lower() in {".pem", ".key", ".p12", ".pfx", ".sqlite3"}
        or path.name in {"ssh-config", "authorized_keys"}
        or path.name.startswith(("id_rsa", "id_ed25519", "known_hosts"))
    )


def literals(root):
    path = root / ".ops-private/forbidden-values.json"
    if not path.exists():
        return []
    if path.is_symlink() or path.stat().st_size > MAX_BYTES:
        raise ValueError("Unsafe private signature file")
    values = json.loads(path.read_text())
    if not isinstance(values, list) or any(not isinstance(v, str) or len(v) < 5 for v in values):
        raise ValueError("Private signatures must be a JSON array of strings of at least five characters")
    return [value.lower().encode() for value in values]


def inspect(name, payload, signatures):
    reasons = []
    if private_path(name):
        reasons.append("private operations path")
    for label, pattern in PATTERNS.items():
        if pattern.search(payload):
            reasons.append(label)
    if any(value in payload.lower() or value in name.lower().encode() for value in signatures):
        reasons.append("registered private connection detail")
    return reasons


def workspace_files(root):
    # A disposable Git directory applies the actual .gitignore without changing
    # the checkout or relying on its index. Include tracked files when available:
    # adding an ignore rule does not untrack previously committed data.
    with tempfile.TemporaryDirectory(prefix="crispr-privacy-") as directory:
        subprocess.run(["git", "init", "--bare", "--quiet", directory], check=True)
        result = subprocess.run([
            "git", "--git-dir", directory, "--work-tree", str(root),
            "ls-files", "--others", "--exclude-standard", "-z",
        ], check=True, capture_output=True)
    paths = set(result.stdout.split(b"\0")) - {b""}
    try:
        paths.update(git(root, "ls-files", "-z").split(b"\0"))
    except ValueError:
        pass
    return sorted(os.fsdecode(path) for path in paths if path)


def file_entries(root, paths):
    for name in paths:
        path = root / name
        if path.is_symlink():
            yield name, b"", "symlink cannot be published"
        elif not path.is_file():
            yield name, b"", "not a regular file"
        elif path.stat().st_size > MAX_BYTES:
            yield name, b"", "file exceeds publication inspection limit"
        else:
            yield name, path.read_bytes(), None


def tree_entries(root, revision=None):
    rows = git(root, "ls-files", "--stage", "-z") if revision is None else git(root, "ls-tree", "-r", "-z", revision)
    for row in rows.split(b"\0"):
        if not row:
            continue
        meta, raw_name = row.split(b"\t", 1)
        fields = meta.split()
        mode, oid = fields[0], fields[1] if revision is None else fields[2]
        name = os.fsdecode(raw_name)
        if mode not in {b"100644", b"100755"}:
            yield name, b"", "non-file Git entry cannot be published"
        elif int(git(root, "cat-file", "-s", oid.decode())) > MAX_BYTES:
            yield name, b"", "blob exceeds publication inspection limit"
        else:
            yield name, git(root, "cat-file", "blob", oid.decode()), None


def pushed_commits(root, stream):
    commits = set()
    for line in stream:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError("Invalid pre-push input")
        _, local, _, remote = fields
        if not all(re.fullmatch(r"[0-9a-f]{40,64}", value) for value in (local, remote)):
            raise ValueError("Invalid pre-push object")
        if set(local) == {"0"}:
            continue
        # Tags must resolve to commits; a tag pointing straight at a private
        # blob/tree must not evade commit-tree inspection.
        git(root, "rev-parse", "--verify", f"{local}^{{commit}}")
        # Inspect all reachable history for new refs, including commits deleted
        # from the tip. Existing refs inspect every newly introduced commit.
        revisions = [local] if set(remote) == {"0"} else [f"{remote}..{local}"]
        commits.update(git(root, "rev-list", *revisions).decode().splitlines())
    return sorted(commits)


def check(root, mode, export=None, stream=None):
    signatures = literals(root)
    findings = set()
    count = 0
    if mode == "workspace":
        batches = [file_entries(root, workspace_files(root))]
    elif mode == "index":
        batches = [tree_entries(root)]
    elif mode == "pre-push":
        batches = (tree_entries(root, commit) for commit in pushed_commits(root, stream))
    else:
        paths = []
        for directory, names, files in os.walk(export, followlinks=False):
            for name in list(names):
                path = Path(directory) / name
                if path.is_symlink():
                    paths.append(str(path.relative_to(export)))
                if name in GENERATED or path.is_symlink():
                    names.remove(name)
            for name in files:
                # Match only the exporter's explicit source exclusions. The
                # destination scan remains strict, including ignored files.
                if mode == "export-source" and (
                    (name.startswith(".env") and name != ".env.example")
                    or name.endswith((".local", ".log", ".pem", ".key"))
                    or name == ".DS_Store"
                ):
                    continue
                paths.append(str((Path(directory) / name).relative_to(export)))
        batches = [file_entries(export, paths)]
    for entries in batches:
        for name, payload, error in entries:
            count += 1
            for reason in ([error] if error else []) + inspect(name, payload, signatures):
                findings.add((name, reason))
    for name, reason in sorted(findings):
        safe_name = name
        for value in signatures:
            safe_name = re.sub(re.escape(value.decode()), "[private identifier]", safe_name, flags=re.IGNORECASE)
        print(f"BLOCKED {json.dumps(safe_name)}: {reason}", file=sys.stderr)
    print(f"Private-data check: {count} file versions inspected; {len(findings)} finding(s).")
    return not findings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("workspace", "index", "pre-push", "export", "export-source"))
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--export", type=Path)
    args = parser.parse_args()
    if args.mode.startswith("export") and (args.export is None or not args.export.is_dir() or args.export.is_symlink()):
        parser.error("--export must be a normal directory")
    try:
        passed = check(args.root.resolve(), args.mode, args.export, sys.stdin)
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print(f"Private-data check failed: {error}", file=sys.stderr)
        return 1
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
