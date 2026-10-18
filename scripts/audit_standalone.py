#!/usr/bin/env python3
"""Detect executable/configuration references and symlinks to retired source trees."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def audit(retired):
    retired = [str(path.resolve()) for path in retired]
    findings = []
    scanned = 0
    code_roots = [
        "scripts",
        "backend/src",
        "Integration/crispr_integration",
        "frontend/src",
        "deploy",
    ]
    extensions = {
        ".py",
        ".sh",
        ".json",
        ".toml",
        ".js",
        ".vue",
        ".service",
        ".timer",
        ".yml",
        ".template",
    }
    for relative in code_roots:
        for directory, names, files in os.walk(ROOT / relative):
            names[:] = [name for name in names if name != "__pycache__"]
            for name in files:
                path = Path(directory) / name
                if path.suffix not in extensions:
                    continue
                scanned += 1
                payload = path.read_text()
                if any(value in payload for value in retired):
                    findings.append(
                        {"kind": "code_reference", "path": str(path.relative_to(ROOT))}
                    )
    symlinks = 0
    for relative in [".runtime", "backend/.venv", "tools"]:
        for directory, names, files in os.walk(ROOT / relative):
            names[:] = [name for name in names if name not in {".git", "__pycache__"}]
            for name in [*names, *files]:
                path = Path(directory) / name
                if not path.is_symlink():
                    continue
                symlinks += 1
                target = str(path.resolve())
                if any(
                    target == value or target.startswith(value + os.sep)
                    for value in retired
                ):
                    findings.append(
                        {
                            "kind": "retired_symlink",
                            "path": str(path.relative_to(ROOT)),
                            "target": target,
                        }
                    )
    result = {
        "status": "passed" if not findings else "failed",
        "retired_roots": retired,
        "scanned_code_files": scanned,
        "checked_symlinks": symlinks,
        "findings": findings,
        "scope": "Runtime code/configuration and runtime/tool symlinks. Historical provenance documents are intentionally excluded.",
    }
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retired", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = audit(args.retired)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    raise SystemExit(0 if result["status"] == "passed" else 1)
