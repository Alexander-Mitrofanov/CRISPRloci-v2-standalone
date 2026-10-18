#!/usr/bin/env python3
"""Install local publication guards without overwriting another hooks policy."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--allow-source-archive", action="store_true")
    args = parser.parse_args()
    check = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--show-toplevel"], capture_output=True, text=True)
    if check.returncode or Path(check.stdout.strip()).resolve() != ROOT:
        if args.allow_source_archive:
            print("Source archive: Git guards NOT installed. Install them in the real checkout before publishing.", file=sys.stderr)
            return
        raise SystemExit("No accessible Git metadata. Install guards in the real checkout before committing or pushing.")
    setting = subprocess.run(["git", "-C", str(ROOT), "config", "--get", "core.hooksPath"], capture_output=True, text=True)
    if setting.returncode not in (0, 1):
        raise SystemExit("Cannot inspect existing Git hook configuration")
    if setting.stdout.strip() not in ("", ".githooks", str(ROOT / ".githooks")):
        raise SystemExit("Existing hooksPath preserved. Integrate the privacy hooks into that policy before publishing.")
    if not setting.stdout.strip():
        hooks = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--git-path", "hooks"], check=True, capture_output=True, text=True)
        directory = Path(hooks.stdout.strip())
        if not directory.is_absolute():
            directory = ROOT / directory
        if directory.is_dir() and any(
            path.is_file() and not path.name.endswith(".sample") and os.access(path, os.X_OK)
            for path in directory.iterdir()
        ):
            raise SystemExit("Existing default Git hooks preserved. Integrate the privacy hooks before publishing.")
    for name in ("pre-commit", "pre-push"):
        path = ROOT / ".githooks" / name
        path.chmod(0o755)
    subprocess.run(["git", "-C", str(ROOT), "config", "--local", "core.hooksPath", ".githooks"], check=True)
    print("Installed staged-blob and outgoing-history privacy guards.")


if __name__ == "__main__":
    main()
