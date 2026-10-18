"""Fetch frozen public source files, verify hashes, and extract valid supplements.

Run from a checkout containing benchmarks/literature/download_manifest.json.
Use --verify-only for an entirely offline integrity check. The separate atlas
streamer handles the 5.27 GB source without keeping a second full copy on disk.
"""

import argparse
import hashlib
import json
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "benchmarks/literature"


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    rows = json.loads((ROOT / "download_manifest.json").read_text())
    for row in rows:
        path = ROOT / row["path"]
        if not path.exists():
            if args.verify_only:
                raise FileNotFoundError(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + ".partial")
            with (
                urllib.request.urlopen(row["url"], timeout=90) as stream,
                temporary.open("wb") as out,
            ):
                while block := stream.read(1024 * 1024):
                    out.write(block)
            if digest(temporary) != row["sha256"]:
                raise ValueError(f"Source changed or download invalid: {row['url']}")
            temporary.replace(path)
        if digest(path) != row["sha256"]:
            raise ValueError(f"Checksum mismatch: {path}")
        if row["kind"] == "supplement_zip" and not args.verify_only:
            directory = path.with_suffix("")
            with zipfile.ZipFile(path) as archive:
                for member in archive.infolist():
                    target = (directory / member.filename).resolve()
                    if not target.is_relative_to(directory.resolve()):
                        raise ValueError(f"Unsafe archive member: {member.filename}")
                archive.extractall(directory)
    print(f"Verified {len(rows)} frozen source files")


if __name__ == "__main__":
    main()
