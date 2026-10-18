#!/usr/bin/env python3
"""Read-only inventory of committed CrisprCasDb campaign publications.

This standalone standard-library script can run over SSH (``python3 -B -``).
It does not invoke scientific tools, modify results, or import another project.
Bundle hashes are recorded here and must be verified by the fetching process.
Collection-index and raw-document hashes are verified during inventory.
"""

from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys


PROFILE = "5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540"
ROOTS = [
    ("original", "/data/crispr/results-v232-eightvm-20260910", 0),
    ("recovery", "/data/crispr/recovery-20260920/results", 10),
    ("redistribution", "/data/crispr/recovery-20260920/redistribution-20260920-stopped/results", 20),
    ("rebalance", "/data/crispr/recovery-20260920/rebalance-20260921/results", 30),
    ("tail", "/data/crispr/recovery-20260922-tail/results", 40),
    ("retained", "/data/crispr/recovery-20260922-tail/retained-source", 5),
]
ACCESSION = re.compile(r"GC[AF]_\d+\.\d+\Z")
DIGEST = re.compile(r"[0-9a-f]{64}\Z")
MAX_DOCUMENT = 256 * 1024 * 1024


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_bytes(path: Path) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Expected a regular source file: {path}")
    with path.open("rb") as stream:
        data = stream.read(MAX_DOCUMENT + 1)
    if len(data) > MAX_DOCUMENT:
        raise ValueError(f"Source metadata exceeds the bounded read: {path}")
    return data


def document(path: Path) -> tuple[bytes, dict]:
    data = read_bytes(path)
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as stream:
            decoded = stream.read(MAX_DOCUMENT + 1)
        if len(decoded) > MAX_DOCUMENT:
            raise ValueError(f"Expanded source metadata is too large: {path}")
    else:
        decoded = data
    value = json.loads(decoded)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return data, value


def safe_child(root: Path, relative: str) -> Path:
    if not isinstance(relative, str):
        raise ValueError("Source member path must be text")
    pure = PurePosixPath(relative)
    if pure.is_absolute() or ".." in pure.parts or "\\" in relative:
        raise ValueError(f"Unsafe source member path: {relative!r}")
    path = root.joinpath(*pure.parts)
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Source member escapes its root: {relative!r}")
    return path


def checked_entry(entry: dict) -> dict:
    accession, status = entry.get("accession"), entry.get("status")
    if not isinstance(accession, str) or not ACCESSION.fullmatch(accession):
        raise ValueError(f"Invalid accession: {accession!r}")
    if status not in {"complete", "failed", "unavailable"}:
        raise ValueError(f"Nonterminal source record: {accession}: {status!r}")
    sha = entry.get("sha256")
    if not isinstance(sha, str) or not DIGEST.fullmatch(sha):
        raise ValueError(f"Missing source checksum: {accession}")
    result = {"accession": accession, "status": status, "sha256": sha}
    if "summary" in entry:
        result["summary"] = entry["summary"]
    if "error" in entry:
        result["error"] = entry["error"]
    return result


def profile_check(value: dict, path: Path) -> None:
    actual = value.get("analysis_profile_sha256")
    if actual is None:
        actual = value.get("provenance", {}).get("analysis_profile_sha256")
    if actual != PROFILE:
        raise ValueError(f"Unexpected scientific profile at {path}: {actual}")


def inventory(roots: list[tuple[str, str, int]], node: str | None = None) -> dict:
    sources = []
    indexed = set()
    seen_paths = set()
    for name, directory, precedence in roots:
        root = Path(directory)
        if not root.exists():
            continue
        if root.is_symlink():
            raise ValueError(f"Source root must not be a symlink: {root}")
        for index_path in sorted(root.rglob("public/index.json")):
            _, public_index = document(index_path)
            for publication in public_index["publications"]:
                profile_check(publication, index_path)
                shard_name = publication["shard_name"]
                if not re.fullmatch(r"shard-\d+", shard_name):
                    raise ValueError(f"Invalid shard name: {shard_name!r}")
                collection_path = index_path.parent.parent / "shards" / shard_name / "collection_index.json"
                raw, collection = document(collection_path)
                if digest(raw) != publication["collection_index_sha256"]:
                    raise ValueError(f"Collection-index checksum mismatch: {collection_path}")
                profile_check(collection, collection_path)
                bundle = publication["bundle"]
                if not DIGEST.fullmatch(bundle["sha256"]):
                    raise ValueError(f"Invalid bundle checksum: {index_path}")
                path = safe_child(index_path.parent, bundle["download_path"])
                if path.is_symlink() or not path.is_file() or path.stat().st_size != bundle["bytes"]:
                    raise ValueError(f"Missing or wrong-sized bundle: {path}")
                entries = [checked_entry(entry) for entry in collection["entries"]]
                indexed.update(str(safe_child(collection_path.parent, entry["file"])) for entry in collection["entries"])
                if str(path) in seen_paths:
                    continue
                seen_paths.add(str(path))
                sources.append({
                    "kind": "bundle", "path": str(path), "sha256": bundle["sha256"],
                    "bytes": bundle["bytes"], "precedence": precedence,
                    "origin": name, "node": node, "entries": entries,
                    "profile_sha256": PROFILE,
                    "collection_index_sha256": digest(raw),
                    "generated_at": collection.get("generated_at"),
                })
        # Unfinished batches and the retained VM3 subset can contain committed
        # terminal documents without an independently published catalogue shard.
        for pattern in ("assemblies/*.json.gz", "failures/*.json.gz"):
            for path in sorted(root.rglob(pattern)):
                if str(path) in indexed or str(path) in seen_paths:
                    continue
                raw, value = document(path)
                profile_check(value, path)
                entry = checked_entry({
                    "accession": value["assembly"]["accession"],
                    "status": value["status"], "sha256": digest(raw),
                    **({"summary": value["summary"]} if "summary" in value else {}),
                    **({"error": value["error"]} if "error" in value else {}),
                })
                sources.append({
                    "kind": "raw", "path": str(path), "sha256": digest(raw),
                    "bytes": len(raw), "precedence": precedence, "origin": name,
                    "node": node, "entries": [entry], "profile_sha256": PROFILE,
                    "generated_at": value.get("finished_at", value.get("provenance", {}).get("finished_at")),
                })
                seen_paths.add(str(path))
        # NCBI omissions have no scientific result bundle. Their committed job
        # receipt is evidence of unavailability, never evidence of absence.
        for path in sorted(root.rglob("job-receipt.json")):
            raw, receipt = document(path)
            unavailable = receipt.get("unavailable", {})
            if not unavailable:
                continue
            if not isinstance(unavailable, dict):
                raise ValueError(f"Invalid unavailable outcomes: {path}")
            entries = [checked_entry({
                "accession": accession, "status": "unavailable", "sha256": digest(raw),
                "error": {"type": "NCBIUnavailable", "message": outcome.get("reason")},
            }) for accession, outcome in sorted(unavailable.items())]
            sources.append({
                "kind": "unavailable", "path": str(path), "sha256": digest(raw),
                "bytes": len(raw), "precedence": precedence, "origin": name,
                "node": node, "entries": entries, "generated_at": receipt.get("finished_at"),
            })
    counts = Counter(entry["status"] for source in sources for entry in source["entries"])
    return {
        "schema": "crisprloci.source-inventory.v1", "node": node,
        "profile_sha256": PROFILE, "created_at": datetime.now(timezone.utc).isoformat(),
        "roots": [directory for _, directory, _ in roots], "sources": sources,
        "candidate_counts": dict(counts),
        "bundle_hash_verification": "required_at_fetch",
    }


def merge_inventories(inventories: list[dict]) -> dict:
    """Select one authoritative outcome per accession across all machines."""
    sources = []
    selected = {}
    for value in inventories:
        if value.get("profile_sha256") != PROFILE:
            raise ValueError("Cannot combine different scientific profiles")
        for source in value["sources"]:
            source_index = len(sources)
            sources.append({**source, "entries": []})
            for entry in source["entries"]:
                entry = checked_entry(entry)
                accession = entry["accession"]
                previous = selected.get(accession)
                priority = (entry["status"] == "complete", source["precedence"], source.get("generated_at") or "")
                if previous is not None:
                    old_entry, _, old_priority = previous
                    if entry["status"] == old_entry["status"] == "complete" and entry["sha256"] != old_entry["sha256"]:
                        raise ValueError(f"Conflicting successful results: {accession}")
                    if priority <= old_priority:
                        continue
                selected[accession] = (entry, source_index, priority)
    for entry, source_index, _ in selected.values():
        sources[source_index]["entries"].append(entry)
    counts = Counter(entry["status"] for entry, _, _ in selected.values())
    return {
        "schema": "crisprloci.source-plan.v1", "profile_sha256": PROFILE,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "counts": dict(counts), "total": len(selected),
        "sources": [source for source in sources if source["entries"]],
        "selection": "successful outcome first; then recovery precedence and recorded time",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node")
    parser.add_argument("--root", action="append", help="Explicit source root; later roots have higher precedence")
    parser.add_argument("--merge", nargs="+", metavar="INVENTORY", help="Merge saved inventories and select one outcome per accession")
    parser.add_argument("--expect-total", type=int)
    parser.add_argument("--gzip", action="store_true")
    args = parser.parse_args()
    if args.merge:
        result = merge_inventories([document(Path(path))[1] for path in args.merge])
        if args.expect_total is not None and result["total"] != args.expect_total:
            raise ValueError(f"Coverage mismatch: {result['total']} != {args.expect_total}")
    else:
        roots = [(f"source-{index}", path, index * 10) for index, path in enumerate(args.root)] if args.root else ROOTS
        result = inventory(roots, args.node)
    data = json.dumps(result, separators=(",", ":"), sort_keys=True).encode("utf-8") + b"\n"
    sys.stdout.buffer.write(gzip.compress(data, compresslevel=3, mtime=0) if args.gzip else data)


if __name__ == "__main__":
    main()
