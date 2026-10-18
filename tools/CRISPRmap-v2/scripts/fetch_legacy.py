#!/usr/bin/env python3
"""Explicitly download and verify the frozen original CRISPRmap 2013 assets.

The repository does not redistribute these third-party assets. Analysis never
calls this downloader. All URLs and content digests are pinned below.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

NAMESPACE = "crisprmap-2013-v1"
BASE = "https://rna.informatik.uni-freiburg.de"
ASSETS = {
    "reference.fa": {
        "url": BASE + "/CRISPRmap/oridata/CRISPRmap_data.fa",
        "sha256": "87d0c59372648ef1a2fbb8dd0bc06e887dbf1fbc5d4615e3e02eb488801c972e",
        "bytes": 305908,
    },
    "original-input.fa": {
        "url": BASE + "/DownloadFile.jsp?jobID=8888684&toolName=CRISPRmap&fileNameOrType=input.fa",
        "sha256": "c6c99cba8b451212db420604104af566f31d429d49a294482d8e8f1fc09e2cb8",
        "bytes": 3781,
    },
    "original-results.tsv": {
        "url": BASE + "/DownloadFile.jsp?jobID=8888684&toolName=CRISPRmap&fileNameOrType=Results.tab",
        "sha256": "26ae52c701eacc5c32a5ce123eba94951fdacfb3274d4517b3c75fcfa5efbe25",
        "bytes": 1407,
    },
}


def checked_bytes(data: bytes, spec: dict) -> bytes:
    if len(data) != spec["bytes"] or hashlib.sha256(data).hexdigest() != spec["sha256"]:
        raise ValueError("Asset differs from pinned content; refusing to update the historical reference")
    return data


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    name = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".crisprmap-", delete=False) as handle:
            name = handle.name
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if name and os.path.exists(name):
            os.unlink(name)


def metadata_from_legacy(records: list[dict]) -> dict[str, dict]:
    """Parse labels from historical identifiers; X remains missing."""
    metadata = {}
    for record in records:
        match = re.search(r"_(M\d+|X)_(F\d+|X)$", record["id"])
        if match is None:
            raise ValueError(f"Unrecognized original reference identifier: {record['id']}")
        motif, family = match.groups()
        description = record.get("description", "")
        accession = re.search(r"(?:^|\s)((?:NC|NZ)_\d+(?:\.\d+)?|[A-Z]{1,4}\d+(?:\.\d+)?)_", description)
        metadata[record["id"]] = {
            "family": None if family == "X" else family,
            "motif": None if motif == "X" else motif,
            "accession": accession.group(1) if accession else None,
        }
    return metadata


def parse_original_results(text: str) -> list[dict]:
    """Original Results.tab is motif, family, superclass, description, sequence."""
    rows = []
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != 5:
            raise ValueError(f"Original result line {line_number} has {len(fields)} columns, expected 5")
        motif, family, superclass, description, sequence = fields
        match = re.fullmatch(r"ID(\d+) UserId=(\S+)(?: CRISPRmap=(\S+))?", description)
        if not match:
            raise ValueError(f"Original result line {line_number} has an invalid identifier")
        if not re.fullmatch(r"[ACGUTacgut]+", sequence):
            raise ValueError(f"Original result line {line_number} has an invalid sequence")
        if any(value != "-" and not value.isdigit() for value in [motif, family]):
            raise ValueError(f"Original result line {line_number} has an invalid class")
        rows.append({
            "ordinal": int(match.group(1)), "original_user_id": match.group(2),
            "original_reference_id": match.group(3), "sequence": sequence.upper().replace("U", "T"),
            "family": None if family == "-" else "F" + family,
            "motif": None if motif == "-" else "M" + motif,
            "superclass": None if superclass == "-" else superclass,
        })
    if [r["ordinal"] for r in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("Original result ordinals are missing or out of order")
    return rows


def fetch_assets(destination: Path, cache: Path | None = None) -> dict:
    verified = {}
    for name, spec in ASSETS.items():
        existing = cache / name if cache else destination / name
        if existing.is_file():
            data = checked_bytes(existing.read_bytes(), spec)
        elif cache:
            raise ValueError(f"Missing pinned cache asset: {existing}")
        else:
            request = urllib.request.Request(spec["url"], headers={"User-Agent": "CRISPRmap-v2-reference-import/2.0"})
            with urllib.request.urlopen(request, timeout=45) as response:
                if not response.geturl().startswith("https://rna.informatik.uni-freiburg.de/"):
                    raise ValueError("Unexpected download redirect")
                data = checked_bytes(response.read(spec["bytes"] + 1), spec)
        verified[name] = data
    # Verify every download before exposing any of them.
    for name, data in verified.items():
        atomic_write(destination / name, data)
    return verified


def build_legacy_reference(data_dir: Path) -> dict:
    from crisprmap.io import read_fasta
    from crisprmap.reference import build_reference
    raw = checked_bytes((data_dir / "reference.fa").read_bytes(), ASSETS["reference.fa"])
    records, digest = read_fasta(data_dir / "reference.fa")
    if hashlib.sha256(raw).hexdigest() != digest:
        raise ValueError("Reference changed while importing")
    metadata = metadata_from_legacy(records)
    provenance = {
        "namespace": NAMESPACE, "source_url": ASSETS["reference.fa"]["url"],
        "source_sha256": digest, "source_release": "CRISPRmap v1.0-2013 reference data",
        "retrieval_verified_date": "2026-09-30",
        "citation": "Lange et al. 2013; doi:10.1093/nar/gkt606",
        "license_status": "Third-party data license not separately specified; downloaded by operator, not redistributed",
    }
    return build_reference(records, metadata, provenance)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data" / "legacy" / "2013")
    parser.add_argument("--cache-dir", type=Path, help="Read the three pinned filenames locally instead of using network")
    args = parser.parse_args(argv)
    try:
        from crisprmap.io import read_fasta
        assets = fetch_assets(args.output_dir, args.cache_dir)
        records, digest = read_fasta(args.output_dir / "reference.fa")
        if len(records) != 3527:
            raise ValueError("Original reference must contain exactly 3527 records")
        reference = build_legacy_reference(args.output_dir)
        original = parse_original_results(assets["original-results.tsv"].decode("utf-8"))
        if len(original) != 16:
            raise ValueError("Original example must contain exactly 16 records")
        manifest = {
            "namespace": NAMESPACE, "source_assets": ASSETS,
            "original_example_job_id": "8888684", "original_version": "v1.0-2013",
            "original_runtime_seconds": None,
            "baseline_kind": "published_original_output_replay",
            "original_execution_status": "not_rerun_original_service_offline_and_local_source_unavailable",
        }
        for name, doc in [("reference.json", reference), ("original-results.json", original), ("manifest.json", manifest)]:
            atomic_write(args.output_dir / name, (json.dumps(doc, indent=2, sort_keys=True) + "\n").encode())
        print(json.dumps({"reference": str(args.output_dir / "reference.json"), "records": len(records),
                          "original_example_queries": len(original), "namespace": NAMESPACE}, indent=2))
        return 0
    except (OSError, ValueError) as error:
        print(f"Import failed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
