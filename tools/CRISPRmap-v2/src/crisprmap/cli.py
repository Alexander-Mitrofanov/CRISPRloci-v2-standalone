"""Explicit reference construction and offline repeat analysis."""

from __future__ import annotations

import argparse
import csv
from io import StringIO
import json
import math
from pathlib import Path
import sys

from . import __version__
from .core import search
from .exports import write_report
from .io import InputError, atomic_json, read_bytes, read_fasta
from .reference import build_reference, load_reference


def _metadata(raw: bytes) -> dict:
    """Require an unambiguous TSV header and complete rows before binding labels."""
    try:
        reader = csv.reader(StringIO(raw.decode("utf-8"), newline=""), delimiter="\t", strict=True)
        header = next(reader, [])
        if (len(header) != len(set(header)) or not {"id", "family", "motif"} <= set(header)
                or not set(header) <= {"id", "family", "motif", "accession"}):
            raise InputError("metadata requires unique id, family and motif columns, with optional accession")
        metadata = {}
        for values in reader:
            if not values:  # Ignore blank lines, as csv.DictReader did.
                continue
            if len(values) != len(header):
                raise InputError("metadata row width differs from its header")
            row = dict(zip(header, values))
            if row["id"] in metadata:
                raise InputError("duplicate metadata id")
            metadata[row["id"]] = {key: value or None for key, value in row.items() if key != "id"}
        return metadata
    except (csv.Error, UnicodeError) as exc:
        raise InputError("Malformed metadata TSV") from exc


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--version", action="version", version=f"CRISPRmap-v2 {__version__}")
    commands = result.add_subparsers(dest="command", required=True)
    classify = commands.add_parser("classify", help="map repeats against a frozen local reference")
    classify.add_argument("--input", type=Path, required=True)
    classify.add_argument("--reference", type=Path, required=True)
    classify.add_argument("--output", type=Path, required=True, help="new report directory")
    classify.add_argument("--max-distance", type=int, default=3, help="uncalibrated edit cutoff (default: 3)")
    classify.add_argument("--top-k", type=int, default=20, help="display cap; all best ties determine evidence")
    classify.add_argument("--max-comparisons", type=int, default=5_000_000)
    classify.add_argument("--deadline-seconds", type=float, default=300)
    build = commands.add_parser("build-reference", help="build a provenance-bearing reference from FASTA/TSV")
    build.add_argument("--fasta", type=Path, required=True)
    build.add_argument("--metadata", type=Path, required=True, help="TSV: id, family, motif, optional accession")
    build.add_argument("--namespace", required=True)
    build.add_argument("--source-url", required=True)
    build.add_argument("--output", type=Path, required=True)
    inspect = commands.add_parser("inspect-reference", help="verify reference checksums and show summary")
    inspect.add_argument("--reference", type=Path, required=True)
    catalog_build = commands.add_parser("catalog-build", help="import a frozen source plan into a resumable catalogue")
    catalog_build.add_argument("--plan", type=Path, required=True)
    catalog_build.add_argument("--profile-sha256", required=True)
    catalog_build.add_argument("--output", type=Path, required=True)
    catalog_build.add_argument("--min-free-gib", type=float, default=2)
    merge = commands.add_parser("catalog-merge", help="merge disjoint completed catalogue partitions")
    merge.add_argument("--parts", type=Path, nargs="+", required=True)
    merge.add_argument("--output", type=Path, required=True)
    merge.add_argument("--min-free-gib", type=float, default=2)
    index = commands.add_parser("catalog-index", help="build a complete disk-backed edit-distance index")
    index.add_argument("--catalog", type=Path, required=True)
    index.add_argument("--output", type=Path, required=True)
    index.add_argument("--max-distance", type=int, default=6)
    index.add_argument("--min-similarity", type=float, default=.85)
    indexed = commands.add_parser("catalog-search", help="stream catalogue matches and independent cluster candidates")
    indexed.add_argument("--input", type=Path, required=True)
    indexed.add_argument("--catalog", type=Path, required=True)
    indexed.add_argument("--index", type=Path, required=True)
    indexed.add_argument("--clusters", type=Path)
    indexed.add_argument("--cluster-similarity", type=float, default=.9)
    indexed.add_argument("--output", type=Path, required=True)
    indexed.add_argument("--max-distance", type=int, default=3)
    indexed.add_argument("--top-k", type=int, default=20)
    indexed.add_argument("--deadline-seconds", type=float, default=300)
    forest = commands.add_parser("catalog-cluster-forest", help="build/resume one lane of exact sequence components")
    forest.add_argument("--catalog", type=Path, required=True)
    forest.add_argument("--index", type=Path, required=True)
    forest.add_argument("--output", type=Path, required=True)
    forest.add_argument("--lane", type=int, default=0)
    forest.add_argument("--lanes", type=int, default=1)
    forest.add_argument("--thresholds", nargs="+", default=["0.85", "0.90", "0.95"])
    forest.add_argument("--max-core-nodes", type=int, default=2_000_000)
    forest.add_argument("--deadline-seconds", type=float, default=172800)
    forest.add_argument("--checkpoint-queries", type=int, default=16,
                        help="commit every 1–1024 completed queries in this lane (default: 16)")
    assemble = commands.add_parser("catalog-cluster-assemble", help="validate all graph lanes and assign sequence group IDs")
    assemble.add_argument("--catalog", type=Path, required=True)
    assemble.add_argument("--forests", type=Path, nargs="+", required=True)
    assemble.add_argument("--output", type=Path, required=True)
    fold = commands.add_parser("catalog-fold", help="compute separate pinned RNA folding evidence for core nodes")
    fold.add_argument("--catalog", type=Path, required=True)
    fold.add_argument("--output", type=Path, required=True)
    fold.add_argument("--lane", type=int, default=0)
    fold.add_argument("--lanes", type=int, default=1)
    fold.add_argument("--deadline-seconds", type=float, default=172800)
    fold.add_argument("--checkpoint-nodes", type=int, default=100,
                      help="commit every 1–10000 newly folded nodes (default: 100)")
    benchmark = commands.add_parser("catalog-benchmark", help="measure leakage-aware coverage and engineering correctness")
    benchmark.add_argument("--catalog", type=Path, required=True)
    benchmark.add_argument("--index", type=Path, required=True)
    benchmark.add_argument("--legacy-reference", type=Path)
    benchmark.add_argument("--output", type=Path, required=True)
    benchmark.add_argument("--max-queries", type=int, default=1000)
    benchmark.add_argument("--seed", default="20261001")
    cat_inspect = commands.add_parser("catalog-inspect", help="verify catalogue hash and show frozen provenance")
    cat_inspect.add_argument("--catalog", type=Path, required=True)
    return result


def modern_command(args):
    if hasattr(args, "min_free_gib") and (not math.isfinite(args.min_free_gib) or not 0 <= args.min_free_gib <= 1048576):
        raise ValueError("min-free-gib must be finite and between 0 and 1048576")
    if args.command == "catalog-build":
        from .catalog import build_partition
        return build_partition(args.plan, args.output, args.profile_sha256, min_free_bytes=int(args.min_free_gib * 1024**3))
    if args.command == "catalog-merge":
        from .catalog import merge_catalog
        return merge_catalog(args.parts, args.output, min_free_bytes=int(args.min_free_gib * 1024**3))
    if args.command == "catalog-index":
        from .catalog_index import build_index
        return build_index(args.catalog, args.output, max_distance=args.max_distance, min_similarity=args.min_similarity)
    if args.command == "catalog-search":
        from .indexed_search import search_catalog
        records, input_sha = read_fasta(args.input)
        return search_catalog(records, args.catalog, args.index, args.output, max_distance=args.max_distance,
                              top_k=args.top_k, deadline_seconds=args.deadline_seconds, clusters_dir=args.clusters,
                              cluster_similarity=args.cluster_similarity, input_sha256=input_sha)
    if args.command == "catalog-cluster-forest":
        from .clusters import build_forest
        return build_forest(args.catalog, args.index, args.output, lane=args.lane, lanes=args.lanes,
                            thresholds=args.thresholds, max_core_nodes=args.max_core_nodes,
                            deadline_seconds=args.deadline_seconds, checkpoint_queries=args.checkpoint_queries)
    if args.command == "catalog-cluster-assemble":
        from .clusters import assemble_clusters
        return assemble_clusters(args.catalog, args.forests, args.output)
    if args.command == "catalog-fold":
        from .structure import fold_catalog
        return fold_catalog(args.catalog, args.output, lane=args.lane, lanes=args.lanes,
                            deadline_seconds=args.deadline_seconds, checkpoint_nodes=args.checkpoint_nodes)
    if args.command == "catalog-benchmark":
        from .modern_benchmark import run_benchmark
        return run_benchmark(args.catalog, args.index, args.output, legacy_reference=args.legacy_reference,
                             max_queries=args.max_queries, seed=args.seed)
    from .catalog import load_catalog
    connection, manifest = load_catalog(args.catalog)
    connection.close()
    return manifest


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command.startswith("catalog-"):
            print(json.dumps(modern_command(args), indent=2))
        elif args.command == "classify":
            if args.output.exists() or args.output.is_symlink():
                raise ValueError("output path already exists; choose a new directory")
            records, input_sha = read_fasta(args.input)
            reference, reference_sha = load_reference(args.reference)
            report = search(records, reference, max_distance=args.max_distance, top_k=args.top_k,
                            max_comparisons=args.max_comparisons, deadline_seconds=args.deadline_seconds)
            report["software"] = {"name": "CRISPRmap-v2", "version": __version__}
            report["input"] = {"sha256": input_sha, "record_count": len(records)}
            report["reference"]["file_sha256"] = reference_sha
            report["reference"]["provenance"] = reference["provenance"]
            write_report(report, args.output)
            print(f"Wrote {len(records)} repeat results to {args.output / 'report.json'}")
        elif args.command == "build-reference":
            import hashlib

            records, source_sha = read_fasta(args.fasta)
            raw = read_bytes(args.metadata, limit=16 * 1024 * 1024)
            metadata = _metadata(raw)
            document = build_reference(records, metadata, {
                "namespace": args.namespace, "source_url": args.source_url, "source_sha256": source_sha,
                "metadata_sha256": hashlib.sha256(raw).hexdigest(),
            })
            atomic_json(args.output, document)
            print(f"Built reference with {len(document['entries'])} entries")
        else:
            document, raw_sha = load_reference(args.reference)
            print(json.dumps({"namespace": document["provenance"]["namespace"],
                              "entries": len(document["entries"]),
                              "content_sha256": document["content_sha256"], "file_sha256": raw_sha}, indent=2))
    except (ValueError, OSError, TimeoutError, RuntimeError) as error:
        print(f"crisprmap: {error}", file=sys.stderr)
        return 2
    return 0
