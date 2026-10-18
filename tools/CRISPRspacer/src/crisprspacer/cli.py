"""The standalone command-line boundary. Errors never masquerade as zero hits."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path
import sys

from . import __version__
from .artifacts import build_reference, load_inventory, load_reference, prepare_bundle, validate_result, write_result
from .inputs import CATEGORIES, from_identify, from_spacers
from .io import InputError, LimitError, metadata, read_fasta, read_fasta_snapshot
from .search import Policy


def search_options(parser):
    parser.add_argument("--substitutions", type=int, default=0, choices=range(3), help="Full-length substitutions only (default: exact)")
    parser.add_argument("--min-length", type=int, default=1)
    parser.add_argument("--max-length", type=int, default=1000)
    parser.add_argument("--category", action="append", choices=sorted(CATEGORIES), help="Repeat to select categories; default Bona-fide and Possible")
    parser.add_argument("--flank", type=int, default=10, help="Forward-target flank context; no PAM inference")
    parser.add_argument("--array-padding", type=int, default=0, help="Explicit additional exclusion bases around known arrays")
    parser.add_argument("--max-hits", type=int, default=100_000)
    parser.add_argument("--max-candidates", type=int, default=2_000_000)
    parser.add_argument("--max-pairs", type=int, default=100_000)
    parser.add_argument("--max-evidence-bytes", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--max-work", type=int, default=500_000_000)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, required=True, help="New output directory; existing paths are never overwritten")


def source_options(parser):
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--identify", type=Path, nargs="+", required=True, help="Canonical CRISPRidentify report.json files")
    parser.add_argument("--metadata", type=Path, help="Explicit source assembly/topology JSON")


def parser():
    top = argparse.ArgumentParser(prog="crisprspacer", description="Observational CRISPR spacer sequence analysis; offline and reproducible.")
    top.add_argument("--version", action="version", version=f"CRISPRspacer {__version__}")
    commands = top.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare", help="Validate canonical array reports and export an occurrence inventory")
    source_options(prep)
    prep.add_argument("--output", type=Path, required=True)
    self_mode = commands.add_parser("self", help="Search source-bound spacers in their source records/declared assemblies")
    source_options(self_mode)
    search_options(self_mode)
    self_mode.add_argument("--scope", choices=["same-record", "same-assembly", "all"], default="same-record")
    search = commands.add_parser("search", help="Search supplied spacer occurrences against supplied/reference target DNA")
    queries = search.add_mutually_exclusive_group(required=True)
    queries.add_argument("--spacers", type=Path, help="FASTA with unique occurrence IDs; duplicate sequences are retained")
    queries.add_argument("--inventory", type=Path, help="Prepared inventory directory or inventory.json")
    targets = search.add_mutually_exclusive_group(required=True)
    targets.add_argument("--targets", type=Path)
    targets.add_argument("--reference", type=Path)
    search.add_argument("--target-metadata", type=Path)
    search_options(search)
    index = commands.add_parser("index", help="Publish a verified immutable local reference sequence inventory")
    index.add_argument("--targets", type=Path, required=True)
    index.add_argument("--metadata", type=Path)
    index.add_argument("--name", required=True)
    index.add_argument("--source-label", default="user-supplied", help="Source/accession collection label without whitespace")
    index.add_argument("--release", default="unspecified", help="Source release identifier; explicitly unspecified by default")
    index.add_argument("--output", type=Path, required=True)
    validate = commands.add_parser("validate", help="Check hashes, source metrics, exports and deterministic replay")
    validate.add_argument("result", type=Path)
    reference_check = commands.add_parser("validate-reference", help="Verify all records and manifest in a reference bundle")
    reference_check.add_argument("reference", type=Path)
    db_index = commands.add_parser("database-index", help="Build an offline SQLite spacer index (operator command)")
    db_index.add_argument("--spacers", type=Path, required=True)
    db_index.add_argument("--name", required=True)
    db_index.add_argument("--source-url", required=True)
    db_index.add_argument("--source-release", required=True)
    db_index.add_argument("--output", type=Path, required=True)
    db_search = commands.add_parser("database-search", help="Check target DNA against an indexed spacer collection; exact only")
    db_search.add_argument("--database", type=Path, required=True)
    db_search.add_argument("--targets", type=Path, required=True)
    db_search.add_argument("--expected-database-sha256")
    db_search.add_argument("--timeout", type=int, default=300)
    db_search.add_argument("--max-hits", type=int, default=10_000)
    db_search.add_argument("--max-candidates", type=int, default=2_000_000)
    db_search.add_argument("--output", type=Path, required=True)
    db_check = commands.add_parser("validate-database", help="Verify an immutable SQLite spacer collection")
    db_check.add_argument("database", type=Path)
    db_check.add_argument("--expected-database-sha256")
    db_result = commands.add_parser("validate-database-result", help="Replay a database search against its pinned snapshot")
    db_result.add_argument("result", type=Path)
    db_result.add_argument("--database", type=Path, required=True)
    viral_plan = commands.add_parser("viral-plan", help="Operator metadata-only plan for the current NCBI viral RefSeq release")
    viral_plan.add_argument("--output", type=Path, required=True)
    viral_download = commands.add_parser("viral-download", help="Operator-only download of all genomic files in a frozen NCBI plan")
    viral_download.add_argument("--plan", type=Path, required=True)
    viral_download.add_argument("--output", type=Path, required=True)
    viral_download.add_argument("--min-free-bytes", type=int, default=20 * 1024**3)
    viral_index = commands.add_parser("viral-index", help="Build an immutable local viral BLAST reference from supplied FASTA")
    viral_index.add_argument("--sources", type=Path, nargs="+", required=True)
    viral_index.add_argument("--source-manifest", type=Path)
    viral_index.add_argument("--name", required=True)
    viral_index.add_argument("--source-release", required=True)
    viral_index.add_argument("--source-url", required=True)
    viral_index.add_argument("--makeblastdb", default="makeblastdb")
    viral_index.add_argument("--blastdbcmd", default="blastdbcmd")
    viral_index.add_argument("--timeout", type=int, default=3600)
    viral_index.add_argument("--output", type=Path, required=True)
    viral_check = commands.add_parser("validate-viral-reference", help="Verify viral reference checksums and native sequence roundtrip")
    viral_check.add_argument("reference", type=Path)
    viral_check.add_argument("--expected-reference-sha256")
    viral_check.add_argument("--blastdbcmd", help="Override the recorded blastdbcmd path with an identical executable")
    viral_check.add_argument("--timeout", type=int, default=3600)
    viral_search = commands.add_parser("viral-search", help="Search supplied spacers for viral reference matches using local short-query BLAST")
    viral_search.add_argument("--spacers", type=Path, required=True)
    viral_search.add_argument("--reference", type=Path, required=True)
    viral_search.add_argument("--expected-reference-sha256")
    viral_search.add_argument("--blastn", default="blastn")
    viral_search.add_argument("--max-mismatches", type=int, choices=range(3), default=2)
    viral_search.add_argument("--timeout", type=int, default=300)
    viral_search.add_argument("--threads", type=int, default=2)
    viral_search.add_argument("--max-hits", type=int, default=10_000)
    viral_search.add_argument("--max-hsps", type=int, default=200_000, help="Raw HSP budget; exceeding it fails the result")
    viral_search.add_argument("--output", type=Path, required=True)
    viral_result = commands.add_parser("validate-viral-result", help="Validate retained BLAST evidence against the pinned reference")
    viral_result.add_argument("result", type=Path)
    viral_result.add_argument("--reference", type=Path, required=True)
    viral_result.add_argument("--expected-reference-sha256")
    return top


def main(argv=None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "viral-plan":
            from .ncbi_viral import plan_ncbi_viral
            result = plan_ncbi_viral(args.output)
        elif args.command == "viral-download":
            from .ncbi_viral import download_ncbi_viral
            result = download_ncbi_viral(args.plan, args.output, min_free_bytes=args.min_free_bytes)
        elif args.command == "viral-index":
            from .viral_reference import build_viral_reference
            result = build_viral_reference(args.sources, args.output, name=args.name, source_release=args.source_release,
                                          source_url=args.source_url, source_manifest=args.source_manifest,
                                          makeblastdb=args.makeblastdb, blastdbcmd=args.blastdbcmd, timeout=args.timeout)
        elif args.command == "validate-viral-reference":
            from .viral_reference import validate_viral_reference
            result = validate_viral_reference(args.reference, expected=args.expected_reference_sha256,
                                              blastdbcmd=args.blastdbcmd, timeout=args.timeout)
        elif args.command == "viral-search":
            from .viral_search import search_viruses
            report = search_viruses(args.spacers, args.reference, args.output, blastn=args.blastn,
                                   expected_reference_sha256=args.expected_reference_sha256,
                                   max_mismatches=args.max_mismatches, timeout=args.timeout, threads=args.threads,
                                   max_hits=args.max_hits, max_hsps=args.max_hsps)
            result = {"status": "completed", "outcome": report["outcome"], "counts": report["counts"]}
        elif args.command == "validate-viral-result":
            from .viral_search import validate_viral_search
            result = validate_viral_search(args.result, reference=args.reference,
                                          expected_reference_sha256=args.expected_reference_sha256)
        elif args.command in {"database-index", "database-search", "validate-database", "validate-database-result"}:
            from .database import build_database, search_database, validate_database, validate_database_result
            if args.command == "database-index":
                result = build_database(args.spacers, args.output, args.name, args.source_url, args.source_release)
            elif args.command == "database-search":
                report = search_database(args.database, args.targets, args.output, args.expected_database_sha256,
                                         timeout=args.timeout, max_hits=args.max_hits, max_candidates=args.max_candidates)
                result = {"status": "completed", "outcome": report["outcome"], "counts": report["counts"]}
            elif args.command == "validate-database":
                ref = validate_database(args.database, args.expected_database_sha256)
                result = {"status": "valid", "reference": ref}
            else:
                result = validate_database_result(args.result, args.database)
        elif args.command == "validate":
            result = validate_result(args.result)
        elif args.command == "validate-reference":
            records, reference = load_reference(args.reference)
            result = {"status": "valid", "records": len(records), "reference": reference}
        elif args.command == "index":
            targets, raw_hash = read_fasta_snapshot(args.targets)
            targets = metadata(targets, args.metadata)
            build_reference(args.output, targets, args.name, raw_hash, args.source_label, args.release)
            result = {"status": "completed", "records": len(targets)}
        else:
            reference = None
            if args.command in {"self", "prepare"}:
                sources = metadata(read_fasta(args.source), args.metadata)
                inventory = from_identify(sources, args.identify)
                targets = sources
            else:
                if args.inventory:
                    inventory = load_inventory(args.inventory)
                else:
                    records, raw_hash = read_fasta_snapshot(args.spacers)
                    inventory = from_spacers(records, raw_hash)
                if args.reference:
                    if args.target_metadata:
                        raise InputError("Reference metadata is immutable; rebuild a new reference to change it")
                    targets, reference = load_reference(args.reference)
                else:
                    targets = metadata(read_fasta(args.targets), args.target_metadata)
            if args.command == "prepare":
                prepare_bundle(args.output, inventory)
                result = {"status": "completed", "query_occurrences": len(inventory["queries"]), "arrays": len(inventory["arrays"])}
            else:
                policy = Policy(substitutions=args.substitutions, min_length=args.min_length, max_length=args.max_length,
                    categories=tuple(args.category or ("Bona-fide", "Possible")), scope=getattr(args, "scope", "all"),
                    flank=args.flank, array_padding=args.array_padding, max_hits=args.max_hits,
                    max_candidates=args.max_candidates, max_pairs=args.max_pairs, max_evidence_bytes=args.max_evidence_bytes,
                    max_work=args.max_work, timeout=args.timeout)
                report = write_result(args.output, inventory, targets, policy, args.command == "self", reference)
                result = {"status": "completed", "outcome": report["outcome"], "counts": report["counts"]}
        print(json.dumps(result, sort_keys=True, allow_nan=False))
        return 0
    except LimitError as exc:
        print(json.dumps({"status": "resource_limit", "error": str(exc), "search_complete": False}), file=sys.stderr)
        return 3
    except (InputError, OSError, TypeError, KeyError, AttributeError, RecursionError, sqlite3.Error) as exc:
        message = str(exc) if isinstance(exc, InputError) else "Malformed input, inaccessible file, or invalid contract structure"
        print(json.dumps({"status": "error", "error": message, "search_complete": False}), file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print(json.dumps({"status": "cancelled", "search_complete": False}), file=sys.stderr)
        return 130
