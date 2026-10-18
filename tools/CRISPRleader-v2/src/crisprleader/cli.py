"""Standalone offline CLI; errors are distinct from completed scientific no-calls."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

from . import __version__
from .errors import InputError, LeaderError
from .models import AnalysisConfig


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="crisprleader", description="Reference-supported CRISPR leader evidence")
    root.add_argument("--version", action="version", version=f"CRisprleader v2 {__version__}")
    commands = root.add_subparsers(dest="command", required=True)
    extract = commands.add_parser("extract", help="Extract both context windows; no leader prediction")
    extract.add_argument("--fasta", required=True, type=Path)
    extract.add_argument("--arrays", required=True, type=Path)
    extract.add_argument("--output", required=True, type=Path, help="New output directory")
    extract.add_argument("--flank-length", type=int, default=500)
    predict = commands.add_parser("predict", help="Unavailable until scientific implementation/evaluation is complete")
    predict.add_argument("--fasta", required=True, type=Path)
    predict.add_argument("--arrays", required=True, type=Path)
    predict.add_argument("--references", type=Path, help="Reserved for validated biological references")
    predict.add_argument("--output", required=True, type=Path, help="New output directory")
    predict.add_argument("--flank-length", type=int, default=500)
    predict.add_argument("--threads", type=int, default=1)
    adapt = commands.add_parser("adapt", help="Convert source-bound scientific annotations")
    adapters = adapt.add_subparsers(dest="adapter", required=True)
    identify = adapters.add_parser("identify", help="Read canonical CRISPRidentify reports")
    identify.add_argument("--fasta", required=True, type=Path)
    identify.add_argument("--reports", nargs="+", required=True, type=Path)
    identify.add_argument("--output", required=True, type=Path)
    validate = commands.add_parser("validate", help="Independently validate annotations or a result")
    validate.add_argument("--fasta", required=True, type=Path)
    choice = validate.add_mutually_exclusive_group(required=True)
    choice.add_argument("--arrays", type=Path)
    choice.add_argument("--result", type=Path)
    references = commands.add_parser("references", help="Inspect immutable local references")
    actions = references.add_subparsers(dest="reference_command", required=True)
    verify = actions.add_parser("verify", help="Verify schema, provenance and file hashes")
    verify.add_argument("--bundle", type=Path)
    benchmark = commands.add_parser("benchmark", help="Describe supplied orientation labels/calls; no superiority claim")
    benchmark.add_argument("--dataset", required=True, type=Path)
    benchmark.add_argument("--predictions", required=True, type=Path)
    benchmark.add_argument("--output", required=True, type=Path, help="New JSON report filename")
    benchmark.add_argument("--bootstrap-replicates", type=int, default=2000)
    benchmark.add_argument("--seed", type=int, default=0)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "benchmark":
            from .evaluation import evaluate
            from .export import write_json_report
            from .validation import read_json
            dataset = read_json(args.dataset, limit=16 * 1024 * 1024)
            predictions = read_json(args.predictions, limit=32 * 1024 * 1024)
            report = evaluate(dataset, predictions, bootstrap_replicates=args.bootstrap_replicates,
                              seed=args.seed)
            report["evaluator"] = {"name": "crisprleader benchmark", "version": __version__,
                                   "analysis_kind": "descriptive_orientation_report"}
            # Receipt covers the parsed content actually evaluated, without a second file read.
            report["input_document_receipts"] = {
                name: {"canonical_json_sha256": hashlib.sha256(json.dumps(
                    value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
                ).encode("utf-8")).hexdigest()}
                for name, value in (("dataset", dataset), ("predictions", predictions))
            }
            target = write_json_report(report, args.output)
            print(json.dumps({"execution_status": "completed", "output": str(target),
                              "scientific_superiority_established": False}))
        elif args.command == "extract":
            from .pipeline import extract
            from .export import write_bundle
            report = extract(args.fasta, args.arrays, config=AnalysisConfig(flank_length=args.flank_length))
            target = write_bundle(report, args.output)
            print(json.dumps({"execution_status": "completed", "analysis_mode": "context_extraction",
                              "output": str(target), "arrays": len(report["arrays"]),
                              "candidates": len(report["candidates"]), "prediction_performed": False}))
        elif args.command == "predict":
            from .pipeline import predict
            from .export import write_bundle
            report = predict(args.fasta, args.arrays, references=args.references,
                             config=AnalysisConfig(flank_length=args.flank_length, threads=args.threads))
            target = write_bundle(report, args.output)
            counts: dict[str, int] = {}
            for row in report["predictions"]:
                counts[row["status"]] = counts.get(row["status"], 0) + 1
            print(json.dumps({"execution_status": "completed", "output": str(target), "decisions": counts}))
        elif args.command == "adapt":
            from .adapters.identify import adapt_identify
            from .export import write_json_report
            document = adapt_identify(args.fasta, args.reports)
            write_json_report(document, args.output)
            print(json.dumps({"status": "valid", "arrays": len(document["arrays"])}))
        elif args.command == "references":
            from .references import load_bundle
            bundle = load_bundle(args.bundle)
            print(json.dumps({"status": "valid", "bundle_id": bundle["bundle_id"],
                              "version": bundle["version"], "families": len(bundle["families"]),
                              "motifs": len(bundle["motifs"])}))
        elif args.command == "validate":
            from .sequence_io import load_fasta, load_inputs
            from .validation import read_json, validate_result
            from dataclasses import replace
            if args.arrays:
                _, arrays = load_inputs(args.fasta, args.arrays)
                print(json.dumps({"status": "valid", "arrays": len(arrays)}))
            else:
                document = read_json(args.result)
                sources = load_fasta(args.fasta)
                rows = document.get("sources")
                if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
                    raise InputError("Result sources must be a list of objects")
                for row in rows:
                    if isinstance(row.get("id"), str) and row["id"] in sources:
                        sources[row["id"]] = replace(sources[row["id"]], topology=row.get("topology", "unknown"),
                                                     topology_provenance=row.get("topology_provenance"))
                validate_result(document, sources)
                print(json.dumps({"status": "valid", "predictions": len(document["predictions"])}))
        return 0
    except (LeaderError, OSError, ValueError) as error:
        print(json.dumps({"execution_status": "failed", "error": type(error).__name__, "message": str(error)}),
              file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print('{"execution_status":"cancelled","message":"Interrupted"}', file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
