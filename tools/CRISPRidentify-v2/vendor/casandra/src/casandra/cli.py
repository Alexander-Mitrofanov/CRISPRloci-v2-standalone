"""Command-line interface for the final CasAndra predictor."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from casandra.predictor import build_production_model, predict_proteins


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="casandra", description="Local CPU-only Cas predictor")
    commands = parser.add_subparsers(dest="command", required=True)
    build = commands.add_parser("build-model", help="freeze the post-benchmark production model")
    build.add_argument("--dataset", required=True, type=Path)
    build.add_argument("--profile-catalogue", required=True, type=Path)
    build.add_argument("--alignments", required=True, type=Path)
    build.add_argument("--benchmark-metrics", required=True, type=Path)
    build.add_argument("--output", required=True, type=Path)
    predict = commands.add_parser("predict-proteins", help="classify amino-acid FASTA records")
    predict.add_argument("--model", required=True, type=Path)
    predict.add_argument("--input", required=True, type=Path)
    predict.add_argument("--output", required=True, type=Path)
    predict.add_argument("--threads", type=int, default=6)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "build-model":
        result = build_production_model(
            args.dataset,
            args.profile_catalogue,
            args.alignments,
            args.benchmark_metrics,
            args.output,
        )
    elif args.command == "predict-proteins":
        result = predict_proteins(args.model, args.input, args.output, threads=args.threads)
    else:
        raise AssertionError(args.command)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
