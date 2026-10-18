"""CasAndra production command line."""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from casandra import __version__
from casandra.api import DEFAULT_MODEL_BUNDLE, AnalysisOptions, Analyzer, GenomeInput


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="casandra",
        description="Offline CPU-only Cas protein and cassette predictor",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    inspect = commands.add_parser("inspect-model", help="verify and summarize a model bundle")
    inspect.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path)

    genome = commands.add_parser("predict-genome", help="call/reuse CDSs and predict Cas cassettes")
    genome.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path)
    inputs = genome.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--genome", type=Path, help="nucleotide FASTA")
    inputs.add_argument("--genbank", type=Path, help="GenBank with supplied CDS annotations")
    genome.add_argument("--gff3", type=Path, help="GFF3 CDS annotations for --genome")
    genome.add_argument("--output", required=True, type=Path)
    genome.add_argument("--gene-mode", choices=["auto", "single", "meta"], default="auto")
    genome.add_argument("--translation-table", type=int, default=11)
    genome.add_argument("--threads", type=int, default=6)

    proteins = commands.add_parser("predict-proteins", help="classify amino-acid FASTA records")
    proteins.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path)
    proteins.add_argument("--input", required=True, type=Path)
    proteins.add_argument("--output", required=True, type=Path)
    proteins.add_argument("--threads", type=int, default=6)

    annotation = commands.add_parser(
        "annotate-proteins",
        help="annotate amino-acid FASTA records with atomic provenance-bound outputs",
    )
    annotation.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path)
    annotation.add_argument("--input", required=True, type=Path)
    annotation.add_argument("--output", required=True, type=Path)
    annotation.add_argument("--threads", type=int, default=6)

    cassette = commands.add_parser(
        "classify-cassette",
        help="classify one ordered amino-acid FASTA as a CRISPR-Cas cassette",
    )
    cassette.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path)
    cassette.add_argument("--input", required=True, type=Path)
    cassette.add_argument("--output", required=True, type=Path)
    cassette.add_argument("--threads", type=int, default=6)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "predict-genome" and args.gff3 is not None and args.genome is None:
        raise SystemExit("--gff3 requires --genome")
    analyzer = Analyzer.from_bundle(args.model)
    if args.command == "inspect-model":
        result = analyzer.inspect()
    elif args.command == "predict-proteins":
        result = analyzer.predict_proteins(args.input, args.output, threads=args.threads)
    elif args.command == "annotate-proteins":
        result = analyzer.annotate_proteins(
            args.input, args.output, threads=args.threads
        ).summary
    elif args.command == "classify-cassette":
        result = analyzer.classify_cassette(
            args.input, args.output, threads=args.threads
        ).summary
    elif args.command == "predict-genome":
        result = analyzer.analyze_genome(
            GenomeInput(genome_fasta=args.genome, genbank=args.genbank, gff3=args.gff3),
            args.output,
            AnalysisOptions(
                gene_mode=args.gene_mode,
                translation_table=args.translation_table,
                threads=args.threads,
            ),
        ).summary
    else:
        raise AssertionError(f"Unhandled command: {args.command}")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
