"""CasAndra Cas-only production command line."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from casandra.pipeline import predict_genome
from casandra.production_v3 import predict_proteins_v3

DEFAULT_MODEL_BUNDLE = Path(__file__).resolve().parent / "models"


def build_parser():
    parser = argparse.ArgumentParser(prog="casandra", description="Offline CPU-only Cas protein and cassette predictor")
    commands = parser.add_subparsers(dest="command", required=True)
    genome = commands.add_parser("predict-genome", help="call/reuse CDSs and predict Cas cassettes")
    genome.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path,
                        help="CasAndra model bundle (default: bundled model)")
    inputs = genome.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--genome", type=Path, help="nucleotide FASTA")
    inputs.add_argument("--genbank", type=Path, help="GenBank with supplied CDS annotations")
    genome.add_argument("--gff3", type=Path, help="GFF3 CDS annotations for --genome")
    genome.add_argument("--output", required=True, type=Path)
    genome.add_argument("--gene-mode", choices=["auto", "single", "meta"], default="auto")
    genome.add_argument("--translation-table", type=int, default=11)
    genome.add_argument("--threads", type=int, default=6)
    proteins = commands.add_parser("predict-proteins", help="classify amino-acid FASTA records")
    proteins.add_argument("--model", default=DEFAULT_MODEL_BUNDLE, type=Path,
                          help="CasAndra model bundle (default: bundled model)")
    proteins.add_argument("--input", required=True, type=Path)
    proteins.add_argument("--output", required=True, type=Path)
    proteins.add_argument("--threads", type=int, default=6)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "predict-genome":
        if args.gff3 is not None and args.genome is None:
            raise SystemExit("--gff3 requires --genome")
        result = predict_genome(
            args.model,
            args.output,
            genome_fasta=args.genome,
            genbank=args.genbank,
            gff3=args.gff3,
            gene_mode=args.gene_mode,
            translation_table=args.translation_table,
            threads=args.threads,
        )
    elif args.command == "predict-proteins":
        config = json.loads((args.model / "config.json").read_text(encoding="utf-8"))
        result = predict_proteins_v3(
            args.model / config["protein_model_directory"],
            args.input,
            args.output,
            threads=args.threads,
        )
    else:
        raise AssertionError(args.command)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0
