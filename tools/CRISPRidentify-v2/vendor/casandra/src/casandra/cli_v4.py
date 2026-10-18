"""Final fast CPU-only CasAndra command line."""

from __future__ import annotations

import json

from casandra.cli_v2 import build_parser
from casandra.pipeline_v3 import predict_genome_v3
from casandra.production_v4 import predict_proteins_v4


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "predict-genome":
        result = predict_genome_v3(args.model, args.output, genome_fasta=args.genome, genbank=args.genbank, gff3=args.gff3, gene_mode=args.gene_mode, translation_table=args.translation_table, threads=args.threads)
    elif args.command == "predict-proteins":
        config = json.loads((args.model / "config.json").read_text(encoding="utf-8"))
        result = predict_proteins_v4(args.model / config["protein_model_directory"], args.input, args.output, threads=args.threads)
    else:
        raise AssertionError(args.command)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

