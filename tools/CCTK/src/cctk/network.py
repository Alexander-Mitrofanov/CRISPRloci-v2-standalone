from __future__ import annotations

#!/usr/bin/env python3

import sys
import argparse

from . import sequence_operations, file_handling

description = """
usage: cctk network [-h] -i <path> [-t <path>] [-o <path>] [--min-shared <int>]

required argument:
  -i, --input   Array_IDs.txt or Array_seqs.txt

optional arguments:
  -h, --help    show this help message and exit
  -t, --types   array CRISPR subtypes file. 2 columns: Array Type
  -o, --outdir  output directory path. Default ./
  --min-shared  minimum spacers shared to draw an edge in network
"""


def build_parser(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "-i",
        "--input",
        metavar="VALUE",
        required=True,
        help="Array_IDs.txt or Array_seqs.txt",
    )
    parser.add_argument(
        "-t",
        "--types",
        metavar="VALUE",
        help="what CRISPR subtype is each array. 2 columns: Array\tType",
    )
    parser.add_argument(
        "-o", "--outdir", metavar="VALUE", default="./", help="output directory path"
    )
    parser.add_argument(
        "--min-shared",
        metavar="VALUE",
        default=1,
        type=int,
        help="Minimum number of spacers shared to draw an edge",
    )

    return parser


def main(args: argparse.Namespace) -> None:
    from pathlib import Path
    from .policy import validate_command
    from .comparison import components

    validate_command(args, "network")
    arrays = file_handling.read_array_file(args.input)
    types = file_handling.read_array_types_file(args.types) if args.types else {}
    if args.types and set(types) != set(arrays):
        raise ValueError("array type annotations must cover exactly the input arrays")
    records = []
    for ident, spacers in arrays.items():
        array = file_handling.FoundArray()
        array.id, array.spacers = ident, spacers
        array.repeat_id = types.get(ident, "")
        records.append(array)
    edges = [
        edge
        for edge in sequence_operations.build_network(records)
        if edge.nshared >= args.min_shared
    ]
    output = Path(args.outdir)
    output.mkdir(parents=True, exist_ok=True)
    file_handling.write_network_file(edges, str(output / "Array_network.txt"), args.min_shared)
    clusters = components(list(arrays), [{"a": edge.a, "b": edge.b} for edge in edges])
    file_handling.write_array_cluster_file(clusters, str(output / "Array_clusters.txt"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a network representation of the relationships \
        between arrays."
    )
    parser = build_parser(parser)

    if len(sys.argv) == 1:
        parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    main(args)
