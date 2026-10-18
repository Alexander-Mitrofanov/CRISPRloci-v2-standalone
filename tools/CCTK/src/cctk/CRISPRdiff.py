from __future__ import annotations
import sys
from collections import defaultdict
from collections.abc import Sequence
from itertools import combinations, permutations
from random import Random
import argparse
from . import colour_schemes, file_handling, plotting
from .io import MAX_ARRAYS


#!/usr/bin/env python3

_rng = Random(2)
sample, randrange, seed = _rng.sample, _rng.randrange, _rng.seed


description = """
usage: cctk crisprdiff [-h] -a <path> -o <path> [--iterations <int>] \
[--preordered] [--approx-ordered] [--seed <int>] [--colour-file <path>] \
[--colour-scheme-outfile <path>] [--colour-scheme-infile <path>] \
[--force-colour-unique] [--line-width <float>] [--dpi <int>] \
[--connection-outline] [--plot-width <float>] [--plot-height <float>] \
[--font-size <float>] [arrays_to_align <list of names>]

optional arguments:
  -h, --help         show this help message and exit

required arguments:
  -a, --array-file   Array_IDs.txt or Array_seqs.txt
  -o, --out-file     output plot file name and path

positional arguments:
  must be at the end of your command
  arrays_to_align    IDs of the arrays you want to analyse. Default: all

running parameters:
  control run behaviour

  --iterations      number of iterations of order determination. Default = 100
  --preordered      array order you provided is the one you want plotted
  --approx-ordered  array order you provided should be optimized slightly
  --seed            set seed for random processes

colour scheme files:
  set inputs and outputs for optional colour scheme files

  --colour-file     file with custom colour list
  --colour-scheme-outfile
                    output file to store json format colour schemes
  --colour-scheme-infile
                    input file json format colour scheme
  --force-colour-unique
                    override black colour of unique spacers. Instead \
use specified colour scheme

plotting parameters:
  control elements of the produced plot

  --line-width      scale factor for lines between identical spacers
  --dpi             resolution of output image
  --connection-outline 
                    add outline colour to lines connecting identical spacers
  --plot-width      width of plot in inches. Default = 3
  --plot-height     height of plot in inches. Default = 3
  --font-size       font size. Defualt automatically scaled.
"""


def _pair_overlaps(arrays_dict: dict[str, list[str]]) -> dict[tuple[str, str], int]:
    """Cache distinct-spacer overlaps once per ordering search."""
    spacers = {ident: set(values) for ident, values in arrays_dict.items()}
    return {(a, b): len(spacers[a] & spacers[b]) for a, b in combinations(spacers, 2)}


def _order_score(overlaps: dict[tuple[str, str], int], order: Sequence[str]) -> int:
    score = 0
    for a, b in zip(order, order[1:]):
        value = overlaps.get((a, b), overlaps.get((b, a)))
        if value is None:
            print("can't find combination %s, %s" % (b, a))
        else:
            score += value
    return score


def get_list_score(arrays_dict: dict[str, list[str]], arrays_order: Sequence[str]) -> int:
    """Sum set intersections of neighbouring arrays; repeated occurrences count once."""
    return _order_score(_pair_overlaps(arrays_dict), arrays_order)


def decide_array_order_global_best(arrays_dict: dict[str, list[str]]) -> list[str]:
    """Exhaustive optimum for at most eight arrays; keep the first tied order."""
    if not 1 <= len(arrays_dict) <= 8:
        raise ValueError("exhaustive ordering supports one to eight arrays")
    elements = list(arrays_dict.keys())
    possible_orders = permutations(elements, len(elements))
    overlaps = _pair_overlaps(arrays_dict)
    score = 0
    best_order = list(elements)
    for i in possible_orders:
        new_score = _order_score(overlaps, i)
        if new_score > score:
            score = int(new_score)
            best_order = list(i)
    return best_order


def jiggle_list_to_local_max(
    arrays_dict: dict[str, list[str]], order: list[str]
) -> tuple[list[str], int]:
    """Retain the original neighbour-swap walk and strict improvement rule."""
    overlaps = _pair_overlaps(arrays_dict)
    repeat = True
    best_order = list(order)
    best_score = _order_score(overlaps, order)
    while repeat:
        repeat = False
        for i in range(len(order)):
            try:
                # switch list element with the element before it
                order[i], order[i - 1] = order[i - 1], order[i]
                new_score = _order_score(overlaps, order)
                if new_score > best_score:
                    repeat = True
                    best_score = int(new_score)
                    best_order = list(order)
            except IndexError:
                pass
            try:
                # switch list element with the element after it
                order[i], order[i + 1] = order[i + 1], order[i]
                new_score = _order_score(overlaps, order)
                if new_score > best_score:
                    repeat = True
                    best_score = int(new_score)
                    best_order = list(order)
            except IndexError:
                pass
    return best_order, best_score


def shuffle_random_arrays(order: list[str]) -> list[str]:
    """Swap two distinct positions in place; a singleton has no possible swap."""
    if len(order) < 2:
        return order
    ran = len(order)
    a = randrange(ran)
    b = randrange(ran)
    while b == a:
        b = randrange(ran)
    order[a], order[b] = order[b], order[a]

    return order


def decide_array_order_local_best(
    arrays_dict: dict[str, list[str]], reps: int, nits: int
) -> tuple[list[str], int]:
    """Seeded swap walk with restarts; this heuristic does not prove optimality.

    Preserve the original random draws, accepted improvements and tie handling.
    Supplied order is the fallback when every candidate has zero sharing.
    """
    if not 1 <= len(arrays_dict) <= MAX_ARRAYS:
        raise ValueError(f"local ordering supports one to {MAX_ARRAYS} arrays")
    if any(type(value) is not int or not 1 <= value <= 1000 for value in (reps, nits)):
        raise ValueError("local ordering repetitions and iterations must be 1–1000")
    elements = list(arrays_dict.keys())
    overall_best_score = 0
    overall_best_order = list(elements)
    overlaps = _pair_overlaps(arrays_dict)
    if len(elements) == 1:
        return overall_best_order, overall_best_score
    for i in range(nits):
        best_score = 0
        rep = 0
        order = sample(elements, len(elements))
        best_order = list(order)
        while rep < reps:
            rep += 1
            order = shuffle_random_arrays(order)
            score = _order_score(overlaps, order)
            if score > best_score:
                rep = 0
                best_score = int(score)
                best_order = list(order)
        if best_score > overall_best_score:
            overall_best_score = int(best_score)
            overall_best_order = list(best_order)
    return overall_best_order, overall_best_score


def build_parser(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        "-a", "--array-file", required=True, help="Specify array representatives file."
    )
    parser.add_argument("-o", "--out-file", required=True, help="Specify output file.")

    run_params = parser.add_argument_group("Running parameters", "Control run behaviour.")
    run_params.add_argument(
        "--iterations",
        required=False,
        default=100,
        metavar="VALUE",
        type=int,
        help="Number of attempts to find the best order of arrays for \
        plotting. Higher number can improve the array order. Default = 100",
    )
    run_params.add_argument(
        "--optimize-order", action="store_true", help="Explicitly run bounded order optimization"
    )
    run_params.add_argument(
        "--preordered",
        required=False,
        action="store_true",
        help="Declare that the array order you provided is the one you \
        want plotted.",
    )
    run_params.add_argument(
        "--approx-ordered",
        required=False,
        action="store_true",
        help="Declare that the array order you provided should be \
        optimized slightly before plotting. Optimization involves \
        switching order of adjacent arrays in list as long as that \
        increases the total number of shared spacers among all \
        neighbours.",
    )
    run_params.add_argument(
        "--seed",
        metavar="VALUE",
        required=False,
        type=int,
        default=2,
        help="The order of outline and fill colours assigned to \
        spacers is semi-random. Change it by providing a number here \
        to change which colours are assigned to each spacer.",
    )

    cs_files = parser.add_argument_group(
        "Colour scheme files",
        "Set inputs and outputs for optional colour scheme files.",
    )
    cs_files.add_argument(
        "--colour-file",
        metavar="VALUE",
        required=False,
        help="Specify file with custom colour list (Optional). \
        Colours must be hex codes. One colour per line with no header \
        line in file. e.g. #fd5925.",
    )
    cs_files.add_argument(
        "--colour-scheme-outfile",
        metavar="VALUE",
        required=False,
        help="Specify output file to store json format dictionary of \
        the colour schemes used for spacers in this run.",
    )
    cs_files.add_argument(
        "--colour-scheme-infile",
        metavar="VALUE",
        required=False,
        help="Specify input file containing json format dictionary of \
        the colour scheme to be used for spacers in this run. Any \
        spacers not in the input file will be coloured according to \
        the normal process.",
    )
    cs_files.add_argument(
        "--force-colour-unique",
        required=False,
        action="store_true",
        help="override black colour of unique spacers. Instead use \
        specified colour scheme",
    )

    plot_params = parser.add_argument_group(
        "Plotting parameters", "Control elements of the produced plot."
    )
    plot_params.add_argument(
        "--line-width",
        metavar="VALUE",
        required=False,
        default=1,
        type=float,
        help="Control the width of lines connecting shared spacers. Line width \
        will be multiplied by the given number. Default = 1.0",
    )
    plot_params.add_argument(
        "--dpi",
        metavar="VALUE",
        required=False,
        type=int,
        default=300,
        help="Resolution of output image. Only relevant for bitmap \
        formats such as PNG. Has no effect on SVG outputs.",
    )
    plot_params.add_argument(
        "--connection-outline",
        required=False,
        action="store_true",
        help="Add outline colour to lines connecting identical spacers.",
    )
    plot_params.add_argument(
        "--plot-width",
        type=float,
        default=3,
        metavar="VALUE",
        help="Width of plot in inches. Default = 3",
    )
    plot_params.add_argument(
        "--plot-height",
        type=float,
        default=3,
        metavar="VALUE",
        help="Height of plot in inches. Default = 3",
    )
    plot_params.add_argument(
        "--font-size",
        type=float,
        metavar="VALUE",
        default=0,
        help="Set font size. Defualt automatically scaled.",
    )

    parser.add_argument(
        "arrays_to_align",
        nargs="*",
        help="Specify the arrays for which you want to plot an \
        alignment. **Must come at the end of your command after all \
        other arguments.**",
    )

    return parser


def main(args: argparse.Namespace) -> None:
    from .policy import validate_command

    validate_command(args, "crisprdiff")

    array_dict = file_handling.read_array_file(args.array_file)

    if len(args.arrays_to_align) == 0:
        arrays_of_interest_dict = array_dict
        arrays_of_interest = [a for a in array_dict.keys()]
    else:
        arrays_of_interest_dict = {}
        arrays_of_interest = args.arrays_to_align
        for array in arrays_of_interest:
            arrays_of_interest_dict[array] = array_dict[array]

    seed(args.seed)
    if args.preordered or (not args.approx_ordered and not getattr(args, "optimize_order", False)):
        # Reverse order so order of list plots top to bottom
        array_order = arrays_of_interest[::-1]
    elif args.approx_ordered:
        array_order, _ = jiggle_list_to_local_max(arrays_of_interest_dict, arrays_of_interest[::-1])
    else:
        if len(arrays_of_interest) < 9:
            array_order = decide_array_order_global_best(arrays_of_interest_dict)
        else:
            array_order, local_best_score = decide_array_order_local_best(
                arrays_of_interest_dict, 100, args.iterations
            )
            print(
                "The score (sum of spacers shared between all neighbouring "
                "arrays) of the best ordering found was: {}".format(local_best_score)
            )

    ## Find spacers present in more than one array
    all_spacers = list(arrays_of_interest_dict.values())
    occurrences = defaultdict(int)
    for array in all_spacers:
        for spacer in dict.fromkeys(array):
            occurrences[spacer] += 1

    imp_spacers = []

    for k, v in occurrences.items():
        if v > 1:
            imp_spacers.append(k)

    sys.stderr.write(
        "Identified {} spacers present in more than one array.\n".format(len(imp_spacers))
    )

    # Process colour file related args
    spacer_colours = colour_schemes.process_colour_args(args, imp_spacers, args.force_colour_unique)

    plotting.plot_diffplot(
        arrays_of_interest_dict,
        array_order,
        imp_spacers,
        spacer_colours,
        text_size=args.font_size,
        plot_width=args.plot_width,
        plot_height=args.plot_height,
        dpi=args.dpi,
        outfile=args.out_file,
        line_width=args.line_width,
        connection_outline=args.connection_outline,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Given an array file and a list of arrays to "
        "align, produces a plot showing shared spacers among arrays "
        "and their location within each array. Array file format is "
        "whitespace separated columns where column 1 is the array ID "
        "and columns 3 onwards are spacer IDs, names, or sequences."
        "\ne.g. CRISPRdiff.py -a Array_IDs.txt -o "
        "test.png -c cols.txt -iter 100 155 9 204 73 97",
    )
    parser = build_parser(parser)

    if len(sys.argv) == 1:
        parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    main(args)
