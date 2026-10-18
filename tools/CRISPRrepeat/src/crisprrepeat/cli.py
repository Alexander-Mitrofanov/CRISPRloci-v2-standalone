"""Standalone CLI; all scientific calculations run in the caller's isolated worker."""

import argparse
from pathlib import Path
import signal
import sys

from . import __version__
from .analysis import analyze
from .exports import write_outputs
from .folding import Settings


def main(argv=None):
    parser = argparse.ArgumentParser(description="CRISPR repeat structure and array-context evidence")
    parser.add_argument("--version", action="version", version=f"CRISPRrepeat {__version__}")
    parser.add_argument("--input", type=Path, required=True, help="FASTA file or folder")
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument(
        "--reports", type=Path, help="Canonical CRISPRidentify reports; input is original genome FASTA"
    )
    parser.add_argument("--molecule", choices=["DNA", "RNA"], default="DNA")
    parser.add_argument("--orientation", choices=["both", "forward", "reverse"])
    parser.add_argument("--categories", nargs="+", default=["Bona-fide", "Possible"])
    parser.add_argument("--temperature", type=float, default=37)
    parser.add_argument("--window", type=int, default=150)
    parser.add_argument("--span", type=int, default=80)
    parser.add_argument("--no-lonely-pairs", action="store_true")
    parser.add_argument("--reference", type=Path, help="Optional operator-owned frozen reference JSON")
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args(argv)
    try:
        if not 1 <= args.timeout <= 3600:
            raise ValueError("Timeout must be 1–3600 seconds")
        if args.output.exists() or args.output.is_symlink():
            raise ValueError("Output must not exist")
        # Default SIGALRM terminates even inside native computation. No partial
        # report is ever mistaken for a completed job; worker containment remains outermost.
        signal.alarm(args.timeout)
        report = analyze(
            args.input,
            reports=args.reports,
            molecule=args.molecule,
            orientation=args.orientation or ("forward" if args.molecule == "RNA" else "both"),
            settings=Settings(args.temperature, args.window, args.span, args.no_lonely_pairs),
            categories=args.categories,
            reference=args.reference,
        )
        write_outputs(args.output, report)
        print(f"Completed: {report['counts']['instances']} repeat/orientation instances")
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(f"CRISPRrepeat: {error}", file=sys.stderr)
        return 1
    finally:
        signal.alarm(0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
