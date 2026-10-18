"""Small stable command-line interface; no network or environment mutation."""
import argparse
import sys
from pathlib import Path

from . import __version__
from .errors import InputError, ModelError, OutputError
from .inference import classify
from .input import read_fasta
from .output import publish


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Classify repeat FASTA using the frozen CRISPRCasTyper model")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("classify")
    command.add_argument("--input", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--threshold", type=float, default=0.75)
    args = parser.parse_args(argv)
    try:
        records, digest = read_fasta(args.input.absolute())
        report = classify(records, input_sha256=digest, threshold=args.threshold)
        publish(report, args.output.absolute())
    except (InputError, ModelError, OutputError, OSError) as error:
        print(f"repeattyper: {error}", file=sys.stderr)
        return 2 if isinstance(error, InputError) else 3 if isinstance(error, ModelError) else 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
