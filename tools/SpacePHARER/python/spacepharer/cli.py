"""Stable offline entry point. Reference preparation is an operator action."""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

from . import __version__
from .engine import EngineError
from .io import InputError, json_value, read
from .pipeline import associate
from .reference import genome_list, prepare, verify


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    ref = sub.add_parser("prepare-reference")
    inputs = ref.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--genome", type=Path, action="append")
    inputs.add_argument("--genome-list", type=Path)
    ref.add_argument("--name", required=True)
    ref.add_argument("--output", type=Path, required=True)
    ref.add_argument("--source-metadata", type=Path)
    ref.add_argument("--scope", choices=("targeted_panel", "broad_phage_reference"), default="targeted_panel")
    ref.add_argument("--build-timeout", type=int, default=14_400)
    check = sub.add_parser("verify-reference")
    check.add_argument("--reference", type=Path, required=True)
    run = sub.add_parser("associate")
    for name in ("input", "groups", "reference", "output"):
        run.add_argument("--" + name, type=Path, required=True)
    run.add_argument("--fdr", type=float, default=0.05)
    run.add_argument("--memory-mb", type=int, default=2048)
    for command in (ref, run):
        command.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare-reference":
            if not 1 <= args.threads <= 8:
                raise InputError("Invalid thread count")
            paths = genome_list(args.genome_list) if args.genome_list else args.genome
            metadata = json_value(read(args.source_metadata, 8192)) if args.source_metadata else None
            prepare(
                paths,
                args.output,
                args.name,
                args.threads,
                source=metadata,
                scope=args.scope,
                build_timeout=args.build_timeout,
            )
        elif args.command == "verify-reference":
            print(json.dumps(verify(args.reference), sort_keys=True))
        else:
            associate(
                args.input, args.groups, args.reference, args.output, args.fdr, args.threads, args.memory_mb
            )
    except (InputError, ValueError, KeyError, TypeError, json.JSONDecodeError, sqlite3.Error) as error:
        print(f"Invalid input: {error}", file=sys.stderr)
        return 2
    except EngineError as error:
        print(str(error), file=sys.stderr)
        return 3
    except OSError as error:
        print(f"I/O failure: {error.strerror}", file=sys.stderr)
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
