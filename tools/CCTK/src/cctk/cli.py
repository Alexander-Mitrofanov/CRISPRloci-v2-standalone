"""Installed CCTK entry point with explicit subcommands and failure codes."""

from __future__ import annotations
import argparse
import importlib
from pathlib import Path
import sys
import signal
import os
import tempfile
from . import __version__
from .errors import CCTKError

COMMANDS = {
    "network": "network",
    "crisprdiff": "CRISPRdiff",
    "crisprtree": "CRISPRtree",
    "constrain": "constrain",
    "evolve": "evolve",
    "blast": "blast",
    "spacerblast": "spacerblast",
    "minced": "minced",
    "quickrun": "quickrun",
}


def parser() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser(
        prog="cctk", description="Observed CRISPR comparisons and explicit evolutionary analyses."
    )
    cli.add_argument("--version", action="version", version=f"CCTK {__version__}")
    commands = cli.add_subparsers(dest="command", required=True)
    compare = commands.add_parser(
        "compare",
        help="Structured DNA-occurrence comparison with gaps and explicit identity policy",
    )
    compare.add_argument("-i", "--input", type=Path, required=True)
    compare.add_argument("-o", "--outdir", type=Path, required=True)
    compare.add_argument("--identity", choices=["exact", "reverse-complement"], default="exact")
    compare.add_argument("--min-shared", type=int, default=1)
    compare.add_argument("--order", choices=["supplied", "id"], default="supplied")
    compare.add_argument("--plot", action="store_true")
    compare.add_argument("--seed", type=int, default=2)
    for command, module in COMMANDS.items():
        implementation = importlib.import_module("." + module, __package__)
        child = commands.add_parser(command, help=implementation.__doc__)
        implementation.build_parser(child)
        child.set_defaults(run=implementation.main)
    return cli


def main(argv: list[str] | None = None) -> int:
    def terminate(signum: int, frame: object) -> None:
        raise CCTKError("analysis terminated")

    previous = signal.signal(signal.SIGTERM, terminate)
    cache = None
    previous_cache = os.environ.get("MPLCONFIGDIR")
    try:
        # Apply only to the CLI, before command imports reach pyplot. A worker
        # may supply its own private cache; otherwise each invocation owns one.
        if previous_cache is None:
            cache = tempfile.TemporaryDirectory(prefix="cctk-matplotlib-")
            os.environ["MPLCONFIGDIR"] = cache.name
        import matplotlib

        matplotlib.use("Agg", force=True)
        args = parser().parse_args(argv)
        if args.command == "compare":
            from .comparison import run_compare

            run_compare(
                args.input,
                args.outdir,
                identity=args.identity,
                min_shared=args.min_shared,
                order=args.order,
                plot=args.plot,
                seed=args.seed,
            )
        else:
            args.run(args)
        return 0
    except (CCTKError, OSError, KeyError, ValueError, RecursionError) as error:
        print(f"CCTK error: {error}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130
    finally:
        signal.signal(signal.SIGTERM, previous)
        if cache is not None:
            cache.cleanup()
            os.environ.pop("MPLCONFIGDIR", None)
