"""Command-line entry points; importing the package never starts analysis."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from . import __version__
from .config import Settings, load_settings
from .errors import TracrError
from .pipeline import doctor, run


def boolean(value: str) -> bool:
    normalized = value.lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise argparse.ArgumentTypeError("Expected true or false")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(
        description="CRISPRtracrRNA v3: reproducible tracrRNA candidates"
    )
    root.add_argument("--version", action="version", version=f"CRISPRtracrRNA {__version__}")
    commands = root.add_subparsers(dest="command", required=True)
    for name in ("run", "doctor"):
        command = commands.add_parser(
            name, help="Run analysis" if name == "run" else "Check dependencies"
        )
        command.add_argument("--config", type=Path, help="TOML configuration file")
        command.add_argument(
            "--mode", "--run_type", choices=["complete", "model", "complete_run", "model_run"]
        )
        command.add_argument("--model-type", "--model_type", type=str.upper, choices=["II", "V"])
        command.add_argument("--threads", type=int)
        command.add_argument("--timeout", type=float)
        command.add_argument("--intarna-batch-size", type=int)
        command.add_argument(
            "--short-anti-repeat-rescue",
            action=argparse.BooleanOptionalAction,
            default=None,
            help="Experimental Type II short-seed rescue with anchored independent support",
        )
        command.add_argument("--search-backend", choices=["fasta36", "blastn"])
        command.add_argument("--model-rescue", action=argparse.BooleanOptionalAction, default=None)
        command.add_argument("--model-rescue-evalue", type=float)
        command.add_argument("--mask-arrays", action=argparse.BooleanOptionalAction, default=None)
        command.add_argument("--identify-preset", choices=["balanced", "region-screen", "full"])
        command.add_argument("--gene-mode", choices=["auto", "single", "meta"])
        command.add_argument("--model", dest="models", type=str, action="append")
        command.add_argument(
            "--type-v-motifs",
            "--perform_type_v_anti_repeat_analysis",
            type=boolean,
            nargs="?",
            const=True,
            default=None,
        )
        command.add_argument("--similarity", "--anti_repeat_similarity_threshold", type=float)
        command.add_argument("--coverage", "--anti_repeat_coverage_threshold", type=float)
        command.add_argument("--tail-length", type=int)
        command.add_argument("--upstream-length", type=int)
        command.add_argument("--merge-gap", type=int)
        if name == "run":
            inputs = command.add_mutually_exclusive_group(required=True)
            inputs.add_argument("--input", "--file", dest="input", type=Path)
            inputs.add_argument("--input-folder", "--input_folder", dest="input", type=Path)
            command.add_argument(
                "--output", "--output_folder", type=Path, default=Path("CRISPRtracrRNA_output")
            )
            command.add_argument("--verbose", action="store_true")
    return root


def main(argv: list[str] | None = None) -> int:
    args_list = list(sys.argv[1:] if argv is None else argv)
    # Preserve the v2 direct flag invocation while keeping explicit subcommands.
    if args_list and args_list[0].startswith("--") and args_list[0] not in {"--help", "--version"}:
        args_list.insert(0, "run")
    cli = parser()
    args = cli.parse_args(args_list)
    logging.basicConfig(
        level=logging.DEBUG if getattr(args, "verbose", False) else logging.INFO,
        format="%(levelname)s: %(message)s",
    )
    overrides = {
        key: value for key, value in vars(args).items() if key in Settings.__dataclass_fields__
    }
    if overrides.get("mode") in {"complete_run", "model_run"}:
        overrides["mode"] = overrides["mode"].removesuffix("_run")
    if overrides.get("models"):
        overrides["models"] = tuple(str(Path(p).resolve()) for p in overrides["models"])
    try:
        settings = load_settings(args.config, overrides)
        if args.command == "doctor":
            result = doctor(settings)
            print(json.dumps(result, indent=2))
            return 0 if result["ready"] else 1
        report = run(args.input, args.output, settings)
        print(
            f"Completed {report['summary']['sequence_count']} sequences; "
            f"{report['summary']['prediction_count']} candidates. "
            f"Report: {args.output.resolve() / 'report.json'}"
        )
        return 0
    except (TracrError, OSError, ValueError, TypeError) as exc:
        logging.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        logging.error("Cancelled; external processes terminated and evidence retained.")
        return 130
