"""Command-line interface for CRISPR tool integration."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

from .adapter import adapt_reports
from .annotations import run_annotation_pipeline
from .runner import load_runner_config, run_pipeline
from .scientific_options import load_tool_options_environment


def _positive_int(value: str) -> int:
    converted = int(value)
    if converted < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return converted


def _non_negative_int(value: str) -> int:
    converted = int(value)
    if converted < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return converted


def _positive_float(value: str) -> float:
    converted = float(value)
    if converted <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return converted


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="crispr-tools",
        description="CRISPRloci v2: arrays, Cas systems, tracrRNA and CRISPR evolution.",
    )
    parser.add_argument("--version", action="version", version="%(prog)s 0.14.0")
    subparsers = parser.add_subparsers(dest="action", required=True)

    adapt = subparsers.add_parser("adapt", help="convert existing CRISPRidentify JSON output to CCDB input")
    adapt.add_argument("input", help="report.json, Complete_json_summary.json, or CRISPRidentify result folder")
    adapt.add_argument("output", help="new or empty adapter output directory")
    adapt.add_argument(
        "--source-fasta",
        type=Path,
        help="original FASTA file/folder; require report ID, length, and sequence hash to match",
    )
    adapt.add_argument(
        "--categories", nargs="+", default=["Bona-fide", "Possible"],
        help="CRISPRidentify categories to retain (default: Bona-fide Possible)",
    )
    adapt.add_argument("--minimum-arrays", type=_positive_int, default=2)

    run = subparsers.add_parser("run", help="run the connected three-tool workflow")
    run.add_argument("input", help="input FASTA/folder or existing CRISPRidentify result path")
    run.add_argument("output", help="new or empty job output directory")
    run.add_argument(
        "--input-kind", choices=["fasta", "identify-results"], default="fasta",
        help="reuse existing CRISPRidentify output with identify-results",
    )
    run.add_argument("--runner-config", type=Path, help="JSON arrays defining isolated worker commands")
    run.add_argument(
        "--tools", choices=["identify", "both", "spacerplacer", "crispr-evor", "casandra", "tracrrna", "loci", "leader", "repeats", "repeat-context", "protospacer", "viral-search", "repeat-type", "repeat-map", "spacer-association", "array-compare", "locus-compare", "prophage-context", "insertion-context", "defense-context", "rna-family"], default="both",
        help=(
            "analysis selection: identify stops after detection/normalization; both uses CRISPR-evOr's "
            "orientation-selected SpacerPlacer reconstruction (default: both)"
        ),
    )
    run.add_argument("--array-map", type=Path)
    run.add_argument("--association-reference", type=Path)
    run.add_argument("--association-reference-sha256")
    run.add_argument("--association-groups", type=Path)
    run.add_argument("--repeat-reference", type=Path)
    run.add_argument("--repeat-reference-sha256")
    run.add_argument("--viral-reference", type=Path)
    run.add_argument("--viral-reference-sha256")
    run.add_argument("--viral-max-mismatches", type=int, choices=range(3), default=2)
    run.add_argument("--spacer-database", type=Path, help="Operator-installed immutable spacer collection")
    run.add_argument("--spacer-database-sha256", help="Operator-pinned SHA-256 of SQLite database bytes")
    run.add_argument("--leader-flank-length", type=_positive_int, default=500)
    run.add_argument("--molecule", choices=["DNA", "RNA"], default="DNA")
    run.add_argument("--tracr-model-type", choices=["II", "V"], default="II")
    run.add_argument("--categories", nargs="+", default=["Bona-fide", "Possible"])
    run.add_argument("--minimum-arrays", type=_positive_int, default=2)
    run.add_argument("--spacer-distance", type=_non_negative_int, default=1)
    run.add_argument("--cpu", type=_positive_int)
    run.add_argument(
        "--identify-folder-jobs", type=_positive_int,
        help=(
            "CRISPRidentify lanes for independent FASTA records; the worker clips the effective "
            "count to the input-record and --cpu limits (default: worker default of 1)"
        ),
    )
    run.add_argument("--stage-timeout", type=_positive_float, default=3600.0)
    run.add_argument(
        "--disable-bias-corrections", action="store_true",
        help="disable optional alpha/rho bias corrections if they are numerically undefined for the dataset",
    )
    run.add_argument(
        "--rho-bias-failure-policy",
        choices=["fallback_uncorrected", "fail"],
        default="fallback_uncorrected",
        help=(
            "when a requested rho correction has no finite fit, retry only that group/orientation without it "
            "and record a warning, or fail strictly (default: fallback_uncorrected)"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config = None
        if args.action == "run":
            mode = {"identify": "detection", "both": "orientation", "spacerplacer": "reconstruction", "crispr-evor": "orientation", "casandra": "cas", "repeat-context": "repeat_context", "viral-search": "viral_search", "repeat-map": "repeat_map", "repeat-type": "repeat_type", "spacer-association": "spacer_association", "array-compare": "array_compare", "locus-compare": "locus_compare", "prophage-context": "prophage_context", "insertion-context": "insertion_context", "defense-context": "defense_context", "rna-family": "rna_family"}.get(args.tools, args.tools)
            config = replace(load_runner_config(args.runner_config), tool_options=load_tool_options_environment(mode))
        if args.action == "adapt":
            result = adapt_reports(
                args.input,
                args.output,
                categories=args.categories,
                minimum_arrays=args.minimum_arrays,
                expected_fasta_path=args.source_fasta,
            )
            display_result = {
                "status": "completed",
                "manifest": str(Path(args.output).expanduser().resolve() / "manifest.json"),
                "input_array_count": result["input_array_count"],
                "selected_array_count": result["selected_array_count"],
                "emitted_array_count": result["emitted_array_count"],
                "emitted_group_count": result["emitted_group_count"],
                "skipped_array_count": len(result["skipped_arrays"]),
            }
        elif args.tools == "locus-compare":
            from .locus_compare import run_locus_compare_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Locus compare requires supplied genomic DNA")
            display_result = run_locus_compare_pipeline(
                args.input, args.output, runner_config=config, cpu=args.cpu or 1,
                timeout_seconds=args.stage_timeout)
        elif args.tools == "prophage-context":
            from .prophage_context import run_prophage_context_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Prophage context requires supplied genomic DNA")
            display_result = run_prophage_context_pipeline(
                args.input, args.output, runner_config=config, cpu=args.cpu or 1,
                timeout_seconds=args.stage_timeout)
        elif args.tools == "insertion-context":
            from .insertion_context import run_insertion_context_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Insertion context requires supplied genomic DNA")
            display_result = run_insertion_context_pipeline(
                args.input, args.output, runner_config=config, cpu=args.cpu or 1,
                timeout_seconds=args.stage_timeout)
        elif args.tools == "defense-context":
            from .defense_context import run_defense_context_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Defense context requires supplied genomic DNA")
            display_result = run_defense_context_pipeline(
                args.input, args.output, runner_config=config, cpu=args.cpu or 1,
                timeout_seconds=args.stage_timeout)
        elif args.tools == "rna-family":
            from .rna_family import run_rna_family_pipeline
            if args.input_kind != "fasta":
                raise ValueError("Rna family requires supplied sequences")
            display_result = run_rna_family_pipeline(
                args.input, args.output, runner_config=config, cpu=args.cpu or 1,
                timeout_seconds=args.stage_timeout, molecule=args.molecule)
        elif args.tools == "array-compare":
            from .array_compare import run_array_compare_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA" or args.array_map is None:
                raise ValueError("Array comparison requires observed DNA and a structured array map")
            display_result = run_array_compare_pipeline(args.input,args.output,runner_config=config,map_path=args.array_map,timeout_seconds=args.stage_timeout)
        elif args.tools == "spacer-association":
            from .association import run_association_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA" or args.association_groups is None or args.association_reference is None:
                raise ValueError("Spacer association requires DNA spacers, explicit groups and an operator panel")
            display_result = run_association_pipeline(args.input,args.output,runner_config=config,
                reference_path=args.association_reference,reference_sha256=args.association_reference_sha256,
                groups_path=args.association_groups,cpu=args.cpu or 1,timeout_seconds=args.stage_timeout)
        elif args.tools == "repeat-type":
            from .repeat_type import run_repeat_type_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Repeat typing requires DNA repeat FASTA")
            display_result = run_repeat_type_pipeline(args.input, args.output, runner_config=config, timeout_seconds=args.stage_timeout)
        elif args.tools == "repeat-map":
            from .repeat_map import run_repeat_map_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Repeat mapping requires DNA repeat FASTA")
            display_result = run_repeat_map_pipeline(args.input, args.output,
                runner_config=config,
                reference_path=args.repeat_reference, reference_sha256=args.repeat_reference_sha256,
                timeout_seconds=args.stage_timeout)
        elif args.tools == "viral-search":
            from .viral import run_viral_pipeline
            if args.input_kind != "fasta":
                raise ValueError("Viral search requires supplied spacer FASTA")
            display_result = run_viral_pipeline(args.input, args.output, runner_config=config,
                reference_path=args.viral_reference, reference_sha256=args.viral_reference_sha256,
                max_mismatches=args.viral_max_mismatches, threads=args.cpu or 2, timeout_seconds=args.stage_timeout)
        elif args.tools == "protospacer":
            from .protospacer import run_protospacer_pipeline
            if args.input_kind != "fasta" or args.molecule != "DNA":
                raise ValueError("Protospacer search requires original DNA FASTA")
            display_result = run_protospacer_pipeline(
                args.input, args.output, runner_config=config,
                database_path=args.spacer_database, database_sha256=args.spacer_database_sha256,
                timeout_seconds=args.stage_timeout,
            )
        elif args.tools in {"repeats", "repeat-context"}:
            from .repeats import run_repeat_pipeline
            if args.input_kind != "fasta":
                raise ValueError("Repeat workflows require original FASTA input")
            display_result = run_repeat_pipeline(
                args.input, args.output, runner_config=config,
                context=args.tools == "repeat-context", molecule=args.molecule,
                categories=args.categories, cpu=args.cpu, timeout_seconds=args.stage_timeout,
            )
        elif args.tools == "leader":
            from .leader import run_leader_pipeline
            if args.input_kind != "fasta":
                raise ValueError("Leader context requires FASTA input")
            display_result = run_leader_pipeline(
                args.input, args.output, runner_config=config,
                cpu=args.cpu, timeout_seconds=args.stage_timeout, flank_length=args.leader_flank_length,
            )
        elif args.tools in {"casandra", "tracrrna", "loci"}:
            if args.input_kind != "fasta":
                raise ValueError("Annotation workflows require FASTA input")
            display_result = run_annotation_pipeline(
                args.input, args.output, runner_config=config,
                selection=args.tools, cpu=args.cpu, timeout_seconds=args.stage_timeout,
                model_type=args.tracr_model_type, flank_length=args.leader_flank_length,
            )
        else:
            tools = {
                "identify": (),
                "both": ("spacerplacer", "crispr_evor"),
                "spacerplacer": ("spacerplacer",),
                "crispr-evor": ("crispr_evor",),
            }[args.tools]
            result = run_pipeline(
                args.input,
                args.output,
                input_kind=args.input_kind,
                runner_config=config,
                tools=tools,
                categories=args.categories,
                minimum_arrays=args.minimum_arrays,
                spacer_distance=args.spacer_distance,
                cpu=args.cpu,
                identify_folder_jobs=args.identify_folder_jobs,
                timeout_seconds=args.stage_timeout,
                bias_corrections=not args.disable_bias_corrections,
                rho_bias_failure_policy=args.rho_bias_failure_policy,
            )
            display_result = result
    except (OSError, ValueError, RuntimeError) as error:
        print(f"crispr-tools: error: {error}", file=sys.stderr)
        return 1
    print(json.dumps(display_result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
