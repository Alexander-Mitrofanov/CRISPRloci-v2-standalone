"""Shared standalone resource and scientific input preconditions."""

from __future__ import annotations
import argparse
import math
from pathlib import Path
from .errors import CCTKError
from .io import safe_open, read_arrays


def plot_limits(width: float, height: float, dpi: float, path: str) -> None:
    if (
        any(not math.isfinite(x) or x <= 0 for x in (width, height, dpi))
        or width * height * dpi * dpi > 20_000_000
        or max(width, height) * dpi > 12000
    ):
        raise CCTKError("plot exceeds finite dimension or 20-million-pixel limit")
    target = Path(path)
    if target.suffix.lower() not in {".png", ".pdf", ".svg"} or any(
        p.is_symlink() for p in (target, *target.parents)
    ):
        raise CCTKError("plots require a regular PNG, PDF or SVG output path")


def validate_command(args: argparse.Namespace, command: str) -> None:
    """Apply to direct Python command APIs as well as the installed CLI."""
    for name, value in vars(args).items():
        if isinstance(value, float) and not math.isfinite(value):
            raise CCTKError(f"{name} must be finite")
    for name, maximum in [
        ("iterations", 1000),
        ("replicates", 1000),
        ("num_threads", 64),
        ("num_events", 10000),
        ("initial_length", 500),
        ("max_cluster", 200),
        ("min_shared", 500),
        ("batch_size", 10000),
        ("max_target_seqs", 100000),
    ]:
        value = getattr(args, name, None)
        if value is not None and (type(value) is not int or not 1 <= value <= maximum):
            raise CCTKError(f"{name} is outside supported limits")
    if hasattr(args, "seed") and args.seed is None:
        args.seed = 2
    for name in [
        "acquisition",
        "deletion",
        "insertion",
        "rep_indel",
        "duplication",
        "trailer_loss",
        "no_ident",
    ]:
        value = getattr(args, name, None)
        if value is not None and (not 0 <= value <= 100000):
            raise CCTKError("event cost or frequency outside supported limits")
    if command == "evolve":
        if args.acquisition == 0 and args.trailer_loss == 0:
            raise CCTKError("deletion-only simulation can reach an impossible one-spacer state")
        if (
            not 0 <= args.loss_rate <= 100
            or args.acquisition + args.deletion + args.trailer_loss <= 0
        ):
            raise CCTKError("invalid simulation frequencies")
    if command == "quickrun" and (not args.indir or not Path(args.indir).is_dir()):
        raise CCTKError("quickrun requires a genome directory")
    if hasattr(args, "plot_width"):
        path = getattr(args, "out_file", None) or getattr(args, "plot_file", None) or "bounded.png"
        plot_limits(args.plot_width, args.plot_height, args.dpi, path)
    for name in (
        "array_file",
        "input_tree",
        "genome_array_file",
        "types",
        "repeats_file",
        "spacer_file",
        "colour_file",
        "colour_scheme_infile",
    ):
        path = getattr(args, name, None)
        if path:
            with safe_open(path):
                pass
    if hasattr(args, "array_file") and args.array_file:
        arrays = read_arrays(args.array_file)
        selected = (
            getattr(args, "arrays_to_align", None)
            or getattr(args, "arrays_to_join", None)
            or list(arrays)
        )
        if len(set(selected)) != len(selected) or any(x not in arrays for x in selected):
            raise CCTKError("selected arrays must be unique input identifiers")
        if command in {"crisprtree", "constrain"} and len(arrays) > 50:
            raise CCTKError("parsimony methods support at most 50 arrays per run")
