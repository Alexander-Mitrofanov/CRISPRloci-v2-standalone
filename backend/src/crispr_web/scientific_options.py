"""Bounded overrides for flags shipped by the pinned scientific CLIs.

This allowlist is deliberately duplicated at the API and isolated Integration
boundary. Neither accepts arbitrary argv, executable names, or filesystem paths.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass

MAX_TOOL_OPTIONS_BYTES = 8192
TOOL_OPTIONS_ENV = "CRISPRLOCI_TOOL_OPTIONS"


@dataclass(frozen=True)
class Option:
    kind: str
    default: object
    minimum: float | None = None
    maximum: float | None = None
    choices: tuple = ()
    modes: tuple[str, ...] = ()
    boolean_style: str = "value"


IDENTIFY_MODES = ("detection", "orientation", "reconstruction", "leader", "repeat_context")
TOOL_MODES = {
    "cctk": ("array_compare",),
    "crispridentify": IDENTIFY_MODES,
    "spacerplacer": ("reconstruction",),
    "crispr_evor": ("orientation",),
    "casandra": ("cas",),
    "crispr_tracrrna": ("tracrrna", "loci"),
    "crisprrepeat": ("repeats", "repeat_context"),
    "crisprmap": ("repeat_map",),
    "repeattyper": ("repeat_type",),
    "spacepharer": ("spacer_association",),
}
_EVOLUTION = {
    "--seed": Option("integer", 2357, 0, 2147483647),
    "--tree_construction_method": Option("select", "upgma", choices=("upgma", "nj")),
    "--insertion_rate": Option("number", 1e-5, 1e-10, 100),
    "--deletion_rate": Option("number", 0.1, 1e-10, 100),
    "--extend_branches": Option("number", 0, 0, 100),
    "--significance_level": Option("number", 0.05, 1e-6, 0.999999),
}
OPTIONS = {
    "crispridentify": {
        "--fast_run": Option("boolean", True),
        "--fast_run_seed_profile": Option("select", "balanced", choices=("balanced", "legacy", "sensitive")),
        "--degenerated": Option("boolean", True),
        "--min_len_rep": Option("integer", 21, 10, 100),
        "--max_len_rep": Option("integer", 55, 10, 100),
        "--min_len_spacer": Option("integer", 18, 1, 500),
        "--max_len_spacer": Option("integer", 78, 1, 500),
        "--min_repeats": Option("integer", 3, 2, 20),
    },
    "spacerplacer": dict(_EVOLUTION),
    "crispr_evor": {
        **_EVOLUTION,
        "--orientation_decision_boundary": Option("number", 5, 0, 1000),
    },
    "casandra": {
        "--gene-mode": Option("select", "auto", choices=("auto", "single", "meta")),
    },
    "crispr_tracrrna": {
        "--merge-gap": Option("integer", 150, 0, 1000, modes=("tracrrna",)),
        "--identify-preset": Option("select", "balanced", choices=("balanced", "region-screen", "full"), modes=("loci",)),
        "--gene-mode": Option("select", "auto", choices=("auto", "single", "meta"), modes=("loci",)),
        "--model-rescue": Option("boolean", True, modes=("loci",), boolean_style="optional"),
        "--model-rescue-evalue": Option("number", 1e-5, 1e-20, 1, modes=("loci",)),
        "--similarity": Option("number", 0.7, 0, 1, modes=("loci",)),
        "--coverage": Option("number", 0.6, 0, 1, modes=("loci",)),
        "--tail-length": Option("integer", 150, 1, 1000, modes=("loci",)),
        "--upstream-length": Option("integer", 5, 0, 100, modes=("loci",)),
    },
    "crisprrepeat": {
        "--temperature": Option("number", 37, 0, 100),
        "--window": Option("integer", 150, 10, 300, modes=("repeat_context",)),
        "--span": Option("integer", 80, 5, 300, modes=("repeat_context",)),
        "--no-lonely-pairs": Option("boolean", False, boolean_style="switch"),
    },
    "cctk": {
        "--identity": Option("select", "exact", choices=("exact", "reverse-complement")),
        "--min-shared": Option("integer", 1, 1, 500),
    },
    "spacepharer": {"--fdr": Option("number", 0.05, 0.001, 0.5)},
    "repeattyper": {"--threshold": Option("number", 0.75, 0, 1)},
    "crisprmap": {
        "--max-distance": Option("integer", 3, 0, 3),
        "--top-k": Option("integer", 20, 1, 100),
    },
}


def validate_tool_options(value: object, mode: str | None = None) -> dict:
    if not isinstance(value, dict) or len(value) > len(OPTIONS):
        raise ValueError("tool_options must be an object of supported tools")
    try:
        serialized = json.dumps(value, allow_nan=False, separators=(",", ":"))
    except (TypeError, ValueError) as error:
        raise ValueError("tool_options must contain finite JSON scalar values") from error
    if len(serialized.encode("utf-8")) > MAX_TOOL_OPTIONS_BYTES:
        raise ValueError("tool_options exceeds its byte limit")
    validated = {}
    for tool, overrides in value.items():
        if tool not in OPTIONS or not isinstance(overrides, dict):
            raise ValueError(f"Unsupported scientific tool options: {tool}")
        if mode is not None and mode not in TOOL_MODES[tool]:
            raise ValueError(f"{tool} options do not apply to {mode}")
        values = {}
        for flag, scalar in overrides.items():
            spec = OPTIONS[tool].get(flag)
            if spec is None:
                raise ValueError(f"Unsupported {tool} flag: {flag}")
            if mode is not None and spec.modes and mode not in spec.modes:
                raise ValueError(f"{flag} does not apply to {mode}")
            if spec.kind == "boolean":
                valid = type(scalar) is bool
            elif spec.kind == "integer":
                valid = type(scalar) is int
            elif spec.kind == "number":
                valid = type(scalar) is int or (type(scalar) is float and math.isfinite(scalar))
            else:
                valid = type(scalar) is str and scalar in spec.choices
            if not valid:
                raise ValueError(f"Invalid value for {tool} {flag}")
            if spec.minimum is not None and not spec.minimum <= scalar <= spec.maximum:
                raise ValueError(f"{tool} {flag} must be between {spec.minimum} and {spec.maximum}")
            values[flag] = scalar
        if tool == "crispridentify":
            for lower, upper in (("--min_len_rep", "--max_len_rep"), ("--min_len_spacer", "--max_len_spacer")):
                if values.get(lower, OPTIONS[tool][lower].default) > values.get(upper, OPTIONS[tool][upper].default):
                    raise ValueError(f"{lower} must not exceed {upper}")
            if values.get("--fast_run") is False and "--fast_run_seed_profile" in values:
                raise ValueError("--fast_run_seed_profile requires --fast_run True")
        if tool == "crisprrepeat" and values.get("--span", 80) > values.get("--window", 150):
            raise ValueError("--span must not exceed --window")
        if tool == "crispr_tracrrna" and values.get("--model-rescue") is False and "--model-rescue-evalue" in values:
            raise ValueError("--model-rescue-evalue requires model rescue")
        if values:
            validated[tool] = values
    return validated


def option_arguments(options: dict, tool: str) -> list[str]:
    """Serialize trusted, validated overrides using each real CLI's bool syntax."""
    validated = validate_tool_options(options)
    arguments = []
    for flag, value in validated.get(tool, {}).items():
        spec = OPTIONS[tool][flag]
        if spec.kind == "boolean" and spec.boolean_style == "optional":
            arguments.append(flag if value else "--no-" + flag[2:])
        elif spec.kind == "boolean" and spec.boolean_style == "switch":
            if value:
                arguments.append(flag)
        else:
            arguments.extend((flag, str(value)))
    return arguments
