"""Frozen, user-facing runtime presets and their provenance contract."""

import hashlib
import json


PRESET_SCHEMA_VERSION = 1

# These are the settings whose scientific or output semantics are controlled
# by a named preset.  Resource controls such as CPU count and folder lanes are
# deliberately excluded so operators can tune them for their scheduler.
BALANCED_DEFAULTS = {
    "output_profile": "standard",
    "accepted_only": False,
    "strand": True,
    "is_element": False,
    "cas": False,
    "fast_run": True,
    "fast_run_seed_profile": "balanced",
    "fast_run_rescue": True,
    "fast_run_phase_repair": True,
    "fast_two_repeat_terminal_repair": False,
    "fast_accepted_family_rescue": False,
    "split_competitive_components": True,
    "degenerated": True,
    "architecture_rescue": False,
    "reconstruction_policy": "conservative",
    "fuzzy_search_mode": "native",
    "fuzzy_search_flank": 100,
    "fuzzy_family_padding": 15,
    "fuzzy_boundary_paths": True,
}

REGION_SCREEN_CONFIG = {
    "output_profile": "machine",
    "accepted_only": True,
    "strand": False,
    "is_element": False,
    "cas": False,
    "fast_run": True,
    "fast_run_seed_profile": "balanced",
    "fast_run_rescue": False,
    "fast_run_phase_repair": False,
    "fast_two_repeat_terminal_repair": False,
    "fast_accepted_family_rescue": False,
    "split_competitive_components": False,
    "degenerated": False,
    "architecture_rescue": False,
    "reconstruction_policy": "legacy",
    "fuzzy_search_mode": "native",
    "fuzzy_search_flank": 100,
    "fuzzy_family_padding": 15,
    "fuzzy_boundary_paths": False,
}

# The historical full-array search is selected by fast_run=False.  Optional
# downstream annotations remain opt-in: enabling Cas or IS annotation would
# change the requested analysis rather than the array-search preset.
FULL_CONFIG = {
    destination: value
    for destination, value in dict(
        BALANCED_DEFAULTS,
        fast_run=False,
        fast_run_rescue=False,
        fast_run_phase_repair=False,
        split_competitive_components=False,
        reconstruction_policy="legacy",
    ).items()
    if destination not in ("cas", "is_element")
}

PRESET_CONFIGS = {
    "balanced": BALANCED_DEFAULTS,
    "region-screen": REGION_SCREEN_CONFIG,
    "full": FULL_CONFIG,
}

PRESET_CONTRACTS = {
    "balanced": (
        "Default high-recall fast analysis with core-supported reconstruction. "
        "Direct scientific flags may customize this preset."
    ),
    "region-screen": (
        "Validated high-throughput biological-region recall screen. It may "
        "change candidate representations, categories, and exact borders."
    ),
    "full": (
        "Historical full array search with standard output; optional Cas and "
        "IS-element annotations remain independently opt-in."
    ),
}


class PresetConflictError(ValueError):
    """An explicit scientific option contradicts an atomic preset."""


def explicit_option_destinations(parser, argv):
    """Return argparse destinations explicitly present in an argument vector.

    Both ``--option value`` and ``--option=value`` forms are recognized.  The
    parser remains the authority for aliases and unknown-option errors.
    """

    destinations = set()
    for token in argv:
        if token == "--":
            break
        if not isinstance(token, str) or not token.startswith("--"):
            continue
        option = token.split("=", 1)[0]
        action = parser._option_string_actions.get(option)
        if action is not None:
            destinations.add(action.dest)
    return frozenset(destinations)


def resolve_preset(preset, values, explicit_destinations=()):
    """Resolve one preset without silently overriding explicit conflicts.

    ``balanced`` preserves the historical direct-flag behavior.  The two
    non-default presets are atomic: a controlled option may be repeated only
    when its parsed value agrees with the preset.
    """

    if preset not in PRESET_CONFIGS:
        raise ValueError("unknown runtime preset: {}".format(preset))
    missing = sorted(set(BALANCED_DEFAULTS).difference(values))
    if missing:
        raise ValueError(
            "runtime configuration is missing: {}".format(
                ", ".join(missing)
            )
        )

    resolved = dict(values)
    explicit_destinations = set(explicit_destinations)
    if preset == "balanced":
        # Architecture rescue is an independently selected workflow. Preserve
        # it when the reconstruction option is omitted, but never silently
        # discard an explicit request for both incompatible policies.
        if values["architecture_rescue"]:
            if ("reconstruction_policy" in explicit_destinations
                    and values["reconstruction_policy"] != "legacy"):
                raise PresetConflictError(
                    "conservative reconstruction requires --preset balanced "
                    "and --architecture_rescue False"
                )
            resolved["reconstruction_policy"] = "legacy"
        return resolved

    expected = PRESET_CONFIGS[preset]
    conflicts = []
    for destination, required in expected.items():
        observed = values[destination]
        if destination in explicit_destinations and observed != required:
            conflicts.append((destination, observed, required))
    if conflicts:
        details = "; ".join(
            "--{}={} (requires {})".format(
                destination, json.dumps(observed), json.dumps(required)
            )
            for destination, observed, required in conflicts
        )
        raise PresetConflictError(
            "explicit option conflicts with --preset {}: {}".format(
                preset, details
            )
        )

    resolved.update(expected)
    return resolved


def preset_provenance(
        preset, resolved_config, explicit_destinations=(), preset_explicit=False):
    """Return deterministic provenance for a resolved runtime configuration."""

    canonical = json.dumps(
        resolved_config,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    controlled_explicit = sorted(
        set(explicit_destinations).intersection(BALANCED_DEFAULTS)
    )
    return {
        "schema_version": PRESET_SCHEMA_VERSION,
        "name": preset,
        "selection": "explicit" if preset_explicit else "default",
        "contract": PRESET_CONTRACTS[preset],
        "explicit_controlled_options": controlled_explicit,
        "effective_config_sha256": hashlib.sha256(canonical).hexdigest(),
    }
