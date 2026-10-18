"""Discovery and bounded publication of ancestral reconstruction sources."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Mapping

from .common import (
    _MAX_RECONSTRUCTION_TOTAL_INPUT_BYTES,
    _MAX_RECONSTRUCTION_WARNINGS,
    _checkpoint,
    _safe_directory,
)
from .reconstruction_events import (
    _branch_reconstruction_events,
    _canonical_spacer_ids,
    _special_reconstruction_events,
)
from .reconstruction_trees import _load_reconstruction_source, _reconstruction_name


def _public_reconstruction_source(
    output_root: Path,
    source: Mapping[str, object],
    *,
    canonical: Mapping[str, int],
    budget: list[int],
    warn: Callable[[str], None],
    checkpoint: Callable[[], None] | None,
) -> dict[str, object]:
    gains, losses = _branch_reconstruction_events(
        output_root,
        source,
        canonical=canonical,
        budget=budget,
        warn=warn,
        checkpoint=checkpoint,
    )
    special = _special_reconstruction_events(
        output_root,
        source,
        canonical=canonical,
        budget=budget,
        warn=warn,
        checkpoint=checkpoint,
    )
    raw_nodes = source["nodes"]
    top_order = source["top_order"]
    if not isinstance(raw_nodes, list) or not isinstance(top_order, list):
        raise ValueError("reconstruction source is incomplete")
    nodes: list[dict[str, object]] = []
    for raw_node in raw_nodes:
        _checkpoint(checkpoint)
        if not isinstance(raw_node, Mapping):
            raise ValueError("reconstruction node is invalid")
        name = raw_node.get("name")
        spacer_identities = raw_node.get("spacer_identities")
        if not isinstance(name, str) or not isinstance(spacer_identities, list):
            raise ValueError("reconstruction node is incomplete")
        nodes.append(
            {
                "name": name,
                "spacers": [canonical[identity] for identity in spacer_identities],
                "gains": gains.get(name, []),
                "loss_blocks": losses.get(name, []),
                **{
                    public_name: event_map.get(name, [])
                    for public_name, event_map in special.items()
                },
            }
        )
    return {
        "group": source["group"],
        "hypothesis": source["hypothesis"],
        "newick": source["newick"],
        "spacer_order": [canonical[identity] for identity in top_order],
        "nodes": nodes,
        "acquisition_count": sum(len(node["gains"]) for node in nodes),
        "deletion_count": sum(
            len(block) for node in nodes for block in node["loss_blocks"]
        ),
    }


def _discover_reconstruction_directories(
    output_root: Path,
    *,
    max_files: int,
    checkpoint: Callable[[], None] | None,
) -> tuple[dict[str, dict[str, Path]], bool]:
    found: dict[str, dict[str, Path]] = {}
    truncated = False
    for folder, hypothesis in (("0_forward", "input"), ("0_reversed", "reverse")):
        events_root = _safe_directory(
            output_root,
            output_root / "crispr_evor" / folder / "reconstructed_events",
        )
        if events_root is None:
            continue
        candidates: list[tuple[str, Path]] = []
        scanned = 0
        try:
            entries = events_root.iterdir()
            for entry in entries:
                _checkpoint(checkpoint)
                scanned += 1
                if scanned > max_files:
                    truncated = True
                    break
                name = entry.name
                if not name.endswith("_rec") or name.endswith("_guide_rec"):
                    continue
                group = name[: -len("_rec")]
                try:
                    group = _reconstruction_name(group, kind="group")
                except ValueError:
                    continue
                safe = _safe_directory(output_root, entry)
                if safe is not None:
                    candidates.append((group, safe))
        except OSError:
            continue
        for group, directory in sorted(candidates, key=lambda item: item[0]):
            found.setdefault(group, {})[hypothesis] = directory
    return found, truncated


def _orientation_reconstruction_rows(
    output_root: Path,
    max_rows: int,
    max_files: int,
    warnings: list[str],
    checkpoint: Callable[[], None] | None,
) -> tuple[list[dict[str, object]], bool]:
    directories, discovery_truncated = _discover_reconstruction_directories(
        output_root,
        max_files=max_files,
        checkpoint=checkpoint,
    )
    groups = sorted(directories)
    selected_groups = groups[:max_rows]
    budget = [_MAX_RECONSTRUCTION_TOTAL_INPUT_BYTES]
    sources_by_group: dict[str, list[dict[str, object]]] = {}
    local_warning_count = 0
    suppressed_warning_count = 0

    def add_warning(message: str) -> None:
        nonlocal local_warning_count, suppressed_warning_count
        if local_warning_count < _MAX_RECONSTRUCTION_WARNINGS:
            warnings.append(message)
            local_warning_count += 1
        else:
            suppressed_warning_count += 1

    for group in selected_groups:
        for hypothesis in ("input", "reverse"):
            directory = directories[group].get(hypothesis)
            if directory is None:
                continue
            try:
                source = _load_reconstruction_source(
                    output_root,
                    directory,
                    group=group,
                    hypothesis=hypothesis,
                    budget=budget,
                    checkpoint=checkpoint,
                )
            except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
                add_warning(
                    "Could not summarize CRISPR-evOr reconstruction for "
                    f"{group} ({hypothesis}): {type(error).__name__}"
                )
                continue
            sources_by_group.setdefault(group, []).append(source)

    rows: list[dict[str, object]] = []
    for group in selected_groups:
        sources = sources_by_group.get(group, [])
        if not sources:
            continue
        try:
            canonical = _canonical_spacer_ids(sources)
        except ValueError as error:
            add_warning(
                "Could not canonicalize CRISPR-evOr reconstruction for "
                f"{group}: {type(error).__name__}"
            )
            continue
        identity_sets = [set(source["identities"]) for source in sources]
        if len(identity_sets) == 2 and identity_sets[0] != identity_sets[1]:
            add_warning(
                "CRISPR-evOr reconstruction hypotheses contain different spacer "
                f"identity sets for {group}; canonical IDs cover their union."
            )
        for source in sorted(
            sources,
            key=lambda item: 0 if item["hypothesis"] == "input" else 1,
        ):
            hypothesis = str(source["hypothesis"])

            def source_warning(detail: str) -> None:
                add_warning(
                    f"CRISPR-evOr reconstruction for {group} ({hypothesis}): {detail}."
                )

            try:
                rows.append(
                    _public_reconstruction_source(
                        output_root,
                        source,
                        canonical=canonical,
                        budget=budget,
                        warn=source_warning,
                        checkpoint=checkpoint,
                    )
                )
            except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
                add_warning(
                    "Could not normalize CRISPR-evOr reconstruction for "
                    f"{group} ({hypothesis}): {type(error).__name__}"
                )
    if suppressed_warning_count:
        warnings.append(
            f"{suppressed_warning_count} additional CRISPR-evOr reconstruction "
            "warning(s) were omitted."
        )
    return rows, discovery_truncated or len(groups) > len(selected_groups)


