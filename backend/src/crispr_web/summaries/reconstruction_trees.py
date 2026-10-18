"""Bounded reconstruction-tree and ancestral-event normalization."""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping

from .common import (
    _MAX_RECONSTRUCTION_IDENTITY_CHARACTERS,
    _MAX_RECONSTRUCTION_NAME_CHARACTERS,
    _MAX_RECONSTRUCTION_NODES,
    _MAX_RECONSTRUCTION_SPACERS,
    _MAX_RECONSTRUCTION_STATE_CELLS,
    _checkpoint,
    _json,
    _mapping,
    _read_reconstruction_bytes,
    _reconstruction_json,
    _safe_file,
)

def _public_newick(value: object) -> str:
    if not isinstance(value, str):
        raise ValueError("tree is not a string")
    tree = value.strip()
    if not tree or len(tree) > 10_000 or any(ord(character) < 32 for character in tree):
        raise ValueError("tree is empty, oversized, or contains control characters")
    return tree


def _orientation_tree_rows(
    output_root: Path,
    max_rows: int,
    checkpoint: Callable[[], None] | None,
) -> tuple[list[dict[str, str]], bool]:
    path = output_root / "crispr_evor" / "additional_data" / "orientation_trees.json"
    if _safe_file(output_root, path) is None:
        return [], False
    document = _mapping(_json(output_root, path, checkpoint=checkpoint))
    tree_maps = {
        policy: document.get(policy) for policy in ("forward", "reverse", "selected")
    }
    if any(not isinstance(tree_map, Mapping) for tree_map in tree_maps.values()):
        raise ValueError("orientation tree document is missing a policy map")
    group_names = set(tree_maps["selected"])
    if (
        set(tree_maps["forward"]) != group_names
        or set(tree_maps["reverse"]) != group_names
        or any(
            not isinstance(group, str) or not group or len(group) > 240
            for group in group_names
        )
    ):
        raise ValueError("orientation tree policy maps do not contain the same groups")
    ordered_groups = sorted(group_names)
    rows = [
        {
            "group": group,
            "forward_newick": _public_newick(tree_maps["forward"][group]),
            "reverse_newick": _public_newick(tree_maps["reverse"][group]),
            "selected_newick": _public_newick(tree_maps["selected"][group]),
        }
        for group in ordered_groups[:max_rows]
    ]
    return rows, len(ordered_groups) > len(rows)


def _reconstruction_tree_rows(
    output_root: Path,
    max_rows: int,
    checkpoint: Callable[[], None] | None,
) -> tuple[list[dict[str, str]], bool]:
    path = output_root / "spacerplacer" / "additional_data" / "dict_nwk_trees.json"
    if _safe_file(output_root, path) is None:
        return [], False
    document = _mapping(_json(output_root, path, checkpoint=checkpoint))
    group_names = set(document)
    if any(
        not isinstance(group, str) or not group or len(group) > 240
        for group in group_names
    ):
        raise ValueError("reconstruction tree document contains an invalid group")
    ordered_groups = sorted(group_names)
    rows = [
        {
            "group": group,
            "selected_newick": _public_newick(document[group]),
        }
        for group in ordered_groups[:max_rows]
    ]
    return rows, len(ordered_groups) > len(rows)


def _reconstruction_name(value: object, *, kind: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > _MAX_RECONSTRUCTION_NAME_CHARACTERS
        or any(ord(character) < 32 for character in value)
    ):
        raise ValueError(f"reconstruction {kind} name is invalid")
    return value


def _reconstruction_identity_map(document: object) -> dict[str, int]:
    if (
        not isinstance(document, Mapping)
        or not document
        or len(document) > _MAX_RECONSTRUCTION_SPACERS
    ):
        raise ValueError("reconstruction spacer map is invalid or oversized")
    identities: dict[str, int] = {}
    numbers: set[int] = set()
    for identity, number in document.items():
        if (
            not isinstance(identity, str)
            or not identity
            or len(identity) > _MAX_RECONSTRUCTION_IDENTITY_CHARACTERS
            or any(ord(character) < 32 for character in identity)
            or isinstance(number, bool)
            or not isinstance(number, int)
            or number < 1
            or number > _MAX_RECONSTRUCTION_SPACERS
            or number in numbers
        ):
            raise ValueError("reconstruction spacer map contains an invalid entry")
        identities[identity] = number
        numbers.add(number)
    if numbers != set(range(1, len(identities) + 1)):
        raise ValueError("reconstruction spacer numbering is not contiguous")
    return identities


def _reconstruction_top_order(
    document: object,
    identities: Mapping[str, int],
) -> list[str]:
    if (
        not isinstance(document, list)
        or len(document) != len(identities)
        or len(document) > _MAX_RECONSTRUCTION_SPACERS
    ):
        raise ValueError("reconstruction spacer order is invalid or oversized")
    order: list[str] = []
    seen: set[str] = set()
    for identity in document:
        if (
            not isinstance(identity, str)
            or identity not in identities
            or identity in seen
        ):
            raise ValueError("reconstruction spacer order contains an invalid entry")
        order.append(identity)
        seen.add(identity)
    if seen != set(identities):
        raise ValueError("reconstruction spacer order and numbering do not agree")
    return order


def _reconstruction_node_states(
    document: object,
    top_order: list[str],
    checkpoint: Callable[[], None] | None,
) -> list[dict[str, object]]:
    if not isinstance(document, Mapping):
        raise ValueError("reconstruction state document is not an object")
    raw_states = document.get("rec_spacers")
    if (
        not isinstance(raw_states, Mapping)
        or not raw_states
        or len(raw_states) > _MAX_RECONSTRUCTION_NODES
        or len(raw_states) * len(top_order) > _MAX_RECONSTRUCTION_STATE_CELLS
    ):
        raise ValueError("reconstruction node states are invalid or oversized")
    nodes: list[dict[str, object]] = []
    for raw_name, profile in raw_states.items():
        _checkpoint(checkpoint)
        name = _reconstruction_name(raw_name, kind="node")
        if not isinstance(profile, list) or len(profile) != len(top_order):
            raise ValueError("reconstruction node profile has the wrong width")
        present: list[str] = []
        for identity, state in zip(top_order, profile, strict=True):
            if (
                isinstance(state, bool)
                or not isinstance(state, int)
                or state not in (0, 1)
            ):
                raise ValueError("reconstruction node profile is not binary")
            if state == 1:
                present.append(identity)
        nodes.append({"name": name, "spacer_identities": present})
    return nodes


def _load_reconstruction_source(
    output_root: Path,
    directory: Path,
    *,
    group: str,
    hypothesis: str,
    budget: list[int],
    checkpoint: Callable[[], None] | None,
) -> dict[str, object]:
    prefix = directory.name
    identities = _reconstruction_identity_map(
        _reconstruction_json(
            output_root,
            directory / f"{prefix}_spacer_names_to_numbers.json",
            budget=budget,
            checkpoint=checkpoint,
        )
    )
    top_order = _reconstruction_top_order(
        _reconstruction_json(
            output_root,
            directory / f"{prefix}_top_order.json",
            budget=budget,
            checkpoint=checkpoint,
        ),
        identities,
    )
    nodes = _reconstruction_node_states(
        _reconstruction_json(
            output_root,
            directory / f"{prefix}_rec_spacers.json",
            budget=budget,
            checkpoint=checkpoint,
        ),
        top_order,
        checkpoint,
    )
    newick = _public_newick(
        _read_reconstruction_bytes(
            output_root,
            directory / f"{prefix}_tree.nwk",
            budget=budget,
            checkpoint=checkpoint,
        ).decode("utf-8")
    )
    return {
        "group": group,
        "hypothesis": hypothesis,
        "directory": directory,
        "prefix": prefix,
        "identities": identities,
        "top_order": top_order,
        "nodes": nodes,
        "newick": newick,
    }


