"""Canonical spacer identities and reconstructed-event normalization."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Mapping

from .common import (
    _MAX_RECONSTRUCTION_EVENT_REFERENCES,
    _MAX_RECONSTRUCTION_NODES,
    _MAX_RECONSTRUCTION_SPACERS,
    _RECONSTRUCTION_EVENT_FIELDS,
    _reconstruction_json,
    _safe_file,
)
from .reconstruction_trees import _reconstruction_name


def _canonical_spacer_ids(sources: list[dict[str, object]]) -> dict[str, int]:
    canonical: dict[str, int] = {}
    ordered_sources = sorted(
        sources,
        key=lambda source: 0 if source["hypothesis"] == "input" else 1,
    )
    for source in ordered_sources:
        identities = source["identities"]
        if not isinstance(identities, Mapping):
            raise ValueError("reconstruction spacer map was not retained")
        for identity, _number in sorted(identities.items(), key=lambda item: item[1]):
            if not isinstance(identity, str):
                raise ValueError("reconstruction spacer identity is invalid")
            if identity not in canonical:
                canonical[identity] = len(canonical) + 1
    if len(canonical) > _MAX_RECONSTRUCTION_SPACERS:
        raise ValueError("combined reconstruction spacer set is oversized")
    return canonical


def _canonical_event_ids(
    value: object,
    *,
    local_to_identity: Mapping[int, str],
    canonical: Mapping[str, int],
    reference_count: list[int],
) -> list[int]:
    if not isinstance(value, list) or len(value) > _MAX_RECONSTRUCTION_SPACERS:
        raise ValueError("reconstruction event list is invalid or oversized")
    reference_count[0] += len(value)
    if reference_count[0] > _MAX_RECONSTRUCTION_EVENT_REFERENCES:
        raise ValueError("reconstruction event data are oversized")
    result: list[int] = []
    for local_number in value:
        if (
            isinstance(local_number, bool)
            or not isinstance(local_number, int)
            or local_number not in local_to_identity
        ):
            raise ValueError("reconstruction event contains an invalid spacer number")
        result.append(canonical[local_to_identity[local_number]])
    if len(result) != len(set(result)):
        raise ValueError("reconstruction event contains duplicate spacer numbers")
    return result


def _canonical_event_map(
    value: object,
    *,
    node_names: frozenset[str],
    local_to_identity: Mapping[int, str],
    canonical: Mapping[str, int],
    reference_count: list[int],
) -> dict[str, list[int]]:
    if not isinstance(value, Mapping) or len(value) > _MAX_RECONSTRUCTION_NODES:
        raise ValueError("reconstruction event map is invalid or oversized")
    result: dict[str, list[int]] = {}
    for raw_name, raw_ids in value.items():
        name = _reconstruction_name(raw_name, kind="event node")
        if name not in node_names:
            raise ValueError("reconstruction event names an unknown node")
        result[name] = _canonical_event_ids(
            raw_ids,
            local_to_identity=local_to_identity,
            canonical=canonical,
            reference_count=reference_count,
        )
    return result


def _canonical_loss_blocks(
    value: object,
    *,
    local_to_identity: Mapping[int, str],
    canonical: Mapping[str, int],
    reference_count: list[int],
) -> list[list[int]]:
    if not isinstance(value, list) or len(value) > _MAX_RECONSTRUCTION_SPACERS:
        raise ValueError("reconstruction loss blocks are invalid or oversized")
    if not value:
        return []
    if all(isinstance(item, int) and not isinstance(item, bool) for item in value):
        losses = _canonical_event_ids(
            value,
            local_to_identity=local_to_identity,
            canonical=canonical,
            reference_count=reference_count,
        )
        return [[loss] for loss in losses]
    if not all(isinstance(item, list) for item in value):
        raise ValueError("reconstruction loss blocks mix incompatible values")
    blocks = [
        _canonical_event_ids(
            block,
            local_to_identity=local_to_identity,
            canonical=canonical,
            reference_count=reference_count,
        )
        for block in value
    ]
    if any(not block for block in blocks):
        raise ValueError("reconstruction loss blocks contain an empty block")
    flattened = [spacer for block in blocks for spacer in block]
    if len(flattened) != len(set(flattened)):
        raise ValueError("reconstruction loss blocks overlap")
    return blocks


def _canonical_loss_map(
    value: object,
    *,
    node_names: frozenset[str],
    local_to_identity: Mapping[int, str],
    canonical: Mapping[str, int],
    reference_count: list[int],
) -> dict[str, list[list[int]]]:
    if not isinstance(value, Mapping) or len(value) > _MAX_RECONSTRUCTION_NODES:
        raise ValueError("reconstruction loss map is invalid or oversized")
    result: dict[str, list[list[int]]] = {}
    for raw_name, raw_blocks in value.items():
        name = _reconstruction_name(raw_name, kind="loss node")
        if name not in node_names:
            raise ValueError("reconstruction loss data name an unknown node")
        result[name] = _canonical_loss_blocks(
            raw_blocks,
            local_to_identity=local_to_identity,
            canonical=canonical,
            reference_count=reference_count,
        )
    return result


def _branch_reconstruction_events(
    output_root: Path,
    source: Mapping[str, object],
    *,
    canonical: Mapping[str, int],
    budget: list[int],
    warn: Callable[[str], None],
    checkpoint: Callable[[], None] | None,
) -> tuple[dict[str, list[int]], dict[str, list[list[int]]]]:
    directory = source["directory"]
    prefix = source["prefix"]
    identities = source["identities"]
    raw_nodes = source["nodes"]
    if (
        not isinstance(directory, Path)
        or not isinstance(prefix, str)
        or not isinstance(identities, Mapping)
        or not isinstance(raw_nodes, list)
    ):
        raise ValueError("reconstruction source is incomplete")
    path = directory / f"{prefix}_rec_gains_losses.json"
    if _safe_file(output_root, path) is None:
        return {}, {}
    node_names = frozenset(
        node["name"]
        for node in raw_nodes
        if isinstance(node, Mapping) and isinstance(node.get("name"), str)
    )
    local_to_identity = {number: identity for identity, number in identities.items()}
    reference_count = [0]
    try:
        document = _reconstruction_json(
            output_root,
            path,
            budget=budget,
            checkpoint=checkpoint,
        )
        if not isinstance(document, Mapping):
            raise ValueError("reconstruction branch events are not an object")
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        warn("branch gains/losses were malformed and were omitted")
        return {}, {}
    try:
        gains = _canonical_event_map(
            document.get("rec_gains"),
            node_names=node_names,
            local_to_identity=local_to_identity,
            canonical=canonical,
            reference_count=reference_count,
        )
    except ValueError:
        warn("branch gains were malformed and were omitted")
        gains = {}
    try:
        losses = _canonical_loss_map(
            document.get("rec_losses"),
            node_names=node_names,
            local_to_identity=local_to_identity,
            canonical=canonical,
            reference_count=reference_count,
        )
    except ValueError:
        warn("branch losses were malformed and were omitted")
        losses = {}
    return gains, losses


def _special_reconstruction_events(
    output_root: Path,
    source: Mapping[str, object],
    *,
    canonical: Mapping[str, int],
    budget: list[int],
    warn: Callable[[str], None],
    checkpoint: Callable[[], None] | None,
) -> dict[str, dict[str, list[int]]]:
    result: dict[str, dict[str, list[int]]] = {
        public_name: {} for public_name in _RECONSTRUCTION_EVENT_FIELDS.values()
    }
    directory = source["directory"]
    prefix = source["prefix"]
    identities = source["identities"]
    raw_nodes = source["nodes"]
    if (
        not isinstance(directory, Path)
        or not isinstance(prefix, str)
        or not isinstance(identities, Mapping)
        or not isinstance(raw_nodes, list)
    ):
        raise ValueError("reconstruction source is incomplete")
    path = directory / f"{prefix}_other_events.json"
    if _safe_file(output_root, path) is None:
        return result
    try:
        document = _reconstruction_json(
            output_root,
            path,
            budget=budget,
            checkpoint=checkpoint,
        )
        if not isinstance(document, Mapping):
            raise ValueError("reconstruction special events are not an object")
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        warn("special events were malformed and were omitted")
        return result
    node_names = frozenset(
        node["name"]
        for node in raw_nodes
        if isinstance(node, Mapping) and isinstance(node.get("name"), str)
    )
    local_to_identity = {number: identity for identity, number in identities.items()}
    reference_count = [0]
    for raw_name, public_name in _RECONSTRUCTION_EVENT_FIELDS.items():
        raw_events = document.get(raw_name, {})
        try:
            result[public_name] = _canonical_event_map(
                raw_events,
                node_names=node_names,
                local_to_identity=local_to_identity,
                canonical=canonical,
                reference_count=reference_count,
            )
        except ValueError:
            warn(f"{public_name} events were malformed and were omitted")
    return result


