"""Integration adapter manifest validation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path

from .files import _read_json, _safe_nonempty_output
from .validation_adapter import _adapter_group_leaves
from .validation_common import (
    _MAX_RESULT_GROUPS,
    _SUPPORTED_ADAPTER_SCHEMA_MAJOR,
    _supported_schema_major,
)


def _validated_adapter_groups(
    output_root: Path,
    adapter_stage: Mapping[str, object],
    checkpoint: Callable[[], None] | None = None,
) -> dict[str, frozenset[str]] | None:
    manifest_path = output_root / "adapter" / "manifest.json"
    if not _safe_nonempty_output(output_root, manifest_path):
        return None
    manifest = _read_json(output_root, manifest_path, checkpoint=checkpoint)
    if not _supported_schema_major(
        manifest,
        "crispr-tool-integration-adapter",
        _SUPPORTED_ADAPTER_SCHEMA_MAJOR,
    ):
        return None
    group_count = manifest.get("emitted_group_count")
    if isinstance(group_count, bool) or not isinstance(group_count, int):
        return None
    if (
        group_count < 0
        or group_count > _MAX_RESULT_GROUPS
        or adapter_stage.get("emitted_group_count") != group_count
    ):
        return None
    expected_group_leaves = _adapter_group_leaves(manifest, group_count)
    emitted_array_count = manifest.get("emitted_array_count")
    if (
        isinstance(emitted_array_count, bool)
        or not isinstance(emitted_array_count, int)
        or emitted_array_count < 0
        or expected_group_leaves is None
        or emitted_array_count
        != sum(len(leaves) for leaves in expected_group_leaves.values())
        or adapter_stage.get("emitted_array_count") != emitted_array_count
    ):
        return None
    return expected_group_leaves
