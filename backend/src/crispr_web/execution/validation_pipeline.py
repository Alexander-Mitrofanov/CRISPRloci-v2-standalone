"""End-to-end scientific pipeline output validation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from pathlib import Path

from ..db import ClaimedJob
from .contracts import SUCCESS_STATUSES
from .files import _read_json
from .validation_common import _nonempty_string, _schema_version_is
from .validation_identify import _validate_identify_reports
from .validation_manifest import _validated_adapter_groups
from .validation_orientation import _validate_orientation_outputs
from .validation_reconstruction import _validate_reconstruction_outputs


def _validated_pipeline_status(
    claimed: ClaimedJob,
    output_root: Path,
    expected_source_bindings: Mapping[str, tuple[int, str]],
    checkpoint: Callable[[], None] | None = None,
    *, spacer_reference: Path | None = None, association_reference: Path | None = None, input_root: Path | None = None, repeat_reference: Path | None = None, viral_reference: Path | None = None,
) -> str | None:
    """Validate cross-file scientific semantics instead of trusting exit status."""

    integration = _read_json(
        output_root,
        output_root / "integration_result.json",
        checkpoint=checkpoint,
    )
    if not _schema_version_is(integration, "crispr-tool-integration-run", "1.0.0"):
        return None
    # Older runtimes must not silently ignore per-job overrides. Default-only
    # jobs retain compatibility with reports produced before this optional field.
    if integration.get("tool_options", {}) != claimed.tool_options:
        return None
    status_value = integration.get("status")
    if status_value not in SUCCESS_STATUSES:
        return None
    stages = integration.get("stages")
    if not isinstance(stages, Mapping):
        return None
    if claimed.mode == "locus_compare":
        from .validation_locus_compare import validate_locus_compare_pipeline
        return "completed" if validate_locus_compare_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    if claimed.mode == "prophage_context":
        from .validation_prophage_context import validate_prophage_context_pipeline
        return "completed" if validate_prophage_context_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    if claimed.mode == "insertion_context":
        from .validation_insertion_context import validate_insertion_context_pipeline
        return "completed" if validate_insertion_context_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    if claimed.mode == "defense_context":
        from .validation_defense_context import validate_defense_context_pipeline
        return "completed" if validate_defense_context_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    if claimed.mode == "rna_family":
        from .validation_rna_family import validate_rna_family_pipeline
        return "completed" if validate_rna_family_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    if claimed.mode == "array_compare":
        from .validation_array_compare import validate_array_compare_pipeline
        return "completed" if validate_array_compare_pipeline(claimed,output_root,expected_source_bindings,integration,input_root,checkpoint) else None
    if claimed.mode == "spacer_association":
        from .validation_association import validate_association_pipeline
        return "completed" if validate_association_pipeline(claimed,output_root,expected_source_bindings,integration,input_root,association_reference,checkpoint) else None
    if claimed.mode == "repeat_type":
        from .validation_repeat_type import validate_repeat_type_pipeline
        return "completed" if validate_repeat_type_pipeline(claimed, output_root, expected_source_bindings, integration, input_root, checkpoint) else None
    if claimed.mode == "repeat_map":
        from .validation_repeat_map import validate_repeat_map_pipeline
        return "completed" if validate_repeat_map_pipeline(claimed, output_root, expected_source_bindings,
            integration, input_root, repeat_reference, checkpoint) else None
    if claimed.mode == "viral_search":
        from .validation_viral import validate_viral_pipeline
        return "completed" if validate_viral_pipeline(claimed, output_root, expected_source_bindings,
            integration, input_root, viral_reference, checkpoint) else None
    if claimed.mode == "protospacer":
        from .validation_protospacer import validate_protospacer_pipeline
        return "completed" if validate_protospacer_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint, database_path=spacer_reference
        ) else None
    if claimed.mode in {"repeats", "repeat_context"}:
        from .validation_repeats import validate_repeat_pipeline
        return "completed" if validate_repeat_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    if claimed.mode in {"leader", "loci"}:
        from .validation_leader import validate_leader_pipeline
        if not validate_leader_pipeline(claimed, output_root, expected_source_bindings, integration, input_root, checkpoint):
            return None
        if claimed.mode == "leader":
            return "completed"
    if claimed.mode in {"cas", "tracrrna", "loci"}:
        from .validation_annotations import validate_annotation_pipeline
        return "completed" if validate_annotation_pipeline(
            claimed, output_root, expected_source_bindings, integration, input_root, checkpoint
        ) else None
    identify = stages.get("crispridentify")
    adapter = stages.get("adapter")
    if not isinstance(identify, Mapping) or identify.get("status") != "completed":
        return None
    if not isinstance(adapter, Mapping) or adapter.get("status") != "completed":
        return None
    if not _validate_identify_reports(
        output_root,
        expected_source_bindings,
        checkpoint,
    ):
        return None

    expected_group_leaves = _validated_adapter_groups(
        output_root,
        adapter,
        checkpoint,
    )
    if expected_group_leaves is None:
        return None
    group_count = len(expected_group_leaves)
    expected_groups = set(expected_group_leaves)

    if claimed.mode == "detection":
        return "completed" if status_value == "completed" else None

    stage_name = "spacerplacer" if claimed.mode == "reconstruction" else "crispr_evor"
    downstream = stages.get(stage_name)
    if not isinstance(downstream, Mapping):
        return None
    if group_count == 0:
        if status_value != "completed_no_eligible_groups":
            return None
        return (
            "completed_no_eligible_groups"
            if downstream.get("status") == "skipped"
            else None
        )
    if status_value != "completed" or downstream.get("status") != "completed":
        return None
    result_groups = downstream.get("result_groups")
    if (
        downstream.get("result_count") != group_count
        or not isinstance(result_groups, list)
        or len(result_groups) != group_count
        or any(not _nonempty_string(group) for group in result_groups)
        or set(result_groups) != expected_groups
    ):
        return None

    valid_outputs = (
        _validate_reconstruction_outputs(
            output_root,
            expected_group_leaves,
            checkpoint,
            expected_bias_corrections=claimed.bias_corrections_effective,
            downstream_stage=downstream,
        )
        if claimed.mode == "reconstruction"
        else _validate_orientation_outputs(
            output_root,
            downstream,
            expected_group_leaves,
            checkpoint,
            expected_bias_corrections=claimed.bias_corrections_effective,
            expected_orientation_boundary=claimed.tool_options.get("crispr_evor", {}).get("--orientation_decision_boundary"),
        )
    )
    return "completed" if valid_outputs else None
