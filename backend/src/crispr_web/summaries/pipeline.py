"""Sanitized Integration pipeline-stage summary extraction."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable, Mapping

from .common import (
    _BIAS_STAGE_NAMES,
    _MAX_BIAS_WARNINGS,
    _bias_fallback_metadata,
    _bias_policy,
    _json,
    _public_bias_warning,
    _safe_file,
)
from .reconstruction import (
    _mapping,
)

def _pipeline_summary(
    output_root: Path,
    max_rows: int,
    warnings: list[str],
    checkpoint: Callable[[], None] | None,
) -> tuple[str | None, dict[str, Any], list[dict[str, str]]]:
    path = output_root / "integration_result.json"
    if _safe_file(output_root, path) is None:
        warnings.append("Integration did not produce integration_result.json")
        return None, {}, []
    try:
        document = _mapping(_json(output_root, path, checkpoint=checkpoint))
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
        warnings.append(
            f"Could not summarize Integration result: {type(error).__name__}"
        )
        return None, {}, []
    stages: dict[str, Any] = {}
    public_bias_warnings: list[dict[str, str]] = []
    bias_warning_limit = min(max_rows * 2, _MAX_BIAS_WARNINGS)
    for name, stage in _mapping(document.get("stages")).items():
        if not isinstance(stage, Mapping):
            continue
        stage_name = str(name)
        public_stage = {
            key: stage.get(key)
            for key in (
                "status",
                "duration_seconds",
                "emitted_array_count",
                "emitted_group_count",
                "tree_source",
            )
            if key in stage
        }
        if stage_name in _BIAS_STAGE_NAMES:
            policy = _bias_policy(stage_name, stage)
            if policy is not None:
                public_stage["bias_correction_policy"] = policy
            elif "bias_correction_policy" in stage:
                warnings.append(
                    f"Malformed {stage_name} bias-correction policy was omitted."
                )

            has_fallback_metadata = (
                "bias_fallback_count" in stage or "bias_warnings" in stage
            )
            if has_fallback_metadata:
                fallback_metadata = (
                    _bias_fallback_metadata(
                        stage_name,
                        stage,
                        policy,
                        max_warnings=bias_warning_limit,
                    )
                    if policy is not None
                    else None
                )
                if fallback_metadata is None:
                    warnings.append(
                        f"Malformed {stage_name} bias-fallback provenance was omitted."
                    )
                else:
                    fallback_count, stage_bias_warnings = fallback_metadata
                    public_stage["bias_fallback_count"] = fallback_count
                    public_stage["bias_warnings"] = stage_bias_warnings
                    public_bias_warnings.extend(
                        _public_bias_warning(stage_name, warning)
                        for warning in stage_bias_warnings
                    )
        stages[stage_name] = public_stage
    return (
        str(document.get("status")) if document.get("status") is not None else None,
        {"stages": stages},
        public_bias_warnings,
    )


