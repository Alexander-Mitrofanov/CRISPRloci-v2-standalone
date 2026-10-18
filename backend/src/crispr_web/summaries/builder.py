"""Top-level public-summary composition."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Iterable

from .annotations import annotation_summary
from .leader import leader_summary
from .repeats import repeat_summary
from .protospacer import protospacer_summary
from .viral import viral_summary
from .repeat_type import repeat_type_summary
from .association import association_summary
from .array_compare import array_compare_summary
from .locus_compare import locus_compare_summary
from .prophage_context import prophage_context_summary
from .insertion_context import insertion_context_summary
from .defense_context import defense_context_summary
from .rna_family import rna_family_summary
from .repeat_map import repeat_map_summary
from .stages import (
    _adapter_summary,
    _detection_summary,
    _orientation_summary,
    _pipeline_summary,
    _reconstruction_summary,
)

def build_summary(
    output_root: Path,
    *,
    max_arrays: int,
    max_rows: int,
    max_files: int,
    accepted_categories: Iterable[str],
    checkpoint: Callable[[], None] | None = None,
) -> dict[str, Any]:
    warnings: list[str] = []
    if checkpoint is not None:
        checkpoint()
    pipeline_status, pipeline, bias_warnings = _pipeline_summary(
        output_root, max_rows, warnings, checkpoint
    )
    summary: dict[str, Any] = {
        "schema_version": "1.0.0",
        "pipeline_status": pipeline_status,
        "pipeline": pipeline,
        "detection": _detection_summary(
            output_root,
            max_arrays,
            max_files,
            frozenset(accepted_categories),
            warnings,
            checkpoint,
        ),
        "adapter": _adapter_summary(output_root, warnings, max_rows, checkpoint),
        "orientation": _orientation_summary(
            output_root, max_rows, max_files, warnings, checkpoint
        ),
        "reconstruction": _reconstruction_summary(
            output_root, max_rows, warnings, checkpoint
        ),
        "warnings": [*warnings, *bias_warnings],
    }
    cas, tracr, detection = annotation_summary(
        output_root, max_rows=max_rows, accepted_categories=frozenset(accepted_categories),
        warnings=warnings, checkpoint=checkpoint, max_arrays=max_arrays,
    )
    summary.update(casandra=cas, tracrrna=tracr)
    summary["crisprleader"] = leader_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["repeats"] = repeat_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["array_compare"] = array_compare_summary(output_root,max_rows=max_rows,checkpoint=checkpoint)
    summary["locus_compare"] = locus_compare_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["prophage_context"] = prophage_context_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["insertion_context"] = insertion_context_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["defense_context"] = defense_context_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["rna_family"] = rna_family_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["spacer_association"] = association_summary(output_root,max_rows=max_rows,checkpoint=checkpoint)
    summary["repeat_type"] = repeat_type_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["repeat_map"] = repeat_map_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["viral_search"] = viral_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    summary["protospacer"] = protospacer_summary(output_root, max_rows=max_rows, checkpoint=checkpoint)
    if detection is not None:
        summary["detection"] = detection
    summary["warnings"] = [*warnings, *bias_warnings]
    return summary
