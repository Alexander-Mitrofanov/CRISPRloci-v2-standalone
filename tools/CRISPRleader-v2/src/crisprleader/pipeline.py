"""Validated context extraction; scientific prediction has a separate readiness gate."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from . import __version__
from .errors import ExecutionError, InputError
from .models import AnalysisConfig, json_record
from .validation import validate_result


def extract(fasta: Path, annotations: Path, *,
            config: AnalysisConfig | None = None) -> dict[str, Any]:
    """Extract source-bound context without asserting leader presence or boundaries."""
    from .coordinates import extract_candidates
    from .sequence_io import load_inputs

    config = config or AnalysisConfig()
    if type(config.flank_length) is not int or not 1 <= config.flank_length <= 5000:
        raise InputError("flank_length must be between 1 and 5000")
    sources, arrays = load_inputs(fasta, annotations, config=config)
    candidates: list[dict[str, Any]] = []
    for array in arrays:
        flanks = extract_candidates(sources[array.source_id], array, config=config)
        if len(candidates) + len(flanks) > config.max_candidates:
            raise InputError("Candidate count exceeds the configured bound")
        candidates.extend(json_record(flanks))
    report: dict[str, Any] = {
        "schema_version": "1.0.0",
        "tool": {"name": "CRisprleader v2", "version": __version__},
        "execution_status": "completed",
        "analysis_mode": "context_extraction",
        "coordinate_system": "zero_based_half_open",
        "sources": [{"id": source.id, "length": source.length, "sha256": source.sha256,
                     "topology": source.topology, "topology_provenance": source.topology_provenance}
                    for source in sources.values()],
        "arrays": json_record(arrays),
        "candidates": candidates,
        "predictions": [],
        "reference": None,
        "effective_config": json_record(config),
        "provenance": {"created_at": datetime.now(timezone.utc).isoformat(),
                       "input_filename": fasta.name, "annotation_filename": annotations.name},
        "warnings": ["Context extraction only: no scientific inference was performed.",
                     "Extracted windows are not inferred leader lengths or functional annotations."],
    }
    validate_result(report, sources)
    return report


def predict(fasta: Path, annotations: Path, *, references: Path | None = None,
            config: AnalysisConfig | None = None) -> dict[str, Any]:
    """Fail explicitly until the biological implementation and evaluation are complete."""
    raise ExecutionError(
        "Scientific prediction is unavailable in this development snapshot: the inference "
        "policy, default biological reference bundle and independent evaluation are incomplete. "
        "Use 'extract' for context windows only; see docs/STATUS.md."
    )
