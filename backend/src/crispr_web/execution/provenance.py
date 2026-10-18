"""Sanitized worker and scientific release provenance."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .. import __version__
from ..db import ClaimedJob, utc_iso
from ..service import CATEGORY_POLICY_CATEGORIES
from .contracts import MODE_TO_TOOL
from .context_provenance import CONTEXT_CITATIONS
from .files import _read_json


def _provenance(
    *,
    job: Mapping[str, Any],
    claimed: ClaimedJob,
    output_root: Path,
    summary: Mapping[str, Any],
    command: list[str],
    backend_release_id: str,
    backend_wheel_sha256: str,
    scientific_release_manifest_sha256: str,
    scientific_artifacts: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    integration = _read_json(output_root, output_root / "integration_result.json")
    stages: dict[str, Any] = {}
    raw_stages = integration.get("stages")
    if isinstance(raw_stages, Mapping):
        for name, raw_stage in raw_stages.items():
            if not isinstance(raw_stage, Mapping):
                continue
            stages[str(name)] = {
                key: raw_stage[key]
                for key in (
                    "status",
                    "duration_seconds",
                    "return_code",
                    "tree_source",
                    "tree_policy",
                )
                if key in raw_stage
            }
    orientation = _read_json(
        output_root,
        output_root / "crispr_evor" / "crispr_evor_results.json",
    )
    warnings = list(claimed.option_warnings)
    summary_warnings = summary.get("warnings")
    if isinstance(summary_warnings, list):
        warnings.extend(str(item) for item in summary_warnings)
    return {
        "job": {
            "job_id": claimed.job_id,
            "mode": claimed.mode,
            "attempt": claimed.attempt,
            "created_at": utc_iso(job.get("created_at")),
        },
        "input": {
            "normalized_sha256": job.get("input_sha256"),
            "source_ids": job.get("source_ids", []),
            "record_count": job.get("record_count"),
            "base_count": job.get("base_count"),
        },
        "options": {
            "tool_options": claimed.tool_options,
            "molecule": claimed.molecule,
            "spacer_database_sha256": claimed.spacer_database_sha256,
            "association_reference_sha256": claimed.association_reference_sha256,
            "association_groups_sha256": claimed.association_groups_sha256,
            "array_map_sha256": claimed.array_map_sha256,
            "association_grouping": claimed.association_grouping,
            "repeat_reference_sha256": claimed.repeat_reference_sha256,
            "viral_reference_sha256": claimed.viral_reference_sha256,
            "viral_max_mismatches": claimed.viral_max_mismatches,
            "leader_flank_length": claimed.leader_flank_length, "tracr_model_type": claimed.tracr_model_type,
            "category_policy": claimed.category_policy,
            "categories": list(CATEGORY_POLICY_CATEGORIES[claimed.category_policy]),
            "spacer_distance": claimed.spacer_distance,
            "bias_corrections_requested": claimed.bias_corrections_requested,
            "bias_corrections_effective": claimed.bias_corrections_effective,
        },
        "software": {
            "backend_version": __version__,
            "backend_release_id": backend_release_id,
            "backend_wheel_sha256": backend_wheel_sha256,
            "scientific_release_manifest_sha256": (scientific_release_manifest_sha256),
            "scientific_artifacts": {
                key: dict(identity) for key, identity in scientific_artifacts.items()
            },
            "integration_schema": integration.get("schema"),
            "crispridentify_report_schema_versions": (
                summary.get("detection", {}).get("report_schema_versions", {})
                if isinstance(summary.get("detection"), Mapping)
                else {}
            ),
            "crispr_evor_version": orientation.get("version"),
        },
        "execution": {
            "logical_invocation": {
                "tool_options": claimed.tool_options,
                "action": "run",
                "analysis": MODE_TO_TOOL[claimed.mode],
                "categories": list(CATEGORY_POLICY_CATEGORIES[claimed.category_policy]),
                "spacer_distance": claimed.spacer_distance,
                "bias_corrections": claimed.bias_corrections_effective,
                "cpu": int(command[command.index("--cpu") + 1]),
                "identify_folder_jobs": min(
                    int(command[command.index("--identify-folder-jobs") + 1]),
                    int(job.get("record_count") or 1),
                    int(command[command.index("--cpu") + 1]),
                ),
                "identify_folder_jobs_configured": int(
                    command[command.index("--identify-folder-jobs") + 1]
                ),
                "stage_timeout_seconds": int(
                    command[command.index("--stage-timeout") + 1]
                ),
            },
            "pipeline_status": integration.get("status"),
            "stages": stages,
        },
        "citations": CONTEXT_CITATIONS[claimed.mode] if claimed.mode in CONTEXT_CITATIONS else ([{"tool":"CCTK", "role":"Observed spacer sharing and ordered occurrence comparison", "doi":"10.1089/crispr.2022.0080", "url":"https://doi.org/10.1089/crispr.2022.0080"}] if claimed.mode == "array_compare" else [{"tool":"SpacePHARER", "role":"Offline grouped spacer-to-target association hypothesis with empirical controls", "doi":"10.1093/bioinformatics/btab222", "url":"https://doi.org/10.1093/bioinformatics/btab222"}] if claimed.mode == "spacer_association" else [{"tool": "RepeatTyper / CRISPRCasTyper", "role": "Repeat-derived subtype model evidence; closed-set scores are not calibrated confidence", "doi": "10.1089/crispr.2020.0059", "url": "https://doi.org/10.1089/crispr.2020.0059"}] if claimed.mode == "repeat_type" else [{"tool": "CRISPRmap", "role": "Historical repeat reference mapping; near hits are candidates", "doi": "10.1093/nar/gkt606", "url": "https://doi.org/10.1093/nar/gkt606"}] if claimed.mode == "repeat_map" else [{"tool": "NCBI BLAST+ and viral RefSeq", "role": "Heuristic short-spacer alignment against viral references", "url": "https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/"}] if claimed.mode == "viral_search" else [{"tool": "CRISPRspacer reference collection", "role": "Observed reference spacer sequences",
                        "url": summary.get("protospacer", {}).get("reference", {}).get("source_url")} ]
                      if claimed.mode == "protospacer" else ([
            {"tool": "ViennaRNA", "role": "RNA folding and local pair/accessibility calculations",
             "doi": "10.1186/1748-7188-6-26", "url": "https://doi.org/10.1186/1748-7188-6-26"},
            {"tool": "CRISPRrepeat method background", "role": "Experimental repeat-context study, not a validation of predicted cleavage",
             "doi": "10.1093/nar/gkw786", "url": "https://doi.org/10.1093/nar/gkw786"},
        ] if claimed.mode in {"repeats", "repeat_context"} else []) + [
            {
                "tool": "CRISPRidentify",
                "role": "CRISPR array detection and classification",
                "doi": "10.1093/nar/gkaa1158",
                "url": "https://doi.org/10.1093/nar/gkaa1158",
            },
            {
                "tool": "SpacerPlacer",
                "role": "ancestral CRISPR spacer-array reconstruction",
                "doi": "10.1093/nar/gkae772",
                "url": "https://doi.org/10.1093/nar/gkae772",
            },
            {
                "tool": "CRISPR-evOr",
                "role": "evolutionary acquisition-orientation comparison",
                "doi": "10.1371/journal.pcbi.1013706",
                "url": "https://doi.org/10.1371/journal.pcbi.1013706",
            },
        ]),
        "viral_reference": summary.get("viral_search", {}).get("reference") if claimed.mode == "viral_search" else None,
        "protospacer_reference": summary.get("protospacer", {}).get("reference") if claimed.mode == "protospacer" else None,
        "repeat_analysis": summary.get("repeats", {}).get("model") if claimed.mode in {"repeats", "repeat_context"} else None,
        "warnings": warnings,
    }
