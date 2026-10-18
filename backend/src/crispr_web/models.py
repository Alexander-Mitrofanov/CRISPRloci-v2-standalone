"""Typed public API models."""

from __future__ import annotations

from enum import StrEnum
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from .scientific_options import validate_tool_options


class ScientificMode(StrEnum):
    detection = "detection"
    orientation = "orientation"
    reconstruction = "reconstruction"
    cas = "cas"
    tracrrna = "tracrrna"
    loci = "loci"
    leader = "leader"
    repeats = "repeats"
    repeat_context = "repeat_context"
    protospacer = "protospacer"
    viral_search = "viral_search"
    repeat_map = "repeat_map"
    repeat_type = "repeat_type"
    spacer_association = "spacer_association"
    array_compare = "array_compare"
    locus_compare = "locus_compare"
    prophage_context = "prophage_context"
    insertion_context = "insertion_context"
    defense_context = "defense_context"
    rna_family = "rna_family"


class CategoryPolicy(StrEnum):
    bona_fide_possible = "bona_fide_possible"
    bona_fide_only = "bona_fide_only"


class JobStatus(StrEnum):
    queued = "queued"
    running = "running"
    completed = "completed"
    completed_no_eligible_groups = "completed_no_eligible_groups"
    failed = "failed"
    cancelled = "cancelled"


class JobSubmission(BaseModel):
    model_config = ConfigDict(extra="forbid")

    molecule: Literal["DNA", "RNA"] = Field(default="DNA", description="RNA is supported for supplied repeats and homologous RNA families, in transcribed 5-prime to 3-prime orientation.")

    sequence: str = Field(
        description="IUPAC FASTA, raw sequence, or a serialized crisprloci.sequence-set v1 JSON file. JSON kind must match the mode's input port. RNA requires repeats or rna_family mode and molecule=RNA."
    )
    filename: str | None = Field(
        default=None,
        max_length=255,
        description="Display-only source filename; paths are discarded.",
    )
    mode: ScientificMode = Field(
        description=(
            "detection runs CRISPRidentify; reconstruction additionally runs "
            "SpacerPlacer; orientation runs evOr and its selected reconstruction; "
            "cas runs CasAndra; tracrrna screens covariance models; loci combines "
            "CRISPRidentify, CasAndra, complete tracrRNA and flank context evidence; leader extracts both adjacent windows without predicting a leader; repeats folds supplied repeat sequences; repeat_context detects arrays and compares isolated/local repeat structure; protospacer searches submitted DNA against an immutable observed-spacer collection for exact full-length matches; repeat_type predicts repeat-derived subtype evidence using a frozen 37-class model; repeat_map maps supplied repeats to a frozen annotated reference with at most three global edits; viral_search searches supplied spacers against an immutable viral RefSeq BLAST reference; locus_compare aligns source-verified protein neighborhoods; prophage_context joins source-spacer matches to candidate prophage intervals; insertion_context maps insertion-sequence evidence beside Cas annotations; defense_context identifies systems within pinned model coverage; rna_family aligns user-justified homologous RNA families and predicts a consensus structure."
        )
    )
    association_grouping: Literal["single_set", "per_record", "per_source"] | None = Field(default=None, description="Required for spacer association: explicit grouping assertion; source records are not automatically hosts.")
    viral_max_mismatches: int = Field(default=2, ge=0, le=2, strict=True, description="Maximum substitutions in full-query viral BLAST matches; retrieval remains heuristic.")
    leader_flank_length: int = Field(default=500, ge=1, le=5000, strict=True, description="Bases extracted on each side of accepted arrays; not a predicted leader boundary.")
    tracr_model_type: Literal["II", "V"] = Field(default="II", description="Type II or Type V-K covariance models for tracrRNA workflows.")
    tool_options: dict[str, dict[str, Any]] = Field(default_factory=dict, description="Overrides keyed by scientific tool and its exact shipped CLI flag. Only allowlisted flags applicable to this workflow are accepted.")

    @model_validator(mode="after")
    def validated_tool_options(self):
        if self.mode.value == "spacer_association" and self.association_grouping is None:
            raise ValueError("Choose an explicit spacer association grouping")
        self.tool_options = validate_tool_options(self.tool_options, self.mode.value)
        return self

    category_policy: CategoryPolicy = CategoryPolicy.bona_fide_possible
    spacer_distance: int = Field(default=1, ge=0, le=2, strict=True)
    bias_corrections: bool = Field(
        default=True,
        strict=True,
        description="Request bias corrections; operator safety policy may disable them.",
    )


class InputSummary(BaseModel):
    filename: str
    record_count: int
    base_count: int
    sha256: str
    source_ids: list[str]


class JobOptions(BaseModel):
    tool_options: dict[str, dict[str, Any]] = Field(default_factory=dict)
    molecule: Literal["DNA", "RNA"] = "DNA"
    association_grouping: Literal["single_set", "per_record", "per_source"] | None = None
    viral_max_mismatches: int = 2
    leader_flank_length: int = 500
    tracr_model_type: Literal["II", "V"] = "II"
    category_policy: CategoryPolicy
    categories: list[str]
    spacer_distance: int
    bias_corrections_requested: bool
    bias_corrections_effective: bool
    warnings: list[str] = Field(default_factory=list)


class ArtifactView(BaseModel):
    artifact_id: str
    name: str
    size_bytes: int
    sha256: str
    media_type: str
    download_url: str


class ErrorView(BaseModel):
    code: str
    message: str


class JobView(BaseModel):
    job_id: str
    status: JobStatus
    mode: ScientificMode
    created_at: str
    updated_at: str
    started_at: str | None = None
    finished_at: str | None = None
    expires_at: str | None = None
    queue_position: int | None = None
    cancel_requested: bool = False
    input: InputSummary
    options: JobOptions
    summary: dict[str, Any] | None = None
    artifacts: list[ArtifactView] = Field(default_factory=list)
    error: ErrorView | None = None


class JobCreated(BaseModel):
    access_token: str = Field(
        description="Shown once. Supply as `Authorization: Bearer …` for this job."
    )
    token_type: Literal["Bearer"] = "Bearer"
    job: JobView


class CancelResponse(BaseModel):
    job: JobView


class PublicConfig(BaseModel):
    service: str
    api_version: str
    modes: list[ScientificMode]
    category_policies: list[CategoryPolicy]
    max_request_bytes: int
    max_archive_bytes: int
    max_total_bases: int
    max_record_bases: int
    max_records: int
    max_header_characters: int
    max_queued_jobs: int
    max_active_jobs: int
    max_active_jobs_per_client: int
    identify_folder_jobs: int
    retention_seconds: int
    accepted_categories: list[str]
    model_score_is_probability: Literal[False] = False


class ScientificToolIdentity(BaseModel):
    name: str
    display_name: str
    version: str
    kind: Literal["source_tree", "wheel", "git_source"]
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    sha256_scope: Literal["source_tree_v1", "file_bytes", "source_inventory_v1"]
    repository: str | None = Field(default=None, exclude_if=lambda value: value is None)
    commit: str | None = Field(default=None, exclude_if=lambda value: value is None)


class HealthResponse(BaseModel):
    status: Literal["ok", "degraded"]
    database: Literal["ok", "error"]
    worker: Literal["ok", "stale", "not_started", "storage_low", "error"]
    storage: Literal["ok", "low", "error"]
    scientific_release: Literal["ok", "error"]
    service: str
    version: str
    release_id: str
    backend_wheel_sha256: str
    scientific_release_manifest_sha256: str
    scientific_artifacts: dict[str, ScientificToolIdentity]


class VersionResponse(BaseModel):
    service: str
    version: str
    api_version: str
    integration_contract: str
    release_id: str
    backend_wheel_sha256: str
    scientific_release_manifest_sha256: str
    scientific_artifacts: dict[str, ScientificToolIdentity]
