"""Public, source-preserving catalogue records, separate from private jobs."""

from __future__ import annotations

from typing import Generic, TypeVar

from pydantic import BaseModel, ConfigDict, Field


class CatalogModel(BaseModel):
    model_config = ConfigDict(extra="forbid")


class CatalogCounts(CatalogModel):
    assemblies: int = Field(ge=0)
    arrays: int = Field(ge=0)
    repeats: int = Field(ge=0)
    spacers: int = Field(ge=0)
    cas_systems: int = Field(ge=0)
    cas_genes: int = Field(ge=0)


class CatalogCoverage(CatalogModel):
    catalog_assemblies: int | None = Field(default=None, ge=0)
    successful_assemblies: int | None = Field(default=None, ge=0)
    failed_assemblies: int | None = Field(default=None, ge=0)
    unavailable_assemblies: int | None = Field(default=None, ge=0)
    imported_assemblies: int | None = Field(default=None, ge=0)


class CatalogSummary(CatalogModel):
    available: bool
    snapshot_id: str | None = Field(default=None, max_length=200)
    title: str | None = Field(default=None, max_length=300)
    source_url: str | None = Field(default=None, max_length=1000)
    generated_at: str | None = Field(default=None, max_length=100)
    coordinate_system: str = "1-based inclusive; forward-source coordinates"
    counts: CatalogCounts | None = None
    coverage: CatalogCoverage | None = None
    limitations: list[str] = Field(default_factory=list, max_length=100)
    provenance: dict[str, str] = Field(default_factory=dict, max_length=100)


class CatalogAssembly(CatalogModel):
    id: int
    accession: str
    status: str
    validation_status: str | None
    sequence_count: int | None
    array_count: int | None
    cas_system_count: int | None
    cas_gene_count: int | None


class CatalogCall(CatalogModel):
    id: int
    accession: str
    sequence_id: str
    tool: str
    tool_version: str | None
    call_id: str
    start: int | None
    end: int | None
    length: int | None


class CatalogArray(CatalogCall):
    strand: str | None
    repeat_consensus: str | None
    repeat_count: int | None
    spacer_count: int | None
    evidence_level: str | None
    category: str | None
    certainty_score: float | None


class CatalogRepeat(CatalogModel):
    id: int
    array_id: int
    accession: str
    sequence_id: str
    ordinal: int
    start: int | None
    end: int | None
    length: int | None
    sequence: str | None


class CatalogSpacer(CatalogRepeat):
    kind: str
    boundary_position: int | None


class CatalogCasSystem(CatalogCall):
    cas_class: str | None
    cas_type: str | None
    subtype: str | None
    confidence: float | None
    confidence_kind: str | None
    experimental: bool | None
    missing_gene_count: int | None


class CatalogCasGene(CatalogCall):
    cas_system_id: int | None
    strand: str | None
    annotation: str | None
    profile: str | None
    cas_class: str | None
    cas_type: str | None
    subtype: str | None
    score: float | None
    score_kind: str | None
    protein_length: int | None
    translation_table: int | None
    is_cas: bool | None


CatalogRecord = TypeVar("CatalogRecord", bound=CatalogModel)


class CatalogPage(CatalogModel, Generic[CatalogRecord]):
    items: list[CatalogRecord]
    next_cursor: str | None
    limit: int
