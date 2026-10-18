"""Public service metadata and readiness routes."""

from __future__ import annotations

from ...context_modes import context_mode_available
from typing import Protocol

from fastapi import APIRouter, HTTPException, Request, Response

from ... import __version__
from ...config import Settings
from ...models import (
    CategoryPolicy,
    HealthResponse,
    PublicConfig,
    ScientificMode,
    VersionResponse,
)
from ..middleware import _is_direct_loopback_request


API_VERSION = "v1"


class ReadinessResponder(Protocol):
    def __call__(self, response: Response, *, force: bool) -> HealthResponse: ...


def build_service_router(
    resolved: Settings,
    readiness_response: ReadinessResponder,
) -> APIRouter:
    router = APIRouter()

    @router.get("/api/v1/health", response_model=HealthResponse, tags=["service"])
    def health(response: Response) -> HealthResponse:
        return readiness_response(response, force=False)

    @router.get(
        "/api/v1/health/fresh",
        response_model=HealthResponse,
        include_in_schema=False,
    )
    def fresh_health(request: Request, response: Response) -> HealthResponse:
        if not _is_direct_loopback_request(request):
            raise HTTPException(status_code=404, detail="not found")
        return readiness_response(response, force=True)

    @router.get("/api/v1/config", response_model=PublicConfig, tags=["service"])
    def public_config() -> PublicConfig:
        return PublicConfig(
            service="crisprloci-v2",
            api_version=API_VERSION,
            modes=[mode for mode in ScientificMode if context_mode_available(mode.value, resolved.scientific_release_identities) and (mode != ScientificMode.array_compare or "cctk" in resolved.scientific_release_identities) and (mode != ScientificMode.spacer_association or resolved.association_reference_path is not None) and (mode != ScientificMode.repeat_type or "repeattyper" in resolved.scientific_release_identities) and (mode != ScientificMode.protospacer or resolved.spacer_database_path is not None) and (mode != ScientificMode.repeat_map or resolved.repeat_reference_path is not None) and (mode != ScientificMode.viral_search or resolved.viral_reference_path is not None)],
            category_policies=list(CategoryPolicy),
            max_request_bytes=resolved.max_request_bytes,
            max_archive_bytes=resolved.max_archive_bytes,
            max_total_bases=resolved.max_total_bases,
            max_record_bases=resolved.max_record_bases,
            max_records=resolved.max_records,
            max_header_characters=resolved.max_header_characters,
            max_queued_jobs=resolved.max_queued_jobs,
            max_active_jobs=resolved.max_active_jobs,
            max_active_jobs_per_client=resolved.max_active_jobs_per_client,
            identify_folder_jobs=resolved.identify_folder_jobs,
            retention_seconds=resolved.retention_seconds,
            accepted_categories=["Bona-fide", "Possible"],
        )

    @router.get("/api/v1/version", response_model=VersionResponse, tags=["service"])
    def version() -> VersionResponse:
        return VersionResponse(
            service="crisprloci-v2",
            version=resolved.service_version or __version__,
            api_version=API_VERSION,
            integration_contract="crispr-tool-integration-run/1.0.0",
            release_id=resolved.release_id,
            backend_wheel_sha256=resolved.backend_wheel_sha256,
            scientific_release_manifest_sha256=(
                resolved.scientific_release_manifest_sha256
            ),
            scientific_artifacts=resolved.scientific_release_identities,
        )

    return router
