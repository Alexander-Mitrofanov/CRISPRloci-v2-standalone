"""FastAPI application for anonymous, capability-protected scientific jobs."""

from __future__ import annotations

from contextlib import asynccontextmanager
import logging
import os
import threading
import time
from typing import AsyncIterator

from fastapi import Depends, FastAPI, HTTPException, Request, Response, status
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer

from .capacity import available_capacity
from .config import Settings, origin_strings
from .db import Store
from .http.middleware import (
    JsonErrorMiddleware,
    RequestLimitMiddleware,
    SecurityHeadersMiddleware,
)
from .http.routes.jobs import build_job_router
from .http.routes.catalog import build_catalog_router
from .http.routes.service import API_VERSION, build_service_router
from .models import HealthResponse
from .service import JobService, require_unprivileged_process


LOGGER = logging.getLogger("crispr_web.api")
bearer = HTTPBearer(auto_error=False, scheme_name="JobAccessToken")

__all__ = [
    "API_VERSION",
    "JsonErrorMiddleware",
    "RequestLimitMiddleware",
    "SecurityHeadersMiddleware",
    "create_app",
]


def create_app(settings: Settings | None = None) -> FastAPI:
    require_unprivileged_process()
    os.umask(0o077)
    resolved = settings or Settings.from_env()
    store = Store(resolved)
    service = JobService(resolved, store)
    readiness_lock = threading.Lock()
    readiness_cached: tuple[float, dict[str, str]] | None = None

    def probe_readiness(*, force: bool) -> dict[str, str]:
        nonlocal readiness_cached
        now = time.monotonic()
        cached = readiness_cached
        if not force and cached is not None and now < cached[0]:
            return cached[1]
        with readiness_lock:
            now = time.monotonic()
            cached = readiness_cached
            if not force and cached is not None and now < cached[0]:
                return cached[1]
            try:
                database = "ok" if store.database_ok() else "error"
            except Exception:
                LOGGER.exception("database readiness probe failed")
                database = "error"
            try:
                worker = store.worker_health()
            except Exception:
                LOGGER.exception("worker readiness probe failed")
                worker = "error"
            try:
                free_bytes, free_inodes = available_capacity(resolved.data_root)
                required_bytes = max(
                    resolved.submission_admission_free_bytes,
                    resolved.worker_admission_free_bytes,
                )
                required_inodes = max(
                    resolved.submission_admission_free_inodes,
                    resolved.worker_admission_free_inodes,
                )
                storage_state = (
                    "ok"
                    if (free_bytes >= required_bytes and free_inodes >= required_inodes)
                    else "low"
                )
            except OSError:
                storage_state = "error"
            try:
                resolved.verify_scientific_release()
                scientific_release = "ok"
            except Exception:
                LOGGER.exception("scientific release readiness probe failed")
                scientific_release = "error"
            values = {
                "database": database,
                "worker": worker,
                "storage": storage_state,
                "scientific_release": scientific_release,
            }
            readiness_cached = (
                time.monotonic() + resolved.readiness_cache_seconds,
                values,
            )
            return values

    @asynccontextmanager
    async def lifespan(_app: FastAPI) -> AsyncIterator[None]:
        service.cleanup_expired()
        yield

    app = FastAPI(
        title="CRISPR evolutionary workflow API",
        summary="CRISPRidentify, SpacerPlacer, and CRISPR-evOr as durable jobs",
        description=(
            "Submit bounded IUPAC DNA FASTA. A per-job 256-bit bearer capability "
            "is required to inspect, cancel, or download that job. Raw CRISPRidentify "
            "model scores are not calibrated probabilities; category is the primary result."
        ),
        version=resolved.service_version,
        lifespan=lifespan,
        docs_url=None,
        openapi_url="/api/v1/openapi.json",
        redoc_url=None,
    )
    app.state.settings = resolved
    app.state.store = store
    app.state.service = service

    @app.exception_handler(RequestValidationError)
    async def bounded_validation_error(
        _request: Request, error: RequestValidationError
    ) -> JSONResponse:
        details = []
        for item in error.errors()[:50]:
            location = [
                part if isinstance(part, int) else str(part)[:120]
                for part in item.get("loc", ())[:20]
            ]
            details.append(
                {
                    "loc": location,
                    "msg": str(item.get("msg") or "invalid request")[:240],
                    "type": str(item.get("type") or "validation_error")[:120],
                }
            )
        return JSONResponse(
            status_code=status.HTTP_422_UNPROCESSABLE_CONTENT,
            content={"detail": details},
        )

    app.add_middleware(RequestLimitMiddleware, max_bytes=resolved.max_request_bytes)
    app.add_middleware(JsonErrorMiddleware)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origin_strings(resolved.cors_origins),
        allow_credentials=False,
        allow_methods=["GET", "POST", "DELETE", "OPTIONS"],
        allow_headers=["Authorization", "Content-Type"],
        expose_headers=["Content-Disposition", "Content-Length", "Retry-After"],
        max_age=600,
    )
    # Starlette prepends middleware registrations. Register this last so every
    # response, including CORS preflights and JSON error fallbacks, is protected.
    app.add_middleware(SecurityHeadersMiddleware, hsts=resolved.hsts)

    def capability(
        job_id: str,
        credentials: HTTPAuthorizationCredentials | None = Depends(bearer),
    ) -> str:
        if credentials is None or credentials.scheme.lower() != "bearer":
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="a job bearer token is required",
                headers={"WWW-Authenticate": "Bearer"},
            )
        if not service.authorize(job_id, credentials.credentials):
            # Do not reveal whether the locator or the capability was wrong.
            raise HTTPException(status_code=404, detail="job not found")
        return credentials.credentials

    def readiness_response(response: Response, *, force: bool) -> HealthResponse:
        values = probe_readiness(force=force)
        overall = (
            "ok"
            if values["database"] == "ok"
            and values["worker"] == "ok"
            and values["storage"] == "ok"
            and values["scientific_release"] == "ok"
            else "degraded"
        )
        if overall == "degraded":
            response.status_code = status.HTTP_503_SERVICE_UNAVAILABLE
        return HealthResponse(
            status=overall,
            database=values["database"],
            worker=values["worker"],
            storage=values["storage"],
            service="crisprloci-v2",
            version=resolved.service_version,
            release_id=resolved.release_id,
            backend_wheel_sha256=resolved.backend_wheel_sha256,
            scientific_release=values["scientific_release"],
            scientific_release_manifest_sha256=(
                resolved.scientific_release_manifest_sha256
            ),
            scientific_artifacts=resolved.scientific_release_identities,
        )

    app.include_router(build_service_router(resolved, readiness_response))
    app.include_router(build_job_router(resolved, service, capability))
    app.include_router(build_catalog_router(resolved))

    return app
