"""Capability-protected job lifecycle and download routes."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, Request, status
from fastapi.responses import FileResponse

from ...config import Settings
from ...db import ClientCapacityError, QueueCapacityError, SubmissionRateLimitError
from ...fasta import FastaValidationError
from ...models import CancelResponse, JobCreated, JobSubmission, JobView
from ...security import client_ip
from ...service import JobService, ReferenceUnavailableError, StorageCapacityError


def build_job_router(
    resolved: Settings,
    service: JobService,
    capability: Callable[..., str],
) -> APIRouter:
    router = APIRouter()

    @router.post(
        "/api/v1/jobs",
        response_model=JobCreated,
        status_code=status.HTTP_202_ACCEPTED,
        tags=["jobs"],
    )
    def submit_job(submission: JobSubmission, request: Request) -> JobCreated:
        address = client_ip(request, resolved.trusted_proxy_networks)
        try:
            access_token, job = service.submit(submission, address)
        except ReferenceUnavailableError as error:
            raise HTTPException(status_code=503, detail=str(error)) from error
        except FastaValidationError as error:
            raise HTTPException(status_code=422, detail=str(error)) from error
        except SubmissionRateLimitError as error:
            raise HTTPException(
                status_code=429,
                detail="submission rate limit exceeded; try again later",
                headers={"Retry-After": str(error.retry_after_seconds)},
            ) from error
        except ClientCapacityError as error:
            raise HTTPException(
                status_code=429,
                detail=(
                    "this network already has the maximum number of active jobs; "
                    "wait for one to finish or cancel it"
                ),
                headers={"Retry-After": "60"},
            ) from error
        except StorageCapacityError as error:
            raise HTTPException(
                status_code=503,
                detail="the service is temporarily below its free-space reserve",
                headers={"Retry-After": "300"},
            ) from error
        except QueueCapacityError as error:
            raise HTTPException(
                status_code=503,
                detail="the compute queue is full; try again later",
                headers={"Retry-After": "30"},
            ) from error
        return JobCreated(access_token=access_token, job=JobView.model_validate(job))

    @router.get("/api/v1/jobs/{job_id}", response_model=JobView, tags=["jobs"])
    def get_job(job_id: str, _token: str = Depends(capability)) -> JobView:
        try:
            return JobView.model_validate(service.get(job_id))
        except KeyError as error:
            raise HTTPException(status_code=404, detail="job not found") from error

    @router.delete("/api/v1/jobs/{job_id}", response_model=CancelResponse, tags=["jobs"])
    def cancel_job(job_id: str, _token: str = Depends(capability)) -> CancelResponse:
        try:
            return CancelResponse(job=JobView.model_validate(service.cancel(job_id)))
        except KeyError as error:
            raise HTTPException(status_code=404, detail="job not found") from error

    @router.get(
        "/api/v1/jobs/{job_id}/result",
        response_class=FileResponse,
        responses={
            200: {
                "description": "Complete reproducibility bundle; always an attachment.",
                "content": {"application/zip": {}},
            }
        },
        tags=["jobs"],
    )
    def download_result(
        job_id: str,
        _token: str = Depends(capability),
    ) -> FileResponse:
        try:
            path, artifact = service.result_bundle(job_id)
        except (KeyError, FileNotFoundError) as error:
            raise HTTPException(
                status_code=404, detail="result bundle not found"
            ) from error
        return FileResponse(
            path,
            media_type="application/zip",
            filename="results.zip",
            content_disposition_type="attachment",
        )

    @router.get(
        "/api/v1/jobs/{job_id}/artifacts/{artifact_id}",
        response_class=FileResponse,
        responses={
            200: {
                "description": "Allowlisted result file; always an attachment.",
                "content": {"application/octet-stream": {}},
            }
        },
        tags=["jobs"],
    )
    def download_artifact(
        job_id: str,
        artifact_id: str,
        _token: str = Depends(capability),
    ) -> FileResponse:
        try:
            path, artifact = service.artifact_file(job_id, artifact_id)
        except (KeyError, FileNotFoundError) as error:
            raise HTTPException(status_code=404, detail="artifact not found") from error
        # FileResponse always emits Content-Disposition: attachment when a
        # filename is supplied. HTML is not allowlisted at all.
        return FileResponse(
            path,
            media_type=artifact["media_type"],
            filename=Path(artifact["name"]).name,
            content_disposition_type="attachment",
        )

    return router
