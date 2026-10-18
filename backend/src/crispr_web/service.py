"""Application service joining FASTA storage, queue state, and capability views."""

from __future__ import annotations

import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import stat
import time
from typing import Any

from .artifacts import file_sha256, resolve_registered_file
from .capacity import has_capacity
from .context_modes import context_mode_available
from .config import Settings
from .db import Store, utc_iso
from .exchange import sequence_input
from .fasta import (
    materialize_fasta,
    normalize_fasta,
    safe_upload_filename,
)
from .models import JobSubmission
from .security import new_access_token, token_digest


JOB_ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")
ARTIFACT_ID_PATTERN = re.compile(r"^[0-9a-f]{24}$")
ORPHAN_QUARANTINE_PATTERN = re.compile(r"^[0-9a-f]{32}[.][0-9a-f]{16}$")
_ORPHAN_SUBMISSION_GRACE_SECONDS = 3_600
_SUBMISSION_LOCK_NAME = ".submission.lock"


class StorageCapacityError(RuntimeError):
    pass


class ReferenceUnavailableError(RuntimeError):
    pass


CATEGORY_POLICY_CATEGORIES = {
    "bona_fide_possible": ("Bona-fide", "Possible"),
    "bona_fide_only": ("Bona-fide",),
}


def _remove_private_job_tree(root: Path, *, ignore_errors: bool = False) -> bool:
    """Remove one real directory tree without following a top-level symlink."""

    try:
        metadata = root.lstat()
    except FileNotFoundError:
        return True
    except OSError:
        if ignore_errors:
            return False
        raise
    if not stat.S_ISDIR(metadata.st_mode):
        if ignore_errors:
            return False
        raise OSError(f"private job tree is not a directory: {root}")
    try:
        for current, _directories, _files in os.walk(root, followlinks=False):
            Path(current).chmod(0o700)
        shutil.rmtree(root)
    except OSError:
        if not ignore_errors:
            raise
        return False
    return True


def require_unprivileged_process() -> None:
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        raise RuntimeError(
            "the CRISPR web API and worker must run as an unprivileged service account"
        )


class JobService:
    def __init__(self, settings: Settings, store: Store | None = None):
        self.settings = settings
        self.store = store or Store(settings)

    def job_root(self, job_id: str) -> Path:
        if not JOB_ID_PATTERN.fullmatch(job_id):
            raise KeyError(job_id)
        return self.settings.jobs_root / job_id

    def _purge_cancelled_job_data(self, job_id: str) -> bool:
        """Quarantine then remove data for a cancellation that won while queued."""

        root = self.job_root(job_id)
        quarantine = self.settings.data_root / ".orphan-quarantine"
        try:
            quarantine.mkdir(mode=0o700, exist_ok=True)
            metadata = quarantine.lstat()
        except OSError:
            return False
        if not stat.S_ISDIR(metadata.st_mode):
            return False
        os.chmod(quarantine, 0o700)
        target = quarantine / f"{job_id}.{os.urandom(8).hex()}"
        try:
            os.rename(root, target)
        except FileNotFoundError:
            return True
        except OSError:
            return False
        return _remove_private_job_tree(target, ignore_errors=True)

    def submit(
        self, submission: JobSubmission, client_ip: str
    ) -> tuple[str, dict[str, Any]]:
        if not context_mode_available(submission.mode.value, self.settings.scientific_release_identities):
            raise ReferenceUnavailableError("The requested context analysis is not provisioned in this release")
        if submission.mode.value == "spacer_association" and self.settings.association_reference_path is None:
            raise ReferenceUnavailableError("Association target/control panel is unavailable")
        if submission.mode.value == "repeat_type" and "repeattyper" not in self.settings.scientific_release_identities:
            raise ReferenceUnavailableError("RepeatTyper is not provisioned in this release")
        if submission.mode.value == "repeat_map" and self.settings.repeat_reference_path is None:
            raise ReferenceUnavailableError("Repeat mapping reference is unavailable")
        if submission.mode.value == "protospacer" and self.settings.spacer_database_path is None:
            raise ReferenceUnavailableError("Protospacer reference database is unavailable")
        if submission.mode.value == "viral_search" and self.settings.viral_reference_path is None:
            raise ReferenceUnavailableError("Viral reference database is unavailable")
        self.store.preflight_submission(client_ip)
        if not has_capacity(
            self.settings.data_root,
            required_bytes=self.settings.submission_admission_free_bytes,
            required_inodes=self.settings.submission_admission_free_inodes,
        ):
            raise StorageCapacityError("the service is below its byte or inode reserve")
        display_filename = safe_upload_filename(submission.filename)
        if submission.molecule == "RNA" and submission.mode.value not in {"repeats", "rna_family"}:
            from .fasta import FastaValidationError
            raise FastaValidationError("RNA input is supported only in repeat and RNA-family analyses")
        array_map = None
        if submission.mode.value == "array_compare":
            if "cctk" not in self.settings.scientific_release_identities:
                raise ReferenceUnavailableError("CCTK is not provisioned in this release")
            from .array_exchange import array_input
            input_text, array_map = array_input(submission.sequence, max_bytes=self.settings.max_request_bytes)
            handoff = None
        else:
            input_text, handoff = sequence_input(
                submission.sequence, mode=submission.mode.value, molecule=submission.molecule,
                max_bytes=self.settings.max_request_bytes,
            )
        context_limits = {"locus_compare": (8, 50_000), "prophage_context": (20, 2_000_000),
                          "insertion_context": (20, 2_000_000), "defense_context": (20, 2_000_000),
                          "rna_family": (50, 500)}
        record_limit, length_limit = context_limits.get(submission.mode.value, (self.settings.max_records, self.settings.max_record_bases))
        normalized = normalize_fasta(
            input_text,
            max_total_bases=min(self.settings.max_total_bases, 5_000_000) if submission.mode.value == "protospacer" else min(self.settings.max_total_bases, 2_000_000) if submission.mode.value in context_limits else self.settings.max_total_bases,
            max_record_bases=min(length_limit, (min(self.settings.max_record_bases, 80) if submission.mode.value == "viral_search" else min(self.settings.max_record_bases, 200) if submission.mode.value in {"repeats", "repeat_map", "repeat_type", "spacer_association"} else min(self.settings.max_record_bases, 5_000_000) if submission.mode.value == "protospacer" else self.settings.max_record_bases)),
            max_records=min(record_limit, min(self.settings.max_records, 100) if submission.mode.value == "repeat_map" else min(self.settings.max_records, 1000) if submission.mode.value in {"repeats", "repeat_context", "repeat_type", "spacer_association"} else self.settings.max_records),
            max_header_characters=self.settings.max_header_characters,
            fallback_identifier=Path(display_filename).stem,
            molecule=submission.molecule,
        )
        if submission.mode.value in {"locus_compare", "rna_family"} and normalized.record_count < 2:
            from .fasta import FastaValidationError
            raise FastaValidationError("This comparison requires at least two sequences")
        if submission.mode.value == "rna_family":
            from .fasta import FastaValidationError
            if any(len(record.sequence) < 30 or set(record.sequence) - set("ACGU" if submission.molecule == "RNA" else "ACGT") for record in normalized.records):
                raise FastaValidationError("RNA families require unambiguous sequences of 30–500 bases")
            if len({record.sequence for record in normalized.records}) < 2:
                raise FastaValidationError("RNA families require at least two distinct sequences")
        from .exchange import association_groups
        group_document = association_groups(normalized.source_ids,submission.association_grouping,handoff) if submission.mode.value == "spacer_association" else None
        effective_bias_corrections = (
            submission.bias_corrections and not self.settings.disable_bias_corrections
        )
        option_warnings = []
        if submission.bias_corrections and not effective_bias_corrections:
            option_warnings.append(
                "Bias corrections were disabled by the operator safety policy."
            )
        access_token = new_access_token()
        # A random 128-bit identifier is only a locator; the independent
        # 256-bit capability token is the authorization secret.
        job_id = os.urandom(16).hex()
        root = self.job_root(job_id)
        root.mkdir(mode=0o700)
        lock_path = root / _SUBMISSION_LOCK_NAME
        lock_flags = os.O_CREAT | os.O_EXCL | os.O_RDWR
        lock_flags |= getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        lock_fd = os.open(lock_path, lock_flags, 0o600)
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
        committed = False
        try:
            materialize_fasta(root, normalized)
            if handoff is not None:
                with (root / "handoff.json").open("x", encoding="utf-8") as stream:
                    json.dump(handoff, stream, allow_nan=False)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.chmod(root / "handoff.json", 0o600)
            if group_document is not None:
                with (root / "association-groups.json").open("x",encoding="utf8") as stream:
                    json.dump(group_document,stream,allow_nan=False)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.chmod(root / "association-groups.json",0o600)
            if array_map is not None:
                with (root / "array-map.json").open("x",encoding="utf8") as stream:
                    json.dump(array_map,stream,allow_nan=False)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.chmod(root / "array-map.json",0o600)
            self.store.create_job(
                job_id=job_id,
                access_token_digest=token_digest(access_token),
                mode=submission.mode.value,
                molecule=submission.molecule,
                tool_options=submission.tool_options,
                array_map_sha256=file_sha256(root/"array-map.json") if array_map is not None else None,
                association_reference_sha256=self.settings.association_reference_sha256 if submission.mode.value == "spacer_association" else None,
                association_grouping=submission.association_grouping if submission.mode.value == "spacer_association" else None,
                association_groups_sha256=file_sha256(root/"association-groups.json") if group_document is not None else None,
                repeat_reference_sha256=self.settings.repeat_reference_sha256 if submission.mode.value == "repeat_map" else None,
                viral_max_mismatches=submission.viral_max_mismatches,
                viral_reference_sha256=self.settings.viral_reference_sha256 if submission.mode.value == "viral_search" else None,
                spacer_database_sha256=self.settings.spacer_database_sha256 if submission.mode.value == "protospacer" else None,
                tracr_model_type=submission.tracr_model_type,
                leader_flank_length=submission.leader_flank_length,
                filename=display_filename,
                record_count=normalized.record_count,
                base_count=normalized.total_bases,
                input_sha256=normalized.sha256,
                source_ids=normalized.source_ids,
                category_policy=submission.category_policy.value,
                spacer_distance=submission.spacer_distance,
                bias_corrections_requested=submission.bias_corrections,
                bias_corrections_effective=effective_bias_corrections,
                option_warnings=option_warnings,
                client_ip=client_ip,
            )
            committed = True
        except BaseException:
            # If the database commit did not happen, the locked private tree is
            # safe to remove. If commit state cannot be read, leave it for the
            # conservative reconciliation pass instead of risking a live row.
            try:
                persisted = self.store.get_job(job_id)
            except Exception:
                persisted = object()
            if persisted is None:
                _remove_private_job_tree(root, ignore_errors=True)
            raise
        finally:
            if committed:
                try:
                    lock_path.unlink()
                except OSError:
                    pass
            fcntl.flock(lock_fd, fcntl.LOCK_UN)
            os.close(lock_fd)
        job = self.store.get_job(job_id)
        if job is None:
            raise RuntimeError("newly committed job is missing")
        return access_token, self.public_view(job)

    def authorize(self, job_id: str, token: str) -> bool:
        return bool(JOB_ID_PATTERN.fullmatch(job_id)) and self.store.authorize(
            job_id, token
        )

    def get(self, job_id: str) -> dict[str, Any]:
        if not JOB_ID_PATTERN.fullmatch(job_id):
            raise KeyError(job_id)
        job = self.store.get_job(job_id)
        if (
            job is None
            or job.get("expires_at") is None
            or float(job["expires_at"]) <= time.time()
        ):
            raise KeyError(job_id)
        return self.public_view(job)

    def cancel(self, job_id: str) -> dict[str, Any]:
        purge_data = self.store.request_cancel(job_id)
        if purge_data:
            self._purge_cancelled_job_data(job_id)
        return self.get(job_id)

    def artifact_file(
        self, job_id: str, artifact_id: str
    ) -> tuple[Path, dict[str, Any]]:
        if not JOB_ID_PATTERN.fullmatch(job_id) or not ARTIFACT_ID_PATTERN.fullmatch(
            artifact_id
        ):
            raise FileNotFoundError(artifact_id)
        job = self.store.get_job(job_id)
        if (
            job is None
            or job.get("expires_at") is None
            or float(job["expires_at"]) <= time.time()
            or job["status"]
            not in {"completed", "completed_no_eligible_groups", "failed"}
        ):
            raise FileNotFoundError(artifact_id)
        artifact = self.store.artifact(job_id, artifact_id)
        if artifact is None:
            raise FileNotFoundError(artifact_id)
        path = resolve_registered_file(self.job_root(job_id), artifact["relative_path"])
        if (
            path.stat().st_size != artifact["size_bytes"]
            or file_sha256(path) != artifact["sha256"]
        ):
            raise FileNotFoundError(artifact_id)
        return path, artifact

    def result_bundle(self, job_id: str) -> tuple[Path, dict[str, Any]]:
        job = self.store.get_job(job_id)
        if (
            job is None
            or job.get("expires_at") is None
            or float(job["expires_at"]) <= time.time()
        ):
            raise FileNotFoundError(job_id)
        for artifact in job.get("artifacts", []):
            if artifact["name"] == "results.zip":
                return self.artifact_file(job_id, artifact["artifact_id"])
        raise FileNotFoundError("results.zip")

    def reconcile_orphaned_job_trees(self, now: float | None = None) -> int:
        """Quarantine old job directories absent from the durable queue."""

        timestamp = time.time() if now is None else now
        cutoff = timestamp - _ORPHAN_SUBMISSION_GRACE_SECONDS
        quarantine = self.settings.data_root / ".orphan-quarantine"
        try:
            quarantine.mkdir(mode=0o700, exist_ok=True)
            quarantine_metadata = quarantine.lstat()
        except OSError:
            return 0
        if not stat.S_ISDIR(quarantine_metadata.st_mode):
            return 0
        os.chmod(quarantine, 0o700)

        # Finish deletions interrupted after an earlier atomic quarantine move.
        try:
            quarantined = tuple(quarantine.iterdir())
        except OSError:
            quarantined = ()
        for candidate in quarantined:
            if not ORPHAN_QUARANTINE_PATTERN.fullmatch(candidate.name):
                continue
            try:
                metadata = candidate.lstat()
            except OSError:
                continue
            if stat.S_ISDIR(metadata.st_mode):
                _remove_private_job_tree(candidate, ignore_errors=True)

        removed = 0
        try:
            candidates = tuple(self.settings.jobs_root.iterdir())
            # Avoid opening SQLite and deserializing every job once per retained
            # directory. A final point query below still closes commit races for
            # the small set of apparent orphans.
            registered_job_ids = self.store.registered_job_ids()
        except OSError:
            return 0
        for root in candidates:
            if not JOB_ID_PATTERN.fullmatch(root.name):
                continue
            try:
                metadata = root.lstat()
            except OSError:
                continue
            if not stat.S_ISDIR(metadata.st_mode) or metadata.st_mtime > cutoff:
                continue
            if root.name in registered_job_ids:
                continue

            lock_fd: int | None = None
            lock_path = root / _SUBMISSION_LOCK_NAME
            try:
                flags = os.O_RDWR | getattr(os, "O_CLOEXEC", 0)
                flags |= getattr(os, "O_NOFOLLOW", 0)
                try:
                    lock_fd = os.open(lock_path, flags)
                except FileNotFoundError:
                    lock_fd = None
                if lock_fd is not None:
                    try:
                        fcntl.flock(lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    except BlockingIOError:
                        continue
                try:
                    refreshed = root.lstat()
                except OSError:
                    continue
                if (
                    not stat.S_ISDIR(refreshed.st_mode)
                    or refreshed.st_mtime > cutoff
                    or self.store.get_job(root.name) is not None
                ):
                    continue
                target = quarantine / f"{root.name}.{os.urandom(8).hex()}"
                try:
                    os.rename(root, target)
                except OSError:
                    continue
                removed += 1
                _remove_private_job_tree(target, ignore_errors=True)
            except OSError:
                continue
            finally:
                if lock_fd is not None:
                    try:
                        fcntl.flock(lock_fd, fcntl.LOCK_UN)
                    except OSError:
                        pass
                    os.close(lock_fd)
        return removed

    def cleanup_expired(self) -> int:
        now = time.time()
        self.store.purge_expired_rate_events(now)
        # Retention must not depend on the worker being healthy. Recover leases
        # here as well as in the worker so a persistently unavailable worker
        # cannot strand a running row (and its private tree) past expiry.
        self.store.recover_stale(now)
        removed = self.reconcile_orphaned_job_trees(now)
        for job_id in self.store.expired_job_ids(now):
            if not JOB_ID_PATTERN.fullmatch(job_id):
                continue
            if not self.store.freeze_expired_job(job_id, now):
                continue
            root = self.job_root(job_id)
            if root.exists():
                try:
                    _remove_private_job_tree(root)
                except OSError:
                    continue
            if self.store.delete_expired_job(job_id, now):
                removed += 1
        return removed

    @staticmethod
    def public_view(job: dict[str, Any]) -> dict[str, Any]:
        artifacts = [
            {
                "artifact_id": artifact["artifact_id"],
                "name": artifact["name"],
                "size_bytes": artifact["size_bytes"],
                "sha256": artifact["sha256"],
                "media_type": artifact["media_type"],
                "download_url": (
                    f"/api/v1/jobs/{job['job_id']}/artifacts/{artifact['artifact_id']}"
                ),
            }
            for artifact in job.get("artifacts", [])
        ]
        error = None
        if job.get("error_code"):
            error = {
                "code": job["error_code"],
                "message": job.get("error_message")
                or "The scientific workflow failed.",
            }
        return {
            "job_id": job["job_id"],
            "status": job["status"],
            "mode": job["mode"],
            "created_at": utc_iso(job["created_at"]),
            "updated_at": utc_iso(job["updated_at"]),
            "started_at": utc_iso(job.get("started_at")),
            "finished_at": utc_iso(job.get("finished_at")),
            "expires_at": utc_iso(job.get("expires_at")),
            "queue_position": job.get("queue_position"),
            "cancel_requested": bool(job.get("cancel_requested")),
            "input": {
                "filename": job["filename"],
                "record_count": job["record_count"],
                "base_count": job["base_count"],
                "sha256": job["input_sha256"],
                "source_ids": job.get("source_ids", []),
            },
            "options": {
                "tool_options": job.get("tool_options", {}),
                "molecule": job.get("molecule", "DNA"),
                "association_grouping": job.get("association_grouping"),
                "viral_max_mismatches": job.get("viral_max_mismatches", 2),
                "leader_flank_length": job.get("leader_flank_length", 500),
                "tracr_model_type": job.get("tracr_model_type", "II"),
                "category_policy": job["category_policy"],
                "categories": list(CATEGORY_POLICY_CATEGORIES[job["category_policy"]]),
                "spacer_distance": job["spacer_distance"],
                "bias_corrections_requested": bool(job["bias_corrections_requested"]),
                "bias_corrections_effective": bool(job["bias_corrections_effective"]),
                "warnings": job.get("option_warnings", []),
            },
            "summary": job.get("summary"),
            "artifacts": artifacts,
            "error": error,
        }
