"""Canonical FASTA input binding verification."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

from ..config import Settings
from ..fasta import FastaValidationError, normalize_fasta
from .files import _read_bounded_output, _safe_nonempty_output


def _report_source_bindings(
    job_root: Path,
    expected_source_ids: tuple[str, ...],
    *,
    max_bytes: int = 20_000_000,
) -> dict[str, tuple[int, str]] | None:
    records_root = job_root / "input" / "records"
    bindings: dict[str, tuple[int, str]] = {}
    try:
        for path in sorted(records_root.iterdir()):
            payload = _read_bounded_output(job_root, path, max_bytes=max_bytes)
            if payload is None:
                return None
            lines = payload.decode("ascii").splitlines()
            if not lines or not lines[0].startswith(">"):
                return None
            source_id = lines[0][1:]
            sequence = "".join(lines[1:]).upper()
            if (
                not source_id
                or not sequence
                or source_id in bindings
                or source_id not in expected_source_ids
            ):
                return None
            bindings[source_id] = (
                len(sequence),
                hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            )
    except (OSError, UnicodeError):
        return None
    return bindings if tuple(bindings) == expected_source_ids else None


def _verify_job_inputs(
    settings: Settings,
    job_root: Path,
    job: Mapping[str, Any],
) -> bool:
    """Verify the exact canonical FASTA files before and after untrusted tools."""

    normalized_path = job_root / "input" / "normalized.fasta"
    records_root = job_root / "input" / "records"
    if (
        not _safe_nonempty_output(job_root, normalized_path)
        or not records_root.is_dir()
        or records_root.is_symlink()
    ):
        return False
    try:
        raw = _read_bounded_output(
            job_root, normalized_path, max_bytes=settings.max_request_bytes
        )
        if raw is None:
            return False
        expected_digest = job.get("input_sha256")
        if not isinstance(expected_digest, str):
            return False
        if hashlib.sha256(raw).hexdigest() != expected_digest:
            return False
        text = raw.decode("ascii")
        normalized = normalize_fasta(
            text,
            max_total_bases=settings.max_total_bases,
            max_record_bases=settings.max_record_bases,
            max_records=settings.max_records,
            max_header_characters=settings.max_header_characters,
            molecule=job.get("molecule", "DNA"),
        )
        if normalized.combined != text or normalized.sha256 != expected_digest:
            return False
        source_ids = tuple(job.get("source_ids", ()))
        if normalized.source_ids != source_ids:
            return False
        actual_names = {path.name for path in records_root.iterdir()}
        expected_names = {record.filename for record in normalized.records}
        if actual_names != expected_names:
            return False
        for record in normalized.records:
            path = records_root / record.filename
            if not _safe_nonempty_output(job_root, path):
                return False
            wrapped = "\n".join(
                record.sequence[offset : offset + 80]
                for offset in range(0, len(record.sequence), 80)
            )
            expected = f">{record.identifier}\n{wrapped}\n".encode("ascii")
            if _read_bounded_output(job_root, path, max_bytes=len(expected)) != expected:
                return False
    except (OSError, UnicodeError, FastaValidationError, TypeError):
        return False
    return True

