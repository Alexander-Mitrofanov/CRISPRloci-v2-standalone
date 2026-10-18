"""Composition of sanitized provenance, downloads, and a bounded result ZIP."""

from __future__ import annotations

import copy
from collections.abc import Callable, Iterable, Mapping
import json
import os
from pathlib import Path, PurePosixPath
from typing import Any

from .archive import _ArchiveLimitExceeded, _write_zip
from .discovery import _bundle_outputs, _public_text_safe
from .safety import (
    _allowed_individual,
    _checkpoint,
    _contains_private_path,
    _fsync_published,
    _private_markers,
    _record,
    _safe_file,
    _write_sanitized_manifest,
    file_sha256,
)


def _rewrite_provenance(path: Path, document: Mapping[str, Any]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def collect_artifacts(
    job_root: Path,
    output_root: Path,
    *,
    provenance: Mapping[str, Any],
    max_files: int,
    max_file_bytes: int,
    max_archive_bytes: int,
    private_roots: Iterable[str | Path] = (),
    checkpoint: Callable[[], None] | None = None,
) -> tuple[list[dict[str, object]], list[str]]:
    """Generate sanitized provenance, a bounded ZIP, and narrow downloads."""

    _checkpoint(checkpoint)
    markers = _private_markers((*private_roots, job_root, output_root))
    warnings: list[str] = []
    attempt_name = output_root.name
    exports = job_root / "exports" / attempt_name
    exports.mkdir(parents=True, mode=0o700, exist_ok=False)

    bundle_files: list[tuple[str, Path]] = []
    sanitized_manifest: Path | None = None
    raw_manifest = _safe_file(output_root, output_root / "adapter" / "manifest.json")
    if raw_manifest is not None:
        sanitized_manifest, redactions = _write_sanitized_manifest(
            raw_manifest,
            exports / "adapter-manifest.sanitized.json",
            markers=markers,
            max_bytes=max_file_bytes,
            checkpoint=checkpoint,
        )
        if sanitized_manifest is None:
            warnings.append(
                "Adapter manifest was omitted because it could not be safely sanitized."
            )
        else:
            bundle_files.append(("results/adapter/manifest.json", sanitized_manifest))
            if redactions:
                warnings.append(
                    f"Adapter manifest published with {redactions} private field(s) removed."
                )

    normalized_input = _safe_file(job_root, job_root / "input" / "normalized.fasta")
    if normalized_input is None:
        warnings.append("Normalized input was unavailable for the result bundle.")
    elif normalized_input.stat().st_size == 0:
        warnings.append("Normalized input was empty and omitted from the result bundle.")
    elif normalized_input.stat().st_size > max_file_bytes:
        warnings.append("Normalized input exceeds the bundle file limit.")
    else:
        bundle_files.insert(0, ("input/normalized.fasta", normalized_input))

    used_bytes = sum(path.stat().st_size for _, path in bundle_files)
    output_files, discovery_warnings = _bundle_outputs(
        output_root,
        max_file_bytes=max_file_bytes,
        max_files=max_files,
        max_total_bytes=max(0, max_archive_bytes - used_bytes),
        markers=markers,
        checkpoint=checkpoint,
    )
    bundle_files.extend(output_files)
    warnings.extend(discovery_warnings)

    used_bytes = sum(path.stat().st_size for _, path in bundle_files)
    launcher_log_root = job_root / "worker-logs"
    for suffix in ("stdout.log", "stderr.log"):
        log_path = launcher_log_root / f"{output_root.name}.{suffix}"
        safe_log = _safe_file(job_root, log_path)
        if safe_log is None:
            continue
        size = safe_log.stat().st_size
        if size == 0:
            warnings.append(f"Worker log omitted because it was empty: {log_path.name}")
            continue
        if size > max_file_bytes or used_bytes + size > max_archive_bytes:
            warnings.append(f"Worker log exceeds a bundle limit: {log_path.name}")
            continue
        try:
            log_is_safe = _public_text_safe(
                safe_log,
                markers=markers,
                max_bytes=max_file_bytes,
                checkpoint=checkpoint,
            )
        except OSError:
            log_is_safe = False
        if not log_is_safe:
            warnings.append(
                f"Worker log omitted because it contained a private path: {log_path.name}"
            )
            continue
        bundle_files.append((f"logs/{log_path.name}", safe_log))
        used_bytes += size

    file_manifest = [
        {
            "path": archive_name,
            "size_bytes": path.stat().st_size,
            "sha256": file_sha256(path, checkpoint),
        }
        for archive_name, path in bundle_files
    ]
    document = copy.deepcopy(dict(provenance))
    document["schema"] = {
        "name": "crisprloci-v2-provenance",
        "version": "1.0.0",
    }
    document["bundle_files"] = file_manifest
    if warnings:
        document["warnings"] = [*document.get("warnings", []), *warnings]
    if _contains_private_path(document, markers):
        raise ValueError("generated provenance contains a private filesystem path")

    provenance_path = exports / "provenance.json"
    with provenance_path.open("x", encoding="utf-8") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.chmod(provenance_path, 0o600)

    bundle_files.append(("provenance.json", provenance_path))
    total_size = sum(path.stat().st_size for _, path in bundle_files)
    archive_path: Path | None = None
    if total_size <= max_archive_bytes:
        temporary = exports / ".results.zip.tmp"
        archive_path = exports / "results.zip"
        try:
            try:
                _write_zip(
                    temporary,
                    bundle_files,
                    max_bytes=max_archive_bytes,
                    checkpoint=checkpoint,
                )
            except _ArchiveLimitExceeded:
                archive_path = None
                warning = (
                    "Result ZIP omitted because its final size exceeds the "
                    "configured archive limit."
                )
                warnings.append(warning)
                document["warnings"] = [*document.get("warnings", []), warning]
                _rewrite_provenance(provenance_path, document)
            else:
                os.chmod(temporary, 0o600)
                os.replace(temporary, archive_path)
        finally:
            temporary.unlink(missing_ok=True)
    else:
        warning = "Result ZIP omitted because safe outputs exceed the configured archive limit."
        warnings.append(warning)
        document["warnings"] = [*document.get("warnings", []), warning]
        _rewrite_provenance(provenance_path, document)

    records: list[dict[str, object]] = []
    published_paths: list[Path] = []
    for archive_name, safe in output_files:
        relative = PurePosixPath(archive_name.removeprefix("results/"))
        if _allowed_individual(relative):
            published_paths.append(safe)
            records.append(
                _record(
                    job_root,
                    safe,
                    relative.as_posix(),
                    checkpoint=checkpoint,
                )
            )
    if sanitized_manifest is not None:
        published_paths.append(sanitized_manifest)
        records.append(
            _record(
                job_root,
                sanitized_manifest,
                "adapter/manifest.json",
                checkpoint=checkpoint,
            )
        )
    published_paths.append(provenance_path)
    records.append(
        _record(job_root, provenance_path, "provenance.json", checkpoint=checkpoint)
    )
    if archive_path is not None:
        published_paths.append(archive_path)
        records.append(
            _record(job_root, archive_path, "results.zip", checkpoint=checkpoint)
        )
    _fsync_published(job_root, published_paths)
    return records, warnings
