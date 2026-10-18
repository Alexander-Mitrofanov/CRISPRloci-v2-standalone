"""Stable façade for safe public artifact discovery and publication."""

from __future__ import annotations

from .publication.archive import (
    _ArchiveLimitExceeded,
    _CappedArchiveWriter,
    _write_zip,
)
from .publication.collector import collect_artifacts
from .publication.contracts import (
    _BUNDLE_SUFFIXES,
    _DENIED_PARTS,
    _EMBEDDED_ABSOLUTE_PATH,
    _IDENTIFY_FILENAMES,
    _INDIVIDUAL_EXACT,
    _INTERNAL_JSON,
    _MAX_PUBLIC_JSON_BYTES,
    _MEDIA_TYPES,
    _OMIT,
    _PRIVATE_MANIFEST_KEYS,
    _STATIC_PRIVATE_ROOTS,
)
from .publication.discovery import _bundle_outputs, _public_text_safe
from .publication.safety import (
    _allowed_individual,
    _checkpoint,
    _contains_private_path,
    _fsync_published,
    _private_markers,
    _public_json_safe,
    _read_json_bounded,
    _record,
    _safe_file,
    _sanitize_manifest_value,
    _write_sanitized_manifest,
    file_sha256,
    resolve_registered_file,
)


__all__ = [
    "_ArchiveLimitExceeded",
    "_BUNDLE_SUFFIXES",
    "_CappedArchiveWriter",
    "_DENIED_PARTS",
    "_EMBEDDED_ABSOLUTE_PATH",
    "_IDENTIFY_FILENAMES",
    "_INDIVIDUAL_EXACT",
    "_INTERNAL_JSON",
    "_MAX_PUBLIC_JSON_BYTES",
    "_MEDIA_TYPES",
    "_OMIT",
    "_PRIVATE_MANIFEST_KEYS",
    "_STATIC_PRIVATE_ROOTS",
    "_allowed_individual",
    "_bundle_outputs",
    "_checkpoint",
    "_contains_private_path",
    "_fsync_published",
    "_private_markers",
    "_public_json_safe",
    "_public_text_safe",
    "_read_json_bounded",
    "_record",
    "_safe_file",
    "_sanitize_manifest_value",
    "_write_sanitized_manifest",
    "_write_zip",
    "collect_artifacts",
    "file_sha256",
    "resolve_registered_file",
]
