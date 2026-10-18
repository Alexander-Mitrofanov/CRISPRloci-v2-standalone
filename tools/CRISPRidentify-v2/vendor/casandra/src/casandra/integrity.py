"""Fail-closed validation for packaged CasAndra model bundles."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


class ModelIntegrityError(RuntimeError):
    """A model bundle is missing, corrupt, or internally inconsistent."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json_object(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ModelIntegrityError(
            "cannot read model metadata {}: {}".format(path, error)
        ) from error
    if not isinstance(value, dict):
        raise ModelIntegrityError("model metadata is not an object: {}".format(path))
    return value


def verify_artifact_manifest(directory: Path) -> dict:
    """Verify every runtime artifact declared by ``directory/manifest.json``."""

    directory = Path(directory).resolve()
    manifest_path = directory / "manifest.json"
    manifest = _read_json_object(manifest_path)
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        raise ModelIntegrityError(
            "model manifest has no runtime artifacts: {}".format(manifest_path)
        )
    for relative_name, expected in sorted(artifacts.items()):
        if not isinstance(relative_name, str) or not isinstance(expected, dict):
            raise ModelIntegrityError(
                "malformed artifact entry in {}".format(manifest_path)
            )
        candidate = (directory / relative_name).resolve()
        try:
            candidate.relative_to(directory)
        except ValueError as error:
            raise ModelIntegrityError(
                "artifact escapes model directory: {!r}".format(relative_name)
            ) from error
        if not candidate.is_file():
            raise ModelIntegrityError("missing model artifact: {}".format(candidate))
        if candidate.stat().st_size != expected.get("size"):
            raise ModelIntegrityError("model artifact size mismatch: {}".format(candidate))
        if _sha256(candidate) != expected.get("sha256"):
            raise ModelIntegrityError("model artifact checksum mismatch: {}".format(candidate))
    return manifest


def validate_model_bundle(model_bundle: Path) -> dict:
    """Validate the complete deployable bundle and its two nested models."""

    model_bundle = Path(model_bundle).resolve()
    manifest = verify_artifact_manifest(model_bundle)
    config = _read_json_object(model_bundle / "config.json")
    if config.get("schema_version") != 4 or manifest.get("schema_version") != 4:
        raise ModelIntegrityError("unsupported CasAndra model schema")
    if config.get("bundle_id") != manifest.get("bundle_id"):
        raise ModelIntegrityError("CasAndra bundle identifiers do not match")
    for key in ("protein_model_directory", "cassette_architecture_directory"):
        relative_name = config.get(key)
        if not isinstance(relative_name, str) or not relative_name:
            raise ModelIntegrityError("invalid model directory setting: {}".format(key))
        verify_artifact_manifest(model_bundle / relative_name)
    return manifest


__all__ = ("ModelIntegrityError", "validate_model_bundle", "verify_artifact_manifest")
