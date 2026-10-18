"""Trusted model-bundle validation and preloading for inference."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyhmmer

from casandra.utils import sha256_file


class ModelIntegrityError(ValueError):
    """Raised when a model bundle does not match its signed-by-hash manifest."""


def _json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ModelIntegrityError(f"Cannot read model metadata {path}: {error}") from error
    if not isinstance(value, dict):
        raise ModelIntegrityError(f"Model metadata must be a JSON object: {path}")
    return value


def _artifact_path(root: Path, relative: str) -> Path:
    candidate = (root / relative).resolve()
    resolved_root = root.resolve()
    if candidate != resolved_root and resolved_root not in candidate.parents:
        raise ModelIntegrityError(f"Artifact path escapes model directory: {relative!r}")
    return candidate


def verify_artifacts(
    root: Path,
    manifest: dict[str, Any],
    *,
    field: str = "artifacts",
) -> dict[str, Any]:
    """Verify every file declared by a model manifest."""
    entries = manifest.get(field)
    if not isinstance(entries, dict) or not entries:
        raise ModelIntegrityError(f"Model manifest has no {field!r} mapping: {root}")
    for relative, metadata in sorted(entries.items()):
        if not isinstance(relative, str) or not isinstance(metadata, dict):
            raise ModelIntegrityError(f"Invalid artifact entry in {root / 'manifest.json'}")
        path = _artifact_path(root, relative)
        if not path.is_file():
            raise ModelIntegrityError(f"Missing model artifact: {path}")
        expected_size = metadata.get("size")
        expected_hash = metadata.get("sha256")
        if not isinstance(expected_size, int) or not isinstance(expected_hash, str):
            raise ModelIntegrityError(f"Incomplete artifact metadata: {relative}")
        if path.stat().st_size != expected_size:
            raise ModelIntegrityError(
                f"Model artifact size mismatch for {relative}: "
                f"expected {expected_size}, found {path.stat().st_size}"
            )
        actual_hash = sha256_file(path)
        if actual_hash != expected_hash:
            raise ModelIntegrityError(
                f"Model artifact checksum mismatch for {relative}: "
                f"expected {expected_hash}, found {actual_hash}"
            )
    return manifest


def verify_bundle(bundle_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate the complete deployment bundle before any pickle is loaded."""
    bundle_dir = Path(bundle_dir)
    manifest = _json_object(bundle_dir / "manifest.json")
    verify_artifacts(bundle_dir, manifest)
    config = normalize_bundle_config(_json_object(bundle_dir / "config.json"))
    for key in ("protein_model_directory", "cassette_architecture_directory"):
        relative = config.get(key)
        if not isinstance(relative, str) or not _artifact_path(bundle_dir, relative).is_dir():
            raise ModelIntegrityError(f"Invalid bundle directory setting {key!r}")
    return manifest, config


def normalize_bundle_config(config: dict[str, Any]) -> dict[str, Any]:
    """Read the fail-closed successor contract without changing legacy bundles.

    Legacy runtimes require protein_model_directory at the top level. Its
    absence from a v2 bundle forces them to reject it before model loading.
    """
    if "specificity_runtime_contract" not in config:
        return config
    contract = config["specificity_runtime_contract"]
    if (not isinstance(contract, dict) or type(contract.get("schema_version")) is not int
            or contract["schema_version"] != 2 or contract.get("specificity_required") is not True):
        raise ModelIntegrityError("Unsupported required specificity runtime contract")
    if "protein_model_directory" in config:
        raise ModelIntegrityError("Ambiguous legacy and required-specificity model paths")
    relative = contract.get("protein_model_directory")
    if not isinstance(relative, str):
        raise ModelIntegrityError("Required specificity protein directory is missing")
    return {**config, "protein_model_directory": relative}


def normalize_protein_card(card: dict[str, Any]) -> dict[str, Any]:
    """Normalize required core settings only for a recognized overlay contract."""
    if "required_specificity_core" not in card:
        return card
    core = card["required_specificity_core"]
    if (not isinstance(core, dict) or type(core.get("schema_version")) is not int
            or core["schema_version"] != 2 or not isinstance(card.get("specificity_overlay"), dict)):
        raise ModelIntegrityError("Required specificity core contract is incomplete")
    if "report_evalue" in card or "detection_threshold" in card:
        raise ModelIntegrityError("Ambiguous legacy and required-specificity scan settings")
    if (type(core.get("report_evalue")) not in (int, float) or core["report_evalue"] != 1000.0
            or type(core.get("detection_threshold")) not in (int, float) or core["detection_threshold"] != 0.0):
        raise ModelIntegrityError("Required specificity core scan settings differ from full361")
    return {**card, "report_evalue": core["report_evalue"],
            "detection_threshold": core["detection_threshold"]}


def validate_threads(threads: int) -> int:
    if isinstance(threads, bool) or not isinstance(threads, int) or not 1 <= threads <= 64:
        raise ValueError("threads must be an integer between 1 and 64")
    return threads


@dataclass(frozen=True)
class ProteinModel:
    """A verified protein model kept in memory across analyses."""

    directory: Path
    manifest: dict[str, Any]
    card: dict[str, Any]
    metadata: dict[str, dict[str, Any]]
    hmms: tuple[pyhmmer.plan7.HMM, ...]
    specificity: Any = None

    @property
    def model_id(self) -> str:
        return str(self.manifest["model_id"])

    @classmethod
    def from_directory(cls, model_dir: Path, *, verify: bool = True) -> "ProteinModel":
        model_dir = Path(model_dir)
        manifest = _json_object(model_dir / "manifest.json")
        if verify:
            verify_artifacts(model_dir, manifest)
        card = normalize_protein_card(_json_object(model_dir / "model_card.json"))
        metadata = _json_object(model_dir / "profile_metadata.json")
        hmms: list[pyhmmer.plan7.HMM] = []
        with pyhmmer.plan7.HMMFile(model_dir / "profiles.hmm") as handle:
            while (hmm := handle.read()) is not None:
                hmms.append(hmm)
        if not hmms:
            raise ModelIntegrityError(f"Protein model contains no HMMs: {model_dir}")
        overlay = None
        if 'specificity_overlay' in card:
            from casandra.specificity import SpecificityOverlay
            configuration = card['specificity_overlay']
            if not isinstance(configuration, dict):
                raise ModelIntegrityError('Specificity overlay configuration must be an object')
            overlay = SpecificityOverlay.from_protein_directory(
                model_dir, configuration, metadata, card, hmms)
        return cls(model_dir, manifest, card, metadata, tuple(hmms), overlay)
