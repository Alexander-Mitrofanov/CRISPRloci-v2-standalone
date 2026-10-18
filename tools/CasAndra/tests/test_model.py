import hashlib
from pathlib import Path

import pytest

from casandra.api import DEFAULT_MODEL_BUNDLE
from casandra.model import ModelIntegrityError, verify_artifacts, verify_bundle


def test_bundled_model_manifest_is_complete_and_valid():
    manifest, config = verify_bundle(DEFAULT_MODEL_BUNDLE)
    assert manifest["bundle_id"] == config["bundle_id"]
    assert config["offline_inference"] is True
    assert config["hybrid_routing"] == (
        "Type II and Type III ordered-profile architecture ExtraTrees; "
        "direct profiles otherwise"
    )
    assert config["routing_policy_evidence"] == {
        "aggregate_accuracy_reestimated": False,
        "type_ii_basis": (
            "targeted canonical Cas9-Cas1-Cas2 regression using the packaged "
            "deployment-refit architecture model"
        ),
    }
    assert config["minimum_runtime_version"] == "0.3.0.dev2"
    assert config["specificity_runtime_contract"]["specificity_required"] is True
    assert manifest["bundle_id"] == "casandra-full361-v3-cap15-2026-09-20-v1"
    assert hashlib.sha256((DEFAULT_MODEL_BUNDLE / "manifest.json").read_bytes()).hexdigest() == (
        "966c4522d38585db6c3dbca14566aa25d520b7c165f487abe5f35845441eb1b3"
    )


def test_manifest_verification_rejects_tampering(tmp_path: Path):
    artifact = tmp_path / "model.bin"
    artifact.write_bytes(b"trusted")
    manifest = {
        "artifacts": {
            "model.bin": {
                "size": artifact.stat().st_size,
                "sha256": hashlib.sha256(b"trusted").hexdigest(),
            }
        }
    }
    verify_artifacts(tmp_path, manifest)
    artifact.write_bytes(b"corrupt")
    with pytest.raises(ModelIntegrityError, match="checksum mismatch"):
        verify_artifacts(tmp_path, manifest)


def test_manifest_verification_rejects_path_escape(tmp_path: Path):
    manifest = {"artifacts": {"../outside": {"size": 0, "sha256": "0" * 64}}}
    with pytest.raises(ModelIntegrityError, match="escapes"):
        verify_artifacts(tmp_path, manifest)
