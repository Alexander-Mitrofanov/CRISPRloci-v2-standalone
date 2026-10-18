"""Synthetic software checks; these are not biological validation."""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path

import pytest

from crisprleader.errors import ReferenceError
from crisprleader.references import load_bundle, reference_receipt, validate_manifest


def synthetic_manifest() -> dict:
    return {
        "schema_version": "1.0.0", "bundle_id": "synthetic-test-only", "version": "0.0.1",
        "families": [{
            "id": "synthetic-family", "repeat_sequences": ["ACGT"], "core_sequences": ["AACCGG"],
            "boundary_kind": "unresolved", "training_locus_ids": ["synthetic-locus"],
            "citations": ["synthetic:software-test"], "metadata": {"synthetic": True},
        }],
        "motifs": [],
        "provenance": {"build": {"purpose": "software test"}, "sources": [], "rights": [],
                       "limitations": ["No biological validation"]},
        "files": {}, "policy": {"production": False},
    }


def write_bundle(root: Path, manifest: dict | None = None) -> dict:
    manifest = synthetic_manifest() if manifest is None else manifest
    root.mkdir(parents=True, exist_ok=True)
    (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return manifest


def test_verified_member_receipt_and_retained_content(tmp_path: Path) -> None:
    root = tmp_path / "bundle"
    (root / "models").mkdir(parents=True)
    payload = b'{"engineering": true}\n'
    (root / "models/example.json").write_bytes(payload)
    manifest = synthetic_manifest()
    manifest["files"] = {"models/example.json": hashlib.sha256(payload).hexdigest()}
    manifest["families"][0].update({"model_file": "models/example.json", "model_format": "json"})
    write_bundle(root, manifest)
    bundle = load_bundle(root / "manifest.json")
    receipt = reference_receipt(bundle)
    assert receipt["manifest_sha256"] == hashlib.sha256((root / "manifest.json").read_bytes()).hexdigest()
    assert receipt["files"] == manifest["files"]
    assert bundle["_members"]["models/example.json"] == payload
    (root / "models/example.json").write_bytes(b"changed")
    assert bundle["_members"]["models/example.json"] == payload
    receipt["files"].clear()
    assert reference_receipt(bundle)["files"]


@pytest.mark.parametrize("missing", ["schema_version", "families", "provenance", "policy", "files"])
def test_missing_required_manifest_field(missing: str) -> None:
    manifest = synthetic_manifest()
    del manifest[missing]
    with pytest.raises(ReferenceError, match="keys invalid"):
        validate_manifest(manifest)


def test_unknown_fields_and_incompatible_schema_rejected() -> None:
    manifest = synthetic_manifest()
    manifest["unknown"] = True
    with pytest.raises(ReferenceError, match="unknown"):
        validate_manifest(manifest)
    del manifest["unknown"]
    manifest["schema_version"] = "99.0"
    with pytest.raises(ReferenceError, match="unsupported reference schema"):
        validate_manifest(manifest)


def test_duplicate_json_keys_and_nonfinite_values_rejected(tmp_path: Path) -> None:
    (tmp_path / "manifest.json").write_text('{"files":{}, "files":{}}')
    with pytest.raises(ReferenceError, match="duplicate JSON key"):
        load_bundle(tmp_path)
    (tmp_path / "manifest.json").write_text('{"x":NaN}')
    with pytest.raises(ReferenceError, match="nonfinite"):
        load_bundle(tmp_path)
    manifest = synthetic_manifest()
    manifest["provenance"]["build"]["unsafe_number"] = float("inf")
    with pytest.raises(ReferenceError, match="nonfinite"):
        validate_manifest(manifest)
    (tmp_path / "manifest.json").write_text('{"x":1e999}')
    with pytest.raises(ReferenceError, match="nonfinite"):
        load_bundle(tmp_path)


@pytest.mark.parametrize("member", ["../escape", "/absolute", "a/../b", "a//b", "a/./b",
                                    "a\\b", "a:b", "manifest.json", "code.pkl", "code.py"])
def test_unsafe_member_rejected(member: str) -> None:
    manifest = synthetic_manifest()
    manifest["files"][member] = "0" * 64
    with pytest.raises(ReferenceError, match="unsafe|executable"):
        validate_manifest(manifest)


def test_checksum_mismatch_fails_execution(tmp_path: Path) -> None:
    manifest = synthetic_manifest()
    manifest["files"] = {"model.json": "0" * 64}
    (tmp_path / "model.json").write_text("{}")
    write_bundle(tmp_path, manifest)
    with pytest.raises(ReferenceError, match="checksum mismatch"):
        load_bundle(tmp_path)


def test_missing_and_unlisted_members_rejected(tmp_path: Path) -> None:
    manifest = synthetic_manifest()
    manifest["files"] = {"missing.json": "0" * 64}
    write_bundle(tmp_path, manifest)
    with pytest.raises(ReferenceError, match="inventory mismatch"):
        load_bundle(tmp_path)
    write_bundle(tmp_path)
    (tmp_path / "unlisted.txt").write_text("unlisted")
    with pytest.raises(ReferenceError, match="inventory mismatch"):
        load_bundle(tmp_path)


def test_symlinked_file_directory_and_root_rejected(tmp_path: Path) -> None:
    root = tmp_path / "bundle"
    write_bundle(root)
    target = tmp_path / "outside"
    target.write_text("payload")
    (root / "link.txt").symlink_to(target)
    with pytest.raises(ReferenceError, match="symlink"):
        load_bundle(root)
    (root / "link.txt").unlink()
    (root / "linkdir").symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(ReferenceError, match="symlink"):
        load_bundle(root)
    alias = tmp_path / "alias"
    alias.symlink_to(root, target_is_directory=True)
    with pytest.raises(ReferenceError, match="unavailable"):
        load_bundle(alias)


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="FIFO unavailable")
def test_fifo_is_rejected_without_blocking(tmp_path: Path) -> None:
    write_bundle(tmp_path)
    os.mkfifo(tmp_path / "pipe")
    with pytest.raises(ReferenceError, match="nonregular"):
        load_bundle(tmp_path)


def test_size_limits_apply_before_unbounded_read(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import crisprleader.references as references

    write_bundle(tmp_path)
    monkeypatch.setattr(references, "MAX_MANIFEST_BYTES", 10)
    with pytest.raises(ReferenceError, match="byte limit"):
        load_bundle(tmp_path)


def test_missing_reference_is_setup_error(tmp_path: Path) -> None:
    with pytest.raises(ReferenceError, match="install a verified biological reference"):
        load_bundle(tmp_path / "missing")


def test_duplicate_family_and_bad_models_fail_closed() -> None:
    manifest = synthetic_manifest()
    manifest["families"].append(copy.deepcopy(manifest["families"][0]))
    with pytest.raises(ReferenceError, match="duplicate family"):
        validate_manifest(manifest)
    manifest["families"].pop()
    manifest["families"][0]["model_file"] = "model.json"
    with pytest.raises(ReferenceError, match="supplied together"):
        validate_manifest(manifest)
    manifest["families"][0]["model_format"] = "json"
    with pytest.raises(ReferenceError, match="not checksum-listed"):
        validate_manifest(manifest)


def test_malformed_json_is_reference_error(tmp_path: Path) -> None:
    (tmp_path / "manifest.json").write_bytes(b"\xff")
    with pytest.raises(ReferenceError, match="invalid reference manifest"):
        load_bundle(tmp_path)
