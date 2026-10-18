from __future__ import annotations

import hashlib
import json

import pytest

from crispr_web.config import _origins


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("https://example.org", ("https://example.org",)),
        ("https://example.org:8443", ("https://example.org:8443",)),
        ("http://localhost:5173", ("http://localhost:5173",)),
        ("http://127.0.0.1:8000", ("http://127.0.0.1:8000",)),
        ("http://[::1]:8000", ("http://[::1]:8000",)),
    ],
)
def test_exact_canonical_cors_origins(raw, expected):
    assert _origins(raw) == expected


@pytest.mark.parametrize(
    "raw",
    [
        "*",
        "https://*.example.org",
        "http://example.org",
        "https://example.org/",
        "https://example.org/path",
        "https://example.org?query=yes",
        "https://user@example.org",
        "https://EXAMPLE.org",
        "https://example.org:443",
        " https://example.org",
        "https://example.org,https://example.org",
        "https://exa_mple.org",
    ],
)
def test_ambiguous_or_unsafe_cors_origins_are_rejected(raw):
    with pytest.raises(ValueError):
        _origins(raw)


def test_release_identifiers_are_required_exact_lowercase_sha256(
    monkeypatch, scientific_manifest
):
    from crispr_web.config import Settings

    manifest_path, manifest_digest, _document = scientific_manifest
    monkeypatch.setenv("CRISPR_RELEASE_MANIFEST", str(manifest_path))
    monkeypatch.delenv("CRISPR_WEB_RELEASE_ID", raising=False)
    monkeypatch.delenv("CRISPR_WEB_BACKEND_WHEEL_SHA256", raising=False)
    monkeypatch.delenv("CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256", raising=False)
    with pytest.raises(ValueError, match="CRISPR_WEB_RELEASE_ID is required"):
        Settings.from_env()

    monkeypatch.setenv("CRISPR_WEB_RELEASE_ID", "a" * 64)
    with pytest.raises(ValueError, match="CRISPR_WEB_BACKEND_WHEEL_SHA256 is required"):
        Settings.from_env()

    monkeypatch.setenv("CRISPR_WEB_BACKEND_WHEEL_SHA256", "B" * 64)
    with pytest.raises(ValueError, match="64 lowercase hexadecimal"):
        Settings.from_env()

    monkeypatch.setenv("CRISPR_WEB_BACKEND_WHEEL_SHA256", "b" * 64)
    with pytest.raises(
        ValueError, match="CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256 is required"
    ):
        Settings.from_env()

    monkeypatch.setenv("CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256", "C" * 64)
    with pytest.raises(ValueError, match="64 lowercase hexadecimal"):
        Settings.from_env()

    monkeypatch.setenv("CRISPR_WEB_SCIENTIFIC_RELEASE_MANIFEST_SHA256", manifest_digest)
    settings = Settings.from_env()
    assert settings.release_id == "a" * 64
    assert settings.backend_wheel_sha256 == "b" * 64
    assert settings.scientific_release_manifest_sha256 == manifest_digest
    assert set(settings.scientific_release_identities) == {
        "crispridentify",
        "spacerplacer",
        "crispr_evor",
        "integration",
    }
    assert all(
        "path" not in identity
        for identity in settings.scientific_release_identities.values()
    )


def test_manifest_exact_bytes_and_schema_are_fail_closed(
    settings_factory, scientific_manifest, tmp_path
):
    manifest_path, manifest_digest, document = scientific_manifest
    settings = settings_factory(
        scientific_release_manifest=manifest_path,
        scientific_release_manifest_sha256=manifest_digest,
    )
    assert settings.scientific_release_ok() is True

    manifest_path.write_bytes(manifest_path.read_bytes() + b" ")
    assert settings.scientific_release_ok() is False
    with pytest.raises(ValueError, match="SHA-256"):
        settings.verify_scientific_release()

    document.pop("scientific_artifacts")
    payload = (
        json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    malformed_path = tmp_path / "malformed-release.json"
    malformed_path.write_bytes(payload)
    with pytest.raises(ValueError, match="missing required"):
        settings_factory(
            scientific_release_manifest=malformed_path,
            scientific_release_manifest_sha256=hashlib.sha256(payload).hexdigest(),
        )


def test_runtime_file_layout_is_not_part_of_backend_identity_contract(
    settings_factory, scientific_manifest, tmp_path
):
    _manifest_path, _manifest_digest, document = scientific_manifest
    document["runtime_files"] = {"future-layout": {"opaque": True}}
    document["environment_locks"] = {"future": True}
    document["mafft"] = {"future": True}
    document["native_runtime"] = {"future": True}
    payload = (
        json.dumps(document, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")
    path = tmp_path / "future-runtime-layout.json"
    path.write_bytes(payload)
    settings = settings_factory(
        scientific_release_manifest=path,
        scientific_release_manifest_sha256=hashlib.sha256(payload).hexdigest(),
    )
    assert set(settings.scientific_release_identities) == {
        "crispridentify",
        "spacerplacer",
        "crispr_evor",
        "integration",
    }


def test_root_managed_stable_manifest_symlink_is_supported(
    settings_factory, scientific_manifest, tmp_path
):
    manifest_path, manifest_digest, _document = scientific_manifest
    stable_path = tmp_path / "stable-release-manifest.json"
    try:
        stable_path.symlink_to(manifest_path)
    except OSError:
        pytest.skip("symlinks unavailable on this filesystem")
    settings = settings_factory(
        scientific_release_manifest=stable_path,
        scientific_release_manifest_sha256=manifest_digest,
    )
    assert settings.scientific_release_ok() is True
