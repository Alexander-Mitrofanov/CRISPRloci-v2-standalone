from __future__ import annotations

import pytest

from crispr_web.config import _validated_scientific_manifest


def manifest(*optional, version="4.0.0"):
    names = ["crispridentify", "spacerplacer", "crispr_evor", "casandra", "crispr_tracrrna", *optional]
    identities = {
        name: {
            "name": name,
            "display_name": name,
            "version": "1.0.0",
            "kind": "git_source",
            "path": f"/private/sources/{name}",
            "sha256": "a" * 64,
            "sha256_scope": "source_inventory_v1",
            "repository": f"https://github.com/Alexander-Mitrofanov/{name}.git",
            "commit": "b" * 40,
        }
        for name in names
    }
    identities["integration"] = {
        "name": "integration",
        "display_name": "CRISPRloci integration",
        "version": "0.5.0",
        "kind": "wheel",
        "path": "/private/integration.whl",
        "sha256": "c" * 64,
        "sha256_scope": "file_bytes",
    }
    return {
        "schema": "crispr-web-scientific-release",
        "schema_version": version,
        "release_id": "d" * 64,
        "scientific_artifacts": identities,
    }


@pytest.mark.parametrize("optional", [("crisprrepeat",), ("crisprrepeat", "crisprleader")])
def test_repeat_manifest_with_or_without_leader_is_private_path_free(optional):
    identities = _validated_scientific_manifest(manifest(*optional))
    assert identities["crisprrepeat"]["commit"] == "b" * 40
    assert identities["integration"]["version"] == "0.5.0"
    assert ("crisprleader" in identities) is ("crisprleader" in optional)
    assert all("path" not in item for item in identities.values())


@pytest.mark.parametrize("optional", [(), ("crisprleader",), ("crisprrepeat", "unknown")])
def test_repeat_manifest_requires_repeat_and_rejects_unrecognized_tools(optional):
    with pytest.raises(ValueError, match="artifact set"):
        _validated_scientific_manifest(manifest(*optional))


@pytest.mark.parametrize(("version", "optional"), [("2.0.0", ()), ("3.0.0", ("crisprleader",))])
def test_existing_manifest_schemas_remain_valid(version, optional):
    assert _validated_scientific_manifest(manifest(*optional, version=version))


@pytest.mark.parametrize("version", ["2.0.0", "3.0.0"])
def test_repeat_requires_new_manifest_version(version):
    with pytest.raises(ValueError, match="artifact set"):
        _validated_scientific_manifest(manifest("crisprrepeat", version=version))


@pytest.mark.parametrize(("field", "value"), [("commit", "main"), ("sha256", "invalid"), ("kind", "wheel"), ("path", "relative/path")])
def test_repeat_identity_is_validated_fail_closed(field, value):
    document = manifest("crisprrepeat")
    document["scientific_artifacts"]["crisprrepeat"][field] = value
    with pytest.raises(ValueError):
        _validated_scientific_manifest(document)
