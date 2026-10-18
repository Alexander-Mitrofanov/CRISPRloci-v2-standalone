"""Compatibility with a frozen predecessor contract; no old checkout required."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def assert_additive(before, after, path="contract"):
    if isinstance(before, dict):
        assert isinstance(after, dict), path
        for key, value in before.items():
            assert key in after, f"Removed {path}.{key}"
            assert_additive(value, after[key], f"{path}.{key}")
    elif isinstance(before, list):
        assert isinstance(after, list), path
        for value in before:
            assert value in after, f"Removed value in {path}: {value}"
    else:
        assert before == after, f"Changed {path}"


def test_predecessor_public_contract_is_preserved_as_an_additive_extension():
    before = json.loads((ROOT / "tests/fixtures/predecessor-api-v1.json").read_text())
    after = json.loads((ROOT / "frontend/src/contracts/public-api-v1.json").read_text())
    assert_additive(before, after)


def test_worker_keeps_predecessor_isolation_and_writable_temporary_storage():
    unit = (ROOT / "deploy/systemd/crisprloci-worker.service").read_text().splitlines()
    for directive in (
        "User=crisprsvc", "KillMode=control-group", "UMask=0077",
        "ProtectSystem=strict", "NoNewPrivileges=yes",
        "RestrictAddressFamilies=AF_UNIX", "IPAddressDeny=any",
        "ReadWritePaths=/srv/crispr/jobs/web",
    ):
        assert unit.count(directive) == 1
    mounts = [line for line in unit if line.startswith("TemporaryFileSystem=")]
    assert len(mounts) == 2
    assert all("mode=1777" in line and "nr_inodes=" in line and "size=" in line for line in mounts)
