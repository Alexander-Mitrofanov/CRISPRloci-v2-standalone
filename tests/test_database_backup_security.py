from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import stat
import subprocess
import sys
from types import ModuleType, SimpleNamespace

import pytest


DEPLOY = Path(__file__).resolve().parents[1] / "deploy"


def load_module(name: str, filename: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, DEPLOY / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


backup = load_module("crispr_database_backup_security", "backup-database.py")


@pytest.fixture(autouse=True)
def local_root_policy(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(backup, "ROOT_UID", os.getuid())
    monkeypatch.setattr(backup, "ROOT_GID", os.getgid())
    # Tiny synthetic databases should not require production-sized disk reserves.
    # Capacity rejection tests inject their own statvfs evidence.
    monkeypatch.setattr(backup, "MIN_BACKUP_DESTINATION_FREE_BYTES", 1_000_000)


def private_directory(path: Path) -> Path:
    path.mkdir(mode=0o700)
    path.chmod(0o700)
    return path


def database(path: Path, value: str = "retained") -> Path:
    connection = sqlite3.connect(path)
    connection.execute("CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT)")
    connection.execute(
        "INSERT INTO metadata(name, value) VALUES ('schema_version', '7')"
    )
    connection.execute("CREATE TABLE jobs (id INTEGER PRIMARY KEY, value TEXT)")
    connection.execute("INSERT INTO jobs(value) VALUES (?)", (value,))
    connection.commit()
    connection.close()
    path.chmod(0o600)
    return path


def crash_left_wal_database(path: Path, value: str) -> Path:
    script = """
import os
import sqlite3
import sys

os.umask(0o077)
connection = sqlite3.connect(sys.argv[1])
assert connection.execute("PRAGMA journal_mode = WAL").fetchone() == ("wal",)
connection.execute("PRAGMA wal_autocheckpoint = 0")
connection.execute("CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT)")
connection.execute(
    "INSERT INTO metadata(name, value) VALUES ('schema_version', '7')"
)
connection.execute("CREATE TABLE jobs (value TEXT)")
connection.execute("INSERT INTO jobs(value) VALUES (?)", (sys.argv[2],))
connection.commit()
os.chmod(sys.argv[1], 0o600)
os._exit(0)
"""
    subprocess.run(
        [sys.executable, "-B", "-c", script, str(path), value],
        check=True,
    )
    assert path.is_file()
    assert Path(f"{path}-wal").is_file()
    return path


def backup_fixture(tmp_path: Path) -> tuple[Path, Path, dict[str, object]]:
    source_root = private_directory(tmp_path / "service")
    source = database(source_root / "queue.sqlite3")
    transaction = private_directory(tmp_path / "transaction-42")
    artifact = transaction / "queue.sqlite3"
    backup.create_backup(source, artifact, "release-before")
    document = json.loads(
        artifact.with_suffix(".sqlite3.json").read_text(encoding="utf-8")
    )
    return artifact, transaction, document


def restore_temporary_name(
    source: Path,
    destination: Path,
    staging: Path,
) -> str:
    return backup._transaction_name(
        "restore",
        staging,
        destination,
        str(source),
    )


def test_backup_manifest_schema_stays_compatible_and_restore_is_atomic(
    tmp_path: Path,
) -> None:
    artifact, staging, document = backup_fixture(tmp_path)
    assert set(document) == {
        "schema",
        "release_before_switch",
        "database_schema_version",
        "created_at_unix",
        "size_bytes",
        "sha256",
    }
    assert document["schema"] == {
        "name": "crispr-web-database-backup",
        "version": "1.0.0",
    }

    destination_root = private_directory(tmp_path / "live")
    destination = destination_root / "queue.sqlite3"
    backup.restore_backup(
        artifact,
        destination,
        staging_directory=staging,
        expected_uid=os.getuid(),
        expected_gid=os.getgid(),
        mode=0o600,
    )

    connection = sqlite3.connect(f"file:{destination}?mode=ro", uri=True)
    assert connection.execute("SELECT value FROM jobs").fetchall() == [("retained",)]
    connection.close()
    metadata = destination.stat()
    assert stat.S_IMODE(metadata.st_mode) == 0o600
    assert metadata.st_uid == os.getuid()
    assert metadata.st_gid == os.getgid()
    assert not any(
        Path(f"{destination}{suffix}").exists()
        for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )
    assert not tuple(staging.glob(".*.restore-*"))


def test_restore_cleans_stale_transaction_temp_and_all_sidecars(
    tmp_path: Path,
) -> None:
    artifact, staging, _ = backup_fixture(tmp_path)
    destination_root = private_directory(tmp_path / "live")
    destination = destination_root / "queue.sqlite3"
    temporary_name = restore_temporary_name(artifact, destination, staging)
    for path in (
        staging / temporary_name,
        *(
            staging / f"{temporary_name}{suffix}"
            for suffix in backup.SQLITE_SIDECAR_SUFFIXES
        ),
        *(Path(f"{destination}{suffix}") for suffix in backup.SQLITE_SIDECAR_SUFFIXES),
    ):
        path.write_bytes(b"stale")

    backup.restore_backup(
        artifact,
        destination,
        staging_directory=staging,
        expected_uid=os.getuid(),
        expected_gid=os.getgid(),
        mode=0o600,
    )

    assert destination.is_file()
    assert not (staging / temporary_name).exists()
    assert not any(
        (staging / f"{temporary_name}{suffix}").exists()
        for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )
    assert not any(
        Path(f"{destination}{suffix}").exists()
        for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )


def test_restore_cleans_deterministic_temp_even_when_manifest_is_malformed(
    tmp_path: Path,
) -> None:
    artifact, staging, _ = backup_fixture(tmp_path)
    destination_root = private_directory(tmp_path / "live")
    destination = destination_root / "queue.sqlite3"
    temporary_name = restore_temporary_name(artifact, destination, staging)
    staging_stale = (
        staging / temporary_name,
        *(
            staging / f"{temporary_name}{suffix}"
            for suffix in backup.SQLITE_SIDECAR_SUFFIXES
        ),
    )
    destination_sidecars = tuple(
        Path(f"{destination}{suffix}") for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )
    for path in (*staging_stale, *destination_sidecars):
        path.write_bytes(b"stale")
    manifest = artifact.with_suffix(".sqlite3.json")
    manifest.write_text("{malformed", encoding="utf-8")
    manifest.chmod(0o600)

    with pytest.raises(SystemExit, match="backup manifest is malformed"):
        backup.restore_backup(
            artifact,
            destination,
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )

    assert not any(path.exists() for path in staging_stale)
    assert not any(path.exists() for path in destination_sidecars)
    assert not destination.exists()


def test_malformed_restore_preserves_crash_left_live_wal_byte_for_byte(
    tmp_path: Path,
) -> None:
    artifact, staging, _ = backup_fixture(tmp_path)
    manifest = artifact.with_suffix(".sqlite3.json")
    manifest.write_text("{malformed", encoding="utf-8")
    manifest.chmod(0o600)
    destination_root = private_directory(tmp_path / "live")
    destination = crash_left_wal_database(
        destination_root / "queue.sqlite3", "live-wal-only"
    )
    live_paths = tuple(
        path
        for path in destination_root.iterdir()
        if path.name.startswith(destination.name)
    )
    before = {path.name: path.read_bytes() for path in live_paths}

    with pytest.raises(SystemExit, match="backup manifest is malformed"):
        backup.restore_backup(
            artifact,
            destination,
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )

    assert {path.name: path.read_bytes() for path in live_paths} == before
    connection = sqlite3.connect(f"file:{destination}?mode=ro", uri=True)
    assert connection.execute("SELECT value FROM jobs").fetchall() == [
        ("live-wal-only",)
    ]
    connection.close()


def test_restore_rejects_symlink_sidecar_without_touching_target(
    tmp_path: Path,
) -> None:
    artifact, staging, _ = backup_fixture(tmp_path)
    destination_root = private_directory(tmp_path / "live")
    destination = destination_root / "queue.sqlite3"
    victim = tmp_path / "victim"
    victim.write_text("do not touch", encoding="utf-8")
    sidecar = Path(f"{destination}-wal")
    sidecar.symlink_to(victim)

    with pytest.raises(SystemExit, match="unsafe .*SQLite sidecar"):
        backup.restore_backup(
            artifact,
            destination,
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )

    assert sidecar.is_symlink()
    assert victim.read_text(encoding="utf-8") == "do not touch"
    assert not destination.exists()


def test_restore_failure_removes_generated_temp_and_orphan_sidecars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    artifact, staging, _ = backup_fixture(tmp_path)
    destination_root = private_directory(tmp_path / "live")
    destination = destination_root / "queue.sqlite3"
    temporary_name = restore_temporary_name(artifact, destination, staging)
    original_integrity = backup.integrity
    integrity_calls = 0
    destination_sidecars = tuple(
        Path(f"{destination}{suffix}") for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )
    for path in destination_sidecars:
        path.write_bytes(b"preexisting live state")

    def fail_after_destination_creation(connection: sqlite3.Connection) -> None:
        nonlocal integrity_calls
        original_integrity(connection)
        integrity_calls += 1
        if integrity_calls == 2:
            for path in (
                *(
                    staging / f"{temporary_name}{suffix}"
                    for suffix in backup.SQLITE_SIDECAR_SUFFIXES
                ),
            ):
                path.write_bytes(b"generated during failure")
            backup.fail("injected restore integrity failure")

    monkeypatch.setattr(backup, "integrity", fail_after_destination_creation)
    with pytest.raises(SystemExit, match="injected restore integrity failure"):
        backup.restore_backup(
            artifact,
            destination,
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )

    assert not destination.exists()
    assert not (staging / temporary_name).exists()
    assert not any(
        (staging / f"{temporary_name}{suffix}").exists()
        for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )
    assert not any(path.exists() for path in destination_sidecars)


def test_restore_never_traverses_symlink_destination_parent(tmp_path: Path) -> None:
    artifact, staging, _ = backup_fixture(tmp_path)
    adversary = private_directory(tmp_path / "adversary")
    victim = adversary / "queue.sqlite3"
    victim.write_bytes(b"adversary bytes")
    alias = tmp_path / "live-alias"
    alias.symlink_to(adversary, target_is_directory=True)

    with pytest.raises(SystemExit, match="contains a symlink"):
        backup.restore_backup(
            artifact,
            alias / "queue.sqlite3",
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )

    assert victim.read_bytes() == b"adversary bytes"


def test_restore_rejects_service_parent_as_staging_area(tmp_path: Path) -> None:
    artifact, _, _ = backup_fixture(tmp_path)
    destination_root = private_directory(tmp_path / "live")
    staging = private_directory(destination_root / "attacker-controlled-parent")

    with pytest.raises(SystemExit, match="must not be beneath"):
        backup.restore_backup(
            artifact,
            destination_root / "queue.sqlite3",
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )


def test_restore_requires_exact_private_staging_metadata(tmp_path: Path) -> None:
    artifact, _, _ = backup_fixture(tmp_path)
    destination_root = private_directory(tmp_path / "live")
    staging = private_directory(tmp_path / "restore-stage")
    staging.chmod(0o750)

    with pytest.raises(SystemExit, match="exactly root:root mode 0700"):
        backup.restore_backup(
            artifact,
            destination_root / "queue.sqlite3",
            staging_directory=staging,
            expected_uid=os.getuid(),
            expected_gid=os.getgid(),
            mode=0o600,
        )


def test_backup_failure_removes_stale_and_generated_sidecars(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source_root = private_directory(tmp_path / "service")
    source = database(source_root / "queue.sqlite3")
    transaction = private_directory(tmp_path / "transaction-42")
    destination = transaction / "queue.sqlite3"
    temporary_name = backup._transaction_name(
        "backup", transaction, destination, "release-before"
    )
    stale = (
        transaction / temporary_name,
        transaction / f"{temporary_name}.json",
        *(
            transaction / f"{temporary_name}{suffix}"
            for suffix in backup.SQLITE_SIDECAR_SUFFIXES
        ),
        *(Path(f"{destination}{suffix}") for suffix in backup.SQLITE_SIDECAR_SUFFIXES),
    )
    for path in stale:
        path.write_bytes(b"stale")

    original_integrity = backup.integrity
    integrity_calls = 0

    def fail_after_destination_creation(connection: sqlite3.Connection) -> None:
        nonlocal integrity_calls
        original_integrity(connection)
        integrity_calls += 1
        if integrity_calls == 2:
            for path in (
                *(
                    transaction / f"{temporary_name}{suffix}"
                    for suffix in backup.SQLITE_SIDECAR_SUFFIXES
                ),
                *(
                    Path(f"{destination}{suffix}")
                    for suffix in backup.SQLITE_SIDECAR_SUFFIXES
                ),
            ):
                path.write_bytes(b"generated during failure")
            backup.fail("injected integrity failure")

    monkeypatch.setattr(backup, "integrity", fail_after_destination_creation)
    with pytest.raises(SystemExit, match="injected integrity failure"):
        backup.create_backup(source, destination, "release-before")

    assert not destination.exists()
    assert not destination.with_suffix(".sqlite3.json").exists()
    assert not any(path.exists() for path in stale)


def test_backup_checkpoints_live_wal_without_publishing_sidecars(
    tmp_path: Path,
) -> None:
    source_root = private_directory(tmp_path / "service")
    source = source_root / "queue.sqlite3"
    connection = sqlite3.connect(source)
    assert connection.execute("PRAGMA journal_mode = WAL").fetchone() == ("wal",)
    connection.execute("CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT)")
    connection.execute(
        "INSERT INTO metadata(name, value) VALUES ('schema_version', '7')"
    )
    connection.execute("CREATE TABLE jobs (value TEXT)")
    connection.execute("INSERT INTO jobs(value) VALUES ('from-wal')")
    connection.commit()
    source.chmod(0o600)
    assert Path(f"{source}-wal").exists()

    transaction = private_directory(tmp_path / "transaction-42")
    destination = transaction / "queue.sqlite3"
    backup.create_backup(source, destination, "release-before")
    connection.close()
    assert not any(
        Path(f"{destination}{suffix}").exists()
        for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )
    restored = sqlite3.connect(f"file:{destination}?mode=ro&immutable=1", uri=True)
    assert restored.execute("SELECT value FROM jobs").fetchall() == [("from-wal",)]
    restored.close()
    assert not any(
        Path(f"{destination}{suffix}").exists()
        for suffix in backup.SQLITE_SIDECAR_SUFFIXES
    )


def test_backup_rejects_symlink_sidecar_and_parent(tmp_path: Path) -> None:
    source_root = private_directory(tmp_path / "service")
    source = database(source_root / "queue.sqlite3")
    transaction = private_directory(tmp_path / "transaction-42")
    destination = transaction / "queue.sqlite3"
    victim = tmp_path / "victim"
    victim.write_text("safe", encoding="utf-8")
    Path(f"{destination}-journal").symlink_to(victim)

    with pytest.raises(SystemExit, match="unsafe transient or SQLite sidecar"):
        backup.create_backup(source, destination, "release-before")
    assert victim.read_text(encoding="utf-8") == "safe"

    Path(f"{destination}-journal").unlink()
    alias = tmp_path / "transaction-alias"
    alias.symlink_to(transaction, target_is_directory=True)
    with pytest.raises(SystemExit, match="contains a symlink"):
        backup.create_backup(source, alias / "queue.sqlite3", "release-before")
    assert not destination.exists()


def test_existing_backup_destination_keeps_crash_left_wal_on_rejection(
    tmp_path: Path,
) -> None:
    source_root = private_directory(tmp_path / "service")
    source = database(source_root / "queue.sqlite3", "new-source")
    transaction = private_directory(tmp_path / "transaction-42")
    destination = crash_left_wal_database(
        transaction / "queue.sqlite3", "existing-wal-only"
    )
    existing_paths = tuple(
        path for path in transaction.iterdir() if path.name.startswith(destination.name)
    )
    before = {path.name: path.read_bytes() for path in existing_paths}

    with pytest.raises(SystemExit, match="backup destination already exists"):
        backup.create_backup(source, destination, "release-before")

    assert {path.name: path.read_bytes() for path in existing_paths} == before
    connection = sqlite3.connect(f"file:{destination}?mode=ro", uri=True)
    assert connection.execute("SELECT value FROM jobs").fetchall() == [
        ("existing-wal-only",)
    ]
    connection.close()


def test_backup_refuses_source_above_byte_ceiling(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source_root = private_directory(tmp_path / "service")
    source = database(source_root / "queue.sqlite3")
    transaction = private_directory(tmp_path / "transaction-42")
    destination = transaction / "queue.sqlite3"
    monkeypatch.setattr(backup, "MAX_DATABASE_BACKUP_BYTES", 1)

    with pytest.raises(SystemExit, match="exceeds the backup byte ceiling"):
        backup.create_backup(source, destination, "release-before")

    assert not destination.exists()
    assert tuple(transaction.iterdir()) == ()


def test_backup_refuses_low_destination_bytes_and_inodes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source_root = private_directory(tmp_path / "service")
    source = database(source_root / "queue.sqlite3")
    transaction = private_directory(tmp_path / "transaction-42")
    destination = transaction / "queue.sqlite3"
    monkeypatch.setattr(
        backup.os,
        "fstatvfs",
        lambda _descriptor: SimpleNamespace(
            f_bavail=0,
            f_frsize=4_096,
            f_favail=0,
        ),
    )

    with pytest.raises(SystemExit, match="free-space reserve"):
        backup.create_backup(source, destination, "release-before")

    assert not destination.exists()
    assert tuple(transaction.iterdir()) == ()
