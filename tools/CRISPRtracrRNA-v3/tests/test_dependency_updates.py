"""Exercise the update workflow using local repositories, without network access."""

import importlib.util
import json
import subprocess
from copy import deepcopy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "dependency_manager", ROOT / "scripts/dependencies.py"
)
dependencies = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dependencies)


def git(path, *args):
    return subprocess.check_output(["git", "-C", str(path), *args], text=True).strip()


def commit(path, text):
    (path / "runtime.py").write_text(text)
    git(path, "add", "runtime.py")
    git(path, "-c", "user.name=Test", "-c", "user.email=test@example.invalid", "commit", "-m", text)
    return git(path, "rev-parse", "HEAD")


@pytest.fixture
def checkout(tmp_path, monkeypatch):
    upstream = tmp_path / "upstream"
    upstream.mkdir()
    git(upstream, "init", "-b", "main")
    first = commit(upstream, "v1")
    root = tmp_path / "project"
    root.mkdir()
    monkeypatch.setattr(dependencies, "ROOT", root)
    lock = {"CasAndra": {"url": str(upstream), "commit": first}}
    (root / "upstream-lock.json").write_text(json.dumps(lock) + "\n")
    dependencies.sync("CasAndra", lock["CasAndra"])
    return upstream, root, root / "external/CasAndra", lock


def test_independent_update_preserves_unmodified_sources_and_lock(checkout):
    upstream, root, target, lock = checkout
    first = lock["CasAndra"]["commit"]
    second = commit(upstream, "v2")
    dependencies.update("CasAndra", "main", lock)
    assert git(target, "rev-parse", "HEAD") == second
    assert not git(target, "status", "--porcelain")
    saved = json.loads((root / "upstream-lock.json").read_text())
    assert saved["CasAndra"]["previous_commit"] == first
    assert saved["CasAndra"]["commit"] == second
    (target / "runtime.py").write_text("local modification")
    with pytest.raises(RuntimeError, match="local changes"):
        dependencies.update("CasAndra", "main", lock)
    assert (target / "runtime.py").read_text() == "local modification"


def test_annotated_tag_records_commit_and_allows_later_sync_and_update(checkout):
    upstream, root, target, lock = checkout
    second = commit(upstream, "v2")
    git(
        upstream,
        "-c",
        "user.name=Test",
        "-c",
        "user.email=test@example.invalid",
        "tag",
        "-a",
        "release",
        "-m",
        "Release",
    )
    assert git(upstream, "rev-parse", "release") != second
    dependencies.update("CasAndra", "release", lock)
    assert lock["CasAndra"]["commit"] == git(target, "rev-parse", "HEAD") == second
    assert json.loads((root / "upstream-lock.json").read_text()) == lock
    dependencies.sync("CasAndra", lock["CasAndra"])
    third = commit(upstream, "v3")
    dependencies.update("CasAndra", "main", lock)
    assert lock["CasAndra"]["commit"] == git(target, "rev-parse", "HEAD") == third


@pytest.mark.parametrize("attached", [False, True])
@pytest.mark.parametrize("failure", ["staging", "replacement", "checkout"])
def test_failed_update_preserves_lock_checkout_and_owned_temporary_files(
    checkout, monkeypatch, attached, failure
):
    upstream, root, target, lock = checkout
    first = lock["CasAndra"]["commit"]
    second = commit(upstream, "v2")
    if attached:
        git(target, "checkout", "main")
    original = deepcopy(lock)
    lock_bytes = (root / "upstream-lock.json").read_bytes()
    unrelated = root / ".upstream-lock-unrelated.json"
    unrelated.write_text("belongs to someone else")

    if failure == "staging":

        def fail_write(value, handle, **kwargs):
            handle.write("partial JSON")
            raise OSError("injected staging failure")

        monkeypatch.setattr(dependencies.json, "dump", fail_write)
    elif failure == "replacement":

        def fail_replace(self, target):
            raise OSError("injected replacement failure")

        monkeypatch.setattr(Path, "replace", fail_replace)
    else:
        original_git = dependencies.git

        def fail_checkout(*args, **kwargs):
            result = original_git(*args, **kwargs)
            if args == ("checkout", "--detach", second):
                raise OSError("injected checkout failure after HEAD changed")
            return result

        monkeypatch.setattr(dependencies, "git", fail_checkout)

    with pytest.raises(OSError, match=f"injected {failure} failure"):
        dependencies.update("CasAndra", "main", lock)
    assert lock == original
    assert (root / "upstream-lock.json").read_bytes() == lock_bytes
    assert git(target, "rev-parse", "HEAD") == first
    assert git(target, "rev-parse", "--abbrev-ref", "HEAD") == ("main" if attached else "HEAD")
    assert not git(target, "status", "--porcelain")
    assert sorted(root.glob(".upstream-lock-*.json")) == [unrelated]
    assert unrelated.read_text() == "belongs to someone else"


def test_restoration_failure_reports_both_errors_and_recovery_state(checkout, monkeypatch):
    upstream, root, target, lock = checkout
    first = lock["CasAndra"]["commit"]
    second = commit(upstream, "v2")
    original = deepcopy(lock)
    lock_bytes = (root / "upstream-lock.json").read_bytes()
    original_git = dependencies.git

    def fail_restore(*args, **kwargs):
        if args == ("checkout", "--detach", first):
            raise OSError("injected restoration failure")
        return original_git(*args, **kwargs)

    def fail_replace(self, target):
        raise OSError("injected replacement failure")

    monkeypatch.setattr(dependencies, "git", fail_restore)
    monkeypatch.setattr(Path, "replace", fail_replace)
    with pytest.raises(RuntimeError) as error:
        dependencies.update("CasAndra", "main", lock)
    message = str(error.value)
    for detail in (
        "injected replacement failure",
        "injected restoration failure",
        str(target),
        first,
    ):
        assert detail in message
    assert "lock was not published" in message
    assert lock == original
    assert (root / "upstream-lock.json").read_bytes() == lock_bytes
    assert git(target, "rev-parse", "HEAD") == second
    assert not list(root.glob(".upstream-lock-*.json"))


@pytest.mark.parametrize("attached", [False, True])
@pytest.mark.parametrize("restoration_fails", [False, True])
def test_cleanup_failure_preserves_primary_error_and_recovery_diagnostics(
    checkout, monkeypatch, attached, restoration_fails
):
    upstream, root, target, lock = checkout
    first = lock["CasAndra"]["commit"]
    second = commit(upstream, "v2")
    if attached:
        git(target, "checkout", "main")
    original = deepcopy(lock)
    lock_bytes = (root / "upstream-lock.json").read_bytes()
    original_git = dependencies.git
    original_unlink = Path.unlink
    publication_error = OSError("injected publication failure")
    restore_command = (
        ("checkout", "--no-guess", "main")
        if attached
        else (
            "checkout",
            "--detach",
            first,
        )
    )

    def fail_replace(self, target):
        raise publication_error

    def fail_cleanup(self, *args, **kwargs):
        if self.name.startswith(".upstream-lock-"):
            raise PermissionError("injected cleanup failure")
        return original_unlink(self, *args, **kwargs)

    def conditional_restore(*args, **kwargs):
        if restoration_fails and args == restore_command:
            raise OSError("injected restoration failure")
        return original_git(*args, **kwargs)

    monkeypatch.setattr(Path, "replace", fail_replace)
    monkeypatch.setattr(Path, "unlink", fail_cleanup)
    monkeypatch.setattr(dependencies, "git", conditional_restore)
    with pytest.raises(RuntimeError if restoration_fails else OSError) as error:
        dependencies.update("CasAndra", "main", lock)
    if restoration_fails:
        assert error.value.__cause__ is publication_error
        for detail in ("injected publication failure", "injected restoration failure", first):
            assert detail in str(error.value)
    else:
        assert error.value is publication_error
    residue = list(root.glob(".upstream-lock-*.json"))
    assert len(residue) == 1
    assert error.value.__notes__ == [
        f"Could not remove staged lock {residue[0]}: injected cleanup failure"
    ]
    assert lock == original
    assert (root / "upstream-lock.json").read_bytes() == lock_bytes
    assert git(target, "rev-parse", "HEAD") == (second if restoration_fails else first)
    expected_branch = "main" if attached and not restoration_fails else "HEAD"
    assert git(target, "rev-parse", "--abbrev-ref", "HEAD") == expected_branch
    assert not git(target, "status", "--porcelain")


def test_successful_publication_does_not_remove_the_consumed_staging_path(checkout, monkeypatch):
    upstream, root, target, lock = checkout
    second = commit(upstream, "v2")
    original_unlink = Path.unlink
    cleanup_attempts = []

    def fail_cleanup(self, *args, **kwargs):
        if self.name.startswith(".upstream-lock-"):
            cleanup_attempts.append(self)
            raise PermissionError("injected cleanup failure")
        return original_unlink(self, *args, **kwargs)

    monkeypatch.setattr(Path, "unlink", fail_cleanup)
    dependencies.update("CasAndra", "main", lock)
    assert cleanup_attempts == []
    assert json.loads((root / "upstream-lock.json").read_text()) == lock
    assert lock["CasAndra"]["commit"] == git(target, "rev-parse", "HEAD") == second
    assert not list(root.glob(".upstream-lock-*.json"))


@pytest.mark.parametrize("ref", ["", " ", "--all", "main:refs/heads/main"])
def test_unsupported_ref_is_rejected_before_git(ref, monkeypatch):
    def unexpected_git(*args, **kwargs):
        pytest.fail("invalid refs must not invoke Git")

    monkeypatch.setattr(dependencies, "git", unexpected_git)
    with pytest.raises(ValueError, match="Git tag, branch, or commit"):
        dependencies.update("CasAndra", ref, {})


@pytest.mark.parametrize("mismatch", ["remote", "commit"])
def test_update_refuses_remote_or_head_mismatch_before_fetch(checkout, monkeypatch, mismatch):
    upstream, root, target, lock = checkout
    if mismatch == "remote":
        lock["CasAndra"]["url"] = str(upstream / "other")
    else:
        lock["CasAndra"]["commit"] = "0" * 40
    original = deepcopy(lock)
    lock_bytes = (root / "upstream-lock.json").read_bytes()
    original_git = dependencies.git

    def no_fetch(*args, **kwargs):
        assert args[0] != "fetch"
        return original_git(*args, **kwargs)

    monkeypatch.setattr(dependencies, "git", no_fetch)
    with pytest.raises(RuntimeError, match="remote|differs from lock"):
        dependencies.update("CasAndra", "main", lock)
    assert lock == original
    assert (root / "upstream-lock.json").read_bytes() == lock_bytes
    assert not git(target, "status", "--porcelain")


@pytest.mark.parametrize("failure", ["clone", "checkout"])
def test_failed_sync_leaves_no_destination_and_corrected_retry_succeeds(
    tmp_path, monkeypatch, failure
):
    upstream = tmp_path / "upstream"
    upstream.mkdir()
    git(upstream, "init", "-b", "main")
    first = commit(upstream, "v1")
    root = tmp_path / "project"
    monkeypatch.setattr(dependencies, "ROOT", root)
    good = {"url": str(upstream), "commit": first}
    bad = dict(good)
    bad["url" if failure == "clone" else "commit"] = "nonexistent"
    with pytest.raises(subprocess.CalledProcessError):
        dependencies.sync("CasAndra", bad)
    assert not list((root / "external").iterdir())
    dependencies.sync("CasAndra", good)
    assert git(root / "external/CasAndra", "rev-parse", "HEAD") == first
    assert list((root / "external").iterdir()) == [root / "external/CasAndra"]


def test_sync_preserves_preexisting_directory(tmp_path, monkeypatch):
    monkeypatch.setattr(dependencies, "ROOT", tmp_path)
    target = tmp_path / "external/CasAndra"
    target.mkdir(parents=True)
    sentinel = target / "keep.txt"
    sentinel.write_text("original contents")
    with pytest.raises(subprocess.CalledProcessError):
        dependencies.sync("CasAndra", {"url": "unused", "commit": "unused"})
    assert list(target.iterdir()) == [sentinel]
    assert sentinel.read_text() == "original contents"


def test_sync_keeps_untracked_build_products_in_existing_pinned_checkout(checkout):
    _, _, target, lock = checkout
    product = target / "build/native.so"
    product.parent.mkdir()
    product.write_bytes(b"build product")
    dependencies.sync("CasAndra", lock["CasAndra"])
    assert product.read_bytes() == b"build product"
    assert git(target, "rev-parse", "HEAD") == lock["CasAndra"]["commit"]
