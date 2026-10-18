from __future__ import annotations

import io
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import check_private_data as privacy


def git(root, *args):
    return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True, text=True).stdout.strip()


@pytest.fixture
def repo(tmp_path):
    git(tmp_path, "init", "-q")
    git(tmp_path, "config", "user.name", "Guard test")
    git(tmp_path, "config", "user.email", "test@example.org")
    (tmp_path / ".gitignore").write_text(".ops-private/\noutput/\n.env*\n!.env.example\n")
    git(tmp_path, "add", ".gitignore")
    git(tmp_path, "commit", "-qm", "Initial")
    return tmp_path


def test_force_added_private_record_is_blocked_even_when_worktree_is_clean(repo, capsys):
    private = repo / ".ops-private"
    private.mkdir()
    (private / "connections.md").write_text("Restricted notes")
    assert privacy.check(repo, "workspace")
    git(repo, "add", "-f", ".ops-private/connections.md")
    assert not privacy.check(repo, "index")
    assert not privacy.check(repo, "workspace")
    assert "Restricted notes" not in capsys.readouterr().err


def test_index_uses_staged_blob_and_new_branch_checks_deleted_history(repo):
    text = "Host" + "Name private-host.example\n"
    (repo / "notes.md").write_text(text)
    git(repo, "add", "notes.md")
    (repo / "notes.md").write_text("Sanitized working copy")
    assert not privacy.check(repo, "index")
    git(repo, "commit", "-qm", "Unsafe staged copy")
    git(repo, "rm", "-f", "notes.md")
    git(repo, "commit", "-qm", "Deleted later")
    tip = git(repo, "rev-parse", "HEAD")
    stream = io.StringIO(f"refs/heads/main {tip} refs/heads/main {'0' * 40}\n")
    assert not privacy.check(repo, "pre-push", stream=stream)


def test_registered_vm_name_cannot_be_copied_into_frontend(repo, capsys):
    private = repo / ".ops-private"
    private.mkdir()
    (private / "forbidden-values.json").write_text('["synthetic-vm-marker"]')
    frontend = repo / "frontend"
    frontend.mkdir()
    (frontend / "synthetic-vm-marker.md").write_text("synthetic-vm-marker")
    assert not privacy.check(repo, "export", frontend)
    assert "synthetic-vm-marker" not in capsys.readouterr().err


def test_push_cannot_hide_data_in_a_tag_pointing_to_a_blob(repo):
    (repo / "notes.txt").write_text("private notes")
    oid = git(repo, "hash-object", "-w", "notes.txt")
    with pytest.raises(ValueError):
        privacy.check(repo, "pre-push", stream=io.StringIO(f"refs/tags/notes {oid} refs/tags/notes {'0' * 40}\n"))


def test_symlink_and_nested_private_files_are_blocked_from_export(repo):
    frontend = repo / "frontend"
    frontend.mkdir()
    (frontend / "notes").symlink_to(repo / ".git/config")
    assert not privacy.check(repo, "export", frontend)
    (frontend / "notes").unlink()
    (frontend / ".ops-private").mkdir()
    (frontend / ".ops-private/connections.md").write_text("private")
    assert not privacy.check(repo, "export", frontend)


def test_clean_examples_and_deleted_refs_are_allowed(repo):
    (repo / ".env.example").write_text("VITE_API_BASE_URL=https://analysis.example.org\n")
    git(repo, "add", ".env.example")
    assert privacy.check(repo, "index")
    tip = git(repo, "rev-parse", "HEAD")
    assert privacy.check(repo, "pre-push", stream=io.StringIO(f"(delete) {'0' * 40} refs/heads/main {tip}\n"))


def test_hooks_reject_commit_and_push_to_local_remote(repo, tmp_path):
    source = Path(__file__).resolve().parents[1]
    (repo / "scripts").mkdir()
    for name in ("check_private_data.py", "install_git_guards.py"):
        shutil.copyfile(source / "scripts" / name, repo / "scripts" / name)
    shutil.copytree(source / ".githooks", repo / ".githooks")
    subprocess.run([sys.executable, str(repo / "scripts/install_git_guards.py")], check=True)
    (repo / "notes.md").write_text("Identity" + "File ~/.ssh/synthetic-test-key\n")
    git(repo, "add", ".")
    result = subprocess.run(["git", "-C", str(repo), "commit", "-qm", "Blocked"], capture_output=True)
    assert result.returncode != 0
    git(repo, "-c", "core.hooksPath=/dev/null", "commit", "-qm", "Bypassed commit for push test")
    remote = tmp_path / "remote.git"
    subprocess.run(["git", "init", "--bare", "-q", str(remote)], check=True)
    git(repo, "remote", "add", "guard-test", str(remote))
    pushed = subprocess.run(["git", "-C", str(repo), "push", "guard-test", "HEAD:main"], capture_output=True)
    assert pushed.returncode != 0
    assert not (remote / "refs/heads/main").exists()


def test_inspection_failure_blocks_unknown_push_object(repo):
    with pytest.raises(ValueError):
        privacy.check(repo, "pre-push", stream=io.StringIO(f"refs/heads/main {'f' * 40} refs/heads/main {'0' * 40}\n"))


def test_export_ignores_local_source_config_but_rejects_it_at_destination(repo):
    frontend = repo / "frontend"
    frontend.mkdir()
    (frontend / ".env.local").write_text("LOCAL_ONLY=1\n")
    assert privacy.check(repo, "export-source", frontend)
    assert not privacy.check(repo, "export", frontend)


def test_actual_frontend_export_rejects_copied_vm_before_changing_target(tmp_path):
    source = Path(__file__).resolve().parents[1]
    project = tmp_path / "source"
    for name in ("scripts", "deploy", "frontend", ".ops-private"):
        (project / name).mkdir(parents=True)
    shutil.copyfile(source / "scripts/check_private_data.py", project / "scripts/check_private_data.py")
    exporter = project / "deploy/export-pages-frontend.sh"
    shutil.copyfile(source / "deploy/export-pages-frontend.sh", exporter)
    (project / ".ops-private/forbidden-values.json").write_text('["synthetic-export-vm"]')
    (project / "frontend/README.md").write_text("synthetic-export-vm")
    (project / "frontend/.env.local").write_text("PRIVATE_CONFIG=1\n")
    target = tmp_path / "target"
    target.mkdir()
    git(target, "init", "-q")
    git(target, "config", "user.name", "Export test")
    git(target, "config", "user.email", "test@example.org")
    origin = "https://github.com/example/frontend.git"
    git(target, "remote", "add", "origin", origin)
    (target / "keep-until-checked").write_text("unchanged")
    git(target, "add", ".")
    git(target, "commit", "-qm", "Initial")
    command = ["bash", str(exporter), str(target), origin]
    blocked = subprocess.run(command, capture_output=True)
    assert blocked.returncode != 0
    assert (target / "keep-until-checked").read_text() == "unchanged"
    (project / "frontend/README.md").write_text("Safe public frontend")
    passed = subprocess.run(command, capture_output=True)
    assert passed.returncode == 0, passed.stderr.decode()
    assert (target / "README.md").read_text() == "Safe public frontend"
    assert not (target / ".env.local").exists()
    assert not (target / "keep-until-checked").exists()


def test_installer_preserves_existing_default_hooks(repo):
    source = Path(__file__).resolve().parents[1]
    (repo / "scripts").mkdir()
    installer = repo / "scripts/install_git_guards.py"
    shutil.copyfile(source / "scripts/install_git_guards.py", installer)
    existing = repo / ".git/hooks/pre-commit"
    existing.write_text("#!/bin/sh\nexit 0\n")
    existing.chmod(0o755)
    result = subprocess.run([sys.executable, str(installer)], capture_output=True)
    assert result.returncode != 0
    assert b"Existing default Git hooks preserved" in result.stderr
    assert existing.read_text() == "#!/bin/sh\nexit 0\n"
    assert subprocess.run(["git", "-C", str(repo), "config", "--get", "core.hooksPath"], capture_output=True).returncode == 1


def test_source_archive_setup_reports_uninstalled_guards(tmp_path):
    source = Path(__file__).resolve().parents[1]
    (tmp_path / "scripts").mkdir()
    installer = tmp_path / "scripts/install_git_guards.py"
    shutil.copyfile(source / "scripts/install_git_guards.py", installer)
    strict = subprocess.run([sys.executable, str(installer)], capture_output=True)
    assert strict.returncode != 0
    archive = subprocess.run([sys.executable, str(installer), "--allow-source-archive"], capture_output=True)
    assert archive.returncode == 0
    assert b"NOT installed" in archive.stderr
