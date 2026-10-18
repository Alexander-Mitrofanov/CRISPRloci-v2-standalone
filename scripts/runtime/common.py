"""Local operator helpers for immutable source and runtime releases."""

from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[2]
RUNTIME = ROOT / ".runtime"


def run(command, *, cwd=None, env=None):
    print("+", " ".join(str(part) for part in command), flush=True)
    subprocess.run([str(part) for part in command], cwd=cwd, env=env, check=True)


def digest(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(1024 * 1024):
            result.update(block)
    return result.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    os.replace(temporary, path)


def load_lock(path):
    data = json.loads(Path(path).read_text())
    import sys

    sys.path.insert(0, str(ROOT / "scripts"))
    from sync_tools import validate_entry

    if data.get("schema_version") != 1 or len(data.get("tools", [])) not in set(range(5, 18)):
        raise ValueError("Expected a supported five- through seventeen-tool lock schema")
    names = {entry["name"] for entry in data["tools"]}
    baseline = {
        "CRISPRidentify-v2",
        "CasAndra",
        "CRISPRtracrRNA-v3",
        "SpacerPlacer",
        "evOr",
    }
    optional = {"CRISPRleader-v2", "CRISPRrepeat", "CRISPRspacer", "CRISPRmap-v2", "RepeatTyper", "SpacePHARER", "CCTK", "clinker-v2", "Phigaro-v2", "ISEScan-v2", "DefenseFinder-v2", "RNAalifold-v2"}
    if not baseline <= names <= baseline | optional or len(names) != len(data["tools"]):
        raise ValueError("Tool lock must contain the supported tool set without duplicates")
    for entry in data["tools"]:
        validate_entry(entry)
        if (
            entry["repository"]
            != f"https://github.com/Alexander-Mitrofanov/{entry['name']}.git"
        ):
            raise ValueError("Unexpected repository identity")
    return data


def lock_id(lock):
    return hashlib.sha256(
        json.dumps(lock, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def export_repository(entry, destination):
    from bundled_tools import export_bundled

    if export_bundled(ROOT, entry, destination):
        return
    repository = ROOT / entry["path"]
    if (
        subprocess.check_output(
            ["git", "-C", str(repository), "config", "--get", "remote.origin.url"],
            text=True,
        ).strip()
        != entry["repository"]
    ):
        raise ValueError("Tool checkout origin does not match the lock")
    with tempfile.TemporaryFile() as archive:
        subprocess.run(
            ["git", "-C", str(repository), "archive", "--format=tar", entry["commit"]],
            stdout=archive,
            check=True,
        )
        archive.seek(0)
        with tarfile.open(fileobj=archive) as tar:
            for member in tar.getmembers():
                if (
                    not (member.isfile() or member.isdir())
                    or Path(member.name).is_absolute()
                    or ".." in Path(member.name).parts
                ):
                    raise ValueError("Git export contains a link or unsafe entry")
            tar.extractall(destination, filter="data")


def source_inventory(path):
    excluded = {
        ".git",
        ".venv",
        ".envs",
        "node_modules",
        "__pycache__",
        ".pytest_cache",
        ".ruff_cache",
        "build",
        "dist",
    }
    files = {}
    for file in sorted(Path(path).rglob("*")):
        relative = file.relative_to(path)
        # zstd ships essential source CMake recipes under build/. Include only
        # this pinned SpacePHARER vendored source, preserving old tool receipts.
        vendored_recipe = Path(path).name == "SpacePHARER" and relative.parts[:5] == ("lib", "mmseqs", "lib", "zstd", "build")
        if any(
            (part in excluded and not (part == "build" and vendored_recipe)) or part.endswith(".egg-info") for part in relative.parts
        ):
            continue
        if file.name.startswith(".env") and file.name != ".env.example":
            continue
        if (
            file.is_file()
            and not file.is_symlink()
            and file.suffix not in {".pyc", ".pyo"}
        ):
            files[str(relative)] = digest(file)
    return files


def environment_inventory(path):
    """Hash installed Python/native code and models; bytecode caches are expendable.

    Symlink destinations and executable permissions are part of the receipt.
    This detects local environment edits as well as changes to source snapshots.
    """
    files = {}
    for directory, names, filenames in os.walk(path):
        names[:] = sorted(name for name in names if name != "__pycache__")
        for name in sorted(
            [
                *filenames,
                *(name for name in names if (Path(directory) / name).is_symlink()),
            ]
        ):
            file = Path(directory) / name
            if file.suffix in {".pyc", ".pyo"}:
                continue
            entry = {"executable": bool(file.stat().st_mode & 0o111)}
            if file.is_symlink():
                entry["link"] = os.readlink(file)
            if file.is_file():
                entry["sha256"] = digest(file)
            files[str(file.relative_to(path))] = entry
    return files
