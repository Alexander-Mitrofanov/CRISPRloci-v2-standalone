"""Operator-only pinned native/model provisioning for locus-context tools."""

from __future__ import annotations

import json
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen

from .common import atomic_json, digest
from .progress import emit

DEFENSE_COMMIT = "5bdb41257fc7234fb8c3b70e3d362c52541860dd"
DEFENSE_MANIFEST = "bb4dc6de8dfce9179e5d14b1649d1af1e27ca5b0e14840787e883aa46a8d7348"
PVOG_URL = "https://ftp.ncbi.nlm.nih.gov/pub/kristensen/pVOGs/downloads/All/AllvogHMMprofiles.tar.gz"
PVOG_ARCHIVE = "1c5bd9ec035fcae466760ab58c3e2dbd088e12cf1de608e8a1074b8b98bda5a9"
PVOG_HMM = "0ddeafff975f1cfd86686e95f52c7a9980f9ad10fa082797313097ca42be97ed"
PVOG_COUNT = 9518
PVOG_BYTES = 767855419


def download(url, destination, *, limit, sha256=None):
    """Bounded download to a new operator staging file; never called by a job."""
    destination = Path(destination)
    temporary = destination.with_name(destination.name + ".partial")
    try:
        with urlopen(url, timeout=120) as response, temporary.open("xb") as stream:
            total = 0
            expected = int(response.headers.get("Content-Length", 0)) or None
            emit("Downloading " + destination.name, 0, expected, "bytes")
            while block := response.read(1024 * 1024):
                total += len(block)
                if total > limit:
                    raise ValueError("Reference archive exceeds provisioning limit")
                stream.write(block)
                emit("Downloading " + destination.name, total, expected, "bytes")
        if sha256 and digest(temporary) != sha256:
            raise ValueError("Reference archive differs from pinned identity")
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def safe_member(member):
    path = PurePosixPath(member.name)
    if (
        path.is_absolute()
        or ".." in path.parts
        or not (member.isfile() or member.isdir())
    ):
        raise ValueError("Unsafe reference archive member")
    return path


def install_defense_models(source, prefix):
    manifest = source / "defense_context/models.lock.json"
    if digest(manifest) != DEFENSE_MANIFEST:
        raise ValueError("DefenseFinder model manifest is not the reviewed release")
    expected = json.loads(manifest.read_text())["files"]
    parent = prefix / "share/defensefinder"
    destination = parent / "defense-finder-models"
    if not destination.exists():
        parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".defense-models-", dir=parent) as tmp:
            tmp = Path(tmp)
            archive = tmp / "models.tar.gz"
            staging = tmp / "defense-finder-models"
            staging.mkdir()
            download(
                f"https://codeload.github.com/mdmparis/defense-finder-models/tar.gz/{DEFENSE_COMMIT}",
                archive,
                limit=200_000_000,
            )
            seen = set()
            total = 0
            with tarfile.open(archive, "r|gz") as package:
                for member in package:
                    name = safe_member(member)
                    relative = "/".join(name.parts[1:])
                    if not member.isfile() or relative not in expected:
                        continue
                    if relative in seen or not 0 < member.size <= 20_000_000:
                        raise ValueError("Invalid DefenseFinder model member")
                    total += member.size
                    if total > 1_000_000_000:
                        raise ValueError(
                            "DefenseFinder models exceed provisioning budget"
                        )
                    file = staging / relative
                    file.parent.mkdir(parents=True, exist_ok=True)
                    with package.extractfile(member) as data, file.open("xb") as output:
                        shutil.copyfileobj(data, output, length=1024 * 1024)
                    if digest(file) != expected[relative]:
                        raise ValueError(
                            "DefenseFinder model bytes differ from pinned manifest"
                        )
                    seen.add(relative)
            if seen != set(expected):
                raise ValueError("Incomplete DefenseFinder archive")
            staging.replace(destination)
    actual = {
        p.relative_to(destination).as_posix(): digest(p)
        for p in destination.rglob("*")
        if p.is_file() and not p.is_symlink()
    }
    if actual != expected or any(p.is_symlink() for p in destination.rglob("*")):
        raise ValueError("Installed DefenseFinder models differ from pinned manifest")
    return parent, {
        "version": "1.2.3",
        "commit": DEFENSE_COMMIT,
        "manifest_sha256": DEFENSE_MANIFEST,
        "path": str(destination),
        "file_count": len(expected),
    }


def concatenate_pvogs(archive, destination):
    """Read gzip once, then concatenate sorted safe members without quadratic seeks."""
    if digest(archive) != PVOG_ARCHIVE:
        raise ValueError("pVOG archive checksum mismatch")
    with tempfile.TemporaryDirectory(
        prefix=".pvogs-", dir=destination.parent
    ) as temporary:
        temporary = Path(temporary)
        names = set()
        total = 0
        with tarfile.open(archive, "r|gz") as package:
            for member in package:
                path = safe_member(member)
                if member.isdir():
                    continue
                if not re.fullmatch(
                    r"AllvogHMMprofiles/VOG[0-9]+[a-z]?\.hmm", str(path)
                ):
                    raise ValueError("Unexpected pVOG archive member")
                if member.name in names or not 0 < member.size <= 20_000_000:
                    raise ValueError("Duplicate or oversized pVOG profile")
                total += member.size
                if total > PVOG_BYTES or len(names) >= PVOG_COUNT:
                    raise ValueError("pVOG model inventory exceeds pinned limits")
                names.add(member.name)
                target = temporary / path.name
                with package.extractfile(member) as data, target.open("xb") as output:
                    shutil.copyfileobj(data, output, length=1024 * 1024)
        if len(names) != PVOG_COUNT or total != PVOG_BYTES:
            raise ValueError("pVOG model inventory differs from pinned release")
        combined = temporary / "pvogs.hmm"
        with combined.open("xb") as output:
            for name in sorted(names):
                with (temporary / PurePosixPath(name).name).open("rb") as data:
                    shutil.copyfileobj(data, output, length=1024 * 1024)
        if digest(combined) != PVOG_HMM:
            raise ValueError("pVOG concatenation identity mismatch")
        combined.replace(destination)


def install_pvogs(prefix):
    folder = prefix / "share/phigaro"
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / "pvogs.hmm"
    if not target.exists():
        archive = folder / ".profiles.tar.gz"
        try:
            download(PVOG_URL, archive, limit=200_000_000, sha256=PVOG_ARCHIVE)
            concatenate_pvogs(archive, target)
        finally:
            archive.unlink(missing_ok=True)
    if target.is_symlink() or digest(target) != PVOG_HMM:
        raise ValueError("Installed pVOG identity mismatch")
    receipt = {
        "name": "pVOGs",
        "url": PVOG_URL,
        "archive_sha256": PVOG_ARCHIVE,
        "sha256": PVOG_HMM,
        "profiles": PVOG_COUNT,
        "bytes": PVOG_BYTES,
        "path": str(target),
    }
    atomic_json(folder / "reference.json", receipt)
    return target, receipt


def native_identity(path, version, arguments, marker, env):
    path = path.resolve(strict=True)
    if not path.is_file():
        raise ValueError("Native executable must be an installed regular file")
    process = subprocess.run(
        [str(path), *arguments],
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=30,
    )
    if process.returncode not in (0, 1) or marker not in process.stdout:
        raise ValueError("Native executable version check failed: " + path.name)
    return {"version": version, "path": str(path), "sha256": digest(path)}


def verify_context_native_receipt(tool, prefix):
    """Verify recorded native/model identities without executing or downloading."""
    native = tool.get("native_tools")
    expected = {
        "clinker-v2": set(),
        "DefenseFinder-v2": {"hmmsearch", "defense-finder-models"},
        "Phigaro-v2": {"hmmsearch", "prodigal", "pVOGs"},
        "ISEScan-v2": {
            "hmmsearch",
            "phmmer",
            "FragGeneScan",
            "blastn",
            "blastp",
            "makeblastdb",
        },
        "RNAalifold-v2": {"ViennaRNA", "mafft"},
    }[tool["name"]]
    if not isinstance(native, dict) or set(native) != expected:
        raise ValueError("Context native/model receipt missing: " + tool["name"])
    for name, item in native.items():
        path = Path(item["path"]).resolve(strict=True)
        path.relative_to(prefix.resolve(strict=True))
        if name == "defense-finder-models":
            if (
                item.get("commit") != DEFENSE_COMMIT
                or item.get("manifest_sha256") != DEFENSE_MANIFEST
            ):
                raise ValueError("DefenseFinder model receipt changed")
        elif digest(path) != item["sha256"]:
            raise ValueError("Context native/model bytes changed: " + name)
        if name in {"hmmsearch", "phmmer"} and item["version"] != "3.4":
            raise ValueError("HMMER version differs")
        if name == "FragGeneScan" and item["version"] != "1.32":
            raise ValueError("FragGeneScan version differs")
        if name == "prodigal" and item["version"] != "2.6.3":
            raise ValueError("Prodigal version differs")
        if name in {"blastn", "blastp", "makeblastdb"} and item["version"] != "2.16.0":
            raise ValueError("BLAST version differs")
        if name == "pVOGs" and (
            item["sha256"] != PVOG_HMM
            or item["profiles"] != PVOG_COUNT
            or item["archive_sha256"] != PVOG_ARCHIVE
        ):
            raise ValueError("pVOG reference receipt changed")
        if name == "mafft" and item["version"] != "7.526":
            raise ValueError("MAFFT version differs")
        if name == "ViennaRNA" and (
            item["version"] != "2.7.0" or item["binding_version"] != "2.7.0"
        ):
            raise ValueError("RNAalifold ViennaRNA identity differs")
