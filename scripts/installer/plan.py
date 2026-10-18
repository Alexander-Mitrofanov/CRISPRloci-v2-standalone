"""Validate explicit user choices and describe the actual installation stages."""
from __future__ import annotations

import os
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
GIB = 1024**3
MIN_FREE = 35 * GIB
REFERENCES = {
    "spacer": {"title": "Observed spacer search", "modes": ["sample", "download", "local"],
               "description": "Sample: 3 spacers. Download: the pinned 2017 spacerome, 720,391 source records (about 13 MiB compressed). Local: your spacer FASTA."},
    "viral": {"title": "Viral reference search", "modes": ["sample", "download", "local"],
              "description": "Sample: synthetic targets. Download: a frozen current NCBI viral RefSeq release, potentially many GiB. Local: your viral DNA FASTA."},
    "association": {"title": "Spacer-group associations", "modes": ["sample", "download", "local"],
                    "description": "Sample: bundled phages. Download: NCBI complete DNA phage/prokaryotic-virus records under the documented selection policy. Local: a directory of genome FASTA files."},
    "repeat": {"title": "Repeat reference mapping", "modes": ["sample", "historical", "local"],
               "description": "Sample: bundled repeat reference. Historical: checksum-pinned 2013 CRISPRmap reference. Local: a compatible reference.json."},
    "catalog": {"title": "Searchable reference catalogue", "modes": ["sample", "local", "import"],
                "description": "Sample: a clearly labelled synthetic catalogue. Local: a complete compatible catalogue snapshot directory. Import: a checksum-pinned CRISPRcasDb export bundle. The hosted production catalogue is not bundled."},
}


def defaults():
    return {"cpu": min(4, os.cpu_count() or 2), "developer": False,
            "references": {key: {"mode": "sample", "path": "", "sha256": ""} for key in REFERENCES}}


def choices(value):
    if not isinstance(value, dict) or set(value) - {"cpu", "developer", "references"}:
        raise ValueError("Unknown installer options")
    result = defaults()
    cpu = value.get("cpu", result["cpu"])
    if type(cpu) is not int or not 1 <= cpu <= 32:
        raise ValueError("Choose between 1 and 32 CPU threads")
    if type(value.get("developer", False)) is not bool:
        raise ValueError("Developer checks must be true or false")
    result.update(cpu=cpu, developer=value.get("developer", False))
    references = value.get("references", {})
    if not isinstance(references, dict) or set(references) - set(REFERENCES):
        raise ValueError("Unknown reference selection")
    for key, raw in references.items():
        if not isinstance(raw, dict) or set(raw) - {"mode", "path", "sha256"}:
            raise ValueError("Invalid reference options: " + key)
        item = {**result["references"][key], **raw}
        if item["mode"] not in REFERENCES[key]["modes"]:
            raise ValueError("Unsupported reference choice: " + key)
        for field in ("path", "sha256"):
            if not isinstance(item[field], str) or len(item[field]) > 4096 or "\0" in item[field]:
                raise ValueError("Invalid reference " + field)
        if item["mode"] in {"local", "import"}:
            path = Path(item["path"]).expanduser()
            if not item["path"] or not path.is_absolute() or not path.exists():
                raise ValueError(f"{REFERENCES[key]['title']}: select an existing absolute local path")
            if any(p.is_symlink() for p in (path, *path.parents)):
                raise ValueError("Reference paths must not traverse symlinks")
            directory = key in {"association", "catalog"} and item["mode"] == "local"
            if directory != path.is_dir():
                raise ValueError("Choose a reference directory" if directory else "Choose a reference file")
            item["path"] = str(path)
        if item["mode"] == "import" and not re.fullmatch(r"[0-9a-f]{64}", item["sha256"]):
            raise ValueError("Catalogue import requires the source bundle's SHA-256")
        result["references"][key] = item
    return result


def work_directory(root, value):
    lock = root / "tools/tools.lock.json"
    if lock.is_file() and lock.stat().st_size > 1024**2:
        raise ValueError("Scientific source lock exceeds its metadata limit")
    identity = {"choices": value, "sources": hashlib.sha256(lock.read_bytes() if lock.is_file() else b"").hexdigest()}
    return root / ".installer/work" / hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:20]


def system(root=ROOT, value=None):
    free = shutil.disk_usage(root).free
    minimum = MIN_FREE
    resuming = False
    if value is not None:
        checkpoint = work_directory(root, value) / "completed.json"
        if checkpoint.is_file() and not checkpoint.is_symlink() and checkpoint.stat().st_size < 16384:
            completed = json.loads(checkpoint.read_text())
            if isinstance(completed, list) and "bootstrap" in completed and (root / ".installer/bootstrap.json").is_file():
                minimum = 5 * GIB
                resuming = True
    errors = []
    if platform.system() != "Linux" or platform.machine() not in {"x86_64", "amd64"}:
        errors.append("Scientific installation requires Linux x86-64. Windows users can use an x86-64 Linux installation in WSL2; macOS is not supported by these pinned native environments.")
    if free < minimum:
        errors.append("At least 5 GiB free space is required to resume this installation; unfinished downloads may need more." if resuming else "At least 35 GiB free disk space is required for a fresh software installation. Large reference downloads need additional space.")
    if not (root / "tools/tools.lock.json").is_file():
        errors.append("The complete CRISPRloci source directory is missing.")
    return {"platform": platform.system(), "architecture": platform.machine(), "free_bytes": free,
            "minimum_bytes": minimum, "resuming": resuming, "cpu": os.cpu_count() or 1, "project": str(root), "errors": errors}


def stages(value):
    result = [
        {"id": "preflight", "title": "Check source, platform and disk space"},
        {"id": "bootstrap", "title": "Install the private Python, Node and build environment"},
        {"id": "application", "title": "Install the application and interface dependencies"},
        {"id": "tools", "title": "Build 17 isolated scientific environments and required models"},
        {"id": "release", "title": "Build and verify the immutable application release"},
    ]
    for key, item in value["references"].items():
        result.append({"id": "reference-" + key,
                       "title": REFERENCES[key]["title"] + ": " + item["mode"]})
    result += [{"id": "gate", "title": "Validate all analysis modes with real scientific jobs"}]
    if value["developer"]:
        result.append({"id": "tests", "title": "Run developer regression checks"})
    result += [{"id": "verify-local", "title": "Check selected databases and all modes through the local API"}]
    result += [{"id": "activate", "title": "Activate the validated local installation"}]
    return result
