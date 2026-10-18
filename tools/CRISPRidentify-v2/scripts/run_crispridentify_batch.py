#!/usr/bin/env python3
"""Bounded, restartable CRISPRidentify execution for large FASTA cohorts.

The detector remains a single-input program.  This supervisor gives each input
an isolated process group and immutable attempt directory, validates its
machine report against the source FASTA, and publishes ``completion.json`` only
after validation.  A disk-backed manifest index and a bounded Future window
keep memory use independent of cohort size.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import csv
import dataclasses
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import sqlite3
import stat
import subprocess
import sys
import tempfile
import threading
import time
from typing import Dict, Iterable, Iterator, Mapping, Optional, Sequence, Tuple


BATCH_SCHEMA = "CRISPRidentify-batch"
BATCH_VERSION = "1.0.0"
COMPLETION_SCHEMA = "CRISPRidentify-batch-completion"
COMPLETION_VERSION = 1
MANIFEST_FIELDS = ("sample_id", "fasta", "bytes", "sha256")
STATUS_FIELDS = (
    "ordinal",
    "sample_id",
    "status",
    "attempts",
    "elapsed_seconds",
    "error_code",
    "error_message",
    "completion",
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_COMMIT_RE = re.compile(r"^[0-9a-f]{40,64}$")
SAFE_SLUG_RE = re.compile(r"[^A-Za-z0-9_.-]+")
ATTEMPT_RE = re.compile(r"^attempt-([0-9a-f]{12})-([0-9]{4})$")
REPORT_SCHEMA = "CRISPRidentify-report"
CURRENT_REPORT_VERSION = "1.2.0"
ACCEPTED_CATEGORIES = frozenset(("Bona-fide", "Possible"))
CANONICAL_REPORT_CATEGORIES = (
    "Bona-fide",
    "Possible",
    "Low score",
)
LIKELY_TANDEM_SCHEMA = "CRISPRidentify-likely-tandem-quarantine"
LIKELY_TANDEM_SCHEMA_VERSION = "1.1.0"
LIKELY_TANDEM_RULE = "complete-link-low-spacer-diversity"
LIKELY_TANDEM_RULE_VERSION = "1.0.0"
LIKELY_TANDEM_BASIS_METHOD = "low_spacer_diversity_tandem_demotion"
LIKELY_TANDEM_TRANSACTION_RULE = "likely-tandem-category-quarantine"
LIKELY_TANDEM_TRANSACTION_RULE_VERSION = "1.1.0"
LIKELY_TANDEM_POLICY = {
    "minimum_spacer_count": 2,
    "maximum_spacer_count": 8,
    "minimum_pairwise_edit_identity": 0.95,
    "canonical_alphabet": "ACGT",
}
KNOWN_CATEGORIES = frozenset(
    (
        "Bona-fide",
        "Alternative",
        "Possible",
        "Possible Discarded",
        "Low score",
    )
)
PRESETS = ("region-screen", "balanced", "full")
FORBIDDEN_EXTRA_OPTIONS = frozenset(
    (
        "--file",
        "--input_folder",
        "--input_folder_multifasta",
        "--result_folder",
        "--cpu",
        "--folder_jobs",
        "--preset",
    )
)


class BatchError(RuntimeError):
    """A cohort-level contract was violated."""


class InputError(BatchError):
    """An input or manifest cannot be processed safely."""


class CompletionError(BatchError):
    """A completion marker or report failed validation."""


@dataclasses.dataclass(frozen=True)
class Job:
    ordinal: int
    sample_id: str
    fasta: Path
    source_bytes: int
    source_sha256: str


@dataclasses.dataclass(frozen=True)
class ToolConfiguration:
    runner: Path
    preset: str
    preset_mode: str
    analysis_arguments: Tuple[str, ...]
    extra_arguments: Tuple[str, ...]
    cpus_per_worker: int
    cas: bool
    cas_gene_mode: str
    cas_timeout_seconds: int
    topology_manifest: Optional[Path]
    runtime_tree_sha256: str
    runtime_version: str
    runtime_commit: Optional[str]


@dataclasses.dataclass(frozen=True)
class WorkerConfiguration:
    run_dir: Path
    scratch_root: Path
    tool: ToolConfiguration
    config_sha256: str
    timeout_seconds: float
    termination_grace_seconds: float
    max_attempts: int
    resume_verify_input: bool


@dataclasses.dataclass(frozen=True)
class JobResult:
    ordinal: int
    sample_id: str
    status: str
    attempts: int
    elapsed_seconds: float
    error_code: str = ""
    error_message: str = ""
    completion: str = ""


class ProcessRegistry:
    """Track active process groups so SIGINT/SIGTERM can stop all children."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._processes: Dict[int, subprocess.Popen] = {}
        self.stop_requested = threading.Event()

    def add(self, process: subprocess.Popen) -> None:
        with self._lock:
            self._processes[process.pid] = process

    def remove(self, process: subprocess.Popen) -> None:
        with self._lock:
            self._processes.pop(process.pid, None)

    def terminate_all(self, grace_seconds: float) -> None:
        self.stop_requested.set()
        with self._lock:
            processes = tuple(self._processes.values())
        for process in processes:
            terminate_process_group(process, grace_seconds)


ACTIVE_PROCESSES = ProcessRegistry()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def _canonical_sha256(value: object) -> str:
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


def _strict_json_equal(left: object, right: object) -> bool:
    """Compare JSON values without Python's bool/int numeric coercion."""

    try:
        return _canonical_json_bytes(left) == _canonical_json_bytes(right)
    except (TypeError, ValueError):
        return False


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(str(path), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write_bytes(path: Path, content: bytes) -> None:
    """Durably replace one file without exposing partially written content."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=".{}-".format(path.name)
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(path))
        _fsync_directory(path.parent)
    except BaseException:
        with contextlib.suppress(FileNotFoundError):
            temporary.unlink()
        raise


def atomic_write_json(path: Path, value: object) -> None:
    atomic_write_bytes(path, _canonical_json_bytes(value) + b"\n")


def _sample_slug(sample_id: str) -> str:
    readable = SAFE_SLUG_RE.sub("_", sample_id).strip("_.-")[:80] or "sample"
    suffix = hashlib.sha256(sample_id.encode("utf-8")).hexdigest()[:12]
    return "{}--{}".format(readable, suffix)


def shard_for_sample(sample_id: str, shard_count: int) -> int:
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    digest = hashlib.sha256(sample_id.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % shard_count


def _validate_manifest_text(value: str, label: str, line_number: int) -> str:
    if not isinstance(value, str) or not value:
        raise InputError("manifest line {} has empty {}".format(line_number, label))
    if value != value.strip() or any(ord(char) < 32 for char in value):
        raise InputError(
            "manifest line {} has unsafe {} whitespace/control characters".format(
                line_number, label
            )
        )
    return value


def index_manifest(
    manifest: Path,
    state_db: Path,
    shard_count: int,
    shard_index: int,
) -> Tuple[str, int, int]:
    """Validate and disk-index one selected shard without retaining all rows."""

    if not manifest.is_file():
        raise InputError("manifest is not a readable file: {}".format(manifest))
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise InputError("shard index must satisfy 0 <= index < count")

    manifest_sha256 = _sha256_file(manifest)
    temporary = state_db.with_name(".{}-{}".format(state_db.name, os.getpid()))
    with contextlib.suppress(FileNotFoundError):
        temporary.unlink()
    connection = sqlite3.connect(str(temporary))
    total = 0
    selected = 0
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute(
            """CREATE TABLE jobs (
                ordinal INTEGER PRIMARY KEY,
                sample_id TEXT NOT NULL UNIQUE,
                fasta TEXT NOT NULL,
                source_bytes INTEGER NOT NULL,
                source_sha256 TEXT NOT NULL,
                status TEXT NOT NULL DEFAULT 'pending',
                attempts INTEGER NOT NULL DEFAULT 0,
                elapsed_seconds REAL NOT NULL DEFAULT 0,
                error_code TEXT NOT NULL DEFAULT '',
                error_message TEXT NOT NULL DEFAULT '',
                completion TEXT NOT NULL DEFAULT ''
            )"""
        )
        with manifest.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != MANIFEST_FIELDS:
                raise InputError(
                    "manifest header must be exactly: {}".format(
                        "\t".join(MANIFEST_FIELDS)
                    )
                )
            for line_number, row in enumerate(reader, 2):
                total += 1
                sample_id = _validate_manifest_text(
                    row.get("sample_id", ""), "sample_id", line_number
                )
                fasta_text = _validate_manifest_text(
                    row.get("fasta", ""), "fasta", line_number
                )
                try:
                    source_bytes = int(row.get("bytes", ""))
                except (TypeError, ValueError):
                    raise InputError(
                        "manifest line {} bytes is not an integer".format(line_number)
                    )
                if source_bytes < 1:
                    raise InputError(
                        "manifest line {} bytes must be positive".format(line_number)
                    )
                source_sha256 = row.get("sha256", "")
                if SHA256_RE.fullmatch(source_sha256) is None:
                    raise InputError(
                        "manifest line {} sha256 must be lowercase SHA-256".format(
                            line_number
                        )
                    )
                if shard_for_sample(sample_id, shard_count) != shard_index:
                    continue
                fasta = Path(fasta_text).expanduser()
                if not fasta.is_absolute():
                    fasta = manifest.parent / fasta
                fasta = fasta.absolute()
                try:
                    connection.execute(
                        "INSERT INTO jobs(ordinal,sample_id,fasta,source_bytes,source_sha256) "
                        "VALUES(?,?,?,?,?)",
                        (total, sample_id, str(fasta), source_bytes, source_sha256),
                    )
                except sqlite3.IntegrityError as error:
                    raise InputError(
                        "duplicate selected-shard sample_id {!r}".format(sample_id)
                    ) from error
                selected += 1
                if selected % 10000 == 0:
                    connection.commit()
        if total == 0:
            raise InputError("manifest contains no data rows")
        connection.commit()
        connection.close()
        state_db.parent.mkdir(parents=True, exist_ok=True)
        os.replace(str(temporary), str(state_db))
        _fsync_directory(state_db.parent)
    except BaseException:
        with contextlib.suppress(Exception):
            connection.close()
        with contextlib.suppress(FileNotFoundError):
            temporary.unlink()
        raise
    return manifest_sha256, total, selected


def iter_indexed_jobs(state_db: Path, batch_size: int = 1024) -> Iterator[Job]:
    last_ordinal = 0
    connection = sqlite3.connect(str(state_db))
    try:
        while True:
            rows = connection.execute(
                "SELECT ordinal,sample_id,fasta,source_bytes,source_sha256 "
                "FROM jobs WHERE ordinal > ? ORDER BY ordinal LIMIT ?",
                (last_ordinal, batch_size),
            ).fetchall()
            if not rows:
                break
            for row in rows:
                last_ordinal = row[0]
                yield Job(row[0], row[1], Path(row[2]), row[3], row[4])
    finally:
        connection.close()


def _walk_regular_files(root: Path) -> Iterator[Path]:
    """Yield regular files using bounded directory-stack memory."""

    stack = [root]
    while stack:
        directory = stack.pop()
        try:
            entries = os.scandir(str(directory))
        except OSError as error:
            raise InputError("cannot scan {}: {}".format(directory, error)) from error
        with entries:
            for entry in entries:
                try:
                    if entry.is_dir(follow_symlinks=False):
                        stack.append(Path(entry.path))
                    elif entry.is_file(follow_symlinks=False):
                        yield Path(entry.path)
                except OSError as error:
                    raise InputError(
                        "cannot inspect {}: {}".format(entry.path, error)
                    ) from error


def create_manifest(
    input_root: Path,
    output: Path,
    suffixes: Sequence[str],
) -> int:
    """Create a sorted, checksummed manifest via a disk-backed index."""

    input_root = input_root.expanduser().absolute()
    output = output.expanduser().absolute()
    if not input_root.is_dir():
        raise InputError("input root is not a directory: {}".format(input_root))
    normalized_suffixes = tuple(
        suffix.lower() if suffix.startswith(".") else "." + suffix.lower()
        for suffix in suffixes
    )
    if not normalized_suffixes:
        raise InputError("at least one FASTA suffix is required")

    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, database_name = tempfile.mkstemp(
        dir=str(output.parent), prefix=".crispridentify-manifest-", suffix=".sqlite"
    )
    os.close(descriptor)
    database = Path(database_name)
    connection = sqlite3.connect(str(database))
    count = 0
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute(
            "CREATE TABLE entries(sample_id TEXT PRIMARY KEY, fasta TEXT UNIQUE, "
            "source_bytes INTEGER NOT NULL, source_sha256 TEXT NOT NULL)"
        )
        for path in _walk_regular_files(input_root):
            if path.suffix.lower() not in normalized_suffixes:
                continue
            relative = path.relative_to(input_root)
            sample_id = relative.with_suffix("").as_posix()
            source_bytes = path.stat().st_size
            if source_bytes < 1:
                raise InputError("FASTA is empty: {}".format(path))
            try:
                connection.execute(
                    "INSERT INTO entries VALUES(?,?,?,?)",
                    (sample_id, str(path), source_bytes, _sha256_file(path)),
                )
            except sqlite3.IntegrityError as error:
                raise InputError(
                    "input files produce duplicate sample_id {!r}".format(sample_id)
                ) from error
            count += 1
            if count % 1000 == 0:
                connection.commit()
        if count == 0:
            raise InputError("no FASTA files matched the requested suffixes")
        connection.commit()

        temporary = output.with_name(".{}-{}".format(output.name, os.getpid()))
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(MANIFEST_FIELDS)
            cursor = connection.execute(
                "SELECT sample_id,fasta,source_bytes,source_sha256 "
                "FROM entries ORDER BY sample_id"
            )
            for sample_id, fasta, source_bytes, source_sha256 in cursor:
                relative_fasta = os.path.relpath(fasta, str(output.parent))
                writer.writerow(
                    (sample_id, relative_fasta, source_bytes, source_sha256)
                )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(output))
        _fsync_directory(output.parent)
    finally:
        connection.close()
        with contextlib.suppress(FileNotFoundError):
            database.unlink()
    return count


def _runtime_version(project_root: Path) -> str:
    try:
        return (project_root / "VERSION").read_text(encoding="utf-8").strip()
    except OSError:
        return "unknown"


def _runtime_commit(project_root: Path) -> Optional[str]:
    try:
        completed = subprocess.run(
            ("git", "-C", str(project_root), "rev-parse", "HEAD"),
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    value = completed.stdout.strip()
    return (
        value if completed.returncode == 0 and GIT_COMMIT_RE.fullmatch(value) else None
    )


def _runtime_files(project_root: Path, runner: Path) -> Iterator[Tuple[str, Path]]:
    roots = (
        "VERSION",
        "CRISPRidentify.py",
        "crispridentify",
        "scripts/run_crispridentify.sh",
        "components",
        "trained_models",
        "tools/blasting",
        "vendor/casandra",
    )
    seen = set()
    for relative in roots:
        path = project_root / relative
        if path.is_file() and not path.is_symlink():
            seen.add(path.resolve())
            yield relative, path
        elif path.is_dir():
            for directory, names, files in os.walk(str(path)):
                names[:] = sorted(name for name in names if name != "__pycache__")
                for name in sorted(files):
                    child = Path(directory) / name
                    if (
                        child.is_symlink()
                        or not child.is_file()
                        or name.endswith(".pyc")
                    ):
                        continue
                    seen.add(child.resolve())
                    yield child.relative_to(project_root).as_posix(), child
    resolved_runner = runner.resolve()
    if resolved_runner not in seen:
        yield "external-runner", resolved_runner


def runtime_tree_sha256(project_root: Path, runner: Path) -> str:
    digest = hashlib.sha256()
    for relative, path in _runtime_files(project_root, runner):
        encoded = relative.encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "big"))
        digest.update(encoded)
        size = path.stat().st_size
        digest.update(size.to_bytes(8, "big"))
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def _probe_preset_support(runner: Path) -> bool:
    try:
        completed = subprocess.run(
            (str(runner), "--help"),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=60,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise BatchError("cannot inspect detector CLI: {}".format(error)) from error
    if completed.returncode != 0:
        raise BatchError(
            "detector --help failed with exit {}: {}".format(
                completed.returncode, completed.stdout[-1000:]
            )
        )
    return "--preset" in completed.stdout


def _fallback_preset_arguments(preset: str) -> Tuple[str, ...]:
    common = (
        "--fasta_report",
        "False",
        "--html_report",
        "False",
        "--is_element",
        "False",
        "--architecture_rescue",
        "False",
        "--fuzzy_search_mode",
        "native",
        "--fuzzy_search_flank",
        "100",
        "--fuzzy_family_padding",
        "15",
    )
    if preset == "region-screen":
        return common + (
            "--output_profile",
            "machine",
            "--accepted_only",
            "True",
            "--strand",
            "False",
            "--cas",
            "False",
            "--fast_run",
            "True",
            "--fast_run_seed_profile",
            "balanced",
            "--fast_run_rescue",
            "False",
            "--fast_run_phase_repair",
            "False",
            "--fast_two_repeat_terminal_repair",
            "False",
            "--fast_accepted_family_rescue",
            "False",
            "--split_competitive_components",
            "False",
            "--degenerated",
            "False",
            "--fuzzy_boundary_paths",
            "False",
        )
    if preset == "balanced":
        return common + (
            "--output_profile",
            "standard",
            "--accepted_only",
            "False",
            "--strand",
            "True",
            "--cas",
            "False",
            "--fast_run",
            "True",
            "--fast_run_seed_profile",
            "balanced",
            "--fast_run_rescue",
            "True",
            "--fast_run_phase_repair",
            "True",
            "--fast_two_repeat_terminal_repair",
            "False",
            "--fast_accepted_family_rescue",
            "False",
            "--split_competitive_components",
            "True",
            "--degenerated",
            "True",
            "--fuzzy_boundary_paths",
            "True",
        )
    if preset == "full":
        return common + (
            "--output_profile",
            "standard",
            "--accepted_only",
            "False",
            "--strand",
            "True",
            "--cas",
            "False",
            "--fast_run",
            "False",
            "--fast_run_rescue",
            "False",
            "--fast_run_phase_repair",
            "False",
            "--fast_two_repeat_terminal_repair",
            "False",
            "--fast_accepted_family_rescue",
            "False",
            "--split_competitive_components",
            "False",
            "--degenerated",
            "True",
            "--fuzzy_boundary_paths",
            "True",
        )
    raise ValueError("unknown preset: {}".format(preset))


def _validate_extra_arguments(arguments: Sequence[str]) -> Tuple[str, ...]:
    arguments = tuple(arguments)
    for token in arguments:
        if not isinstance(token, str) or "\x00" in token:
            raise InputError("extra tool arguments must be safe strings")
        option = token.split("=", 1)[0]
        if option in FORBIDDEN_EXTRA_OPTIONS:
            raise InputError("{} is controlled by crispridentify-batch".format(option))
    return arguments


def build_tool_configuration(
    args: argparse.Namespace, project_root: Path
) -> ToolConfiguration:
    runner = args.runner.expanduser().absolute()
    if not runner.is_file() or not os.access(str(runner), os.X_OK):
        raise InputError("detector launcher is not executable: {}".format(runner))
    has_preset = _probe_preset_support(runner)
    if has_preset:
        preset_mode = "detector-preset"
        analysis_arguments = ("--preset", args.preset)
    else:
        preset_mode = "compatibility-expansion"
        analysis_arguments = _fallback_preset_arguments(args.preset)
    if args.cas and args.preset == "region-screen":
        raise InputError(
            "--cas cannot be combined with region-screen; run CasAndra as a "
            "separately budgeted balanced/full batch"
        )
    if args.cas:
        analysis_arguments += (
            "--cas",
            "True",
            "--cas_gene_mode",
            args.cas_gene_mode,
            "--cas_timeout_seconds",
            str(args.cas_timeout_seconds),
        )
    topology_manifest = None
    if args.sequence_topology_manifest is not None:
        topology_manifest = args.sequence_topology_manifest.expanduser().absolute()
        if not topology_manifest.is_file():
            raise InputError(
                "sequence topology manifest is not a file: {}".format(topology_manifest)
            )
        analysis_arguments += ("--sequence_topology_manifest", str(topology_manifest))
    extra_arguments = _validate_extra_arguments(args.tool_arg)
    return ToolConfiguration(
        runner=runner,
        preset=args.preset,
        preset_mode=preset_mode,
        analysis_arguments=analysis_arguments,
        extra_arguments=extra_arguments,
        cpus_per_worker=args.cpus_per_worker,
        cas=args.cas,
        cas_gene_mode=args.cas_gene_mode,
        cas_timeout_seconds=args.cas_timeout_seconds,
        topology_manifest=topology_manifest,
        runtime_tree_sha256=runtime_tree_sha256(project_root, runner),
        runtime_version=_runtime_version(project_root),
        runtime_commit=_runtime_commit(project_root),
    )


def _tool_config_record(tool: ToolConfiguration) -> Mapping[str, object]:
    topology = None
    if tool.topology_manifest is not None:
        topology = {
            "path": str(tool.topology_manifest),
            "sha256": _sha256_file(tool.topology_manifest),
        }
    return {
        "preset": tool.preset,
        "preset_mode": tool.preset_mode,
        "analysis_arguments": list(tool.analysis_arguments),
        "extra_arguments": list(tool.extra_arguments),
        "cpus_per_worker": tool.cpus_per_worker,
        "cas": tool.cas,
        "cas_gene_mode": tool.cas_gene_mode if tool.cas else None,
        "cas_timeout_seconds": tool.cas_timeout_seconds if tool.cas else None,
        "sequence_topology_manifest": topology,
        "runtime": {
            "version": tool.runtime_version,
            "git_commit": tool.runtime_commit,
            "tree_sha256": tool.runtime_tree_sha256,
        },
    }


def _read_fasta(path: Path) -> Tuple[str, Dict[str, str]]:
    """Return the byte hash and unique source-ID to uppercase sequence map."""

    digest = hashlib.sha256()
    sources: Dict[str, list] = {}
    current_id = None
    try:
        handle = path.open("rb")
    except OSError as error:
        raise InputError("cannot read FASTA {}: {}".format(path, error)) from error
    with handle:
        for line_number, raw_line in enumerate(handle, 1):
            digest.update(raw_line)
            try:
                line = raw_line.decode("ascii").strip()
            except UnicodeDecodeError as error:
                raise InputError(
                    "FASTA {} line {} is not ASCII".format(path, line_number)
                ) from error
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                source_id = header.split()[0] if header else ""
                if not source_id:
                    raise InputError(
                        "FASTA {} line {} has an empty header".format(path, line_number)
                    )
                if source_id in sources:
                    raise InputError(
                        "FASTA {} repeats source ID {!r}".format(path, source_id)
                    )
                sources[source_id] = []
                current_id = source_id
            else:
                if current_id is None:
                    raise InputError(
                        "FASTA {} has sequence before its first header".format(path)
                    )
                sequence = "".join(line.split()).upper()
                sources[current_id].append(sequence)
    if not sources:
        raise InputError("FASTA has no records: {}".format(path))
    joined = {source_id: "".join(parts) for source_id, parts in sources.items()}
    for source_id, sequence in joined.items():
        if not sequence:
            raise InputError(
                "FASTA source {!r} has an empty sequence".format(source_id)
            )
    return digest.hexdigest(), joined


def _safe_relative_file(root: Path, relative_text: str) -> Path:
    relative = PurePosixPath(relative_text)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise CompletionError("completion contains an unsafe report path")
    target = root.joinpath(*relative.parts)
    try:
        metadata = target.lstat()
    except OSError as error:
        raise CompletionError(
            "completion report is missing: {}".format(target)
        ) from error
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise CompletionError("completion report is not a regular file")
    return target


def _is_plain_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _validated_interval(
    value: object,
    label: str,
    source_length: int,
    containing: Optional[Mapping[str, int]] = None,
) -> Mapping[str, int]:
    if not isinstance(value, dict):
        raise CompletionError("{} must be an object".format(label))
    start = value.get("start")
    end = value.get("end")
    length = value.get("length")
    if (
        not _is_plain_int(start)
        or not _is_plain_int(end)
        or not _is_plain_int(length)
        or start < 1
        or end < start
        or length != end - start + 1
        or end > source_length
    ):
        raise CompletionError(
            "{} is not a valid source interval".format(label)
        )
    interval = {"start": start, "end": end, "length": length}
    if containing is not None and (
        start < containing["start"] or end > containing["end"]
    ):
        raise CompletionError("{} lies outside its array".format(label))
    return interval


def _validated_tags(value: object, label: str) -> Tuple[str, ...]:
    if value is None:
        value = []
    if (
        not isinstance(value, list)
        or any(
            not isinstance(tag, str)
            or re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", tag) is None
            for tag in value
        )
        or value != sorted(set(value))
    ):
        raise CompletionError(
            "{} tags must be sorted unique identifiers".format(label)
        )
    return tuple(value)


def _global_edit_distance(first: str, second: str) -> int:
    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for first_index, first_base in enumerate(first, 1):
        current = [first_index]
        for second_index, second_base in enumerate(second, 1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[second_index] + 1,
                    previous[second_index - 1]
                    + (first_base != second_base),
                )
            )
        previous = current
    return previous[-1]


def _likely_tandem_evidence(
    spacers: Sequence[str],
) -> Tuple[bool, Optional[Mapping[str, object]]]:
    """Replay the fixed complete-link rule without importing detector code."""

    normalized = tuple(spacer.upper() for spacer in spacers)
    if not 2 <= len(normalized) <= 8:
        return False, None
    if any(
        not spacer or not set(spacer).issubset({"A", "C", "G", "T"})
        for spacer in normalized
    ):
        return False, None

    identities = []
    distances = []
    minimum_pair = None
    for left_index, left in enumerate(normalized[:-1]):
        for right_index in range(left_index + 1, len(normalized)):
            right = normalized[right_index]
            distance = _global_edit_distance(left, right)
            maximum_length = max(len(left), len(right))
            numerator = maximum_length - distance
            identity = numerator / float(maximum_length)
            identities.append(identity)
            distances.append(distance)
            if (
                minimum_pair is None
                or numerator * minimum_pair[1]
                < minimum_pair[0] * maximum_length
                or (
                    numerator * minimum_pair[1]
                    == minimum_pair[0] * maximum_length
                    and (left_index, right_index) < minimum_pair[2]
                )
            ):
                minimum_pair = (
                    numerator,
                    maximum_length,
                    (left_index, right_index),
                    distance,
                )

    threshold_satisfied = (
        (minimum_pair[1] - minimum_pair[3]) * 20
        >= minimum_pair[1] * 19
    )
    if not threshold_satisfied:
        return False, None
    lengths = tuple(len(spacer) for spacer in normalized)
    distinct_count = len(set(normalized))
    evidence = {
        "decision": "likely_tandem",
        "category": "likely_tandem",
        "reason_code": "all_spacer_pairs_at_or_above_identity_threshold",
        "rule_name": LIKELY_TANDEM_RULE,
        "rule_version": LIKELY_TANDEM_RULE_VERSION,
        "coordinate_system": "source_forward_strand",
        "policy": dict(LIKELY_TANDEM_POLICY),
        "spacer_count": len(normalized),
        "nonempty_spacer_count": len(normalized),
        "canonical_spacer_count": len(normalized),
        "distinct_spacer_count": distinct_count,
        "exact_duplicate_excess_count": len(normalized) - distinct_count,
        "shortest_spacer_length": min(lengths),
        "longest_spacer_length": max(lengths),
        "pairwise_comparison_count": len(identities),
        "minimum_pairwise_edit_identity": (
            minimum_pair[0] / float(minimum_pair[1])
        ),
        "maximum_pairwise_edit_identity": max(identities),
        "maximum_pairwise_edit_distance": max(distances),
        "minimum_identity_pair_1based": [
            minimum_pair[2][0] + 1,
            minimum_pair[2][1] + 1,
        ],
        "minimum_identity_pair_edit_distance": minimum_pair[3],
        "minimum_identity_pair_maximum_length": minimum_pair[1],
        "complete_link_threshold_satisfied": True,
    }
    return True, evidence


def _validate_array(
    array: object,
    sequence: Optional[str],
    source_length: int,
    region_screen: bool,
    strict: bool,
) -> Optional[Mapping[str, object]]:
    if not isinstance(array, dict):
        raise CompletionError("report array must be an object")
    array_id = array.get("id")
    category = array.get("category")
    if category not in KNOWN_CATEGORIES:
        raise CompletionError("report array has unknown category")
    if region_screen and category not in ACCEPTED_CATEGORIES:
        raise CompletionError("region-screen report contains a non-accepted array")
    if not strict:
        return None
    if category not in CANONICAL_REPORT_CATEGORIES:
        raise CompletionError("current report contains an internal category")

    tags = _validated_tags(array.get("tags"), "array {!r}".format(array_id))
    if "tandem" in tags and category != "Low score":
        raise CompletionError("tandem tag is restricted to Low score arrays")
    validation = array.get("validation")
    if (
        not isinstance(validation, dict)
        or validation.get("status") != "valid"
        or validation.get("source_reconstructed") is not True
        or validation.get("method")
        != "independent_original_sequence_slicing"
    ):
        raise CompletionError(
            "array {!r} lacks canonical source validation".format(array_id)
        )

    outer = _validated_interval(
        array.get("source_interval"),
        "array {!r} source_interval".format(array_id),
        source_length,
    )
    source_sequence = array.get("source_sequence")
    if (
        not isinstance(source_sequence, str)
        or not source_sequence.isascii()
        or source_sequence != source_sequence.upper()
        or len(source_sequence) != outer["length"]
    ):
        raise CompletionError(
            "array {!r} source_sequence length differs".format(array_id)
        )
    if sequence is not None and source_sequence != sequence[
        outer["start"] - 1 : outer["end"]
    ]:
        raise CompletionError("report array sequence does not match its source")

    def expected_slice(interval: Mapping[str, int]) -> str:
        offset_start = interval["start"] - outer["start"]
        offset_end = interval["end"] - outer["start"] + 1
        return source_sequence[offset_start:offset_end]

    repeat_count = array.get("repeat_count")
    repeats = array.get("repeats")
    if (
        not _is_plain_int(repeat_count)
        or repeat_count < 1
        or not isinstance(repeats, list)
        or len(repeats) != repeat_count
    ):
        raise CompletionError(
            "array {!r} repeat table/count differs".format(array_id)
        )
    validated_repeats = []
    previous_end = None
    for position, repeat in enumerate(repeats, 1):
        if (
            not isinstance(repeat, dict)
            or not _is_plain_int(repeat.get("index"))
            or repeat.get("index") != position
        ):
            raise CompletionError(
                "array {!r} repeat indexes differ".format(array_id)
            )
        interval = _validated_interval(
            repeat.get("source_interval"),
            "array {!r} repeat {}".format(array_id, position),
            source_length,
            containing=outer,
        )
        repeat_sequence = repeat.get("sequence")
        if (
            not isinstance(repeat_sequence, str)
            or not repeat_sequence
            or not repeat_sequence.isascii()
            or repeat_sequence != repeat_sequence.upper()
            or len(repeat_sequence) != interval["length"]
            or repeat_sequence != expected_slice(interval)
        ):
            raise CompletionError(
                "array {!r} repeat {} differs from source".format(
                    array_id, position
                )
            )
        if previous_end is not None and interval["start"] <= previous_end:
            raise CompletionError(
                "array {!r} repeats overlap or are unordered".format(array_id)
            )
        validated_repeats.append(
            {
                "index": position,
                "source_interval": interval,
                "sequence": repeat_sequence,
            }
        )
        previous_end = interval["end"]
    if (
        validated_repeats[0]["source_interval"]["start"] != outer["start"]
        or validated_repeats[-1]["source_interval"]["end"] != outer["end"]
    ):
        raise CompletionError(
            "array {!r} interval does not span its repeats".format(array_id)
        )

    spacers = array.get("spacers")
    deletions = array.get("spacer_deletions")
    architecture = array.get("architecture")
    spacer_count = array.get("spacer_count")
    non_deleted_count = array.get("non_deleted_spacer_count")
    deletion_count = array.get("spacer_deletion_count")
    if (
        not isinstance(spacers, list)
        or not isinstance(deletions, list)
        or not isinstance(architecture, list)
        or not _is_plain_int(spacer_count)
        or not _is_plain_int(non_deleted_count)
        or not _is_plain_int(deletion_count)
        or spacer_count != repeat_count - 1
        or non_deleted_count != len(spacers)
        or deletion_count != len(deletions)
        or non_deleted_count + deletion_count != spacer_count
    ):
        raise CompletionError(
            "array {!r} spacer table/count differs".format(array_id)
        )

    expected_spacers = []
    expected_deletions = []
    expected_architecture = []
    rule_spacers = []
    for position, repeat in enumerate(validated_repeats, 1):
        expected_architecture.append(dict(repeat, type="repeat"))
        if position == repeat_count:
            continue
        right = validated_repeats[position]["source_interval"]
        left = repeat["source_interval"]
        gap_length = right["start"] - left["end"] - 1
        if gap_length == 0:
            boundary = {
                "position": left["end"],
                "coordinate_system": "0-based interbase",
                "reference": "input_forward_strand",
            }
            deletion = {
                "index": position,
                "source_boundary": boundary,
                "flanking_repeat_indices": [position, position + 1],
                "length": 0,
                "sequence": "",
            }
            expected_deletions.append(deletion)
            expected_architecture.append(dict(deletion, type="spacer-deletion"))
            rule_spacers.append("")
        else:
            interval = {
                "start": left["end"] + 1,
                "end": right["start"] - 1,
                "length": gap_length,
            }
            spacer = {
                "index": position,
                "source_interval": interval,
                "sequence": expected_slice(interval),
            }
            expected_spacers.append(spacer)
            expected_architecture.append(dict(spacer, type="spacer"))
            rule_spacers.append(spacer["sequence"])
    if (
        not _strict_json_equal(spacers, expected_spacers)
        or not _strict_json_equal(deletions, expected_deletions)
        or not _strict_json_equal(architecture, expected_architecture)
        or source_sequence
        != "".join(element.get("sequence", "") for element in architecture)
    ):
        raise CompletionError(
            "array {!r} architecture/table differs from source".format(array_id)
        )

    if category in ACCEPTED_CATEGORIES:
        likely_tandem, _evidence = _likely_tandem_evidence(rule_spacers)
        if likely_tandem:
            raise CompletionError(
                "accepted report array {!r} still qualifies as likely "
                "tandem".format(array_id)
            )
    return {
        "source_interval": outer,
        "source_sequence": source_sequence,
        "repeats": validated_repeats,
        "spacers": expected_spacers,
        "rule_spacers": tuple(rule_spacers),
        "tags": tags,
    }


def _validate_ledger_architecture(
    record: Mapping[str, object],
    outer: Mapping[str, int],
    sequence: Optional[str],
) -> Mapping[str, object]:
    architecture = record.get("source_architecture")
    if (
        not isinstance(architecture, dict)
        or architecture.get("coordinate_system")
        != "zero_based_half_open"
        or architecture.get("reference") != "input_forward_strand"
    ):
        raise CompletionError("likely-tandem ledger architecture is invalid")
    repeats = architecture.get("repeats")
    spacers = architecture.get("spacers")
    spacer_copies = record.get("spacers_source_forward")
    if (
        not isinstance(repeats, list)
        or not isinstance(spacers, list)
        or not isinstance(spacer_copies, list)
        or len(repeats) != len(spacers) + 1
        or len(spacers) != len(spacer_copies)
    ):
        raise CompletionError(
            "likely-tandem ledger architecture cardinality differs"
        )

    outer_start = outer["start"] - 1
    outer_end = outer["end"]

    def validate_elements(
        elements: Sequence[object], label: str
    ) -> Tuple[Mapping[str, object], ...]:
        validated = []
        for position, element in enumerate(elements, 1):
            if (
                not isinstance(element, dict)
                or not _is_plain_int(element.get("index"))
                or element.get("index") != position
            ):
                raise CompletionError(
                    "likely-tandem ledger {} indexes differ".format(label)
                )
            interval = element.get("source_interval")
            if not isinstance(interval, dict) or set(interval) != {"start", "end"}:
                raise CompletionError(
                    "likely-tandem ledger {} interval is invalid".format(label)
                )
            start = interval.get("start")
            end = interval.get("end")
            element_sequence = element.get("sequence")
            if (
                not _is_plain_int(start)
                or not _is_plain_int(end)
                or start < outer_start
                or end <= start
                or end > outer_end
                or not isinstance(element_sequence, str)
                or not element_sequence
                or not element_sequence.isascii()
                or element_sequence != element_sequence.upper()
                or len(element_sequence) != end - start
                or element.get("sequence_sha256")
                != hashlib.sha256(element_sequence.encode("ascii")).hexdigest()
            ):
                raise CompletionError(
                    "likely-tandem ledger {} source element is invalid".format(
                        label
                    )
                )
            if sequence is not None and element_sequence != sequence[start:end]:
                raise CompletionError(
                    "likely-tandem ledger {} differs from source".format(label)
                )
            validated.append(
                {
                    "index": position,
                    "start": start,
                    "end": end,
                    "sequence": element_sequence,
                }
            )
        return tuple(validated)

    validated_repeats = validate_elements(repeats, "repeat")
    validated_spacers = validate_elements(spacers, "spacer")
    if tuple(item["sequence"] for item in validated_spacers) != tuple(
        spacer_copies
    ):
        raise CompletionError("likely-tandem ledger spacer copies differ")
    if (
        validated_repeats[0]["start"] != outer_start
        or validated_repeats[-1]["end"] != outer_end
    ):
        raise CompletionError(
            "likely-tandem ledger architecture does not span its interval"
        )

    reconstructed = []
    for position, repeat in enumerate(validated_repeats):
        reconstructed.append(repeat["sequence"])
        if position < len(validated_spacers):
            spacer = validated_spacers[position]
            if (
                repeat["end"] != spacer["start"]
                or spacer["end"] != validated_repeats[position + 1]["start"]
            ):
                raise CompletionError(
                    "likely-tandem ledger architecture is not contiguous"
                )
            reconstructed.append(spacer["sequence"])
    reconstructed_sequence = "".join(reconstructed)
    if (
        len(reconstructed_sequence) != outer["length"]
        or architecture.get("reconstructed_source_span_sha256")
        != hashlib.sha256(
            reconstructed_sequence.encode("ascii")
        ).hexdigest()
        or (
            sequence is not None
            and reconstructed_sequence != sequence[outer_start:outer_end]
        )
    ):
        raise CompletionError(
            "likely-tandem ledger reconstructed source span differs"
        )
    return {
        "repeats": validated_repeats,
        "spacers": validated_spacers,
        "reconstructed_sequence": reconstructed_sequence,
    }


def _core_array_matches_ledger(
    array: Mapping[str, object],
    geometry: Mapping[str, object],
    record: Mapping[str, object],
    ledger_geometry: Mapping[str, object],
) -> bool:
    basis = array.get("classification_basis")
    transition = basis.get("transition") if isinstance(basis, dict) else None
    if (
        array.get("category") != "Low score"
        or "tandem" not in geometry["tags"]
        or not _strict_json_equal(
            array.get("source_interval"), record.get("source_interval")
        )
        or not _strict_json_equal(
            array.get("legacy_interval_key"), record.get("locus_key")
        )
        or isinstance(array.get("certainty_score"), bool)
        or not isinstance(array.get("certainty_score"), (int, float))
        or not math.isfinite(float(array["certainty_score"]))
        or not _strict_json_equal(
            array.get("certainty_score"), record.get("score")
        )
        or not isinstance(basis, dict)
        or basis.get("method") != LIKELY_TANDEM_BASIS_METHOD
        or basis.get("disposition") != "likely_tandem"
        or basis.get("tag") != "tandem"
        or not _strict_json_equal(
            basis.get("likely_tandem"), record.get("evidence")
        )
        or not isinstance(transition, dict)
        or transition.get("rule_name") != LIKELY_TANDEM_TRANSACTION_RULE
        or transition.get("rule_version")
        != LIKELY_TANDEM_TRANSACTION_RULE_VERSION
        or not _strict_json_equal(
            transition.get("source_category_index"),
            record.get("source_category_index"),
        )
        or transition.get("source_category_name")
        != record.get("source_category_name")
        or transition.get("target_category_index") != 4
        or transition.get("target_category_name") != "Low score"
        or not _strict_json_equal(
            transition.get("row_ordinal"), record.get("row_ordinal")
        )
        or not _strict_json_equal(
            transition.get("locus_key"), record.get("locus_key")
        )
        or transition.get("reported_disposition") != "likely_tandem"
        or transition.get("reported_tag") != "tandem"
    ):
        return False

    core_repeats = geometry["repeats"]
    core_spacers = geometry["spacers"]
    ledger_repeats = ledger_geometry["repeats"]
    ledger_spacers = ledger_geometry["spacers"]
    if (
        len(core_repeats) != len(ledger_repeats)
        or len(core_spacers) != len(ledger_spacers)
    ):
        return False
    for core, ledger in zip(core_repeats, ledger_repeats):
        interval = core["source_interval"]
        if (
            core["sequence"] != ledger["sequence"]
            or interval["start"] - 1 != ledger["start"]
            or interval["end"] != ledger["end"]
        ):
            return False
    for core, ledger in zip(core_spacers, ledger_spacers):
        interval = core["source_interval"]
        if (
            core["sequence"] != ledger["sequence"]
            or interval["start"] - 1 != ledger["start"]
            or interval["end"] != ledger["end"]
        ):
            return False
    return geometry["source_sequence"] == ledger_geometry[
        "reconstructed_sequence"
    ]


def _validate_likely_tandem_extension(
    report: Mapping[str, object],
    arrays: Sequence[Mapping[str, object]],
    geometries: Sequence[Mapping[str, object]],
    sequence: Optional[str],
    source_sha256: str,
    source_length: int,
    accepted_only: bool,
) -> None:
    extensions = report.get("extensions")
    if not isinstance(extensions, dict):
        raise CompletionError(
            "current report lacks the likely-tandem audit ledger"
        )
    extension = extensions.get("likely_tandem")
    if not isinstance(extension, dict):
        raise CompletionError(
            "current report lacks the likely-tandem audit ledger"
        )
    if (
        extension.get("schema")
        != {
            "name": LIKELY_TANDEM_SCHEMA,
            "version": LIKELY_TANDEM_SCHEMA_VERSION,
        }
        or extension.get("source_sha256") != source_sha256
        or extension.get("disposition") != "likely_tandem"
        or extension.get("tag") != "tandem"
        or extension.get("reported_category") != "Low score"
        or extension.get("status") != "applied"
        or extension.get("external_tool_support_used_for_decision") is not False
    ):
        raise CompletionError("likely-tandem audit ledger envelope differs")
    rule = extension.get("rule")
    if (
        not isinstance(rule, dict)
        or rule.get("name") != LIKELY_TANDEM_RULE
        or rule.get("version") != LIKELY_TANDEM_RULE_VERSION
        or rule.get("linkage") != "complete_link_all_spacer_pairs"
        or rule.get("identity")
        != "1 - global_levenshtein/max(sequence_lengths)"
        or not _strict_json_equal(rule.get("policy"), LIKELY_TANDEM_POLICY)
    ):
        raise CompletionError("likely-tandem audit rule differs")

    records = extension.get("arrays")
    array_count = extension.get("array_count")
    stats = extension.get("transaction_stats")
    required_stats = {
        "total_rows_before",
        "total_rows_after",
        "accepted_rows_reviewed",
        "accepted_rows_before",
        "accepted_rows_after",
        "likely_tandem_count",
        "retained_count",
        "bona_fide_quarantined",
        "possible_quarantined",
        "diagnostic_rows_before",
        "diagnostic_rows_after",
    }
    if (
        not isinstance(records, list)
        or not _is_plain_int(array_count)
        or array_count != len(records)
        or not isinstance(stats, dict)
        or set(stats) != required_stats
        or any(not _is_plain_int(value) or value < 0 for value in stats.values())
        or stats["likely_tandem_count"] != len(records)
    ):
        raise CompletionError("likely-tandem audit counts are invalid")

    core_by_id = {
        array["id"]: (array, geometry)
        for array, geometry in zip(arrays, geometries)
    }
    mapped_core_ids = set()
    observed_record_keys = set()
    observed_record_ids = set()
    source_category_counts = {"Bona-fide": 0, "Possible": 0}
    for position, record in enumerate(records, 1):
        if not isinstance(record, dict):
            raise CompletionError("likely-tandem ledger row must be an object")
        source_category_index = record.get("source_category_index")
        expected_source_name = {0: "Bona-fide", 2: "Possible"}.get(
            source_category_index
        )
        score = record.get("score")
        if (
            record.get("id") != "LT{}".format(position)
            or record.get("action") != "quarantined"
            or record.get("disposition") != "likely_tandem"
            or record.get("tag") != "tandem"
            or not _is_plain_int(source_category_index)
            or expected_source_name is None
            or record.get("source_category_name") != expected_source_name
            or record.get("target_category_index") != 4
            or record.get("target_category_name") != "Low score"
            or not _is_plain_int(record.get("row_ordinal"))
            or record["row_ordinal"] < 0
            or isinstance(score, bool)
            or not isinstance(score, (int, float))
            or not math.isfinite(float(score))
        ):
            raise CompletionError("likely-tandem ledger row is invalid")
        if record["id"] in observed_record_ids:
            raise CompletionError("likely-tandem ledger IDs are not unique")
        observed_record_ids.add(record["id"])
        outer = _validated_interval(
            record.get("source_interval"),
            "likely-tandem ledger source_interval",
            source_length,
        )
        ledger_geometry = _validate_ledger_architecture(
            record, outer, sequence
        )
        spacers = tuple(
            item["sequence"] for item in ledger_geometry["spacers"]
        )
        qualifies, evidence = _likely_tandem_evidence(spacers)
        if (
            not qualifies
            or not _strict_json_equal(record.get("evidence"), evidence)
            or record.get("repeat_count") != len(spacers) + 1
        ):
            raise CompletionError(
                "likely-tandem ledger evidence does not replay"
            )
        try:
            locus_key_json = json.dumps(
                record.get("locus_key"),
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
        except (TypeError, ValueError) as error:
            raise CompletionError(
                "likely-tandem ledger locus key is invalid"
            ) from error
        record_key = (
            outer["start"],
            outer["end"],
            source_category_index,
            locus_key_json,
            record["row_ordinal"],
        )
        if record_key in observed_record_keys:
            raise CompletionError("likely-tandem ledger rows are not unique")
        observed_record_keys.add(record_key)
        source_category_counts[expected_source_name] += 1

        if accepted_only:
            if (
                record.get("reported_array_id") is not None
                or record.get("core_projection") != "omitted_by_accepted_only"
            ):
                raise CompletionError(
                    "accepted-only likely-tandem ledger projection differs"
                )
        else:
            reported_id = record.get("reported_array_id")
            if (
                not isinstance(reported_id, str)
                or not reported_id
                or reported_id in mapped_core_ids
                or record.get("core_projection") != "reported_as_low_score"
                or reported_id not in core_by_id
            ):
                raise CompletionError(
                    "likely-tandem ledger/core projection is invalid"
                )
            core_array, core_geometry = core_by_id[reported_id]
            if not _core_array_matches_ledger(
                core_array, core_geometry, record, ledger_geometry
            ):
                raise CompletionError(
                    "likely-tandem ledger does not match its core array"
                )
            mapped_core_ids.add(reported_id)

    tandem_core_ids = {
        array["id"]
        for array, geometry in zip(arrays, geometries)
        if "tandem" in geometry["tags"]
    }
    if accepted_only:
        if tandem_core_ids:
            raise CompletionError(
                "accepted-only report contains a tandem-tagged core array"
            )
    elif tandem_core_ids != mapped_core_ids:
        raise CompletionError(
            "tandem-tagged core arrays and audit ledger differ"
        )

    summary = report["summary"]
    category_counts = summary["category_counts"]
    if (
        stats["accepted_rows_reviewed"] != stats["accepted_rows_before"]
        or stats["accepted_rows_after"]
        != stats["accepted_rows_before"] - len(records)
        or stats["retained_count"] != stats["accepted_rows_after"]
        or stats["bona_fide_quarantined"]
        != source_category_counts["Bona-fide"]
        or stats["possible_quarantined"]
        != source_category_counts["Possible"]
        or stats["bona_fide_quarantined"]
        + stats["possible_quarantined"]
        != len(records)
        or stats["total_rows_before"] != stats["total_rows_after"]
        or stats["diagnostic_rows_after"]
        != stats["diagnostic_rows_before"] + len(records)
        or stats["accepted_rows_after"] != summary["accepted_array_count"]
        or (
            not accepted_only
            and stats["diagnostic_rows_after"]
            != category_counts["Low score"]
        )
        or (accepted_only and category_counts["Low score"] != 0)
    ):
        raise CompletionError(
            "likely-tandem audit transaction stats differ from report"
        )


def validate_report(
    path: Path,
    expected_sources: Optional[Mapping[str, str]],
    expected_record: Optional[Mapping[str, object]],
    region_screen: bool,
    expected_preset: Optional[str] = None,
    require_preset_provenance: bool = False,
) -> Mapping[str, object]:
    try:
        report = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CompletionError(
            "cannot parse report {}: {}".format(path, error)
        ) from error
    if not isinstance(report, dict):
        raise CompletionError("report root must be an object")
    strict_completion = expected_sources is not None or expected_record is not None
    run = report.get("run")
    configuration = run.get("configuration") if isinstance(run, dict) else None
    if require_preset_provenance and not isinstance(configuration, dict):
        raise CompletionError("report lacks detector preset provenance")
    if isinstance(configuration, dict) and expected_preset is not None:
        if configuration.get("preset") != expected_preset:
            raise CompletionError("report detector preset differs")
    schema = report.get("schema")
    if not isinstance(schema, dict) or schema.get("name") != REPORT_SCHEMA:
        raise CompletionError("report schema name differs")
    if not isinstance(schema.get("version"), str) or not schema["version"]:
        raise CompletionError("report schema version is missing")
    if strict_completion and schema["version"] != CURRENT_REPORT_VERSION:
        raise CompletionError("completion requires the current report schema")
    source = report.get("source")
    if not isinstance(source, dict):
        raise CompletionError("report source is missing")
    source_id = source.get("id")
    source_sha256 = source.get("sha256")
    source_length = source.get("length")
    if (
        not isinstance(source_id, str)
        or not source_id
        or SHA256_RE.fullmatch(str(source_sha256)) is None
        or isinstance(source_length, bool)
        or not isinstance(source_length, int)
        or source_length < 1
    ):
        raise CompletionError("report source identity is invalid")

    sequence = None
    if expected_sources is not None:
        sequence = expected_sources.get(source_id)
        if sequence is None:
            raise CompletionError("report source ID is absent from input FASTA")
        observed_source_hash = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if source_sha256 != observed_source_hash or source_length != len(sequence):
            raise CompletionError("report source hash/length differs from input FASTA")
    if expected_record is not None:
        if (
            source_id != expected_record.get("source_id")
            or source_sha256 != expected_record.get("source_sha256")
            or source_length != expected_record.get("source_length")
        ):
            raise CompletionError("report source differs from completion record")

    arrays = report.get("arrays")
    summary = report.get("summary")
    if not isinstance(arrays, list) or not isinstance(summary, dict):
        raise CompletionError("report arrays/summary contract is missing")
    if summary.get("array_count") != len(arrays):
        raise CompletionError("report summary array count differs")
    identifiers = []
    geometries = []
    for array in arrays:
        geometry = _validate_array(
            array,
            sequence,
            source_length,
            region_screen,
            strict_completion,
        )
        geometries.append(geometry)
        identifiers.append(array.get("id"))
    if any(not isinstance(value, str) or not value for value in identifiers):
        raise CompletionError("report array ID is missing")
    if len(identifiers) != len(set(identifiers)):
        raise CompletionError("report array IDs are not unique")
    if strict_completion:
        parameters = run.get("parameters") if isinstance(run, dict) else None
        accepted_only = (
            parameters.get("param_accepted_only")
            if isinstance(parameters, dict)
            else None
        )
        if not isinstance(accepted_only, bool):
            raise CompletionError(
                "current report lacks accepted-only projection provenance"
            )
        if region_screen and accepted_only is not True:
            raise CompletionError(
                "region-screen report is not an accepted-only projection"
            )
        observed_counts = {
            category: sum(
                array.get("category") == category for array in arrays
            )
            for category in CANONICAL_REPORT_CATEGORIES
        }
        report_validation = report.get("validation")
        if (
            summary.get("accepted_array_count")
            != observed_counts["Bona-fide"] + observed_counts["Possible"]
            or summary.get("category_counts") != observed_counts
        ):
            raise CompletionError("report summary category counts differ")
        if (
            not isinstance(report_validation, dict)
            or report_validation.get("status") != "valid"
            or report_validation.get("method")
            != "independent_original_sequence_slicing"
            or report_validation.get("validated_array_count") != len(arrays)
            or report_validation.get("invalid_array_count") != 0
        ):
            raise CompletionError("report validation summary differs")
        _validate_likely_tandem_extension(
            report,
            arrays,
            geometries,
            sequence,
            source_sha256,
            source_length,
            accepted_only,
        )
    return {
        "source_id": source_id,
        "source_sha256": source_sha256,
        "source_length": source_length,
        "schema_version": schema["version"],
        "array_count": len(arrays),
        "report_sha256": _sha256_file(path),
    }


def collect_and_validate_reports(
    result_root: Path,
    sources: Mapping[str, str],
    preset: str,
    require_preset_provenance: bool,
) -> Tuple[Mapping[str, object], ...]:
    report_paths = sorted(result_root.rglob("report.json"))
    if not report_paths:
        raise CompletionError("detector produced no report.json")
    records = []
    seen_sources = set()
    for report_path in report_paths:
        if report_path.is_symlink() or not report_path.is_file():
            raise CompletionError("report.json must be a regular file")
        record = dict(
            validate_report(
                report_path,
                sources,
                None,
                preset == "region-screen",
                expected_preset=preset,
                require_preset_provenance=require_preset_provenance,
            )
        )
        source_id = record["source_id"]
        if source_id in seen_sources:
            raise CompletionError("detector produced duplicate source reports")
        seen_sources.add(source_id)
        record["path"] = report_path.relative_to(result_root).as_posix()
        records.append(record)
    if seen_sources != set(sources):
        missing = sorted(set(sources) - seen_sources)
        extra = sorted(seen_sources - set(sources))
        raise CompletionError(
            "detector report/source set differs: missing={}, extra={}".format(
                missing[:5], extra[:5]
            )
        )
    return tuple(records)


def _completion_path(worker: WorkerConfiguration, job: Job) -> Path:
    return worker.run_dir / "samples" / _sample_slug(job.sample_id) / "completion.json"


def _validate_completion(
    completion_path: Path,
    job: Job,
    worker: WorkerConfiguration,
) -> Mapping[str, object]:
    try:
        completion = json.loads(completion_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CompletionError(
            "cannot parse completion marker: {}".format(error)
        ) from error
    if not isinstance(completion, dict):
        raise CompletionError("completion root must be an object")
    expected = {
        "schema": COMPLETION_SCHEMA,
        "version": COMPLETION_VERSION,
        "sample_id": job.sample_id,
        "source_file_sha256": job.source_sha256,
        "source_file_bytes": job.source_bytes,
        "config_sha256": worker.config_sha256,
    }
    for key, value in expected.items():
        if completion.get(key) != value:
            raise CompletionError("completion {} differs".format(key))
    reports = completion.get("reports")
    if not isinstance(reports, list) or not reports:
        raise CompletionError("completion reports must be non-empty")
    item_root = completion_path.parent
    source_ids = set()
    for report_record in reports:
        if not isinstance(report_record, dict):
            raise CompletionError("completion report record must be an object")
        target = _safe_relative_file(item_root, report_record.get("path", ""))
        if _sha256_file(target) != report_record.get("report_sha256"):
            raise CompletionError("completion report hash differs")
        validated = validate_report(
            target,
            expected_sources=None,
            expected_record=report_record,
            region_screen=worker.tool.preset == "region-screen",
            expected_preset=worker.tool.preset,
            require_preset_provenance=(worker.tool.preset_mode == "detector-preset"),
        )
        if validated["array_count"] != report_record.get("array_count"):
            raise CompletionError("completion array count differs")
        source_ids.add(validated["source_id"])
    if len(source_ids) != len(reports):
        raise CompletionError("completion contains duplicate source reports")
    if worker.resume_verify_input:
        observed_sha256, _sources = _read_fasta(job.fasta)
        if observed_sha256 != job.source_sha256:
            raise CompletionError("input FASTA hash differs during resume verification")
    return completion


def _archive_invalid_completion(path: Path) -> None:
    if not path.exists():
        return
    archive = path.parent / "superseded"
    archive.mkdir(parents=True, exist_ok=True)
    try:
        digest = _sha256_file(path)[:12]
    except OSError:
        digest = "unreadable"
    target = archive / "completion-{}-{}.json".format(time.time_ns(), digest)
    os.replace(str(path), str(target))
    _fsync_directory(path.parent)


def _existing_attempt_count(item_root: Path, config_sha256: str) -> int:
    attempts_root = item_root / "attempts"
    if not attempts_root.is_dir():
        return 0
    prefix = config_sha256[:12]
    numbers = []
    for path in attempts_root.iterdir():
        match = ATTEMPT_RE.fullmatch(path.name)
        if match and match.group(1) == prefix:
            numbers.append(int(match.group(2)))
    return max(numbers, default=0)


def terminate_process_group(process: subprocess.Popen, grace_seconds: float) -> None:
    with contextlib.suppress(ProcessLookupError):
        os.killpg(process.pid, signal.SIGTERM)
    if process.poll() is None:
        with contextlib.suppress(subprocess.TimeoutExpired):
            process.wait(timeout=grace_seconds)
    # A grandchild may ignore SIGTERM and outlive its direct parent.  Target
    # the original process group once more even when the parent has exited.
    with contextlib.suppress(ProcessLookupError):
        os.killpg(process.pid, signal.SIGKILL)
    if process.poll() is None:
        with contextlib.suppress(subprocess.TimeoutExpired):
            process.wait(timeout=max(1.0, grace_seconds))


def run_process(
    command: Sequence[str],
    cwd: Path,
    environment: Mapping[str, str],
    stdout_path: Path,
    stderr_path: Path,
    timeout_seconds: float,
    grace_seconds: float,
) -> Tuple[str, Optional[int], float]:
    start = time.monotonic()
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(
            tuple(command),
            cwd=str(cwd),
            env=dict(environment),
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        ACTIVE_PROCESSES.add(process)
        try:
            deadline = start + timeout_seconds
            status = None
            return_code = None
            while status is None:
                if ACTIVE_PROCESSES.stop_requested.is_set():
                    terminate_process_group(process, grace_seconds)
                    return_code = process.returncode
                    status = "interrupted"
                    break
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    terminate_process_group(process, grace_seconds)
                    return_code = process.returncode
                    status = "timeout"
                    break
                try:
                    return_code = process.wait(timeout=min(1.0, remaining))
                except subprocess.TimeoutExpired:
                    continue
                status = "success" if return_code == 0 else "exit_nonzero"
        finally:
            ACTIVE_PROCESSES.remove(process)
    return status, return_code, time.monotonic() - start


def _attempt_command(
    worker: WorkerConfiguration,
    job: Job,
    result_root: Path,
) -> Tuple[str, ...]:
    command = (
        str(worker.tool.runner),
        "--file",
        str(job.fasta),
        "--result_folder",
        str(result_root),
        "--cpu",
        str(worker.tool.cpus_per_worker),
        "--folder_jobs",
        "1",
        "--parallel",
        "True",
    )
    return command + worker.tool.analysis_arguments + worker.tool.extra_arguments


def _attempt_environment(worker: WorkerConfiguration) -> Dict[str, str]:
    environment = os.environ.copy()
    threads = str(worker.tool.cpus_per_worker)
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[name] = threads
    environment["TMPDIR"] = str(worker.scratch_root)
    environment["PYTHONHASHSEED"] = "0"
    environment["PYTHONNOUSERSITE"] = "1"
    environment["LC_ALL"] = "C"
    return environment


def execute_job(job: Job, worker: WorkerConfiguration) -> JobResult:
    started = time.monotonic()
    item_root = worker.run_dir / "samples" / _sample_slug(job.sample_id)
    item_root.mkdir(parents=True, exist_ok=True)
    completion_path = item_root / "completion.json"
    if completion_path.exists():
        try:
            _validate_completion(completion_path, job, worker)
        except CompletionError:
            _archive_invalid_completion(completion_path)
        else:
            return JobResult(
                job.ordinal,
                job.sample_id,
                "skipped",
                0,
                time.monotonic() - started,
                completion=completion_path.relative_to(worker.run_dir).as_posix(),
            )

    if ACTIVE_PROCESSES.stop_requested.is_set():
        return JobResult(
            job.ordinal,
            job.sample_id,
            "failed",
            0,
            time.monotonic() - started,
            "interrupted",
            "batch shutdown was requested",
        )

    existing_attempts = _existing_attempt_count(item_root, worker.config_sha256)
    if existing_attempts >= worker.max_attempts:
        return JobResult(
            job.ordinal,
            job.sample_id,
            "failed",
            existing_attempts,
            time.monotonic() - started,
            "retry_cap_exhausted",
            "{} prior attempt(s) reached --max-attempts".format(existing_attempts),
        )

    try:
        metadata = job.fasta.stat()
    except OSError as error:
        return JobResult(
            job.ordinal,
            job.sample_id,
            "failed",
            existing_attempts,
            time.monotonic() - started,
            "input_missing",
            str(error),
        )
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_size != job.source_bytes:
        return JobResult(
            job.ordinal,
            job.sample_id,
            "failed",
            existing_attempts,
            time.monotonic() - started,
            "input_size_mismatch",
            "input is not a regular file or byte count differs from manifest",
        )
    try:
        source_file_sha256, sources = _read_fasta(job.fasta)
    except InputError as error:
        return JobResult(
            job.ordinal,
            job.sample_id,
            "failed",
            existing_attempts,
            time.monotonic() - started,
            "invalid_fasta",
            str(error),
        )
    if source_file_sha256 != job.source_sha256:
        return JobResult(
            job.ordinal,
            job.sample_id,
            "failed",
            existing_attempts,
            time.monotonic() - started,
            "input_hash_mismatch",
            "input FASTA SHA-256 differs from manifest",
        )

    attempts_used = existing_attempts
    final_code = ""
    final_message = ""
    for attempt_number in range(existing_attempts + 1, worker.max_attempts + 1):
        if ACTIVE_PROCESSES.stop_requested.is_set():
            final_code = "interrupted"
            final_message = "batch shutdown was requested"
            break
        attempts_used = attempt_number
        attempt_root = (
            item_root
            / "attempts"
            / "attempt-{}-{:04d}".format(worker.config_sha256[:12], attempt_number)
        )
        attempt_root.mkdir(parents=True, exist_ok=False)
        result_root = attempt_root / "results"
        result_root.mkdir()
        command = _attempt_command(worker, job, result_root)
        running_record = {
            "schema": "CRISPRidentify-batch-attempt",
            "version": 1,
            "sample_id": job.sample_id,
            "attempt": attempt_number,
            "config_sha256": worker.config_sha256,
            "source_file_sha256": job.source_sha256,
            "command": list(command),
            "status": "running",
            "started_unix_ns": time.time_ns(),
        }
        atomic_write_json(attempt_root / "attempt.json", running_record)
        try:
            process_status, return_code, duration = run_process(
                command,
                cwd=worker.run_dir,
                environment=_attempt_environment(worker),
                stdout_path=attempt_root / "stdout.log",
                stderr_path=attempt_root / "stderr.log",
                timeout_seconds=worker.timeout_seconds,
                grace_seconds=worker.termination_grace_seconds,
            )
            reports = None
            if process_status == "success":
                try:
                    raw_reports = collect_and_validate_reports(
                        result_root,
                        sources,
                        worker.tool.preset,
                        worker.tool.preset_mode == "detector-preset",
                    )
                    attempt_relative = attempt_root.relative_to(item_root)
                    reports = tuple(
                        dict(
                            record,
                            path=(
                                attempt_relative / "results" / record["path"]
                            ).as_posix(),
                        )
                        for record in raw_reports
                    )
                except CompletionError as error:
                    process_status = "invalid_output"
                    final_code = "invalid_output"
                    final_message = str(error)
            elif process_status == "timeout":
                final_code = "timeout"
                final_message = "detector exceeded {:.3f} seconds".format(
                    worker.timeout_seconds
                )
            elif process_status == "interrupted":
                final_code = "interrupted"
                final_message = "batch shutdown was requested"
            else:
                final_code = "exit_nonzero"
                final_message = "detector exited with status {}".format(return_code)

            finished_record = dict(running_record)
            finished_record.update(
                {
                    "status": process_status,
                    "return_code": return_code,
                    "elapsed_seconds": duration,
                    "finished_unix_ns": time.time_ns(),
                    "stdout_sha256": _sha256_file(attempt_root / "stdout.log"),
                    "stderr_sha256": _sha256_file(attempt_root / "stderr.log"),
                    "error_code": final_code,
                    "error_message": final_message,
                }
            )
            atomic_write_json(attempt_root / "attempt.json", finished_record)
            if process_status == "success" and reports is not None:
                completion = {
                    "schema": COMPLETION_SCHEMA,
                    "version": COMPLETION_VERSION,
                    "sample_id": job.sample_id,
                    "source_file": str(job.fasta),
                    "source_file_bytes": job.source_bytes,
                    "source_file_sha256": job.source_sha256,
                    "config_sha256": worker.config_sha256,
                    "successful_attempt": attempt_number,
                    "reports": list(reports),
                    "completed_unix_ns": time.time_ns(),
                }
                atomic_write_json(completion_path, completion)
                _validate_completion(completion_path, job, worker)
                return JobResult(
                    job.ordinal,
                    job.sample_id,
                    "succeeded",
                    attempts_used,
                    time.monotonic() - started,
                    completion=completion_path.relative_to(worker.run_dir).as_posix(),
                )
        except (OSError, subprocess.SubprocessError, BatchError) as error:
            final_code = "runner_error"
            final_message = str(error)
            failed_record = dict(running_record)
            failed_record.update(
                {
                    "status": "runner_error",
                    "finished_unix_ns": time.time_ns(),
                    "error_code": final_code,
                    "error_message": final_message,
                }
            )
            atomic_write_json(attempt_root / "attempt.json", failed_record)

    return JobResult(
        job.ordinal,
        job.sample_id,
        "failed",
        attempts_used,
        time.monotonic() - started,
        final_code or "retry_cap_exhausted",
        final_message or "all attempts failed",
    )


def _update_job(connection: sqlite3.Connection, result: JobResult) -> None:
    connection.execute(
        "UPDATE jobs SET status=?,attempts=?,elapsed_seconds=?,error_code=?,"
        "error_message=?,completion=? WHERE ordinal=?",
        (
            result.status,
            result.attempts,
            result.elapsed_seconds,
            result.error_code,
            result.error_message,
            result.completion,
            result.ordinal,
        ),
    )
    connection.commit()


def run_bounded(
    jobs: Iterable[Job],
    worker: WorkerConfiguration,
    workers: int,
    state_db: Path,
    max_in_flight: Optional[int] = None,
) -> None:
    """Execute a streaming job iterator with at most a fixed Future window."""

    if workers < 1:
        raise ValueError("workers must be positive")
    if max_in_flight is None:
        max_in_flight = workers * 2
    if max_in_flight < workers:
        raise ValueError("max_in_flight must be at least workers")
    connection = sqlite3.connect(str(state_db))
    iterator = iter(jobs)
    exhausted = False
    pending: Dict[concurrent.futures.Future, Job] = {}
    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            while pending or not exhausted:
                while not exhausted and len(pending) < max_in_flight:
                    if ACTIVE_PROCESSES.stop_requested.is_set():
                        exhausted = True
                        break
                    try:
                        job = next(iterator)
                    except StopIteration:
                        exhausted = True
                        break
                    future = executor.submit(execute_job, job, worker)
                    pending[future] = job
                if not pending:
                    continue
                done, _not_done = concurrent.futures.wait(
                    pending, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for future in sorted(done, key=lambda value: pending[value].ordinal):
                    job = pending.pop(future)
                    try:
                        result = future.result()
                    except BaseException as error:
                        result = JobResult(
                            job.ordinal,
                            job.sample_id,
                            "failed",
                            0,
                            0.0,
                            "worker_exception",
                            "{}: {}".format(type(error).__name__, error),
                        )
                    _update_job(connection, result)
    finally:
        connection.close()


def _write_ledgers(
    state_db: Path,
    run_dir: Path,
    summary_base: Mapping[str, object],
    interrupted: bool,
) -> Mapping[str, object]:
    connection = sqlite3.connect(str(state_db))
    counts = {
        row[0]: row[1]
        for row in connection.execute(
            "SELECT status,COUNT(*) FROM jobs GROUP BY status"
        )
    }
    status_tmp = run_dir / ".status.tsv.tmp"
    failures_tmp = run_dir / ".failures.jsonl.tmp"
    with (
        status_tmp.open("w", encoding="utf-8", newline="") as status_handle,
        failures_tmp.open("w", encoding="utf-8", newline="") as failure_handle,
    ):
        status_writer = csv.writer(status_handle, delimiter="\t", lineterminator="\n")
        status_writer.writerow(STATUS_FIELDS)
        rows = connection.execute(
            "SELECT ordinal,sample_id,status,attempts,elapsed_seconds,error_code,"
            "error_message,completion FROM jobs ORDER BY ordinal"
        )
        for row in rows:
            status_writer.writerow(row)
            if row[2] == "failed":
                failure_handle.write(
                    json.dumps(
                        {
                            "ordinal": row[0],
                            "sample_id": row[1],
                            "attempts": row[3],
                            "elapsed_seconds": row[4],
                            "error_code": row[5],
                            "error_message": row[6],
                        },
                        allow_nan=False,
                        sort_keys=True,
                    )
                    + "\n"
                )
        status_handle.flush()
        failure_handle.flush()
        os.fsync(status_handle.fileno())
        os.fsync(failure_handle.fileno())
    connection.close()
    os.replace(str(status_tmp), str(run_dir / "status.tsv"))
    os.replace(str(failures_tmp), str(run_dir / "failures.jsonl"))
    _fsync_directory(run_dir)

    finished_ns = time.time_ns()
    started_ns = summary_base.get("started_unix_ns")
    elapsed_seconds = (
        max(0.0, (finished_ns - started_ns) / 1_000_000_000.0)
        if isinstance(started_ns, int) and not isinstance(started_ns, bool)
        else None
    )
    selected = sum(counts.values())
    fresh_succeeded = counts.get("succeeded", 0)
    completed = fresh_succeeded + counts.get("skipped", 0)
    fresh_rate = (
        fresh_succeeded * 3600.0 / elapsed_seconds
        if elapsed_seconds is not None and elapsed_seconds > 0
        else None
    )
    summary = dict(summary_base)
    summary.update(
        {
            "status": "interrupted"
            if interrupted
            else (
                "failed"
                if counts.get("failed", 0) or counts.get("pending", 0)
                else "complete"
            ),
            "counts": {
                "selected": selected,
                "succeeded": fresh_succeeded,
                "skipped": counts.get("skipped", 0),
                "failed": counts.get("failed", 0),
                "pending": counts.get("pending", 0),
            },
            "finished_unix_ns": finished_ns,
            "elapsed_seconds": elapsed_seconds,
            # Resume skips prove completeness but consume almost no detector time;
            # excluding them prevents a resume from reporting fictitious throughput.
            "completed_genomes_per_hour": fresh_rate,
            "throughput_basis": "fresh_succeeded_only_excludes_resume_skips",
            "completion_fraction": completed / selected if selected else 1.0,
            "status_ledger_sha256": _sha256_file(run_dir / "status.tsv"),
            "failure_ledger_sha256": _sha256_file(run_dir / "failures.jsonl"),
        }
    )
    atomic_write_json(run_dir / "summary.json", summary)
    return summary


@contextlib.contextmanager
def exclusive_run_lock(run_dir: Path) -> Iterator[None]:
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / ".batch.lock"
    handle = path.open("a+")
    try:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise BatchError("another batch process holds {}".format(path)) from error
        yield
    finally:
        with contextlib.suppress(OSError):
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()


def _available_cpus() -> int:
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def _check_resources(args: argparse.Namespace) -> None:
    requested = args.workers * args.cpus_per_worker
    available = _available_cpus()
    if requested > available and not args.allow_oversubscribe:
        raise InputError(
            "workers * cpus-per-worker requests {} CPUs but affinity exposes {}; "
            "reduce concurrency or pass --allow-oversubscribe explicitly".format(
                requested, available
            )
        )
    for path, label in (
        (args.run_dir, "run directory"),
        (args.scratch_root, "scratch"),
    ):
        path.mkdir(parents=True, exist_ok=True)
        free_gib = shutil.disk_usage(str(path)).free / (1024**3)
        if free_gib < args.min_free_gib:
            raise InputError(
                "{} has {:.2f} GiB free; --min-free-gib requires {:.2f}".format(
                    label, free_gib, args.min_free_gib
                )
            )


def _configuration_record(
    args: argparse.Namespace,
    tool: ToolConfiguration,
    manifest_sha256: str,
    total_rows: int,
    selected_rows: int,
) -> Mapping[str, object]:
    return {
        "schema": BATCH_SCHEMA,
        "version": BATCH_VERSION,
        "manifest": {
            "path": str(args.manifest),
            "sha256": manifest_sha256,
            "total_rows": total_rows,
            "selected_rows": selected_rows,
        },
        "shard": {"count": args.shard_count, "index": args.shard_index},
        "execution": {
            "workers": args.workers,
            "cpus_per_worker": args.cpus_per_worker,
            "timeout_seconds": args.timeout_seconds,
            "termination_grace_seconds": args.termination_grace_seconds,
            "max_attempts": args.max_attempts,
            "resume_verify_input": args.resume_verify_input,
            "scratch_root": str(args.scratch_root),
        },
        "tool": _tool_config_record(tool),
    }


def _install_signal_handlers(grace_seconds: float):
    previous = {}

    def handler(signum, _frame):
        ACTIVE_PROCESSES.stop_requested.set()
        # Signal-handler work stays intentionally bounded.  Worker waits observe
        # terminated children and the scheduler stops submitting new jobs.
        with ACTIVE_PROCESSES._lock:
            processes = tuple(ACTIVE_PROCESSES._processes.values())
        for process in processes:
            with contextlib.suppress(ProcessLookupError):
                os.killpg(process.pid, signal.SIGTERM)

    for signum in (signal.SIGINT, signal.SIGTERM):
        previous[signum] = signal.getsignal(signum)
        signal.signal(signum, handler)
    return previous


def _restore_signal_handlers(previous: Mapping[int, object]) -> None:
    for signum, handler in previous.items():
        signal.signal(signum, handler)


def command_run(args: argparse.Namespace, project_root: Path) -> int:
    args.manifest = args.manifest.expanduser().absolute()
    args.run_dir = args.run_dir.expanduser().absolute()
    args.scratch_root = args.scratch_root.expanduser().absolute()
    _check_resources(args)
    started_ns = time.time_ns()
    ACTIVE_PROCESSES.stop_requested.clear()
    previous_handlers = _install_signal_handlers(args.termination_grace_seconds)
    try:
        with exclusive_run_lock(args.run_dir):
            state_db = args.run_dir / "state.sqlite"
            candidate_state = args.run_dir / ".state-{}.sqlite".format(os.getpid())
            with contextlib.suppress(FileNotFoundError):
                candidate_state.unlink()
            try:
                manifest_sha256, total_rows, selected_rows = index_manifest(
                    args.manifest,
                    candidate_state,
                    args.shard_count,
                    args.shard_index,
                )
                tool = build_tool_configuration(args, project_root)
                configuration = _configuration_record(
                    args, tool, manifest_sha256, total_rows, selected_rows
                )
                config_sha256 = _canonical_sha256(configuration)
                configuration = dict(configuration)
                configuration["config_sha256"] = config_sha256
                config_path = args.run_dir / "config.json"
                if config_path.exists():
                    try:
                        existing = json.loads(config_path.read_text(encoding="utf-8"))
                    except (OSError, json.JSONDecodeError) as error:
                        raise BatchError(
                            "existing config.json is invalid: {}".format(error)
                        )
                    if existing != configuration:
                        raise BatchError(
                            "run directory belongs to a different configuration; "
                            "use a new run directory instead of mixing results"
                        )
                else:
                    atomic_write_json(config_path, configuration)
                os.replace(str(candidate_state), str(state_db))
                _fsync_directory(args.run_dir)
            finally:
                with contextlib.suppress(FileNotFoundError):
                    candidate_state.unlink()

            worker = WorkerConfiguration(
                run_dir=args.run_dir,
                scratch_root=args.scratch_root,
                tool=tool,
                config_sha256=config_sha256,
                timeout_seconds=args.timeout_seconds,
                termination_grace_seconds=args.termination_grace_seconds,
                max_attempts=args.max_attempts,
                resume_verify_input=args.resume_verify_input,
            )
            run_bounded(
                iter_indexed_jobs(state_db),
                worker,
                args.workers,
                state_db,
                max_in_flight=args.max_in_flight,
            )
            interrupted = ACTIVE_PROCESSES.stop_requested.is_set()
            summary = _write_ledgers(
                state_db,
                args.run_dir,
                {
                    "schema": BATCH_SCHEMA + "-summary",
                    "version": 1,
                    "config_sha256": config_sha256,
                    "manifest_sha256": manifest_sha256,
                    "shard_count": args.shard_count,
                    "shard_index": args.shard_index,
                    "started_unix_ns": started_ns,
                },
                interrupted,
            )
            return 0 if summary["status"] == "complete" else 1
    finally:
        if ACTIVE_PROCESSES.stop_requested.is_set():
            ACTIVE_PROCESSES.terminate_all(args.termination_grace_seconds)
        _restore_signal_handlers(previous_handlers)


def build_parser(project_root: Path) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="restartable, bounded CRISPRidentify cohort execution"
    )
    parser.add_argument(
        "--version",
        action="version",
        version="crispridentify-batch {} (CRISPRidentify {})".format(
            BATCH_VERSION, _runtime_version(project_root)
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    manifest = subparsers.add_parser(
        "manifest", help="create a sorted checksummed FASTA manifest"
    )
    manifest.add_argument("--input-root", type=Path, required=True)
    manifest.add_argument("--output", type=Path, required=True)
    manifest.add_argument(
        "--suffix",
        action="append",
        default=None,
        help="accepted suffix; repeat as needed (default: .fna/.fa/.fasta/.fas)",
    )

    run = subparsers.add_parser("run", help="run or safely resume one manifest shard")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--run-dir", type=Path, required=True)
    run.add_argument("--runner", type=Path, default=project_root / "crispridentify")
    run.add_argument("--preset", choices=PRESETS, default="region-screen")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--cpus-per-worker", type=int, default=2)
    run.add_argument("--max-in-flight", type=int, default=None)
    run.add_argument("--timeout-seconds", type=float, default=1800.0)
    run.add_argument("--termination-grace-seconds", type=float, default=10.0)
    run.add_argument("--max-attempts", type=int, default=3)
    run.add_argument("--shard-count", type=int, default=1)
    run.add_argument("--shard-index", type=int, default=0)
    run.add_argument(
        "--scratch-root",
        type=Path,
        default=Path(os.environ.get("TMPDIR", tempfile.gettempdir())),
    )
    run.add_argument("--min-free-gib", type=float, default=10.0)
    run.add_argument("--allow-oversubscribe", action="store_true")
    run.add_argument("--resume-verify-input", action="store_true")
    run.add_argument("--cas", action="store_true")
    run.add_argument(
        "--cas-gene-mode", choices=("auto", "single", "meta"), default="auto"
    )
    run.add_argument("--cas-timeout-seconds", type=int, default=3600)
    run.add_argument("--sequence-topology-manifest", type=Path)
    run.add_argument(
        "--tool-arg",
        action="append",
        default=[],
        help="extra detector argv token; repeat and use --tool-arg=--option",
    )
    return parser


def _positive(value: object, name: str, allow_zero: bool = False) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise InputError(
            "--{} must be {}".format(name, "non-negative" if allow_zero else "positive")
        )
    if (allow_zero and value < 0) or (not allow_zero and value <= 0):
        raise InputError(
            "--{} must be {}".format(name, "non-negative" if allow_zero else "positive")
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    project_root = Path(__file__).resolve().parents[1]
    parser = build_parser(project_root)
    args = parser.parse_args(argv)
    try:
        if args.command == "manifest":
            suffixes = args.suffix or (".fna", ".fa", ".fasta", ".fas")
            count = create_manifest(args.input_root, args.output, suffixes)
            print("wrote {} FASTA file(s) to {}".format(count, args.output))
            return 0

        for value, name in (
            (args.workers, "workers"),
            (args.cpus_per_worker, "cpus-per-worker"),
            (args.timeout_seconds, "timeout-seconds"),
            (args.termination_grace_seconds, "termination-grace-seconds"),
            (args.max_attempts, "max-attempts"),
            (args.shard_count, "shard-count"),
            (args.cas_timeout_seconds, "cas-timeout-seconds"),
        ):
            _positive(value, name)
        _positive(args.shard_index, "shard-index", allow_zero=True)
        _positive(args.min_free_gib, "min-free-gib", allow_zero=True)
        if args.max_in_flight is not None:
            _positive(args.max_in_flight, "max-in-flight")
            if args.max_in_flight < args.workers:
                raise InputError("--max-in-flight must be at least --workers")
        if args.shard_index >= args.shard_count:
            raise InputError("--shard-index must be smaller than --shard-count")
        return command_run(args, project_root)
    except (BatchError, OSError, sqlite3.Error) as error:
        print("crispridentify-batch: error: {}".format(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
