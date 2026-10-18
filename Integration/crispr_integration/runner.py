"""Fail-fast, shell-free orchestration for the three isolated tool workers."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field as dataclass_field
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import subprocess
import time
from typing import Iterable, Mapping, Sequence

from .adapter import adapt_reports
from .parser import normalize_category
from .scientific_options import tool_arguments


class StageError(RuntimeError):
    def __init__(self, stage: str, message: str):
        self.stage = stage
        super().__init__(f"{stage}: {message}")


@dataclass(frozen=True)
class RunnerConfig:
    crispridentify: tuple[str, ...]
    spacerplacer: tuple[str, ...]
    crispr_evor: tuple[str, ...]
    casandra: tuple[str, ...] = ("casandra",)
    crispr_tracrrna: tuple[str, ...] = ("crisprtracrrna",)
    crisprleader: tuple[str, ...] = ("crisprleader",)
    crisprrepeat: tuple[str, ...] = ("crisprrepeat",)
    crisprspacer: tuple[str, ...] = ("crisprspacer",)
    cctk: tuple[str, ...] = ("cctk",)
    clinker: tuple[str, ...] = ("clinker-v2",)
    phigaro: tuple[str, ...] = ("phigaro-v2", "predict")
    isescan: tuple[str, ...] = ("isescan-v2",)
    defensefinder: tuple[str, ...] = ("defense-context",)
    rnaalifold: tuple[str, ...] = ("rnaalifold-v2",)
    spacepharer: tuple[str, ...] = ("spacepharer",)
    repeattyper: tuple[str, ...] = ("repeattyper",)
    crisprmap: tuple[str, ...] = ("crisprmap",)
    blastn: tuple[str, ...] = ("blastn",)
    tool_options: dict = dataclass_field(default_factory=dict)


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_runner_config() -> RunnerConfig:
    root = _project_root()
    return RunnerConfig(
        crispridentify=(str(root / "tools" / "CRISPRidentify-v2" / "crispridentify"),),
        spacerplacer=("spacerplacer",),
        crispr_evor=("crispr-evor",),
    )


def _command(value: object, name: str, base: Path) -> tuple[str, ...]:
    if not isinstance(value, list) or not value or not all(isinstance(item, str) and item for item in value):
        raise ValueError(f"runner config field {name!r} must be a non-empty JSON string array")
    if any("\0" in item for item in value):
        raise ValueError(f"runner config field {name!r} contains a NUL byte")
    command = list(value)
    first = Path(command[0]).expanduser()
    if not first.is_absolute() and ("/" in command[0] or "\\" in command[0]):
        command[0] = str((base / first).resolve())
    return tuple(command)


def load_runner_config(path: str | Path | None) -> RunnerConfig:
    if path is None:
        return default_runner_config()
    config_path = Path(path).expanduser().resolve()
    try:
        with config_path.open("r", encoding="utf-8") as stream:
            document = json.load(stream)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read runner config {config_path}: {error}") from error
    if not isinstance(document, Mapping):
        raise ValueError("runner config must be a JSON object")
    return RunnerConfig(
        crispridentify=_command(document.get("crispridentify"), "crispridentify", config_path.parent),
        spacerplacer=_command(document.get("spacerplacer"), "spacerplacer", config_path.parent),
        crispr_evor=_command(document.get("crispr_evor"), "crispr_evor", config_path.parent),
        crisprleader=_command(document.get("crisprleader", ["crisprleader"]), "crisprleader", config_path.parent),
        crisprrepeat=_command(document.get("crisprrepeat", ["crisprrepeat"]), "crisprrepeat", config_path.parent),
        cctk=_command(document.get("cctk", ["cctk"]), "cctk", config_path.parent),
        clinker=_command(document.get("clinker", ["clinker-v2"]), "clinker", config_path.parent),
        phigaro=_command(document.get("phigaro", ["phigaro-v2", "predict"]), "phigaro", config_path.parent),
        isescan=_command(document.get("isescan", ["isescan-v2"]), "isescan", config_path.parent),
        defensefinder=_command(document.get("defensefinder", ["defense-context"]), "defensefinder", config_path.parent),
        rnaalifold=_command(document.get("rnaalifold", ["rnaalifold-v2"]), "rnaalifold", config_path.parent),
        spacepharer=_command(document.get("spacepharer", ["spacepharer"]), "spacepharer", config_path.parent),
        repeattyper=_command(document.get("repeattyper", ["repeattyper"]), "repeattyper", config_path.parent),
        crisprmap=_command(document.get("crisprmap", ["crisprmap"]), "crisprmap", config_path.parent),
        blastn=_command(document.get("blastn", ["blastn"]), "blastn", config_path.parent),
        crisprspacer=_command(document.get("crisprspacer", ["crisprspacer"]), "crisprspacer", config_path.parent),
        casandra=_command(document.get("casandra", ["casandra"]), "casandra", config_path.parent),
        crispr_tracrrna=_command(document.get("crispr_tracrrna", ["crisprtracrrna"]), "crispr_tracrrna", config_path.parent),
    )


def _atomic_json(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    os.replace(temporary, path)


def _empty_job_root(path: str | Path) -> Path:
    root = Path(path).expanduser().resolve()
    if root.exists() and not root.is_dir():
        raise FileExistsError(f"job output exists and is not a directory: {root}")
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"job output directory is not empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    return root


def _resolve_program(command: Sequence[str]) -> None:
    program = command[0]
    if "/" in program or "\\" in program:
        path = Path(program)
        if not path.is_file():
            raise FileNotFoundError(f"worker executable does not exist: {path}")
        if not os.access(path, os.X_OK):
            raise PermissionError(f"worker executable is not executable: {path}")
    elif shutil.which(program) is None:
        raise FileNotFoundError(f"worker executable is not on PATH: {program}")


def _run_process(
    stage: str,
    command: Sequence[str],
    logs: Path,
    *,
    timeout_seconds: float,
) -> dict[str, object]:
    _resolve_program(command)
    stdout_path = logs / f"{stage}.stdout.log"
    stderr_path = logs / f"{stage}.stderr.log"
    started = time.monotonic()
    with stdout_path.open("x", encoding="utf-8") as stdout, stderr_path.open("x", encoding="utf-8") as stderr:
        try:
            # Commands come only from the operator-owned runner config and are
            # passed as an argv array with no shell interpretation.
            process = subprocess.Popen(  # noqa: S603
                list(command),
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
        except OSError as error:
            raise StageError(stage, f"could not start worker: {error}") from error
        timed_out = False
        try:
            return_code = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                return_code = process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                return_code = process.wait()
    result: dict[str, object] = {
        "command": list(command),
        "duration_seconds": round(time.monotonic() - started, 6),
        "return_code": return_code,
        "status": "timed_out" if timed_out else ("completed" if return_code == 0 else "failed"),
        "stderr": str(stderr_path),
        "stdout": str(stdout_path),
    }
    if timed_out:
        raise StageError(stage, f"worker exceeded {timeout_seconds:g} seconds")
    if return_code != 0:
        raise StageError(stage, f"worker exited with status {return_code}; see {stderr_path}")
    return result


def _require_nonempty_regular_file(stage: str, path: Path, description: str) -> None:
    if path.is_symlink():
        raise StageError(stage, f"required {description} must not be a symbolic link: {path}")
    try:
        is_file = path.is_file()
        size = path.stat().st_size if is_file else None
    except OSError as error:
        raise StageError(stage, f"cannot inspect required {description} {path}: {error}") from error
    if not is_file:
        raise StageError(stage, f"required {description} is missing or not a regular file: {path}")
    if size == 0:
        raise StageError(stage, f"required {description} is empty: {path}")


def _read_json_object(stage: str, path: Path, description: str) -> Mapping[str, object]:
    _require_nonempty_regular_file(stage, path, description)
    try:
        with path.open("r", encoding="utf-8") as stream:
            document = json.load(stream)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise StageError(stage, f"required {description} is not valid JSON: {path}: {error}") from error
    if not isinstance(document, Mapping):
        raise StageError(stage, f"required {description} must contain a JSON object: {path}")
    return document


def _nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _finite_number(value: object) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value)


def _read_result_csv(
    stage: str,
    path: Path,
    description: str,
    *,
    required_fields: set[str],
) -> dict[str, dict[str, str]]:
    """Read a worker CSV and return one validated row per unique group."""
    _require_nonempty_regular_file(stage, path, description)
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            fieldnames = reader.fieldnames or []
            if len(fieldnames) != len(set(fieldnames)):
                raise StageError(stage, f"required {description} contains duplicate columns")
            missing_fields = required_fields.difference(fieldnames)
            if missing_fields:
                raise StageError(
                    stage,
                    f"required {description} is missing columns: {sorted(missing_fields)}",
                )

            rows: dict[str, dict[str, str]] = {}
            for line_number, row in enumerate(reader, start=2):
                if None in row:
                    raise StageError(stage, f"required {description} has extra fields on line {line_number}")
                group = (row.get("name") or "").strip()
                if not group:
                    raise StageError(stage, f"required {description} has an empty group on line {line_number}")
                if group in rows:
                    raise StageError(stage, f"required {description} has duplicate group {group!r}")
                rows[group] = {key: value or "" for key, value in row.items()}
    except StageError:
        raise
    except (OSError, UnicodeError, csv.Error) as error:
        raise StageError(stage, f"required {description} is not valid CSV: {path}: {error}") from error
    return rows


def _require_exact_groups(
    stage: str,
    actual: set[str],
    expected: set[str],
    description: str,
) -> None:
    missing = expected.difference(actual)
    unexpected = actual.difference(expected)
    if missing or unexpected:
        details = []
        if missing:
            details.append(f"missing {sorted(missing)}")
        if unexpected:
            details.append(f"unexpected {sorted(unexpected)}")
        raise StageError(stage, f"{description} groups disagree with the adapter manifest: {'; '.join(details)}")


def _adapter_group_leaf_names(manifest: Mapping[str, object]) -> dict[str, set[str]]:
    stage = "adapter"
    groups = manifest.get("groups")
    emitted_count = manifest.get("emitted_group_count")
    if not isinstance(groups, list) or isinstance(emitted_count, bool) or not isinstance(emitted_count, int):
        raise StageError(stage, "adapter manifest has invalid group metadata")

    group_leaves: dict[str, set[str]] = {}
    for index, entry in enumerate(groups):
        if not isinstance(entry, Mapping) or not _nonempty_string(entry.get("name")):
            raise StageError(stage, f"adapter manifest groups[{index}].name must be a non-empty string")
        name = str(entry["name"]).strip()
        if name in group_leaves:
            raise StageError(stage, f"adapter manifest contains duplicate group {name!r}")
        arrays = entry.get("arrays")
        array_count = entry.get("array_count")
        if (
            not isinstance(arrays, list)
            or isinstance(array_count, bool)
            or not isinstance(array_count, int)
            or array_count < 2
        ):
            raise StageError(stage, f"adapter manifest group {name!r} has invalid array metadata")
        leaves: set[str] = set()
        for array_index, array in enumerate(arrays):
            if not isinstance(array, Mapping) or not _nonempty_string(array.get("ccdb_file")):
                raise StageError(
                    stage,
                    f"adapter manifest group {name!r} arrays[{array_index}].ccdb_file is invalid",
                )
            relative = PurePosixPath(str(array["ccdb_file"]))
            if (
                relative.is_absolute()
                or len(relative.parts) != 4
                or relative.parts[:2] != ("ccdb", name)
                or relative.parts[2] not in {"pos_strand", "neg_strand"}
                or relative.suffix != ".fasta"
            ):
                raise StageError(
                    stage,
                    f"adapter manifest group {name!r} arrays[{array_index}].ccdb_file is unsafe",
                )
            leaf = relative.stem
            if _NEWICK_LABEL.fullmatch(leaf) is None:
                raise StageError(stage, f"adapter manifest group {name!r} has an unsafe tree leaf name")
            if leaf in leaves:
                raise StageError(stage, f"adapter manifest group {name!r} has duplicate tree leaf {leaf!r}")
            leaves.add(leaf)
        if array_count != len(arrays) or array_count != len(leaves):
            raise StageError(
                stage,
                f"adapter manifest group {name!r} array_count disagrees with its array entries",
            )
        group_leaves[name] = leaves
    if emitted_count != len(group_leaves):
        raise StageError(
            stage,
            f"adapter emitted_group_count={emitted_count} disagrees with its {len(group_leaves)} group entries",
        )
    return group_leaves


def _csv_finite_float(stage: str, row: Mapping[str, str], group: str, field: str) -> float:
    raw_value = row.get(field)
    try:
        value = float(raw_value) if raw_value is not None else math.nan
    except (TypeError, ValueError) as error:
        raise StageError(stage, f"result CSV field {field!r} is invalid for group {group!r}") from error
    if not math.isfinite(value):
        raise StageError(stage, f"result CSV field {field!r} is non-finite for group {group!r}")
    return value


def _csv_nonnegative_float(stage: str, row: Mapping[str, str], group: str, field: str) -> float:
    value = _csv_finite_float(stage, row, group, field)
    if value < 0:
        raise StageError(
            stage,
            f"result CSV field {field!r} must be non-negative for group {group!r}",
        )
    return value


def _csv_nonnegative_int(stage: str, row: Mapping[str, str], group: str, field: str) -> int:
    raw_value = (row.get(field) or "").strip()
    if re.fullmatch(r"[0-9]+", raw_value) is None:
        raise StageError(
            stage,
            f"result CSV field {field!r} must be a non-negative integer for group {group!r}",
        )
    return int(raw_value)


_NEWICK_LABEL = re.compile(r"[A-Za-z0-9_.+|-]+")
_NEWICK_NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?")
_NEWICK_BRANCH_ROUNDING_ERROR = 5.0001e-6


@dataclass(frozen=True)
class _NewickSummary:
    leaves: tuple[str, ...]
    total_branch_length: float
    explicit_branch_count: int


def _validate_generated_newick(stage: str, value: object, description: str) -> _NewickSummary:
    """Validate the restricted Newick emitted by the bundled tree writer."""
    if not _nonempty_string(value):
        raise StageError(stage, f"{description} must be a non-empty Newick string")
    tree = str(value).strip()
    position = 0
    branch_lengths: list[float] = []

    def skip_whitespace() -> None:
        nonlocal position
        while position < len(tree) and tree[position].isspace():
            position += 1

    def parse_label(*, required: bool) -> str | None:
        nonlocal position
        skip_whitespace()
        match = _NEWICK_LABEL.match(tree, position)
        if match is None:
            if required:
                raise StageError(stage, f"{description} has a missing or unsafe node label")
            return None
        position = match.end()
        return match.group(0)

    def parse_branch_length() -> None:
        nonlocal position
        skip_whitespace()
        if position >= len(tree) or tree[position] != ":":
            return
        position += 1
        skip_whitespace()
        match = _NEWICK_NUMBER.match(tree, position)
        if match is None:
            raise StageError(stage, f"{description} has an invalid branch length")
        branch_length = float(match.group(0))
        if not math.isfinite(branch_length):
            raise StageError(stage, f"{description} has a non-finite branch length")
        if branch_length < 0:
            raise StageError(stage, f"{description} has a negative branch length")
        branch_lengths.append(branch_length)
        position = match.end()

    def parse_subtree() -> list[str]:
        nonlocal position
        skip_whitespace()
        if position >= len(tree):
            raise StageError(stage, f"{description} ends inside a subtree")
        if tree[position] == "(":
            position += 1
            leaves = parse_subtree()
            child_count = 1
            while True:
                skip_whitespace()
                if position < len(tree) and tree[position] == ",":
                    position += 1
                    leaves.extend(parse_subtree())
                    child_count += 1
                    continue
                break
            if child_count < 2 or position >= len(tree) or tree[position] != ")":
                raise StageError(stage, f"{description} has an invalid internal node")
            position += 1
            parse_label(required=False)
        else:
            leaf = parse_label(required=True)
            leaves = [str(leaf)]
        parse_branch_length()
        return leaves

    leaves = parse_subtree()
    skip_whitespace()
    if position >= len(tree) or tree[position] != ";":
        raise StageError(stage, f"{description} is missing its terminal semicolon")
    position += 1
    skip_whitespace()
    if position != len(tree):
        raise StageError(stage, f"{description} has trailing content")
    seen_leaves: set[str] = set()
    duplicate_leaves: set[str] = set()
    for leaf in leaves:
        if leaf in seen_leaves:
            duplicate_leaves.add(leaf)
        seen_leaves.add(leaf)
    if duplicate_leaves:
        raise StageError(stage, f"{description} has duplicate leaf labels: {sorted(duplicate_leaves)}")
    return _NewickSummary(
        leaves=tuple(leaves),
        total_branch_length=math.fsum(branch_lengths),
        explicit_branch_count=len(branch_lengths),
    )


def _validate_tree_leaf_identity(
    stage: str,
    value: object,
    description: str,
    expected_leaves: set[str],
) -> _NewickSummary:
    summary = _validate_generated_newick(stage, value, description)
    actual_leaves = set(summary.leaves)
    missing = expected_leaves.difference(actual_leaves)
    unexpected = actual_leaves.difference(expected_leaves)
    if missing or unexpected:
        details = []
        if missing:
            details.append(f"missing {sorted(missing)}")
        if unexpected:
            details.append(f"unexpected {sorted(unexpected)}")
        raise StageError(
            stage,
            f"{description} leaf identities disagree with the adapter manifest: {'; '.join(details)}",
        )
    if len(summary.leaves) < 2:
        raise StageError(stage, f"{description} must contain at least two leaves")
    return summary


def _validate_reported_tree_length(
    stage: str,
    reported: float,
    summary: _NewickSummary,
    group: str,
) -> None:
    # The bundled ETE writer serializes each branch to five decimal places,
    # while the CSV retains the full-precision pre-serialization sum.  Bound
    # the accepted discrepancy by the worst-case cumulative rounding error.
    tolerance = _NEWICK_BRANCH_ROUNDING_ERROR * summary.explicit_branch_count + 1e-12
    if abs(reported - summary.total_branch_length) > tolerance:
        raise StageError(
            stage,
            f"result CSV field 'tree length' disagrees with the Newick branch sum "
            f"for group {group!r}",
        )


def _validate_orientation_result(entry: object, index: int) -> str:
    stage = "crispr_evor"
    location = f"results[{index}]"
    if not isinstance(entry, Mapping):
        raise StageError(stage, f"CRISPR-evOr {location} must be a JSON object")

    group = entry.get("group")
    if not _nonempty_string(group):
        raise StageError(stage, f"CRISPR-evOr {location}.group must be a non-empty string")

    prediction = entry.get("prediction")
    if not isinstance(prediction, str) or prediction not in ("Forward", "Reverse", "ND"):
        raise StageError(stage, f"CRISPR-evOr {location}.prediction is invalid: {prediction!r}")

    expected_reverse = {"Forward": False, "Reverse": True, "ND": None}[prediction]
    if entry.get("recommended_reverse") is not expected_reverse:
        raise StageError(
            stage,
            f"CRISPR-evOr {location}.recommended_reverse disagrees with prediction {prediction!r}",
        )

    decisive = entry.get("decisive")
    if not isinstance(decisive, bool) or decisive != (prediction != "ND"):
        raise StageError(stage, f"CRISPR-evOr {location}.decisive disagrees with prediction {prediction!r}")

    threshold = entry.get("confidence_threshold")
    if not _finite_number(threshold) or threshold < 0:
        raise StageError(stage, f"CRISPR-evOr {location}.confidence_threshold must be finite and non-negative")

    for field in (
        "forward_ln_likelihood_bdm",
        "reverse_ln_likelihood_bdm",
        "forward_minus_reverse_ln_likelihood_bdm",
    ):
        if not _finite_number(entry.get(field)):
            raise StageError(stage, f"CRISPR-evOr {location}.{field} must be a finite number")

    forward = float(entry["forward_ln_likelihood_bdm"])
    reverse = float(entry["reverse_ln_likelihood_bdm"])
    difference = float(entry["forward_minus_reverse_ln_likelihood_bdm"])
    if not math.isclose(difference, forward - reverse, rel_tol=1e-9, abs_tol=1e-9):
        raise StageError(stage, f"CRISPR-evOr {location} has an inconsistent likelihood difference")
    expected_prediction = (
        "Reverse" if difference < -float(threshold)
        else "Forward" if difference > float(threshold)
        else "ND"
    )
    if prediction != expected_prediction:
        raise StageError(
            stage,
            f"CRISPR-evOr {location}.prediction disagrees with its likelihood difference and threshold",
        )
    return group


def _validate_bias_correction_policy(
    document: Mapping[str, object],
    *,
    expected_requested: bool,
    expected_rho_failure: str,
) -> Mapping[str, object]:
    stage = "crispr_evor"
    policy = document.get("bias_correction_policy")
    if not isinstance(policy, Mapping):
        raise StageError(stage, "CRISPR-evOr result metadata is missing bias_correction_policy")
    for field in ("alpha_requested", "rho_requested"):
        if policy.get(field) is not expected_requested:
            raise StageError(
                stage,
                f"CRISPR-evOr bias policy {field} disagrees with the requested correction policy",
            )
    if policy.get("rho_failure") != expected_rho_failure:
        raise StageError(
            stage,
            "CRISPR-evOr bias policy rho_failure disagrees with the requested failure policy",
        )
    if policy.get("fallback_scope") != "per_group_per_orientation":
        raise StageError(stage, "CRISPR-evOr bias policy has an invalid fallback scope")
    if not _nonempty_string(policy.get("orientation_likelihood_effect")):
        raise StageError(stage, "CRISPR-evOr bias policy is missing its orientation likelihood effect")
    return policy


def _validate_result_bias_corrections(
    entry: Mapping[str, object],
    *,
    index: int,
    requested: bool,
    rho_failure_policy: str,
) -> list[dict[str, str]]:
    stage = "crispr_evor"
    location = f"results[{index}].bias_corrections"
    corrections = entry.get("bias_corrections")
    if not isinstance(corrections, Mapping):
        raise StageError(stage, f"CRISPR-evOr {location} must be a JSON object")

    warnings: list[dict[str, str]] = []
    for orientation in ("forward", "reverse"):
        effective = corrections.get(orientation)
        if not isinstance(effective, Mapping):
            raise StageError(stage, f"CRISPR-evOr {location}.{orientation} must be a JSON object")
        alpha_status = effective.get("alpha")
        rho_status = effective.get("rho")
        allowed_alpha = {"applied", "not_applicable"} if requested else {"disabled", "not_applicable"}
        allowed_rho = (
            {"applied", "fallback_uncorrected", "not_applicable"}
            if requested
            else {"disabled", "not_applicable"}
        )
        if alpha_status not in allowed_alpha:
            raise StageError(stage, f"CRISPR-evOr {location}.{orientation}.alpha is invalid")
        if rho_status not in allowed_rho:
            raise StageError(stage, f"CRISPR-evOr {location}.{orientation}.rho is invalid")

        warning_code = effective.get("warning_code")
        warning = effective.get("warning")
        if rho_status == "fallback_uncorrected":
            if rho_failure_policy != "fallback_uncorrected":
                raise StageError(stage, "CRISPR-evOr reported a rho fallback under a strict policy")
            if warning_code != "rho_bias_fit_unavailable_used_uncorrected" or not _nonempty_string(warning):
                raise StageError(stage, f"CRISPR-evOr {location}.{orientation} fallback lacks warning provenance")
            warnings.append({
                "group": str(entry["group"]),
                "orientation": orientation,
                "code": str(warning_code),
                "message": str(warning),
            })
        elif warning_code is not None or warning is not None:
            raise StageError(stage, f"CRISPR-evOr {location}.{orientation} has a spurious fallback warning")
    return warnings


def _validate_orientation_trees(
    output: Path,
    reference: object,
    group_leaves: Mapping[str, set[str]],
    predictions: Mapping[str, str],
) -> tuple[Path, dict[str, _NewickSummary]]:
    stage = "crispr_evor"
    groups = set(group_leaves)
    if not _nonempty_string(reference):
        raise StageError(stage, "CRISPR-evOr orientation_trees must be a non-empty relative path")
    relative = Path(reference)
    if relative.is_absolute() or ".." in relative.parts:
        raise StageError(stage, f"CRISPR-evOr orientation_trees path is unsafe: {reference!r}")

    candidate = output / relative
    try:
        resolved_output = output.resolve(strict=True)
        resolved_candidate = candidate.resolve(strict=True)
    except OSError as error:
        raise StageError(stage, f"cannot resolve CRISPR-evOr orientation tree output {candidate}: {error}") from error
    if resolved_output not in resolved_candidate.parents:
        raise StageError(stage, f"CRISPR-evOr orientation_trees path escapes its output directory: {reference!r}")

    tree_document = _read_json_object(stage, resolved_candidate, "CRISPR-evOr orientation tree output")
    tree_summaries: dict[str, dict[str, _NewickSummary]] = {}
    for policy in ("forward", "reverse", "selected"):
        tree_map = tree_document.get(policy)
        if not isinstance(tree_map, Mapping):
            raise StageError(stage, f"CRISPR-evOr orientation tree output is missing the {policy!r} tree map")
        invalid = [
            name
            for name, tree in tree_map.items()
            if not _nonempty_string(name) or not _nonempty_string(tree)
        ]
        if invalid:
            raise StageError(stage, f"CRISPR-evOr {policy!r} tree map contains invalid entries: {invalid}")
        missing = groups.difference(tree_map)
        if missing:
            raise StageError(stage, f"CRISPR-evOr {policy!r} tree map is missing groups: {sorted(missing)}")
        unexpected = set(tree_map).difference(groups)
        if unexpected:
            raise StageError(stage, f"CRISPR-evOr {policy!r} tree map has unexpected groups: {sorted(unexpected)}")
        policy_summaries: dict[str, _NewickSummary] = {}
        for group in sorted(groups):
            policy_summaries[group] = _validate_tree_leaf_identity(
                stage,
                tree_map[group],
                f"CRISPR-evOr {policy!r} tree for group {group!r}",
                group_leaves[group],
            )
        tree_summaries[policy] = policy_summaries
    for group in sorted(groups):
        selected_policy = "reverse" if predictions[group] == "Reverse" else "forward"
        if tree_document["selected"][group] != tree_document[selected_policy][group]:
            raise StageError(
                stage,
                f"CRISPR-evOr selected tree for group {group!r} disagrees with "
                f"prediction {predictions[group]!r}",
            )
    return resolved_candidate, tree_summaries["selected"]


def _validate_spacerplacer_outputs(
    output: Path,
    *,
    expected_group_leaves: Mapping[str, set[str]],
    expected_bias_corrections: bool,
    expected_rho_failure_policy: str,
) -> dict[str, object]:
    stage = "spacerplacer"
    expected_groups = set(expected_group_leaves)
    results_path = output / "0_results.csv"
    rows = _read_result_csv(
        stage,
        results_path,
        "SpacerPlacer result CSV",
        required_fields={
            "name",
            "deletion_rate_bdm",
            "insertion_rate_bdm",
            "ln_lh_bdm",
            "nb of spacers in alignment",
            "nb of reconstructed insertions",
            "nb of reconstructed deletions",
            "Deletion model preferred by LRT",
            "tree length",
            "run_time",
            "rho_bias_failure_policy",
            "alpha_bias_correction_status",
            "rho_bias_correction_status",
            "rho_bias_correction_warning_code",
            "rho_bias_correction_warning",
        },
    )
    _require_exact_groups(stage, set(rows), expected_groups, "SpacerPlacer result CSV")
    bias_warnings: list[dict[str, str]] = []
    from .evolution_validation import validate_metrics
    for group, row in rows.items():
        try:
            validate_metrics(row)
        except ValueError as error:
            raise StageError(stage, str(error)) from error
        for field in (
            "deletion_rate_bdm",
            "insertion_rate_bdm",
            "tree length",
            "run_time",
        ):
            _csv_nonnegative_float(stage, row, group, field)
        _csv_finite_float(stage, row, group, "ln_lh_bdm")
        for field in (
            "nb of spacers in alignment",
            "nb of reconstructed insertions",
            "nb of reconstructed deletions",
        ):
            _csv_nonnegative_int(stage, row, group, field)
        preferred_model = (row.get("Deletion model preferred by LRT") or "").strip()
        if preferred_model not in {"BDM", "IDM"}:
            raise StageError(
                stage,
                f"result CSV field 'Deletion model preferred by LRT' is invalid for group {group!r}",
            )
        reported_policy = (row.get("rho_bias_failure_policy") or "").strip()
        if reported_policy != expected_rho_failure_policy:
            raise StageError(
                stage,
                f"SpacerPlacer rho failure policy for group {group!r} disagrees with the requested policy",
            )
        alpha_status = (row.get("alpha_bias_correction_status") or "").strip()
        expected_alpha_status = "applied" if expected_bias_corrections else "disabled"
        if alpha_status != expected_alpha_status:
            raise StageError(
                stage,
                f"SpacerPlacer alpha correction status for group {group!r} disagrees with the request",
            )
        rho_status = (row.get("rho_bias_correction_status") or "").strip()
        allowed_rho_statuses = (
            {"applied", "fallback_uncorrected"}
            if expected_bias_corrections
            else {"disabled"}
        )
        if rho_status not in allowed_rho_statuses:
            raise StageError(
                stage,
                f"SpacerPlacer rho correction status for group {group!r} disagrees with the request",
            )
        warning_code = (row.get("rho_bias_correction_warning_code") or "").strip()
        warning = (row.get("rho_bias_correction_warning") or "").strip()
        if rho_status == "fallback_uncorrected":
            if expected_rho_failure_policy != "fallback_uncorrected":
                raise StageError(stage, "SpacerPlacer reported a rho fallback under a strict policy")
            if warning_code != "rho_bias_fit_unavailable_used_uncorrected" or not warning:
                raise StageError(
                    stage,
                    f"SpacerPlacer rho fallback for group {group!r} lacks warning provenance",
                )
            bias_warnings.append({
                "group": group,
                "orientation": "input_order",
                "code": warning_code,
                "message": warning,
            })
        elif warning_code or warning:
            raise StageError(
                stage,
                f"SpacerPlacer group {group!r} has a spurious rho fallback warning",
            )
    tree_path = output / "additional_data" / "dict_nwk_trees.json"
    tree_map = _read_json_object(stage, tree_path, "SpacerPlacer tree output")
    invalid_names = [name for name in tree_map if not _nonempty_string(name)]
    if invalid_names:
        raise StageError(stage, f"SpacerPlacer tree output has invalid group names: {invalid_names}")
    _require_exact_groups(stage, set(tree_map), expected_groups, "SpacerPlacer tree output")
    for group in sorted(expected_groups):
        summary = _validate_tree_leaf_identity(
            stage,
            tree_map[group],
            f"SpacerPlacer tree for group {group!r}",
            expected_group_leaves[group],
        )
        _validate_reported_tree_length(
            stage,
            _csv_nonnegative_float(stage, rows[group], group, "tree length"),
            summary,
            group,
        )
    return {
        "bias_correction_policy": {
            "alpha_requested": expected_bias_corrections,
            "rho_requested": expected_bias_corrections,
            "rho_failure": expected_rho_failure_policy,
            "fallback_scope": "per_group",
            "orientation_likelihood_effect": (
                "none; standalone SpacerPlacer preserves the supplied input orientation"
            ),
        },
        "bias_fallback_count": len(bias_warnings),
        "bias_warnings": bias_warnings,
        "result_count": len(rows),
        "result_groups": sorted(rows),
        "results": str(results_path),
        "trees": str(tree_path),
    }


def _validate_crispr_evor_outputs(
    output: Path,
    *,
    expected_group_leaves: Mapping[str, set[str]],
    expected_tree_policy: str,
    expected_bias_corrections: bool,
    expected_rho_failure_policy: str,
    expected_orientation_boundary: float | None = None,
) -> dict[str, object]:
    stage = "crispr_evor"
    expected_groups = set(expected_group_leaves)
    results_path = output / "0_results.csv"
    metadata_path = output / "crispr_evor_results.json"
    csv_rows = _read_result_csv(
        stage,
        results_path,
        "CRISPR-evOr result CSV",
        required_fields={
            "name",
            "predicted orientation",
            "recommend reversing array",
            "confidence threshold",
            "ln_lh_bdm",
            "reversed_ln_lh_bdm",
            "ln_lh_bdm - reversed_ln_lh_bdm",
            "nb of reconstructed insertions",
            "Deletion model preferred by LRT",
            "tree length",
        },
    )
    document = _read_json_object(stage, metadata_path, "CRISPR-evOr result metadata")

    if document.get("tool") != "CRISPR-evOr":
        raise StageError(stage, "CRISPR-evOr result metadata has an invalid or missing tool identifier")
    if not _nonempty_string(document.get("version")):
        raise StageError(stage, "CRISPR-evOr result metadata has an invalid or missing version")

    results = document.get("results")
    if not isinstance(results, list):
        raise StageError(stage, "CRISPR-evOr result metadata field 'results' must be a JSON array")
    bias_policy = _validate_bias_correction_policy(
        document,
        expected_requested=expected_bias_corrections,
        expected_rho_failure=expected_rho_failure_policy,
    )
    groups: set[str] = set()
    result_entries: dict[str, Mapping[str, object]] = {}
    bias_warnings: list[dict[str, str]] = []
    for index, entry in enumerate(results):
        group = _validate_orientation_result(entry, index)
        if expected_orientation_boundary is not None and entry.get("confidence_threshold") != expected_orientation_boundary:
            raise StageError(stage, "CRISPR-evOr did not apply the requested orientation boundary")
        if group in groups:
            raise StageError(stage, f"CRISPR-evOr result metadata contains duplicate group {group!r}")
        groups.add(group)
        result_entries[group] = entry
        bias_warnings.extend(_validate_result_bias_corrections(
            entry,
            index=index,
            requested=expected_bias_corrections,
            rho_failure_policy=expected_rho_failure_policy,
        ))

    _require_exact_groups(stage, groups, expected_groups, "CRISPR-evOr result metadata")
    if set(csv_rows) != groups:
        missing = groups.difference(csv_rows)
        unexpected = set(csv_rows).difference(groups)
        raise StageError(
            stage,
            "CRISPR-evOr result CSV groups disagree with result metadata: "
            f"missing {sorted(missing)}; unexpected {sorted(unexpected)}",
        )
    for group, entry in result_entries.items():
        row = csv_rows[group]
        from .evolution_validation import validate_metrics
        try:
            validate_metrics(row)
        except ValueError as error:
            raise StageError(stage, str(error)) from error
        for field in ("deletion_rate_bdm", "insertion_rate_bdm"):
            _csv_nonnegative_float(stage, row, group, field)
        _csv_nonnegative_int(stage, row, group, "nb of reconstructed deletions")
        selected = entry["bias_corrections"]["reverse" if entry["prediction"] == "Reverse" else "forward"]
        for column, key in (("alpha_bias_correction_status", "alpha"), ("rho_bias_correction_status", "rho"),
                            ("rho_bias_correction_warning_code", "warning_code"), ("rho_bias_correction_warning", "warning")):
            if column not in row or (row[column].strip() or None) != selected.get(key):
                raise StageError(stage, "Selected CSV correction differs from hypothesis metadata")
        if row.get("rho_bias_failure_policy") != expected_rho_failure_policy:
            raise StageError(stage, "Selected CSV rho failure policy mismatch")
        if (row.get("predicted orientation") or "").strip() != entry["prediction"]:
            raise StageError(stage, f"CRISPR-evOr CSV prediction disagrees with metadata for group {group!r}")
        recommendation = (row.get("recommend reversing array") or "").strip().lower()
        expected_recommendation = {False: "false", True: "true", None: "nd"}[entry["recommended_reverse"]]
        if recommendation != expected_recommendation:
            raise StageError(stage, f"CRISPR-evOr CSV recommendation disagrees with metadata for group {group!r}")
        selected_metadata_field = (
            "reverse_ln_likelihood_bdm"
            if entry["prediction"] == "Reverse"
            else "forward_ln_likelihood_bdm"
        )
        opposite_metadata_field = (
            "forward_ln_likelihood_bdm"
            if entry["prediction"] == "Reverse"
            else "reverse_ln_likelihood_bdm"
        )
        comparisons = (
            ("confidence threshold", float(entry["confidence_threshold"])),
            ("ln_lh_bdm", float(entry[selected_metadata_field])),
            ("reversed_ln_lh_bdm", float(entry[opposite_metadata_field])),
            (
                "ln_lh_bdm - reversed_ln_lh_bdm",
                float(entry["forward_minus_reverse_ln_likelihood_bdm"]),
            ),
        )
        for csv_field, expected_value in comparisons:
            csv_value = _csv_finite_float(stage, row, group, csv_field)
            if not math.isclose(csv_value, expected_value, rel_tol=1e-9, abs_tol=1e-9):
                raise StageError(
                    stage,
                    f"CRISPR-evOr CSV field {csv_field!r} disagrees with metadata for group {group!r}",
                )
        _csv_nonnegative_int(stage, row, group, "nb of reconstructed insertions")
        preferred_model = (row.get("Deletion model preferred by LRT") or "").strip()
        if preferred_model not in {"BDM", "IDM"}:
            raise StageError(
                stage,
                f"result CSV field 'Deletion model preferred by LRT' is invalid for group {group!r}",
            )

    status = document.get("status")
    expected_status = "completed" if results else "completed_no_results"
    if status != expected_status:
        raise StageError(
            stage,
            f"CRISPR-evOr result status {status!r} does not match its results; expected {expected_status!r}",
        )

    semantics = document.get("orientation_semantics")
    if not isinstance(semantics, Mapping) or any(
        not _nonempty_string(semantics.get(prediction)) for prediction in ("Forward", "Reverse", "ND")
    ):
        raise StageError(stage, "CRISPR-evOr result metadata has incomplete orientation_semantics")

    tree_policy = document.get("tree_policy")
    if not isinstance(tree_policy, str) or tree_policy not in ("estimated_separately", "provided_shared"):
        raise StageError(stage, f"CRISPR-evOr result metadata has invalid tree_policy {tree_policy!r}")
    if tree_policy != expected_tree_policy:
        raise StageError(
            stage,
            f"CRISPR-evOr reported tree_policy {tree_policy!r}; expected {expected_tree_policy!r}",
        )
    orientation_trees_path, selected_tree_summaries = _validate_orientation_trees(
        output,
        document.get("orientation_trees"),
        expected_group_leaves,
        {group: str(entry["prediction"]) for group, entry in result_entries.items()},
    )
    for group in sorted(groups):
        _validate_reported_tree_length(
            stage,
            _csv_nonnegative_float(stage, csv_rows[group], group, "tree length"),
            selected_tree_summaries[group],
            group,
        )
    return {
        "bias_correction_policy": dict(bias_policy),
        "bias_fallback_count": len(bias_warnings),
        "bias_warnings": bias_warnings,
        "orientation_trees": str(orientation_trees_path),
        "result_count": len(results),
        "result_groups": sorted(groups),
        "result_json": str(metadata_path),
        "result_status": status,
        "results": str(results_path),
        "tree_policy": tree_policy,
    }


def _downstream_arguments(
    ccdb_path: Path,
    output_path: Path,
    *,
    spacer_distance: int,
    tree_path: Path | None = None,
    bias_corrections: bool = True,
    alignment_workers: int | None = None,
    aligned_groups_path: Path | None = None,
    rho_bias_failure_policy: str | None = None,
) -> list[str]:
    arguments = [
        str(ccdb_path),
        str(output_path),
        "--input_type", "ccdb",
        "--cluster_spacers",
        "--cluster_spacers_max_distance", str(spacer_distance),
        "--no_plot_reconstruction",
        "--save_reconstructed_events",
        "--no_plot_order_graph",
        "--verbosity", "0",
    ]
    if tree_path is not None:
        arguments.extend(("--tree_path", str(tree_path)))
    if not bias_corrections:
        arguments.extend(("--no_alpha_bias_correction", "--no_rho_bias_correction"))
    if rho_bias_failure_policy is not None:
        arguments.extend(("--rho_bias_failure_policy", rho_bias_failure_policy))
    if alignment_workers is not None:
        arguments.extend(("--alignment_workers", str(alignment_workers)))
    if aligned_groups_path is not None:
        arguments.extend(("--aligned_groups_path", str(aligned_groups_path)))
    return arguments


def run_pipeline(
    input_path: str | Path,
    output_path: str | Path,
    *,
    input_kind: str = "fasta",
    runner_config: RunnerConfig | None = None,
    tools: Iterable[str] = ("spacerplacer", "crispr_evor"),
    categories: Iterable[str] = ("Bona-fide", "Possible"),
    minimum_arrays: int = 2,
    spacer_distance: int = 1,
    cpu: int | None = None,
    identify_folder_jobs: int | None = None,
    timeout_seconds: float = 3600,
    bias_corrections: bool = True,
    rho_bias_failure_policy: str = "fallback_uncorrected",
) -> dict[str, object]:
    """Run identification/adaptation/reconstruction as one isolated job."""

    if input_kind not in {"fasta", "identify-results"}:
        raise ValueError("input_kind must be 'fasta' or 'identify-results'")
    if spacer_distance < 0:
        raise ValueError("spacer_distance must be non-negative")
    if cpu is not None and cpu < 1:
        raise ValueError("cpu must be a positive integer")
    if identify_folder_jobs is not None and (
        isinstance(identify_folder_jobs, bool) or identify_folder_jobs < 1
    ):
        raise ValueError("identify_folder_jobs must be a positive integer")
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")
    if rho_bias_failure_policy not in {"fail", "fallback_uncorrected"}:
        raise ValueError("rho_bias_failure_policy must be 'fail' or 'fallback_uncorrected'")
    selected_tools = tuple(dict.fromkeys(tools))
    selected_categories = tuple(categories)
    unsupported = set(selected_tools) - {"spacerplacer", "crispr_evor"}
    if unsupported:
        raise ValueError(f"unsupported downstream tool(s): {sorted(unsupported)}")

    source = Path(input_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"pipeline input does not exist: {source}")
    job = _empty_job_root(output_path)
    logs = job / "logs"
    logs.mkdir()
    config = runner_config or default_runner_config()
    result_path = job / "integration_result.json"
    result: dict[str, object] = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "status": "running",
        "input": str(source),
        "input_kind": input_kind,
        "tool_options": config.tool_options,
        "stages": {},
    }
    _atomic_json(result_path, result)

    try:
        if input_kind == "fasta":
            identify_output = job / "crispridentify"
            input_flag = "--input_folder" if source.is_dir() else "--file"
            command = [
                *config.crispridentify,
                input_flag, str(source),
                "--result_folder", str(identify_output),
                "--strand", "True",
                "--html_report", "False",
                "--fasta_report", "False",
                "--cas", "False",
                "--accepted_only", (
                    "False"
                    if any(normalize_category(category) == "Low score" for category in selected_categories)
                    else "True"
                ),
            ]
            if cpu is not None:
                command.extend(("--cpu", str(cpu)))
            if identify_folder_jobs is not None:
                # CRISPRidentify bounds the effective lane count by the number
                # of input records and CPUs in its global CPU budget. This is
                # an operator-selected concurrency shape, not an extra budget.
                command.extend(("--folder_jobs", str(identify_folder_jobs)))
            command.extend(tool_arguments(config, "crispridentify"))
            result["stages"]["crispridentify"] = _run_process(
                "crispridentify", command, logs, timeout_seconds=timeout_seconds
            )
            reports = identify_output
        else:
            reports = source
            result["stages"]["crispridentify"] = {
                "status": "reused_existing_results",
                "path": str(reports),
            }

        adapter_output = job / "adapter"
        adapter_started = time.monotonic()
        adapter_manifest = adapt_reports(
            reports,
            adapter_output,
            categories=selected_categories,
            minimum_arrays=minimum_arrays,
            expected_fasta_path=source if input_kind == "fasta" else None,
        )
        result["stages"]["adapter"] = {
            "status": "completed",
            "duration_seconds": round(time.monotonic() - adapter_started, 6),
            "manifest": str(adapter_output / "manifest.json"),
            "emitted_array_count": adapter_manifest["emitted_array_count"],
            "emitted_group_count": adapter_manifest["emitted_group_count"],
            "source_binding": adapter_manifest["source_binding"],
        }
        expected_group_leaves = _adapter_group_leaf_names(adapter_manifest)

        if adapter_manifest["emitted_group_count"] == 0:
            if not selected_tools:
                result["status"] = "completed"
                _atomic_json(result_path, result)
                return result
            for tool in selected_tools:
                result["stages"][tool] = {
                    "status": "skipped",
                    "reason": "adapter produced no group with enough related arrays",
                }
            result["status"] = "completed_no_eligible_groups"
            _atomic_json(result_path, result)
            return result

        ccdb_path = adapter_output / "ccdb"
        # CRISPR-evOr already runs the complete SpacerPlacer reconstruction in
        # both directions and writes the selected-orientation reconstruction to
        # its own 0_results.csv. A standalone input-order result is misleading
        # whenever CRISPR-evOr selects Reverse, so retain that worker only for
        # the explicit SpacerPlacer-only workflow.
        run_standalone_spacerplacer = (
            "spacerplacer" in selected_tools and "crispr_evor" not in selected_tools
        )
        if run_standalone_spacerplacer:
            spacerplacer_output = job / "spacerplacer"
            command = [
                *config.spacerplacer,
                *_downstream_arguments(
                    ccdb_path,
                    spacerplacer_output,
                    spacer_distance=spacer_distance,
                    bias_corrections=bias_corrections,
                    alignment_workers=cpu,
                    rho_bias_failure_policy=rho_bias_failure_policy,
                ),
            ]
            command.extend(tool_arguments(config, "spacerplacer"))
            spacerplacer_stage = _run_process(
                "spacerplacer", command, logs, timeout_seconds=timeout_seconds
            )
            result["stages"]["spacerplacer"] = spacerplacer_stage
            try:
                spacerplacer_stage.update(
                    _validate_spacerplacer_outputs(
                        spacerplacer_output,
                        expected_group_leaves=expected_group_leaves,
                        expected_bias_corrections=bias_corrections,
                        expected_rho_failure_policy=rho_bias_failure_policy,
                    )
                )
            except StageError:
                spacerplacer_stage["status"] = "failed_output_validation"
                raise

        if "crispr_evor" in selected_tools:
            crispr_evor_output = job / "crispr_evor"
            command = [
                *config.crispr_evor,
                *_downstream_arguments(
                    ccdb_path,
                    crispr_evor_output,
                    spacer_distance=spacer_distance,
                    bias_corrections=bias_corrections,
                    alignment_workers=cpu,
                    rho_bias_failure_policy=rho_bias_failure_policy,
                ),
            ]
            command.extend(tool_arguments(config, "crispr_evor"))
            crispr_evor_stage = _run_process(
                "crispr_evor", command, logs, timeout_seconds=timeout_seconds
            )
            result["stages"]["crispr_evor"] = crispr_evor_stage
            try:
                crispr_evor_stage.update(
                    _validate_crispr_evor_outputs(
                        crispr_evor_output,
                        expected_group_leaves=expected_group_leaves,
                        expected_tree_policy="estimated_separately",
                        expected_orientation_boundary=config.tool_options.get("crispr_evor", {}).get("--orientation_decision_boundary"),
                        expected_bias_corrections=bias_corrections,
                        expected_rho_failure_policy=rho_bias_failure_policy,
                    )
                )
            except StageError:
                crispr_evor_stage["status"] = "failed_output_validation"
                raise
            crispr_evor_stage["tree_source"] = "crispr_evor_orientation_specific_estimates"
            if "spacerplacer" in selected_tools:
                result["stages"]["spacerplacer"] = {
                    "status": "fulfilled_by_crispr_evor",
                    "reason": (
                        "CRISPR-evOr includes SpacerPlacer and reports the "
                        "reconstruction selected after orientation comparison"
                    ),
                    "results": str(crispr_evor_output / "0_results.csv"),
                }

        result["status"] = "completed"
        _atomic_json(result_path, result)
        return result
    except Exception as error:
        result["status"] = "failed"
        result["error"] = {
            "type": type(error).__name__,
            "message": str(error),
            "stage": error.stage if isinstance(error, StageError) else "adapter_or_validation",
        }
        _atomic_json(result_path, result)
        raise
