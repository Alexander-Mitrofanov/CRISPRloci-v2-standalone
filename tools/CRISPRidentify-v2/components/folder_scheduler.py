"""Pure helpers for scheduling independent input files across CPU lanes.

The module deliberately contains no process-launching or filesystem mutation.
Callers can therefore validate a folder and build child-process argument
vectors before starting any work.
"""

import math
import os
from pathlib import Path
from typing import Callable, Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union


CPURequest = Union[int, str, None]
FileWeight = Union[int, float]
WeightSource = Union[Mapping[Path, FileWeight], Callable[[Path], FileWeight]]

HIDDEN_LANE_FLAG = "--_folder_lane"
_LANE_VALUE_OPTIONS = (
    "--input_folder",
    "--file",
    "--input_folder_multifasta",
    "--result_folder",
    "--cpu",
    "--folder_jobs",
)


def _normalized_cpu_ids(cpu_ids: Iterable[int]) -> Tuple[int, ...]:
    """Return sorted, unique non-negative CPU identifiers."""
    normalized = []
    for cpu_id in cpu_ids:
        if isinstance(cpu_id, bool):
            raise ValueError("CPU identifiers must be integers")
        try:
            cpu_id = int(cpu_id)
        except (TypeError, ValueError):
            raise ValueError("CPU identifiers must be integers")
        if cpu_id < 0:
            raise ValueError("CPU identifiers must be non-negative")
        normalized.append(cpu_id)

    if not normalized:
        raise ValueError("at least one CPU must be available")
    if len(normalized) != len(set(normalized)):
        raise ValueError("CPU identifiers must be unique")
    return tuple(sorted(normalized))


def get_allowed_cpu_ids(
        affinity_getter: Optional[Callable[[int], Iterable[int]]] = None,
        cpu_count_getter: Optional[Callable[[], Optional[int]]] = None,
) -> Tuple[int, ...]:
    """Return CPUs available to this process, respecting scheduler affinity.

    ``os.sched_getaffinity(0)`` is preferred where available.  Platforms that
    do not expose affinity, or reject the query, fall back to
    ``os.cpu_count()``.  A conservative single CPU is used if the platform
    reports no count.

    The injectable getters keep the behavior deterministic in unit tests.
    """
    if affinity_getter is None:
        affinity_getter = getattr(os, "sched_getaffinity", None)
    if cpu_count_getter is None:
        cpu_count_getter = os.cpu_count

    if affinity_getter is not None:
        try:
            affinity = tuple(affinity_getter(0))
            if affinity:
                return _normalized_cpu_ids(affinity)
        except (AttributeError, NotImplementedError, OSError):
            pass

    cpu_count = cpu_count_getter()
    if not isinstance(cpu_count, int) or isinstance(cpu_count, bool) or cpu_count < 1:
        cpu_count = 1
    return tuple(range(cpu_count))


def resolve_cpu_budget(
        requested: CPURequest,
        allowed_cpus: Optional[Iterable[int]] = None,
        affinity_getter: Optional[Callable[[int], Iterable[int]]] = None,
        cpu_count_getter: Optional[Callable[[], Optional[int]]] = None,
) -> int:
    """Resolve ``requested`` to a positive budget clipped to allowed CPUs.

    ``None`` and the case-insensitive string ``"ALL"`` request every allowed
    CPU.  Positive integer strings are accepted for compatibility with the
    existing CLI.
    """
    if allowed_cpus is None:
        available = get_allowed_cpu_ids(affinity_getter, cpu_count_getter)
    else:
        available = _normalized_cpu_ids(allowed_cpus)

    if requested is None:
        return len(available)
    if isinstance(requested, str):
        normalized = requested.strip()
        if normalized.upper() == "ALL":
            return len(available)
        try:
            requested = int(normalized)
        except ValueError:
            raise ValueError("CPU budget must be ALL or a positive integer")
    if isinstance(requested, bool) or not isinstance(requested, int) or requested < 1:
        raise ValueError("CPU budget must be ALL or a positive integer")
    return min(requested, len(available))


def resolve_cpu_ids(
        requested: CPURequest,
        allowed_cpus: Optional[Iterable[int]] = None,
        affinity_getter: Optional[Callable[[int], Iterable[int]]] = None,
        cpu_count_getter: Optional[Callable[[], Optional[int]]] = None,
) -> Tuple[int, ...]:
    """Return the deterministic subset of allowed CPU IDs for a request."""
    if allowed_cpus is None:
        available = get_allowed_cpu_ids(affinity_getter, cpu_count_getter)
    else:
        available = _normalized_cpu_ids(allowed_cpus)
    budget = resolve_cpu_budget(requested, allowed_cpus=available)
    return available[:budget]


def partition_cpu_lanes(
        allowed_cpus: Iterable[int],
        lane_count: int,
) -> Tuple[Tuple[int, ...], ...]:
    """Split CPU IDs into balanced, disjoint, non-empty lanes.

    More requested lanes than CPUs are clipped to the number of CPUs.
    """
    cpus = _normalized_cpu_ids(allowed_cpus)
    if isinstance(lane_count, bool) or not isinstance(lane_count, int) or lane_count < 1:
        raise ValueError("lane_count must be a positive integer")

    lane_count = min(lane_count, len(cpus))
    small_size, larger_lane_count = divmod(len(cpus), lane_count)
    lanes = []
    start = 0
    for lane_index in range(lane_count):
        size = small_size + (1 if lane_index < larger_lane_count else 0)
        lanes.append(cpus[start:start + size])
        start += size
    return tuple(lanes)


def pipeline_result_key(file_path: Union[str, os.PathLike]) -> str:
    """Return the result key used by :class:`components.pipeline.Pipeline`."""
    name = Path(file_path).name
    if "." not in name:
        key = name
    else:
        # Only the final suffix is a file-format extension. Truncating at the
        # first dot aliases distinct inputs such as strain.v1.fasta and
        # strain.v2.fasta to the same output directory.
        key = name.rsplit(".", 1)[0]

    if not key:
        raise ValueError(
            "input file {!r} has an empty Pipeline result key".format(name)
        )
    if key in (".", ".."):
        raise ValueError(
            "input file {!r} has an unsafe Pipeline result key {!r}".format(
                name, key
            )
        )
    return key


def validate_unique_result_keys(
        files: Iterable[Union[str, os.PathLike]],
) -> Dict[str, Path]:
    """Index files by Pipeline result key, rejecting collisions.

    Validation is intended to run before lanes or child processes are
    launched, preventing two inputs from writing to the same result folder.
    """
    indexed = {}
    for file_path in files:
        path = Path(file_path)
        key = pipeline_result_key(path)
        if not key:
            raise ValueError(
                "input file {!r} has an empty Pipeline result key".format(str(path))
            )
        if key in indexed:
            raise ValueError(
                "input files {!r} and {!r} share Pipeline result key {!r}".format(
                    str(indexed[key]), str(path), key
                )
            )
        indexed[key] = path
    return indexed


def _validated_weight(value: FileWeight, path: Path) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("weight for {!r} must be numeric".format(str(path)))
    value = float(value)
    if not math.isfinite(value) or value < 0:
        raise ValueError(
            "weight for {!r} must be finite and non-negative".format(str(path))
        )
    return value




def rewrite_lane_argv(
        argv: Sequence[str],
        input_folder: Union[str, os.PathLike],
        result_folder: Union[str, os.PathLike],
        cpu: Union[int, str],
        folder_jobs: int = 1,
        add_hidden_lane_flag: bool = True,
) -> Tuple[str, ...]:
    """Return child argv with lane-specific options safely replaced.

    Tokens are handled as an argument vector; no quoting, shell parsing, or
    string command construction occurs.  Both ``--option value`` and
    ``--option=value`` forms are removed before canonical replacements are
    appended.  The hidden lane marker is a standalone flag.
    """
    if isinstance(folder_jobs, bool) or not isinstance(folder_jobs, int) or folder_jobs < 1:
        raise ValueError("folder_jobs must be a positive integer")
    if isinstance(cpu, bool):
        raise ValueError("cpu must be ALL or a positive integer")
    cpu_text = str(cpu)
    if cpu_text.upper() != "ALL":
        try:
            if int(cpu_text) < 1:
                raise ValueError
        except ValueError:
            raise ValueError("cpu must be ALL or a positive integer")

    tokens = list(argv)
    rewritten = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if not isinstance(token, str):
            raise TypeError("argv tokens must be strings")

        if token == HIDDEN_LANE_FLAG or token.startswith(HIDDEN_LANE_FLAG + "="):
            index += 1
            continue

        option = token.split("=", 1)[0]
        if option in _LANE_VALUE_OPTIONS:
            if "=" in token:
                index += 1
                continue
            if index + 1 >= len(tokens):
                raise ValueError("missing value for {}".format(option))
            index += 2
            continue

        rewritten.append(token)
        index += 1

    rewritten.extend((
        "--input_folder", os.fspath(input_folder),
        "--result_folder", os.fspath(result_folder),
        "--cpu", cpu_text,
        "--folder_jobs", str(folder_jobs),
    ))
    if add_hidden_lane_flag:
        rewritten.append(HIDDEN_LANE_FLAG)
    return tuple(rewritten)
