"""Validated settings and independently replaceable external commands."""

from __future__ import annotations

import math
import tomllib
from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

from .errors import TracrError

DATA = Path(__file__).parent / "data"


def _finite_number(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and (isinstance(value, int) or math.isfinite(value))
    )


def _nonempty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value)


def _string_array(value: object, *, allow_empty: bool = False) -> bool:
    return (
        isinstance(value, (tuple, list))
        and (allow_empty or bool(value))
        and all(_nonempty_string(part) for part in value)
    )


@dataclass(frozen=True)
class Settings:
    mode: str = "complete"
    model_type: str = "II"
    threads: int = 1
    timeout: float = 3600
    similarity: float = 0.7
    coverage: float = 0.6
    tail_length: int = 150
    upstream_length: int = 5
    merge_gap: int = 150
    type_v_motifs: bool = False
    search_backend: str = "fasta36"
    model_rescue: bool = True
    model_rescue_evalue: float = 1e-5
    mask_arrays: bool = True
    intarna_batch_size: int = 128
    short_anti_repeat_rescue: bool = False
    identify_command: tuple[str, ...] = ("crispridentify",)
    casandra_command: tuple[str, ...] = ("casandra",)
    identify_preset: str = "balanced"
    gene_mode: str = "auto"
    cmscan: str = "cmscan"
    fasta36: str = "fasta36"
    blastn: str = "blastn"
    intarna: str = "IntaRNA"
    erpin: str = str(DATA / "erpin/erpin")
    erpin_model: str = str(DATA / "erpin/rho-indep.epn")
    models: tuple[str, ...] = ()
    environment: dict[str, str] = field(default_factory=dict)
    weights: dict[str, float] = field(
        default_factory=lambda: {
            "array": 0.5,
            "identity": 0.5,
            "coverage": 0.5,
            "identity_coverage": 0.5,
            "interaction": 0.6,
            "model": 0.9,
            "terminator": 0.9,
            "orientation": 0.1,
            "anti_tail": 0.1,
            "tail_terminator": 0.1,
        }
    )

    def __post_init__(self) -> None:
        choices = {
            "mode": ("complete", "model"),
            "model_type": ("II", "V"),
            "search_backend": ("fasta36", "blastn"),
            "identify_preset": ("balanced", "region-screen", "full"),
            "gene_mode": ("auto", "single", "meta"),
        }
        for name, options in choices.items():
            if getattr(self, name) not in options:
                raise TracrError(f"{name} must be one of {options}")
        for name in ("threads", "tail_length", "intarna_batch_size"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise TracrError(f"{name} must be a positive integer")
        for name in ("upstream_length", "merge_gap"):
            value = getattr(self, name)
            if type(value) is not int or value < 0:
                raise TracrError(f"{name} must be a nonnegative integer")
        for name in ("similarity", "coverage"):
            value = getattr(self, name)
            if not _finite_number(value) or not 0 <= value <= 1:
                raise TracrError(f"{name} must be finite and between 0 and 1")
        for name in ("timeout", "model_rescue_evalue"):
            value = getattr(self, name)
            if not _finite_number(value) or value <= 0:
                raise TracrError(f"{name} must be finite and positive")
        if type(self.type_v_motifs) is not bool:
            raise TracrError("type_v_motifs must be a boolean")
        if type(self.model_rescue) is not bool:
            raise TracrError("model_rescue must be a boolean")
        if type(self.short_anti_repeat_rescue) is not bool:
            raise TracrError("short_anti_repeat_rescue must be a boolean")
        if self.short_anti_repeat_rescue and (self.mode != "complete" or self.model_type != "II"):
            raise TracrError("Short anti-repeat rescue requires Type II complete mode")
        if type(self.mask_arrays) is not bool:
            raise TracrError("mask_arrays must be a boolean")
        if self.type_v_motifs and self.model_type != "V":
            raise TracrError("Type V motif analysis requires model_type V")
        for name in ("cmscan", "fasta36", "blastn", "intarna", "erpin", "erpin_model"):
            if not _nonempty_string(getattr(self, name)):
                raise TracrError(f"{name} must be a nonempty string")
        if not _string_array(self.models, allow_empty=True):
            raise TracrError("models must be an array of paths")
        for name in ("identify_command", "casandra_command"):
            if not _string_array(getattr(self, name)):
                raise TracrError(f"{name} must be a nonempty array of nonempty strings")
        if not isinstance(self.environment, dict) or any(
            not isinstance(k, str) or not isinstance(v, str) for k, v in self.environment.items()
        ):
            raise TracrError("environment must map names to strings")
        expected = {
            "array",
            "identity",
            "coverage",
            "identity_coverage",
            "interaction",
            "model",
            "terminator",
            "orientation",
            "anti_tail",
            "tail_terminator",
        }
        if (
            not isinstance(self.weights, Mapping)
            or set(self.weights) != expected
            or any(not _finite_number(w) or w < 0 for w in self.weights.values())
        ):
            raise TracrError(
                f"weights must be finite, nonnegative, and use keys {sorted(expected)}"
            )

    def model_paths(self) -> tuple[Path, ...]:
        names = ["II.cm"] if self.model_type == "II" else [f"V-{i}.cm" for i in range(1, 5)]
        return (
            tuple(Path(p).resolve() for p in self.models)
            if self.models
            else tuple(DATA / "models" / name for name in names)
        )

    def public(self) -> dict[str, Any]:
        result = asdict(self)
        # Environment may include credentials: keep values out of reports and logs.
        result["environment"] = sorted(self.environment)
        result["models"] = [str(p) for p in self.model_paths()]
        return result


def load_settings(path: Path | None, overrides: dict[str, Any]) -> Settings:
    values: dict[str, Any] = {}
    if path:
        with path.open("rb") as handle:
            values = tomllib.load(handle)
        # File-valued settings resolve relative to the config, not the caller's CWD.
        for key in ("models",):
            if key in values:
                if not _string_array(values[key], allow_empty=True):
                    raise TracrError("models must be an array of paths")
                values[key] = tuple(str((path.parent / p).resolve()) for p in values[key])
        for key in ("erpin_model",):
            if key in values:
                if not _nonempty_string(values[key]):
                    raise TracrError(f"{key} must be a nonempty string")
                values[key] = str((path.parent / values[key]).resolve())
        for key in ("cmscan", "fasta36", "blastn", "intarna", "erpin"):
            if key not in values:
                continue
            if not _nonempty_string(values[key]):
                raise TracrError(f"{key} must be a nonempty string")
            if "/" in values[key]:
                # Resolving a venv's Python symlink would bypass its environment.
                values[key] = str((path.parent / values[key]).absolute())
        for key in ("identify_command", "casandra_command"):
            if key in values:
                if not _string_array(values[key]):
                    raise TracrError(f"{key} must be a nonempty array of nonempty strings")
                values[key] = tuple(
                    str((path.parent / part).absolute())
                    if isinstance(part, str)
                    and "/" in part
                    and not part.startswith("-")
                    and "=" not in part
                    else part
                    for part in values[key]
                )
    values.update({k: v for k, v in overrides.items() if v is not None})
    unknown = values.keys() - Settings.__dataclass_fields__.keys()
    if unknown:
        raise TracrError(f"Unknown configuration keys: {sorted(unknown)}")
    return Settings(**values)
