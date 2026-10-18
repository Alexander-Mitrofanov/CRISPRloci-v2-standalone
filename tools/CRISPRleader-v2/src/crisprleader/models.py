"""Frozen shared records. Coordinates are zero-based half-open source-forward."""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from functools import cached_property
import hashlib
from typing import Any, Literal

Strand = Literal["+", "-", "?"]


@dataclass(frozen=True)
class Segment:
    start: int
    end: int

    def __post_init__(self) -> None:
        if type(self.start) is not int or type(self.end) is not int:
            raise ValueError("Segment endpoints must be integers")
        if not 0 <= self.start < self.end:
            raise ValueError("Segment must have 0 <= start < end")


@dataclass(frozen=True)
class SourceRecord:
    id: str
    sequence: str
    topology: Literal["linear", "circular", "unknown"] = "unknown"
    topology_provenance: str | None = None
    header: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @cached_property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.upper().encode("ascii")).hexdigest()

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class RepeatAnchor:
    segments: tuple[Segment, ...]
    sequence: str  # Source-forward sequence, always checked against original genome.


@dataclass(frozen=True)
class ArrayAnnotation:
    id: str
    source_id: str
    segments: tuple[Segment, ...]
    terminal_repeats: tuple[RepeatAnchor, RepeatAnchor]  # First/last on forward arc.
    consensus_repeat: str
    strand: Strand = "?"
    strand_constraint: bool = False
    category: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class FlankCandidate:
    id: str
    source_id: str
    array_id: str
    side: Literal["left", "right"]
    strand: Literal["+", "-"]
    segments: tuple[Segment, ...]  # Source-forward arc. Empty at a clipped edge.
    sequence: str  # Oriented leader-to-array, excludes terminal repeat.
    requested_length: int
    truncated_reason: str | None
    anchor: RepeatAnchor
    observed_repeat: str  # Terminal repeat oriented into this hypothesis.
    consensus_repeat: str  # Consensus oriented into this hypothesis.
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def available_length(self) -> int:
        return len(self.sequence)

    @property
    def junction_sequence(self) -> str:
        return self.sequence + self.observed_repeat


@dataclass(frozen=True)
class RawEvidence:
    candidate_id: str
    family_id: str
    method: str
    score: float | None
    score_units: str
    repeat_identity: float | None
    coverage: float | None
    start: int | None  # Oriented candidate-local core interval, excludes repeat.
    end: int | None
    passed: bool
    reasons: tuple[str, ...] = ()
    details: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class MotifHit:
    candidate_id: str
    model_id: str
    domain: Literal["leader", "junction", "repeat"]
    start: int  # Oriented junction_sequence offsets; can include terminal repeat.
    end: int
    strand: Literal["+", "-"]
    score: float
    pvalue: float | None = None
    qvalue: float | None = None
    details: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class LeaderDecision:
    source_id: str
    array_id: str
    status: Literal["supported", "ambiguous", "no_supported_hit", "unsupported"]
    selected_candidate_id: str | None
    core_start: int | None
    core_end: int | None
    reasons: tuple[str, ...]
    evidence: tuple[RawEvidence, ...]
    motifs: tuple[MotifHit, ...]
    details: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class AnalysisConfig:
    flank_length: int = 500
    max_flank_length: int = 5000
    max_arrays: int = 10000
    max_bases: int = 50000000
    max_candidates: int = 20000
    threads: int = 1
    timeout_seconds: float = 120.0
    policy_version: str = "leader-policy-v1"


def json_record(value: Any) -> Any:
    """Copy record containers in one traversal, preserving missing values and scores."""
    if hasattr(value, "__dataclass_fields__"):
        return {item.name: json_record(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, dict):
        return {key: json_record(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_record(item) for item in value]
    return value
