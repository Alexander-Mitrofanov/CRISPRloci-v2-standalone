"""Immutable domain values; internal intervals are zero-based, half-open."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from typing import Any, Literal

Strand = Literal["+", "-"]


@dataclass(frozen=True, order=True)
class Interval:
    start: int
    end: int

    def __post_init__(self) -> None:
        if type(self.start) is not int or type(self.end) is not int:
            raise ValueError("Interval coordinates must be integers")
        if self.start < 0 or self.end <= self.start:
            raise ValueError(f"Invalid half-open interval [{self.start}, {self.end})")

    @classmethod
    def from_public(cls, start: int, end: int) -> Interval:
        if type(start) is not int or type(end) is not int:
            raise ValueError("Public coordinates must be integers")
        if start > end:
            raise ValueError("Public coordinates must be ordered; strand is separate")
        return cls(start - 1, end)

    @property
    def length(self) -> int:
        return self.end - self.start

    def overlap(self, other: Interval) -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))

    def distance(self, other: Interval) -> int:
        """Number of intervening bases (zero for touching or overlapping intervals)."""
        return max(0, other.start - self.end, self.start - other.end)

    def public(self) -> dict[str, int]:
        return {"start": self.start + 1, "end": self.end, "length": self.length}


@dataclass(frozen=True)
class Sequence:
    id: str
    description: str
    bases: str

    @property
    def sha256(self) -> str:
        return sha256(self.bases.encode("ascii")).hexdigest()


@dataclass(frozen=True)
class Array:
    id: str
    interval: Interval
    repeat: str
    strand: Strand | None
    category: str
    score: float


@dataclass(frozen=True)
class Hit:
    query: str
    interval: Interval
    strand: Strand
    score: float
    evalue: float
    model: str = ""

    def public(self) -> dict[str, Any]:
        return {
            "interval": self.interval.public(),
            "strand": self.strand,
            "score": self.score,
            "evalue": self.evalue,
            "model": self.model,
        }


@dataclass(frozen=True)
class AntiRepeat:
    array: Array
    repeat: str
    repeat_strand: Strand
    interval: Interval
    strand: Strand
    identity: float
    coverage: float
    complementary: bool


@dataclass(frozen=True)
class Candidate:
    id: str
    anti: AntiRepeat
    interval: Interval
    sequence: str
    anti_local: Interval

    def to_genome(self, local: Interval) -> Interval:
        if local.end > self.interval.length:
            raise ValueError("Local interval exceeds candidate length")
        if self.anti.strand == "+":
            return Interval(self.interval.start + local.start, self.interval.start + local.end)
        return Interval(self.interval.end - local.end, self.interval.end - local.start)
