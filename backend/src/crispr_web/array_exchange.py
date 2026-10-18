"""Typed observed-array submission contract; uploaded origins remain annotations."""

from __future__ import annotations

import hashlib
from typing import Annotated, Literal
from pydantic import Field, model_validator
from .exchange import ExchangeModel, strict_json, MAX_EXCHANGE_BYTES
from .fasta import FastaValidationError

Identifier = Annotated[
    str, Field(min_length=1, max_length=240, pattern=r"^[^\x00-\x1f\x7f]+$")
]


class ArraySchema(ExchangeModel):
    name: Literal["crisprloci.arrays"]
    version: Literal["1.0.0"]


class Source(ExchangeModel):
    id: Identifier
    length: int = Field(ge=1, le=1_000_000_000)
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class Interval(ExchangeModel):
    start: int = Field(ge=1)
    end: int = Field(ge=1)

    @model_validator(mode="after")
    def ordered(self):
        if self.end < self.start:
            raise ValueError("reversed interval")
        return self


class Boundary(ExchangeModel):
    exact: bool | None = None
    occurrence_path: str | None = Field(default=None, max_length=80)
    repeat_phase: str | None = Field(default=None, max_length=80)


class Unit(ExchangeModel):
    id: Identifier | None = None
    ordinal: int = Field(ge=1, le=501)
    sequence: str | None = Field(
        default=None, min_length=1, max_length=1000, pattern=r"^[ACGTRYSWKMBDHVN]+$"
    )
    status: Literal["deletion", "not_reported", "not_source_verified"] | None = None
    start: int | None = Field(default=None, ge=1)
    end: int | None = Field(default=None, ge=1)
    boundary_position: int | None = Field(default=None, ge=0)


class Array(ExchangeModel):
    id: Identifier
    source_id: Identifier
    category: str = Field(min_length=1, max_length=80)
    model_score: float | None = Field(allow_inf_nan=False)
    strand: Literal["+", "-"] | None
    interval: Interval
    boundary: Boundary
    repeats: list[Unit] | None = Field(max_length=501)
    spacers: list[Unit] = Field(min_length=1, max_length=500)


class ArraySet(ExchangeModel):
    schema_: ArraySchema = Field(alias="schema")
    coordinate_system: Literal["1-based-inclusive-source-forward"]
    score_semantics: Literal["uncalibrated-model-score"]
    categories: list[str] = Field(min_length=1, max_length=10)
    sources: list[Source] = Field(min_length=1, max_length=200)
    arrays: list[Array] = Field(min_length=1, max_length=200)

    @model_validator(mode="after")
    def consistency(self):
        sources = {s.id: s for s in self.sources}
        if len(sources) != len(self.sources) or len(set(self.categories)) != len(
            self.categories
        ):
            raise ValueError("duplicate source or category")
        if sum(len(a.spacers) for a in self.arrays) > 10000:
            raise ValueError("too many spacer occurrences")
        seen = set()
        for array in self.arrays:
            scope = (array.source_id, array.id)
            if (
                scope in seen
                or array.source_id not in sources
                or array.category not in self.categories
            ):
                raise ValueError("invalid scoped array identity or category")
            seen.add(scope)
            if array.interval.end > sources[array.source_id].length:
                raise ValueError("array exceeds source")
            occupied = []
            for kind, units in [
                ("repeat", array.repeats or []),
                ("spacer", array.spacers),
            ]:
                ids = set()
                previous = array.interval.start - 1
                for ordinal, unit in enumerate(units, 1):
                    if unit.ordinal != ordinal or (
                        unit.id is not None and unit.id in ids
                    ):
                        raise ValueError("duplicate occurrence or invalid order")
                    if unit.id is not None:
                        ids.add(unit.id)
                    if unit.start is not None or unit.end is not None:
                        if (
                            unit.start is None
                            or unit.end is None
                            or not (
                                array.interval.start
                                <= unit.start
                                <= unit.end
                                <= array.interval.end
                            )
                        ):
                            raise ValueError("unit interval exceeds array bounds")
                    if unit.boundary_position is not None and unit.status != "deletion":
                        raise ValueError("only deletion slots carry boundary positions")
                    if unit.status == "deletion":
                        if (
                            unit.sequence is not None
                            or unit.start is not None
                            or unit.end is not None
                        ):
                            raise ValueError(
                                "deletion cannot contain DNA or an interval"
                            )
                        if (
                            unit.boundary_position is None
                            or not previous
                            <= unit.boundary_position
                            <= array.interval.end
                        ):
                            raise ValueError(
                                "deletion needs a consistent source boundary"
                            )
                        previous = unit.boundary_position
                    elif unit.status is not None:
                        if kind == "spacer":
                            raise ValueError(
                                "missing spacer evidence cannot be omitted"
                            )
                        if unit.sequence is not None:
                            raise ValueError("unverified repeat cannot carry sequence")
                    else:
                        if (
                            unit.sequence is None
                            or unit.id is None
                            or unit.start is None
                            or unit.end is None
                            or unit.boundary_position is not None
                        ):
                            raise ValueError(
                                "observed units require DNA, identity and source interval"
                            )
                        if (
                            unit.end - unit.start + 1 != len(unit.sequence)
                            or not array.interval.start
                            <= unit.start
                            <= unit.end
                            <= array.interval.end
                            or unit.start <= previous
                        ):
                            raise ValueError(
                                "inconsistent observed unit interval/order"
                            )
                        if kind == "spacer" and set(unit.sequence) - set("ACGT"):
                            raise ValueError(
                                "ambiguous spacers are unsupported for exact sharing"
                            )
                        occupied.append((unit.start, unit.end))
                        previous = unit.end
            occupied.sort()
            if any(
                first[1] >= second[0] for first, second in zip(occupied, occupied[1:])
            ):
                raise ValueError("observed intervals overlap")
        return self


def native_input(document):
    """Revalidate the map and return deterministic native arrays plus real DNA FASTA."""
    packet = ArraySet.model_validate(document)
    native = []
    records = {}
    for index, array in enumerate(packet.arrays):
        slots = []
        for unit in array.spacers:
            if unit.status == "deletion":
                ident = f"gap{unit.ordinal:06d}"
                while ident in {u.id for u in array.spacers if u.id is not None}:
                    ident += "_"
                slots.append(
                    {
                        "id": unit.id or ident,
                        "ordinal": unit.ordinal,
                        "sequence": None,
                        "status": "deletion",
                    }
                )
            else:
                ident = f"occ{index:06d}_{unit.ordinal:06d}"
                records[ident] = unit.sequence
                slots.append(
                    {
                        "id": unit.id,
                        "ordinal": unit.ordinal,
                        "sequence": unit.sequence,
                        "status": "observed",
                    }
                )
        native.append(
            {
                "id": f"array{index:06d}",
                "source_id": array.source_id,
                "array_id": array.id,
                "strand": array.strand,
                "spacers": slots,
            }
        )
    if not records:
        raise ValueError(
            "at least one observed spacer is required; deletion-only arrays may accompany it"
        )
    return {
        "schema": {"name": "cctk.arrays", "version": "1.0.0"},
        "arrays": native,
    }, records


def array_input(text, *, max_bytes):
    if len(text.encode()) > min(max_bytes, MAX_EXCHANGE_BYTES):
        raise FastaValidationError("Array JSON exceeds the input byte limit")
    try:
        document = strict_json(text.lstrip("\ufeff"))
        native, records = native_input(document)
    except (ValueError, TypeError, KeyError, RecursionError):
        raise FastaValidationError(
            "Supply consistent crisprloci.arrays v1 JSON with source-bound A/C/G/T spacers; missing spacer sequences are not accepted"
        ) from None
    return "".join(f">{key}\n{seq}\n" for key, seq in records.items()), {
        "schema": {"name": "crisprloci.array-map", "version": "1.0.0"},
        "submitted_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "document": document,
    }
