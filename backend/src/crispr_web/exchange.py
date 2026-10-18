"""Small, typed sequence ports shared by publication and submission.

These are data contracts, never instructions to launch tools or fetch files.
Native scientific validation remains independent of Integration.
"""

from __future__ import annotations

import hashlib
import json
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from .fasta import FastaValidationError


INPUT_KINDS = {
    "locus_compare": "target",
    "prophage_context": "target",
    "insertion_context": "target",
    "defense_context": "target",
    "rna_family": "rna",
    "array_compare": "array",
    "detection": "target",
    "orientation": "target",
    "reconstruction": "target",
    "cas": "target",
    "tracrrna": "target",
    "loci": "target",
    "leader": "target",
    "repeat_context": "target",
    "protospacer": "target",
    "repeats": "repeat",
    "repeat_map": "repeat",
    "repeat_type": "repeat",
    "viral_search": "spacer",
    "spacer_association": "spacer",
}
MAX_EXCHANGE_BYTES = 20_000_000


class ExchangeModel(BaseModel):
    model_config = ConfigDict(extra="forbid", strict=True)


class SequenceSchema(ExchangeModel):
    name: Literal["crisprloci.sequence-set"]
    version: Literal["1.0.0"]


class SequenceOrigin(ExchangeModel):
    source_id: str = Field(min_length=1, max_length=240, pattern=r"^[^\x00-\x1f\x7f]+$")
    source_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    array_id: str = Field(min_length=1, max_length=240, pattern=r"^[^\x00-\x1f\x7f]+$")
    ordinal: int = Field(ge=1)
    start: int = Field(ge=1)
    end: int = Field(ge=1)


class SequenceRecord(ExchangeModel):
    id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_-]{0,63}$")
    sequence: str = Field(
        min_length=1, max_length=MAX_EXCHANGE_BYTES, pattern=r"^[ACGTURYSWKMBDHVN]+$"
    )
    origin: SequenceOrigin | None = None

    @model_validator(mode="after")
    def check_interval(self):
        if self.origin and self.origin.end - self.origin.start + 1 != len(
            self.sequence
        ):
            raise ValueError("origin interval must match the sequence length")
        return self


class SequenceSet(ExchangeModel):
    schema_: SequenceSchema = Field(alias="schema")
    kind: Literal["target", "repeat", "spacer", "rna"]
    molecule: Literal["DNA", "RNA"]
    complete: bool
    records: list[SequenceRecord] = Field(max_length=100_000)

    @model_validator(mode="after")
    def check_records(self):
        if self.molecule == "RNA" and self.kind not in {"repeat", "rna"}:
            raise ValueError("only repeat and RNA-family ports accept RNA")
        ids = [row.id for row in self.records]
        if len(ids) != len(set(ids)):
            raise ValueError("sequence identifiers must be unique")
        forbidden = "U" if self.molecule == "DNA" else "T"
        if any(forbidden in row.sequence for row in self.records):
            raise ValueError("sequence alphabet does not match the molecule")
        return self

    def document(self):
        return self.model_dump(by_alias=True, exclude_none=True)


def strict_json(text):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result

    def invalid_constant(_value):
        raise ValueError("non-finite JSON number")

    return json.loads(text, object_pairs_hook=unique, parse_constant=invalid_constant)


def sequence_input(text: str, *, mode: str, molecule: str, max_bytes: int):
    """Return FASTA and optional validated handoff metadata; never truncate."""
    if not text.lstrip("\ufeff \t\r\n").startswith("{"):
        return text, None
    if len(text.encode("utf-8")) > min(max_bytes, MAX_EXCHANGE_BYTES):
        raise FastaValidationError("Sequence JSON exceeds the input byte limit")
    try:
        packet = SequenceSet.model_validate(strict_json(text.lstrip("\ufeff")))
    except (ValueError, ValidationError, RecursionError):
        raise FastaValidationError(
            "Invalid sequence JSON; upload a CRISPRloci sequence-set v1 file"
        ) from None
    if not packet.complete:
        raise FastaValidationError(
            "Incomplete sequence sets cannot be used as module inputs"
        )
    if not packet.records:
        raise FastaValidationError("This sequence set has no sequences to analyze")
    if packet.kind != INPUT_KINDS[mode]:
        raise FastaValidationError(
            f"{mode} requires {INPUT_KINDS[mode]} sequences, not {packet.kind} sequences"
        )
    if packet.molecule != molecule:
        raise FastaValidationError(
            f"Select {packet.molecule} input to use this sequence set"
        )
    fasta = "".join(f">{row.id}\n{row.sequence}\n" for row in packet.records)
    # Preserve immediate biological origins without copying prior run logs or
    # allowing caller-supplied evidence to become a claim of server validation.
    metadata = {
        "schema": packet.schema_.model_dump(),
        "kind": packet.kind,
        "submitted_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "origins": {
            row.id: row.origin.model_dump() for row in packet.records if row.origin
        },
    }
    return fasta, metadata


def sequence_set(kind, molecule, records, *, complete=True):
    return SequenceSet.model_validate(
        {
            "schema": {"name": "crisprloci.sequence-set", "version": "1.0.0"},
            "kind": kind,
            "molecule": molecule,
            "complete": complete,
            "records": records,
        }
    ).document()


def association_groups(source_ids, grouping, handoff=None):
    """An explicit user grouping assertion; source annotations are not verification."""
    keys = list(source_ids)
    if grouping == "single_set":
        rows = [{"id":"set1", "members":keys}]
    elif grouping == "per_record":
        rows = [{"id":key, "members":[key]} for key in keys]
    elif grouping == "per_source":
        origins = (handoff or {}).get("origins", {})
        if set(origins) != set(keys):
            raise FastaValidationError("Grouping by source requires an origin for every spacer in a sequence-set")
        buckets = {}
        for key in keys:
            origin = origins[key]
            buckets.setdefault((origin["source_id"],origin["source_sha256"]),[]).append(key)
        rows = [{"id":f"source{i+1}","members":members} for i,members in enumerate(buckets.values())]
    else:
        raise FastaValidationError("Choose an explicit spacer association grouping")
    if not 1 <= len(rows) <= 100:
        raise FastaValidationError("Spacer association supports at most 100 groups")
    return {"schema_version":"spacepharer.groups/1.0.0", "groups":rows}
