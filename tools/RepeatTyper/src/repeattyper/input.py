"""Bounded FASTA parsing; IDs identify occurrences, never unique sequences."""
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

from .errors import InputError

MAX_BYTES = 2 * 1024 * 1024
MAX_RECORDS = 10000
MAX_LENGTH = 200
DNA = frozenset("ACGT")
IUPAC = frozenset("ACGTRYSWKMBDHVN")


@dataclass(frozen=True, slots=True)
class Record:
    id: str
    sequence: str


def validate_records(records: list[Record]) -> None:
    if not isinstance(records, (list, tuple)):
        raise InputError("Repeats must be an ordered collection of Record values")
    if len(records) > MAX_RECORDS:
        raise InputError("At most 10000 repeats are supported")
    ids: set[str] = set()
    for record in records:
        if not isinstance(record, Record):
            raise InputError("Repeats must contain Record values")
        if (not isinstance(record.id, str)
                or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,200}", record.id) or record.id in ids):
            raise InputError("Repeat IDs must be unique, bounded ASCII identifiers")
        ids.add(record.id)
        if (not isinstance(record.sequence, str)
                or not 1 <= len(record.sequence) <= MAX_LENGTH or set(record.sequence) - IUPAC):
            raise InputError("Repeats must contain 1–200 uppercase IUPAC DNA bases")


def read_fasta(path: Path) -> tuple[list[Record], str]:
    if any(part.is_symlink() for part in (path, *path.parents)) or not path.is_file():
        raise InputError("Input must be a regular, non-symlink FASTA file")
    with path.open("rb") as stream:
        payload = stream.read(MAX_BYTES + 1)
    if len(payload) > MAX_BYTES:
        raise InputError("Input exceeds 2 MiB")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError as exc:
        raise InputError("Input must be ASCII FASTA") from exc
    records: list[Record] = []
    identifier: str | None = None
    fragments: list[str] = []
    length = 0
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                records.append(Record(identifier, "".join(fragments).upper()))
            if len(records) == MAX_RECORDS:
                raise InputError("At most 10000 repeats are supported")
            header = line[1:].strip()
            identifier = header.split(maxsplit=1)[0] if header else ""
            fragments, length = [], 0
        else:
            if identifier is None:
                raise InputError("Sequence before FASTA header")
            length += len(line)
            if length > MAX_LENGTH:
                raise InputError("Repeats must contain 1–200 uppercase IUPAC DNA bases")
            fragments.append(line)
    if identifier is not None:
        records.append(Record(identifier, "".join(fragments).upper()))
    validate_records(records)
    return records, hashlib.sha256(payload).hexdigest()
