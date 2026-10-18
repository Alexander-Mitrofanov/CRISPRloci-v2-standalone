"""Strict multi-record FASTA I/O and strand-safe coordinate operations."""

from __future__ import annotations

import gzip
from collections.abc import Iterable, Iterator
from pathlib import Path

from .errors import TracrError
from .models import Interval, Sequence, Strand

IUPAC = "ACGTRYSWKMBDHVN"
_COMPLEMENT = str.maketrans(IUPAC, "TGCAYRSWMKVHDBN")


def normalize(bases: str) -> str:
    bases = bases.upper().replace("U", "T")
    invalid = set(bases) - set(IUPAC)
    if invalid or not bases:
        raise TracrError(f"Empty sequence or unsupported nucleotide symbols: {sorted(invalid)}")
    return bases


def reverse_complement(bases: str) -> str:
    return normalize(bases).translate(_COMPLEMENT)[::-1]


def subsequence(bases: str, interval: Interval, strand: Strand = "+") -> str:
    if strand not in ("+", "-"):
        raise TracrError("Sequence strand must be '+' or '-'")
    if interval.end > len(bases):
        raise TracrError("Sequence interval exceeds contig length")
    part = bases[interval.start : interval.end]
    return part if strand == "+" else reverse_complement(part)


def read_fasta(path: Path) -> Iterator[Sequence]:
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    seen: set[str] = set()
    header: str | None = None
    chunks: list[str] = []

    def record() -> Sequence:
        assert header is not None
        identifier = header.split()[0]
        if identifier in seen:
            raise TracrError(f"Duplicate FASTA identifier {identifier!r} in {path}")
        seen.add(identifier)
        return Sequence(identifier, header, normalize("".join(chunks)))

    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield record()
                header, chunks = line[1:].strip(), []
                if not header:
                    raise TracrError(f"Empty FASTA header at {path}:{number}")
            elif header is None:
                raise TracrError(f"Sequence before FASTA header at {path}:{number}")
            else:
                chunks.append("".join(line.split()))
    if header is None:
        raise TracrError(f"No FASTA records in {path}")
    yield record()


def write_fasta(path: Path, records: Iterable[Sequence]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.id}\n")
            for start in range(0, len(record.bases), 80):
                handle.write(record.bases[start : start + 80] + "\n")


def input_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path.resolve()]
    if not path.is_dir():
        raise TracrError(f"Input does not exist: {path}")
    extensions = (".fa", ".fna", ".fasta", ".fas", ".fa.gz", ".fna.gz", ".fasta.gz")
    files = sorted(
        p.resolve() for p in path.iterdir() if p.is_file() and p.name.lower().endswith(extensions)
    )
    if not files:
        raise TracrError(f"No FASTA files in {path}")
    return files
