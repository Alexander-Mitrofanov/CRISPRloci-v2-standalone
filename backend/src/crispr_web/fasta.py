"""Bounded FASTA validation and deterministic per-record materialization."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
import os
from pathlib import Path
import re
import unicodedata


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_IUPAC_DNA_INPUT = IUPAC_DNA | frozenset(symbol.lower() for symbol in IUPAC_DNA)
_UNSAFE_NAME = re.compile(r"[^A-Za-z0-9_-]+")
_SEQUENCE_WHITESPACE = re.compile(r"\s+")
_LINE_BOUNDARIES = frozenset(
    {
        "\n",
        "\r",
        "\v",
        "\f",
        "\x1c",
        "\x1d",
        "\x1e",
        "\x85",
        "\u2028",
        "\u2029",
    }
)


class FastaValidationError(ValueError):
    """The submitted sequence is not bounded, normalized IUPAC DNA FASTA."""


@dataclass(frozen=True, slots=True)
class FastaRecord:
    identifier: str
    original_header: str
    sequence: str
    filename: str


@dataclass(frozen=True, slots=True)
class NormalizedFasta:
    records: tuple[FastaRecord, ...]
    total_bases: int
    combined: str

    @property
    def record_count(self) -> int:
        return len(self.records)

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.combined.encode("ascii")).hexdigest()

    @property
    def source_ids(self) -> tuple[str, ...]:
        return tuple(record.identifier for record in self.records)


def _safe_identifier(value: str, fallback: str, max_length: int = 64) -> str:
    # Canonical identifiers must survive queue verification and JSON handoffs.
    # In particular, a truncation can leave a trailing '_' or '-'.
    if len(value) <= max_length and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", value):
        return value
    ascii_value = (
        unicodedata.normalize("NFKD", value)
        .encode("ascii", errors="ignore")
        .decode("ascii")
    )
    cleaned = _UNSAFE_NAME.sub("_", ascii_value).strip("_-")[:max_length]
    return cleaned or fallback


def safe_upload_filename(value: str | None) -> str:
    """Reduce a browser-provided display name to one harmless basename."""

    raw = Path(value or "web_input.fasta").name
    stem = Path(raw).stem
    return f"{_safe_identifier(stem, 'web_input', 80)}.fasta"


def _iter_text_lines(value: str):
    """Yield ``str.splitlines()``-compatible lines without materializing a list."""

    start = 0
    offset = 0
    length = len(value)
    while offset < length:
        character = value[offset]
        if character not in _LINE_BOUNDARIES:
            offset += 1
            continue
        yield value[start:offset]
        if character == "\r" and offset + 1 < length and value[offset + 1] == "\n":
            offset += 1
        offset += 1
        start = offset
    if start < length:
        yield value[start:]


def _invalid_iupac_symbols(value: str, limit: int = 12, alphabet=IUPAC_DNA) -> tuple[str, ...]:
    """Return a bounded sample of literal non-ASCII-IUPAC input symbols."""

    invalid: list[str] = []
    seen: set[str] = set()
    for symbol in value:
        if symbol in alphabet or symbol in {base.lower() for base in alphabet} or symbol in seen:
            continue
        seen.add(symbol)
        invalid.append(symbol)
        if len(invalid) >= limit:
            break
    return tuple(invalid)


def normalize_fasta(
    text: str,
    *,
    max_total_bases: int,
    max_record_bases: int,
    max_records: int,
    max_header_characters: int,
    fallback_identifier: str = "web_input",
    molecule: str = "DNA",
) -> NormalizedFasta:
    """Validate FASTA and rewrite headers/line wrapping to a safe contract.

    Raw DNA text is accepted as a one-record convenience.  Submitted header
    text is retained only as metadata in memory; scientific workers receive a
    unique sanitized identifier, preventing headers from becoming paths.
    """

    if molecule not in {"DNA", "RNA"}:
        raise FastaValidationError("molecule must be DNA or RNA")
    alphabet = IUPAC_DNA if molecule == "DNA" else (IUPAC_DNA - {"T"}) | {"U"}
    if not isinstance(text, str):
        raise FastaValidationError("sequence must be text")
    document = text.lstrip("\ufeff").strip()
    if not document:
        raise FastaValidationError("enter or upload a DNA sequence")
    if not document.startswith(">"):
        document = f">{fallback_identifier}\n{document}"

    raw_records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_stream = io.StringIO()
    sequence_bases = 0
    accepted_bases = 0

    def finish_record() -> None:
        nonlocal sequence_stream, sequence_bases, accepted_bases
        if header is None:
            return
        if sequence_bases == 0:
            raise FastaValidationError(f"FASTA record {header!r} is empty")
        sequence = sequence_stream.getvalue()
        invalid = sorted(set(sequence).difference(alphabet))
        if invalid:
            symbols = ", ".join(repr(value) for value in invalid[:12])
            raise FastaValidationError(
                f"FASTA record {header!r} contains non-IUPAC DNA symbols: {symbols}"
            )
        if len(raw_records) >= max_records:
            raise FastaValidationError(f"input exceeds the {max_records}-record limit")
        raw_records.append((header, sequence))
        accepted_bases += sequence_bases
        sequence_stream = io.StringIO()
        sequence_bases = 0

    for line_number, raw_line in enumerate(_iter_text_lines(document), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish_record()
            header = line[1:].strip()
            if not header:
                raise FastaValidationError(
                    f"line {line_number} has an empty FASTA header"
                )
            if len(header) > max_header_characters:
                raise FastaValidationError(
                    f"line {line_number} FASTA header exceeds "
                    f"{max_header_characters} characters"
                )
            if any(
                ord(character) < 32 or ord(character) == 127 for character in header
            ):
                raise FastaValidationError(
                    f"line {line_number} FASTA header contains control characters"
                )
            continue
        if header is None:
            raise FastaValidationError(
                f"line {line_number} contains sequence before a FASTA header"
            )
        raw_sequence = _SEQUENCE_WHITESPACE.sub("", line)
        invalid = _invalid_iupac_symbols(raw_sequence, alphabet=alphabet)
        if invalid:
            symbols = ", ".join(repr(value) for value in invalid)
            raise FastaValidationError(
                f"FASTA record {header!r} contains non-IUPAC DNA symbols: {symbols}"
            )
        # Case conversion is deliberately performed only after literal ASCII
        # validation. Some Unicode characters (for example, long s) uppercase
        # to valid-looking IUPAC letters and must never become biological input.
        sequence = raw_sequence.upper()
        next_record_bases = sequence_bases + len(sequence)
        if next_record_bases > max_record_bases:
            raise FastaValidationError(
                f"FASTA record {header!r} exceeds the {max_record_bases:,}-base limit"
            )
        if accepted_bases + next_record_bases > max_total_bases:
            raise FastaValidationError(
                f"input exceeds the {max_total_bases:,}-base total limit"
            )
        sequence_stream.write(sequence)
        sequence_bases = next_record_bases
    finish_record()

    if not raw_records:
        raise FastaValidationError("input contains no FASTA records")
    total_bases = accepted_bases

    records: list[FastaRecord] = []
    original_ids: set[str] = set()
    sanitized_ids: set[str] = set()
    for index, (original_header, sequence) in enumerate(raw_records, 1):
        original_identifier = original_header.split()[0]
        if original_identifier in original_ids:
            raise FastaValidationError(
                f"duplicate FASTA record identifier: {original_identifier!r}"
            )
        original_ids.add(original_identifier)
        identifier = _safe_identifier(original_identifier, f"sequence_{index}")
        if identifier in sanitized_ids:
            raise FastaValidationError(
                "FASTA record identifiers collide after safe filename "
                f"normalization: {identifier!r}"
            )
        sanitized_ids.add(identifier)
        records.append(
            FastaRecord(
                identifier=identifier,
                original_header=original_header,
                sequence=sequence,
                filename=f"{index:03d}_{identifier}.fasta",
            )
        )

    lines: list[str] = []
    for record in records:
        lines.append(f">{record.identifier}")
        lines.extend(
            record.sequence[offset : offset + 80]
            for offset in range(0, len(record.sequence), 80)
        )
    return NormalizedFasta(
        records=tuple(records),
        total_bases=total_bases,
        combined="\n".join(lines) + "\n",
    )


def fsync_directory(path: Path) -> None:
    """Persist directory entries before the queue can reference them."""

    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_private(path: Path, value: str) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o400)
    try:
        with os.fdopen(descriptor, "w", encoding="ascii", newline="\n") as stream:
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
    except BaseException:
        path.unlink(missing_ok=True)
        raise


def materialize_fasta(job_root: Path, normalized: NormalizedFasta) -> None:
    """Write combined and per-record FASTA beneath a fresh private job root."""

    input_root = job_root / "input"
    records_root = input_root / "records"
    input_root.mkdir(mode=0o700)
    records_root.mkdir(mode=0o700)
    _write_private(input_root / "normalized.fasta", normalized.combined)
    for record in normalized.records:
        wrapped = "\n".join(
            record.sequence[offset : offset + 80]
            for offset in range(0, len(record.sequence), 80)
        )
        _write_private(
            records_root / record.filename,
            f">{record.identifier}\n{wrapped}\n",
        )
    os.chmod(records_root, 0o500)
    fsync_directory(records_root)
    os.chmod(input_root, 0o500)
    fsync_directory(input_root)
    fsync_directory(job_root)
    fsync_directory(job_root.parent)
