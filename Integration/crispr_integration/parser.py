"""Version-tolerant readers for CRISPRidentify machine reports.

The current canonical ``report.json`` format is preferred because it contains
source-validated spacer sequences and explicit orientation metadata.  The
legacy per-sequence and aggregate JSON formats are supported as a bounded
fallback.  Human-readable CSV/HTML output is intentionally not parsed.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping, Sequence


DNA_ALPHABET = frozenset("ACGTURYSWKMBDHVN")
FASTA_DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")
CANONICAL_SCHEMA_NAME = "CRISPRidentify-report"
SUPPORTED_CANONICAL_SCHEMA_MAJORS = frozenset({1})
_SEMANTIC_VERSION = re.compile(r"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
DEFAULT_MAX_REPORT_BYTES = 256 * 1024 * 1024
DEFAULT_MAX_REPORTS = 10_000
DEFAULT_MAX_ARRAYS = 200_000
DEFAULT_MAX_SPACERS_PER_ARRAY = 10_000


class ReportFormatError(ValueError):
    """A report cannot be interpreted without guessing at scientific data."""


@dataclass(frozen=True, order=True)
class SourceBinding:
    """Identity of one normalized FASTA record used by CRISPRidentify."""

    source_id: str
    length: int
    sha256: str


@dataclass(frozen=True)
class ArrayRecord:
    """One CRISPR array represented in the source report's sequence frame."""

    report_path: Path
    source_id: str
    array_id: str
    category: str
    consensus_repeat: str
    spacers: tuple[str, ...]
    strand: str | None
    sequence_orientation: str
    omitted_empty_spacers: int = 0


def _read_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Read FASTA records with CRISPRidentify's uppercase-sequence semantics."""

    if path.is_symlink():
        raise ReportFormatError(f"symbolic-link FASTA inputs are not accepted: {path}")
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []

    def finish_record(line_number: int) -> None:
        if header is None:
            return
        if not sequence_parts:
            raise ReportFormatError(
                f"{path}: empty FASTA record {header!r} before line {line_number}"
            )
        records.append((header, "".join(sequence_parts).upper()))

    try:
        with path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    finish_record(line_number)
                    header = line[1:].strip()
                    sequence_parts = []
                    if not header:
                        raise ReportFormatError(
                            f"{path}: line {line_number} has an empty FASTA header"
                        )
                    continue
                if header is None:
                    raise ReportFormatError(
                        f"{path}: line {line_number} contains sequence before its FASTA header"
                    )
                invalid = sorted(set(line) - FASTA_DNA_ALPHABET)
                if invalid:
                    raise ReportFormatError(
                        f"{path}: line {line_number} contains non-IUPAC DNA symbols: {invalid}"
                    )
                sequence_parts.append(line)
    except (OSError, UnicodeError) as error:
        raise ReportFormatError(f"cannot read FASTA input {path}: {error}") from error

    finish_record(line_number + 1 if "line_number" in locals() else 1)
    if not records:
        raise ReportFormatError(f"{path}: contains no FASTA records")
    return records


def load_fasta_source_bindings(path: str | Path) -> tuple[SourceBinding, ...]:
    """Return source identities for a FASTA file or flat folder of FASTA files."""

    root = Path(path).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"FASTA input does not exist: {root}")
    if root.is_file():
        files = [root]
        require_single_record = False
    elif root.is_dir():
        files = sorted(candidate for candidate in root.iterdir() if candidate.is_file())
        if not files:
            raise ReportFormatError(f"{root}: contains no regular input files")
        require_single_record = True
    else:
        raise ReportFormatError(f"FASTA input is not a regular file or directory: {root}")

    bindings: list[SourceBinding] = []
    seen_ids: set[str] = set()
    for file_path in files:
        records = _read_fasta_records(file_path)
        if require_single_record and len(records) != 1:
            raise ReportFormatError(
                f"{file_path}: folder inputs must contain exactly one FASTA record per file"
            )
        for header, sequence in records:
            source_id = _required_text(header.split()[0], "FASTA source ID", str(file_path))
            if source_id in seen_ids:
                raise ReportFormatError(f"duplicate FASTA source ID: {source_id!r}")
            seen_ids.add(source_id)
            bindings.append(SourceBinding(
                source_id=source_id,
                length=len(sequence),
                sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            ))
    return tuple(bindings)


def _normalized_key(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower()).strip("_")


def normalize_category(value: object) -> str:
    normalized = _normalized_key(value)
    aliases = {
        "bona_fide": "Bona-fide",
        "bonafide": "Bona-fide",
        "high_confidence": "Bona-fide",
        "possible": "Possible",
        "medium_confidence": "Possible",
        "low_score": "Low score",
        "lowscore": "Low score",
        "low_confidence": "Low score",
    }
    return aliases.get(normalized, str(value).strip())


def _get(mapping: Mapping[str, Any], *aliases: str, default: Any = None) -> Any:
    by_normalized_key = {_normalized_key(key): value for key, value in mapping.items()}
    for alias in aliases:
        key = _normalized_key(alias)
        if key in by_normalized_key:
            return by_normalized_key[key]
    return default


def _required_text(value: object, label: str, context: str) -> str:
    if value is None:
        raise ReportFormatError(f"{context}: missing {label}")
    text = str(value).strip()
    if not text:
        raise ReportFormatError(f"{context}: empty {label}")
    if any(character in text for character in "\r\n\0"):
        raise ReportFormatError(f"{context}: {label} contains a control character")
    return text


def _dna(value: object, label: str, context: str, *, allow_empty: bool = False) -> str:
    if not isinstance(value, str):
        raise ReportFormatError(f"{context}: {label} must be a DNA string")
    sequence = re.sub(r"\s+", "", value).upper()
    if not sequence and not allow_empty:
        raise ReportFormatError(f"{context}: {label} is empty")
    invalid = set(sequence) - DNA_ALPHABET
    if invalid:
        raise ReportFormatError(
            f"{context}: {label} contains invalid DNA character(s) {sorted(invalid)}"
        )
    return sequence


def _sequence_value(value: object, label: str, context: str) -> tuple[str, str]:
    """Return sequence plus ``source`` or ``biological`` frame."""

    if isinstance(value, str):
        return _dna(value, label, context, allow_empty=True), "unspecified"
    if not isinstance(value, Mapping):
        raise ReportFormatError(f"{context}: {label} must be a string or object")

    for key in ("source_sequence", "sequence", "dna", "value"):
        candidate = _get(value, key)
        if candidate is not None:
            frame = "source" if _normalized_key(key) == "source_sequence" else "unspecified"
            return _dna(candidate, label, context, allow_empty=True), frame
    candidate = _get(value, "oriented_sequence", "biological_sequence")
    if candidate is not None:
        return _dna(candidate, label, context, allow_empty=True), "biological"
    raise ReportFormatError(f"{context}: {label} object has no recognized sequence field")


def _strand(value: object) -> str | None:
    if isinstance(value, Mapping):
        nested = _get(value, "strand", "label", "direction", "orientation")
        status = _normalized_key(_get(value, "status", default=""))
        if nested is None and status in {"not_computed", "unknown", "nd"}:
            return None
        value = nested
    if value is None:
        return None
    raw = str(value).strip()
    normalized = _normalized_key(value)
    if raw in {"+", "+1", "1"} or normalized in {
        "forward", "positive", "pos", "plus",
    }:
        return "+"
    if raw in {"-", "-1"} or normalized in {
        "reverse", "reversed", "negative", "neg", "minus",
    }:
        return "-"
    if normalized in {
        "", "none", "unknown", "nd", "not_computed",
        "orientation_was_not_computed", "forward_orientation_was_not_computed",
    }:
        return None
    # A changed producer label is preserved as unknown rather than guessed.
    return None


def _array_sequence_frame(
    report: Mapping[str, Any], observed_frames: set[str]
) -> str:
    if "biological" in observed_frames:
        if observed_frames - {"biological"}:
            raise ReportFormatError("an array mixes source and biologically oriented spacer fields")
        return "biological"

    coordinate_system = _get(report, "coordinate_system", default={})
    if isinstance(coordinate_system, Mapping):
        reference = _normalized_key(_get(coordinate_system, "reference", default=""))
        separate = _get(coordinate_system, "orientation_is_separate")
        if reference == "input_forward_strand" or separate is True:
            return "source"
    # Current and legacy unlabelled ``sequence`` fields are source-oriented.
    return "source"


def _canonical_arrays(document: Mapping[str, Any]) -> Sequence[Any] | None:
    arrays = _get(document, "arrays", "crispr_arrays", "detected_arrays")
    if isinstance(arrays, list):
        return arrays
    results = _get(document, "results", "data")
    if isinstance(results, Mapping):
        arrays = _get(results, "arrays", "crispr_arrays", "detected_arrays")
        if isinstance(arrays, list):
            return arrays
    return None


def _canonical_schema_version(document: Mapping[str, Any], report_path: Path) -> str:
    schema = document.get("schema")
    if not isinstance(schema, Mapping):
        raise ReportFormatError(f"{report_path}: canonical report is missing a schema object")
    name = schema.get("name")
    if name != CANONICAL_SCHEMA_NAME:
        raise ReportFormatError(
            f"{report_path}: unsupported canonical schema name {name!r}; "
            f"expected {CANONICAL_SCHEMA_NAME!r}"
        )
    version = schema.get("version")
    if not isinstance(version, str) or _SEMANTIC_VERSION.fullmatch(version) is None:
        raise ReportFormatError(
            f"{report_path}: canonical schema version must be a three-part semantic version"
        )
    major = int(version.split(".", 1)[0])
    if major not in SUPPORTED_CANONICAL_SCHEMA_MAJORS:
        raise ReportFormatError(
            f"{report_path}: unsupported {CANONICAL_SCHEMA_NAME} schema major {major}; "
            f"supported major(s): {sorted(SUPPORTED_CANONICAL_SCHEMA_MAJORS)}"
        )
    return version


def _canonical_source_binding(
    document: Mapping[str, Any], report_path: Path
) -> SourceBinding:
    source = document.get("source")
    if not isinstance(source, Mapping):
        raise ReportFormatError(
            f"{report_path}: canonical report is missing its source provenance object"
        )
    source_id = _required_text(source.get("id"), "source.id", str(report_path))
    length = source.get("length")
    if isinstance(length, bool) or not isinstance(length, int) or length < 1:
        raise ReportFormatError(f"{report_path}: source.length must be a positive integer")
    sha256 = source.get("sha256")
    if not isinstance(sha256, str) or _SHA256.fullmatch(sha256) is None:
        raise ReportFormatError(
            f"{report_path}: source.sha256 must be a lowercase SHA-256 digest"
        )
    scope = source.get("sha256_scope")
    if scope != "uppercase_sequence":
        raise ReportFormatError(
            f"{report_path}: source.sha256_scope must be 'uppercase_sequence'"
        )
    return SourceBinding(source_id=source_id, length=length, sha256=sha256)


def _consume_expected_source(
    expected: Counter[SourceBinding], claim: SourceBinding, report_path: Path
) -> None:
    if expected[claim] > 0:
        expected[claim] -= 1
        return
    same_id = [binding for binding, count in expected.items() if count > 0 and binding.source_id == claim.source_id]
    if not same_id:
        raise ReportFormatError(
            f"{report_path}: report source ID {claim.source_id!r} does not match any submitted FASTA record"
        )
    expected_lengths = sorted({binding.length for binding in same_id})
    expected_hashes = sorted({binding.sha256 for binding in same_id})
    mismatches: list[str] = []
    if claim.length not in expected_lengths:
        mismatches.append(f"length {claim.length} != submitted {expected_lengths}")
    if claim.sha256 not in expected_hashes:
        mismatches.append("uppercase-sequence SHA-256 differs from submitted FASTA")
    raise ReportFormatError(
        f"{report_path}: report source provenance does not match submitted FASTA for "
        f"{claim.source_id!r}: {'; '.join(mismatches) or 'duplicate report'}"
    )


def _source_id(document: Mapping[str, Any], report_path: Path, fallback: str | None = None) -> str:
    source = _get(document, "source", "sequence", default={})
    if isinstance(source, Mapping):
        value = _get(source, "id", "accession", "name", "sequence_id")
        if value is not None:
            return _required_text(value, "source ID", str(report_path))
    value = _get(document, "source_id", "sequence_id", "accession")
    if value is not None:
        return _required_text(value, "source ID", str(report_path))
    inferred = report_path.parent.name if report_path.name == "report.json" else report_path.stem
    return _required_text(fallback or inferred, "source ID", str(report_path))


def _parse_canonical(
    document: Mapping[str, Any], report_path: Path, source_binding: SourceBinding
) -> list[ArrayRecord]:
    arrays = _canonical_arrays(document)
    if arrays is None:
        raise ReportFormatError(f"{report_path}: no recognized array list")
    source_id = source_binding.source_id
    parsed: list[ArrayRecord] = []
    for index, raw in enumerate(arrays, 1):
        context = f"{report_path}: array {index}"
        if not isinstance(raw, Mapping):
            raise ReportFormatError(f"{context} must be an object")
        array_id = _required_text(
            _get(raw, "id", "array_id", "region_id", "name", default=f"array-{index}"),
            "array ID",
            context,
        )
        category = normalize_category(
            _get(raw, "category", "category_name", "classification", "status", default="Unknown")
        )
        repeat_value = _get(raw, "consensus_repeat", "repeat_consensus", "consensus", "direct_repeat")
        repeat_container = _get(raw, "repeat")
        if repeat_value is None and isinstance(repeat_container, Mapping):
            repeat_value = _get(repeat_container, "consensus", "sequence")
        repeat = _dna(repeat_value, "consensus repeat", context)

        spacer_values = _get(raw, "spacers", "list_spacers", "spacer_sequences")
        if isinstance(spacer_values, Mapping):
            spacer_values = _get(spacer_values, "items", "sequences", "records")
        if not isinstance(spacer_values, list):
            raise ReportFormatError(f"{context}: spacers must be a list")
        if len(spacer_values) > DEFAULT_MAX_SPACERS_PER_ARRAY:
            raise ReportFormatError(f"{context}: too many spacers")

        spacers: list[str] = []
        omitted = 0
        observed_frames: set[str] = set()
        for spacer_index, value in enumerate(spacer_values, 1):
            sequence, frame = _sequence_value(value, f"spacer {spacer_index}", context)
            if frame != "unspecified":
                observed_frames.add(frame)
            if sequence:
                spacers.append(sequence)
            else:
                omitted += 1

        orientation = _get(raw, "orientation", "strand", "direction")
        parsed.append(ArrayRecord(
            report_path=report_path,
            source_id=source_id,
            array_id=array_id,
            category=category,
            consensus_repeat=repeat,
            spacers=tuple(spacers),
            strand=_strand(orientation),
            sequence_orientation=_array_sequence_frame(document, observed_frames),
            omitted_empty_spacers=omitted,
        ))
    return parsed


def _legacy_documents(document: Mapping[str, Any]) -> Iterable[tuple[str | None, Mapping[str, Any]]]:
    arrays = _get(document, "Arrays")
    if isinstance(arrays, Mapping):
        yield None, document
        return
    for source_id, child in document.items():
        if isinstance(child, Mapping) and isinstance(_get(child, "Arrays"), Mapping):
            yield str(source_id), child


def _parse_legacy(
    document: Mapping[str, Any], report_path: Path
) -> list[ArrayRecord]:
    parsed: list[ArrayRecord] = []
    found = False
    for fallback_source, child in _legacy_documents(document):
        found = True
        source_id = _source_id(child, report_path, fallback_source)
        categories = _get(child, "Arrays")
        if not isinstance(categories, Mapping):
            raise ReportFormatError(f"{report_path}: legacy Arrays value must be an object")
        running_index = 0
        for raw_category, values in categories.items():
            category = normalize_category(raw_category)
            if isinstance(values, Mapping):
                records = list(values.items())
            elif isinstance(values, list):
                records = list(enumerate(values, 1))
            else:
                raise ReportFormatError(
                    f"{report_path}: legacy category {raw_category!r} must contain objects"
                )
            for legacy_id, raw in records:
                running_index += 1
                context = f"{report_path}: legacy array {legacy_id}"
                if not isinstance(raw, Mapping):
                    raise ReportFormatError(f"{context} must be an object")
                spacer_values = _get(raw, "list_spacers", "spacers", "spacer_sequences")
                if not isinstance(spacer_values, list):
                    raise ReportFormatError(f"{context}: spacers must be a list")
                if len(spacer_values) > DEFAULT_MAX_SPACERS_PER_ARRAY:
                    raise ReportFormatError(f"{context}: too many spacers")
                spacers: list[str] = []
                omitted = 0
                for spacer_index, value in enumerate(spacer_values, 1):
                    sequence, _ = _sequence_value(value, f"spacer {spacer_index}", context)
                    if sequence:
                        spacers.append(sequence)
                    else:
                        omitted += 1
                repeat = _dna(
                    _get(raw, "consensus_repeat", "repeat_consensus", "consensus", "direct_repeat"),
                    "consensus repeat",
                    context,
                )
                array_id = _required_text(
                    _get(raw, "id", "array_id", "region_id", default=f"legacy-{running_index}-{legacy_id}"),
                    "array ID",
                    context,
                )
                parsed.append(ArrayRecord(
                    report_path=report_path,
                    source_id=source_id,
                    array_id=array_id,
                    category=category,
                    consensus_repeat=repeat,
                    spacers=tuple(spacers),
                    strand=_strand(_get(raw, "orientation", "strand", "direction")),
                    sequence_orientation="source",
                    omitted_empty_spacers=omitted,
                ))
    if not found:
        raise ReportFormatError(f"{report_path}: unrecognized CRISPRidentify JSON format")
    return parsed


def discover_report_files(path: str | Path, *, max_reports: int = DEFAULT_MAX_REPORTS) -> list[Path]:
    root = Path(path).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"CRISPRidentify result path does not exist: {root}")
    if root.is_file():
        reports = [root]
    elif root.is_dir():
        reports = sorted(candidate for candidate in root.rglob("report.json") if candidate.is_file())
        if not reports:
            aggregate = root / "Complete_json_summary.json"
            if aggregate.is_file():
                reports = [aggregate]
            else:
                legacy = sorted(
                    candidate for candidate in root.rglob("*.json")
                    if candidate.is_file() and candidate.parent.name.lower() in {"json", "json_report", "legacy_json"}
                )
                reports = legacy
    else:
        raise ReportFormatError(f"CRISPRidentify result path is not a regular file or directory: {root}")
    if not reports:
        raise ReportFormatError(
            f"{root}: no report.json or recognized legacy JSON output was found"
        )
    if len(reports) > max_reports:
        raise ReportFormatError(f"{root}: found {len(reports)} reports; limit is {max_reports}")
    for report in reports:
        if report.is_symlink():
            raise ReportFormatError(f"symbolic-link reports are not accepted: {report}")
    return reports


def load_identified_arrays(
    path: str | Path,
    *,
    expected_fasta_path: str | Path | None = None,
    max_report_bytes: int = DEFAULT_MAX_REPORT_BYTES,
    max_reports: int = DEFAULT_MAX_REPORTS,
    max_arrays: int = DEFAULT_MAX_ARRAYS,
) -> list[ArrayRecord]:
    """Load current or legacy CRISPRidentify reports with bounded resource use."""

    reports = discover_report_files(path, max_reports=max_reports)
    expected_sources = (
        Counter(load_fasta_source_bindings(expected_fasta_path))
        if expected_fasta_path is not None
        else None
    )
    arrays: list[ArrayRecord] = []
    for report in reports:
        size = report.stat().st_size
        if size > max_report_bytes:
            raise ReportFormatError(
                f"{report}: report is {size} bytes; limit is {max_report_bytes}"
            )
        try:
            with report.open("r", encoding="utf-8") as stream:
                document = json.load(stream)
        except (OSError, UnicodeError, json.JSONDecodeError) as error:
            raise ReportFormatError(f"cannot read JSON report {report}: {error}") from error
        if not isinstance(document, Mapping):
            raise ReportFormatError(f"{report}: top-level JSON value must be an object")

        canonical_arrays = _canonical_arrays(document)
        if "schema" in document or canonical_arrays is not None:
            _canonical_schema_version(document, report)
            if canonical_arrays is None:
                raise ReportFormatError(f"{report}: no recognized array list")
            source_binding = _canonical_source_binding(document, report)
            if expected_sources is not None:
                _consume_expected_source(expected_sources, source_binding, report)
            report_arrays = _parse_canonical(document, report, source_binding)
        else:
            if expected_sources is not None:
                raise ReportFormatError(
                    f"{report}: legacy reports cannot be bound to submitted FASTA provenance"
                )
            report_arrays = _parse_legacy(document, report)
        arrays.extend(report_arrays)
        if len(arrays) > max_arrays:
            raise ReportFormatError(
                f"identified array count exceeds the configured limit of {max_arrays}"
            )
    if expected_sources is not None:
        missing = sorted(
            binding.source_id
            for binding, count in expected_sources.items()
            if count > 0
        )
        if missing:
            raise ReportFormatError(
                "CRISPRidentify reports are missing submitted FASTA source(s): "
                + ", ".join(repr(source_id) for source_id in missing)
            )
    return arrays
