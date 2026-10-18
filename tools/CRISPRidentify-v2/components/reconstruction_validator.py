"""Independent validation of CRISPR candidate coordinates and sequences.

The detector's coordinates are deliberately treated as claims to verify, not
as a source of truth.  All coordinates accepted by this module are zero-based
repeat starts.  ``coordinate_offset`` is added to every supplied start so a
candidate found in a local sequence window can be checked against the original
whole-genome sequence.

This module does not import detection, enhancement, evaluation, or output
code.  It can therefore be used at any pipeline boundary and in benchmark
post-processing without reproducing their coordinate calculations.
"""

from dataclasses import dataclass, field
from functools import lru_cache
import os
from typing import Any, Dict, Optional, Tuple


COORDINATE_SYSTEM = "0-based, end-exclusive"
IUPAC_DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVNacgtryswkmbdhvn")


# Validation is deliberately called at several independent pipeline fences.
# Keep only the most recently seen immutable source so those fences do not
# repeatedly allocate and uppercase a complete multi-megabase genome.  The
# cache is identity based: equal strings from different callers are normalized
# independently, and str subclasses retain their potentially custom ``upper``
# semantics.
_NORMALIZED_GENOME_CACHE = (None, None)


def _normalized_genome_sequence(genome_sequence):
    global _NORMALIZED_GENOME_CACHE

    if type(genome_sequence) is not str:
        return genome_sequence.upper()
    cache = _NORMALIZED_GENOME_CACHE
    if genome_sequence is cache[0]:
        return cache[1]

    normalized = genome_sequence.upper()
    _NORMALIZED_GENOME_CACHE = (genome_sequence, normalized)
    return normalized


@dataclass(frozen=True)
class ReconstructionIssue:
    """One failed reconstruction invariant."""

    code: str
    message: str
    component: Optional[str] = None
    index: Optional[int] = None
    details: Dict[str, Any] = field(default_factory=dict)

    def as_dict(self):
        value = {
            "code": self.code,
            "message": self.message,
        }
        if self.component is not None:
            value["component"] = self.component
        if self.index is not None:
            value["index"] = self.index
        if self.details:
            value["details"] = dict(self.details)
        return value


@dataclass(frozen=True)
class ReconstructionResult:
    """Structured result returned by :func:`validate_reconstruction`."""

    genome_length: int
    coordinate_offset: int
    local_repeat_starts: Tuple[Any, ...]
    genomic_repeat_starts: Tuple[Optional[int], ...]
    repeat_count: int
    spacer_count: int
    array_start: Optional[int]
    array_end: Optional[int]
    issues: Tuple[ReconstructionIssue, ...]

    @property
    def valid(self):
        return not self.issues

    def as_dict(self):
        coordinates = {
            "system": COORDINATE_SYSTEM,
            "coordinate_offset": self.coordinate_offset,
            "local_repeat_starts": list(self.local_repeat_starts),
            "genomic_repeat_starts": list(self.genomic_repeat_starts),
            "array_start": self.array_start,
            "array_end": self.array_end,
        }
        if self.array_start is not None and self.array_end is not None:
            coordinates["array_start_1based"] = self.array_start + 1
            coordinates["array_end_1based_inclusive"] = self.array_end
        return {
            "valid": self.valid,
            "genome_length": self.genome_length,
            "repeat_count": self.repeat_count,
            "spacer_count": self.spacer_count,
            "coordinates": coordinates,
            "issues": [issue.as_dict() for issue in self.issues],
        }


class ReconstructionValidationError(ValueError):
    """Raised when an invalid candidate reaches a strict pipeline boundary."""

    def __init__(self, result):
        self.result = result
        codes = ", ".join(issue.code for issue in result.issues)
        super().__init__("invalid CRISPR reconstruction: {}".format(codes))


def _fasta_file_identity(path):
    stat_result = os.stat(path)
    return (
        stat_result.st_dev,
        stat_result.st_ino,
        stat_result.st_size,
        stat_result.st_mtime_ns,
        stat_result.st_ctime_ns,
    )


@lru_cache(maxsize=2)
def _read_single_fasta_cached(path, file_identity):
    """Read and validate one immutable file-identity snapshot."""

    header = None
    sequence_parts = []
    record_count = 0
    with open(path, encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record_count += 1
                if record_count > 1:
                    raise ValueError(
                        "{} contains more than one FASTA record".format(path)
                    )
                header = line[1:].strip()
                if not header:
                    raise ValueError("{} line {} has an empty FASTA header".format(
                        path, line_number
                    ))
                continue
            if header is None:
                raise ValueError(
                    "{} line {} contains sequence before its FASTA header".format(
                        path, line_number
                    )
                )
            symbols = set(line)
            invalid = symbols.difference(IUPAC_DNA_ALPHABET)
            if invalid and any(symbol.isspace() for symbol in invalid):
                raise ValueError("{} line {} contains embedded whitespace".format(
                    path, line_number
                ))
            if invalid:
                invalid = sorted(invalid)
                raise ValueError(
                    "{} line {} contains non-IUPAC DNA symbols: {}".format(
                        path,
                        line_number,
                        ", ".join(repr(symbol) for symbol in invalid),
                    )
                )
            sequence_parts.append(line)
    if record_count != 1 or not sequence_parts:
        raise ValueError("{} must contain one non-empty FASTA record".format(path))
    return header, "".join(sequence_parts).upper()


def read_single_fasta(path):
    """Read one non-empty FASTA record, reusing its validated snapshot.

    Several independent correctness fences consume the same genome during one
    pipeline run. The bounded cache keys the immutable parsed value by the
    resolved path and filesystem identity, so an in-place replacement or edit
    is revalidated while repeated reads avoid rescanning every nucleotide.
    """

    normalized_path = os.path.realpath(os.fspath(path))
    identity = _fasta_file_identity(normalized_path)
    result = _read_single_fasta_cached(normalized_path, identity)
    if _fasta_file_identity(normalized_path) != identity:
        return read_single_fasta(normalized_path)
    return result


def _as_tuple(value, field_name):
    if isinstance(value, (str, bytes)) or value is None:
        raise TypeError("{} must be a sequence, not a string".format(field_name))
    try:
        return tuple(value)
    except TypeError:
        raise TypeError("{} must be an iterable sequence".format(field_name))


def _is_integer(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _first_sequence_difference(expected, observed):
    limit = min(len(expected), len(observed))
    for index in range(limit):
        if expected[index] != observed[index]:
            return index
    if len(expected) != len(observed):
        return limit
    return None


def _sequence_mismatch_details(expected, observed):
    return {
        "expected_length": len(expected),
        "observed_length": len(observed),
        "first_difference": _first_sequence_difference(expected, observed),
        "expected": expected,
        "observed": observed,
    }


def _append_issue(
    issues,
    code,
    message,
    component=None,
    index=None,
    **details
):
    issues.append(
        ReconstructionIssue(
            code=code,
            message=message,
            component=component,
            index=index,
            details=details,
        )
    )


def validate_reconstruction(
    genome_sequence,
    repeat_starts,
    repeats,
    spacers,
    coordinate_offset=0,
):
    """Verify a CRISPR candidate by reconstructing it from the original DNA.

    Args:
        genome_sequence: Original complete DNA sequence.
        repeat_starts: Zero-based repeat starts.  They may be local to a
            search window.
        repeats: Ungapped repeat sequences in increasing genomic order.
        spacers: Spacer sequences in increasing genomic order.  There must be
            exactly one fewer spacer than repeat.
        coordinate_offset: Zero-based genomic start of the local search
            window.  Use zero when ``repeat_starts`` are already genomic.

    Returns:
        :class:`ReconstructionResult`.  ``array_end`` is zero-based and
        end-exclusive.  The corresponding 1-based inclusive end has the same
        numeric value and is included by ``as_dict``.
    """

    if not isinstance(genome_sequence, str):
        raise TypeError("genome_sequence must be a string")
    if not _is_integer(coordinate_offset) or coordinate_offset < 0:
        raise ValueError("coordinate_offset must be a non-negative integer")

    local_starts = _as_tuple(repeat_starts, "repeat_starts")
    repeat_sequences = _as_tuple(repeats, "repeats")
    spacer_sequences = _as_tuple(spacers, "spacers")
    genome = _normalized_genome_sequence(genome_sequence)
    issues = []

    if not repeat_sequences:
        _append_issue(
            issues,
            "no_repeats",
            "a CRISPR reconstruction must contain at least one repeat",
            component="array",
        )
    if len(local_starts) != len(repeat_sequences):
        _append_issue(
            issues,
            "repeat_start_count_mismatch",
            "repeat starts and repeat sequences have different counts",
            component="array",
            repeat_start_count=len(local_starts),
            repeat_count=len(repeat_sequences),
        )
    expected_spacer_count = max(0, len(repeat_sequences) - 1)
    if len(spacer_sequences) != expected_spacer_count:
        _append_issue(
            issues,
            "spacer_count_mismatch",
            "a candidate must contain exactly one spacer between each repeat",
            component="array",
            expected_spacer_count=expected_spacer_count,
            observed_spacer_count=len(spacer_sequences),
        )

    genomic_starts = []
    for index, local_start in enumerate(local_starts):
        if not _is_integer(local_start):
            genomic_starts.append(None)
            _append_issue(
                issues,
                "invalid_repeat_start",
                "repeat start must be an integer",
                component="repeat",
                index=index,
                observed=repr(local_start),
            )
            continue
        genomic_start = coordinate_offset + local_start
        genomic_starts.append(genomic_start)
        if local_start < 0:
            _append_issue(
                issues,
                "negative_local_repeat_start",
                "repeat start is negative in its local coordinate space",
                component="repeat",
                index=index,
                local_start=local_start,
                genomic_start=genomic_start,
            )

    for index in range(1, len(genomic_starts)):
        previous = genomic_starts[index - 1]
        current = genomic_starts[index]
        if previous is not None and current is not None and current <= previous:
            _append_issue(
                issues,
                "non_monotonic_repeat_starts",
                "repeat starts must be strictly increasing in genomic order",
                component="repeat",
                index=index,
                previous_genomic_start=previous,
                genomic_start=current,
            )

    checked_repeat_count = min(len(genomic_starts), len(repeat_sequences))
    for index in range(checked_repeat_count):
        repeat = repeat_sequences[index]
        genomic_start = genomic_starts[index]
        if not isinstance(repeat, str) or not repeat:
            _append_issue(
                issues,
                "invalid_repeat_sequence",
                "repeat sequence must be a non-empty string",
                component="repeat",
                index=index,
                observed=repr(repeat),
            )
            continue
        if genomic_start is None:
            continue
        genomic_end = genomic_start + len(repeat)
        if genomic_start < 0 or genomic_end > len(genome):
            _append_issue(
                issues,
                "repeat_out_of_bounds",
                "repeat interval lies outside the original genome",
                component="repeat",
                index=index,
                genomic_start=genomic_start,
                genomic_end=genomic_end,
                genome_length=len(genome),
            )
            continue
        observed = repeat.upper()
        expected = genome[genomic_start:genomic_end]
        if observed != expected:
            _append_issue(
                issues,
                "repeat_sequence_mismatch",
                "repeat sequence does not match the original genome",
                component="repeat",
                index=index,
                genomic_start=genomic_start,
                genomic_end=genomic_end,
                **_sequence_mismatch_details(expected, observed)
            )

    checked_spacer_count = min(
        len(spacer_sequences),
        max(0, checked_repeat_count - 1),
    )
    for index in range(checked_spacer_count):
        repeat = repeat_sequences[index]
        spacer = spacer_sequences[index]
        genomic_start = genomic_starts[index]
        next_genomic_start = genomic_starts[index + 1]
        if not isinstance(spacer, str):
            _append_issue(
                issues,
                "invalid_spacer_sequence",
                "spacer sequence must be a string",
                component="spacer",
                index=index,
                observed=repr(spacer),
            )
            continue
        if (
            not isinstance(repeat, str)
            or not repeat
            or genomic_start is None
            or next_genomic_start is None
        ):
            continue

        spacer_start = genomic_start + len(repeat)
        spacer_end = next_genomic_start
        expected_next_start = spacer_start + len(spacer)
        if expected_next_start != next_genomic_start:
            _append_issue(
                issues,
                "spacer_boundary_mismatch",
                "repeat start is inconsistent with the preceding repeat and spacer lengths",
                component="spacer",
                index=index,
                spacer_start=spacer_start,
                expected_next_repeat_start=expected_next_start,
                observed_next_repeat_start=next_genomic_start,
            )
        if spacer_end < spacer_start:
            _append_issue(
                issues,
                "overlapping_repeats",
                "the next repeat begins before the preceding repeat ends",
                component="spacer",
                index=index,
                spacer_start=spacer_start,
                spacer_end=spacer_end,
            )
            continue
        if spacer_start < 0 or spacer_end > len(genome):
            _append_issue(
                issues,
                "spacer_out_of_bounds",
                "spacer interval lies outside the original genome",
                component="spacer",
                index=index,
                genomic_start=spacer_start,
                genomic_end=spacer_end,
                genome_length=len(genome),
            )
            continue
        observed = spacer.upper()
        expected = genome[spacer_start:spacer_end]
        if observed != expected:
            _append_issue(
                issues,
                "spacer_sequence_mismatch",
                "spacer sequence does not match its independently derived genome interval",
                component="spacer",
                index=index,
                genomic_start=spacer_start,
                genomic_end=spacer_end,
                **_sequence_mismatch_details(expected, observed)
            )

    array_start = None
    array_end = None
    if (
        repeat_sequences
        and genomic_starts
        and genomic_starts[0] is not None
        and genomic_starts[-1] is not None
        and len(genomic_starts) == len(repeat_sequences)
        and isinstance(repeat_sequences[-1], str)
        and repeat_sequences[-1]
    ):
        array_start = genomic_starts[0]
        array_end = genomic_starts[-1] + len(repeat_sequences[-1])
        if (
            array_start < 0
            or array_end > len(genome)
            or array_end <= array_start
        ):
            _append_issue(
                issues,
                "array_out_of_bounds",
                "whole-array interval lies outside the original genome",
                component="array",
                genomic_start=array_start,
                genomic_end=array_end,
                genome_length=len(genome),
            )

    can_reconstruct = (
        array_start is not None
        and array_end is not None
        and 0 <= array_start < array_end <= len(genome)
        and len(spacer_sequences) == len(repeat_sequences) - 1
        and all(isinstance(sequence, str) for sequence in repeat_sequences)
        and all(isinstance(sequence, str) for sequence in spacer_sequences)
        and all(start is not None for start in genomic_starts)
        and all(
            genomic_starts[index] < genomic_starts[index + 1]
            for index in range(len(genomic_starts) - 1)
        )
    )
    if can_reconstruct:
        pieces = []
        for index, repeat in enumerate(repeat_sequences):
            pieces.append(repeat)
            if index < len(spacer_sequences):
                pieces.append(spacer_sequences[index])
        observed = "".join(pieces).upper()
        expected = genome[array_start:array_end]
        if observed != expected:
            _append_issue(
                issues,
                "array_sequence_mismatch",
                "interleaved repeats and spacers do not reconstruct the genome interval",
                component="array",
                genomic_start=array_start,
                genomic_end=array_end,
                **_sequence_mismatch_details(expected, observed)
            )

    return ReconstructionResult(
        genome_length=len(genome),
        coordinate_offset=coordinate_offset,
        local_repeat_starts=local_starts,
        genomic_repeat_starts=tuple(genomic_starts),
        repeat_count=len(repeat_sequences),
        spacer_count=len(spacer_sequences),
        array_start=array_start,
        array_end=array_end,
        issues=tuple(issues),
    )


def validate_candidate_reconstruction(
    genome_sequence,
    candidate,
    coordinate_offset=0,
):
    """Validate any candidate exposing the standard CRISPRidentify fields."""

    return validate_reconstruction(
        genome_sequence=genome_sequence,
        repeat_starts=candidate.list_repeat_starts,
        repeats=candidate.list_repeats,
        spacers=candidate.list_spacers,
        coordinate_offset=coordinate_offset,
    )
