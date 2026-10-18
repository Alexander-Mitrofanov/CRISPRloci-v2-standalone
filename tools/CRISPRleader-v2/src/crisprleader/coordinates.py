"""Source-forward circular arcs and explicitly oriented flanking windows.

Extraction windows are a sampling policy, never an inferred biological leader.
Every projected interval is returned in source-forward arc order, including
negative-strand hypotheses; apply the candidate strand only after concatenation.
"""
from __future__ import annotations

from collections.abc import Sequence

from crisprleader.errors import InputError
from crisprleader.models import (
    AnalysisConfig, ArrayAnnotation, FlankCandidate, Segment, SourceRecord,
)

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    """Return upper-case IUPAC DNA reverse complement; reject invalid symbols."""
    if not isinstance(sequence, str) or not set(sequence.upper()) <= IUPAC_DNA:
        raise InputError("sequence must contain only IUPAC DNA symbols (no U or gaps)")
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def validate_arc(
    segments: Sequence[Segment], source_length: int, *, circular: bool = True,
    allow_empty: bool = False, label: str = "segments",
) -> tuple[Segment, ...]:
    """Check a contiguous, at-most-once forward traversal of a source record."""
    if type(source_length) is not int or source_length < 1:
        raise InputError("source length must be a positive integer")
    if not isinstance(segments, (list, tuple)):
        raise InputError(f"{label} must be ordered segments")
    result = tuple(segments)
    if not result and allow_empty:
        return result
    if len(result) not in (1, 2) or any(not isinstance(s, Segment) for s in result):
        raise InputError(f"{label} must have one interval or two origin-crossing intervals")
    if any(s.end > source_length for s in result):
        raise InputError(f"{label} exceeds source length {source_length}")
    if len(result) == 2:
        first, last = result
        if not circular:
            raise InputError(f"{label} crosses origin without explicit circular topology")
        if first.end != source_length or last.start != 0 or last.end > first.start:
            raise InputError(f"{label} must form one contiguous origin-crossing arc without overlap")
    return result


def extract_segments(
    source_sequence: str, segments: Sequence[Segment], strand: str = "+",
) -> str:
    """Reconstruct an arc, then orient it; empty flanks reconstruct as empty DNA."""
    if not isinstance(source_sequence, str) or not source_sequence:
        raise InputError("source sequence must be nonempty DNA")
    if strand not in ("+", "-"):
        raise InputError("reconstruction strand must be + or -")
    arc = validate_arc(segments, len(source_sequence), allow_empty=True)
    dna = "".join(source_sequence[s.start:s.end] for s in arc).upper()
    if not set(dna) <= IUPAC_DNA:
        raise InputError("source segments contain non-IUPAC DNA symbols")
    return dna if strand == "+" else reverse_complement(dna)


def slice_arc(segments: Sequence[Segment], start: int, end: int) -> tuple[Segment, ...]:
    """Slice a validated forward arc using its concatenated offsets."""
    total = sum(s.end - s.start for s in segments)
    if type(start) is not int or type(end) is not int or not 0 <= start < end <= total:
        raise InputError("projected interval requires integer offsets 0 <= start < end <= length")
    result: list[Segment] = []
    offset = 0
    for segment in segments:
        size = segment.end - segment.start
        low, high = max(start, offset), min(end, offset + size)
        if low < high:
            result.append(Segment(segment.start + low - offset, segment.start + high - offset))
        offset += size
    return tuple(result)


def _arc(start: int, length: int, source_length: int) -> tuple[Segment, ...]:
    if length == 0:
        return ()
    start %= source_length
    end = start + length
    if end <= source_length:
        return (Segment(start, end),)
    return (Segment(start, source_length), Segment(0, end - source_length))


def validate_array_geometry(source: SourceRecord, array: ArrayAnnotation) -> None:
    """Validate geometry and observed anchors even for programmatic callers."""
    if not isinstance(source, SourceRecord) or not isinstance(array, ArrayAnnotation):
        raise InputError("extraction requires a source record and array annotation")
    if source.id != array.source_id:
        raise InputError(f"array {array.id} is bound to a different source ID")
    if source.topology not in ("linear", "circular", "unknown"):
        raise InputError("source topology must be linear, circular, or unknown")
    if source.topology == "circular" and (
        not isinstance(source.topology_provenance, str) or not source.topology_provenance.strip()
    ):
        raise InputError("circular topology requires explicit topology provenance")
    if not isinstance(source.sequence, str) or not source.sequence:
        raise InputError("source sequence must be nonempty IUPAC DNA")
    arc = validate_arc(array.segments, source.length, circular=source.topology == "circular")
    size = sum(s.end - s.start for s in arc)
    if len(array.terminal_repeats) != 2:
        raise InputError("array requires exactly two terminal-repeat anchors")
    first, last = array.terminal_repeats
    lengths = []
    for anchor in (first, last):
        validate_arc(anchor.segments, source.length, circular=source.topology == "circular")
        observed = extract_segments(source.sequence, anchor.segments)
        if not isinstance(anchor.sequence, str) or anchor.sequence.upper() != observed:
            raise InputError(f"array {array.id} terminal repeat does not match source sequence")
        lengths.append(len(observed))
    if sum(lengths) > size:
        raise InputError(f"array {array.id} terminal-repeat anchors overlap")
    if (first.segments != slice_arc(arc, 0, lengths[0])
            or last.segments != slice_arc(arc, size - lengths[1], size)):
        raise InputError(f"array {array.id} terminal repeats must be flush with array endpoints")
    if not array.consensus_repeat or not set(array.consensus_repeat.upper()) <= IUPAC_DNA:
        raise InputError("consensus repeat must be nonempty IUPAC DNA in source-forward frame")
    if array.strand not in ("+", "-", "?") or type(array.strand_constraint) is not bool:
        raise InputError("array strand or strand_constraint is invalid")
    if array.strand_constraint and array.strand == "?":
        raise InputError("an unknown strand cannot be a hard strand constraint")


def extract_candidates(
    source: SourceRecord, array: ArrayAnnotation, config: AnalysisConfig | None = None,
) -> tuple[FlankCandidate, FlankCandidate]:
    """Extract both hypotheses; circular windows never enter the focal array."""
    config = config or AnalysisConfig()
    if (type(config.flank_length) is not int or type(config.max_flank_length) is not int
            or not 1 <= config.flank_length <= config.max_flank_length):
        raise InputError("flank_length must be a positive integer at most max_flank_length")
    validate_array_geometry(source, array)
    array_size = sum(s.end - s.start for s in array.segments)
    candidates: list[FlankCandidate] = []
    for side, strand, anchor in (("left", "+", array.terminal_repeats[0]),
                                  ("right", "-", array.terminal_repeats[1])):
        if source.topology == "circular":
            available = source.length - array_size
            reason = "circular_array_boundary"
        else:
            available = (array.segments[0].start if side == "left"
                         else source.length - array.segments[-1].end)
            reason = ("unknown_topology_boundary" if source.topology == "unknown"
                      else "contig_boundary")
        length = min(config.flank_length, available)
        start = (array.segments[0].start - length if side == "left"
                 else array.segments[-1].end)
        segments = _arc(start, length, source.length)
        candidate = FlankCandidate(
            id=f"{array.id}:{side}", source_id=source.id, array_id=array.id,
            side=side, strand=strand, segments=segments,  # type: ignore[arg-type]
            sequence=extract_segments(source.sequence, segments, strand),
            requested_length=config.flank_length,
            truncated_reason=reason if length < config.flank_length else None,
            anchor=anchor,
            observed_repeat=extract_segments(source.sequence, anchor.segments, strand),
            consensus_repeat=(array.consensus_repeat.upper() if strand == "+"
                              else reverse_complement(array.consensus_repeat)),
            metadata={"topology": source.topology, "source_sha256": source.sha256, "source_length": source.length,
                      "extraction_policy": "adjacent-window-excluding-focal-array",
                      "window_is_biological_boundary": False},
        )
        candidates.append(candidate)
    return candidates[0], candidates[1]


def project_interval(
    candidate: FlankCandidate, start: int, end: int, include_repeat: bool = False,
) -> tuple[Segment, ...]:
    """Map an oriented flank/junction interval back to source-forward segments."""
    if not isinstance(candidate, FlankCandidate) or candidate.strand not in ("+", "-"):
        raise InputError("projection requires an oriented flank candidate")
    if type(include_repeat) is not bool:
        raise InputError("include_repeat must be boolean")
    if len(candidate.sequence) != sum(s.end - s.start for s in candidate.segments):
        raise InputError("candidate sequence length does not match its segments")
    segments = candidate.segments
    if include_repeat:
        if len(candidate.observed_repeat) != sum(
            s.end - s.start for s in candidate.anchor.segments
        ):
            raise InputError("observed repeat length does not match anchor segments")
        segments = (segments + candidate.anchor.segments if candidate.strand == "+"
                    else candidate.anchor.segments + segments)
    # Adjacent flank and repeat pieces can be consolidated without crossing origin.
    merged: list[Segment] = []
    for segment in segments:
        if merged and merged[-1].end == segment.start:
            merged[-1] = Segment(merged[-1].start, segment.end)
        else:
            merged.append(segment)
    if merged:
        validate_arc(merged, max(s.end for s in merged))
    total = sum(s.end - s.start for s in merged)
    if type(start) is not int or type(end) is not int or not 0 <= start < end <= total:
        raise InputError("projected interval requires integer offsets 0 <= start < end <= length")
    if candidate.strand == "-":
        start, end = total - end, total - start
    return slice_arc(merged, start, end)
