"""Coordinate helpers shared by sequence-backed CRISPR reporting.

Internally, repeat starts are zero-based indexes into the original input DNA.
Public report coordinates are one-based and inclusive.
"""


FORWARD_STRAND_LABELS = {
    "Forward",
    "Forward (Orientation was not computed)",
}

REVERSE_STRAND_LABELS = {
    "Reverse",
    "Reversed",
}


def is_forward_strand(strand):
    """Return whether *strand* should be displayed in input-sequence order."""

    if strand in FORWARD_STRAND_LABELS:
        return True
    if strand in REVERSE_STRAND_LABELS:
        return False
    raise ValueError("unknown CRISPR strand label: {!r}".format(strand))


def crispr_candidate_span(crispr_candidate):
    """Return the candidate's zero-based, half-open span in the input DNA."""

    starts = list(crispr_candidate.list_repeat_starts)
    repeats = list(crispr_candidate.list_repeats)
    spacers = list(crispr_candidate.list_spacers)
    if not starts or not repeats:
        raise ValueError("a CRISPR candidate must contain at least one repeat")
    if len(starts) != len(repeats):
        raise ValueError("repeat starts and repeat sequences must have equal lengths")
    if len(spacers) != len(repeats) - 1:
        raise ValueError("a CRISPR candidate must have one spacer between repeats")
    if any(
            not isinstance(start, int)
            or isinstance(start, bool)
            or start < 0
            for start in starts):
        raise ValueError("repeat starts must be non-negative integer offsets")
    if any(not isinstance(repeat, str) or not repeat for repeat in repeats):
        raise ValueError("repeat sequences must be non-empty strings")
    if any(not isinstance(spacer, str) for spacer in spacers):
        raise ValueError("spacer sequences must be strings")

    cursor = starts[0]
    array_start = cursor
    for index, (start, repeat) in enumerate(zip(starts, repeats)):
        if start != cursor:
            raise ValueError(
                "candidate repeat/spacer geometry does not reconstruct "
                "the supplied starts"
            )
        cursor += len(repeat)
        if index < len(spacers):
            cursor += len(spacers[index])
    return array_start, cursor


def crispr_candidate_output_coordinates(crispr_candidate):
    """Return one-based inclusive start/end and sequence length for a candidate."""

    start_zero_based, end_exclusive = crispr_candidate_span(crispr_candidate)
    return {
        "start": start_zero_based + 1,
        "end": end_exclusive,
        "length": end_exclusive - start_zero_based,
    }
