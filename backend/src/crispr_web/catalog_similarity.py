"""Literal, full-length nucleotide neighbourhoods; no alignment or tool execution."""

from __future__ import annotations

from collections.abc import Iterator

DNA = "ACGTRYSWKMBDHVN"
_COMPLEMENT = str.maketrans(DNA, "TGCAYRSWMKVHDBN")


def orientations(sequence: str) -> tuple[str, str]:
    return sequence, sequence.translate(_COMPLEMENT)[::-1]


def neighbours(sequence: str) -> Iterator[str]:
    """At most 2 * (1 + 14 * length) exact probes, independent of corpus size."""
    forward, reverse = orientations(sequence)
    for oriented in dict.fromkeys((forward, reverse)):
        yield oriented
        for index, base in enumerate(oriented):
            for replacement in DNA:
                if replacement != base:
                    yield oriented[:index] + replacement + oriented[index + 1:]


def sequence_matcher(query: str):
    choices = orientations(query)

    def matches(sequence: object) -> bool:
        if not isinstance(sequence, str) or len(sequence) != len(query):
            return False
        sequence = sequence.upper()
        if any(base not in DNA for base in sequence):
            return False
        for oriented in choices:
            differences = 0
            for left, right in zip(sequence, oriented):
                if left != right:
                    differences += 1
                    if differences > 1:
                        break
            else:
                return True
        return False

    return matches
