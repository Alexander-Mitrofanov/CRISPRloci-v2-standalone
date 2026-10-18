"""Sequence and clustering helpers shared by SpacerPlacer input adapters."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence


IUPAC_DNA_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)


def reverse_complement(sequence: str) -> str:
    """Return the IUPAC-aware reverse complement of a DNA sequence."""
    return sequence.translate(IUPAC_DNA_COMPLEMENT)[::-1]


def levenshtein_distance(seq1: str, seq2: str, max_distance: int | None = None) -> int:
    """Compute edit distance using linear memory and an optional cutoff.

    When ``max_distance`` is supplied, ``max_distance + 1`` is returned as soon
    as it is known that the sequences cannot be within the requested distance.
    This is substantially faster for the small thresholds used for spacer
    clustering.
    """
    if max_distance is not None and max_distance < 0:
        raise ValueError("max_distance must be non-negative")
    if seq1 == seq2:
        return 0
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1
    if max_distance is not None and len(seq2) - len(seq1) > max_distance:
        return max_distance + 1
    if not seq1:
        distance = len(seq2)
        return min(distance, max_distance + 1) if max_distance is not None else distance

    if max_distance is not None:
        outside = max_distance + 1
        previous = [outside] * (len(seq1) + 1)
        for column in range(min(len(seq1), max_distance) + 1):
            previous[column] = column
        for row_index, char2 in enumerate(seq2, start=1):
            current = [outside] * (len(seq1) + 1)
            if row_index <= max_distance:
                current[0] = row_index
            first_column = max(1, row_index - max_distance)
            last_column = min(len(seq1), row_index + max_distance)
            row_minimum = outside
            for column_index in range(first_column, last_column + 1):
                char1 = seq1[column_index - 1]
                value = min(
                    current[column_index - 1] + 1,
                    previous[column_index] + 1,
                    previous[column_index - 1] + (char1 != char2),
                )
                current[column_index] = value
                row_minimum = min(row_minimum, value)
            if row_minimum > max_distance:
                return outside
            previous = current
        return min(previous[-1], outside)

    previous = list(range(len(seq1) + 1))
    for row_index, char2 in enumerate(seq2, start=1):
        current = [row_index]
        for column_index, char1 in enumerate(seq1, start=1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[column_index] + 1,
                    previous[column_index - 1] + (char1 != char2),
                )
            )
        previous = current
    return previous[-1]


class DisjointSet:
    """Small deterministic union-find implementation."""

    def __init__(self, size: int):
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, first: int, second: int) -> None:
        root_first = self.find(first)
        root_second = self.find(second)
        if root_first == root_second:
            return
        if self.rank[root_first] < self.rank[root_second]:
            root_first, root_second = root_second, root_first
        self.parent[root_second] = root_first
        if self.rank[root_first] == self.rank[root_second]:
            self.rank[root_first] += 1

    def components(self) -> list[list[int]]:
        grouped: dict[int, list[int]] = defaultdict(list)
        for item in range(len(self.parent)):
            grouped[self.find(item)].append(item)
        return sorted(grouped.values(), key=lambda component: component[0])


def cluster_by_edit_distance(sequences: Iterable[str], max_distance: int) -> list[list[str]]:
    """Return transitive, order-stable edit-distance connected components."""
    if max_distance < 0:
        raise ValueError("max_distance must be non-negative")
    unique_sequences = list(dict.fromkeys(sequences))
    components = DisjointSet(len(unique_sequences))
    for first in range(len(unique_sequences)):
        for second in range(first + 1, len(unique_sequences)):
            if levenshtein_distance(
                unique_sequences[first], unique_sequences[second], max_distance=max_distance
            ) <= max_distance:
                components.union(first, second)
    return [[unique_sequences[index] for index in component] for component in components.components()]


def cluster_by_overlap(spacer_arrays: Sequence[Sequence[str]]) -> list[list[int]]:
    """Group array indices into transitive spacer-overlap components."""
    components = DisjointSet(len(spacer_arrays))
    first_array_by_spacer: dict[str, int] = {}
    for array_index, spacers in enumerate(spacer_arrays):
        for spacer in set(spacers):
            previous = first_array_by_spacer.setdefault(spacer, array_index)
            components.union(previous, array_index)
    return components.components()
