"""Exact upstream canonical 4-mer counts, GC fraction and length representation."""
from itertools import product

import numpy as np
from numpy.typing import NDArray

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


KMERS = sorted({min(kmer, reverse_complement(kmer))
                for kmer in map("".join, product("ACGT", repeat=4))})
FEATURE_NAMES = sorted([*KMERS, "GC", "Length"])
INDEX = {name: i for i, name in enumerate(FEATURE_NAMES)}
# Resolve reverse-complement equivalence once, instead of for every window.
KMER_INDEX = {kmer: INDEX[min(kmer, reverse_complement(kmer))]
              for kmer in map("".join, product("ACGT", repeat=4))}


def matrix(sequences: list[str]) -> NDArray[np.float32]:
    # Upstream prepends an all-zero row and removes its prediction. Preserve it
    # as part of the frozen inference contract, including exclusive iteration 41.
    values = np.zeros((len(sequences) + 1, len(FEATURE_NAMES)), dtype=np.float32)
    for row, sequence in enumerate(sequences, start=1):
        columns = [KMER_INDEX[sequence[offset:offset + 4]] for offset in range(len(sequence) - 3)]
        values[row] = np.bincount(columns, minlength=len(FEATURE_NAMES))
        values[row, INDEX["Length"]] = len(sequence)
        values[row, INDEX["GC"]] = (sequence.count("G") + sequence.count("C")) / len(sequence)
    return values
