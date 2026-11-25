"""Sequence utilities with explicit IUPAC DNA semantics."""


_UPPERCASE_IUPAC_COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}

IUPAC_COMPLEMENTS = dict(_UPPERCASE_IUPAC_COMPLEMENTS)
IUPAC_COMPLEMENTS.update(
    {
        nucleotide.lower(): complement.lower()
        for nucleotide, complement in _UPPERCASE_IUPAC_COMPLEMENTS.items()
    }
)
IUPAC_COMPLEMENTS.update({" ": " ", "-": "-", ".": ".", "+": "+"})


def reverse_complement(sequence):
    """Return a case-preserving reverse complement of an IUPAC sequence.

    Alignment display characters used by CRISPRidentify are preserved. Any
    other symbol is rejected instead of being silently emitted unchanged.
    """

    invalid = sorted(set(sequence) - set(IUPAC_COMPLEMENTS))
    if invalid:
        raise ValueError(
            "sequence contains unsupported non-IUPAC symbols: {}".format(
                ", ".join(repr(symbol) for symbol in invalid)
            )
        )
    return "".join(IUPAC_COMPLEMENTS[nucleotide] for nucleotide in reversed(sequence))
