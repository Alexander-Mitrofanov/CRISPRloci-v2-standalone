"""Optional native matcher for CRISPRidentify's fixed fuzzy literal search.

The C implementation mirrors the match order and ENHANCEMATCH behavior of
``regex`` 2019.03.09 for the exact constraint used by detection.  Once built,
it implements that frozen contract independently of the installed ``regex``
package. Unsupported patterns or a missing/incompatible shared library retain
the Python fallback.
"""

import ctypes
from pathlib import Path

from components.fuzzy_alignment import snapshot_owned_fuzzy_match_state


ABI_VERSION = 1
EXPECTED_REGEX_VERSION = "2.5.27"
MIN_VALIDATED_PATTERN = 15
MAX_VALIDATED_PATTERN = 63
MAX_PER_TYPE = 3
MAX_C_INT_LENGTH = (2 ** 31 - 1) - MAX_PER_TYPE
_LITERAL_ALPHABET = frozenset("ACGTN")
_LIBRARY_PATH = Path(__file__).with_name("native") / "libliteral_fuzzy.so"


class NativeMatcherUnavailable(RuntimeError):
    pass


class _CLiteralMatch(ctypes.Structure):
    _fields_ = [
        ("start", ctypes.c_int),
        ("end", ctypes.c_int),
        ("substitutions", ctypes.c_int),
        ("insertions", ctypes.c_int),
        ("deletions", ctypes.c_int),
        ("substitution_pos", ctypes.c_int * MAX_PER_TYPE),
        ("insertion_pos", ctypes.c_int * MAX_PER_TYPE),
        ("deletion_pos", ctypes.c_int * MAX_PER_TYPE),
    ]


class LiteralMatch:
    """Small regex.Match-compatible view used by :class:`FuzzySearch`."""

    # The C matcher intentionally reproduces regex 2019.03.09, whose deletion
    # positions identify the text boundary before the deleted pattern base.
    deletion_position_shift = 1

    __slots__ = (
        "_text", "_start", "_end", "fuzzy_counts", "fuzzy_changes"
    )

    def __init__(self, match, text):
        self._text = text
        self._start = match.start
        self._end = match.end
        self.fuzzy_counts = (
            match.substitutions, match.insertions, match.deletions
        )
        self.fuzzy_changes = (
            list(match.substitution_pos[:match.substitutions]),
            list(match.insertion_pos[:match.insertions]),
            list(match.deletion_pos[:match.deletions]),
        )

    def _fuzzy_alignment_cache_snapshot(self):
        return snapshot_owned_fuzzy_match_state(
            self.fuzzy_counts, self.fuzzy_changes
        )

    def start(self):
        return self._start

    def end(self):
        return self._end

    def span(self):
        return self._start, self._end

    def group(self):
        return self._text[self._start:self._end]


def _load_library():
    try:
        library = ctypes.CDLL(str(_LIBRARY_PATH))
    except OSError as error:
        raise NativeMatcherUnavailable(str(error))

    abi_version = library.literal_fuzzy_abi_version
    abi_version.argtypes = []
    abi_version.restype = ctypes.c_int
    observed_abi = abi_version()
    if observed_abi != ABI_VERSION:
        raise NativeMatcherUnavailable(
            "native matcher ABI {} != {}".format(observed_abi, ABI_VERSION)
        )

    finditer_function = library.literal_fuzzy_finditer
    finditer_function.argtypes = [
        ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_int,
        ctypes.POINTER(_CLiteralMatch), ctypes.c_int,
    ]
    finditer_function.restype = ctypes.c_int
    return library, finditer_function


try:
    _LIBRARY, _FINDITER = _load_library()
except NativeMatcherUnavailable as _LOAD_ERROR:
    _LIBRARY = None
    _FINDITER = None
    LOAD_ERROR = str(_LOAD_ERROR)
else:
    LOAD_ERROR = None

AVAILABLE = _FINDITER is not None


def supports(pattern):
    """Return whether ``pattern`` is inside the differential-test contract."""
    return (
        AVAILABLE
        and MIN_VALIDATED_PATTERN <= len(pattern) <= MAX_VALIDATED_PATTERN
        and all(base in _LITERAL_ALPHABET for base in pattern)
    )


def finditer_many(patterns, text):
    """Return one ordered native match list per supported input pattern."""

    return list(iter_match_batches(patterns, text))


def iter_match_batches(patterns, text):
    """Yield copied matches for one pattern at a time, reusing native scratch.

    Consumers may release a pattern's matches after candidate materialization
    instead of retaining all raw matches alongside the completed candidates.
    Each yielded list owns its records independently of the scratch buffer.
    """

    patterns = tuple(patterns)
    if not patterns:
        return

    for pattern in patterns:
        if not supports(pattern):
            raise NativeMatcherUnavailable(
                "pattern is outside the validated native matcher contract"
            )

    try:
        encoded_patterns = tuple(
            pattern.encode("ascii") for pattern in patterns
        )
        text_bytes = text.encode("ascii")
    except UnicodeEncodeError as error:
        raise NativeMatcherUnavailable(str(error))

    if (
            len(text_bytes) > MAX_C_INT_LENGTH
            or any(
                len(pattern_bytes) > MAX_C_INT_LENGTH
                for pattern_bytes in encoded_patterns
            )):
        raise NativeMatcherUnavailable(
            "pattern or text is too large for the native matcher ABI"
        )

    capacity = max(
        max(
            4,
            len(text) // max(1, len(pattern) - MAX_PER_TYPE) + 4,
        )
        for pattern in patterns
    )
    if capacity > MAX_C_INT_LENGTH:
        raise NativeMatcherUnavailable(
            "match capacity is too large for the native matcher ABI"
        )

    output = (_CLiteralMatch * capacity)()
    for pattern_bytes in encoded_patterns:
        while True:
            count = _FINDITER(
                pattern_bytes,
                len(pattern_bytes),
                text_bytes,
                len(text_bytes),
                output,
                capacity,
            )
            if count >= 0:
                break
            capacity = -count
            if capacity > MAX_C_INT_LENGTH:
                raise NativeMatcherUnavailable(
                    "match capacity is too large for the native matcher ABI"
                )
            output = (_CLiteralMatch * capacity)()

        # Copy every C record before the next pattern overwrites the scratch.
        yield [
                LiteralMatch(output[index], text) for index in range(count)
        ]


def finditer(pattern, text):
    """Return ordered native matches for a supported literal DNA pattern."""
    return finditer_many((pattern,), text)[0]


def signature(matches):
    """Return the strong differential signature used by regression tests."""
    return [
        (match.span(), match.group(), match.fuzzy_counts, match.fuzzy_changes)
        for match in matches
    ]


__all__ = [
    "AVAILABLE", "LOAD_ERROR", "LiteralMatch", "NativeMatcherUnavailable",
    "finditer", "finditer_many", "iter_match_batches", "signature", "supports",
]
