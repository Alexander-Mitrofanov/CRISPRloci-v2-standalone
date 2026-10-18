"""Version-independent reproduction of CRISPRidentify\047s Python 3.7 hashes.

The bundled EDeN-derived classifiers were trained from feature vectors whose
hashed feature identifiers were produced by 64-bit CPython 3.7 with
``PYTHONHASHSEED=0``.  CPython later changed both Unicode and tuple hashing, so
using :func:`hash` directly makes otherwise identical DNA depend on the Python
minor version.

Only the value types used by the sequence-path vectorizer are implemented:
ASCII strings, integers, and tuples containing those values.  The algorithms
below reproduce CPython 3.7.6\047s SipHash-2-4 string hash and tuple hash on a
64-bit build with the zero hash secret.
"""

from functools import lru_cache


_MASK_64 = (1 << 64) - 1
_SIGN_BIT_64 = 1 << 63
_PY_HASH_MODULUS = (1 << 61) - 1
_TUPLE_MULTIPLIER = 1000003


def _rotate_left_64(value, shift):
    value &= _MASK_64
    return ((value << shift) | (value >> (64 - shift))) & _MASK_64


def _sip_round(v0, v1, v2, v3):
    v0 = (v0 + v1) & _MASK_64
    v2 = (v2 + v3) & _MASK_64
    v1 = _rotate_left_64(v1, 13) ^ v0
    v3 = _rotate_left_64(v3, 16) ^ v2
    v0 = _rotate_left_64(v0, 32)
    v2 = (v2 + v1) & _MASK_64
    v0 = (v0 + v3) & _MASK_64
    v1 = _rotate_left_64(v1, 17) ^ v2
    v3 = _rotate_left_64(v3, 21) ^ v0
    v2 = _rotate_left_64(v2, 32)
    return v0, v1, v2, v3


def _as_signed_64(value):
    value &= _MASK_64
    if value >= _SIGN_BIT_64:
        value -= 1 << 64
    return -2 if value == -1 else value


def _siphash24_zero_key(data):
    """Return CPython 3.7\047s 64-bit SipHash-2-4 for a zero secret."""

    v0 = 0x736F6D6570736575
    v1 = 0x646F72616E646F6D
    v2 = 0x6C7967656E657261
    v3 = 0x7465646279746573

    offset = 0
    while len(data) - offset >= 8:
        message = int.from_bytes(data[offset:offset + 8], "little")
        offset += 8
        v3 ^= message
        v0, v1, v2, v3 = _sip_round(v0, v1, v2, v3)
        v0, v1, v2, v3 = _sip_round(v0, v1, v2, v3)
        v0 ^= message

    tail = data[offset:]
    final_block = (len(data) & 0xFF) << 56
    for index, byte in enumerate(tail):
        final_block |= byte << (8 * index)

    v3 ^= final_block
    v0, v1, v2, v3 = _sip_round(v0, v1, v2, v3)
    v0, v1, v2, v3 = _sip_round(v0, v1, v2, v3)
    v0 ^= final_block
    v2 ^= 0xFF
    for _ in range(4):
        v0, v1, v2, v3 = _sip_round(v0, v1, v2, v3)
    return _as_signed_64(v0 ^ v1 ^ v2 ^ v3)


def _hash_integer(value):
    if value >= 0:
        result = value % _PY_HASH_MODULUS
    else:
        result = -((-value) % _PY_HASH_MODULUS)
    return -2 if result == -1 else result


def _hash_tuple(values):
    result = 0x345678
    multiplier = _TUPLE_MULTIPLIER
    remaining = len(values)
    for value in values:
        remaining -= 1
        item_hash = py37_hash(value)
        result = ((result ^ (item_hash & _MASK_64)) * multiplier) & _MASK_64
        multiplier = (
            multiplier + 82520 + remaining + remaining
        ) & _MASK_64
    result = (result + 97531) & _MASK_64
    if result == _MASK_64:
        return -2
    return _as_signed_64(result)


@lru_cache(maxsize=65536)
def py37_hash_int_tuple(*values):
    """Hash an integer tuple with the Python 3.7 recurrence.

    Sequence-path vectorization hashes the same small integer tuples millions
    of times. Its tuple elements have already been reduced to integer graph
    labels, so routing every element back through py37_hash only repeats type
    dispatch and integer hashing. Keeping this specialized cache separate also
    prevents arbitrary input strings from occupying it.

    The recurrence and integer hash below are identical to _hash_tuple; only
    the already-known element type and memoization differ.
    """

    result = 0x345678
    multiplier = _TUPLE_MULTIPLIER
    remaining = len(values)
    for value in values:
        if not isinstance(value, int):
            raise TypeError("legacy integer-tuple hashing requires integers")
        remaining -= 1
        item_hash = _hash_integer(value)
        result = ((result ^ (item_hash & _MASK_64)) * multiplier) & _MASK_64
        multiplier = (
            multiplier + 82520 + remaining + remaining
        ) & _MASK_64
    result = (result + 97531) & _MASK_64
    if result == _MASK_64:
        return -2
    return _as_signed_64(result)


def py37_hash(value):
    """Hash a supported value exactly as 64-bit CPython 3.7.6 did.

    Strings are deliberately restricted to ASCII because CRISPRidentify\047s
    path-vectorizer inputs are DNA/alignment symbols.  Rejecting anything else
    prevents a false claim of reproducing CPython\047s variable-width Unicode
    storage hashing.
    """

    if isinstance(value, str):
        try:
            encoded = value.encode("ascii")
        except UnicodeEncodeError:
            raise ValueError("legacy feature hashing requires ASCII strings")
        if not encoded:
            return 0
        return _siphash24_zero_key(encoded)
    if isinstance(value, int):
        return _hash_integer(value)
    if isinstance(value, tuple):
        return _hash_tuple(value)
    raise TypeError(
        "legacy feature hashing does not support {}".format(
            type(value).__name__
        )
    )
