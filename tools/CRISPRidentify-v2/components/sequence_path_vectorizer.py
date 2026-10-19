"""Fast, exact EDeN vectorization for unweighted sequence path graphs.

CRISPRidentify represents repeats and spacers as simple path graphs and always
vectorizes them with ``r=3`` and ``d=4``. Building NetworkX graphs and running
the generic graph preprocessing pipeline dominates similarity scoring. This
module computes the same hashed, normalized features directly from each string.

Inputs outside this narrow contract are rejected; the legacy generic
EDeN/NetworkX vectorizer is not part of the v2 runtime.
"""

import ctypes
import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from components.legacy_hash import py37_hash
from components.legacy_hash import py37_hash_int_tuple


BITMASK_16 = 65535
BITMASK_32 = 4294967295
FEATURE_SIZE = 65537
RADIUS = 3
DISTANCE = 4
_EXPANDED_DISTANCE_ITEMS = ((0, 0), (1, 2), (2, 4), (3, 6), (4, 8))
_RADIUS_ITEMS = (
    (),
    ((0, 0),),
    ((0, 0), (1, 2)),
    ((0, 0), (1, 2), (2, 4)),
    ((0, 0), (1, 2), (2, 4), (3, 6)),
)

# ``_direct_feature_dict`` immediately masks every feature tuple hash to 16
# bits.  All tuple elements on that path are non-negative integers, so the low
# 16 bits of Python 3.7's recurrence depend only on the low 16 bits of each
# multiplication and xor.  Keeping the recurrence at that width is exact and
# avoids generic integer hashing and signed-64 conversion in the hottest loop.
_TUPLE_HASH_INITIAL_16 = 0x345678 & BITMASK_16
_TUPLE_HASH_MULTIPLIER_16 = 1000003 & BITMASK_16
_TUPLE_HASH_INCREMENT = 82520
_TUPLE_HASH_FINAL_16 = 97531 & BITMASK_16
_NATIVE_LIBRARY_PATH = (
    Path(__file__).with_name("native") / "libliteral_fuzzy.so"
)
_LEGACY_ASCII_LABEL_HASHES = (ctypes.c_uint32 * 128)(*(
    (py37_hash(chr(byte_value)) & BITMASK_16) + 1
    for byte_value in range(128)
))


def _load_native_path_features():
    """Load optional exact path-feature functions from the native matcher."""

    try:
        library = ctypes.CDLL(str(_NATIVE_LIBRARY_PATH))
        abi_version = library.path_features_abi_version
        abi_version.argtypes = []
        abi_version.restype = ctypes.c_int
        observed_abi = abi_version()
        if observed_abi not in (1, 2, 3, 4):
            return None, None, None, None, None
        event_function = library.path_feature_events
        event_function.argtypes = [
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.POINTER(ctypes.c_uint32),
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_uint64),
            ctypes.c_int,
        ]
        event_function.restype = ctypes.c_int
        count_function = None
        row_function = None
        dot_function = None
        ascii_row_function = None
        if observed_abi >= 2:
            count_function = library.path_feature_counts
            count_function.argtypes = event_function.argtypes
            count_function.restype = ctypes.c_int
        if observed_abi >= 3:
            row_function = library.path_feature_row
            row_function.argtypes = [
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.c_int,
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_double),
                ctypes.c_int,
            ]
            row_function.restype = ctypes.c_int
            dot_function = library.path_feature_dot
            dot_function.argtypes = [
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_double),
                ctypes.c_int,
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_double),
                ctypes.c_int,
            ]
            dot_function.restype = ctypes.c_double
        if observed_abi >= 4:
            ascii_row_function = library.path_feature_ascii_row
            ascii_row_function.argtypes = [
                ctypes.c_char_p,
                ctypes.c_int,
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_uint32),
                ctypes.POINTER(ctypes.c_double),
                ctypes.c_int,
            ]
            ascii_row_function.restype = ctypes.c_int
        return (
            event_function, count_function, row_function, dot_function,
            ascii_row_function,
        )
    except (AttributeError, OSError):
        return None, None, None, None, None


(
    _NATIVE_PATH_FEATURE_EVENTS,
    _NATIVE_PATH_FEATURE_COUNTS,
    _NATIVE_PATH_FEATURE_ROW,
    _NATIVE_PATH_FEATURE_DOT,
    _NATIVE_PATH_FEATURE_ASCII_ROW,
) = _load_native_path_features()
NATIVE_PATH_FEATURES_AVAILABLE = _NATIVE_PATH_FEATURE_EVENTS is not None
NATIVE_PATH_FEATURE_COUNTS_AVAILABLE = _NATIVE_PATH_FEATURE_COUNTS is not None
NATIVE_PACKED_PATH_FEATURES_AVAILABLE = (
    _NATIVE_PATH_FEATURE_ROW is not None
    and _NATIVE_PATH_FEATURE_DOT is not None
)
NATIVE_ASCII_PATH_FEATURES_AVAILABLE = (
    NATIVE_PACKED_PATH_FEATURES_AVAILABLE
    and _NATIVE_PATH_FEATURE_ASCII_ROW is not None
)


def _feature_hash3(first, second, third):
    """Return ``(py37_hash((a, b, c)) & 0xffff) + 1`` exactly."""

    result = (
        (_TUPLE_HASH_INITIAL_16 ^ (first & BITMASK_16))
        * _TUPLE_HASH_MULTIPLIER_16
    ) & BITMASK_16
    multiplier = (
        _TUPLE_HASH_MULTIPLIER_16 + _TUPLE_HASH_INCREMENT + 4
    ) & BITMASK_16
    result = ((result ^ (second & BITMASK_16)) * multiplier) & BITMASK_16
    multiplier = (multiplier + _TUPLE_HASH_INCREMENT + 2) & BITMASK_16
    result = ((result ^ (third & BITMASK_16)) * multiplier) & BITMASK_16
    return ((result + _TUPLE_HASH_FINAL_16) & BITMASK_16) + 1


def _feature_hash4(first, second, third, fourth):
    """Return ``(py37_hash((a, b, c, d)) & 0xffff) + 1`` exactly."""

    result = (
        (_TUPLE_HASH_INITIAL_16 ^ (first & BITMASK_16))
        * _TUPLE_HASH_MULTIPLIER_16
    ) & BITMASK_16
    multiplier = (
        _TUPLE_HASH_MULTIPLIER_16 + _TUPLE_HASH_INCREMENT + 6
    ) & BITMASK_16
    result = ((result ^ (second & BITMASK_16)) * multiplier) & BITMASK_16
    multiplier = (multiplier + _TUPLE_HASH_INCREMENT + 4) & BITMASK_16
    result = ((result ^ (third & BITMASK_16)) * multiplier) & BITMASK_16
    multiplier = (multiplier + _TUPLE_HASH_INCREMENT + 2) & BITMASK_16
    result = ((result ^ (fourth & BITMASK_16)) * multiplier) & BITMASK_16
    return ((result + _TUPLE_HASH_FINAL_16) & BITMASK_16) + 1


def _legacy_float_sum(values):
    """Accumulate left-to-right as Python 3.7 did for model features.

    Python 3.12 changed :func:`sum` for float iterables to use a more accurate
    algorithm.  Accuracy is normally desirable, but these normalized values
    are inputs to already-trained trees, so preserving the original bits is
    part of the model format.
    """

    total = 0.0
    for value in values:
        total += value
    return total


class _DenseKernel:
    """Small compatibility wrapper for the ``csr.dot(...).toarray()`` API."""

    __slots__ = ("_array",)

    def __init__(self, array):
        self._array = array

    def toarray(self):
        return self._array


class _TransposeFeatureMatrix:
    __slots__ = ("matrix",)

    def __init__(self, matrix):
        self.matrix = matrix


class _SparseDifference:
    """The part of SciPy's subtraction result used by compatibility tests."""

    __slots__ = ("nnz",)

    def __init__(self, nnz):
        self.nnz = nnz


class _PackedFeatureRow:
    """Native sparse row with sorted feature identifiers."""

    __slots__ = ("feature_ids", "feature_values", "length")

    def __init__(self, feature_ids, feature_values, length):
        self.feature_ids = feature_ids
        self.feature_values = feature_values
        self.length = length

    def __len__(self):
        return self.length

    def as_dict(self):
        return {
            int(self.feature_ids[index]): float(self.feature_values[index])
            for index in range(self.length)
        }


def _copy_packed_feature_row(feature_ids, feature_values, count):
    """Retain one exact-size native row without scalar ctypes marshalling."""

    exact_ids = (ctypes.c_uint32 * count)()
    exact_values = (ctypes.c_double * count)()
    ctypes.memmove(
        exact_ids, feature_ids, count * ctypes.sizeof(ctypes.c_uint32)
    )
    ctypes.memmove(
        exact_values, feature_values, count * ctypes.sizeof(ctypes.c_double)
    )
    return _PackedFeatureRow(exact_ids, exact_values, count)


def _sparse_dot(first_row, second_row):
    """Dot sparse rows in the same sorted-feature order as SciPy CSR."""
    if (
            isinstance(first_row, _PackedFeatureRow)
            and isinstance(second_row, _PackedFeatureRow)
            and _NATIVE_PATH_FEATURE_DOT is not None):
        return float(_NATIVE_PATH_FEATURE_DOT(
            first_row.feature_ids,
            first_row.feature_values,
            first_row.length,
            second_row.feature_ids,
            second_row.feature_values,
            second_row.length,
        ))
    if isinstance(first_row, _PackedFeatureRow):
        first_row = first_row.as_dict()
    if isinstance(second_row, _PackedFeatureRow):
        second_row = second_row.as_dict()
    result = 0.0
    for feature_id in sorted(first_row.keys() & second_row.keys()):
        result += first_row[feature_id] * second_row[feature_id]
    return result


class SparseFeatureMatrix:
    """Minimal sparse matrix needed by CRISPR similarity evaluation.

    The production path only needs row selection and a Gram matrix. Keeping
    normalized feature rows as small dictionaries avoids importing the larger
    SciPy runtime solely for those two operations.
    """

    __array_priority__ = 10000
    __slots__ = (
        "_dot_cache",
        "_feature_size",
        "_row_ids",
        "_rows",
    )

    def __init__(self, rows, feature_size=FEATURE_SIZE):
        self._rows = tuple(rows)
        self._feature_size = feature_size
        self._row_ids = tuple(range(len(self._rows)))
        self._dot_cache = {}

    @classmethod
    def _subset(cls, matrix, row_ids):
        subset = object.__new__(cls)
        subset._rows = matrix._rows
        subset._feature_size = matrix._feature_size
        subset._row_ids = tuple(row_ids)
        subset._dot_cache = matrix._dot_cache
        return subset

    @property
    def shape(self):
        return len(self._row_ids), self._feature_size

    @property
    def T(self):
        return _TransposeFeatureMatrix(self)

    @property
    def nnz(self):
        return sum(len(self._rows[row_id]) for row_id in self._row_ids)

    def __getitem__(self, row_selector):
        if isinstance(row_selector, slice):
            return self._subset(self, self._row_ids[row_selector])
        if isinstance(row_selector, (int, np.integer)):
            return self._subset(self, (self._row_ids[int(row_selector)],))
        return self._subset(
            self,
            (self._row_ids[int(index)] for index in row_selector),
        )

    def _dot_rows(self, first_row_id, second_row_id):
        if first_row_id > second_row_id:
            first_row_id, second_row_id = second_row_id, first_row_id
        cache_key = first_row_id, second_row_id
        try:
            return self._dot_cache[cache_key]
        except KeyError:
            result = _sparse_dot(
                self._rows[first_row_id],
                self._rows[second_row_id],
            )
            self._dot_cache[cache_key] = result
            return result

    def dot(self, other):
        if not isinstance(other, _TransposeFeatureMatrix):
            raise TypeError(
                "SparseFeatureMatrix only supports matrix.dot(matrix.T)"
            )
        other_matrix = other.matrix
        if self._rows is not other_matrix._rows:
            raise ValueError("feature matrices do not share the same rows")

        result = np.empty(
            (len(self._row_ids), len(other_matrix._row_ids)),
            dtype=np.float64,
        )
        symmetric = self._row_ids == other_matrix._row_ids
        for first_position, first_row_id in enumerate(self._row_ids):
            second_start = first_position if symmetric else 0
            for second_position in range(
                    second_start, len(other_matrix._row_ids)):
                value = self._dot_rows(
                    first_row_id,
                    other_matrix._row_ids[second_position],
                )
                result[first_position, second_position] = value
                if symmetric and second_position != first_position:
                    result[second_position, first_position] = value
        return _DenseKernel(result)

    def __rsub__(self, other):
        """Compare exactly with a SciPy matrix without importing SciPy here."""
        if getattr(other, "shape", None) != self.shape:
            return NotImplemented
        getrow = getattr(other, "getrow", None)
        if getrow is None:
            return NotImplemented

        different = 0
        for output_row, row_id in enumerate(self._row_ids):
            expected = self._rows[row_id]
            if isinstance(expected, _PackedFeatureRow):
                expected = expected.as_dict()
            actual_row = getrow(output_row)
            actual = dict(zip(actual_row.indices, actual_row.data))
            for feature_id in set(expected).union(actual):
                if expected.get(feature_id, 0.0) != actual.get(
                        feature_id, 0.0):
                    different += 1
        return _SparseDifference(different)


def supports_direct_vectorization(sequences):
    """Return whether all inputs satisfy the specialized path contract."""
    return bool(sequences) and all(
        isinstance(sequence, str) and sequence
        for sequence in sequences
    )


@lru_cache(maxsize=256)
def _candidate_frontiers(length):
    """Cache EDeN-compatible path frontiers, including set traversal order."""
    frontiers = []
    for vertex in range(length):
        vertex_frontiers = []
        for distance in range(DISTANCE + 1):
            if distance == 0:
                vertex_frontiers.append((vertex,))
                continue

            candidates = []
            left = vertex - distance
            right = vertex + distance
            if left >= 0:
                candidates.append(left)
            if right < length:
                candidates.append(right)
            vertex_frontiers.append(tuple(set(candidates)))
        frontiers.append(tuple(vertex_frontiers))
    return tuple(frontiers)


@lru_cache(maxsize=256)
def _native_candidate_frontiers(length):
    """Pack Python's exact frontier order for the native event generator."""

    if length > BITMASK_16:
        raise ValueError("native path-feature coordinates exceed 16 bits")
    packed = []
    for vertex_frontiers in _candidate_frontiers(length):
        for frontier in vertex_frontiers:
            value = 0
            if frontier:
                value = frontier[0] + 1
                if len(frontier) == 2:
                    value |= (frontier[1] + 1) << 16
            packed.append(value)
    return (ctypes.c_uint32 * len(packed))(*packed)


def _native_feature_group_order(neighborhood_hashes):
    """Return exact grouped feature counts, or ``None`` for Python fallback."""

    if not NATIVE_PATH_FEATURES_AVAILABLE:
        return None
    length = len(neighborhood_hashes)
    if length > BITMASK_16:
        return None

    flattened = []
    for radius_hashes in neighborhood_hashes:
        flattened.extend(radius_hashes)
        flattened.extend((0,) * (RADIUS + 1 - len(radius_hashes)))
    native_neighborhoods = (
        ctypes.c_uint32 * len(flattened)
    )(*flattened)
    native_frontiers = _native_candidate_frontiers(length)
    if _NATIVE_PATH_FEATURE_COUNTS is not None:
        count_capacity = max(1, length * 80)
        count_records = (ctypes.c_uint64 * count_capacity)()
        count = _NATIVE_PATH_FEATURE_COUNTS(
            native_neighborhoods,
            native_frontiers,
            length,
            count_records,
            count_capacity,
        )
        if count >= 0:
            grouped_features = [
                {} for _ in range((DISTANCE + 1) * (RADIUS + 1))
            ]
            group_order = []
            for index in range(count):
                packed = count_records[index]
                group = (packed >> 49) & 31
                feature = (packed >> 32) & 0x1FFFF
                feature_count = packed & BITMASK_32
                features = grouped_features[group]
                if not features:
                    group_order.append(features)
                features[feature] = float(feature_count)
            return group_order

    capacity = max(1, length * 40)
    events = (ctypes.c_uint64 * capacity)()
    count = _NATIVE_PATH_FEATURE_EVENTS(
        native_neighborhoods,
        native_frontiers,
        length,
        events,
        capacity,
    )
    if count < 0:
        capacity = -count
        events = (ctypes.c_uint64 * capacity)()
        count = _NATIVE_PATH_FEATURE_EVENTS(
            native_neighborhoods,
            native_frontiers,
            length,
            events,
            capacity,
        )
    if count < 0 or count > capacity:
        raise RuntimeError("native path-feature event count is invalid")

    grouped_features = [{} for _ in range((DISTANCE + 1) * (RADIUS + 1))]
    group_order = []
    for index in range(count):
        packed = events[index]
        group = (packed >> 34) & 31
        feature = (packed >> 17) & 0x1FFFF
        half_feature = packed & 0x1FFFF
        features = grouped_features[group]
        if not features:
            group_order.append(features)
        features[feature] = features.get(feature, 0.0) + 1.0
        features[half_feature] = features.get(half_feature, 0.0) + 1.0
    return group_order


def _neighborhood_radius_hashes(sequence):
    """Return only the neighborhood radii used by EDeN's r=3 features."""
    hash_value = py37_hash
    hash_int_tuple = py37_hash_int_tuple
    length = len(sequence)
    last_coordinate = 2 * length - 2
    expanded_vertex_hashes = [0] * (last_coordinate + 1)
    label_hashes = {}

    if length == 1:
        hashed_label = (hash_value(sequence[0]) & BITMASK_16) + 1
        expanded_vertex_hashes[0] = (
            hash_int_tuple(hashed_label, 0) & BITMASK_32) + 1
    else:
        edge_label_hash = (hash_value('-') & BITMASK_16) + 1
        edge_vertex_hash = (
            hash_int_tuple(edge_label_hash, 2) & BITMASK_32) + 1
        label_hashes['-'] = edge_label_hash
        last_position = length - 1

        for position, label in enumerate(sequence):
            try:
                hashed_label = label_hashes[label]
            except KeyError:
                hashed_label = (hash_value(label) & BITMASK_16) + 1
                label_hashes[label] = hashed_label

            degree = (
                1 if position == 0 or position == last_position else 2)
            coordinate = position * 2
            expanded_vertex_hashes[coordinate] = (
                hash_int_tuple(hashed_label, degree) & BITMASK_32) + 1
            if coordinate < last_coordinate:
                expanded_vertex_hashes[coordinate + 1] = edge_vertex_hash

    neighborhoods = []
    append_neighborhood = neighborhoods.append
    max_hash_depth = RADIUS * 2
    for root_coordinate in range(0, last_coordinate + 1, 2):
        eccentricity = max(
            root_coordinate, last_coordinate - root_coordinate)
        depth_limit = min(max_hash_depth, eccentricity)
        running_hash = 0xAAAAAAAA
        radius_hashes = []
        append_radius = radius_hashes.append

        for distance in range(depth_limit + 1):
            left = root_coordinate - distance
            right = root_coordinate + distance
            if (left >= 0 and right <= last_coordinate and
                    right != left):
                first = expanded_vertex_hashes[left]
                second = expanded_vertex_hashes[right]
                if first > second:
                    first, second = second, first
                distance_hash = (
                    hash_int_tuple(first, second) & BITMASK_32) + 1
            else:
                coordinate = left if left >= 0 else right
                distance_hash = (
                    hash_int_tuple(expanded_vertex_hashes[coordinate]) &
                    BITMASK_32) + 1

            running_hash ^= hash_int_tuple(
                running_hash, distance_hash, distance)
            if distance % 2 == 0:
                append_radius((running_hash & BITMASK_32) + 1)

        append_neighborhood(radius_hashes)

    return neighborhoods


def _native_packed_feature_row(sequence, output_buffers=None):
    """Return one exact normalized native sparse row, or None."""

    if not sequence.isascii():
        raise ValueError("native packed path-feature row requires ASCII input")
    if _NATIVE_PATH_FEATURE_ROW is None:
        return None
    length = len(sequence)
    if length > BITMASK_16:
        return None
    native_frontiers = _native_candidate_frontiers(length)
    required_capacity = max(1, length * 80)
    if output_buffers is None:
        capacity = required_capacity
        feature_ids = (ctypes.c_uint32 * capacity)()
        feature_values = (ctypes.c_double * capacity)()
    else:
        feature_ids, feature_values, capacity = output_buffers
        if capacity < required_capacity:
            raise ValueError("native packed path-feature buffer is too small")
    if _NATIVE_PATH_FEATURE_ASCII_ROW is not None and sequence.isascii():
        count = _NATIVE_PATH_FEATURE_ASCII_ROW(
            sequence.encode("ascii"),
            length,
            _LEGACY_ASCII_LABEL_HASHES,
            native_frontiers,
            feature_ids,
            feature_values,
            capacity,
        )
        if count != -2:
            if count <= 0 or count > capacity:
                raise RuntimeError("native ASCII path-feature row is invalid")
            return _copy_packed_feature_row(
                feature_ids, feature_values, count
            )

    neighborhood_hashes = _neighborhood_radius_hashes(sequence)
    flattened = []
    for radius_hashes in neighborhood_hashes:
        flattened.extend(radius_hashes)
        flattened.extend((0,) * (RADIUS + 1 - len(radius_hashes)))
    native_neighborhoods = (
        ctypes.c_uint32 * len(flattened)
    )(*flattened)
    count = _NATIVE_PATH_FEATURE_ROW(
        native_neighborhoods,
        native_frontiers,
        length,
        feature_ids,
        feature_values,
        capacity,
    )
    if count <= 0 or count > capacity:
        raise RuntimeError("native packed path-feature row is invalid")
    return _copy_packed_feature_row(feature_ids, feature_values, count)


def _direct_feature_dict(sequence):
    length = len(sequence)
    neighborhood_hashes = _neighborhood_radius_hashes(sequence)
    group_order = _native_feature_group_order(neighborhood_hashes)
    if group_order is None:
        frontiers = _candidate_frontiers(length)
        grouped_features = [
            [{} for _ in range(RADIUS + 1)]
            for _ in range(DISTANCE + 1)
        ]
        group_order = []
        append_group = group_order.append
        for vertex_v in range(length):
            labels_v = neighborhood_hashes[vertex_v]
            labels_v_length = len(labels_v)
            for external_distance, expanded_distance in _EXPANDED_DISTANCE_ITEMS:
                distance_groups = grouped_features[external_distance]
                for vertex_u in frontiers[vertex_v][external_distance]:
                    labels_u = neighborhood_hashes[vertex_u]
                    labels_u_length = len(labels_u)
                    limit = min(labels_v_length, labels_u_length)
                    for external_radius, expanded_radius in _RADIUS_ITEMS[limit]:
                        vertex_v_hash = labels_v[external_radius]
                        vertex_u_hash = labels_u[external_radius]
                        if vertex_v_hash < vertex_u_hash:
                            first_hash, second_hash = (
                                vertex_v_hash, vertex_u_hash)
                        else:
                            first_hash, second_hash = (
                                vertex_u_hash, vertex_v_hash)
                        feature = _feature_hash4(
                            first_hash, second_hash, expanded_radius,
                            expanded_distance,
                        )
                        half_feature = _feature_hash3(
                            vertex_u_hash, expanded_radius, expanded_distance,
                        )
                        features = distance_groups[external_radius]
                        if not features:
                            # Cross-group hash collisions make first-encounter
                            # normalization order observable in the final floats.
                            append_group(features)
                        features[feature] = features.get(feature, 0.0) + 1.0
                        features[half_feature] = (
                            features.get(half_feature, 0.0) + 1.0)

    # Match EDeN's per-radius/distance normalization and its overwrite
    # behavior when features from different groups collide after hashing.
    feature_vector = {}
    for features in group_order:
        sqrt_norm = math.sqrt(_legacy_float_sum(
            count * count for count in features.values()))
        for feature_id, count in features.items():
            feature_vector[feature_id] = float(count) / sqrt_norm

    sqrt_total_norm = math.sqrt(_legacy_float_sum(
        value * value for value in feature_vector.values()))
    return {
        feature_id: value / sqrt_total_norm
        for feature_id, value in feature_vector.items()
    }


def vectorize_sequence_paths(sequences):
    """Return exact EDeN ``r=3, d=4`` rows without requiring SciPy."""
    if not supports_direct_vectorization(sequences):
        raise ValueError('direct vectorization requires non-empty strings')

    packed_buffers = None
    if _NATIVE_PATH_FEATURE_ROW is not None:
        packed_lengths = [
            len(sequence) for sequence in sequences
            if len(sequence) <= BITMASK_16
        ]
        if packed_lengths:
            capacity = max(1, max(packed_lengths) * 80)
            packed_buffers = (
                (ctypes.c_uint32 * capacity)(),
                (ctypes.c_double * capacity)(),
                capacity,
            )

    rows = []
    for sequence in sequences:
        packed = _native_packed_feature_row(sequence, packed_buffers)
        rows.append(
            packed if packed is not None else _direct_feature_dict(sequence)
        )

    return SparseFeatureMatrix(
        rows,
        feature_size=FEATURE_SIZE,
    )
