"""Source-bound, scope-local occurrence correspondence (not boundary truth).

The 50% relation is deliberately separate from reconstruction's 3 bp / 80%
core correspondence. Counting complete ordered mappings saturates at two;
there is no truncation of candidate geometry or enumeration of mapping paths.
"""
import hashlib
import json
import math

RELATION = "unique_ordered_reciprocal_half_overlap_v1"


def checked_geometry(value, source_length=None):
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError("occurrence protection geometry is malformed")
    result = []
    for pair in value:
        if (not isinstance(pair, (list, tuple)) or len(pair) != 2
                or any(type(x) is not int for x in pair) or not 0 <= pair[0] < pair[1]
                or (result and pair[0] < result[-1][1])
                or (source_length is not None and pair[1] > source_length)):
            raise ValueError("occurrence protection intervals are invalid")
        result.append(tuple(pair))
    return tuple(result)


def occurrence_relation(reference, proposed):
    """Return exact map status and integer evidence.

    Reciprocal fractions are intersection / each bound interval length; this
    also makes failed edges independently derivable when no full map exists.
    No float rounding participates in the inclusive half-overlap threshold.
    """
    reference, proposed = checked_geometry(reference), checked_geometry(proposed)
    layers, unmatched, partial, inside_spacer = [], [], [], []
    lower = 0
    for ordinal, (start, end) in enumerate(reference):
        while lower < len(proposed) and proposed[lower][1] <= start:
            lower += 1
        options, overlaps = [], []
        for index in range(lower, len(proposed)):
            a, b = proposed[index]
            if a >= end:
                break
            overlap = min(end, b) - max(start, a)
            overlaps.append(index)
            if 2 * overlap >= end - start and 2 * overlap >= b - a:
                options.append(index)
        if not options:
            unmatched.append(ordinal)
            if overlaps:
                partial.append(ordinal)
            if (0 < lower < len(proposed) and proposed[lower - 1][1] <= start
                    and end <= proposed[lower][0]):
                inside_spacer.append(ordinal)
        if not layers:
            current = {index: (1, None) for index in options}
        else:
            previous = layers[-1]
            previous_indexes = sorted(previous)
            position, count, predecessor = 0, 0, None
            current = {}
            for index in options:
                while position < len(previous_indexes) and previous_indexes[position] < index:
                    old = previous_indexes[position]
                    ways = previous[old][0]
                    predecessor = old if count == 0 and ways == 1 else None
                    count = min(2, count + ways)
                    position += 1
                if count:
                    current[index] = (count, predecessor if count == 1 else None)
        layers.append(current)
    ways = min(2, sum(value[0] for value in layers[-1].values()))
    mapping = []
    if ways == 1:
        index = next(iter(layers[-1]))
        for layer in reversed(layers):
            mapping.append(index)
            index = layer[index][1]
        mapping.reverse()
    pairs = []
    for old, new in enumerate(mapping):
        a, b = reference[old], proposed[new]
        overlap = min(a[1], b[1]) - max(a[0], b[0])
        pairs.append({"reference_index": old, "candidate_index": new,
                      "intersection_length": overlap,
                      "boundary_shifts": [b[0] - a[0], b[1] - a[1]]})
    return {"rule": RELATION, "status": ("unique" if ways == 1 else
            "ambiguous_mapping" if ways else "no_complete_mapping"),
            "mapping_count_capped_at_two": ways, "mapping": mapping,
            "pairs": pairs, "reference_indexes_without_edge": unmatched,
            "reference_indexes_with_partial_overlap_only": partial,
            "reference_indexes_inside_candidate_spacers": inside_spacer}


def _hash(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def occurrence_snapshot(row):
    """Compact immutable-value observation; never alias candidate lists."""
    candidate = row[1]
    starts = tuple(candidate.list_repeat_starts)
    repeats = tuple(candidate.list_repeats)
    spacers = tuple(candidate.list_spacers)
    gapped = tuple(candidate.list_repeats_gaped)
    geometry = checked_geometry(tuple((start, start + len(repeat))
                                     for start, repeat in zip(starts, repeats)))
    if len(starts) != len(repeats) or len(spacers) != len(repeats) - 1:
        raise ValueError("occurrence protection candidate cardinality differs")
    digest = hashlib.sha256()
    for index, repeat in enumerate(repeats):
        digest.update(repeat.encode("ascii"))
        if index < len(spacers):
            if len(spacers[index]) != geometry[index + 1][0] - geometry[index][1]:
                raise ValueError("occurrence protection spacer length differs")
            digest.update(spacers[index].encode("ascii"))
    return {"repeat_intervals": [list(pair) for pair in geometry],
            "score": float(row[0]),
            "representation_sha256": _hash((starts, repeats, gapped, spacers)),
            "source_span_sha256": digest.hexdigest()}


def validate_snapshot(snapshot, source=None, source_length=None):
    required = {"repeat_intervals", "score", "representation_sha256", "source_span_sha256"}
    if not isinstance(snapshot, dict) or set(snapshot) != required:
        raise ValueError("occurrence protection snapshot is malformed")
    geometry = checked_geometry(snapshot["repeat_intervals"],
                                len(source) if source is not None else source_length)
    score = snapshot["score"]
    if type(score) not in (float, int) or not math.isfinite(score) or not .5 <= score <= 1:
        raise ValueError("occurrence protection requires a qualified own score")
    for key in ("representation_sha256", "source_span_sha256"):
        value = snapshot[key]
        if (not isinstance(value, str) or len(value) != 64
                or any(char not in "0123456789abcdef" for char in value)):
            raise ValueError("occurrence protection digest is malformed")
    if source is not None:
        span = source[geometry[0][0]:geometry[-1][1]]
        if hashlib.sha256(span.encode("ascii")).hexdigest() != snapshot["source_span_sha256"]:
            raise ValueError("occurrence protection source span differs")
    return geometry
