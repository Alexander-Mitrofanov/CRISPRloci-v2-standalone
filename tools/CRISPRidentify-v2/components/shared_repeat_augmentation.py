"""Conservative source-backed augmentation of overlapping array paths.

This module is deliberately independent of pipeline policy.  It combines
candidate paths only inside an existing locus and only when exact repeat
intervals prove that the paths belong to one shared-repeat graph component.
Coordinates are zero-based and end-exclusive; circular wrapping is not
inferred.
"""

from collections import defaultdict

from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
    validate_candidate_source_reconstruction,
)


MINIMUM_SHARED_PATH_SPACER_LENGTH = 10
MAXIMUM_SHARED_PATH_SPACER_LENGTH = 140


def _candidate_repeat_path(full_dna, candidate):
    validate_candidate_source_reconstruction(
        full_dna,
        candidate.list_repeats,
        candidate.list_repeat_starts,
        candidate.list_spacers,
    )
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts, candidate.list_repeats
        )
    )


def _source_geometry_key(path):
    return tuple(path)


def _shared_interval_components(paths):
    """Return deterministic connected components of candidate indexes."""

    interval_owners = defaultdict(list)
    for candidate_index, path in enumerate(paths):
        for interval in set(path):
            interval_owners[interval].append(candidate_index)

    adjacency = [set() for _ in paths]
    for owners in interval_owners.values():
        if len(owners) < 2:
            continue
        first_owner = owners[0]
        for owner in owners[1:]:
            adjacency[first_owner].add(owner)
            adjacency[owner].add(first_owner)

    components = []
    unseen = set(range(len(paths)))
    while unseen:
        first = min(unseen)
        stack = [first]
        unseen.remove(first)
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            neighbours = sorted(
                adjacency[current].intersection(unseen), reverse=True
            )
            for neighbour in neighbours:
                unseen.remove(neighbour)
                stack.append(neighbour)
        components.append(tuple(sorted(component)))
    return tuple(components)


def _component_union_geometry(paths, component):
    intervals = sorted({
        interval
        for candidate_index in component
        for interval in paths[candidate_index]
    })

    for previous, current in zip(intervals, intervals[1:]):
        if current[0] < previous[1]:
            return None, "inconsistent_overlap"

    spacer_lengths = [
        current[0] - previous[1]
        for previous, current in zip(intervals, intervals[1:])
    ]
    if any(
            spacer_length < MINIMUM_SHARED_PATH_SPACER_LENGTH
            or spacer_length > MAXIMUM_SHARED_PATH_SPACER_LENGTH
            for spacer_length in spacer_lengths):
        return None, "spacer_out_of_bounds"

    input_edges = {
        edge
        for candidate_index in component
        for edge in zip(
            paths[candidate_index], paths[candidate_index][1:]
        )
    }
    if any(
            edge not in input_edges
            for edge in zip(intervals, intervals[1:])):
        return None, "unsupported_consecutive_edge"

    return tuple(intervals), None


def augment_source_backed_shared_repeat_candidates(
        full_dna, candidates_by_locus):
    """Append source-valid union variants for shared-repeat path components.

    Every original candidate is retained in its original order.  Candidates
    are connected only by an identical half-open repeat interval inside the
    same locus.  A connected component contributes at most one union variant,
    and only when all union intervals are mutually consistent, every derived
    spacer is 10..140 bases, and every consecutive union edge occurs in at
    least one original path.  Variants duplicate by source geometry are not
    appended.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    if not isinstance(candidates_by_locus, dict):
        raise ValueError(
            "shared-repeat augmentation requires a candidate mapping"
        )

    augmented = {}
    stats = {
        "locus_count": len(candidates_by_locus),
        "candidate_count": 0,
        "connected_component_count": 0,
        "multi_candidate_component_count": 0,
        "variants_built": 0,
        "variants_added": 0,
        "duplicate_variants_skipped": 0,
        "inconsistent_overlap_rejections": 0,
        "spacer_out_of_bounds_rejections": 0,
        "unsupported_consecutive_edge_rejections": 0,
    }

    for locus_key, candidates in candidates_by_locus.items():
        if not isinstance(candidates, (tuple, list)):
            raise ValueError(
                "each shared-repeat locus must contain a candidate list"
            )
        originals = list(candidates)
        augmented[locus_key] = originals
        paths = tuple(
            _candidate_repeat_path(full_dna, candidate)
            for candidate in originals
        )
        stats["candidate_count"] += len(originals)

        seen_geometries = {
            _source_geometry_key(path) for path in paths
        }
        components = _shared_interval_components(paths)
        stats["connected_component_count"] += len(components)
        for component in components:
            if len(component) < 2:
                continue
            stats["multi_candidate_component_count"] += 1
            union_geometry, rejection = _component_union_geometry(
                paths, component
            )
            if rejection is not None:
                stats[rejection + "_rejections"] += 1
                continue

            union_key = _source_geometry_key(union_geometry)
            if union_key in seen_geometries:
                stats["duplicate_variants_skipped"] += 1
                continue

            seed_index = max(
                component,
                key=lambda candidate_index: (
                    len(paths[candidate_index]), -candidate_index
                ),
            )
            starts = tuple(start for start, _end in union_geometry)
            ends = tuple(end for _start, end in union_geometry)
            variant = build_source_backed_candidate_from_boundaries(
                full_dna, originals[seed_index], starts, ends
            )
            validate_candidate_source_reconstruction(
                full_dna,
                variant.list_repeats,
                variant.list_repeat_starts,
                variant.list_spacers,
            )
            stats["variants_built"] += 1
            augmented[locus_key].append(variant)
            seen_geometries.add(union_key)
            stats["variants_added"] += 1

    return augmented, stats


__all__ = [
    "MAXIMUM_SHARED_PATH_SPACER_LENGTH",
    "MINIMUM_SHARED_PATH_SPACER_LENGTH",
    "augment_source_backed_shared_repeat_candidates",
]
