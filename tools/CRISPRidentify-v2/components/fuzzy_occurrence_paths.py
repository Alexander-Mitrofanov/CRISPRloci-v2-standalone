"""Source-replayed occurrence paths for fuzzy literal detection.

``regex.finditer`` and the native matcher intentionally implement a
leftmost, non-overlapping fuzzy-search contract.  A fuzzy match that begins a
few bases before a source-exact repeat can therefore consume its start and
hide the exact phase from later matches.  This module keeps that historical
path and adds a deliberately narrow exact-occurrence overlay:

* exact occurrences are searched only in bounded windows around baseline
  fuzzy hits;
* every materialized path is source ordered and non-overlapping;
* the baseline fuzzy path is never replaced or ranked against an exact path;
* every maximal exact path is retained when enumeration is complete;
* discovery is bounded at 8,192 exact nodes and 256 maximal paths.  Reaching
  either ceiling marks the lattice incomplete and emits no exact alternative,
  leaving only the independently valid baseline path eligible downstream.

Exact discovery scans at most one ``O(repeat_length)`` neighborhood per
baseline hit, rather than trying fuzzy alignment at every source position.
This is not an all-start fuzzy enumerator and does not claim universal repeat
recall.  Its contract is the proven overlapping-match suppression case.
"""

import hashlib
from bisect import bisect_left
from collections import namedtuple

from components.fuzzy_alignment import snapshot_owned_fuzzy_match_state


OCCURRENCE_PATH_SCHEMA = "crispridentify.fuzzy-occurrence-lattice.v1"
MAX_MATERIALIZED_EXACT_PATHS = 256
MAX_EXACT_OVERLAY_NODES = 8192


OccurrencePathNode = namedtuple(
    "OccurrencePathNode",
    "node_id start end sequence fuzzy_counts origin",
)
OccurrencePathRecord = namedtuple(
    "OccurrencePathRecord",
    "path_id kind node_ids",
)
OccurrencePathLattice = namedtuple(
    "OccurrencePathLattice",
    (
        "schema lattice_id repeat_candidate nodes compatibility_edges paths "
        "exact_paths_complete"
    ),
)
OccurrencePathPlan = namedtuple(
    "OccurrencePathPlan",
    "match_paths path_ids lattice",
)


class ExactLiteralMatch:
    """Minimal ``regex.Match``-compatible view of one exact source slice."""

    deletion_position_shift = 1
    __slots__ = ("_text", "_start", "_end", "fuzzy_counts", "fuzzy_changes")

    def __init__(self, text, start, end):
        self._text = text
        self._start = start
        self._end = end
        self.fuzzy_counts = (0, 0, 0)
        self.fuzzy_changes = ([], [], [])

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


def _match_signature(match):
    return (
        match.start(),
        match.end(),
        match.group(),
        tuple(match.fuzzy_counts),
    )


def _validate_match(match, text):
    start = match.start()
    end = match.end()
    if (
            not isinstance(start, int)
            or isinstance(start, bool)
            or not isinstance(end, int)
            or isinstance(end, bool)
            or start < 0
            or end <= start
            or end > len(text)):
        raise ValueError("fuzzy occurrence has invalid source coordinates")
    if text[start:end] != match.group():
        raise ValueError("fuzzy occurrence does not reconstruct its source")
    counts = tuple(match.fuzzy_counts)
    if (
            len(counts) != 3
            or any(
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 0
                for value in counts
            )):
        raise ValueError("fuzzy occurrence has invalid edit counts")


def _validate_non_overlapping_path(matches, text):
    previous_end = 0
    for match in matches:
        _validate_match(match, text)
        if match.start() < previous_end:
            raise ValueError(
                "overlapping fuzzy occurrences cannot form source spacers"
            )
        previous_end = match.end()


def _bounded_exact_occurrences(
        repeat, text, baseline_matches,
        max_nodes=MAX_EXACT_OVERLAY_NODES):
    """Return locally sampled exact hits and an explicit completeness flag."""

    repeat_length = len(repeat)
    if len(baseline_matches) > max_nodes:
        return (), False

    per_hit_budget = max(1, max_nodes // max(1, len(baseline_matches)))
    exact = {}
    complete = True
    for match in baseline_matches:
        # An exact [start, start + repeat_length) interval overlaps the
        # baseline interval iff start is in this bounded start domain.
        window_start = max(0, match.start() - repeat_length + 1)
        window_end = min(len(text), match.end() + repeat_length - 1)
        if window_end - window_start < repeat_length:
            continue
        local = []
        found = text.find(repeat, window_start, window_end)
        while found != -1:
            end = found + repeat_length
            if found < match.end() and match.start() < end:
                local.append((found, end))
            found = text.find(repeat, found + 1, window_end)
        if len(local) > per_hit_budget:
            complete = False
            # Preserve the phases nearest the baseline geometry throughout
            # the source instead of exhausting a global budget at the left
            # edge of a low-complexity tract.
            local = sorted(
                local,
                key=lambda span: (
                    abs(span[0] - match.start()), span[0], span[1]
                ),
            )[:per_hit_budget]
        for start, end in local:
            exact[(start, end)] = ExactLiteralMatch(text, start, end)
    return tuple(exact[key] for key in sorted(exact)), complete


def _maximal_non_overlapping_paths(matches, limit):
    """Enumerate source-order maximal interval paths, deterministically.

    For a selected occurrence, only occurrences that cannot have another hit
    inserted before them are immediate successors.  Starting from every
    occurrence with no compatible predecessor therefore enumerates maximal
    independent paths directly instead of iterating over every subset.
    """

    matches = tuple(sorted(matches, key=_match_signature))
    if not matches:
        return (), True
    if all(
            left.end() <= right.start()
            for left, right in zip(matches, matches[1:])):
        return (matches,), True

    lengths = {match.end() - match.start() for match in matches}
    if len(lengths) != 1:
        raise ValueError("source-exact overlay intervals must have one length")
    interval_length = next(iter(lengths))
    starts = tuple(match.start() for match in matches)
    compatible_successors = []
    has_predecessor = []
    for index, match in enumerate(matches):
        first = bisect_left(starts, match.end(), index + 1)
        if first == len(matches):
            compatible_successors.append(())
        else:
            stop = bisect_left(
                starts, starts[first] + interval_length, first + 1
            )
            compatible_successors.append(tuple(range(first, stop)))
        has_predecessor.append(
            starts[0] + interval_length <= match.start()
        )

    paths = []
    complete = True
    for start_index, predecessor_present in enumerate(has_predecessor):
        if predecessor_present:
            continue
        if len(paths) >= limit:
            complete = False
            break

        path = [start_index]
        stack = [(start_index, 0)]
        while stack:
            match_index, successor_offset = stack[-1]
            successors = compatible_successors[match_index]
            if not successors:
                if len(paths) >= limit:
                    complete = False
                    break
                paths.append(tuple(matches[index] for index in path))
                stack.pop()
                path.pop()
                continue
            if successor_offset >= len(successors):
                stack.pop()
                path.pop()
                continue
            if len(paths) >= limit:
                complete = False
                break
            successor = successors[successor_offset]
            stack[-1] = (match_index, successor_offset + 1)
            stack.append((successor, 0))
            path.append(successor)
        if not complete:
            break

    unique = {}
    for path in paths:
        unique[tuple(_match_signature(match) for match in path)] = path
    return tuple(unique[key] for key in sorted(unique)), complete


def _node_id(sequence_start, match, origin):
    fields = (
        str(sequence_start + match.start()),
        str(sequence_start + match.end()),
        match.group(),
        ",".join(str(value) for value in match.fuzzy_counts),
        origin,
    )
    return "fuzzy-node:" + hashlib.sha256(
        "\x1f".join(fields).encode("utf-8")
    ).hexdigest()


def _path_id(lattice_seed, kind, node_ids):
    fields = (lattice_seed, kind) + tuple(node_ids)
    return "fuzzy-path:" + hashlib.sha256(
        "\x1f".join(fields).encode("utf-8")
    ).hexdigest()


def build_occurrence_path_plan(
        repeat, text, baseline_matches, sequence_start=0,
        max_exact_paths=MAX_MATERIALIZED_EXACT_PATHS,
        max_exact_nodes=MAX_EXACT_OVERLAY_NODES):
    """Build the baseline path plus source-exact phase alternatives.

    The returned match paths use window-relative coordinates.  Lattice nodes
    use absolute half-open source coordinates so provenance can be compared
    across cluster windows without trusting window offsets.
    """

    if not isinstance(repeat, str) or not repeat:
        raise ValueError("repeat candidate must be a non-empty string")
    if not isinstance(text, str):
        raise ValueError("fuzzy occurrence source must be a string")
    if (
            not isinstance(sequence_start, int)
            or isinstance(sequence_start, bool)
            or sequence_start < 0):
        raise ValueError("sequence start must be a non-negative integer")
    if (
            not isinstance(max_exact_paths, int)
            or isinstance(max_exact_paths, bool)
            or max_exact_paths <= 0):
        raise ValueError("exact path ceiling must be a positive integer")
    if (
            not isinstance(max_exact_nodes, int)
            or isinstance(max_exact_nodes, bool)
            or max_exact_nodes <= 0):
        raise ValueError("exact node ceiling must be a positive integer")

    baseline = tuple(baseline_matches)
    _validate_non_overlapping_path(baseline, text)
    exact_matches, exact_overlay_complete = _bounded_exact_occurrences(
        repeat, text, baseline, max_exact_nodes
    )
    exact_paths, path_enumeration_complete = _maximal_non_overlapping_paths(
        exact_matches, max_exact_paths
    )
    exact_paths_complete = (
        exact_overlay_complete and path_enumeration_complete
    )
    if not exact_paths_complete:
        # Truncated paths are not a candidate universe.  Retain their sampled
        # nodes plus the incompleteness bit in the lattice, but do not create
        # path records that could be mistaken for classification candidates.
        exact_paths = ()

    match_paths = [baseline]
    kinds = ["baseline_fuzzy"]
    baseline_signature = tuple(
        _match_signature(match) for match in baseline
    )
    seen = {baseline_signature}
    for path in exact_paths:
        signature = tuple(_match_signature(match) for match in path)
        if len(path) >= 2 and signature not in seen:
            _validate_non_overlapping_path(path, text)
            match_paths.append(path)
            kinds.append("source_exact_overlay")
            seen.add(signature)

    if len(match_paths) == 1 and exact_paths_complete:
        return OccurrencePathPlan((baseline,), (None,), None)

    node_origins = {}
    matches_by_signature = {}
    for match in baseline:
        signature = _match_signature(match)
        node_origins.setdefault(signature, set()).add("baseline_fuzzy")
        matches_by_signature[signature] = match
    for match in exact_matches:
        signature = _match_signature(match)
        node_origins.setdefault(signature, set()).add("source_exact")
        matches_by_signature[signature] = match

    nodes = []
    node_ids_by_signature = {}
    for signature in sorted(matches_by_signature):
        match = matches_by_signature[signature]
        origin = "+".join(sorted(node_origins[signature]))
        node_id = _node_id(sequence_start, match, origin)
        node_ids_by_signature[signature] = node_id
        nodes.append(OccurrencePathNode(
            node_id,
            sequence_start + match.start(),
            sequence_start + match.end(),
            match.group(),
            tuple(match.fuzzy_counts),
            origin,
        ))

    lattice_payload = [repeat]
    lattice_payload.extend(node.node_id for node in nodes)
    lattice_seed = hashlib.sha256(
        "\x1f".join(lattice_payload).encode("utf-8")
    ).hexdigest()
    lattice_id = "fuzzy-lattice:" + lattice_seed

    records = []
    path_ids = []
    for kind, path in zip(kinds, match_paths):
        node_ids = tuple(
            node_ids_by_signature[_match_signature(match)]
            for match in path
        )
        path_id = _path_id(lattice_seed, kind, node_ids)
        path_ids.append(path_id)
        records.append(OccurrencePathRecord(path_id, kind, node_ids))

    compatibility_edges = tuple(sorted({
        (left, right)
        for record in records
        for left, right in zip(record.node_ids, record.node_ids[1:])
    }))

    lattice = OccurrencePathLattice(
        OCCURRENCE_PATH_SCHEMA,
        lattice_id,
        repeat,
        tuple(nodes),
        compatibility_edges,
        tuple(records),
        exact_paths_complete,
    )
    return OccurrencePathPlan(
        tuple(match_paths), tuple(path_ids), lattice
    )


__all__ = [
    "ExactLiteralMatch",
    "MAX_EXACT_OVERLAY_NODES",
    "MAX_MATERIALIZED_EXACT_PATHS",
    "OCCURRENCE_PATH_SCHEMA",
    "OccurrencePathLattice",
    "OccurrencePathNode",
    "OccurrencePathPlan",
    "OccurrencePathRecord",
    "build_occurrence_path_plan",
]
