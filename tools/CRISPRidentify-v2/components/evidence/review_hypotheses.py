"""Bounded comparisons for human review of repeat boundaries and spacer tilings.

All intervals are zero-based and half-open. Source proposals carry no model
score or acceptance category. Retained observations refer back to the exact
scored row in the research inventory. No comparison changes that inventory.
"""

from collections import Counter, defaultdict
from statistics import median

from components.refinement.enhancement_search import (
    _source_repeat_alignment,
    _source_repeat_edit_distance,
)


def _valid(path, size):
    return (all(0 <= a < b <= size for a, b in path)
            and all(b < c for (_, b), (c, _) in zip(path, path[1:])))


def _gaps(path):
    return [c - b for (_, b), (c, _) in zip(path, path[1:])]


def _modal(repeats):
    counts = Counter(repeats)
    return min(counts, key=lambda rep: (-counts[rep], rep))


def _summary(source, path):
    repeats = [source[a:b] for a, b in path]
    motif = _modal(repeats)
    return {
        "modal_repeat": motif,
        "modal_copies": repeats.count(motif),
        "spacer_lengths": _gaps(path),
    }


def _single_copy(source, path, repeats, motif):
    # Local endpoints one base apart can encode an indel or a shifted repeat.
    # Examine sparse width outliers or a discordant copy in a short array.
    width = len(motif)
    width_counts = Counter(map(len, repeats))
    indexes = [i for i, rep in enumerate(repeats)
               if width_counts[width] >= 2 and abs(len(rep) - width) == 1]
    if len(path) <= 6 and repeats.count(motif) >= 2:
        indexes.extend(i for i, rep in enumerate(repeats)
                       if len(rep) == width and rep != motif and
                       sum(a != b for a, b in zip(rep, motif)) * 3 >= width)
    # Arrays with many outliers lack an isolated boundary comparison.
    if not 1 <= len(indexes) <= 2:
        return
    for i in indexes:
        a, b = path[i]
        for left, right in ((-1, -1), (1, 1), (0, -1), (1, 0), (0, 1), (-1, 0)):
            proposal = path[:i] + ((a + left, b + right),) + path[i + 1:]
            if not _valid(proposal, len(source)):
                continue
            rep = source[a + left:b + right]
            original = repeats[i]
            evidence = {
                "copy": i + 1,
                "primary_edit_distance": _source_repeat_edit_distance(motif, original),
                "proposed_edit_distance": _source_repeat_edit_distance(motif, rep),
                "primary_alignment_to_modal": _source_repeat_alignment((motif, original)),
                "proposed_alignment_to_modal": _source_repeat_alignment((motif, rep)),
            }
            if len(rep) == width:
                evidence["proposed_positional_mismatches"] = sum(
                    a != b for a, b in zip(rep, motif))
            yield proposal, "local_boundary", evidence


def _column(source, path, side, depth):
    indices = [a - depth if side == "left" else b + depth - 1 for a, b in path]
    if min(indices) < 0 or max(indices) >= len(source):
        return ()
    bases = tuple(source[i] for i in indices)
    return bases if all(base in "ACGT" for base in bases) else ()


def _edges(source, path, repeats):
    if len(path) < 4:
        return
    for side in ("left", "right"):
        columns = []
        for depth in range(1, 10):
            bases = _column(source, path, side, depth)
            if not bases:
                break
            counts = Counter(bases)
            columns.append(dict(sorted(counts.items())))
            if max(counts.values()) * 4 < 3 * len(bases):
                # A drop in support brackets a conserved extension on this edge.
                if depth >= 3:
                    delta = depth - 1
                    proposal = tuple((a - delta, b) if side == "left" else (a, b + delta)
                                     for a, b in path)
                    yield proposal, "conserved_flank", {"side": side, "columns": columns}
                break
        bases = _column(source, path, side, 1)
        if not bases or max(Counter(bases).values()) * 5 >= 3 * len(path):
            continue
        # A weak marginal column can still recur with distinct internal alleles.
        groups = defaultdict(Counter)
        for rep, base in zip(repeats, bases):
            groups[rep][base] += 1
        groups = {r: c for r, c in groups.items() if sum(c.values()) >= 3
                  and max(c.values()) * 4 >= 3 * sum(c.values())}
        if len({c.most_common(1)[0][0] for c in groups.values()}) < 2:
            continue
        proposal = tuple((a - 1, b) if side == "left" else (a, b + 1) for a, b in path)
        yield proposal, "allele_associated_flank", {
            "side": side,
            "internal_repeat_base_counts": dict(sorted(groups.items())),
            "next_column": dict(Counter(_column(source, path, side, 2))),
        }


def _retained_reason(path, other):
    if len(other) == len(path):
        offsets = [(c - a, d - b) for (a, b), (c, d) in zip(path, other)]
        if (len(set(d - c for c, d in other)) == 1
                and max(max(abs(x), abs(y)) for x, y in offsets) <= 3
                and (all(x == 0 for x, _ in offsets) or all(y == 0 for _, y in offsets))):
            return "retained_edge"
        return None
    if len(other) != len(path) + 1:
        return None
    missing = set(other) - set(path)
    if len(missing) == 1:
        interval, = missing
        if other[0] != interval and other[-1] != interval:
            idx = other.index(interval)
            gaps = _gaps(other)
            typical = median(_gaps(other[:idx]) + _gaps(other[idx + 1:]) or gaps)
            if all(abs(g - typical) <= max(2, typical * 0.15) for g in gaps[idx - 1:idx + 1]):
                return "retained_internal_repeat"
        return None
    # A displaced occurrence requires comparing both complete segmentations.
    # Require one displaced old occurrence and multiple shared positions.
    displaced = [old for old in path if not any(
        old[0] < end and start < old[1] for start, end in other)]
    old_gaps, new_gaps = _gaps(path), _gaps(other)
    if (len(displaced) == 1 and len(path) >= 4 and
            any(b <= displaced[0][0] < displaced[0][1] <= c
                for (_, b), (c, _) in zip(other, other[1:])) and
            max(new_gaps) - min(new_gaps) < max(old_gaps) - min(old_gaps)):
        return "retained_resegmentation"
    return None


def compare_paths(source, path, retained):
    """Compare one canonical path with source placements and scored observations.

    'retained' contains (geometry, identifier, category, score) for intersecting
    native observations. Limits bound review volume and source alignment work.
    """
    if not 3 <= len(path) <= 64 or not _valid(path, len(source)):
        return None
    if any(b - a > 128 for a, b in path):
        return None
    repeats = [source[a:b] for a, b in path]
    motif = _modal(repeats)
    proposals = {}
    for generator in (_single_copy(source, path, repeats, motif), _edges(source, path, repeats)):
        for other, reason, evidence in generator:
            if other != path and _valid(other, len(source)):
                proposals.setdefault(other, {})[reason] = evidence
    ranked = defaultdict(dict)
    for other, observation, category, score in retained:
        if other == path or not _valid(other, len(source)):
            continue
        reason = _retained_reason(path, other)
        if reason is None:
            continue
        rank = (len({end - start for start, end in other}), -score, other, observation)
        if other not in ranked[reason] or rank < ranked[reason][other][0]:
            ranked[reason][other] = rank, {"observation_id": observation, "category": category,
                                           "observed_score": score}
    for reason, candidates in ranked.items():
        for rank, evidence in sorted(candidates.values(), key=lambda row: row[0])[:4]:
            proposals.setdefault(rank[2], {})[reason] = evidence
    if not proposals:
        return None
    return {
        "primary_repeat_intervals": path,
        "primary_metrics": _summary(source, path),
        "hypotheses": [
            {"repeat_intervals": other, "evidence": evidence, "metrics": _summary(source, other)}
            for other, evidence in sorted(proposals.items())
        ],
    }
