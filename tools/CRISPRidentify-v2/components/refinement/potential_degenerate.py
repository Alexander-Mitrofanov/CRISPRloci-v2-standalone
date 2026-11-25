"""Search one spacer away from confirmed copies for divergent repeat members."""

from statistics import median

from components.refinement.enhancement_search import (
    _minimal_edit_changes, _source_repeat_edit_distance,
    build_source_backed_candidate_from_boundaries,
)
from components.refinement.repeat_boundary_repair import classification_restriction


def geometry(candidate):
    return tuple((s, s + len(r)) for s, r in zip(
        candidate.list_repeat_starts, candidate.list_repeats))


def added_member_evidence(dna, core, old, new, ordinary):
    """Describe divergent copies added to a source-backed periodic path."""
    labels = {}
    motif = core.consensus
    j = 0
    for a, b in new:
        while j < len(old) and old[j][1] <= a:
            j += 1
        if j < len(old) and old[j][0] < b:
            continue
        sub, ins, delete = (len(v) for v in _minimal_edit_changes(motif, dna[a:b]))
        edits = sub + ins + delete
        if edits <= ordinary and max(sub, ins, delete) <= (ordinary + 1) // 2:
            continue
        location = "internal"
        if b <= old[0][0]:
            location = "left"
        elif a >= old[-1][1]:
            location = "right"
        labels[a, b] = dict(
            status="potential_degenerate_repeat",
            core_consensus=motif,
            edit_distance=edits,
            substitutions=sub,
            insertions=ins,
            deletions=delete,
            ordinary_edit_limit=ordinary,
            location=location,
            placement_method="local_spacer_period",
        )
    return labels


def _placement(dna, motif, left, right, bounds, typical, margin, budget, ordinary):
    """Rank all admissible full-length alignments, including short indels."""
    lo, hi = typical - margin, typical + margin
    length = len(motif)
    ranked = []
    for size in range(length - 2, length + 3):
        start = max(bounds[0], left + lo if left is not None else 0)
        stop = min(bounds[1] - size, left + hi if left is not None else bounds[1])
        if right is not None:
            start = max(start, right - hi - size)
            stop = min(stop, right - lo - size)
        for a in range(start, stop + 1):
            end = a + size
            sequence = dna[a:end]
            if set(sequence) - set("ACGT"):
                continue
            distance = _source_repeat_edit_distance(motif, sequence)
            if distance > budget:
                continue
            deviation = (abs(a - left - typical) if left is not None else 0)
            if right is not None:
                deviation += abs(right - end - typical)
            ranked.append(((distance, abs(size - length), deviation), a, end))
    if not ranked:
        return None
    ranked.sort()
    rank, a, end = ranked[0]
    # Multiple equally supported boundaries, or another nearly as good copy
    # at a distinct position, do not determine an individual interval.
    for rival, first, last in ranked[1:]:
        if rival[0] > rank[0] + 1:
            break
        if rival == rank or first >= end or last <= a:
            return None
    s, i, d = (len(changes) for changes in _minimal_edit_changes(motif, dna[a:end]))
    if s + i + d <= ordinary and max(s, i, d) <= (ordinary + 1) // 2:
        return None
    evidence = dict(status="potential_degenerate_repeat", core_consensus=motif,
                    edit_distance=s+i+d, substitutions=s, insertions=i, deletions=d,
                    ordinary_edit_limit=ordinary, relaxed_edit_limit=budget,
                    location="internal" if left is not None and right is not None
                    else "right" if left is not None else "left")
    return (a, end), evidence


def propose(dna, core, bounds, ordinary):
    path = geometry(core)
    motif = core.consensus
    if (len(path) < 3 or not 20 <= len(motif) <= 60
            or set(motif) - set("ACGT") or len(set(motif)) < 3
            or max(motif.count(b) for b in "ACGT") > 0.65 * len(motif)
            or classification_restriction(core)):
        return {}
    if sum(_source_repeat_edit_distance(motif, seq) <= len(motif) // 10
           for seq in core.list_repeats) < 2:
        return {}
    lengths = [b[0] - a[1] for a, b in zip(path, path[1:])]
    typical_gaps = [gap for gap in lengths if 12 <= gap <= 90]
    if not typical_gaps:
        return {}
    # The lower median prevents one missing internal repeat in a short array
    # from turning the doubled spacer into the spacing reference.
    typical = sorted(typical_gaps)[(len(typical_gaps) - 1) // 2]
    mad = median(abs(gap - typical) for gap in typical_gaps)
    margin = min(12, max(4, typical // 5, round(2 * mad)))
    budget = min(ordinary + 4, len(motif) // 3)
    windows = [(None, path[0][0]), (path[-1][1], None)]
    for (left, right), gap in zip(zip(path, path[1:]), lengths):
        if 2 * (typical - margin) + len(motif) - 2 <= gap <= 2 * (typical + margin) + len(motif) + 2:
            windows.append((left[1], right[0]))
    members = {}
    for left, right in windows:
        placement = _placement(dna, motif, left, right, bounds, typical, margin, budget, ordinary)
        if placement:
            interval, evidence = placement
            members[interval] = evidence
    return members


def complete(dna, categories, scored, ordinary, score, categorize):
    """Replace winning cores only if the expanded paths survive fresh competition."""
    accepted = sorted(
        (geometry(row[1])[0][0], geometry(row[1])[-1][1], key, id(row[1]), row[1])
        for index in (0, 2) for key, rows in categories[index].items() for row in rows
    )
    requests = {}
    proposals = {}
    for idx, (start, end, key, token, core) in enumerate(accepted):
        prev_end = max((x[1] for x in accepted[:idx]), default=0)
        next_start = accepted[idx + 1][0] if idx + 1 < len(accepted) else len(dna)
        if prev_end > start or end > next_start:
            continue
        bounds = ((start + prev_end + 1) // 2 if idx else 0,
                  (end + next_start) // 2 if idx + 1 < len(accepted) else len(dna))
        members = propose(dna, core, bounds, ordinary)
        if members:
            path = geometry(core)
            starts, ends = zip(*sorted((*path, *members)))
            extended = build_source_backed_candidate_from_boundaries(dna, core, starts, ends)
            requests.setdefault(key, []).append(extended)
            proposals[token] = (extended, path, members)

    if not requests:
        return categories, scored, []
    fresh = score(requests)
    if fresh.keys() != requests.keys():
        raise ValueError("degenerate repeat scoring changed locus keys")
    updated = dict(scored)
    result = [dict(category) for category in categories]
    committed = []
    for key, variants in requests.items():
        if len(fresh[key]) != len(variants) or any(r[1] is not v for r, v in zip(fresh[key], variants)):
            raise ValueError("degenerate repeat scoring changed candidate identity")
        fresh_rows = {id(row[1]): row for row in fresh[key]}
        replacements = {token: fresh_rows[id(v)] for token, (v, _, _) in proposals.items()
                        if id(v) in fresh_rows}
        trials = [replacements.get(id(row[1]), row) for row in scored[key]]
        tested = categorize({key: trials})
        expected = {id(replacements[t][1]) if t in replacements else t
                    for _, _, k, t, _ in accepted if k == key}
        winners = {id(r[1]) for index in (0, 2) for r in tested[index].get(key, [])}
        if expected != winners:
            continue
        updated[key] = trials
        for index, category in enumerate(tested):
            if key in category:
                result[index][key] = category[key]
            else:
                result[index].pop(key, None)
        committed.extend((path, members) for t, (_, path, members) in proposals.items()
                         if t in replacements)
    return result, updated, committed


def annotate_members(categories, committed):
    """Carry labels through downstream reconstruction only on source-matching copies."""
    if not committed:
        return
    for index in (0, 2):
        for rows in categories[index].values():
            for row in rows:
                candidate = row[1]
                path = set(geometry(candidate))
                labels = dict(getattr(candidate, "potential_degenerate_repeats", {}))
                for anchors, members in committed:
                    if len(path.intersection(anchors)) >= 2:
                        labels.update((loc, evidence) for loc, evidence in members.items() if loc in path)
                if labels:
                    candidate.potential_degenerate_repeats = labels
