"""Resolve repeat intervals from an accepted core and its local spacer spacing."""

from collections import Counter
from functools import lru_cache
from statistics import median

from components.paths import TOOLS_ROOT
from components.refinement.enhancement_search import (
    _source_repeat_edit_distance as distance,
)
from components.refinement.potential_degenerate import geometry
from components.refinement.repeat_boundary_repair import classification_restriction
from components.sequence.utils import reverse_complement


@lru_cache(maxsize=1)
def _reference_repeats():
    panel = TOOLS_ROOT / "blasting" / "Verified_repeats_dataset1.fa"
    repeats = {
        record.partition("\n")[2].replace("\n", "").upper()
        for record in panel.read_text().split(">")
        if record
    }
    return repeats | {reverse_complement(s) for s in repeats}


def _reference_phase(dna, motif, path):
    """Only restore missing reference columns independently corroborated by the core."""
    if len(path) < 4 or motif in _reference_repeats():
        return motif, path
    phases = set()
    for ref in _reference_repeats():
        extra = len(ref) - len(motif)
        if not 1 <= extra <= 6:
            continue
        for left in range(extra + 1):
            if ref[left : left + len(motif)] != motif:
                continue
            right = extra - left
            # The unmodified edge anchors the full reference width even
            # when a terminal source copy was previously clipped.
            if left and right:
                projected = tuple((a - left, b + right) for a, b in path)
            elif left:
                projected = tuple((b - len(ref), b) for a, b in path)
            else:
                projected = tuple((a, a + len(ref)) for a, b in path)
            if any(a < 0 or b > len(dna) or b - a != len(ref) for a, b in projected):
                continue
            edges = (*range(left), *range(len(ref) - right, len(ref)))
            if all(
                sum(dna[a + k] == ref[k] for a, _ in projected) > len(path) / 2
                for k in edges
            ):
                phases.add((ref, projected))
    return phases.pop() if len(phases) == 1 else (motif, path)


def _spacing(path):
    gaps = [b[0] - a[1] for a, b in zip(path, path[1:])]
    ordinary = sorted(g for g in gaps if 15 <= g <= 80)
    if len(ordinary) < 2:
        return None
    typical = ordinary[(len(ordinary) - 1) // 2]
    mad = median(abs(g - typical) for g in ordinary)
    if mad > 3 or (len(path) == 3 and max(ordinary) - min(ordinary) > 2):
        return None
    return typical, min(9, max(4, typical // 7, round(3 * mad + 2)))


def _fit_width(dna, motif, path, *, include_ends=False):
    """Prefer the modal width at irregular internal copies if it explains the bases."""
    fixed = list(path)
    width = len(motif)
    for n, (a, b) in enumerate(path):
        if abs(b - a - width) not in (1, 2) or (
            not include_ends and n in (0, len(path) - 1)
        ):
            continue
        old_distance = distance(motif, dna[a:b])
        options = []
        starts = {a, b - width}
        if n == 0:
            starts = {b - width}
        elif n == len(path) - 1:
            starts = {a}
        for x in starts:
            y = x + width
            if x < (path[n - 1][1] if n else 0) or y > (
                path[n + 1][0] if n + 1 < len(path) else len(dna)
            ):
                continue
            mismatches = sum(
                base != expected for base, expected in zip(dna[x:y], motif)
            )
            # With multiple independently filled gaps, the inner endpoint
            # anchors the terminal repeat. Allow a one-base displacement
            # when its edit alignment is as good as the clipped copy.
            admissible = mismatches <= old_distance
            if include_ends and n in (0, len(path) - 1):
                admissible |= (
                    distance(motif, dna[x:y]) <= old_distance + 1
                    and mismatches <= old_distance + 2
                )
            if admissible:
                options.append((mismatches, x, y))
        options.sort()
        if options and (len(options) == 1 or options[0][0] < options[1][0]):
            fixed[n] = options[0][1:]
    if any(b > c for (_, b), (c, _) in zip(fixed, fixed[1:])):
        return tuple(path)
    return tuple(fixed)


def recurrent_allele_trim(dna, old, candidate, dominance):
    """Check whether independent source copies support a polymorphic repeat edge."""
    path = geometry(candidate)
    repeats = candidate.list_repeats
    if (len(old) != len(path) or len(set(map(len, repeats))) != 1
            or list(candidate.list_gaped_repeats) != list(repeats)):
        return None
    changed = [i for i, (before, after) in enumerate(zip(old, path)) if before != after]
    if not changed or len(changed) * 4 > len(path):
        return None
    # A single-base trim at the same end must explain every abnormal width.
    shifts = {(path[i][0] - old[i][0], old[i][1] - path[i][1]) for i in changed}
    if shifts not in ({(1, 0)}, {(0, 1)}):
        return None
    edge = 0 if shifts == {(1, 0)} else len(repeats[0]) - 1
    if edge == 0 and 0 in changed or edge != 0 and len(path) - 1 in changed:
        return None  # Only transfer a base into an existing spacer.
    width = len(repeats[0])
    weak = [j for j in (*range(3), *range(width - 3, width))
            if max(Counter(s[j] for s in repeats).values()) < dominance * len(path)]
    if weak != [edge]:
        return None
    support = Counter(dna[a:b] for (a, b), target in zip(old, path) if (a, b) == target)
    # Each corrected copy must match a whole, independently recurring allele.
    # Both variants of the disputed base must already recur in intact copies.
    bases = {s[edge] for s in repeats}
    if (len(bases) != 2 or any(support[repeats[i]] < 2 for i in changed)
            or {s[edge] for s, n in support.items() if n >= 2} != bases):
        return None
    return {"trimmed_side": "left" if edge == 0 else "right",
            "changed_copy_indexes_0based": changed,
            "unchanged_full_allele_counts": dict(sorted(support.items()))}


def _placement(dna, motif, left, right, typical, margin, ordinary):
    width = len(motif)
    lo, hi = (0, len(dna) - width)
    if left is not None:
        lo, hi = max(lo, left + typical - margin), min(hi, left + typical + margin)
    if right is not None:
        lo, hi = (
            max(lo, right - typical - margin - width),
            min(hi, right - typical + margin - width),
        )
    hits = []
    for x in range(lo, hi + 1):
        seq = dna[x : x + width]
        if set(seq) - set("ACGT"):
            continue
        edits = distance(motif, seq)
        if edits > min(width // 3, round(ordinary * width / 20)):
            continue
        hamming = sum(a != b for a, b in zip(motif, seq))
        # A placement that only fits after extensive gapping does not
        # establish both edges of this fixed-width copy.
        if hamming > width // 3:
            continue
        # A well-conserved block plus an appropriate spacer can anchor a
        # divergent terminal tail. Weak similarities alone cannot seed a chain.
        if edits > width // 6:
            run = longest = 0
            for observed, expected in zip(seq, motif):
                run = run + 1 if observed == expected else 0
                longest = max(longest, run)
            if longest < max(9, width // 3):
                continue
        deviation = abs(x - left - typical) if left is not None else 0
        if right is not None:
            deviation += abs(right - x - width - typical)
        hits.append(
            ((edits + min(hamming, edits + 5), deviation, hamming), x, x + width, edits)
        )
    hits.sort()
    if not hits or (len(hits) > 1 and hits[0][0] == hits[1][0]):
        return None
    return hits[0][1:]


def plan_path(dna, core, ordinary):
    """Propose a bounded repeat path using the core and the forward source only."""
    original = geometry(core)
    motif = core.consensus
    if (
        len(original) < 3
        or not 20 <= len(motif) <= 60
        or set(motif) - set("ACGT")
        or len(set(motif)) < 3
        or classification_restriction(core)
    ):
        return original
    motif, phased = _reference_phase(dna, motif, original)
    if any(b >= c for (_, b), (c, _) in zip(phased, phased[1:])):
        return original
    path = _fit_width(dna, motif, phased)
    spacing = _spacing(path)
    if spacing is None:
        return path
    typical, margin = spacing
    if max(Counter(dna[a:b] for a, b in path).values()) < 2:
        return path
    extended = [path[0]]
    for previous, following in zip(path, path[1:]):
        hit = _placement(
            dna, motif, previous[1], following[0], typical, margin, ordinary
        )
        if hit:
            extended.append(hit[:2])
        extended.append(following)
    internal = len(extended) - len(path)
    for side in ("left", "right"):
        for added_on_side in range(64):
            left = extended[-1][1] if side == "right" else None
            right = extended[0][0] if side == "left" else None
            hit = _placement(dna, motif, left, right, typical, margin, ordinary)
            if not hit:
                break
            x, y, edits = hit
            if len(path) == 3 and added_on_side == 0 and edits > len(motif) // 6:
                break
            if side == "left":
                extended.insert(0, (x, y))
            else:
                extended.append((x, y))
            if edits > len(motif) // 6:
                break
        else:
            # A bounded search which reaches its cap has not localized the end.
            return original
    if internal >= 2:
        extended = _fit_width(dna, motif, extended, include_ends=True)
    return tuple(extended)
