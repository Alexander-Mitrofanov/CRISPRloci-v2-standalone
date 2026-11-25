"""Resolve retained source placements using local spacing and repeat alleles.

Inputs are forward, zero-based, half-open intervals from the same source. This
selector only proposes a competing path. The transaction owns source validation,
fresh scoring, coding/tandem checks, occurrence authority and collision checks.
"""

from collections import Counter
from statistics import median

from components.evidence.review_hypotheses import _gaps, _valid
from components.refinement.enhancement_search import _source_repeat_edit_distance as distance
from components.refinement.extent_trials import _width_hypothesis
from components.refinement.periodic_completion import _placement


def _motif(dna, path):
    """Require one modal width and a recurrent, unambiguous allele at that width."""
    repeats = [dna[a:b] for a, b in path]
    widths = Counter(map(len, repeats)).most_common()
    width, count = widths[0]
    if (not 20 <= width <= 60 or count * 2 <= len(path)
            or any(set(s) - set("ACGT") for s in repeats)):
        return None
    alleles = Counter(s for s in repeats if len(s) == width).most_common()
    motif, count = alleles[0]
    if (count < 2 or len(alleles) > 1 and count == alleles[1][1]
            or len(set(motif)) < 3):
        return None
    return motif


def _regular(gaps):
    return bool(gaps and 15 <= min(gaps) and max(gaps) <= 80
                and max(gaps) - min(gaps) <= 2)


def _distinct_spacers(dna, path):
    gaps = [dna[b:c] for (_, b), (c, _) in zip(path, path[1:])]
    return (len(gaps) == len(set(gaps))
            and all(len(set(s)) >= 3 and set(s) <= set("ACGT") for s in gaps))


def _internal(dna, old, other, motif):
    if (len(other) != len(old) + 1 or not set(old) < set(other)
            or other[0] != old[0] or other[-1] != old[-1]):
        return None
    (copy,) = set(other) - set(old)
    a, b = copy
    if b - a != len(motif) or any(y - x != len(motif) for x, y in old):
        return None
    gaps = _gaps(other)
    idx = other.index(copy)
    if not _regular(gaps) or not _distinct_spacers(dna, other):
        return None
    left, right = other[idx - 1][1], other[idx + 1][0]
    typical = median(gaps[:idx - 1] + gaps[idx + 1:])
    margin = max(3, round(typical * 0.15))
    width = len(motif)
    scored = sorted((distance(motif, dna[x:x + width]), x) for x in range(
        max(left, int(left + typical - margin), int(right - width - typical - margin)),
        min(right - width, int(left + typical + margin),
            int(right - width - typical + margin)) + 1))
    # Insist on a placement that wins the sequence search independently of
    # its ideal spacer lengths. A source-tied edge remains an alternative.
    if (len(scored) < 2 or scored[0][1] != a or scored[1][0] < scored[0][0] + 2
            or scored[0][0] * 5 > width):
        return None
    return other, "periodic_internal_insertion", {
        "split_spacer_interval": [left, right], "inserted_interval": list(copy),
        "spacer_lengths": gaps, "motif": motif, "edit_distance": scored[0][0],
        "next_placement_distance": scored[1][0], "placements_compared": len(scored),
        "retained_donor": [list(p) for p in other],
    }


def _width_without_extra_edits(dna, donor, motif):
    # An extension that covers an internal deletion by importing spacer bases
    # pays for both an insertion and the original deletion. Retain the shorter
    # copy unless a full-width, substitution-only explanation costs no more.
    corrected = []
    for old, new in zip(donor, _width_hypothesis(dna, donor, motif)):
        before = dna[slice(*old)]
        after = dna[slice(*new)]
        mismatches = sum(a != b for a, b in zip(after, motif))
        if old != new and mismatches <= distance(before, motif):
            corrected.append(new)
        else:
            corrected.append(old)
    return tuple(corrected)


def _extend(dna, old, donor, motif):
    """Use a retained terminal occurrence, preserving supported internal indels."""
    if len(donor) != len(old) + 1:
        return None
    left = donor[1:] == old
    if not left and donor[:-1] != old:
        return None
    old_gaps = _gaps(old)
    if not _regular(old_gaps):
        return None
    width = len(motif)
    fixed = _width_without_extra_edits(dna, donor, motif)
    if not _regular(_gaps(fixed)) or not _distinct_spacers(dna, fixed):
        return None
    comparison = {}
    if len(old) == 3:
        # Two observed clipped copies must independently support the SAME
        # anchored edge and the source alleles must share an internal variant.
        changed = [(before, after) for before, after in zip(donor, fixed) if before != after]
        edge = 0 if left else -1
        if (len(changed) != 2 or donor[edge] == fixed[edge]
                or any(b - a != width for a, b in fixed)
                or not (all(p[0] == q[0] for p, q in changed)
                        or all(p[1] == q[1] for p, q in changed))):
            return None
        variants = [dna[slice(*q)] for _, q in changed]
        if not any(a == b != m for a, b, m in zip(*variants, motif)):
            return None
        result = fixed
        rule = "retained_substitution_width_extension"
    else:
        # A recurrent shorter allele whose gaps are internal is a positive
        # boundary witness. Keep it intact when extending the retained seed.
        counts = Counter(dna[a:b] for a, b in old if b - a < width)
        if not any(n >= 2 and motif.startswith(s[:5]) and motif.endswith(s[-5:])
                   and distance(motif, s) == width - len(s) for s, n in counts.items()):
            return None
        typical = round(median(old_gaps))
        result = list(fixed)
        projected_sites = {}
        for _ in range(8):
            hit = _placement(
                dna, motif, None if left else result[-1][1],
                result[0][0] if left else None, typical, max(4, typical // 7),
                width // 4,
            )
            if not hit:
                break
            a, b, edits = hit
            full_width = a, b
            # An imported outer base can hide an indel. Prefer a one-base
            # trim only when BOTH motif fit and the adjacent spacer improve.
            options = [(a + 1, b), (a, b - 1)]
            preferred = [
                (x, y) for x, y in options
                if distance(motif, dna[x:y]) < edits
                and abs((result[0][0] - y if left else x - result[-1][1]) - typical)
                < abs((result[0][0] - b if left else a - result[-1][1]) - typical)
            ]
            if len(preferred) == 1:
                a, b = preferred[0]
            projected_sites[a, b] = full_width
            if left:
                result.insert(0, (a, b))
            else:
                result.append((a, b))
            if edits > width // 6:
                break
        else:
            return None
        result = tuple(result)
        if len(result) <= len(donor) or not _regular(_gaps(result)):
            return None
        # An equally long, equal-width interpretation uses the full searched
        # terminal windows and corrects clipped edges in one step. Keep its
        # actual boundaries for independent scoring alongside the indel path.
        windows = tuple(projected_sites.get(p, p) for p in result)
        uniform = _width_hypothesis(dna, windows, motif)
        if (uniform != result and _valid(uniform, len(dna))
                and all(b - a == width for a, b in uniform)
                and all(set(dna[a:b]) <= set("ACGT") for a, b in uniform)
                and _distinct_spacers(dna, uniform)):
            comparison = {
                "alternative_repeat_intervals": [list(p) for p in uniform],
                "alternative_metrics": {
                    "interpretation": "uniform_repeat_width",
                    "spacer_lengths": _gaps(uniform),
                    "repeat_edit_distances": [distance(motif, dna[a:b]) for a, b in uniform],
                    "different_copy_indexes_1based": [
                        i for i, (a, b) in enumerate(zip(result, uniform), 1) if a != b
                    ],
                    "endpoint_status": "source_supported_alternative",
                },
            }
        rule = "retained_indel_allele_extension"
    if not _distinct_spacers(dna, result):
        return None
    return result, rule, {
        "motif": motif, "retained_donor": [list(p) for p in donor],
        "spacer_lengths": _gaps(result),
        "repeat_edit_distances": [distance(motif, dna[a:b]) for a, b in result],
        **comparison,
    }


def _duplication(dna, old, other, retained):
    """A partial tandem duplication can mimic a repeat with a deleted spacer."""
    if (len(other) != len(old) + 1 or len(old) < 5 or retained[other] < .75
            or len({b - a for a, b in other}) != 1):
        return None
    motif = _motif(dna, other)
    if not motif or sum(dna[a:b] == motif for a, b in other) * 2 <= len(other):
        return None
    gaps = _gaps(other)
    if min(gaps) < 15 or not _distinct_spacers(dna, other):
        return None
    typical = round(median(gaps))
    idx = gaps.index(max(gaps))
    if not _regular(gaps[:idx] + gaps[idx + 1:]):
        return None
    gap = other[idx][1], other[idx + 1][0]
    missing = [i for i, (a, b) in enumerate(old) if gap[0] <= a < b <= gap[1]]
    if len(missing) != 1:
        return None
    lost = missing[0]
    a, b = old[lost]
    common = []
    for i, (s, e) in enumerate(old):
        if i == lost:
            continue
        mapped = [(x, y) for x, y in other if x <= s < e <= y]
        if len(mapped) != 1:
            return None
        x, y = mapped[0]
        common.append((s - x, y - e))
    if len(set(common)) != 1 or not 2 <= sum(common[0]) <= len(motif) // 4:
        return None
    # Both edges of a longer phase need a consistent source explanation.
    # All shared observations extend the very same base columns on one edge.
    if min(common[0]) != 0:
        return None
    neighbor = [(s, e) for j, (s, e) in enumerate(old) if abs(j - lost) == 1
                and min(abs(e - a), abs(s - b)) <= 2 and dna[s:e] == dna[a:b]]
    if len(neighbor) != 1:
        return None
    n, _ = neighbor[0]
    period = abs(a - n)
    corrected = gaps.copy()
    corrected[idx] -= period
    if not _regular(corrected):
        return None
    left, right = min(a, n), min(a, n) + b - a
    # Bound the observation to this candidate interval; long periodic
    # sequences outside the locus cannot add support for this decision.
    while left > other[0][0] and dna[left - 1] == dna[left + period - 1]:
        left -= 1
    while right + period < other[-1][1] and dna[right] == dna[right + period]:
        right += 1
    if right - left < period + 4:
        return None
    # Demand an actual alternative which labels the duplicated segment as a
    # further copy, sharing every other unmodified boundary. Interpret the
    # missing spacer explicitly instead of treating an extra copy as a win.
    rivals = [q for q in retained if len(q) == len(other) + 1
              and set(other) < set(q) and retained[q] >= .75
              and 0 <= min(_gaps(q)) <= 2
              and any(x <= a and y >= b for x, y in set(q) - set(other))]
    errors = [distance(motif, dna[x:y]) for x, y in other]
    if max(errors) > len(motif) // 6:
        return None
    detail = {
        "motif": motif, "displaced_interval": [a, b], "duplication_period": period,
        "exact_match_length": right - left,
        "matching_intervals": [[left, right], [left + period, right + period]],
        "shared_repeat_extension": list(common[0]),
        "spacer_lengths": gaps, "duplication_adjusted_spacer_lengths": corrected,
        "primary_spacer_lengths": _gaps(old), "repeat_edit_distances": errors,
        "retained_donor": [list(p) for p in other],
        "alternative_with_short_spacer": [list(p) for p in max(
            rivals, key=lambda q: (retained[q], q))] if rivals else None,
    }
    rank = sum(abs(g - typical) for g in corrected), sum(errors)
    return rank, other, "tandem_duplication_resegmentation", detail


def rank_extension(proposal, *, uniform_eligible, compact_score, uniform_score):
    """Rank two qualified placements of the same extended repeat occurrences.

    The search starts from recurrent indel-bearing alleles. If the source
    also supports a unique modal-width projection with regular spacers, use
    that representation and retain the compact alleles. Native array scores
    describe each interpretation and can depend on the input orientation.
    """
    compact, rule, evidence = proposal
    if (rule != "retained_indel_allele_extension" or not uniform_eligible
            or not evidence.get("alternative_repeat_intervals")):
        return proposal
    metrics = evidence["alternative_metrics"]
    if not _regular(metrics["spacer_lengths"]):
        return proposal
    ranking = {
        "criterion": "qualified_modal_width_with_regular_spacers",
        "compact_score": float(compact_score),
        "uniform_score": float(uniform_score),
        "selected_interpretation": "uniform_repeat_width",
    }
    details = dict(evidence, boundary_ranking=ranking)
    uniform = tuple(map(tuple, evidence["alternative_repeat_intervals"]))
    details.update(
        spacer_lengths=metrics["spacer_lengths"],
        repeat_edit_distances=metrics["repeat_edit_distances"],
        alternative_repeat_intervals=[list(p) for p in compact],
        alternative_metrics=dict(
            metrics,
            interpretation="retained_indel_allele",
            spacer_lengths=evidence["spacer_lengths"],
            repeat_edit_distances=evidence["repeat_edit_distances"],
        ),
    )
    return uniform, rule, details


def select(source, primary, retained):
    """Return one unambiguous source-derived selection or abstain."""
    if not 3 <= len(primary) <= 64 or not _valid(primary, len(source)):
        return None
    motif = _motif(source, primary)
    if motif is None:
        return None
    old_gaps = _gaps(primary)
    regular = _regular(old_gaps)
    tiny_gap = min(old_gaps) <= 2 and max(old_gaps) > 2 * median(old_gaps)
    internal = max(old_gaps) >= 2 * min(old_gaps)
    proposals, resegmentations = {}, []
    for donor, score in sorted(retained.items()):
        if len(donor) != len(primary) + 1 or not _valid(donor, len(source)):
            continue
        if internal and score >= .4:
            choice = _internal(source, primary, donor, motif)
            if choice:
                proposals[choice[0]] = choice
        if regular and score >= .4:
            choice = _extend(source, primary, donor, motif)
            if choice:
                proposals[choice[0]] = choice
        if tiny_gap:
            resolved = _duplication(source, primary, donor, retained)
            if resolved:
                resegmentations.append(resolved)
    if resegmentations:
        resegmentations.sort(key=lambda r: r[0])
        if len(resegmentations) > 1 and resegmentations[0][0] == resegmentations[1][0]:
            return None
        _, path, rule, evidence = resegmentations[0]
        if evidence["alternative_with_short_spacer"] is not None:
            proposals[path] = path, rule, evidence
    return next(iter(proposals.values())) if len(proposals) == 1 else None
