"""Compare source-derived extent proposals against the unmodified accepted core."""

from collections import Counter, defaultdict

from components.evidence.candidate_origins import inherit_candidate_origins
from components.evidence.initial_competition import guard_reconstructed_occurrences
from components.reconstruction.reconstruction_validation import (
    needs_coding_context,
    validate_reconstruction_candidate,
)
from components.refinement.enhancement_search import (
    _source_repeat_edit_distance as distance,
    build_source_backed_candidate_from_boundaries as build_candidate,
)
from components.refinement.periodic_completion import _spacing, plan_path
from components.refinement.potential_degenerate import added_member_evidence, geometry
from components.refinement.repeat_boundary_repair import classification_restriction
from components.rescue.likely_tandem_transaction import validated_row_tags
from components.scoring.evaluation import build_advanced_fuzzy_search_filter


def _span_intersects(path, other):
    return path[0][0] < other[-1][1] and other[0][0] < path[-1][1]


def _contains_anchors(path, anchor):
    # Require an ordered, edge-aligned embedding, stronger than mere overlap.
    i = 0
    for start, end in anchor:
        while i < len(path) and path[i][1] <= start:
            i += 1
        if i == len(path):
            return False
        left, right = path[i]
        if max(abs(left - start), abs(right - end)) > 5:
            return False
        i += 1
    return True


def _usable(row):
    return not classification_restriction(
        row[1]
    ) and "tandem" not in validated_row_tags(row)


def _retained(categories, key, old):
    by_path = {}
    for category in (1, 3, 4):
        for row in categories[category].get(key, ()):
            path = geometry(row[1])
            if (
                not len(old) < len(path) <= 128
                or not _contains_anchors(path, old)
                or not _usable(row)
            ):
                continue
            if path not in by_path or row[0] > by_path[path][0]:
                by_path[path] = row
    return by_path


def _assess(dna, candidates, parameters, score, source_cds):
    if not candidates:
        return []
    rows = score({"extent": candidates})["extent"]
    if len(rows) != len(candidates) or any(
        row[1] is not c for c, row in zip(candidates, rows)
    ):
        raise ValueError("extent scoring changed candidate identity or order")
    structural = build_advanced_fuzzy_search_filter(parameters)
    results = []
    for row in rows:
        cds = tuple(source_cds()) if needs_coding_context(row[1]) else None
        check = validate_reconstruction_candidate(
            dna, row, row, structural, source_cds=cds
        )
        results.append((row, check))
    return results


def _convergent_path(dna, old, donors, ordinary):
    """Require distinct observed extents and a majority-defined, locally bounded phase."""
    support = defaultdict(dict)
    for path, row in donors.items():
        if (
            row[0] < 0.75
            or len(path) < 2 * len(old)
            or len({b - a for a, b in path}) != 1
        ):
            continue
        repeats = [dna[a:b] for a, b in path]
        motif, count = Counter(repeats).most_common(1)[0]
        if (
            3 * count < 2 * len(path)
            or set(motif) - set("ACGT")
            or len(set(dna[b:c] for (_, b), (c, _) in zip(path, path[1:])))
            != len(path) - 1
        ):
            continue
        core = build_candidate(dna, row[1], *zip(*path))
        proposed = plan_path(dna, core, ordinary)
        if (
            not len(path) < len(proposed) <= len(path) + 2
            or proposed[0][1] >= old[0][0]
            or proposed[-1][0] <= old[-1][1]
        ):
            continue
        # Restrict this recovery to simple end extensions. A donor with new
        # internal placements does not corroborate the same ordered phase.
        if not any(
            proposed[i : i + len(path)] == path
            for i in range(len(proposed) - len(path) + 1)
        ):
            continue
        if any(a <= 0 or b >= len(dna) for a, b in path):
            continue
        # A uniform trim can also create many exact copies. Require both
        # neighboring spacer columns to be variable across the observed copies.
        if any(
            3 * max(Counter(dna[pos] for pos in positions).values()) >= 2 * len(path)
            for positions in ([a - 1 for a, _ in path], [b for _, b in path])
        ):
            continue
        support[proposed][path] = row
    qualified = [
        (path, rows)
        for path, rows in support.items()
        if len({len(before) for before in rows}) >= 2
    ]
    return qualified[0] if len(qualified) == 1 else None


def recover_embedded(
    dna, categories, parameters, score, source_cds, reconstructed, deferred
):
    """Replace an isolated ordinary three-copy fragment only on convergent source evidence."""
    anchors = [
        (cat, key, i, row, geometry(row[1]))
        for cat in (0, 2)
        for key, rows in categories[cat].items()
        for i, row in enumerate(rows)
    ]
    updates = []
    established = [
        geometry(row[1])
        for cat in (0, 2)
        for rows in reconstructed[cat].values()
        for row in rows
    ]
    for cat, key, i, row, old in anchors:
        if (
            len(old) != 3
            or not _usable(row)
            or len(set(row[1].list_repeats)) != 1
            or len(row) == 4
            and row[3].get("classification_basis", {}).get("method")
            not in (None, "model_score_band")
        ):
            continue
        donors = _retained(categories, key, old)
        match = _convergent_path(
            dna, old, donors, parameters["param_max_edit_distance"]
        )
        if match is None:
            continue
        path, witnesses = match
        # Earlier stages can promote a longer observation after scoring the
        # fragment. Only override that complete reconstruction for a net gain
        # in distinct occurrences, and keep every previously accepted one.
        existing = [g for g in established if _span_intersects(g, path)]
        if len(existing) != 1 or not _contains_anchors(path, existing[0]):
            continue
        current = existing[0]
        if len(path) <= len(current):
            if path != current and deferred is not None:
                deferred.setdefault(current, {})[path] = {
                    "method": "convergent_boundary_alternative",
                    "donor_repeat_counts": sorted({len(g) for g in witnesses}),
                }
            continue
        if (
            any(
                other is not row and _span_intersects(path, g)
                for _, _, _, other, g in anchors
            )
            or any(_span_intersects(path, previously) for *_, previously in updates)
            or not guard_reconstructed_occurrences(
                path, row, *witnesses.values(), source=dna
            )
        ):
            continue
        candidate = build_candidate(dna, row[1], *zip(*path))
        inherit_candidate_origins(
            candidate, (row[1], *(r[1] for r in witnesses.values()))
        )
        ((fresh, check),) = _assess(dna, [candidate], parameters, score, source_cds)
        if not check.qualified or fresh[0] < 0.75:
            continue
        longest = max(witnesses, key=lambda g: (len(g), witnesses[g][0], g))
        donor = build_candidate(dna, witnesses[longest][1], *zip(*longest))
        observed = {interval for extent in witnesses for interval in extent}
        labels = added_member_evidence(
            dna, donor, longest, path, parameters["param_max_edit_distance"]
        )
        candidate.potential_degenerate_repeats = {
            interval: label
            for interval, label in labels.items()
            if interval not in observed
        }
        candidate.periodic_refinement = {
            "method": "convergent_retained_extents",
            "status": "source_supported_geometry",
            "coordinate_system": "zero_based_half_open",
            "parent_repeat_intervals": [old, *sorted(witnesses)],
            "parent_model_scores": [
                float(row[0]),
                *(float(witnesses[g][0]) for g in sorted(witnesses)),
            ],
        }
        updates.append((cat, key, i, fresh, path))
    if not updates:
        return None
    result = [
        {key: list(rows) for key, rows in category.items()} for category in categories
    ]
    for cat, key, i, fresh, path in sorted(updates, key=lambda u: u[:3], reverse=True):
        prior = result[cat][key].pop(i)
        if not result[cat][key]:
            del result[cat][key]
        result[cat + 1].setdefault(key, []).append(prior)
        result[0].setdefault(key, []).append(fresh)
    return result


def _width_hypothesis(dna, path, motif):
    width = len(motif)
    observed = Counter(dna[start:end] for start, end in path)
    fixed = []
    for i, (start, end) in enumerate(path):
        sites = []
        if 0 < abs(end - start - width) <= 3:
            for left in {start, end - width}:
                right = left + width
                if left < (path[i - 1][1] if i else 0) or right > (
                    path[i + 1][0] if i + 1 < len(path) else len(dna)
                ):
                    continue
                seq = dna[left:right]
                if set(seq) - set("ACGT"):
                    continue
                edits = distance(motif, seq)
                hamming = sum(a != b for a, b in zip(motif, seq))
                # Recurrent clipped copies can share an allele with a divergent
                # edge. Compare both extensions at the modal width in the
                # research evidence, including the cost of newly observed bases.
                clipped = dna[start:end]
                allowance = abs(len(clipped) - width) if observed[clipped] > 1 else 1
                if (
                    edits <= distance(motif, clipped) + allowance
                    and hamming <= width // 3
                ):
                    sites.append(((edits, hamming), left, right))
        sites.sort()
        fixed.append(
            sites[0][1:]
            if sites and (len(sites) == 1 or sites[0][0] < sites[1][0])
            else (start, end)
        )
    return (
        tuple(path)
        if any(b >= c for (_, b), (c, _) in zip(fixed, fixed[1:]))
        else tuple(fixed)
    )


def _terminal_hypotheses(dna, core):
    path = geometry(core)
    motif = core.consensus
    if (
        len(path) < 8
        or not 20 <= len(motif) <= 60
        or set(motif) - set("ACGT")
        or len(set(motif)) < 3
        or any(b - a != len(motif) for a, b in path)
        or max(Counter(core.list_repeats).values()) < 2
    ):
        return
    spacing = _spacing(path)
    if not spacing:
        return
    typical, margin = spacing
    width = len(motif)
    for side, center in (
        ("left", path[0][0] - typical - width),
        ("right", path[-1][1] + typical),
    ):
        rivals = []
        for pos in range(
            max(0, center - margin), min(len(dna) - width, center + margin) + 1
        ):
            seq = dna[pos : pos + width]
            if set(seq) - set("ACGT"):
                continue
            edits = distance(motif, seq)
            hamming = sum(a != b for a, b in zip(motif, seq))
            rivals.append(((edits, hamming, abs(pos - center)), pos))
        rivals.sort()
        if not rivals:
            continue
        (edits, mismatches, deviation), pos = rivals[0]
        if (
            edits > width // 3
            or mismatches > width // 3
            or len(rivals) > 1
            and rivals[1][0] == rivals[0][0]
        ):
            continue
        interval = ((pos, pos + width),)
        yield (
            (interval + path if side == "left" else path + interval),
            {
                "method": "single_terminal_window",
                "side": side,
                "edit_distance": edits,
                "positional_mismatches": mismatches,
                "spacer_deviation": deviation,
                "rival_rank": list(rivals[1][0]) if len(rivals) > 1 else None,
            },
        )


def annotate_trials(dna, categories, parameters, score, source_cds, deferred=None):
    """Attach observational alternatives to accepted anchors for the optional research sidecar."""
    anchors = [
        (key, row, geometry(row[1]))
        for cat in (0, 2)
        for key, rows in categories[cat].items()
        for row in rows
    ]
    for key, row, old in anchors:
        if not _usable(row) or len(old) > 128:
            continue
        donors = _retained(categories, key, old)
        hypotheses = dict(deferred.get(old, {})) if deferred else {}
        # Ranking limits research cost in highly repetitive neighborhoods.
        # The accepted core and all ordinary category decisions are unchanged.
        ranked = sorted(donors, key=lambda g: (-donors[g][0], -len(g), g))[:12]
        for path in ranked:
            donor = donors[path]
            core = build_candidate(dna, donor[1], *zip(*path))
            if max(Counter(core.list_repeats).values()) < 2:
                continue
            proposed = plan_path(dna, core, parameters["param_max_edit_distance"])
            if len(proposed) > len(path):
                hypotheses.setdefault(
                    proposed,
                    {"method": "retained_extent", "donor_repeat_intervals": path},
                )
            phase = _width_hypothesis(dna, proposed, core.consensus)
            if phase != proposed:
                hypotheses.setdefault(
                    phase,
                    {"method": "anchor_edge_width", "donor_repeat_intervals": path},
                )
        for path, details in _terminal_hypotheses(dna, row[1]):
            hypotheses.setdefault(path, details)
        paths = [
            g
            for g in hypotheses
            if g != old
            and _contains_anchors(g, old)
            and not any(
                other is not row and _span_intersects(g, q) for _, other, q in anchors
            )
            and guard_reconstructed_occurrences(g, row, source=dna)
        ]
        paths.sort()
        candidates = [build_candidate(dna, row[1], *zip(*g)) for g in paths]
        records = []
        for g, (fresh, check) in zip(
            paths, _assess(dna, candidates, parameters, score, source_cds)
        ):
            records.append(
                dict(
                    hypotheses[g],
                    repeat_intervals=g,
                    model_score=float(fresh[0]),
                    consensus=fresh[1].consensus,
                    qualified=check.qualified,
                    decision=check.evidence["decision"],
                )
            )
        if records:
            row[1].extent_trials = {
                "status": "unselected",
                "coordinate_system": "zero_based_half_open",
                "unchanged_repeat_intervals": old,
                "proposals": records,
            }
