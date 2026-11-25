"""Select a single divergent terminal copy from a long, regular source array."""

from collections import Counter
from statistics import median_low

from components.refinement.candidates import CrisprConsensus
from components.refinement.enhancement_search import (
    _source_repeat_edit_distance as distance,
)


def _mismatches(first, second):
    return sum(a != b for a, b in zip(first, second))


def _core(source, path):
    if not 8 <= len(path) <= 128:
        return None
    width = path[0][1] - path[0][0]
    if not 20 <= width <= 60 or any(
        b - a != width or a < 0 or b > len(source) for a, b in path
    ):
        return None
    gaps = [c - b for (_, b), (c, _) in zip(path, path[1:])]
    if min(gaps) < 15 or max(gaps) > 80 or max(gaps) - min(gaps) > 4:
        return None
    spacers = [source[b:c] for (_, b), (c, _) in zip(path, path[1:])]
    if 5 * len(set(spacers)) < 4 * len(spacers) or any(
        set(s) - set("ACGT") for s in spacers
    ):
        return None
    repeats = Counter(source[a:b] for a, b in path)
    if any(set(r) - set("ACGT") or len(set(r)) < 3 for r in repeats):
        return None
    # Use each common, intact allele, including ties. Rare errors or a
    # single divergent core copy must not become a template for a new end.
    templates = tuple(sorted(r for r, n in repeats.items() if 4 * n >= len(path)))
    if sum(repeats[t] for t in templates) * 4 < 3 * len(path) or any(
        _mismatches(a, b) > 2 for a in templates for b in repeats
    ):
        return None
    return width, Counter(gaps), spacers, templates


def _end(source, path, core, side):
    width, gaps, spacers, templates = core
    anchor = path[0][0] if side == "left" else path[-1][1]
    typical = median_low(gaps.elements())
    # Search both observed and neighboring spacing, even for the runner-up.
    # Only a full-width window at independently recurring spacing may win.
    windows = []
    for gap in range(min(gaps) - 3, max(gaps) + 4):
        start = anchor - gap - width if side == "left" else anchor + gap
        if start < 0 or start + width > len(source):
            continue
        sequence = source[start : start + width]
        if set(sequence) - set("ACGT"):
            continue
        errors = min(_mismatches(t, sequence) for t in templates)
        windows.append((errors, start, gap))
    if len(windows) < 2:
        return None
    windows.sort()
    errors, start, gap = windows[0]
    if (
        gaps[gap] < 2
        or 10 * errors > 3 * width
        or windows[1][0] - errors < max(3, width // 8)
    ):
        return None
    end = start + width
    hit = source[start:end]
    spacer = source[end:anchor] if side == "left" else source[anchor:start]
    if set(spacer) - set("ACGT") or spacer in spacers:
        return None
    edits = min(distance(t, hit) for t in templates)
    if edits <= 3 or 4 * edits > width:
        return None

    rival = min(
        min(distance(t, source[s : s + width]) for t in templates)
        for _, s, _ in windows[1:]
    )
    if rival <= edits:
        return None
    # A shorter or longer alignment can share an end and fit equally well.
    # In a uniformly full-width core, a tie supports membership but makes
    # the exact boundary provisional. A strictly better variable-width fit
    # precludes selecting the full-width interpretation through this rule.
    variable_rival = width + 2
    tied_endpoints = []
    for g in range(min(gaps) - 3, max(gaps) + 4):
        for size in (width - 2, width - 1, width + 1, width + 2):
            a = anchor - g - size if side == "left" else anchor + g
            if a < 0 or a + size > len(source):
                continue
            seq = source[a : a + size]
            if set(seq) - set("ACGT"):
                continue
            cost = min(distance(t, seq) for t in templates)
            if cost < edits:
                return None
            if cost == edits:
                tied_endpoints.append((a, a + size))
            variable_rival = min(variable_rival, cost)

    # An equally repeat-like window inside a known spacer argues against
    # accepting a distant match purely from its separation from the array.
    background = [
        s[i : i + width] for s in (*spacers, spacer) for i in range(len(s) - width + 1)
    ]
    background_edits = min(
        (distance(t, seq) for t in templates for seq in background), default=None
    )
    if background_edits is not None and background_edits <= edits + 1:
        return None
    evidence = {
        "coordinate_system": "zero_based_half_open",
        "side": side,
        "added_repeat_interval": (start, end),
        "core_repeat_count": len(path),
        "core_width": width,
        "core_allele_counts": {
            t: sum(source[a:b] == t for a, b in path) for t in templates
        },
        "spacer_length": gap,
        "core_median_spacer_length": typical,
        "core_same_length_spacers": gaps[gap],
        "edit_distance": edits,
        "positional_mismatches": errors,
        "next_full_width_edits": rival,
        "next_full_width_mismatches": windows[1][0],
        "best_variable_width_edits": variable_rival,
        "equal_edit_alternative_boundary": bool(tied_endpoints),
        "equal_edit_alternative_intervals": sorted(tied_endpoints),
        "best_within_spacer_edits": background_edits,
    }
    interval = ((start, end),)
    return (interval + path if side == "left" else path + interval), evidence


def select(source, primary, retained):
    """Return one source-backed extension, abstaining if both ends qualify."""
    if not 0.75 <= retained.get(primary, 0) <= 1:
        return None
    core = _core(source, primary)
    if core is None:
        return None
    hits = [
        found
        for side in ("left", "right")
        if (found := _end(source, primary, core, side)) is not None
    ]
    if len(hits) != 1:
        return None
    path, evidence = hits[0]
    return path, "regular_core_divergent_terminal", evidence


def terminal_consensus_basis(alignment, classification):
    """Describe a co-modal core allele when a divergent end tips a core tie.

    The two outer intact copies must support the same complete allele. This
    makes the representative independent of which strand is read first.
    All copies, including the divergent one, remain in the scored alignment.
    """
    if not classification or (
        classification.get("method") != "source_supported_competing_geometry"
        or classification.get("rule") != "regular_core_divergent_terminal"
    ):
        return None
    evidence = classification["evidence"]
    core, terminal = (
        (alignment[1:], alignment[0]) if evidence["side"] == "left"
        else (alignment[:-1], alignment[-1])
    )
    if len(core) != evidence["core_repeat_count"] or core[0] != core[-1]:
        return None
    representative = core[0]
    if set(representative) - set("ACGT "):
        return None
    alleles = Counter(core)
    maximum = max(alleles.values())
    if alleles[representative] != maximum or list(alleles.values()).count(maximum) < 2:
        return None
    raw, aligned = CrisprConsensus(alignment).output()
    if aligned == representative:
        return None
    tied = []
    for index, (reported, scored) in enumerate(zip(representative, aligned)):
        counts = Counter(row[index] for row in core)
        if counts[reported] != max(counts.values()):
            return None
        if reported != scored:
            if (reported not in "ACGT" or scored not in "ACGT"
                    or counts[reported] != counts[scored] or terminal[index] != scored):
                return None
            tied.append({
                "alignment_column_1based": index + 1,
                "core_base_counts": dict(sorted(counts.items())),
                "added_terminal_base": scored,
            })
    return {
        "method": "co_modal_allele_at_both_intact_core_ends",
        "reported_consensus": representative.replace(" ", ""),
        "reported_gapped_consensus": representative,
        "scored_full_alignment_consensus": raw,
        "core_repeat_count": len(core),
        "tied_columns": tied,
    }
