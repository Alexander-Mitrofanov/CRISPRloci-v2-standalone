"""Local, same-membership boundary hypotheses supported by source copies.

Coordinates are source-forward, zero-based, half-open. A proposed geometry
still requires fresh scoring, structural validation and occurrence protection.
"""

from collections import Counter
from statistics import median

from components.evidence.review_hypotheses import _column, _edges, _valid
from components.refinement.enhancement_search import (
    build_source_backed_candidate_from_boundaries,
)
from components.sequence.utils import reverse_complement


def _nearby_exemplars(source, path, motif, variable_side=None):
    """Find a separate, locally spaced cluster of three exact copies."""
    typical = median(c - b for (_, b), (c, _) in zip(path, path[1:]))
    for lo, hi in (
        (max(0, path[0][0] - 10000), path[0][0]),
        (path[-1][1], min(len(source), path[-1][1] + 10000)),
    ):
        cluster = []
        pos = source.find(motif, lo, hi)
        while pos != -1:
            if cluster:
                gap = pos - cluster[-1] - len(motif)
                if not typical / 2 <= gap <= 2 * typical:
                    cluster = []
            cluster.append(pos)
            if len(cluster) >= 3:
                x, y, z = cluster[-3:]
                diverse = source[x + len(motif) : y] != source[y + len(motif) : z]
                # All shorter alleles are substrings of the longer copy. A
                # varying outside base distinguishes their source footprints.
                witnesses = tuple((s, s + len(motif)) for s in (x, y, z))
                if diverse and (
                    variable_side is None
                    or len(set(_column(source, witnesses, variable_side, 1))) > 1
                ):
                    return [x, y, z]
            pos = source.find(motif, pos + len(motif), hi)
    return []


def _conserved_flank(source, primary, path, side, columns):
    depth = len(columns) - 1
    if depth < 4:
        return None
    # Actual spacers, excluding the exterior genome, must agree throughout the
    # added block. At the next base the agreement must end.
    internal = path[1:] if side == "left" else path[:-1]
    inner_primary = primary[1:] if side == "left" else primary[:-1]
    for offset in range(1, depth + 1):
        if len(set(_column(source, inner_primary, side, offset))) != 1:
            return None
    if len(set(_column(source, internal, side, 1))) < 2:
        return None
    repeats = Counter(source[a:b] for a, b in path)
    exemplars = []
    for motif, count in repeats.items():
        if count >= 2:
            hits = _nearby_exemplars(source, path, motif)
            if hits:
                exemplars.append(hits)
    if exemplars:
        return {
            "side": side,
            "conserved_spacer_bases": depth,
            "additional_exact_copy_starts_0based": sorted(exemplars),
        }
    return None


def _allele_flanks(source, path, retained, side):
    # Independent retained subpaths must support full-length copies of both
    # edge alleles. The next outside base must break their association.
    alleles = [source[a:b] for a, b in path]
    outside = _column(source, path, side, 1)
    if not outside:
        return None
    allelic_witnesses = {}
    for allele, occurrences in Counter(alleles).items():
        if occurrences < 3:
            continue
        bases = Counter(x for x, rep in zip(outside, alleles) if rep == allele)
        if 4 * max(bases.values()) >= 3 * occurrences:
            continue
        indexes = {pair for pair, rep in zip(path, alleles) if rep == allele}
        corroborators = [
            other
            for other, score in retained.items()
            if score >= 0.75 and len(other) >= 3 and set(other) <= indexes
        ]
        if corroborators:
            allelic_witnesses[allele] = max(corroborators, key=lambda g: (len(g), g))
    edge = 0 if side == "left" else -1
    if len({allele[edge] for allele in allelic_witnesses}) < 2:
        return None
    coverage = set().union(*allelic_witnesses.values())
    if 4 * len(coverage) < 3 * len(path):
        return None
    return {
        "side": side,
        "retained_subpaths": [
            [list(pair) for pair in witness]
            for witness in sorted(allelic_witnesses.values())
        ],
        "corroborated_copies": len(coverage),
    }


def _width_alternatives(source, path, repeats):
    counts = Counter(repeats)
    if len(counts) != 2:
        return
    (motif, common), (outlier, count) = counts.most_common()
    if count != 1 or common < 2 or len(outlier) != len(motif) + 1:
        return
    index = repeats.index(outlier)
    if index == 0 or index == len(path) - 1:
        return
    # Compare both one-base endpoint trims. Require the inferred difference to
    # occupy only the terminal dinucleotide, where a shifted end is plausible.
    # An interior indel or an exact retained end supports a different model.
    if (
        outlier[1:] == motif
        or outlier[:-1] == motif
        or not any(
            outlier[:i] + outlier[i + 1 :] == motif for i in range(1, len(motif))
        )
    ):
        return
    a, b = path[index]
    for side, changed in (("left", (a + 1, b)), ("right", (a, b - 1))):
        rep = source[slice(*changed)]
        mismatches = [i for i, (x, y) in enumerate(zip(rep, motif)) if x != y]
        terminal = [0, 1] if side == "left" else [len(rep) - 2, len(rep) - 1]
        if mismatches and set(mismatches) <= set(terminal):
            yield (
                path[:index] + (changed,) + path[index + 1 :],
                "uniform_edge_substitutions",
                {
                    "side": side,
                    "outlier_index_0based": index,
                    "exact_other_copies": common,
                    "uniform_positional_substitutions": len(mismatches),
                    "alternative_single_insertion": True,
                },
            )


def _corroborated_width(source, path, repeats):
    for proposed, rule, evidence in _width_alternatives(source, path, repeats):
        index = evidence["outlier_index_0based"]
        allele = source[slice(*proposed[index])]
        hits = _nearby_exemplars(source, proposed, allele, evidence["side"])
        if hits:
            yield (
                proposed,
                rule,
                dict(evidence, additional_exact_copy_starts_0based=hits),
            )


def _unanchored_width(source, path, repeats, index, side):
    """Determine whether the two width conventions need an explicit tie-break."""
    outlier = repeats[index]
    motif = repeats[0]
    if set(outlier) - set("ACGT"):
        return None
    gap_positions = [
        i for i in range(len(outlier)) if outlier[:i] + outlier[i + 1 :] == motif
    ]
    if len(gap_positions) < 2:
        return None
    # A long duplicated flank at *every* modal copy leaves this edge without
    # independent boundary support. A missing or ambiguous flank cannot count.
    width = (len(motif) + 1) // 2
    if min(c - b for (_, b), (c, _) in zip(path, path[1:])) < width:
        return None
    flanks = [
        source[max(0, a - width) : a] if side == "left" else source[b : b + width]
        for i, (a, b) in enumerate(path)
        if i != index
    ]
    flank = flanks[0]
    if len(set(flanks)) != 1 or len(flank) != width or set(flank) - set("ACGT"):
        return None
    # An additional exact allele on either strand offers a separate example:
    # in that situation keep the ordinary conservative/corroboration rules.
    shorter = outlier[1:] if side == "left" else outlier[:-1]
    a, _ = path[index]
    trim_start = a + 1 if side == "left" else a
    expected = {
        motif: {x for i, (x, _) in enumerate(path) if i != index},
        outlier: {a},
        shorter: {trim_start},
    }
    for allele in set(expected) | {reverse_complement(seq) for seq in expected}:
        pos = source.find(allele)
        while pos >= 0:
            if pos not in expected.get(allele, ()):
                return None
            pos = source.find(allele, pos + 1)
    return {
        "side": side,
        "identical_external_modal_flank": flank,
        "external_flank_bases_checked": width,
        "indistinguishable_insertion_positions_0based": gap_positions,
        "additional_exact_allele_occurrences_on_either_strand": 0,
    }


def select_alignment(source, path, parent):
    """Compare a modal-consensus insertion and an ungapped equal-width path."""
    if not 3 <= len(path) <= 64 or not _valid(path, len(source)):
        return None
    if any(not 20 <= b - a <= 60 for a, b in path):
        return None
    if any(c - b > 250 for (_, b), (c, _) in zip(path, path[1:])):
        return None
    repeats = [source[a:b] for a, b in path]
    choices = list(_width_alternatives(source, path, repeats))
    if len(choices) != 1 or any(_corroborated_width(source, path, repeats)):
        return None
    proposed, _, evidence = choices[0]
    index = evidence["outlier_index_0based"]
    old = tuple(parent.list_gaped_repeats)
    if (
        len(old) != len(repeats)
        or old[index] != repeats[index]
        or parent.consensus != repeats[index]
    ):
        return None
    short_rows = {s for i, s in enumerate(old) if i != index}
    if len(short_rows) != 1:
        return None
    (short_row,) = short_rows
    if short_row.count("-") != 1 or " " in short_row:
        return None
    gap = short_row.index("-")
    if (
        short_row[:gap] + short_row[gap + 1 :]
        != old[index][:gap] + old[index][gap + 1 :]
    ):
        return None
    rebuilt = build_source_backed_candidate_from_boundaries(source, None, *zip(*path))
    if rebuilt.consensus != short_row.replace("-", "") or rebuilt.total_mismatches != 1:
        return None
    details = {
        "outlier_index_0based": index,
        "previous_deletions_to_long_consensus": len(path) - 1,
        "insertions_to_modal_consensus": 1,
        "alternative_repeat_intervals": [list(pair) for pair in proposed],
        "alternative_metrics": {
            "uniform_positional_substitutions": evidence[
                "uniform_positional_substitutions"
            ],
            "spacer_lengths": [c - b for (_, b), (c, _) in zip(proposed, proposed[1:])],
        },
    }
    support = _unanchored_width(source, path, repeats, index, evidence["side"])
    if support is not None:
        details["uniform_width_preference"] = support
        return path, "ambiguous_edge_width_competition", details
    return path, "modal_insertion_alignment", details


def rank_alignment(choice, *, uniform_eligible, compact_score, uniform_score):
    """Use the ungapped convention if both interpretations qualify and it scores as well."""
    longer, rule, details = choice
    uniform = tuple(map(tuple, details["alternative_repeat_intervals"]))
    prefer_uniform = uniform_eligible and uniform_score >= compact_score
    ranked = dict(
        details,
        width_selection={
            "selected_interpretation": "uniform_edge_substitutions"
            if prefer_uniform
            else "modal_insertion",
            "uniform_qualified_for_selection": uniform_eligible,
            "uniform_score": float(uniform_score),
            "modal_insertion_score": float(compact_score),
        },
    )
    if not prefer_uniform:
        return longer, rule, ranked
    ranked["selected_metrics"] = details["alternative_metrics"]
    ranked["alternative_repeat_intervals"] = [list(pair) for pair in longer]
    ranked["alternative_metrics"] = {
        "insertions_to_modal_consensus": 1,
        "spacer_lengths": [c - b for (_, b), (c, _) in zip(longer, longer[1:])],
    }
    return uniform, rule, ranked


def _retained_trim(source, path, retained):
    if len(path) < 5:
        return
    for side in ("left", "right"):
        trimmed = tuple((a + 1, b) if side == "left" else (a, b - 1) for a, b in path)
        if trimmed not in retained or retained[trimmed] <= retained[path]:
            continue
        removed = Counter(source[a if side == "left" else b - 1] for a, b in path)
        if len(removed) < 3 or 3 * max(removed.values()) > 2 * len(path):
            continue
        copies = Counter(source[a:b] for a, b in trimmed)
        modal, count = copies.most_common(1)[0]
        if (
            count != len(path) - 1
            or sum(x != y for allele in copies for x, y in zip(allele, modal)) > 1
        ):
            continue
        exact = tuple(p for p in trimmed if source[slice(*p)] == modal)
        if exact not in retained or retained[exact] <= retained[path]:
            continue
        yield (
            trimmed,
            "retained_polymorphic_edge_trim",
            {
                "side": side,
                "edge_base_counts": dict(sorted(removed.items())),
                "exact_trimmed_copies": count,
                "exact_subpath_score": retained[exact],
                "complete_path_score": retained[trimmed],
            },
        )


def select(source, primary, retained):
    """Return the unique supported same-copy path, or abstain on competition."""
    if (
        not 3 <= len(primary) <= 64
        or primary not in retained
        or not _valid(primary, len(source))
    ):
        return None
    if any(not 20 <= b - a <= 60 for a, b in primary):
        return None
    if any(c - b > 250 for (_, b), (c, _) in zip(primary, primary[1:])):
        return None
    repeats = [source[a:b] for a, b in primary]
    spacers = [source[b:c] for (_, b), (c, _) in zip(primary, primary[1:])]
    if len(set(spacers)) != len(spacers) or any(
        set(s) - set("ACGT") for s in repeats + spacers
    ):
        return None
    options = list(_corroborated_width(source, primary, repeats))
    if len({len(rep) for rep in repeats}) == 1:
        options.extend(_retained_trim(source, primary, retained))
        for path, reason, evidence in _edges(source, primary, repeats):
            if not _valid(path, len(source)):
                continue
            if reason == "conserved_flank":
                support = _conserved_flank(
                    source, primary, path, evidence["side"], evidence["columns"]
                )
                rule = "conserved_flank_with_adjacent_exemplars"
            else:
                support = _allele_flanks(source, path, retained, evidence["side"])
                rule = "retained_allele_flanks"
            if support:
                options.append((path, rule, support))
    return options[0] if len(options) == 1 else None
