"""Bounded source-only completion of narrow two-repeat fast-run candidates."""

from collections import Counter
import math
import statistics

import numpy as np

from components.components_evaluated_arrays_enhancement import (
    SOURCE_GLOBAL_PHASE_RULE,
    _longest_exact_source_anchor,
    _source_only_alignment_seed,
    _source_phase_background,
    _source_phase_geometry_score,
    _source_repeat_edit_distance,
    build_source_backed_candidate_from_boundaries,
    validate_candidate_source_reconstruction,
)
from components.fast_run_rescue import candidate_source_interval


MAXIMUM_INWARD_TRIM = 14
CORE_PHASE_SCORE_DELTA = 3.0
MINIMUM_CORE_IDENTITY = 0.85
SECOND_SPACER_RADIUS = 10
TERMINAL_LENGTH_RADIUS = 2
MINIMUM_TERMINAL_IDENTITY = 0.60
MINIMUM_TERMINAL_EXACT_ANCHOR = 10
ANCHOR_ALIGNMENT_OFFSET_RADIUS = 2

# The sealed 749-parent audit observed maxima of 1,515 anchor hits, 37,080
# terminal hypotheses, and 3,410 raw proposals. Every cap fails closed.
MAXIMUM_ANCHOR_HITS_PER_SEED = 2048
MAXIMUM_TERMINAL_HYPOTHESES_PER_SEED = 50000
MAXIMUM_RAW_PROPOSALS_PER_SEED = 4096
MAXIMUM_SOURCE_RANKED_MATERIALIZATIONS = 25
MAXIMUM_SEEDS_PER_SOURCE = 1024
MAXIMUM_LOCUS_WINNERS_PER_SOURCE = 64
REPEAT_DERIVED_SPACER_KMER_LENGTH = 6
MINIMUM_TANDEM_REPEAT_DUST_SCORE = 0.60
MINIMUM_TANDEM_REPEAT_DERIVED_SPACER_COVERAGE = 0.35
MINIMUM_ADMISSION_MODEL_SCORE = 0.57
MAXIMUM_ADMISSION_REPEAT_DERIVED_SPACER_COVERAGE = 0.75
FLANK_CONTINUATION_RADIUS = 250
MINIMUM_REPEAT_FLANK_CONTINUATION_IDENTITY = 0.90
MINIMUM_SPACER_FLANK_CONTINUATION_IDENTITY = 0.55


def _edit_identity(first, second):
    denominator = max(len(first), len(second))
    if denominator == 0:
        return 1.0
    distance = _source_repeat_edit_distance(first, second)
    return 1.0 - distance / float(denominator)


def _spacer_statistics(spacers):
    lengths = tuple(map(len, spacers))
    mean = sum(lengths) / float(len(lengths)) if lengths else 0.0
    cv = statistics.pstdev(lengths) / mean if mean else 0.0
    maximum_identity = max((
        _edit_identity(first, second)
        for index, first in enumerate(spacers[:-1])
        for second in spacers[index + 1:]
    ), default=0.0)
    return cv, maximum_identity


def repeat_derived_spacer_coverage(candidate):
    """Maximum spacer fraction covered by exact repeat-derived 6-mers."""

    kmer_length = REPEAT_DERIVED_SPACER_KMER_LENGTH
    repeat_kmers = {
        repeat[offset:offset + kmer_length]
        for repeat in candidate.list_repeats
        for offset in range(len(repeat) - kmer_length + 1)
    }
    coverages = []
    for spacer in candidate.list_spacers:
        covered = [False] * len(spacer)
        for offset in range(len(spacer) - kmer_length + 1):
            if spacer[offset:offset + kmer_length] not in repeat_kmers:
                continue
            for position in range(offset, offset + kmer_length):
                covered[position] = True
        coverages.append(
            sum(covered) / float(len(spacer)) if spacer else 1.0
        )
    return max(coverages, default=1.0)


def repeat_dust_triplet_score(candidate, repeat=None):
    """DUST-style triplet-pair density of the source-only repeat medoid."""

    if repeat is None:
        repeat = _source_only_alignment_seed(tuple(candidate.list_repeats))
    triplets = Counter(
        repeat[offset:offset + 3]
        for offset in range(len(repeat) - 2)
    )
    pairs = sum(
        count * (count - 1) // 2 for count in triplets.values()
    )
    return pairs / float(max(1, len(repeat) - 3))


def _maximum_same_length_flank_identity(query, flanks):
    """Return the best ungapped identity to either separate source flank."""

    if not isinstance(query, str) or not query:
        raise ValueError("flank-continuation queries must be non-empty strings")
    query_length = len(query)
    query_codes = np.fromiter(
        map(ord, query), dtype=np.uint32, count=query_length
    )
    maximum_matches = 0
    for flank in flanks:
        if not isinstance(flank, str):
            raise ValueError("flank-continuation source flanks must be strings")
        if len(flank) < query_length:
            continue
        flank_codes = np.fromiter(
            map(ord, flank), dtype=np.uint32, count=len(flank)
        )
        # This is the exact Hamming calculation performed by the former
        # nested Python loops.  The view is bounded by the separate 250-base
        # flank contract and does not copy source windows; only the boolean
        # equality matrix is materialized.
        windows = np.lib.stride_tricks.sliding_window_view(
            flank_codes, query_length
        )
        matches = int(np.count_nonzero(
            windows == query_codes, axis=1
        ).max(initial=0))
        if matches > maximum_matches:
            maximum_matches = matches
            if maximum_matches == query_length:
                return 1.0
    return maximum_matches / float(query_length)


def flank_continuation_identities(source, candidate):
    """Measure repeat and spacer continuation outside a fixed array envelope.

    The two 250-base source windows remain separate so no synthetic match can
    bridge the candidate.  Queries use literal, same-length, ungapped windows;
    neither reverse complements nor an external CRISPR comparator participate.
    """

    if not isinstance(source, str):
        raise ValueError("flank-continuation source must be a string")
    start, end = candidate_source_interval(candidate)
    if not 0 <= start < end <= len(source):
        raise ValueError("flank-continuation candidate is outside source")
    flanks = (
        source[max(0, start - FLANK_CONTINUATION_RADIUS):start],
        source[end:min(len(source), end + FLANK_CONTINUATION_RADIUS)],
    )
    repeats = tuple(candidate.list_repeats)
    spacers = tuple(candidate.list_spacers)
    if not repeats or not spacers:
        raise ValueError(
            "flank-continuation candidate needs repeats and spacers"
        )
    return {
        "maximum_repeat_flank_identity": max(
            _maximum_same_length_flank_identity(repeat, flanks)
            for repeat in repeats
        ),
        "maximum_spacer_flank_identity": max(
            _maximum_same_length_flank_identity(spacer, flanks)
            for spacer in spacers
        ),
    }


def fast_two_repeat_admission_allowed(
        model_score,
        repeat_derived_spacer_fraction,
        maximum_repeat_flank_identity,
        maximum_spacer_flank_identity):
    """Apply the comparator-blind precision policy for this synthetic stage."""

    values = (
        model_score,
        repeat_derived_spacer_fraction,
        maximum_repeat_flank_identity,
        maximum_spacer_flank_identity,
    )
    if any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            for value in values):
        return False
    if not all(0.0 <= float(value) <= 1.0 for value in values[1:]):
        return False
    if float(model_score) < MINIMUM_ADMISSION_MODEL_SCORE:
        return False
    if (
            float(repeat_derived_spacer_fraction)
            >= MAXIMUM_ADMISSION_REPEAT_DERIVED_SPACER_COVERAGE):
        return False
    source_continues = (
        float(maximum_repeat_flank_identity)
        >= MINIMUM_REPEAT_FLANK_CONTINUATION_IDENTITY
        and float(maximum_spacer_flank_identity)
        >= MINIMUM_SPACER_FLANK_CONTINUATION_IDENTITY
    )
    return not source_continues


def _geometry(candidate):
    return tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            candidate.list_repeat_starts, candidate.list_repeats
        )
    )


def _source_rank(record):
    """Fixed comparator-blind ordering; greater tuples rank first."""

    return (
        round(record["source_phase_score"], 12),
        round(record["minimum_repeat_identity"], 12),
        round(record["mean_repeat_identity"], 12),
        record["terminal_exact_anchor"],
        -round(record["spacer_cv"], 12),
        -round(record["maximum_spacer_identity"], 12),
        -abs(record["second_spacer"] - record["first_spacer"]),
        -(record["left_trim"] + record["right_trim"]),
        -record["left_trim"],
        -record["right_trim"],
        1 if record["side"] == "right" else 0,
        tuple(
            -value
            for interval in record["repeat_intervals"]
            for value in interval
        ),
    )


def _phase_rule(minimum_repeats):
    return dict(
        SOURCE_GLOBAL_PHASE_RULE,
        minimum_repeats=minimum_repeats,
        minimum_repeat_length=21,
        maximum_repeat_length=55,
        minimum_mean_repeat_length=21,
        maximum_mean_repeat_length=55,
        minimum_spacer_length=18,
        maximum_spacer_length=78,
    )


def _core_phases(source, starts, ends, probabilities, counters):
    phases = []
    rule = _phase_rule(2)
    for left_trim in range(MAXIMUM_INWARD_TRIM + 1):
        phase_starts = tuple(start + left_trim for start in starts)
        for right_trim in range(MAXIMUM_INWARD_TRIM + 1):
            counters["phase_hypotheses"] += 1
            phase_ends = tuple(end - right_trim for end in ends)
            repeats = tuple(
                source[start:end]
                for start, end in zip(phase_starts, phase_ends)
            )
            if (
                    not repeats
                    or min(map(len, repeats)) < 21
                    or max(map(len, repeats)) > 55):
                continue
            medoid = _source_only_alignment_seed(repeats)
            if min(
                    _edit_identity(repeat, medoid)
                    for repeat in repeats) < MINIMUM_CORE_IDENTITY:
                continue
            score = _source_phase_geometry_score(
                source, phase_starts, phase_ends, rule, probabilities
            )
            if score is None:
                continue
            phases.append({
                "left_trim": left_trim,
                "right_trim": right_trim,
                "starts": phase_starts,
                "ends": phase_ends,
                "medoid": medoid,
                "core_source_score": float(score),
            })
    if not phases:
        return ()
    best = max(phase["core_source_score"] for phase in phases)
    retained = tuple(
        phase for phase in phases
        if best - phase["core_source_score"] <= CORE_PHASE_SCORE_DELTA
    )
    counters["retained_core_phases"] = len(retained)
    return retained


def _exact_anchor_hits(source, medoid, search_start, search_end):
    hits = set()
    for medoid_offset in range(
            len(medoid) - MINIMUM_TERMINAL_EXACT_ANCHOR + 1):
        anchor = medoid[
            medoid_offset:
            medoid_offset + MINIMUM_TERMINAL_EXACT_ANCHOR
        ]
        hit = source.find(anchor, search_start, search_end)
        while hit >= 0:
            hits.add((hit, medoid_offset))
            hit = source.find(anchor, hit + 1, search_end)
    return tuple(sorted(hits))


def build_source_winner(source, seed, structural_filter):
    """Return one source-ranked candidate for a seed, or None."""

    counters = Counter()
    starts = tuple(seed.repeat_starts)
    ends = tuple(seed.repeat_ends)
    if not (
            len(starts) == len(ends) == 2
            and all(0 <= start < end <= len(source)
                    for start, end in zip(starts, ends))
            and ends[0] <= starts[1]):
        raise ValueError("two-repeat terminal repair seed is outside source")
    parent = build_source_backed_candidate_from_boundaries(
        source, None, starts, ends
    )
    validate_candidate_source_reconstruction(
        source,
        parent.list_repeats,
        parent.list_repeat_starts,
        parent.list_spacers,
    )
    background, probabilities = _source_phase_background(
        source, (parent,), 5000
    )
    phases = _core_phases(source, starts, ends, probabilities, counters)
    final_rule = _phase_rule(3)
    records = []
    seen = set()
    cap_reason = None

    for phase in phases:
        phase_starts = phase["starts"]
        phase_ends = phase["ends"]
        lengths = tuple(
            end - start for start, end in zip(phase_starts, phase_ends)
        )
        first_spacer = phase_starts[1] - phase_ends[0]
        spacer_lower = max(18, first_spacer - SECOND_SPACER_RADIUS)
        spacer_upper = min(78, first_spacer + SECOND_SPACER_RADIUS)
        length_lower = max(21, min(lengths) - TERMINAL_LENGTH_RADIUS)
        length_upper = min(55, max(lengths) + TERMINAL_LENGTH_RADIUS)
        for side in ("left", "right"):
            if side == "right":
                search_start = max(0, phase_ends[-1] + spacer_lower)
                search_end = min(
                    len(source),
                    phase_ends[-1] + spacer_upper + length_upper,
                )
            else:
                maximum_end = phase_starts[0] - spacer_lower
                minimum_end = phase_starts[0] - spacer_upper
                search_start = max(0, minimum_end - length_upper)
                search_end = min(len(source), maximum_end)
            hits = _exact_anchor_hits(
                source, phase["medoid"], search_start, search_end
            )
            counters["exact_anchor_hits"] += len(hits)
            if (
                    counters["exact_anchor_hits"]
                    > MAXIMUM_ANCHOR_HITS_PER_SEED):
                cap_reason = "anchor_hits"
                break
            for hit, medoid_offset in hits:
                for offset_delta in range(
                        -ANCHOR_ALIGNMENT_OFFSET_RADIUS,
                        ANCHOR_ALIGNMENT_OFFSET_RADIUS + 1):
                    terminal_start = hit - medoid_offset + offset_delta
                    for terminal_length in range(
                            length_lower, length_upper + 1):
                        counters["terminal_hypotheses"] += 1
                        if (
                                counters["terminal_hypotheses"]
                                > MAXIMUM_TERMINAL_HYPOTHESES_PER_SEED):
                            cap_reason = "terminal_hypotheses"
                            break
                        terminal_end = terminal_start + terminal_length
                        if not (
                                0 <= terminal_start < terminal_end
                                <= len(source)):
                            continue
                        if side == "right":
                            second_spacer = terminal_start - phase_ends[-1]
                            repeat_starts = phase_starts + (terminal_start,)
                            repeat_ends = phase_ends + (terminal_end,)
                            terminal_index = 2
                        else:
                            second_spacer = phase_starts[0] - terminal_end
                            repeat_starts = (terminal_start,) + phase_starts
                            repeat_ends = (terminal_end,) + phase_ends
                            terminal_index = 0
                        if not spacer_lower <= second_spacer <= spacer_upper:
                            continue
                        intervals = tuple(zip(repeat_starts, repeat_ends))
                        if intervals in seen:
                            continue
                        seen.add(intervals)
                        terminal = source[terminal_start:terminal_end]
                        if (
                                _edit_identity(terminal, phase["medoid"])
                                < MINIMUM_TERMINAL_IDENTITY):
                            continue
                        terminal_anchor = _longest_exact_source_anchor(
                            terminal, phase["medoid"]
                        )
                        if (
                                terminal_anchor
                                < MINIMUM_TERMINAL_EXACT_ANCHOR):
                            continue
                        repeats = tuple(
                            source[start:end] for start, end in intervals
                        )
                        medoid = _source_only_alignment_seed(repeats)
                        identities = tuple(
                            _edit_identity(repeat, medoid)
                            for repeat in repeats
                        )
                        source_score = _source_phase_geometry_score(
                            source,
                            repeat_starts,
                            repeat_ends,
                            final_rule,
                            probabilities,
                        )
                        if source_score is None:
                            continue
                        spacers = tuple(
                            source[left_end:right_start]
                            for left_end, right_start in zip(
                                repeat_ends, repeat_starts[1:]
                            )
                        )
                        spacer_cv, spacer_identity = _spacer_statistics(
                            spacers
                        )
                        records.append({
                            "repeat_intervals": intervals,
                            "side": side,
                            "terminal_index": terminal_index,
                            "left_trim": phase["left_trim"],
                            "right_trim": phase["right_trim"],
                            "first_spacer": first_spacer,
                            "second_spacer": second_spacer,
                            "terminal_length": terminal_length,
                            "terminal_exact_anchor": terminal_anchor,
                            "core_medoid": phase["medoid"],
                            "source_phase_score": float(source_score),
                            "minimum_repeat_identity": min(identities),
                            "mean_repeat_identity": (
                                sum(identities) / float(len(identities))
                            ),
                            "spacer_cv": spacer_cv,
                            "maximum_spacer_identity": spacer_identity,
                            "seed_interval": seed.source_interval,
                        })
                        counters["raw_source_proposals"] += 1
                        if (
                                counters["raw_source_proposals"]
                                > MAXIMUM_RAW_PROPOSALS_PER_SEED):
                            cap_reason = "raw_source_proposals"
                            break
                    if cap_reason is not None:
                        break
                if cap_reason is not None:
                    break
            if cap_reason is not None:
                break
        if cap_reason is not None:
            break

    counters["cap_exhausted"] = int(cap_reason is not None)
    if cap_reason is not None:
        counters["cap_reason_" + cap_reason] = 1
        return None, dict(counters), background

    records.sort(key=_source_rank, reverse=True)
    for raw_rank, record in enumerate(
            records[:MAXIMUM_SOURCE_RANKED_MATERIALIZATIONS], 1):
        intervals = record["repeat_intervals"]
        candidate = build_source_backed_candidate_from_boundaries(
            source,
            None,
            tuple(start for start, _end in intervals),
            tuple(end for _start, end in intervals),
        )
        counters["materialized"] += 1
        if not structural_filter(candidate):
            counters["materialized_structural_rejections"] += 1
            continue
        repeat_spacer_coverage = repeat_derived_spacer_coverage(candidate)
        repeat_dust_score = repeat_dust_triplet_score(
            candidate, record["core_medoid"]
        )
        record["maximum_repeat_derived_spacer_coverage"] = (
            repeat_spacer_coverage
        )
        record["repeat_dust_triplet_score"] = repeat_dust_score
        if (
                repeat_dust_score >= MINIMUM_TANDEM_REPEAT_DUST_SCORE
                and repeat_spacer_coverage >= (
                    MINIMUM_TANDEM_REPEAT_DERIVED_SPACER_COVERAGE
                )):
            counters["source_tandem_guard_rejections"] += 1
            return None, dict(counters), background
        validate_candidate_source_reconstruction(
            source,
            candidate.list_repeats,
            candidate.list_repeat_starts,
            candidate.list_spacers,
        )
        record["candidate"] = candidate
        record["raw_source_rank"] = raw_rank
        return record, dict(counters), background
    return None, dict(counters), background


def _overlap_components(records):
    ordered = sorted(
        records,
        key=lambda record: (
            candidate_source_interval(record["candidate"]),
            _geometry(record["candidate"]),
        ),
    )
    components = []
    current = []
    current_end = None
    for record in ordered:
        start, end = candidate_source_interval(record["candidate"])
        if not current or start < current_end:
            current.append(record)
            current_end = end if current_end is None else max(current_end, end)
        else:
            components.append(current)
            current = [record]
            current_end = end
    if current:
        components.append(current)
    return components


def generate_source_winners(source, seeds, structural_filter):
    """Generate bounded classifier-blind locus winners and telemetry."""

    seeds = tuple(sorted(set(seeds)))
    stats = Counter(
        seeds_input=len(seeds),
        seed_cap=MAXIMUM_SEEDS_PER_SOURCE,
        locus_winner_cap=MAXIMUM_LOCUS_WINNERS_PER_SOURCE,
    )
    if len(seeds) > MAXIMUM_SEEDS_PER_SOURCE:
        stats["source_seed_cap_exhausted"] = 1
        return (), dict(stats)

    winners = []
    seen_geometries = set()
    for seed in seeds:
        winner, counters, _background = build_source_winner(
            source, seed, structural_filter
        )
        for name, value in counters.items():
            stats[name] += value
        if winner is None:
            continue
        geometry = _geometry(winner["candidate"])
        if geometry in seen_geometries:
            stats["duplicate_winners"] += 1
            continue
        seen_geometries.add(geometry)
        winners.append(winner)
    stats["seed_winners"] = len(winners)

    locus_winners = [
        max(component, key=_source_rank)
        for component in _overlap_components(winners)
    ]
    stats["overlap_components"] = len(locus_winners)
    if len(locus_winners) > MAXIMUM_LOCUS_WINNERS_PER_SOURCE:
        locus_winners = sorted(
            locus_winners, key=_source_rank, reverse=True
        )[:MAXIMUM_LOCUS_WINNERS_PER_SOURCE]
        stats["locus_winner_cap_reached"] = 1
    locus_winners.sort(
        key=lambda record: candidate_source_interval(record["candidate"])
    )
    stats["locus_winners"] = len(locus_winners)
    return tuple(locus_winners), dict(stats)
