import math
import statistics
import sys
from collections import Counter, OrderedDict
from functools import lru_cache
from threading import RLock

import regex
from Bio.Align import PairwiseAligner

from components.components_detection_refinement import CrisprCandidate
from components.candidate_origins import inherit_candidate_origins
from components.fuzzy_alignment import normalize_fuzzy_changes
from components.fuzzy_alignment_merge import (
    merge_fuzzy_alignments,
)


_SOURCE_PAIR_CACHE_MAX_LENGTH = 128
_SOURCE_DISTANCE_CACHE_MAXSIZE = 131072
_SOURCE_EDIT_CACHE_MAXSIZE = 32768
_SOURCE_ALIGNMENT_CACHE_MAXSIZE = 2048
_SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES = 64 * 1024 * 1024
_SOURCE_ALIGNMENT_CACHE_MAX_BASES = 65536
_SOURCE_ALIGNMENT_CACHE_MAX_REPEATS = 1024
_SOURCE_ALIGNMENT_CACHE_MAX_VALUE_BASES = 131072
_SOURCE_ALIGNMENT_CACHE = OrderedDict()
_SOURCE_ALIGNMENT_CACHE_LOCK = RLock()
_SOURCE_ALIGNMENT_CACHE_STATS = Counter()
_SOURCE_ALIGNMENT_CACHE_GENERATION = 0


def _minimal_edit_changes_uncached(pattern, matched_text):
    """Return one deterministic minimal edit script in alignment columns."""
    if pattern == matched_text:
        return [[], [], []]
    rows = len(pattern) + 1
    columns = len(matched_text) + 1
    distances = [[0] * columns for _ in range(rows)]

    for row in range(rows):
        distances[row][0] = row
    for column in range(columns):
        distances[0][column] = column

    for row in range(1, rows):
        for column in range(1, columns):
            substitution_cost = pattern[row - 1] != matched_text[column - 1]
            distances[row][column] = min(
                distances[row - 1][column - 1] + substitution_cost,
                distances[row - 1][column] + 1,
                distances[row][column - 1] + 1,
            )

    operations = []
    row = len(pattern)
    column = len(matched_text)
    while row or column:
        if row and column:
            substitution_cost = pattern[row - 1] != matched_text[column - 1]
            if (distances[row][column]
                    == distances[row - 1][column - 1] + substitution_cost):
                operations.append("match" if not substitution_cost else "substitution")
                row -= 1
                column -= 1
                continue
        if row and distances[row][column] == distances[row - 1][column] + 1:
            operations.append("deletion")
            row -= 1
            continue
        operations.append("insertion")
        column -= 1

    operations.reverse()
    changes = [[], [], []]
    for alignment_column, operation in enumerate(operations):
        if operation == "substitution":
            changes[0].append(alignment_column)
        elif operation == "insertion":
            changes[1].append(alignment_column)
        elif operation == "deletion":
            changes[2].append(alignment_column)
    return changes


@lru_cache(maxsize=_SOURCE_EDIT_CACHE_MAXSIZE)
def _cached_minimal_edit_changes(pattern, matched_text):
    return tuple(tuple(positions) for positions in
                 _minimal_edit_changes_uncached(pattern, matched_text))


def _minimal_edit_changes(pattern, matched_text):
    """Return fresh mutable edit lists backed by bounded immutable reuse."""
    if pattern == matched_text:
        return [[], [], []]
    if max(len(pattern), len(matched_text)) > _SOURCE_PAIR_CACHE_MAX_LENGTH:
        _SOURCE_ALIGNMENT_CACHE_STATS["edit_oversize_bypasses"] += 1
        return _minimal_edit_changes_uncached(pattern, matched_text)
    return [list(positions) for positions in
            _cached_minimal_edit_changes(pattern, matched_text)]


_SOURCE_REPEAT_BITPARALLEL_MAX_LENGTH = 64


def _make_source_unit_edit_scorer():
    """Keep one fixed scalar scorer private; never request native traceback.

    Biopython 1.86's score entry point retains the GIL and allocates its score
    buffers per call. With this configuration captured privately and never
    mutated after initialization, sharing the scorer does not share alignment
    paths or mutable candidate state across calls/threads.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = None
    aligner.match_score = 0.0
    aligner.mismatch_score = -1.0
    aligner.wildcard = None
    gap_names = tuple(
        "{}_{}_{}_score".format(action, location, kind)
        for action in ("open", "extend")
        for location in ("internal", "left", "right")
        for kind in ("insertion", "deletion")
    )
    for name in gap_names:
        setattr(aligner, name, -1.0)
    if (aligner.mode != "global" or aligner.substitution_matrix is not None
            or aligner.match_score != 0.0 or aligner.mismatch_score != -1.0
            or aligner.wildcard is not None
            or any(getattr(aligner, name) != -1.0 for name in gap_names)
            or aligner.algorithm != "Needleman-Wunsch"):
        raise RuntimeError("source unit-edit scorer configuration is invalid")
    score = aligner.score

    def scalar_score(first, second):
        return score(first, second, strand="+")

    return scalar_score


_SOURCE_UNIT_EDIT_SCORE = _make_source_unit_edit_scorer()


def _source_repeat_edit_distance_python(first, second):
    """Return exact Levenshtein distance for source-sliced repeats.

    CRISPR repeats are normally short enough to fit in one 64-bit bit vector.
    Myers' recurrence evaluates one text character with a fixed number of
    integer operations, while preserving the literal character equality used
    by ``_minimal_edit_changes`` (including ambiguous nucleotide symbols).
    Unusually long inputs retain the legacy dynamic-programming semantics.
    """
    if first == second:
        return 0
    if not first:
        return len(second)
    if not second:
        return len(first)
    if max(len(first), len(second)) > _SOURCE_REPEAT_BITPARALLEL_MAX_LENGTH:
        return sum(
            len(changes)
            for changes in _minimal_edit_changes(first, second)
        )

    # The shorter sequence is the bit-vector pattern.  Levenshtein distance is
    # symmetric, and using the shorter pattern minimizes the active bit mask.
    if len(first) > len(second):
        first, second = second, first

    pattern_length = len(first)
    active_mask = (1 << pattern_length) - 1
    highest_bit = 1 << (pattern_length - 1)
    character_masks = {}
    for index, character in enumerate(first):
        character_masks[character] = (
            character_masks.get(character, 0) | (1 << index)
        )

    positive = active_mask
    negative = 0
    distance = pattern_length
    for character in second:
        equal = character_masks.get(character, 0)
        vertical = equal | negative
        horizontal = (((equal & positive) + positive) ^ positive) | equal
        positive_horizontal = negative | ~(horizontal | positive)
        negative_horizontal = positive & horizontal

        if positive_horizontal & highest_bit:
            distance += 1
        elif negative_horizontal & highest_bit:
            distance -= 1

        positive_horizontal = ((positive_horizontal << 1) | 1) & active_mask
        negative_horizontal = (negative_horizontal << 1) & active_mask
        positive = (
            negative_horizontal | ~(vertical | positive_horizontal)
        ) & active_mask
        negative = positive_horizontal & vertical

    return distance


def _source_repeat_edit_distance_uncached(first, second):
    """Exact scalar unit Levenshtein; emitted edit paths remain Python-owned."""
    if first == second:
        return 0
    if not first:
        return len(second)
    if not second:
        return len(first)
    maximum_length = max(len(first), len(second))
    if (maximum_length <= _SOURCE_REPEAT_BITPARALLEL_MAX_LENGTH
            and isinstance(first, str) and isinstance(second, str)
            and first.isascii() and second.isascii()):
        score = _SOURCE_UNIT_EDIT_SCORE(first, second)
        if isinstance(score, bool) or not isinstance(score, (int, float)):
            raise RuntimeError("source unit-edit scorer returned a nonnumeric score")
        distance = -score
        if (not math.isfinite(distance) or not float(distance).is_integer()
                or not 0 <= distance <= maximum_length):
            raise RuntimeError("source unit-edit scorer returned an invalid distance")
        return int(distance)
    return _source_repeat_edit_distance_python(first, second)


@lru_cache(maxsize=_SOURCE_DISTANCE_CACHE_MAXSIZE)
def _cached_source_repeat_edit_distance(first, second):
    return _source_repeat_edit_distance_uncached(first, second)


def _source_repeat_edit_distance(first, second):
    """Exact symmetric distance with bounded, source-string-only memoization."""
    if first == second:
        return 0
    if max(len(first), len(second)) > _SOURCE_PAIR_CACHE_MAX_LENGTH:
        _SOURCE_ALIGNMENT_CACHE_STATS["distance_oversize_bypasses"] += 1
        return _source_repeat_edit_distance_uncached(first, second)
    first, second = sorted((first, second))
    return _cached_source_repeat_edit_distance(first, second)


def _source_only_alignment_seed(repeats):
    """Choose a deterministic medoid using only source-sliced repeats.

    A legacy candidate's gapped alignment and consensus are hypotheses rather
    than source facts.  Reusing either as an alignment seed lets two candidates
    with identical repeat intervals rebuild to different consensuses and
    mismatch counts.  The weighted edit-distance medoid below is invariant to
    candidate identity and input representation.  Exact-repeat multiplicity is
    retained so the usual repeat sequence wins a distance tie before the final
    source-sequence tie breakers.
    """

    repeat_counts = Counter(repeats)
    unique_repeats = tuple(sorted(repeat_counts))
    counts = tuple(repeat_counts[repeat] for repeat in unique_repeats)
    totals = [0] * len(unique_repeats)
    # Exact symmetry gives both weighted objective contributions from one
    # unordered-pair distance; integer totals preserve every previous tie.
    for first_index, first in enumerate(unique_repeats):
        for second_index in range(first_index + 1, len(unique_repeats)):
            distance = _source_repeat_edit_distance(first, unique_repeats[second_index])
            totals[first_index] += counts[second_index] * distance
            totals[second_index] += counts[first_index] * distance
    winner = min(
        range(len(unique_repeats)),
        key=lambda index: (
            totals[index], -counts[index], -len(unique_repeats[index]),
            unique_repeats[index],
        ),
    )
    return unique_repeats[winner]


def _source_repeat_alignment_uncached(repeats):
    seed = _source_only_alignment_seed(repeats)
    changes_by_repeat = {
        repeat: tuple(tuple(indexes) for indexes in
                      _minimal_edit_changes(seed, repeat))
        for repeat in dict.fromkeys(repeats)
    }
    _, alignment = merge_fuzzy_alignments(
        seed, repeats, [changes_by_repeat[repeat] for repeat in repeats],
    )
    return tuple(alignment)


def _source_repeat_alignment(repeats):
    """Cache only immutable alignment strings, never geometry or candidates.

    The key preserves complete order and multiplicity because both medoid
    choice and output row order depend on them. Oversized keys and aligned
    values bypass retention so a pathological alignment cannot inflate the
    bounded cache. Eviction changes cost only.
    """
    key = tuple(repeats)
    cacheable = (len(key) <= _SOURCE_ALIGNMENT_CACHE_MAX_REPEATS
                 and sum(map(len, key)) <= _SOURCE_ALIGNMENT_CACHE_MAX_BASES)
    if cacheable:
        with _SOURCE_ALIGNMENT_CACHE_LOCK:
            generation = _SOURCE_ALIGNMENT_CACHE_GENERATION
            if key in _SOURCE_ALIGNMENT_CACHE:
                _SOURCE_ALIGNMENT_CACHE_STATS["alignment_hits"] += 1
                _SOURCE_ALIGNMENT_CACHE.move_to_end(key)
                return _SOURCE_ALIGNMENT_CACHE[key][0]
            _SOURCE_ALIGNMENT_CACHE_STATS["alignment_misses"] += 1
    else:
        with _SOURCE_ALIGNMENT_CACHE_LOCK:
            _SOURCE_ALIGNMENT_CACHE_STATS["alignment_oversize_key_bypasses"] += 1
    alignment = _source_repeat_alignment_uncached(key)
    if cacheable and sum(map(len, alignment)) <= _SOURCE_ALIGNMENT_CACHE_MAX_VALUE_BASES:
        payload = _source_alignment_payload_bytes(key, alignment)
        with _SOURCE_ALIGNMENT_CACHE_LOCK:
            if generation != _SOURCE_ALIGNMENT_CACHE_GENERATION:
                _SOURCE_ALIGNMENT_CACHE_STATS["alignment_reset_race_bypasses"] += 1
                return alignment
            if key in _SOURCE_ALIGNMENT_CACHE:
                existing = _SOURCE_ALIGNMENT_CACHE[key][0]
                if existing != alignment:
                    raise RuntimeError("duplicate source-alignment computations disagree")
                _SOURCE_ALIGNMENT_CACHE.move_to_end(key)
                _SOURCE_ALIGNMENT_CACHE_STATS["alignment_duplicate_insertions"] += 1
                return existing
            if payload > _SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES:
                _SOURCE_ALIGNMENT_CACHE_STATS["alignment_payload_oversize_bypasses"] += 1
                return alignment
            current = _SOURCE_ALIGNMENT_CACHE_STATS["retained_payload_bytes"]
            while (len(_SOURCE_ALIGNMENT_CACHE) >= _SOURCE_ALIGNMENT_CACHE_MAXSIZE
                   or current + payload > _SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES):
                if current + payload > _SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES:
                    _SOURCE_ALIGNMENT_CACHE_STATS["alignment_byte_limit_evictions"] += 1
                if len(_SOURCE_ALIGNMENT_CACHE) >= _SOURCE_ALIGNMENT_CACHE_MAXSIZE:
                    _SOURCE_ALIGNMENT_CACHE_STATS["alignment_entry_limit_evictions"] += 1
                _, (_, evicted_payload) = _SOURCE_ALIGNMENT_CACHE.popitem(last=False)
                current -= evicted_payload
                _SOURCE_ALIGNMENT_CACHE_STATS["alignment_evictions"] += 1
            _SOURCE_ALIGNMENT_CACHE[key] = (alignment, payload)
            current += payload
            _SOURCE_ALIGNMENT_CACHE_STATS["retained_payload_bytes"] = current
            _SOURCE_ALIGNMENT_CACHE_STATS["peak_retained_payload_bytes"] = max(
                _SOURCE_ALIGNMENT_CACHE_STATS["peak_retained_payload_bytes"], current,
            )
    elif cacheable:
        with _SOURCE_ALIGNMENT_CACHE_LOCK:
            _SOURCE_ALIGNMENT_CACHE_STATS["alignment_oversize_value_bypasses"] += 1
    return alignment


def _source_alignment_payload_bytes(key, alignment):
    """Conservative retained tuple/string payload, counting repeated refs.

    Includes the key tuple, alignment tuple, their referenced strings, and
    the immutable (alignment, payload-size) entry tuple/integer. Dictionary,
    allocator overhead and the separate pair caches are outside this measure;
    the cgroup memory bound remains the authority for total process memory.
    """
    return (sys.getsizeof(key) + sum(map(sys.getsizeof, key))
            + sys.getsizeof(alignment) + sum(map(sys.getsizeof, alignment))
            + sys.getsizeof((alignment, 0)) + sys.getsizeof(0))


def clear_source_alignment_caches():
    """Reset pure computation caches; useful for cold/warm equivalence tests."""
    global _SOURCE_ALIGNMENT_CACHE_GENERATION
    _cached_source_repeat_edit_distance.cache_clear()
    _cached_minimal_edit_changes.cache_clear()
    with _SOURCE_ALIGNMENT_CACHE_LOCK:
        _SOURCE_ALIGNMENT_CACHE.clear()
        _SOURCE_ALIGNMENT_CACHE_STATS.clear()
        _SOURCE_ALIGNMENT_CACHE_GENERATION += 1


def source_alignment_cache_info():
    """Return auditable bounded-cache counters without sequence contents."""
    with _SOURCE_ALIGNMENT_CACHE_LOCK:
        return {
            "distance": dict(_cached_source_repeat_edit_distance.cache_info()._asdict()),
            "edit_scripts": dict(_cached_minimal_edit_changes.cache_info()._asdict()),
            "alignment": dict(_SOURCE_ALIGNMENT_CACHE_STATS,
                              currsize=len(_SOURCE_ALIGNMENT_CACHE),
                              maxsize=_SOURCE_ALIGNMENT_CACHE_MAXSIZE,
                              retained_payload_bytes=_SOURCE_ALIGNMENT_CACHE_STATS["retained_payload_bytes"],
                              peak_retained_payload_bytes=_SOURCE_ALIGNMENT_CACHE_STATS["peak_retained_payload_bytes"],
                              maximum_retained_payload_bytes=_SOURCE_ALIGNMENT_CACHE_MAX_PAYLOAD_BYTES,
                              payload_accounting="sys.getsizeof key/alignment/entry tuples, strings counted per reference, entry integer; excludes dict/allocator/pair-cache overhead",
                              maximum_key_bases=_SOURCE_ALIGNMENT_CACHE_MAX_BASES,
                              maximum_key_repeats=_SOURCE_ALIGNMENT_CACHE_MAX_REPEATS,
                              maximum_value_bases=_SOURCE_ALIGNMENT_CACHE_MAX_VALUE_BASES),
        }


def validate_candidate_source_reconstruction(
        full_dna, list_repeats, list_repeats_starts, list_spacers):
    """Verify zero-based, half-open candidate geometry against source DNA."""
    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    if not list_repeats:
        raise ValueError("an enhanced candidate must contain at least one repeat")
    if len(list_repeats_starts) != len(list_repeats):
        raise ValueError("every repeat must have exactly one start coordinate")
    if len(list_spacers) != len(list_repeats) - 1:
        raise ValueError("spacer count must be one less than repeat count")
    if any(not isinstance(repeat, str) or not repeat for repeat in list_repeats):
        raise ValueError("repeat sequences must be non-empty strings")
    if any(not isinstance(spacer, str) for spacer in list_spacers):
        raise ValueError("spacer sequences must be strings")

    for index, (start, repeat) in enumerate(
            zip(list_repeats_starts, list_repeats)):
        if not isinstance(start, int) or isinstance(start, bool):
            raise ValueError("repeat starts must be integer offsets")
        end = start + len(repeat)
        if start < 0 or end > len(full_dna):
            raise ValueError("repeat interval is outside the source sequence")
        if full_dna[start:end] != repeat:
            raise ValueError(
                "repeat does not reconstruct from its zero-based half-open interval"
            )

        if index < len(list_spacers):
            next_start = list_repeats_starts[index + 1]
            if next_start < end:
                raise ValueError("repeat intervals overlap or are not monotonic")
            if full_dna[end:next_start] != list_spacers[index]:
                raise ValueError(
                    "spacer does not reconstruct between adjacent repeat intervals"
                )


SOURCE_GLOBAL_PHASE_RULE = {
    "maximum_boundary_shift": 12,
    "minimum_repeats": 3,
    "minimum_repeat_length": 15,
    "maximum_repeat_length": 65,
    "minimum_mean_repeat_length": 21,
    "maximum_mean_repeat_length": 55,
    "minimum_spacer_length": 10,
    "maximum_spacer_length": 140,
    "minimum_column_coverage_fraction": 0.80,
    "spacer_edge_width": 3,
    "background_flank": 5000,
    "minimum_score_improvement": 30.0,
    "minimum_winner_margin": 3.0,
}

# Cross-family arbitration is deliberately narrower than the ordinary phase
# collapse.  Five bases of a neighbouring repeat can otherwise leak into a
# competing spacer phase while the three-base global edge penalty remains
# blind to the last two conserved columns.  This rule is used only after both
# competing source paths have independently passed final rescue eligibility;
# it is never substituted for SOURCE_GLOBAL_PHASE_RULE.
SOURCE_CROSS_FAMILY_PHASE_RULE = dict(
    SOURCE_GLOBAL_PHASE_RULE,
    spacer_edge_width=5,
)

_CANONICAL_DNA_BASES = "ACGT"
_SOURCE_PHASE_BOUNDARY_MINIMUM_CANONICAL_COVERAGE = 0.8
_SOURCE_PHASE_BOUNDARY_MINIMUM_INCLUDED_DOMINANCE = 0.8
_SOURCE_PHASE_BOUNDARY_MINIMUM_DOMINANCE_IMPROVEMENT = 0.2
_SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_REPEAT_COUNT = 8
_SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_MODAL_FRACTION = 0.75
_SOURCE_PHASE_BOUNDARY_RELAXED_MAXIMUM_DISCORDANT_CANONICAL = 2
_SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_MODAL_COUNT_IMPROVEMENT = 4
_SOURCE_PHASE_BOUNDARY_NUMERIC_TOLERANCE = 1e-12


def _validate_source_global_phase_rule(rule):
    if set(rule) != set(SOURCE_GLOBAL_PHASE_RULE):
        raise ValueError("source-backed phase rule fields are incomplete")

    integer_fields = (
        "maximum_boundary_shift",
        "minimum_repeats",
        "minimum_repeat_length",
        "maximum_repeat_length",
        "minimum_mean_repeat_length",
        "maximum_mean_repeat_length",
        "minimum_spacer_length",
        "maximum_spacer_length",
        "spacer_edge_width",
        "background_flank",
    )
    for field in integer_fields:
        value = rule[field]
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError("{} must be an integer".format(field))
        if value < 0:
            raise ValueError("{} must be non-negative".format(field))

    if rule["maximum_boundary_shift"] < 1:
        raise ValueError("maximum_boundary_shift must be positive")
    if rule["minimum_repeat_length"] < 1:
        raise ValueError("minimum_repeat_length must be positive")
    if rule["minimum_repeats"] < 2:
        raise ValueError("minimum_repeats must be at least two")
    if rule["maximum_repeat_length"] < rule["minimum_repeat_length"]:
        raise ValueError("repeat-length bounds are inconsistent")
    if (
            rule["maximum_mean_repeat_length"]
            < rule["minimum_mean_repeat_length"]):
        raise ValueError("mean repeat-length bounds are inconsistent")
    if rule["maximum_spacer_length"] < rule["minimum_spacer_length"]:
        raise ValueError("spacer-length bounds are inconsistent")
    if rule["spacer_edge_width"] < 1:
        raise ValueError("spacer_edge_width must be positive")

    for field in (
            "minimum_column_coverage_fraction",
            "minimum_score_improvement",
            "minimum_winner_margin"):
        value = rule[field]
        if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0):
            raise ValueError("{} must be a finite non-negative number".format(
                field
            ))
    if not 0 < rule["minimum_column_coverage_fraction"] <= 1:
        raise ValueError(
            "minimum_column_coverage_fraction must be in (0, 1]"
        )


def _canonical_base_counts(sequence):
    counts = _canonical_base_count_tuple(sequence)
    return dict(zip(_CANONICAL_DNA_BASES, counts))


def _canonical_base_count_tuple(sequence):
    """Return A/C/G/T counts without retaining normalized observations."""

    adenine = cytosine = guanine = thymine = 0
    for character in sequence:
        base = character.upper()
        if base == "A":
            adenine += 1
        elif base == "C":
            cytosine += 1
        elif base == "G":
            guanine += 1
        elif base == "T":
            thymine += 1
    return adenine, cytosine, guanine, thymine


@lru_cache(maxsize=4096)
def _jeffreys_column_log_bayes_factor_from_counts(
        ordered_counts, ordered_background_probabilities):
    """Score one canonical count vector with bounded exact-value reuse."""

    observation_count = sum(ordered_counts)
    if observation_count < 2:
        return 0.0

    # Keep the historical A/C/G/T accumulation order. Its floating-point bits
    # are part of the candidate-ranking contract.
    log_column_marginal = (
        math.lgamma(2.0)
        - math.lgamma(observation_count + 2.0)
        + sum(
            math.lgamma(count + 0.5) - math.lgamma(0.5)
            for count in ordered_counts
        )
    )
    log_background_likelihood = sum(
        count * math.log(probability)
        for count, probability in zip(
            ordered_counts, ordered_background_probabilities
        )
    )
    return log_column_marginal - log_background_likelihood


def _jeffreys_column_log_bayes_factor(observations, background_probabilities):
    """Compare one repeat column with a fixed, local nucleotide background.

    The repeat-column model integrates an unknown four-base distribution under
    a Dirichlet(1/2, ..., 1/2) Jeffreys prior.  The fixed background model uses
    add-one nucleotide probabilities from the candidate's +/-5 kb source
    window.  Non-ACGT source symbols carry no evidence.
    """

    ordered_counts = _canonical_base_count_tuple(observations)
    if sum(ordered_counts) < 2:
        return 0.0
    ordered_background_probabilities = tuple(
        background_probabilities[base] for base in _CANONICAL_DNA_BASES
    )
    return _jeffreys_column_log_bayes_factor_from_counts(
        ordered_counts, ordered_background_probabilities
    )


def _source_phase_geometry_is_valid(
        repeat_starts, repeat_ends, source_length, rule):
    if (
            not repeat_starts
            or len(repeat_starts) != len(repeat_ends)
            or len(repeat_starts) < rule["minimum_repeats"]):
        return False
    repeat_lengths = tuple(
        end - start for start, end in zip(repeat_starts, repeat_ends)
    )
    spacer_lengths = tuple(
        next_start - end
        for end, next_start in zip(repeat_ends, repeat_starts[1:])
    )
    if (
            min(repeat_starts) < 0
            or max(repeat_ends) > source_length
            or min(repeat_lengths) < rule["minimum_repeat_length"]
            or max(repeat_lengths) > rule["maximum_repeat_length"]):
        return False
    mean_repeat_length = sum(repeat_lengths) / len(repeat_lengths)
    if not (
            rule["minimum_mean_repeat_length"]
            <= mean_repeat_length
            <= rule["maximum_mean_repeat_length"]):
        return False
    return bool(spacer_lengths) and all(
        rule["minimum_spacer_length"]
        <= spacer_length
        <= rule["maximum_spacer_length"]
        for spacer_length in spacer_lengths
    )


def build_source_backed_candidate_from_boundaries(
        full_dna, crispr_candidate, repeat_starts, repeat_ends):
    """Rebuild one candidate solely from ordered half-open source intervals.

    Repeat/spacer nucleotides, coordinates, and the deterministic pairwise
    alignment all come from ``full_dna``. Only competition ancestry is inherited
    from an explicit parent; it never determines sequence or alignment.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    repeat_starts = tuple(repeat_starts)
    repeat_ends = tuple(repeat_ends)
    if not repeat_starts or len(repeat_starts) != len(repeat_ends):
        raise ValueError(
            "source-backed repeat starts and ends must have equal nonzero length"
        )
    if any(
            not isinstance(value, int) or isinstance(value, bool)
            for value in repeat_starts + repeat_ends):
        raise ValueError("source-backed repeat boundaries must be integers")
    for index, (start, end) in enumerate(zip(repeat_starts, repeat_ends)):
        if start < 0 or end <= start or end > len(full_dna):
            raise ValueError("source-backed repeat interval is invalid")
        if index and start < repeat_ends[index - 1]:
            raise ValueError(
                "source-backed repeat intervals overlap or are not monotonic"
            )

    repeats = [
        full_dna[start:end]
        for start, end in zip(repeat_starts, repeat_ends)
    ]
    spacers = [
        full_dna[end:next_start]
        for end, next_start in zip(repeat_ends, repeat_starts[1:])
    ]

    # Never reuse the legacy alignment or consensus.  They are representation
    # choices that can differ for the same source geometry and therefore cannot
    # be allowed to change score-relevant rebuilt candidate state.
    # Copy the immutable pure alignment into a fresh candidate. Coordinates,
    # spacers, score-relevant mutable state and lineage are never cached.
    gapped_repeats = list(_source_repeat_alignment(repeats))
    resolved = CrisprCandidate(
        list_repeats=repeats,
        list_repeats_gaped=gapped_repeats,
        list_spacers=spacers,
        list_repeat_starts=list(repeat_starts),
    )
    if crispr_candidate is not None:
        inherit_candidate_origins(resolved, (crispr_candidate,))
    validate_candidate_source_reconstruction(
        full_dna,
        resolved.list_repeats,
        resolved.list_repeat_starts,
        resolved.list_spacers,
    )
    return resolved


def resolve_source_backed_global_phase(
        full_dna, crispr_candidate, rule=None):
    """Conservatively resolve a uniform repeat-boundary phase from source DNA.

    A hypothesis ``(a, b)`` moves every repeat's left boundary by ``a`` and
    right boundary by ``b`` bases.  All pairs in the inclusive configured
    range are evaluated without comparator coordinates or classifier output.
    The score is repeat-column log Bayes-factor evidence minus positive
    repeat-like evidence in the first and last spacer-edge columns.

    The original object is returned unless the unique winner improves on the
    current ``(0, 0)`` geometry by at least 30 natural-log units and beats the
    runner-up by at least 3.  Metadata makes every abstention auditable.
    """

    phase_rule = dict(SOURCE_GLOBAL_PHASE_RULE if rule is None else rule)
    _validate_source_global_phase_rule(phase_rule)
    metadata = {
        "status": "abstained",
        "reason": None,
        "left_shift": 0,
        "right_shift": 0,
        "current_score": None,
        "winner_score": None,
        "score_improvement": None,
        "winner_margin": None,
        "hypotheses_scored": 0,
        "background_window": None,
        "background_probabilities": None,
    }
    if full_dna is None or full_dna == "":
        metadata.update(status="unchanged", reason="no_source_dna")
        return crispr_candidate, metadata
    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")

    validate_candidate_source_reconstruction(
        full_dna,
        crispr_candidate.list_repeats,
        crispr_candidate.list_repeat_starts,
        crispr_candidate.list_spacers,
    )
    repeat_count = len(crispr_candidate.list_repeats)
    if repeat_count < phase_rule["minimum_repeats"]:
        metadata["reason"] = "insufficient_repeats"
        return crispr_candidate, metadata

    starts = tuple(crispr_candidate.list_repeat_starts)
    lengths = tuple(len(repeat) for repeat in crispr_candidate.list_repeats)
    ends = tuple(start + length for start, length in zip(starts, lengths))

    flank = phase_rule["background_flank"]
    background_start = max(0, starts[0] - flank)
    background_end = min(len(full_dna), ends[-1] + flank)
    background_counts = _canonical_base_counts(
        full_dna[background_start:background_end]
    )
    background_size = sum(background_counts.values())
    background_probabilities = {
        base: (background_counts[base] + 1.0) / (background_size + 4.0)
        for base in _CANONICAL_DNA_BASES
    }
    metadata.update({
        "background_window": (background_start, background_end),
        "background_probabilities": dict(background_probabilities),
    })

    coverage_fraction = phase_rule["minimum_column_coverage_fraction"]
    repeat_column_coverage = math.ceil(coverage_fraction * repeat_count)
    spacer_count = repeat_count - 1
    spacer_column_coverage = math.ceil(coverage_fraction * spacer_count)

    def covered_column_score(observations, minimum_coverage):
        canonical = [
            character.upper() for character in observations
            if character.upper() in _CANONICAL_DNA_BASES
        ]
        if len(canonical) < minimum_coverage:
            return 0.0
        return _jeffreys_column_log_bayes_factor(
            canonical, background_probabilities
        )

    maximum_shift = phase_rule["maximum_boundary_shift"]
    shifts = tuple(range(-maximum_shift, maximum_shift + 1))
    relative_begin = -maximum_shift
    relative_end = max(lengths) + maximum_shift
    repeat_prefixes = {}
    for right_shift in shifts:
        prefix = [0.0]
        for relative_position in range(relative_begin, relative_end):
            observations = []
            for start, length in zip(starts, lengths):
                position = start + relative_position
                if (
                        relative_position < length + right_shift
                        and 0 <= position < len(full_dna)):
                    observations.append(full_dna[position])
            score = covered_column_score(
                observations, repeat_column_coverage
            )
            prefix.append(prefix[-1] + score)
        repeat_prefixes[right_shift] = prefix

    edge_width = phase_rule["spacer_edge_width"]
    right_edge_penalties = {}
    left_edge_penalties = {}
    for shift in shifts:
        right_edge_penalties[shift] = sum(
            max(0.0, covered_column_score(
                (full_dna[end + shift + depth] for end in ends[:-1]
                 if 0 <= end + shift + depth < len(full_dna)),
                spacer_column_coverage,
            ))
            for depth in range(edge_width)
        )
        left_edge_penalties[shift] = sum(
            max(0.0, covered_column_score(
                (full_dna[start + shift - 1 - depth] for start in starts[1:]
                 if 0 <= start + shift - 1 - depth < len(full_dna)),
                spacer_column_coverage,
            ))
            for depth in range(edge_width)
        )

    hypotheses = []
    for left_shift in shifts:
        proposed_starts = tuple(start + left_shift for start in starts)
        for right_shift in shifts:
            proposed_ends = tuple(end + right_shift for end in ends)
            if not _source_phase_geometry_is_valid(
                    proposed_starts,
                    proposed_ends,
                    len(full_dna),
                    phase_rule):
                continue

            relative_stop = max(lengths) + right_shift
            repeat_prefix = repeat_prefixes[right_shift]
            repeat_evidence = (
                repeat_prefix[relative_stop - relative_begin]
                - repeat_prefix[left_shift - relative_begin]
            )
            score = (
                repeat_evidence
                - right_edge_penalties[right_shift]
                - left_edge_penalties[left_shift]
            )
            hypotheses.append((score, left_shift, right_shift))

    metadata["hypotheses_scored"] = len(hypotheses)
    current = next(
        (hypothesis for hypothesis in hypotheses
         if hypothesis[1:] == (0, 0)),
        None,
    )
    if current is None or len(hypotheses) < 2:
        metadata["reason"] = "insufficient_valid_hypotheses"
        return crispr_candidate, metadata

    hypotheses.sort(key=lambda hypothesis: (
        -hypothesis[0],
        abs(hypothesis[1]) + abs(hypothesis[2]),
        abs(hypothesis[1]),
        abs(hypothesis[2]),
        -hypothesis[1],
        -hypothesis[2],
    ))
    winner = hypotheses[0]
    runner_up = hypotheses[1]
    improvement = winner[0] - current[0]
    margin = winner[0] - runner_up[0]
    metadata.update({
        "left_shift": winner[1],
        "right_shift": winner[2],
        "current_score": current[0],
        "winner_score": winner[0],
        "score_improvement": improvement,
        "winner_margin": margin,
    })
    if winner[1:] == (0, 0):
        metadata["reason"] = "current_geometry_wins"
        return crispr_candidate, metadata
    if improvement < phase_rule["minimum_score_improvement"]:
        metadata["reason"] = "insufficient_score_improvement"
        return crispr_candidate, metadata
    if margin < phase_rule["minimum_winner_margin"]:
        metadata["reason"] = "ambiguous_winner"
        return crispr_candidate, metadata

    proposed_starts = tuple(start + winner[1] for start in starts)
    proposed_ends = tuple(end + winner[2] for end in ends)
    resolved = build_source_backed_candidate_from_boundaries(
        full_dna, crispr_candidate, proposed_starts, proposed_ends
    )
    metadata.update(status="committed", reason="strong_unique_source_evidence")
    return resolved, metadata




def _candidate_source_boundaries(crispr_candidate):
    starts = tuple(crispr_candidate.list_repeat_starts)
    repeats = tuple(crispr_candidate.list_repeats)
    if not starts or len(starts) != len(repeats):
        raise ValueError(
            "every phase-collapse repeat must have one source start"
        )
    ends = tuple(
        start + len(repeat)
        for start, repeat in zip(starts, repeats)
    )
    return starts, ends


def _bounded_phase_equivalent(first, second, maximum_shift):
    """Return whether two paths describe the same bounded repeat occurrences.

    A real fuzzy-search locus often contains locally different boundaries for
    the same ordered occurrences: one repeat can carry an insertion or one
    alignment can leave a few conserved bases in the adjacent spacer.  A
    constant global shift signature cannot compare those hypotheses.  Equal
    repeat count plus a per-occurrence bound does compare them, while still
    keeping a missing terminal repeat or a genuinely different repeat lattice
    separate.
    """

    return _bounded_source_boundaries_equivalent(
        _candidate_source_boundaries(first),
        _candidate_source_boundaries(second),
        maximum_shift,
    )


def _bounded_source_boundaries_equivalent(
        first_boundaries, second_boundaries, maximum_shift):
    """Compare two precomputed source geometries under the phase rule."""

    first_starts, first_ends = first_boundaries
    second_starts, second_ends = second_boundaries
    return (
        len(first_starts) == len(second_starts)
        and all(
            abs(first_start - second_start) <= maximum_shift
            and abs(first_end - second_end) <= maximum_shift
            for first_start, second_start, first_end, second_end in zip(
                first_starts,
                second_starts,
                first_ends,
                second_ends,
            )
        )
    )


def _phase_boundaries_fit_complete_link_cluster(
        boundaries, cluster, maximum_shift):
    """Return whether one geometry is bounded to every cluster member.

    For each coordinate, being within ``maximum_shift`` of every existing
    value is exactly equivalent to lying between ``maximum - shift`` and
    ``minimum + shift``.  The envelope therefore preserves complete-link
    membership without revisiting every member geometry.
    """

    starts, ends = boundaries
    minimum_starts = cluster["minimum_starts"]
    maximum_starts = cluster["maximum_starts"]
    minimum_ends = cluster["minimum_ends"]
    maximum_ends = cluster["maximum_ends"]
    return (
        len(starts) == len(minimum_starts)
        and len(ends) == len(minimum_ends)
        and all(
            maximum_start - maximum_shift
            <= start
            <= minimum_start + maximum_shift
            and maximum_end - maximum_shift
            <= end
            <= minimum_end + maximum_shift
            for start, end, minimum_start, maximum_start,
            minimum_end, maximum_end in zip(
                starts,
                ends,
                minimum_starts,
                maximum_starts,
                minimum_ends,
                maximum_ends,
            )
        )
    )


def _extend_complete_link_phase_cluster(
        cluster, indexed_candidate, boundaries):
    """Append one accepted geometry and update its complete-link envelope."""

    starts, ends = boundaries
    cluster["members"].append(indexed_candidate)
    cluster["minimum_starts"] = tuple(map(
        min, cluster["minimum_starts"], starts
    ))
    cluster["maximum_starts"] = tuple(map(
        max, cluster["maximum_starts"], starts
    ))
    cluster["minimum_ends"] = tuple(map(
        min, cluster["minimum_ends"], ends
    ))
    cluster["maximum_ends"] = tuple(map(
        max, cluster["maximum_ends"], ends
    ))


def _complete_link_phase_clusters(indexed_candidates, maximum_shift):
    """Build deterministic bounded cliques without transitive chaining."""

    ordered = sorted(
        (
            (indexed_candidate, _candidate_source_boundaries(
                indexed_candidate[1]))
            for indexed_candidate in indexed_candidates
        ),
        key=lambda item: (
            item[1],
            item[0][0],
        ),
    )
    clusters = []
    for indexed_candidate, boundaries in ordered:
        for cluster in clusters:
            if _phase_boundaries_fit_complete_link_cluster(
                    boundaries, cluster, maximum_shift):
                _extend_complete_link_phase_cluster(
                    cluster, indexed_candidate, boundaries
                )
                break
        else:
            starts, ends = boundaries
            clusters.append({
                "members": [indexed_candidate],
                "minimum_starts": starts,
                "maximum_starts": starts,
                "minimum_ends": ends,
                "maximum_ends": ends,
            })
    return [cluster["members"] for cluster in clusters]


def _source_phase_background(full_dna, candidates, flank):
    all_starts = []
    all_ends = []
    for candidate in candidates:
        starts, ends = _candidate_source_boundaries(candidate)
        all_starts.extend(starts)
        all_ends.extend(ends)
    background_start = max(0, min(all_starts) - flank)
    background_end = min(len(full_dna), max(all_ends) + flank)
    counts = _canonical_base_counts(
        full_dna[background_start:background_end]
    )
    canonical_count = sum(counts.values())
    probabilities = {
        base: (counts[base] + 1.0) / (canonical_count + 4.0)
        for base in _CANONICAL_DNA_BASES
    }
    return (
        (background_start, background_end),
        probabilities,
    )


def _source_phase_geometry_score(
        full_dna, repeat_starts, repeat_ends, rule,
        background_probabilities):
    """Score one absolute geometry against one already-frozen background."""

    if not _source_phase_geometry_is_valid(
            repeat_starts, repeat_ends, len(full_dna), rule):
        return None

    repeat_count = len(repeat_starts)
    repeat_lengths = tuple(
        end - start for start, end in zip(repeat_starts, repeat_ends)
    )
    repeat_coverage = math.ceil(
        rule["minimum_column_coverage_fraction"] * repeat_count
    )
    spacer_count = repeat_count - 1
    spacer_coverage = math.ceil(
        rule["minimum_column_coverage_fraction"] * spacer_count
    )
    ordered_background_probabilities = tuple(
        background_probabilities[base] for base in _CANONICAL_DNA_BASES
    )

    def covered_column_score(observations, minimum_coverage):
        ordered_counts = _canonical_base_count_tuple(observations)
        if sum(ordered_counts) < minimum_coverage:
            return 0.0
        return _jeffreys_column_log_bayes_factor_from_counts(
            ordered_counts, ordered_background_probabilities
        )

    repeat_evidence = 0.0
    for relative_position in range(max(repeat_lengths)):
        repeat_evidence += covered_column_score(
            (
                full_dna[start + relative_position]
                for start, length in zip(repeat_starts, repeat_lengths)
                if relative_position < length
            ),
            repeat_coverage,
        )

    edge_width = rule["spacer_edge_width"]
    right_edge_penalty = sum(
        max(0.0, covered_column_score(
            (
                full_dna[end + depth] for end in repeat_ends[:-1]
                if 0 <= end + depth < len(full_dna)
            ),
            spacer_coverage,
        ))
        for depth in range(edge_width)
    )
    left_edge_penalty = sum(
        max(0.0, covered_column_score(
            (
                full_dna[start - 1 - depth] for start in repeat_starts[1:]
                if 0 <= start - 1 - depth < len(full_dna)
            ),
            spacer_coverage,
        ))
        for depth in range(edge_width)
    )
    return repeat_evidence - right_edge_penalty - left_edge_penalty


def _source_column_modal_dominance(full_dna, positions):
    """Measure modal A/C/G/T agreement at homologous source positions."""

    positions = tuple(positions)
    observations = [
        full_dna[position].upper()
        for position in positions
        if (
            0 <= position < len(full_dna)
            and full_dna[position].upper() in _CANONICAL_DNA_BASES
        )
    ]
    counts = _canonical_base_counts(observations)
    canonical_count = sum(counts.values())
    modal_count = max(counts.values()) if canonical_count else 0
    dominance = (
        modal_count / canonical_count
        if canonical_count else 0.0
    )
    return {
        "dominance": dominance,
        "modal_count": modal_count,
        "canonical_count": canonical_count,
        "position_count": len(positions),
    }


def _source_phase_boundary_evidence(
        full_dna, winner_geometry, member_geometries):
    """Require a sharp source-column contrast at each selected boundary.

    A side is checked only when any corresponding boundary varies among group
    members.  Both compared columns must contain canonical bases for at least
    80% of repeats and the selected column must exceed the immediately
    adjacent excluded column by at least 0.2 dominance. The selected column
    normally requires at least 0.8 modal agreement. A finite-count branch can
    admit a lower fraction only with at least eight repeats, at least 75%
    modal support, at most two discordant canonical observations, and an
    absolute modal-count advantage of at least four over the excluded column.
    A fixed 1e-12 tolerance prevents an exact decimal threshold from failing
    solely because of binary floating-point representation.
    """

    winner_starts, winner_ends = winner_geometry
    member_geometries = tuple(member_geometries)
    tolerance = _SOURCE_PHASE_BOUNDARY_NUMERIC_TOLERANCE

    def side_evidence(side):
        if side == "left":
            varies = len({
                geometry[0] for geometry in member_geometries
            }) > 1
            included_positions = winner_starts
            excluded_positions = tuple(
                start - 1 for start in winner_starts
            )
        else:
            varies = len({
                geometry[1] for geometry in member_geometries
            }) > 1
            included_positions = tuple(end - 1 for end in winner_ends)
            excluded_positions = winner_ends

        included = _source_column_modal_dominance(
            full_dna, included_positions
        )
        excluded = _source_column_modal_dominance(
            full_dna, excluded_positions
        )
        improvement = included["dominance"] - excluded["dominance"]
        included_coverage = (
            included["canonical_count"] / included["position_count"]
            if included["position_count"] else 0.0
        )
        excluded_coverage = (
            excluded["canonical_count"] / excluded["position_count"]
            if excluded["position_count"] else 0.0
        )
        repeat_count = included["position_count"]
        relaxed_minimum_modal_count = math.ceil(
            _SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_MODAL_FRACTION
            * repeat_count
        )
        included_discordant_canonical_count = (
            included["canonical_count"] - included["modal_count"]
        )
        modal_count_improvement = (
            included["modal_count"] - excluded["modal_count"]
        )
        standard_included_dominance_passes = (
            included["dominance"] + tolerance
            >= _SOURCE_PHASE_BOUNDARY_MINIMUM_INCLUDED_DOMINANCE
        )
        relaxed_finite_count_support_passes = (
            repeat_count
            >= _SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_REPEAT_COUNT
            and included["modal_count"] >= relaxed_minimum_modal_count
            and included_discordant_canonical_count
            <= _SOURCE_PHASE_BOUNDARY_RELAXED_MAXIMUM_DISCORDANT_CANONICAL
            and modal_count_improvement
            >= (
                _SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_MODAL_COUNT_IMPROVEMENT
            )
        )
        canonical_coverage_passes = (
            included_coverage + tolerance
            >= _SOURCE_PHASE_BOUNDARY_MINIMUM_CANONICAL_COVERAGE
            and excluded_coverage + tolerance
            >= _SOURCE_PHASE_BOUNDARY_MINIMUM_CANONICAL_COVERAGE
        )
        dominance_improvement_passes = (
            improvement + tolerance
            >= _SOURCE_PHASE_BOUNDARY_MINIMUM_DOMINANCE_IMPROVEMENT
        )
        passes = (
            not varies
            or (
                canonical_coverage_passes
                and dominance_improvement_passes
                and (
                    standard_included_dominance_passes
                    or relaxed_finite_count_support_passes
                )
            )
        )
        if not varies:
            evidence_branch = "not_required"
        elif (
                canonical_coverage_passes
                and dominance_improvement_passes
                and standard_included_dominance_passes):
            evidence_branch = "standard_dominance"
        elif (
                canonical_coverage_passes
                and dominance_improvement_passes
                and relaxed_finite_count_support_passes):
            evidence_branch = "relaxed_finite_count"
        else:
            evidence_branch = None
        return {
            "varies": varies,
            "included": included,
            "immediately_excluded": excluded,
            "included_canonical_coverage": included_coverage,
            "immediately_excluded_canonical_coverage": excluded_coverage,
            "dominance_improvement": improvement,
            "modal_count_improvement": modal_count_improvement,
            "included_discordant_canonical_count": (
                included_discordant_canonical_count
            ),
            "relaxed_minimum_included_modal_count": (
                relaxed_minimum_modal_count
            ),
            "canonical_coverage_passes": canonical_coverage_passes,
            "dominance_improvement_passes": (
                dominance_improvement_passes
            ),
            "standard_included_dominance_passes": (
                standard_included_dominance_passes
            ),
            "relaxed_finite_count_support_passes": (
                relaxed_finite_count_support_passes
            ),
            "evidence_branch": evidence_branch,
            "passes": passes,
        }

    left = side_evidence("left")
    right = side_evidence("right")
    return {
        "minimum_canonical_coverage": (
            _SOURCE_PHASE_BOUNDARY_MINIMUM_CANONICAL_COVERAGE
        ),
        "minimum_included_dominance": (
            _SOURCE_PHASE_BOUNDARY_MINIMUM_INCLUDED_DOMINANCE
        ),
        "minimum_dominance_improvement": (
            _SOURCE_PHASE_BOUNDARY_MINIMUM_DOMINANCE_IMPROVEMENT
        ),
        "relaxed_finite_count_rule": {
            "minimum_repeat_count": (
                _SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_REPEAT_COUNT
            ),
            "minimum_modal_fraction": (
                _SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_MODAL_FRACTION
            ),
            "maximum_discordant_canonical_observations": (
                _SOURCE_PHASE_BOUNDARY_RELAXED_MAXIMUM_DISCORDANT_CANONICAL
            ),
            "minimum_modal_count_improvement": (
                _SOURCE_PHASE_BOUNDARY_RELAXED_MINIMUM_MODAL_COUNT_IMPROVEMENT
            ),
            "minimum_dominance_improvement": (
                _SOURCE_PHASE_BOUNDARY_MINIMUM_DOMINANCE_IMPROVEMENT
            ),
        },
        "numeric_tolerance": tolerance,
        "left": left,
        "right": right,
        "passes": left["passes"] and right["passes"],
    }


def compare_source_backed_phase_candidates(full_dna, candidates, rule=None):
    """Return a source-only winner for one bounded phase competition.

    This primitive neither scores classifiers nor changes candidate/category
    membership.  Callers must first establish that the candidates represent
    the same ordered source occurrences.  Every candidate is independently
    reconstructed, all geometries are evaluated against one frozen source
    background, and a winner is returned only when its raw source score is
    uniquely best, both source-boundary gates pass, and its margin reaches the
    named rule threshold.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    try:
        candidates = tuple(candidates)
    except TypeError as error:
        raise ValueError("phase competitors must be iterable") from error
    if len(candidates) < 2:
        raise ValueError("phase competition requires at least two candidates")

    phase_rule = dict(
        SOURCE_CROSS_FAMILY_PHASE_RULE if rule is None else rule
    )
    _validate_source_global_phase_rule(phase_rule)
    geometries = []
    for candidate in candidates:
        validate_candidate_source_reconstruction(
            full_dna,
            candidate.list_repeats,
            candidate.list_repeat_starts,
            candidate.list_spacers,
        )
        geometries.append(_candidate_source_boundaries(candidate))

    repeat_counts = {len(geometry[0]) for geometry in geometries}
    if len(repeat_counts) != 1:
        raise ValueError(
            "phase competitors must contain the same matched occurrence count"
        )
    if next(iter(repeat_counts)) < phase_rule["minimum_repeats"]:
        raise ValueError("phase competition has too few matched occurrences")
    maximum_shift = phase_rule["maximum_boundary_shift"]
    for first_index, first in enumerate(candidates[:-1]):
        for second in candidates[first_index + 1:]:
            if not _bounded_phase_equivalent(
                    first, second, maximum_shift):
                raise ValueError(
                    "phase competitors are outside the bounded source rule"
                )

    distinct_geometries = tuple(sorted(set(geometries)))
    if len(distinct_geometries) < 2:
        return {
            "status": "abstained",
            "reason": "identical_source_geometry",
            "winner_candidate": None,
            "winner_candidate_index": None,
            "winner_margin": 0.0,
            "background_window": None,
            "background_probabilities": None,
            "candidate_evidence": (),
            "classifier_used": False,
            "rule": dict(phase_rule),
        }

    background_window, background_probabilities = _source_phase_background(
        full_dna,
        candidates,
        phase_rule["background_flank"],
    )
    boundary_by_geometry = {
        geometry: _source_phase_boundary_evidence(
            full_dna,
            geometry,
            distinct_geometries,
        )
        for geometry in distinct_geometries
    }
    evidence = []
    for candidate_index, (candidate, geometry) in enumerate(zip(
            candidates, geometries)):
        score = _source_phase_geometry_score(
            full_dna,
            geometry[0],
            geometry[1],
            phase_rule,
            background_probabilities,
        )
        evidence.append({
            "candidate": candidate,
            "candidate_index": candidate_index,
            "geometry": geometry,
            "source_score": score,
            "boundary_evidence": boundary_by_geometry[geometry],
        })

    if any(record["source_score"] is None for record in evidence):
        return {
            "status": "abstained",
            "reason": "geometry_outside_source_phase_rule",
            "winner_candidate": None,
            "winner_candidate_index": None,
            "winner_margin": None,
            "background_window": background_window,
            "background_probabilities": dict(background_probabilities),
            "candidate_evidence": tuple(evidence),
            "classifier_used": False,
            "rule": dict(phase_rule),
        }

    ranked = sorted(
        evidence,
        key=lambda record: (
            -record["source_score"],
            record["geometry"],
            record["candidate_index"],
        ),
    )
    winner = ranked[0]
    runner_up = ranked[1]
    margin = winner["source_score"] - runner_up["source_score"]
    if margin <= 1e-9:
        reason = "source_score_tie"
    elif not winner["boundary_evidence"]["passes"]:
        reason = "weak_source_boundary_evidence"
    elif margin < phase_rule["minimum_winner_margin"]:
        reason = "weak_source_score_margin"
    else:
        reason = "unique_source_phase_winner"

    won = reason == "unique_source_phase_winner"
    return {
        "status": "selected" if won else "abstained",
        "reason": reason,
        "winner_candidate": winner["candidate"] if won else None,
        "winner_candidate_index": (
            winner["candidate_index"] if won else None
        ),
        "winner_margin": margin,
        "background_window": background_window,
        "background_probabilities": dict(background_probabilities),
        "candidate_evidence": tuple(evidence),
        "classifier_used": False,
        "rule": dict(phase_rule),
    }


def collapse_source_backed_phase_equivalent_candidates(
        full_dna, candidates_by_locus, rule=None,
        skip_source_scores_when_boundaries_reject=False):
    """Collapse only bounded source-equivalent boundary representations.

    Classifier output is intentionally unavailable to this function. Within
    one legacy locus, candidates compete only when they contain the same number
    of ordered repeat occurrences and every corresponding start/end differs by
    at most the audited boundary-search radius. Every distinct input geometry
    is scored once against one source background frozen from the union of that
    group. Geometries without sharp source-column boundary contrast are first
    excluded; the remaining winner must have a strictly positive,
    non-negligible source-score margin.

    The emitted candidate is rebuilt from zero-based, half-open FASTA slices;
    provenance retains every replaced runtime identity and geometry for later reporting
    or transactional rollback.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    if not isinstance(candidates_by_locus, dict):
        raise ValueError("phase-collapse candidates must be a mapping")
    if not isinstance(skip_source_scores_when_boundaries_reject, bool):
        raise ValueError("phase-collapse score gate must be boolean")

    phase_rule = dict(SOURCE_GLOBAL_PHASE_RULE if rule is None else rule)
    _validate_source_global_phase_rule(phase_rule)
    maximum_shift = phase_rule["maximum_boundary_shift"]
    stats = {
        "enabled": True,
        "groups_considered": 0,
        "groups_collapsed": 0,
        "groups_ambiguous": 0,
        "groups_outside_phase_rule": 0,
        "groups_weak_boundary_evidence": 0,
        "groups_weak_score_margin": 0,
        "candidates_removed": 0,
        "candidates_input": sum(
            len(candidates) for candidates in candidates_by_locus.values()
        ),
        "candidates_output": 0,
    }
    provenance = []
    collapsed_by_locus = {}

    for locus, candidates in candidates_by_locus.items():
        candidates = list(candidates)
        for candidate in candidates:
            validate_candidate_source_reconstruction(
                full_dna,
                candidate.list_repeats,
                candidate.list_repeat_starts,
                candidate.list_spacers,
            )

        by_repeat_count = {}
        for candidate_index, candidate in enumerate(candidates):
            # Rebuilding a phase winner cannot erase a broken-chain or other
            # classification restriction. Restricted candidates remain intact.
            eligible = getattr(candidate, "source_path_classification_eligible", True)
            if type(eligible) is not bool or not eligible or any(
                    getattr(record, "derivation", None) == "local_chain_family_parent"
                    for record in getattr(candidate, "source_path_provenance", ())):
                continue
            repeat_count = len(candidate.list_repeat_starts)
            by_repeat_count.setdefault(repeat_count, []).append(
                (candidate_index, candidate)
            )

        groups = []
        for repeat_count, indexed_candidates in by_repeat_count.items():
            if repeat_count < phase_rule["minimum_repeats"]:
                continue
            for cluster in _complete_link_phase_clusters(
                    indexed_candidates, maximum_shift):
                distinct_geometries = {
                    _candidate_source_boundaries(candidate)
                    for _index, candidate in cluster
                }
                # Exact source geometry cannot distinguish competing gapped
                # alignments. Preserve such duplicates until a later stage
                # has evidence that can choose among them without silently
                # changing classifier features.
                if len(cluster) < 2 or len(distinct_geometries) < 2:
                    continue
                groups.append(((repeat_count,), cluster))

        replacements = {}
        removed_indexes = set()
        for signature, cluster in groups:
            stats["groups_considered"] += 1
            member_candidates = tuple(
                candidate for _index, candidate in cluster
            )
            member_geometries = tuple(dict.fromkeys(
                _candidate_source_boundaries(candidate)
                for candidate in member_candidates
            ))
            boundary_evidence_by_geometry = None
            if skip_source_scores_when_boundaries_reject:
                boundary_evidence_by_geometry = {
                    geometry: _source_phase_boundary_evidence(
                        full_dna,
                        geometry,
                        member_geometries,
                    )
                    for geometry in member_geometries
                }
                if not any(
                        evidence["passes"]
                        for evidence in boundary_evidence_by_geometry.values()):
                    stats["groups_weak_boundary_evidence"] += 1
                    provenance.append({
                        "status": "unchanged",
                        "reason": "weak_source_boundary_evidence",
                        "locus": locus,
                        "architecture_signature": signature,
                        "background_window": None,
                        "background_probabilities": None,
                        "member_count": len(member_candidates),
                        "member_candidate_indexes": tuple(
                            index for index, _candidate in cluster
                        ),
                        "member_candidate_ids": tuple(
                            id(candidate) for candidate in member_candidates
                        ),
                        "member_geometries": member_geometries,
                        "source_scores": (),
                        "winner_margin": None,
                        "boundary_evidence_by_geometry": tuple(
                            (geometry, boundary_evidence_by_geometry[geometry])
                            for geometry in member_geometries
                        ),
                        "source_scoring_skipped": True,
                        "classifier_used": False,
                    })
                    continue
            background_window, background_probabilities = (
                _source_phase_background(
                    full_dna,
                    member_candidates,
                    phase_rule["background_flank"],
                )
            )
            geometry_scores = {}
            for candidate in member_candidates:
                geometry = _candidate_source_boundaries(candidate)
                geometry_scores.setdefault(
                    geometry,
                    _source_phase_geometry_score(
                        full_dna,
                        geometry[0],
                        geometry[1],
                        phase_rule,
                        background_probabilities,
                    ),
                )

            if any(score is None for score in geometry_scores.values()):
                stats["groups_outside_phase_rule"] += 1
                provenance.append({
                    "status": "unchanged",
                    "reason": "geometry_outside_source_phase_rule",
                    "locus": locus,
                    "architecture_signature": signature,
                    "background_window": background_window,
                    "member_count": len(member_candidates),
                    "member_candidate_indexes": tuple(
                        index for index, _candidate in cluster
                    ),
                    "member_candidate_ids": tuple(
                        id(candidate) for candidate in member_candidates
                    ),
                    "member_geometries": tuple(geometry_scores),
                    "classifier_used": False,
                })
                continue

            all_ranked = sorted(
                geometry_scores.items(),
                key=lambda item: (
                    -item[1],
                    item[0][0],
                    item[0][1],
                ),
            )
            all_score_margin = all_ranked[0][1] - all_ranked[1][1]
            if all_score_margin <= 1e-9:
                stats["groups_ambiguous"] += 1
                provenance.append({
                    "status": "unchanged",
                    "reason": "source_score_tie",
                    "locus": locus,
                    "architecture_signature": signature,
                    "background_window": background_window,
                    "background_probabilities": dict(
                        background_probabilities
                    ),
                    "member_count": len(member_candidates),
                    "member_candidate_indexes": tuple(
                        index for index, _candidate in cluster
                    ),
                    "member_candidate_ids": tuple(
                        id(candidate) for candidate in member_candidates
                    ),
                    "member_geometries": tuple(geometry_scores),
                    "source_scores": tuple(all_ranked),
                    "winner_margin": all_score_margin,
                    "classifier_used": False,
                })
                continue

            if boundary_evidence_by_geometry is None:
                boundary_evidence_by_geometry = {
                    geometry: _source_phase_boundary_evidence(
                        full_dna,
                        geometry,
                        tuple(geometry_scores),
                    )
                    for geometry in geometry_scores
                }
            boundary_supported_ranked = [
                item for item in all_ranked
                if boundary_evidence_by_geometry[item[0]]["passes"]
            ]
            winner_geometry, winner_score = all_ranked[0]
            boundary_evidence = boundary_evidence_by_geometry[
                winner_geometry
            ]
            if not boundary_evidence["passes"]:
                stats["groups_weak_boundary_evidence"] += 1
                provenance.append({
                    "status": "unchanged",
                    "reason": "weak_source_boundary_evidence",
                    "locus": locus,
                    "architecture_signature": signature,
                    "background_window": background_window,
                    "background_probabilities": dict(
                        background_probabilities
                    ),
                    "member_count": len(member_candidates),
                    "member_candidate_indexes": tuple(
                        index for index, _candidate in cluster
                    ),
                    "member_candidate_ids": tuple(
                        id(candidate) for candidate in member_candidates
                    ),
                    "member_geometries": tuple(geometry_scores),
                    "source_scores": tuple(all_ranked),
                    "proposed_canonical_geometry": winner_geometry,
                    "proposed_canonical_source_score": winner_score,
                    "winner_margin": all_score_margin,
                    "boundary_evidence": boundary_evidence,
                    "boundary_evidence_by_geometry": tuple(
                        (geometry, boundary_evidence_by_geometry[geometry])
                        for geometry, _score in all_ranked
                    ),
                    "classifier_used": False,
                })
                continue

            winner_margin = all_score_margin
            if winner_margin < phase_rule["minimum_winner_margin"]:
                stats["groups_weak_score_margin"] += 1
                provenance.append({
                    "status": "unchanged",
                    "reason": "weak_source_score_margin",
                    "locus": locus,
                    "architecture_signature": signature,
                    "background_window": background_window,
                    "background_probabilities": dict(
                        background_probabilities
                    ),
                    "member_count": len(member_candidates),
                    "member_candidate_indexes": tuple(
                        index for index, _candidate in cluster
                    ),
                    "member_candidate_ids": tuple(
                        id(candidate) for candidate in member_candidates
                    ),
                    "member_geometries": tuple(geometry_scores),
                    "source_scores": tuple(
                        (geometry, score)
                        for geometry, score in all_ranked
                    ),
                    "source_boundary_supported_scores": tuple(
                        boundary_supported_ranked
                    ),
                    "winner_margin": winner_margin,
                    "boundary_evidence": boundary_evidence,
                    "minimum_winner_margin": phase_rule[
                        "minimum_winner_margin"
                    ],
                    "classifier_used": False,
                })
                continue

            winner_members = [
                (candidate_index, candidate)
                for candidate_index, candidate in cluster
                if _candidate_source_boundaries(candidate) == winner_geometry
            ]
            _winner_index, winner_seed = min(
                winner_members,
                key=lambda item: (
                    item[1].consensus,
                    item[1].consensus_gaped,
                    tuple(item[1].list_repeats),
                    tuple(item[1].list_spacers),
                    item[0],
                ),
            )
            canonical = build_source_backed_candidate_from_boundaries(
                full_dna,
                winner_seed,
                winner_geometry[0],
                winner_geometry[1],
            )
            inherit_candidate_origins(canonical, (member for _, member in winner_members))
            repair_records = []
            for _index, member in winner_members:
                for record in getattr(member, "boundary_repair_records", ()):
                    if record not in repair_records:
                        repair_records.append(record)
            if repair_records:
                canonical.boundary_repair_records = tuple(repair_records)
            validate_candidate_source_reconstruction(
                full_dna,
                canonical.list_repeats,
                canonical.list_repeat_starts,
                canonical.list_spacers,
            )

            output_index = min(index for index, _candidate in cluster)
            replacements[output_index] = canonical
            removed_indexes.update(
                index for index, _candidate in cluster
                if index != output_index
            )
            stats["groups_collapsed"] += 1
            stats["candidates_removed"] += len(cluster) - 1
            provenance.append({
                "status": "collapsed",
                "reason": "unique_source_score_winner",
                "locus": locus,
                "architecture_signature": signature,
                "background_window": background_window,
                "background_probabilities": dict(background_probabilities),
                "member_count": len(member_candidates),
                "member_candidate_indexes": tuple(
                    index for index, _candidate in cluster
                ),
                "member_candidate_ids": tuple(
                    id(candidate) for candidate in member_candidates
                ),
                "member_geometries": tuple(geometry_scores),
                "source_scores": tuple(
                    (geometry, score) for geometry, score in all_ranked
                ),
                "source_boundary_supported_scores": tuple(
                    boundary_supported_ranked
                ),
                "canonical_candidate": canonical,
                "canonical_candidate_id": id(canonical),
                "output_candidate_index": output_index,
                "canonical_geometry": winner_geometry,
                "canonical_source_score": winner_score,
                "winner_margin": winner_margin,
                "boundary_evidence": boundary_evidence,
                "classifier_used": False,
            })

        collapsed = []
        for candidate_index, candidate in enumerate(candidates):
            if candidate_index in removed_indexes:
                continue
            collapsed.append(replacements.get(candidate_index, candidate))
        collapsed_by_locus[locus] = collapsed
        stats["candidates_output"] += len(collapsed)

    return collapsed_by_locus, provenance, stats

SOURCE_TERMINAL_COMPLETION_RULE = {
    "substitution_rate": 0.20,
    "minimum_substitutions": 4,
    "maximum_substitutions": 8,
    "exact_anchor_length": 10,
    "minimum_spacer": 10,
    "maximum_spacer": 110,
    "minimum_spacer_radius": 12,
    "mad_multiplier": 3,
}


def _validate_source_terminal_rule(rule):
    required = set(SOURCE_TERMINAL_COMPLETION_RULE)
    if set(rule) != required:
        raise ValueError(
            "source-backed terminal completion rule fields are incomplete"
        )

    integer_fields = (
        "minimum_substitutions",
        "maximum_substitutions",
        "exact_anchor_length",
        "minimum_spacer",
        "maximum_spacer",
        "minimum_spacer_radius",
    )
    for field in integer_fields:
        value = rule[field]
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError("{} must be an integer".format(field))

    for field in ("substitution_rate", "mad_multiplier"):
        value = rule[field]
        if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)):
            raise ValueError("{} must be a finite number".format(field))

    if not 0 <= rule["substitution_rate"] <= 1:
        raise ValueError("substitution_rate must be between zero and one")
    if rule["minimum_substitutions"] < 0:
        raise ValueError("minimum_substitutions must be non-negative")
    if rule["maximum_substitutions"] < rule["minimum_substitutions"]:
        raise ValueError(
            "maximum_substitutions must not be smaller than the minimum"
        )
    if rule["exact_anchor_length"] < 1:
        raise ValueError("exact_anchor_length must be positive")
    if rule["minimum_spacer"] < 0:
        raise ValueError("minimum_spacer must be non-negative")
    if rule["maximum_spacer"] < rule["minimum_spacer"]:
        raise ValueError("maximum_spacer must not be smaller than the minimum")
    if rule["minimum_spacer_radius"] < 0:
        raise ValueError("minimum_spacer_radius must be non-negative")
    if rule["mad_multiplier"] < 0:
        raise ValueError("mad_multiplier must be non-negative")


def _source_terminal_spacer_bounds(crispr_candidate, rule):
    repeats = crispr_candidate.list_repeats
    starts = crispr_candidate.list_repeat_starts
    gaps = [
        right_start - (left_start + len(left_repeat))
        for left_start, left_repeat, right_start in zip(
            starts, repeats, starts[1:]
        )
    ]
    if not gaps or any(gap < 0 for gap in gaps):
        raise ValueError(
            "source-backed terminal completion requires non-overlapping repeats"
        )
    median = statistics.median(gaps)
    mad = statistics.median(abs(gap - median) for gap in gaps)
    radius = max(
        rule["minimum_spacer_radius"],
        rule["mad_multiplier"] * mad,
    )
    lower = max(rule["minimum_spacer"], math.floor(median - radius))
    upper = min(rule["maximum_spacer"], math.ceil(median + radius))
    return lower, upper, median, mad


def _source_terminal_max_substitutions(consensus_length, rule):
    raw = math.floor(rule["substitution_rate"] * consensus_length)
    return max(
        rule["minimum_substitutions"],
        min(rule["maximum_substitutions"], raw),
    )


def _has_exact_source_terminal_anchor(consensus, sequence, anchor_length):
    if anchor_length <= 0 or len(consensus) < anchor_length:
        return False
    return any(
        set(consensus[offset:offset + anchor_length]) <= set("ACGT")
        and consensus[offset:offset + anchor_length] in sequence
        for offset in range(len(consensus) - anchor_length + 1)
    )


def source_backed_terminal_repeat_proposals(
        full_dna, crispr_candidate, rule=None):
    """Return at most one conservative source-derived repeat per array side.

    The proposal rule is deliberately one-step and substitution-only. It uses
    the candidate's own consensus and source-derived spacer distribution; no
    comparator coordinates or trained-model outcomes participate.
    """

    rule = dict(SOURCE_TERMINAL_COMPLETION_RULE if rule is None else rule)
    _validate_source_terminal_rule(rule)
    validate_candidate_source_reconstruction(
        full_dna,
        crispr_candidate.list_repeats,
        crispr_candidate.list_repeat_starts,
        crispr_candidate.list_spacers,
    )
    consensus = crispr_candidate.consensus
    repeat_length = len(consensus)
    anchor_length = rule["exact_anchor_length"]
    if repeat_length < anchor_length:
        return ()
    lower, upper, median, mad = _source_terminal_spacer_bounds(
        crispr_candidate, rule
    )
    if lower > upper:
        return ()
    allowed_substitutions = _source_terminal_max_substitutions(
        repeat_length, rule
    )
    first_start = crispr_candidate.list_repeat_starts[0]
    last_end = (
        crispr_candidate.list_repeat_starts[-1]
        + len(crispr_candidate.list_repeats[-1])
    )
    proposals = []
    for side, boundary in (("left", first_start), ("right", last_end)):
        candidates = []
        for spacer_length in range(lower, upper + 1):
            if side == "left":
                repeat_start = boundary - spacer_length - repeat_length
            else:
                repeat_start = boundary + spacer_length
            repeat_end = repeat_start + repeat_length
            if repeat_start < 0 or repeat_end > len(full_dna):
                continue
            sequence = full_dna[repeat_start:repeat_end]
            substitutions = sum(
                observed != expected
                for observed, expected in zip(sequence, consensus)
            )
            if substitutions > allowed_substitutions:
                continue
            if not _has_exact_source_terminal_anchor(
                    consensus, sequence, anchor_length):
                continue
            candidates.append({
                "side": side,
                "repeat_interval": (repeat_start, repeat_end),
                "spacer_length": spacer_length,
                "substitutions": substitutions,
                "spacer_median": median,
                "spacer_mad": mad,
                "maximum_substitutions": allowed_substitutions,
            })
        if candidates:
            proposals.append(min(
                candidates,
                key=lambda item: (
                    item["substitutions"],
                    abs(item["spacer_length"] - median),
                    item["spacer_length"],
                    item["repeat_interval"][0],
                ),
            ))
    return tuple(proposals)


SOURCE_WEAK_TERMINAL_COMPLETION_RULE = {
    "minimum_repeats": 5,
    "minimum_repeat_length": 20,
    "maximum_repeat_length": 55,
    "minimum_mean_repeat_length": 20,
    "maximum_mean_repeat_length": 55,
    "minimum_core_repeat_identity": 0.85,
    "minimum_core_mean_repeat_identity": 0.90,
    "minimum_terminal_repeat_identity": 0.70,
    "maximum_terminal_repeat_identity": 0.85,
    "minimum_final_mean_repeat_identity": 0.90,
    "minimum_exact_anchor_length": 6,
    "minimum_spacer_length": 18,
    "maximum_spacer_length": 80,
    "maximum_spacer_length_cv": 0.05,
    "maximum_distinct_spacer_identity": 0.80,
    "terminal_kmer_length": 5,
    "maximum_terminal_kmer_dominance": 0.80,
    "minimum_spacer_radius": 2,
    "mad_multiplier": 3.0,
    "minimum_column_coverage_fraction": 0.80,
    "spacer_edge_width": 3,
    "background_flank": 5000,
    "minimum_phase_winner_margin": 10.0,
}

_STANDARD_GENETIC_CODE = dict(zip(
    (
        first + second + third
        for first in "TCAG"
        for second in "TCAG"
        for third in "TCAG"
    ),
    (
        "FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG"
    ),
))
_DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _validate_source_weak_terminal_rule(rule):
    if set(rule) != set(SOURCE_WEAK_TERMINAL_COMPLETION_RULE):
        raise ValueError(
            "source-backed weak-terminal completion rule fields are incomplete"
        )

    integer_fields = (
        "minimum_repeats",
        "minimum_repeat_length",
        "maximum_repeat_length",
        "minimum_mean_repeat_length",
        "maximum_mean_repeat_length",
        "minimum_exact_anchor_length",
        "minimum_spacer_length",
        "maximum_spacer_length",
        "terminal_kmer_length",
        "minimum_spacer_radius",
        "spacer_edge_width",
        "background_flank",
    )
    for field in integer_fields:
        value = rule[field]
        if (
                not isinstance(value, int)
                or isinstance(value, bool)
                or value < 0):
            raise ValueError("{} must be a non-negative integer".format(field))

    if rule["minimum_repeats"] < 5:
        raise ValueError("minimum_repeats must be at least five")
    if rule["minimum_repeat_length"] < 1:
        raise ValueError("minimum_repeat_length must be positive")
    if rule["maximum_repeat_length"] < rule["minimum_repeat_length"]:
        raise ValueError("repeat-length bounds are inconsistent")
    if (
            rule["maximum_mean_repeat_length"]
            < rule["minimum_mean_repeat_length"]):
        raise ValueError("mean repeat-length bounds are inconsistent")
    if rule["maximum_spacer_length"] < rule["minimum_spacer_length"]:
        raise ValueError("spacer-length bounds are inconsistent")
    if rule["minimum_exact_anchor_length"] < 1:
        raise ValueError("minimum_exact_anchor_length must be positive")
    if rule["terminal_kmer_length"] < 1:
        raise ValueError("terminal_kmer_length must be positive")
    if rule["spacer_edge_width"] < 1:
        raise ValueError("spacer_edge_width must be positive")

    probability_fields = (
        "minimum_core_repeat_identity",
        "minimum_core_mean_repeat_identity",
        "minimum_terminal_repeat_identity",
        "maximum_terminal_repeat_identity",
        "minimum_final_mean_repeat_identity",
        "maximum_distinct_spacer_identity",
        "maximum_terminal_kmer_dominance",
        "minimum_column_coverage_fraction",
    )
    for field in probability_fields:
        value = rule[field]
        if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or not 0.0 <= value <= 1.0):
            raise ValueError("{} must be finite and in [0, 1]".format(field))
    if not 0.0 < rule["minimum_column_coverage_fraction"] <= 1.0:
        raise ValueError(
            "minimum_column_coverage_fraction must be in (0, 1]"
        )
    if (
            rule["minimum_terminal_repeat_identity"]
            >= rule["maximum_terminal_repeat_identity"]
            or rule["maximum_terminal_repeat_identity"]
            > rule["minimum_core_repeat_identity"]):
        raise ValueError("weak and core repeat-identity bounds are inconsistent")

    non_negative_fields = (
        "maximum_spacer_length_cv",
        "mad_multiplier",
        "minimum_phase_winner_margin",
    )
    for field in non_negative_fields:
        value = rule[field]
        if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
                or value < 0.0):
            raise ValueError("{} must be finite and non-negative".format(field))


def _source_edit_distance(first, second):
    return sum(
        len(changes) for changes in _minimal_edit_changes(first, second)
    )


def _source_edit_identity(first, second):
    denominator = max(len(first), len(second))
    if denominator == 0:
        return 1.0
    return 1.0 - _source_edit_distance(first, second) / float(denominator)


def _longest_exact_source_anchor(first, second):
    previous = [0] * (len(second) + 1)
    maximum = 0
    for first_base in first:
        current = [0]
        for column, second_base in enumerate(second, 1):
            length = (
                previous[column - 1] + 1
                if first_base == second_base else 0
            )
            current.append(length)
            maximum = max(maximum, length)
        previous = current
    return maximum


def _source_spacer_statistics(spacers, terminal_kmer_length):
    lengths = tuple(map(len, spacers))
    mean = sum(lengths) / float(len(lengths)) if lengths else 0.0
    variance = (
        sum((length - mean) ** 2 for length in lengths) / len(lengths)
        if lengths else 0.0
    )
    maximum_identity = 0.0
    for index, first in enumerate(spacers[:-1]):
        for second in spacers[index + 1:]:
            maximum_identity = max(
                maximum_identity,
                _source_edit_identity(first, second),
            )
    terminal_kmers = Counter(
        spacer[-terminal_kmer_length:]
        for spacer in spacers
        if len(spacer) >= terminal_kmer_length
    )
    terminal_dominance = (
        terminal_kmers.most_common(1)[0][1] / float(len(spacers))
        if spacers else 0.0
    )
    return {
        "lengths": lengths,
        "cv": math.sqrt(variance) / mean if mean else 0.0,
        "distinct_count": len(set(spacers)),
        "maximum_identity": maximum_identity,
        "terminal_kmer_dominance": terminal_dominance,
    }


def _source_weak_structure_metrics(full_dna, repeat_intervals, rule):
    repeat_intervals = tuple(sorted(repeat_intervals))
    repeats = tuple(
        full_dna[start:end] for start, end in repeat_intervals
    )
    consensus = _source_only_alignment_seed(repeats)
    repeat_identities = tuple(
        _source_edit_identity(repeat, consensus) for repeat in repeats
    )
    spacers = tuple(
        full_dna[left_end:right_start]
        for (_, left_end), (right_start, _) in zip(
            repeat_intervals, repeat_intervals[1:]
        )
    )
    spacer_stats = _source_spacer_statistics(
        spacers, rule["terminal_kmer_length"]
    )
    return {
        "consensus": consensus,
        "repeat_identities": repeat_identities,
        "minimum_repeat_identity": min(repeat_identities),
        "mean_repeat_identity": (
            sum(repeat_identities) / float(len(repeat_identities))
        ),
        "spacers": spacers,
        "spacer_stats": spacer_stats,
    }


def _translate_source_dna(sequence, frame):
    return "".join(
        _STANDARD_GENETIC_CODE.get(sequence[offset:offset + 3], "X")
        for offset in range(frame, len(sequence) - 2, 3)
    )


def source_translated_period_tandem_identity(full_dna, repeat_intervals):
    """Return the strongest six-frame peptide-unit identity, or ``None``.

    Periodic CRISPR-like DNA is protein-compatible only when every repeat-start
    delta preserves codon phase. In that case, source units from one repeat
    start to the next are translated in all frames and both directions. The
    maximum mean adjacent peptide edit identity is a conservative tandem-protein
    warning, not a CRISPR classifier.
    """

    repeat_intervals = tuple(sorted(repeat_intervals))
    starts = tuple(start for start, _end in repeat_intervals)
    if (
            len(starts) < 3
            or any((right - left) % 3 for left, right in zip(
                starts, starts[1:]
            ))):
        return None
    units = tuple(
        full_dna[left:right] for left, right in zip(starts, starts[1:])
    )
    orientations = (
        units,
        tuple(
            unit.translate(_DNA_COMPLEMENT)[::-1]
            for unit in reversed(units)
        ),
    )
    maximum = 0.0
    for oriented_units in orientations:
        for frame in range(3):
            peptides = tuple(
                _translate_source_dna(unit, frame)
                for unit in oriented_units
            )
            identities = tuple(
                _source_edit_identity(first, second)
                for first, second in zip(peptides, peptides[1:])
                if first and second
            )
            if identities:
                maximum = max(
                    maximum,
                    sum(identities) / float(len(identities)),
                )
    return maximum


def _source_weak_terminal_core_is_eligible(metrics, repeat_count, rule):
    spacer_stats = metrics["spacer_stats"]
    repeat_length = len(metrics["consensus"])
    return (
        repeat_count >= rule["minimum_repeats"]
        and rule["minimum_repeat_length"]
        <= repeat_length
        <= rule["maximum_repeat_length"]
        and metrics["minimum_repeat_identity"]
        >= rule["minimum_core_repeat_identity"]
        and metrics["mean_repeat_identity"]
        >= rule["minimum_core_mean_repeat_identity"]
        and spacer_stats["lengths"]
        and all(
            rule["minimum_spacer_length"]
            <= length
            <= rule["maximum_spacer_length"]
            for length in spacer_stats["lengths"]
        )
        and spacer_stats["cv"] <= rule["maximum_spacer_length_cv"]
        and spacer_stats["distinct_count"] == len(spacer_stats["lengths"])
        and spacer_stats["maximum_identity"]
        <= rule["maximum_distinct_spacer_identity"]
        and spacer_stats["terminal_kmer_dominance"]
        < rule["maximum_terminal_kmer_dominance"]
    )


def _source_weak_phase_background(full_dna, geometries, flank):
    starts = tuple(
        start for geometry in geometries for start, _end in geometry
    )
    ends = tuple(
        end for geometry in geometries for _start, end in geometry
    )
    background_start = max(0, min(starts) - flank)
    background_end = min(len(full_dna), max(ends) + flank)
    counts = _canonical_base_counts(full_dna[background_start:background_end])
    canonical_count = sum(counts.values())
    return {
        base: (counts[base] + 1.0) / (canonical_count + 4.0)
        for base in _CANONICAL_DNA_BASES
    }


def source_backed_weak_terminal_repeat_proposals(
        full_dna, crispr_candidate, rule=None):
    """Return at most one source-only, strongly gated weak terminal repeat.

    This tier exists for a single divergent terminal that the ordinary exact
    ten-base-anchor rule cannot recover. It never changes global phase
    thresholds and never uses comparator coordinates, classifier scores, Cas
    context, or annotations. Every cadence-window position competes only in
    one frozen source-phase objective; only a unique winner with a large margin
    may proceed.

    The returned proposal explicitly requires assessed full-source coding
    context before output promotion. The architecture integration enforces
    that fail-closed contract.
    """

    weak_rule = dict(
        SOURCE_WEAK_TERMINAL_COMPLETION_RULE if rule is None else rule
    )
    _validate_source_weak_terminal_rule(weak_rule)
    validate_candidate_source_reconstruction(
        full_dna,
        crispr_candidate.list_repeats,
        crispr_candidate.list_repeat_starts,
        crispr_candidate.list_spacers,
    )
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(
            crispr_candidate.list_repeat_starts,
            crispr_candidate.list_repeats,
        )
    )
    core_metrics = _source_weak_structure_metrics(
        full_dna, repeat_intervals, weak_rule
    )
    if not _source_weak_terminal_core_is_eligible(
            core_metrics, len(repeat_intervals), weak_rule):
        return ()

    gaps = tuple(
        right_start - left_end
        for (_left_start, left_end), (right_start, _right_end)
        in zip(repeat_intervals, repeat_intervals[1:])
    )
    median = statistics.median(gaps)
    mad = statistics.median(abs(gap - median) for gap in gaps)
    radius = max(
        weak_rule["minimum_spacer_radius"],
        math.ceil(weak_rule["mad_multiplier"] * mad),
    )
    lower = max(
        weak_rule["minimum_spacer_length"], math.floor(median - radius)
    )
    upper = min(
        weak_rule["maximum_spacer_length"], math.ceil(median + radius)
    )
    if lower > upper:
        return ()

    repeat_length = len(core_metrics["consensus"])
    first_start = repeat_intervals[0][0]
    last_end = repeat_intervals[-1][1]
    hypotheses = []
    for side, boundary in (("left", first_start), ("right", last_end)):
        for spacer_length in range(lower, upper + 1):
            repeat_start = (
                boundary - spacer_length - repeat_length
                if side == "left" else boundary + spacer_length
            )
            repeat_end = repeat_start + repeat_length
            if repeat_start < 0 or repeat_end > len(full_dna):
                continue
            interval = (repeat_start, repeat_end)
            geometry = tuple(sorted(repeat_intervals + (interval,)))
            hypotheses.append({
                "side": side,
                "repeat_interval": interval,
                "spacer_length": spacer_length,
                "geometry": geometry,
            })
    if len(hypotheses) < 2:
        return ()

    background = _source_weak_phase_background(
        full_dna,
        tuple(hypothesis["geometry"] for hypothesis in hypotheses),
        weak_rule["background_flank"],
    )
    scored = []
    for hypothesis in hypotheses:
        geometry = hypothesis["geometry"]
        score = _source_phase_geometry_score(
            full_dna,
            tuple(start for start, _end in geometry),
            tuple(end for _start, end in geometry),
            weak_rule,
            background,
        )
        if score is not None:
            scored.append((score, hypothesis))
    if len(scored) < 2:
        return ()
    scored.sort(key=lambda item: (
        item[0],
        item[1]["side"] == "right",
        tuple(-value for value in item[1]["repeat_interval"]),
    ), reverse=True)
    winner_score, winner = scored[0]
    winner_margin = winner_score - scored[1][0]
    if winner_margin < weak_rule["minimum_phase_winner_margin"]:
        return ()

    repeat_start, repeat_end = winner["repeat_interval"]
    repeat = full_dna[repeat_start:repeat_end]
    repeat_identity = _source_edit_identity(
        repeat, core_metrics["consensus"]
    )
    if not (
            weak_rule["minimum_terminal_repeat_identity"]
            <= repeat_identity
            < weak_rule["maximum_terminal_repeat_identity"]):
        return ()
    exact_anchor_length = _longest_exact_source_anchor(
        repeat, core_metrics["consensus"]
    )
    if exact_anchor_length < weak_rule["minimum_exact_anchor_length"]:
        return ()

    final_metrics = _source_weak_structure_metrics(
        full_dna, winner["geometry"], weak_rule
    )
    final_spacers = final_metrics["spacer_stats"]
    if (
            final_metrics["mean_repeat_identity"]
            < weak_rule["minimum_final_mean_repeat_identity"]
            or final_metrics["minimum_repeat_identity"]
            < weak_rule["minimum_terminal_repeat_identity"]
            or final_spacers["cv"] > weak_rule["maximum_spacer_length_cv"]
            or final_spacers["distinct_count"] != len(
                final_spacers["lengths"]
            )
            or final_spacers["maximum_identity"]
            > weak_rule["maximum_distinct_spacer_identity"]
            or final_spacers["terminal_kmer_dominance"]
            >= weak_rule["maximum_terminal_kmer_dominance"]):
        return ()

    # Peptide-period similarity is evidence, not an independently calibrated
    # CRISPR veto. Carry it into the lineage record so downstream review can
    # combine it with the already-batched full-source CDS/HMM context.
    translated_tandem_identity = source_translated_period_tandem_identity(
        full_dna, winner["geometry"]
    )
    return ({
        "tier": "weak",
        "side": winner["side"],
        "repeat_interval": winner["repeat_interval"],
        "spacer_length": winner["spacer_length"],
        "spacer_median": median,
        "spacer_mad": mad,
        "repeat_identity": repeat_identity,
        "exact_anchor_length": exact_anchor_length,
        "source_phase_score": winner_score,
        "source_phase_winner_margin": winner_margin,
        "translated_tandem_identity": translated_tandem_identity,
        "requires_source_coding_context": True,
    },)


def apply_source_backed_weak_terminal_repeat_proposal(
        full_dna, crispr_candidate, proposal, rule=None):
    """Build one independently authorized weak-terminal candidate."""

    authorized = source_backed_weak_terminal_repeat_proposals(
        full_dna, crispr_candidate, rule=rule
    )
    if len(authorized) != 1:
        raise ValueError("the candidate has no authorized weak terminal")
    expected = authorized[0]
    if (
            not isinstance(proposal, dict)
            or proposal.get("side") != expected["side"]
            or tuple(proposal.get("repeat_interval", ()))
            != tuple(expected["repeat_interval"])):
        raise ValueError(
            "weak terminal proposal was not generated by the source-backed rule"
        )
    interval = tuple(expected["repeat_interval"])
    repeat_intervals = [
        (start, start + len(repeat))
        for start, repeat in zip(
            crispr_candidate.list_repeat_starts,
            crispr_candidate.list_repeats,
        )
    ]
    repeat_intervals.append(interval)
    repeat_intervals.sort()
    completed = build_source_backed_candidate_from_boundaries(
        full_dna,
        crispr_candidate,
        tuple(start for start, _end in repeat_intervals),
        tuple(end for _start, end in repeat_intervals),
    )
    validate_candidate_source_reconstruction(
        full_dna,
        completed.list_repeats,
        completed.list_repeat_starts,
        completed.list_spacers,
    )
    return completed


def apply_source_backed_terminal_repeat_proposals(
        full_dna, crispr_candidate, proposals, rule=None):
    """Build and independently reconstruct a candidate plus terminal repeats."""

    proposals = tuple(proposals)
    if not proposals:
        return crispr_candidate
    authorized = {
        (proposal["side"], tuple(proposal["repeat_interval"]))
        for proposal in source_backed_terminal_repeat_proposals(
            full_dna, crispr_candidate, rule=rule
        )
    }
    by_side = {}
    for proposal in proposals:
        side = proposal.get("side")
        interval = proposal.get("repeat_interval")
        if side not in ("left", "right") or side in by_side:
            raise ValueError("terminal proposals must have unique valid sides")
        if (
                not isinstance(interval, (tuple, list))
                or len(interval) != 2
                or not all(isinstance(value, int) and not isinstance(value, bool)
                           for value in interval)
                or interval[0] < 0
                or interval[1] <= interval[0]
                or interval[1] > len(full_dna)):
            raise ValueError("terminal proposal interval is invalid")
        if interval[1] - interval[0] != len(crispr_candidate.consensus):
            raise ValueError(
                "terminal proposal length must equal the input consensus length"
            )
        if (side, tuple(interval)) not in authorized:
            raise ValueError(
                "terminal proposal was not generated by the source-backed rule"
            )
        by_side[side] = (interval[0], interval[1])

    starts = list(crispr_candidate.list_repeat_starts)
    repeats = list(crispr_candidate.list_repeats)
    first_start = starts[0]
    last_end = starts[-1] + len(repeats[-1])
    if "left" in by_side:
        start, end = by_side["left"]
        if end > first_start:
            raise ValueError("left terminal proposal overlaps the input array")
        starts.insert(0, start)
        repeats.insert(0, full_dna[start:end])
    if "right" in by_side:
        start, end = by_side["right"]
        if start < last_end:
            raise ValueError("right terminal proposal overlaps the input array")
        starts.append(start)
        repeats.append(full_dna[start:end])

    spacers = [
        full_dna[start + len(repeat):next_start]
        for start, repeat, next_start in zip(starts, repeats, starts[1:])
    ]
    consensus = crispr_candidate.consensus
    error_indexes = [
        _minimal_edit_changes(consensus, repeat) for repeat in repeats
    ]
    fuzzy_counts = [
        tuple(len(error_positions) for error_positions in changes)
        for changes in error_indexes
    ]
    representation = DotRepresentationMaker(
        list_repeats=repeats,
        list_repeats_starts=starts,
        list_spacers=spacers,
        consensus=consensus,
        relative_error_indexes=error_indexes,
        fuzzy_counts=fuzzy_counts,
    )
    completed = CrisprCandidate(
        list_repeats=repeats,
        list_repeats_gaped=representation.list_gaped_repeats,
        list_spacers=spacers,
        list_repeat_starts=starts,
    )
    inherit_candidate_origins(completed, (crispr_candidate,))
    validate_candidate_source_reconstruction(
        full_dna,
        completed.list_repeats,
        completed.list_repeat_starts,
        completed.list_spacers,
    )
    return completed


class IterativeDegeneratedSearch:
    def __init__(self, full_dna, repeat_seq_candidate,
                 spacer_margin, repeat_seq_candidate_gaped,
                 list_repeats_starts, list_repeats, list_spacers,
                 start_flanking_region_left, end_flanking_region_right,
                 allowed_max_editing_distance,
                 iterative_size_flanking_region, prevent_long_spacers=True,
                 attempt_to_improve_initial_array=True):

        self.final_coordinates = None
        self.full_dna = full_dna
        self.repeat_seq_candidate = repeat_seq_candidate
        self.spacer_margin = spacer_margin
        self.repeats_seq_candidate_gaped = repeat_seq_candidate_gaped
        self.list_repeats_starts = list_repeats_starts
        self.list_repeats = list_repeats
        self.list_spacers = list_spacers
        self.start_flanking_region_left = start_flanking_region_left
        self.end_flanking_region_right = end_flanking_region_right
        self.iterative_size_flanking_region = iterative_size_flanking_region
        self.allowed_max_editing_distance = allowed_max_editing_distance
        self.editing_distances = self._build_editing_distance_schedule(
            allowed_max_editing_distance
        )
        self.prevent_long_spacers = prevent_long_spacers
        self.attempt_to_improve_initial_array = attempt_to_improve_initial_array

        self.flag_make_step_left = True
        self.flag_make_step_right = True

        validate_candidate_source_reconstruction(
            self.full_dna,
            self.list_repeats,
            self.list_repeats_starts,
            self.list_spacers,
        )
        initial_array_start = self.list_repeats_starts[0]
        initial_array_end = (
            self.list_repeats_starts[-1] + len(self.list_repeats[-1])
        )
        if (
                not isinstance(self.start_flanking_region_left, int)
                or isinstance(self.start_flanking_region_left, bool)
                or not isinstance(self.end_flanking_region_right, int)
                or isinstance(self.end_flanking_region_right, bool)
                or self.start_flanking_region_left < 0
                or self.start_flanking_region_left > initial_array_start
                or self.end_flanking_region_right < initial_array_end
                or self.end_flanking_region_right > len(self.full_dna)):
            raise ValueError(
                "enhancement boundaries must contain the source-backed array"
            )
        self._calculate_representation_original_array()
        self._calculate_max_possible_new_spacer_length()

        self._left_flank_iterative_search()
        self._right_flank_iterative_search()
        self._build_new_representation()

    def _calculate_representation_original_array(self):
        self.list_repeats_original = list(self.list_repeats)
        self.list_repeats_starts_original = list(self.list_repeats_starts)
        self.list_spacers_original = list(self.list_spacers)
        self.list_repeats_error_indexes_original = [
            _minimal_edit_changes(self.repeat_seq_candidate, repeat)
            for repeat in self.list_repeats_original
        ]
        self.list_fuzzy_counts_original = [
            tuple(len(error_positions) for error_positions in changes)
            for changes in self.list_repeats_error_indexes_original
        ]

    def _calculate_max_possible_new_spacer_length(self):
        if self.list_spacers:
            self.max_possible_new_spacer_len = max((len(x) for x in self.list_spacers)) + self.spacer_margin
        else:
            self.max_possible_new_spacer_len = 60

    @staticmethod
    def _build_editing_distance_schedule(allowed_max_editing_distance):
        if (not isinstance(allowed_max_editing_distance, int)
                or isinstance(allowed_max_editing_distance, bool)
                or allowed_max_editing_distance < 0):
            raise ValueError("maximum edit distance must be a nonnegative integer")

        first_distance = min(6, allowed_max_editing_distance)
        distances = list(range(first_distance, allowed_max_editing_distance + 1, 2))
        if not distances or distances[-1] != allowed_max_editing_distance:
            distances.append(allowed_max_editing_distance)
        return tuple(distances)

    def _make_new_step_left(self, interval_start):
        if interval_start is not None:
            interval_end = interval_start
            interval_start = max(self.start_flanking_region_left,
                                 interval_start - self.iterative_size_flanking_region)

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = self.editing_distances[0]
            return interval_start, interval_end, interval_dna, distance
        else:
            interval_start = max(self.start_flanking_region_left,
                                 self.list_repeats_starts[0] - self.iterative_size_flanking_region)
            interval_end = self.list_repeats_starts[0]

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = self.editing_distances[0]
            return interval_start, interval_end, interval_dna, distance

    def _make_new_step_right(self, interval_end):
        if interval_end is not None:
            interval_start = interval_end
            interval_end = min(self.end_flanking_region_right,
                               interval_end + self.iterative_size_flanking_region)

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = self.editing_distances[0]
            return interval_start, interval_end, interval_dna, distance
        else:
            interval_start = self.list_repeats_starts[-1] + len(self.list_repeats[-1])
            interval_end = min(self.end_flanking_region_right,
                               interval_start + self.iterative_size_flanking_region)

            interval_dna = self.full_dna[interval_start:interval_end]
            distance = self.editing_distances[0]
            return interval_start, interval_end, interval_dna, distance

    def increase_editing_distance(self, distance):
        for scheduled_distance in self.editing_distances:
            if scheduled_distance > distance:
                return scheduled_distance
        return self.editing_distances[-1]

    def _create_pattern(self, distance):
        if distance not in self.editing_distances:
            raise ValueError("edit distance is outside the configured schedule")
        half_distance = (distance + 1) // 2
        pattern = f"(?be)({self.repeat_seq_candidate})" + "{" + \
                  f"s<={half_distance},i<={half_distance},d<={half_distance},s+i+d<={distance}" + "}"
        return pattern

    def _left_flank_iterative_search(self):
        self.list_repeats_in_intervals_left = []
        self.list_repeats_starts_in_intervals_left = []
        self.list_repeats_errors_in_intervals_left = []
        self.list_fuzzy_counts_in_intervals_left = []
        self.list_spacers_in_intervals_left = []
        self.list_bridge_spacers_left = []

        flag_stop_iterations = False
        most_left_repeat_match_start = None

        while not flag_stop_iterations:
            if self.flag_make_step_left:
                new_step_left = self._make_new_step_left(most_left_repeat_match_start)
                interval_start, interval_end, interval_dna, distance = new_step_left
            else:
                distance = self.increase_editing_distance(distance)

            search_pattern = self._create_pattern(distance)
            ibfs = IntervalBasedFuzzySearch(interval_start, interval_end, interval_dna, search_pattern)

            if ibfs.flag_exist_match:
                bridge_spacer = ibfs.compute_left_bridge_spacer()
                len_bridge_spacer_to_the_left = len(bridge_spacer)
                if len_bridge_spacer_to_the_left < self.max_possible_new_spacer_len:
                    list_repeats = ibfs.get_repeats()
                    list_spacers = ibfs.get_internal_spacers()
                    list_repeats_errors = ibfs.get_relative_error_indexes()

                    list_repeats_relative_starts = ibfs.get_relative_starts()
                    list_repeats_starts = [(x + interval_start) for x in list_repeats_relative_starts]
                    most_left_repeat_match_start = list_repeats_starts[0]
                    list_fuzzy_counts = ibfs.get_errors()

                    self.list_repeats_in_intervals_left.append(list_repeats)
                    self.list_repeats_starts_in_intervals_left.append(list_repeats_starts)
                    self.list_repeats_errors_in_intervals_left.append(list_repeats_errors)
                    self.list_fuzzy_counts_in_intervals_left.append(list_fuzzy_counts)
                    self.list_spacers_in_intervals_left.append(list_spacers)
                    self.list_bridge_spacers_left.append(bridge_spacer)

                    self.flag_make_step_left = True
                else:
                    self.flag_make_step_left = False
                    if distance == self.allowed_max_editing_distance:
                        flag_stop_iterations = True
            else:
                self.flag_make_step_left = False
                if distance == self.allowed_max_editing_distance:
                    flag_stop_iterations = True

    def _right_flank_iterative_search(self):
        self.list_repeats_in_intervals_right = []
        self.list_repeats_starts_in_intervals_right = []
        self.list_repeats_errors_in_intervals_right = []
        self.list_fuzzy_counts_in_intervals_right = []
        self.list_spacers_in_intervals_right = []
        self.list_bridge_spacers_right = []

        flag_stop_iterations = False
        most_right_repeat_match_end = None

        while not flag_stop_iterations:
            if self.flag_make_step_right:
                new_step_right = self._make_new_step_right(most_right_repeat_match_end)
                interval_start, interval_end, interval_dna, distance = new_step_right
            else:
                distance = self.increase_editing_distance(distance)

            search_pattern = self._create_pattern(distance)
            ibfs = IntervalBasedFuzzySearch(interval_start, interval_end, interval_dna, search_pattern)

            if ibfs.flag_exist_match:
                bridge_spacer = ibfs.compute_right_bridge_spacer()
                len_bridge_spacer_to_the_right = len(bridge_spacer)
                if len_bridge_spacer_to_the_right < self.max_possible_new_spacer_len:
                    list_repeats = ibfs.get_repeats()
                    list_spacers = ibfs.get_internal_spacers()
                    list_repeats_errors = ibfs.get_relative_error_indexes()

                    list_repeats_relative_starts = ibfs.get_relative_starts()
                    list_repeats_starts = [(x + interval_start) for x in list_repeats_relative_starts]
                    most_right_repeat_match_end = list_repeats_starts[-1] + len(list_repeats[-1])
                    list_fuzzy_counts = ibfs.get_errors()

                    self.list_repeats_in_intervals_right.append(list_repeats)
                    self.list_repeats_starts_in_intervals_right.append(list_repeats_starts)
                    self.list_repeats_errors_in_intervals_right.append(list_repeats_errors)
                    self.list_fuzzy_counts_in_intervals_right.append(list_fuzzy_counts)
                    self.list_spacers_in_intervals_right.append(list_spacers)
                    self.list_bridge_spacers_right.append(bridge_spacer)

                    self.flag_make_step_right = True
                else:
                    self.flag_make_step_right = False
                    if distance == self.allowed_max_editing_distance:
                        flag_stop_iterations = True
            else:
                self.flag_make_step_right = False
                if distance == self.allowed_max_editing_distance:
                    flag_stop_iterations = True

    def _build_new_representation(self):
        self.list_repeats_in_intervals_left = self.list_repeats_in_intervals_left[::-1]
        self.list_repeats_starts_in_intervals_left = self.list_repeats_starts_in_intervals_left[::-1]
        self.list_repeats_errors_in_intervals_left = self.list_repeats_errors_in_intervals_left[::-1]
        self.list_fuzzy_counts_in_intervals_left = self.list_fuzzy_counts_in_intervals_left[::-1]
        self.list_spacers_in_intervals_left = self.list_spacers_in_intervals_left[::-1]
        self.list_bridge_spacers_left = self.list_bridge_spacers_left[::-1]

        list_repeats_left = [r for list_rep in self.list_repeats_in_intervals_left for r in list_rep]
        list_repeats_starts_left = [s for list_starts in self.list_repeats_starts_in_intervals_left for s in
                                    list_starts]
        list_repeats_errors_in_interval_left = [l_err for list_all_errors in self.list_repeats_errors_in_intervals_left
                                                for l_err in list_all_errors]

        list_fuzzy_counts_left = [fuzzy_count for list_fc in self.list_fuzzy_counts_in_intervals_left
                                  for fuzzy_count in list_fc]

        list_spacers_left = []
        for internal_spacers, bridge_spacer in zip(self.list_spacers_in_intervals_left, self.list_bridge_spacers_left):
            list_spacers_left += internal_spacers
            list_spacers_left.append(bridge_spacer)

        list_repeats_right = [r for list_rep in self.list_repeats_in_intervals_right for r in list_rep]
        list_repeats_starts_right = [s for list_starts in self.list_repeats_starts_in_intervals_right for s in
                                     list_starts]
        list_repeats_errors_in_interval_right = [l_err for list_all_errors in
                                                 self.list_repeats_errors_in_intervals_right
                                                 for l_err in list_all_errors]

        list_fuzzy_counts_right = [fuzzy_count for list_fc in self.list_fuzzy_counts_in_intervals_right
                                   for fuzzy_count in list_fc]

        list_spacers_right = []
        for internal_spacers, bridge_spacer in zip(self.list_spacers_in_intervals_right,
                                                   self.list_bridge_spacers_right):
            list_spacers_right.append(bridge_spacer)
            list_spacers_right += internal_spacers

        list_repeats = list_repeats_left + self.list_repeats_original + list_repeats_right
        list_repeats_starts = list_repeats_starts_left + self.list_repeats_starts_original + list_repeats_starts_right
        list_spacers = list_spacers_left + self.list_spacers_original + list_spacers_right

        list_repeats_error_indexes = list_repeats_errors_in_interval_left + self.list_repeats_error_indexes_original + \
                                     list_repeats_errors_in_interval_right

        list_fuzzy_counts = list_fuzzy_counts_left + self.list_fuzzy_counts_original + list_fuzzy_counts_right

        self.final_coordinates = list_repeats_starts[0], list_repeats_starts[-1] + len(list_repeats[-1])

        dot_representation_maker = DotRepresentationMaker(list_repeats=list_repeats,
                                                          list_repeats_starts=list_repeats_starts,
                                                          list_spacers=list_spacers,
                                                          consensus=self.repeat_seq_candidate,
                                                          relative_error_indexes=list_repeats_error_indexes,
                                                          fuzzy_counts=list_fuzzy_counts)

        list_repeats_gaped = dot_representation_maker.list_gaped_repeats

        self.crispr_candidate = CrisprCandidate(list_repeats, list_repeats_gaped,
                                                list_spacers, list_repeats_starts)
        validate_candidate_source_reconstruction(
            self.full_dna,
            self.crispr_candidate.list_repeats,
            self.crispr_candidate.list_repeat_starts,
            self.crispr_candidate.list_spacers,
        )

    def output(self):
        return self.crispr_candidate


class IntervalBasedFuzzySearch:
    def __init__(self, interval_start, interval_end, interval_dna, search_pattern):
        if (
                not isinstance(interval_start, int)
                or isinstance(interval_start, bool)
                or not isinstance(interval_end, int)
                or isinstance(interval_end, bool)
                or interval_start < 0
                or interval_end < interval_start):
            raise ValueError("fuzzy-search interval is invalid")
        if (
                not isinstance(interval_dna, str)
                or len(interval_dna) != interval_end - interval_start):
            raise ValueError(
                "fuzzy-search interval DNA does not reconstruct its coordinates"
            )
        self.interval_start = interval_start
        self.interval_end = interval_end
        self.interval_dna = interval_dna
        self.search_pattern = search_pattern

        self.len_bridge_spacer_left = None
        self.len_bridge_spacer_right = None

        self._find_match()

    def _find_match(self):
        self.matches = list(regex.finditer(self.search_pattern, self.interval_dna))
        previous_end = 0
        for match_index, match in enumerate(self.matches):
            start = match.start()
            end = match.end()
            if (
                    start < previous_end
                    or end <= start
                    or end > len(self.interval_dna)):
                raise ValueError(
                    "enhancement fuzzy match has invalid source coordinates"
                )
            if self.interval_dna[start:end] != match.group():
                raise ValueError(
                    "enhancement fuzzy match does not reconstruct its source slice"
                )
            previous_end = end
        if self.matches:
            self.flag_exist_match = True
        else:
            self.flag_exist_match = False

    def compute_left_bridge_spacer(self):
        bridge_spacer_start = self.matches[-1].end()
        bridge_spacer = self.interval_dna[bridge_spacer_start:]
        return bridge_spacer

    def compute_right_bridge_spacer(self):
        bridge_spacer_end = self.matches[0].start()
        bridge_spacer = self.interval_dna[:bridge_spacer_end]
        return bridge_spacer

    def get_repeats(self):
        list_repeats = [match.group() for match in self.matches]
        return list_repeats

    def get_internal_spacers(self):
        if len(self.matches) > 1:
            internal_spacer_starts = [match.end() for match in self.matches][:-1]
            internal_spacer_end = [match.start() for match in self.matches][1:]
            internal_spacer_coordinates = [(start, end) for start, end in
                                           zip(internal_spacer_starts, internal_spacer_end)]
            list_internal_spacers = [self.interval_dna[start_end[0]:start_end[1]]
                                     for start_end in internal_spacer_coordinates]
        else:
            list_internal_spacers = []

        return list_internal_spacers

    def get_errors(self):
        return [match.fuzzy_counts for match in self.matches]

    def get_relative_starts(self):
        return [match.start() for match in self.matches]

    def get_relative_error_indexes(self):
        return [normalize_fuzzy_changes(match) for match in self.matches]






class DotRepresentationMaker:
    def __init__(self, list_repeats, list_repeats_starts, list_spacers, consensus,
                 relative_error_indexes, fuzzy_counts):
        self.list_repeats = list_repeats
        self.list_spacers = list_spacers
        self.list_relative_error_indexes = relative_error_indexes
        self.consensus = consensus
        self.list_repeat_absolute_starts = list_repeats_starts
        self.fuzzy_counts = fuzzy_counts

        self.list_all_insertions = []
        self.list_all_deletions = []

        self._create_representation()
        self._get_transitional_repeat()

    def _create_representation(self):
        self.gaped_matching_pattern, self.list_gaped_repeats = (
            merge_fuzzy_alignments(
                self.consensus,
                self.list_repeats,
                self.list_relative_error_indexes,
            )
        )
        self.list_dotted_gaped_repeats = [self.apply_dots(x, self.gaped_matching_pattern)
                                          for x in self.list_gaped_repeats]

    def _get_transitional_repeat(self):
        dict_transitions = {frozenset(("A", "T")): "W",
                            frozenset(("C", "G")): "S",
                            frozenset(("A", "C")): "M",
                            frozenset(("G", "T")): "K",
                            frozenset(("A", "G")): "R",
                            frozenset(("C", "T")): "Y"}

        transition_repeat = []
        for index, char in enumerate(self.gaped_matching_pattern):
            column = [r[index] for r in self.list_gaped_repeats]
            column = set(column)
            difference = column - {char}
            if len(difference) == 1:
                alternative_char = list(difference)[0]
                if frozenset((alternative_char, char)) in dict_transitions:
                    transition_char = dict_transitions[frozenset((alternative_char, char))]
                    transition_repeat.append(transition_char)
                else:
                    transition_repeat.append(char)
            else:
                transition_repeat.append(char)

        transition_repeat = "".join(transition_repeat)
        transition_repeat = transition_repeat.replace(" ", "")
        self.transition_consensus = transition_repeat

    @staticmethod
    def apply_dots(repeat_with_gaps, matching_pattern_with_gaps):
        dotted_repeat = ''
        if len(repeat_with_gaps) != len(matching_pattern_with_gaps):
            raise ValueError
        for char_repeat, char_pattern in zip(repeat_with_gaps, matching_pattern_with_gaps):
            if char_repeat == char_pattern:
                if char_repeat == ' ':
                    dotted_repeat += char_repeat
                else:
                    dotted_repeat += '.'
            else:
                dotted_repeat += char_repeat

        return dotted_repeat

    def create_dot_representation(self):
        dot_repr = ""
        max_len_repeat_start = max(len(str(x)) for x in self.list_repeat_absolute_starts)
        max_len_spacer = max(len(x) for x in self.list_spacers)
        padding_size_repeat = max_len_repeat_start + 3
        padding_size_spacer = max_len_spacer + 3

        for repeat_start, repeat, spacer, errors in zip(self.list_repeat_absolute_starts,
                                                        self.list_dotted_gaped_repeats,
                                                        self.list_spacers,
                                                        self.fuzzy_counts):
            repeat_start += 1
            line_padding_repeat = " " * (padding_size_repeat - len(str(repeat_start)))
            line_padding_spacer = " " * (padding_size_spacer - len(str(spacer)))
            substitutions, insertions, deletions = errors
            error_line = f"s:{substitutions} i:{insertions} d:{deletions}"
            line = f"{repeat_start}{line_padding_repeat}{repeat}   {spacer}{line_padding_spacer}{error_line}\n"
            dot_repr += line

        line_padding_last_repeat = " " * (padding_size_repeat - len(str(self.list_repeat_absolute_starts[-1])))
        last_repeat_start = self.list_repeat_absolute_starts[-1]
        last_repeat = self.list_dotted_gaped_repeats[-1]
        gap_instead_spacer = " " * padding_size_spacer
        l_substitutions, l_insertions, l_deletions = self.fuzzy_counts[-1]
        last_errors = f"s:{l_substitutions} i:{l_insertions} d:{l_deletions}"
        last_repeat_start += 1
        last_line = f"{last_repeat_start}{line_padding_last_repeat}{last_repeat}   {gap_instead_spacer}{last_errors}\n"
        dot_repr += last_line

        dot_repr += "_" * 100 + '\n'
        padding_consensus = " " * padding_size_repeat
        dot_repr += f"{padding_consensus}{self.gaped_matching_pattern}\n"

        dot_repr += "_" * 100 + '\n'
        padding_consensus_no_gap = " " * padding_size_repeat
        dot_repr += f"{padding_consensus_no_gap}{self.consensus}\n"

        dot_repr += "_" * 100 + '\n'
        padding_consensus_transitions = " " * padding_size_repeat
        dot_repr += f"{padding_consensus_transitions}{self.transition_consensus}\n"
        return dot_repr












def create_boundaries_for_intervals(intervals, flank_size, source_length=None):
    if not intervals:
        return []
    if (
            not isinstance(flank_size, int)
            or isinstance(flank_size, bool)
            or flank_size < 0):
        raise ValueError("flank size must be a non-negative integer")
    if source_length is not None and (
            not isinstance(source_length, int)
            or isinstance(source_length, bool)
            or source_length < 0):
        raise ValueError("source length must be a non-negative integer")
    for interval in intervals:
        if (
                not isinstance(interval, (tuple, list))
                or len(interval) != 2
                or not isinstance(interval[0], int)
                or isinstance(interval[0], bool)
                or not isinstance(interval[1], int)
                or isinstance(interval[1], bool)
                or interval[0] < 0
                or interval[1] <= interval[0]
                or (
                    source_length is not None
                    and interval[1] > source_length
                )):
            raise ValueError(
                "array intervals must be valid source-backed half-open spans"
            )

    indexed_intervals = sorted(
        enumerate(intervals),
        key=lambda item: (item[1][0], item[1][1], item[0]),
    )
    sorted_intervals = [interval for _, interval in indexed_intervals]

    starts = [max(0, sorted_intervals[0][0] - flank_size)]
    ends = []
    for first, second in zip(sorted_intervals[:-1], sorted_intervals[1:]):
        end_first = first[1]
        begin_second = second[0]
        if end_first >= begin_second:
            ends.append(end_first)
            starts.append(begin_second)
        elif (end_first + 2 * flank_size) < begin_second:
            ends.append(end_first + flank_size)
            starts.append(begin_second - flank_size)
        else:
            possible_flank = int((begin_second - end_first) / 2)
            ends.append(end_first + possible_flank)
            starts.append(begin_second - possible_flank)

    ends.append(sorted_intervals[-1][1] + flank_size)
    sorted_boundaries = list(zip(starts, ends))

    if source_length is not None:
        sorted_boundaries = [
            (min(start, source_length), min(end, source_length))
            for start, end in sorted_boundaries
        ]

    boundaries = [None] * len(intervals)
    for (original_index, _), boundary in zip(
            indexed_intervals, sorted_boundaries):
        boundaries[original_index] = boundary
    return boundaries

#      Refine regex result (insertion instead of substitution)
#######################################################################
#######################################################################


class ArrayRefinerInsertionsDeletions:
    def __init__(self, crispr_candidate):
        self.crispr_candidate = crispr_candidate
        self.refined_candidate = crispr_candidate

        self._refine_insertions()

    def _refine_insertions(self):
        ari = ArrayRefinerInsertions(self.crispr_candidate)
        self.refined_candidate = ari.refined_candidate

    def _refine_deletions(self):
        pass

    def output(self):
        return self.refined_candidate


class ArrayRefinerInsertions:
    def __init__(self, crispr_candidate):
        self.crispr_candidate = crispr_candidate
        self.refined_candidate = crispr_candidate

        self.original_repeats = self.crispr_candidate.list_repeats
        self.original_repeats_gaped = self.crispr_candidate.list_repeats_gaped
        self.original_spacers = self.crispr_candidate.list_spacers
        self.original_starts = self.crispr_candidate.list_repeat_starts

        self.problematic_spacers_indexes = []

        self._search_potential_problems()
        self._refine_insertions()

    def _search_potential_problems(self):
        len_spacers = [len(s) for s in self.original_spacers]
        for index, len_spacer in enumerate(len_spacers):

            if index == 0:
                if len(len_spacers) > 2:
                    window_spacers_len = len_spacers[1:3]
                    if len(set(window_spacers_len)) == 1:
                        if window_spacers_len[0] != len_spacer:
                            self.problematic_spacers_indexes.append(index)
                else:
                    self.problematic_spacers_indexes.append(index)

            if index == 1:
                if len(len_spacers) > 2:
                    window_spacers_len = len_spacers[0:1] + len_spacers[2:3]
                    if len(set(window_spacers_len)) == 1:
                        if window_spacers_len[0] != len_spacer:
                            self.problematic_spacers_indexes.append(index)
                else:
                    self.problematic_spacers_indexes.append(index)

            if index > 1:
                if index not in self.problematic_spacers_indexes:
                    if index < (len(len_spacers) - 2):
                        window_spacers_len = len_spacers[index - 2:index] + len_spacers[index + 1:index + 3]
                        if len(set(window_spacers_len)) == 1:
                            if window_spacers_len[0] != len_spacer:
                                self.problematic_spacers_indexes.append(index)

            if index == (len(len_spacers) - 2):
                if index not in self.problematic_spacers_indexes:
                    if len(len_spacers) > 2:
                        window_spacers_len = len_spacers[len(len_spacers) - 3:len(len_spacers) - 2] + [len_spacers[-1]]
                        if len(set(window_spacers_len)) == 1:
                            if window_spacers_len[0] != len_spacer:
                                self.problematic_spacers_indexes.append(index)
                    else:
                        self.problematic_spacers_indexes.append(index)

            if index == (len(len_spacers) - 1):
                if index not in self.problematic_spacers_indexes:
                    if len(len_spacers) > 2:
                        window_spacers_len = len_spacers[len(len_spacers) - 2:len(len_spacers) - 1]
                        if len(set(window_spacers_len)) == 1:
                            if window_spacers_len[0] != len_spacer:
                                self.problematic_spacers_indexes.append(index)
                    else:
                        self.problematic_spacers_indexes.append(index)

    def _refine_insertions(self):
        for spacer_index in self.problematic_spacers_indexes:
            osir = OneSpacerInsertionRefiner(self.refined_candidate, spacer_index)
            self.refined_candidate = osir.refined_candidate




class OneSpacerInsertionRefiner:
    def __init__(self, crispr_candidate, spacer_index):
        self.crispr_candidate = crispr_candidate
        self.spacer_index = spacer_index
        self.refined_candidate = crispr_candidate

        self.repeats = self.crispr_candidate.list_repeats
        self.repeats_gaped = self.crispr_candidate.list_repeats_gaped
        self.spacers = self.crispr_candidate.list_spacers
        self.starts = self.crispr_candidate.list_repeat_starts
        self.original_errors = self.crispr_candidate.total_mismatches

        self.indexes_insertions = None
        self.repeat_index_to_refine = None

        self._get_indexes_insertions()
        self._find_repeat_to_refine()

    def _get_indexes_insertions(self):
        self.indexes_insertions = list(set([index for repeat in self.repeats_gaped
                                            for index, char in enumerate(repeat) if char == ' ']))
        self.indexes_insertions = sorted(self.indexes_insertions)

    def _find_repeat_to_refine(self):
        self.left_repeat = self.repeats[self.spacer_index]
        self.right_repeat = self.repeats[self.spacer_index+1]

        if len(self.left_repeat) > len(self.right_repeat):
            self.repeat_index_to_refine = self.spacer_index
            self._refine_left_repeat()
        if len(self.left_repeat) < len(self.right_repeat):
            self.repeat_index_to_refine = self.spacer_index + 1
            self._refine_right_repeat()

    def _candidate_with_recomputed_alignment(
            self, new_repeats, new_spacers, new_starts):
        """Build a structurally valid alignment from source-backed repeats.

        Alignment columns are presentation coordinates, not DNA coordinates.
        In particular, column zero or the final column can be an alignment
        gap.  After moving a DNA base between a repeat and its spacer, rebuild
        every alignment row from the ungapped repeats instead of splicing the
        old rows at either end.
        """
        relative_error_indexes = [
            _minimal_edit_changes(self.crispr_candidate.consensus, repeat)
            for repeat in new_repeats
        ]
        _, new_repeats_gaped = merge_fuzzy_alignments(
            self.crispr_candidate.consensus,
            new_repeats,
            relative_error_indexes,
        )
        refined = CrisprCandidate(
            list_repeats=new_repeats,
            list_repeats_gaped=new_repeats_gaped,
            list_spacers=new_spacers,
            list_repeat_starts=new_starts,
        )
        inherit_candidate_origins(refined, (self.crispr_candidate,))
        return refined

    def _refine_left_repeat(self):
        gaped_repeat = self.repeats_gaped[self.repeat_index_to_refine]
        for global_insertion_index in self.indexes_insertions[::-1]:
            if gaped_repeat[global_insertion_index] != " ":
                new_repeats = self.repeats[:]
                repeat_to_refine = new_repeats[self.repeat_index_to_refine]
                if len(repeat_to_refine) <= 1:
                    break
                moved_base = repeat_to_refine[-1]
                new_repeats[self.repeat_index_to_refine] = repeat_to_refine[:-1]

                new_spacers = self.spacers[:]
                new_spacers[self.spacer_index] = (
                    moved_base + new_spacers[self.spacer_index]
                )

                new_starts = self.starts[:]

                new_candidate = self._candidate_with_recomputed_alignment(
                    new_repeats, new_spacers, new_starts
                )

                new_candidate_mismatches = new_candidate.total_mismatches
                if new_candidate_mismatches <= self.original_errors:
                    self.refined_candidate = new_candidate

                break

    def _refine_right_repeat(self):
        gaped_repeat = self.repeats_gaped[self.repeat_index_to_refine]
        for global_insertion_index in self.indexes_insertions:
            if gaped_repeat[global_insertion_index] != " ":
                new_repeats = self.repeats[:]
                repeat_to_refine = new_repeats[self.repeat_index_to_refine]
                if len(repeat_to_refine) <= 1:
                    break
                moved_base = repeat_to_refine[0]
                new_repeats[self.repeat_index_to_refine] = repeat_to_refine[1:]

                new_spacers = self.spacers[:]
                new_spacers[self.spacer_index] = (
                    new_spacers[self.spacer_index] + moved_base
                )

                new_starts = self.starts[:]
                new_starts[self.repeat_index_to_refine] = new_starts[self.repeat_index_to_refine] + 1

                new_candidate = self._candidate_with_recomputed_alignment(
                    new_repeats, new_spacers, new_starts
                )

                new_candidate_mismatches = new_candidate.total_mismatches
                if new_candidate_mismatches <= self.original_errors:
                    self.refined_candidate = new_candidate

                break
