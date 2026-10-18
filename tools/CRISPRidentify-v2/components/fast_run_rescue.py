"""Bounded, classifier-routed search expansion for ``fast_run``.

The routing score is deliberately not an acceptance threshold.  It only
decides which unresolved source components receive the legacy 4x4 terminal
trim search.  Every resulting candidate still passes through the ordinary
structural filter and compact classifier.
"""

from collections import deque

from components.feature_schema import COMPACT_FEATURE_SCHEMA


ROUTING_POLICY_VERSION = "unresolved_compact_v2"
ROUTING_MIN_SCORE = 0.365
ROUTING_MIN_REPEAT_SIMILARITY = 0.80
ROUTING_MAX_MISMATCHES = 3.0
ROUTING_MIN_AVG_SPACER = 20.0
ROUTING_MAX_AVG_SPACER = 54.0
MAX_TRIGGERED_COMPONENTS = 32
MAX_DELTA_PATTERNS = 4096


def candidate_source_interval(candidate):
    """Return one validated zero-based, half-open candidate envelope."""

    starts = getattr(candidate, "list_repeat_starts", None)
    repeats = getattr(candidate, "list_repeats", None)
    if not isinstance(starts, (tuple, list)) or not isinstance(
            repeats, (tuple, list)):
        raise ValueError("candidate source geometry must use sequences")
    if not starts or len(starts) != len(repeats):
        raise ValueError("candidate repeats and starts must be non-empty peers")
    previous_end = None
    for start, repeat in zip(starts, repeats):
        if (
                not isinstance(start, int)
                or isinstance(start, bool)
                or start < 0
                or not isinstance(repeat, str)
                or not repeat):
            raise ValueError("candidate contains invalid source geometry")
        end = start + len(repeat)
        if previous_end is not None and start < previous_end:
            raise ValueError("candidate repeats must be monotonic and disjoint")
        previous_end = end
    return starts[0], previous_end


def _classification_eligible(candidate):
    eligible = getattr(candidate, "source_path_classification_eligible", True)
    if not isinstance(eligible, bool):
        raise ValueError("source path classification eligibility must be boolean")
    return eligible


def compact_feature_values(scored_row):
    """Map a compact score row to the public, ordered feature schema."""

    if not isinstance(scored_row, (tuple, list)) or len(scored_row) < 3:
        raise ValueError("fast rescue received a malformed score row")
    vectors = scored_row[2]
    if not isinstance(vectors, (tuple, list)) or len(vectors) != 1:
        raise ValueError("compact score row must contain exactly one vector")
    vector = vectors[0]
    try:
        values = tuple(float(value) for value in vector.reshape(-1))
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError("compact score vector is malformed") from error
    if len(values) != len(COMPACT_FEATURE_SCHEMA):
        raise ValueError("compact score vector has the wrong feature width")
    return dict(zip(COMPACT_FEATURE_SCHEMA, values))


def _strict_overlap_components(indexed_rows):
    """Partition rows by transitive, strict source-envelope overlap."""

    ordered = sorted(
        indexed_rows,
        key=lambda item: (item[0][0], item[0][1], item[1]),
    )
    if not ordered:
        return ()
    components = []
    current = [ordered[0]]
    current_end = ordered[0][0][1]
    for record in ordered[1:]:
        interval = record[0]
        if interval[0] < current_end:
            current.append(record)
            current_end = max(current_end, interval[1])
        else:
            components.append(tuple(current))
            current = [record]
            current_end = interval[1]
    components.append(tuple(current))
    return tuple(components)


def select_unresolved_components(scored_rows, structural_filter):
    """Select unresolved components eligible for bounded search expansion.

    Rows below the routing cutoff cannot trigger or bridge a component.  This
    keeps weak background geometry from joining otherwise independent loci.
    The returned records are ordered by descending trigger score, then source
    coordinates, and capped deterministically.
    """

    if not isinstance(scored_rows, dict):
        raise ValueError("scored rows must be a mapping")
    if not callable(structural_filter):
        raise ValueError("structural filter must be callable")

    records = []
    stats = {
        "policy": ROUTING_POLICY_VERSION,
        "scored_rows": 0,
        "classification_ineligible": 0,
        "structurally_invalid": 0,
        "below_routing_score": 0,
        "resolved_components": 0,
        "unresolved_components": 0,
        "failed_repeat_similarity": 0,
        "failed_mismatches": 0,
        "failed_min_avg_spacer": 0,
        "failed_avg_spacer": 0,
        "triggered_before_cap": 0,
        "triggered_components": 0,
        "component_cap_skipped": 0,
    }
    for locus_key, rows in scored_rows.items():
        indexed = []
        for row_index, row in enumerate(rows):
            stats["scored_rows"] += 1
            if not isinstance(row, (tuple, list)) or len(row) < 2:
                raise ValueError("fast rescue received a malformed score row")
            candidate = row[1]
            if not _classification_eligible(candidate):
                stats["classification_ineligible"] += 1
                continue
            if not structural_filter(candidate):
                stats["structurally_invalid"] += 1
                continue
            score = float(row[0])
            if score < ROUTING_MIN_SCORE:
                stats["below_routing_score"] += 1
                continue
            indexed.append((candidate_source_interval(candidate), row_index, row))

        for component in _strict_overlap_components(indexed):
            if any(float(record[2][0]) >= 0.5 for record in component):
                stats["resolved_components"] += 1
                continue
            stats["unresolved_components"] += 1
            top = max(component, key=lambda record: float(record[2][0]))
            features = compact_feature_values(top[2])
            if features["repeat_similarity"] < ROUTING_MIN_REPEAT_SIMILARITY:
                stats["failed_repeat_similarity"] += 1
                continue
            if features["number_mismatches"] > ROUTING_MAX_MISMATCHES:
                stats["failed_mismatches"] += 1
                continue
            if features["avg_spacer_len"] < ROUTING_MIN_AVG_SPACER:
                stats["failed_min_avg_spacer"] += 1
                continue
            if features["avg_spacer_len"] > ROUTING_MAX_AVG_SPACER:
                stats["failed_avg_spacer"] += 1
                continue
            records.append({
                "locus_key": locus_key,
                "source_interval": (
                    min(record[0][0] for record in component),
                    max(record[0][1] for record in component),
                ),
                "trigger_score": float(top[2][0]),
                "trigger_candidate": top[2][1],
                "features": features,
            })

    records.sort(key=lambda record: (
        -record["trigger_score"],
        record["source_interval"][0],
        record["source_interval"][1],
        repr(record["locus_key"]),
    ))
    stats["triggered_before_cap"] = len(records)
    selected = tuple(records[:MAX_TRIGGERED_COMPONENTS])
    stats["triggered_components"] = len(selected)
    stats["component_cap_skipped"] = len(records) - len(selected)
    return selected, stats


def _validated_repeat_group(repeats, trim_profile):
    if not isinstance(repeats, (tuple, list)):
        raise ValueError("repeat seeds must be a sequence")
    if any(not isinstance(repeat, str) or not repeat for repeat in repeats):
        raise ValueError("repeat seeds must be non-empty strings")
    if trim_profile not in ("full", "total3", "one_sided"):
        raise ValueError("fast rescue trim profile is invalid")
    return tuple(repeats)


def _iter_start_end_repeat_variants(repeats, trim_profile):
    originals = set(repeats)
    seen = set(originals)
    offsets = sorted(
        (
            (left, right)
            for left in range(4)
            for right in range(4)
            if left or right
            if trim_profile != "total3" or left + right <= 3
            if trim_profile != "one_sided" or not (left and right)
        ),
        key=lambda pair: (sum(pair), pair[0], pair[1]),
    )
    # Offset-first traversal gives every seed its smallest phase correction
    # before spending budget on deeper trims of an earlier seed.
    for left, right in offsets:
        for repeat in repeats:
            end = len(repeat) - right
            if left >= end:
                continue
            variant = repeat[left:end]
            if variant in seen:
                continue
            seen.add(variant)
            yield variant


def start_end_repeat_variants(
        repeats, maximum=MAX_DELTA_PATTERNS, trim_profile="full"):
    """Return stable unique 4x4 trim patterns absent from the seed domain."""

    if not isinstance(maximum, int) or isinstance(maximum, bool) or maximum < 1:
        raise ValueError("delta-pattern maximum must be a positive integer")
    repeats = _validated_repeat_group(repeats, trim_profile)

    variants = []
    for variant in _iter_start_end_repeat_variants(repeats, trim_profile):
        if len(variants) >= maximum:
            return tuple(variants), True
        variants.append(variant)
    return tuple(variants), False


def round_robin_start_end_repeat_variants(
        repeat_groups, maximum=MAX_DELTA_PATTERNS, trim_profile="full"):
    """Allocate one bounded terminal-variant domain fairly across clusters."""

    if not isinstance(maximum, int) or isinstance(maximum, bool) or maximum < 1:
        raise ValueError("delta-pattern maximum must be a positive integer")
    if not isinstance(repeat_groups, (tuple, list)):
        raise ValueError("repeat-seed groups must be a sequence")
    groups = tuple(
        _validated_repeat_group(repeats, trim_profile)
        for repeats in repeat_groups
    )
    iterators = tuple(
        iter(_iter_start_end_repeat_variants(repeats, trim_profile))
        for repeats in groups
    )
    allocations = [[] for _group in groups]
    active = deque(range(len(groups)))
    allocated = 0
    while active and allocated < maximum:
        group_index = active.popleft()
        try:
            variant = next(iterators[group_index])
        except StopIteration:
            continue
        allocations[group_index].append(variant)
        allocated += 1
        active.append(group_index)
    truncated = False
    if allocated >= maximum:
        # A queued iterator may have yielded its final item exactly at the cap.
        # Look ahead until a real surplus variant is found or every group is
        # proved exhausted; the consumed surplus is intentionally not admitted.
        while active:
            group_index = active.popleft()
            try:
                next(iterators[group_index])
            except StopIteration:
                continue
            truncated = True
            break
    return tuple(tuple(values) for values in allocations), truncated


def candidate_score_fingerprint(candidate):
    """Return a fail-closed structural key for exact first-pass score reuse."""

    eligible = _classification_eligible(candidate)
    try:
        fields = (
            tuple(candidate.list_repeat_starts),
            tuple(candidate.list_repeats),
            tuple(candidate.list_repeats_gaped),
            tuple(candidate.list_spacers),
        )
    except (AttributeError, TypeError) as error:
        raise ValueError("candidate cannot be fingerprinted safely") from error
    return fields + (eligible,)
