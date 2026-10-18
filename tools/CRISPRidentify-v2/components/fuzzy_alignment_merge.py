"""Validated construction and canonical merging of fuzzy alignments."""


def _introduce_alignment_gaps(sequence, positions, gap_character):
    positions = sorted(positions)
    if len(positions) != len(set(positions)):
        raise ValueError("fuzzy edit positions must be unique")

    characters = list(sequence)
    for position in positions:
        if position < 0 or position > len(characters):
            raise ValueError("fuzzy edit position is outside the alignment")
        characters.insert(position, gap_character)
    return "".join(characters)


def build_pairwise_fuzzy_alignment(pattern, matched_text, normalized_changes):
    """Build and validate one gapped pattern/match pair.

    Removing pattern spaces and match ``-`` characters reconstructs the two
    original inputs exactly.  Invalid or mixed coordinate conventions fail
    loudly instead of producing plausible-looking scientific output.
    """
    substitutions, insertions, deletions = normalized_changes
    gapped_pattern = _introduce_alignment_gaps(pattern, insertions, " ")
    gapped_match = _introduce_alignment_gaps(matched_text, deletions, "-")

    if len(gapped_pattern) != len(gapped_match):
        raise ValueError("fuzzy edit coordinates produce unequal alignments")

    observed = [[], [], []]
    for column, (pattern_base, match_base) in enumerate(
            zip(gapped_pattern, gapped_match)):
        if pattern_base == " " and match_base == "-":
            raise ValueError("fuzzy alignment contains a gap-gap column")
        if pattern_base == " ":
            observed[1].append(column)
        elif match_base == "-":
            observed[2].append(column)
        elif pattern_base != match_base:
            observed[0].append(column)

    expected = [
        sorted(substitutions), sorted(insertions), sorted(deletions)
    ]
    if observed != expected:
        raise ValueError("fuzzy edit coordinates do not reconstruct edits")

    return gapped_pattern, gapped_match


def merge_pairwise_fuzzy_alignments(pattern, pairwise_alignments):
    """Merge already validated pairwise alignments canonically.

    The accepted-only detector validates each distinct source match once and
    reuses that pair across sibling occurrence/boundary paths.  This helper is
    deliberately still defensive about the pattern projection, so a stale or
    malformed cached pair fails at the same alignment boundary as the ordinary
    construction path.
    """
    local_columns = []
    insertion_counts = {}
    for gapped_pattern, gapped_match in pairwise_alignments:
        if len(gapped_pattern) != len(gapped_match):
            raise ValueError("fuzzy edit coordinates produce unequal alignments")
        pattern_boundary = 0
        boundary_ordinals = {}
        columns = {}

        for pattern_base, match_base in zip(gapped_pattern, gapped_match):
            if pattern_base == " ":
                ordinal = boundary_ordinals.get(pattern_boundary, 0)
                token = ("insertion", pattern_boundary, ordinal)
                boundary_ordinals[pattern_boundary] = ordinal + 1
                insertion_counts[pattern_boundary] = max(
                    insertion_counts.get(pattern_boundary, 0), ordinal + 1
                )
            else:
                if (pattern_boundary >= len(pattern)
                        or pattern_base != pattern[pattern_boundary]):
                    raise ValueError("gapped alignment changed the pattern")
                token = ("base", pattern_boundary)
                pattern_boundary += 1
            columns[token] = match_base

        if pattern_boundary != len(pattern):
            raise ValueError("gapped alignment did not consume the pattern")
        local_columns.append(columns)

    global_tokens = []
    global_pattern = []
    for boundary in range(len(pattern) + 1):
        for ordinal in range(insertion_counts.get(boundary, 0)):
            global_tokens.append(("insertion", boundary, ordinal))
            global_pattern.append(" ")
        if boundary < len(pattern):
            global_tokens.append(("base", boundary))
            global_pattern.append(pattern[boundary])

    global_repeats = [
        "".join(columns.get(token, " ") for token in global_tokens)
        for columns in local_columns
    ]
    return "".join(global_pattern), global_repeats


def merge_fuzzy_alignments(pattern, matched_texts, normalized_changes):
    """Merge pairwise alignments using canonical pattern-boundary slots.

    A pairwise alignment column is not stable across matches: a deletion in
    one match shifts every later column.  Each inserted text base is therefore
    identified by ``(consumed pattern bases, ordinal at that boundary)``.
    The union of those slots defines the global alignment, while absent slots
    in a particular match are represented by spaces.
    """
    matched_texts = list(matched_texts)
    normalized_changes = list(normalized_changes)
    if len(matched_texts) != len(normalized_changes):
        raise ValueError("every fuzzy match must have one edit-coordinate set")
    pairwise_alignments = [
        build_pairwise_fuzzy_alignment(pattern, matched_text, changes)
        for matched_text, changes in zip(matched_texts, normalized_changes)
    ]
    return merge_pairwise_fuzzy_alignments(pattern, pairwise_alignments)


__all__ = [
    "build_pairwise_fuzzy_alignment",
    "merge_fuzzy_alignments",
    "merge_pairwise_fuzzy_alignments",
]
