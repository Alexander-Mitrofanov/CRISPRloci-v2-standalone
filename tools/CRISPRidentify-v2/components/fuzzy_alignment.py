"""Coordinate normalization for fuzzy literal matches.

The frozen native matcher and regex 2019 expose ``fuzzy_changes`` in the
legacy coordinate convention:

* substitutions and insertions use absolute offsets into the matched text;
* deletions use the text boundary immediately before the deleted pattern
  character (so a leading deletion is reported at ``match.start() - 1``).

Modern regex 2024 reports the deletion boundary one position later.  The
runtime pins and explicitly handles both tested conventions before converting
them to the one alignment-column system used by the rest of the pipeline.

The rest of CRISPRidentify needs one coordinate system instead: zero-based
columns in the final pairwise alignment, after deletion gaps have been added
to the matched text.  This module is the sole conversion point between those
conventions and also merges pairwise alignments without conflating insertion
columns from different pattern boundaries.
"""

import regex


LEGACY_REGEX_VERSION = "2.5.27"
MODERN_REGEX_VERSION = "2.5.148"


def snapshot_owned_fuzzy_match_state(fuzzy_counts, fuzzy_changes):
    """Return the direct immutable owned-match state, or ``None``.

    The fast snapshot is intentionally restricted to the exact containers
    constructed by CRISPRidentify's three match adapters. Any foreign or
    malformed shape retains the recursive cache snapshot path.
    """

    if (
            type(fuzzy_counts) is not tuple
            or len(fuzzy_counts) != 3
            or any(type(value) is not int for value in fuzzy_counts)
            or type(fuzzy_changes) is not tuple
            or len(fuzzy_changes) != 3
            or any(type(positions) is not list for positions in fuzzy_changes)
            or any(
                type(position) is not int
                for positions in fuzzy_changes
                for position in positions
            )):
        return None
    return (
        fuzzy_counts,
        tuple(tuple(positions) for positions in fuzzy_changes),
    )


def _deletion_position_shift(match):
    """Return the tested raw-deletion-to-boundary adjustment."""

    explicit_shift = getattr(match, "deletion_position_shift", None)
    if explicit_shift is not None:
        return explicit_shift
    if regex.__version__ == LEGACY_REGEX_VERSION:
        return 1
    if regex.__version__ == MODERN_REGEX_VERSION:
        return 0
    raise RuntimeError(
        "unsupported regex fuzzy-coordinate convention for version {}"
        .format(regex.__version__)
    )


def _alignment_column(text_offset, deletion_columns):
    """Map a match-relative text offset to its final alignment column."""
    column = text_offset
    for deletion_column in deletion_columns:
        if deletion_column > column:
            break
        column += 1
    return column


def normalize_fuzzy_changes(match):
    """Return substitutions, insertions, and deletions as alignment columns.

    ``match`` may be a ``regex.Match`` or the native match adapter; both
    provide ``start()`` and ``fuzzy_changes``.  Returned positions are
    match-relative, non-negative, and refer to columns in the same final
    gapped alignment:

    * an insertion column receives a space in the pattern;
    * a deletion column receives ``-`` in the matched text;
    * a substitution column contains two different bases.
    """
    match_start = match.start()
    raw_substitutions, raw_insertions, raw_deletions = match.fuzzy_changes

    deletion_shift = _deletion_position_shift(match)
    deletion_columns = sorted(
        position - match_start + deletion_shift
        for position in raw_deletions
    )

    def normalize_text_positions(positions):
        return [
            _alignment_column(position - match_start, deletion_columns)
            for position in positions
        ]

    return [
        normalize_text_positions(raw_substitutions),
        normalize_text_positions(raw_insertions),
        deletion_columns,
    ]

