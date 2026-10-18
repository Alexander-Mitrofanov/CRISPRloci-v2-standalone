"""Pure helpers for planning bounded fuzzy-repeat searches.

The helpers in this module do not execute a fuzzy search.  They only group
repeat candidates and calculate half-open ``(start, end)`` search regions.
Keeping this logic independent from the detection pipeline makes its ordering
and coverage invariants directly testable.

Upstream repeat collections are expected to contain unique, non-empty strings.
The detection layer must reconstruct final results in its original repeat
order after processing the families returned here.
"""


DEFAULT_PADDING = 15


def merge_regions(regions):
    """Return sorted, disjoint half-open regions.

    Overlapping and touching regions are merged.  Touching regions are kept in
    one search call because restarting ``regex.finditer`` can affect fuzzy
    alignment tie-breaking.  Input order and input region objects are not
    modified.
    """
    merged = []
    for start, end in sorted(regions):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(region) for region in merged]


def substring_families(repeats):
    """Group repeats under a longer repeat that contains every member.

    Candidates are considered in a stable, length-descending order.  Thus
    equal-length candidates retain their input order, families retain their
    creation order, and members retain their order from that stable pass.  The
    first member of every returned family contains every other member.
    """
    families = []
    for repeat in sorted(list(repeats), key=len, reverse=True):
        for family in families:
            if repeat in family[0]:
                family.append(repeat)
                break
        else:
            families.append([repeat])
    return families


def common_anchor(family):
    """Return the longest existing member contained in all family members.

    Ties retain family order.  ``None`` is returned for an empty family or
    when no existing family member is common to the entire family.
    """
    members = list(family)
    common = [
        candidate
        for candidate in members
        if all(candidate in repeat for repeat in members)
    ]
    return max(common, key=len) if common else None


def anchored_subfamilies(family):
    """Partition a family so every member contains its selected anchor.

    If one existing member anchors the complete family, a single group is
    returned.  Otherwise a deterministic greedy set cover is used: select the
    candidate covering the most remaining members, then prefer the longer
    candidate and finally its earlier family position.

    Invariants:

    * anchors are existing family members, never synthesized strings;
    * every input member is assigned exactly once;
    * every emitted member contains its group's anchor;
    * member order within a group follows the input family order;
    * the input family is not modified.

    The partition can interleave positions from the original family.  Callers
    that expose results must therefore restore the original detection order.
    """
    members = list(family)
    anchor = common_anchor(members)
    if anchor is not None:
        return [(anchor, list(members))]

    remaining = list(members)
    groups = []
    family_order = {
        pattern: index for index, pattern in enumerate(members)
    }

    while remaining:
        anchor = max(
            members,
            key=lambda candidate: (
                sum(candidate in pattern for pattern in remaining),
                len(candidate),
                -family_order[candidate],
            ),
        )
        anchored_members = [
            pattern for pattern in remaining if anchor in pattern
        ]
        groups.append((anchor, anchored_members))

        covered = set(anchored_members)
        remaining = [
            pattern for pattern in remaining if pattern not in covered
        ]

    return groups


def all_offsets(container, substring):
    """Return all substring offsets, including overlapping occurrences."""
    offsets = []
    offset = container.find(substring)
    while offset >= 0:
        offsets.append(offset)
        offset = container.find(substring, offset + 1)
    return offsets


def aligned_regions_from_spans(pattern, anchor, anchor_spans, text_length,
                               padding=DEFAULT_PADDING):
    """Map ``(start, end)`` anchor spans to regional pattern domains."""
    offsets = all_offsets(pattern, anchor)
    if not offsets:
        return [(0, text_length)]

    regions = []
    for match_start, match_end in anchor_spans:
        for offset in offsets:
            prefix_length = offset
            suffix_length = len(pattern) - offset - len(anchor)
            regions.append((
                max(0, match_start - prefix_length - padding),
                min(text_length, match_end + suffix_length + padding),
            ))
    return merge_regions(regions)




__all__ = [
    "DEFAULT_PADDING",
    "aligned_regions_from_spans",
    "all_offsets",
    "anchored_subfamilies",
    "common_anchor",
    "merge_regions",
    "substring_families",
]
