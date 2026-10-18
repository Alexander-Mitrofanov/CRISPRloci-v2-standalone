"""Deterministic biological-locus grouping for canonical report arrays.

The detector's raw arrays remain the primary observations.  This module only
builds an optional, self-versioned extension that links sufficiently coherent
accepted arrays into biological loci; it never removes, rewrites, or invents
an array.

Consensus similarity uses an *ungapped containment* comparison.  The shorter
consensus is placed, without gaps and without overhangs, at every possible
offset in the longer consensus.  Identity is the largest matching-base
fraction over those placements, with the leftmost offset winning a tie.  The
comparison is deliberately not a local, circular, reverse-complement, or
indel-tolerant alignment.

Grouping is deterministic complete-linkage in source order: an array joins the
first existing group for which it is pairwise coherent with every member.
Consequently every multi-member locus is an all-pair-coherent clique, and a
transitive A--B--C chain cannot silently merge A with C.
"""

import hashlib
import json
import re
import statistics


ACCEPTED_CATEGORIES = frozenset(("Bona-fide", "Possible"))

MIN_REPEAT_COUNT = 3
MIN_CONSENSUS_LENGTH = 18
MIN_CONSENSUS_LENGTH_RATIO = 2.0 / 3.0
MIN_CONSENSUS_IDENTITY = 0.80
MIN_OPPOSITE_STRAND_IDENTITY = 0.95
MIN_INTERVAL_OVERLAP_CONSENSUS_MULTIPLIER = 1
MAX_PERIOD_DIFFERENCE_BP = 3
MAX_PERIOD_DIFFERENCE_FRACTION = 0.10

EXTENSION_SCHEMA_NAME = "CRISPRidentify-biological-loci"
EXTENSION_SCHEMA_VERSION = "1.0.0"

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_ACGT_PATTERN = re.compile(r"^[ACGT]+$")
_KNOWN_STRANDS = frozenset(("+", "-", None))


class BiologicalLocusError(ValueError):
    """Canonical report data is malformed or ambiguous."""


def _is_plain_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _source_metadata(source):
    if not isinstance(source, dict):
        raise BiologicalLocusError("source metadata must be a mapping")

    source_sha256 = source.get("sha256")
    if (
        not isinstance(source_sha256, str)
        or not _SHA256_PATTERN.fullmatch(source_sha256)
    ):
        raise BiologicalLocusError(
            "source.sha256 must be a lowercase hexadecimal SHA-256 digest"
        )

    source_length = source.get("length")
    if not _is_plain_int(source_length) or source_length < 0:
        raise BiologicalLocusError(
            "source.length must be a non-negative integer"
        )

    return source_sha256, source_length


def _interval(value, label, source_length, containing=None):
    if not isinstance(value, dict):
        raise BiologicalLocusError("{} must be a mapping".format(label))

    start = value.get("start")
    end = value.get("end")
    length = value.get("length")
    if (
        not _is_plain_int(start)
        or not _is_plain_int(end)
        or not _is_plain_int(length)
        or start < 1
        or end < start
        or length != end - start + 1
        or end > source_length
    ):
        raise BiologicalLocusError(
            "{} must be a valid one-based inclusive source interval".format(
                label
            )
        )

    if containing is not None:
        if start < containing["start"] or end > containing["end"]:
            raise BiologicalLocusError(
                "{} must be contained by its array interval".format(label)
            )

    return {"start": start, "end": end, "length": length}


def _validated_array(raw, source_length, position):
    if not isinstance(raw, dict):
        raise BiologicalLocusError(
            "array {} must be a mapping".format(position)
        )

    array_id = raw.get("id")
    if not isinstance(array_id, str) or not array_id:
        raise BiologicalLocusError(
            "array {} has no non-empty id".format(position)
        )

    category = raw.get("category")
    if not isinstance(category, str) or not category:
        raise BiologicalLocusError(
            "array {!r} has no category".format(array_id)
        )

    source_interval = _interval(
        raw.get("source_interval"),
        "array {!r} source_interval".format(array_id),
        source_length,
    )

    orientation = raw.get("orientation")
    if not isinstance(orientation, dict):
        raise BiologicalLocusError(
            "array {!r} orientation must be a mapping".format(array_id)
        )
    strand = orientation.get("strand")
    if strand not in _KNOWN_STRANDS:
        raise BiologicalLocusError(
            "array {!r} has an invalid orientation strand".format(array_id)
        )

    consensus = raw.get("consensus_repeat")
    if not isinstance(consensus, str):
        raise BiologicalLocusError(
            "array {!r} consensus_repeat must be a string".format(array_id)
        )
    consensus = consensus.upper()

    repeat_count = raw.get("repeat_count")
    repeats = raw.get("repeats")
    if (
        not _is_plain_int(repeat_count)
        or repeat_count < 1
        or not isinstance(repeats, list)
        or repeat_count != len(repeats)
    ):
        raise BiologicalLocusError(
            "array {!r} repeat_count disagrees with repeats".format(array_id)
        )

    repeat_starts = []
    previous_end = None
    for repeat_position, repeat in enumerate(repeats, 1):
        if not isinstance(repeat, dict):
            raise BiologicalLocusError(
                "array {!r} repeat {} must be a mapping".format(
                    array_id,
                    repeat_position,
                )
            )
        repeat_interval = _interval(
            repeat.get("source_interval"),
            "array {!r} repeat {} source_interval".format(
                array_id,
                repeat_position,
            ),
            source_length,
            containing=source_interval,
        )
        sequence = repeat.get("sequence")
        if (
            not isinstance(sequence, str)
            or len(sequence) != repeat_interval["length"]
        ):
            raise BiologicalLocusError(
                "array {!r} repeat {} sequence length disagrees with its "
                "source interval".format(array_id, repeat_position)
            )
        if previous_end is not None and repeat_interval["start"] <= previous_end:
            raise BiologicalLocusError(
                "array {!r} repeats must be ordered and non-overlapping".format(
                    array_id
                )
            )
        repeat_starts.append(repeat_interval["start"])
        previous_end = repeat_interval["end"]

    periods = [
        right - left for left, right in zip(repeat_starts, repeat_starts[1:])
    ]
    median_period = statistics.median(periods) if periods else None

    return {
        "id": array_id,
        "category": category,
        "source_interval": source_interval,
        "strand": strand,
        "consensus": consensus,
        "repeat_count": repeat_count,
        "repeat_starts": tuple(repeat_starts),
        "repeat_start_set": frozenset(repeat_starts),
        "median_period": median_period,
    }


def _ungapped_containment_identity(first, second):
    """Return ``(identity, offset)`` for shorter fully contained in longer.

    Both inputs must be non-empty strings.  On equal-length input ``first`` is
    treated as the shorter sequence.  Ties resolve to the smallest (leftmost)
    offset, which makes the result independent of dictionary/set iteration.
    """

    if len(first) <= len(second):
        shorter, longer = first, second
    else:
        shorter, longer = second, first

    best_matches = -1
    best_offset = 0
    for offset in range(len(longer) - len(shorter) + 1):
        matches = sum(
            left == right
            for left, right in zip(
                shorter,
                longer[offset:offset + len(shorter)],
            )
        )
        if matches > best_matches:
            best_matches = matches
            best_offset = offset

    return best_matches / float(len(shorter)), best_offset


def _interval_overlap(first, second):
    return max(
        0,
        min(first["end"], second["end"])
        - max(first["start"], second["start"])
        + 1,
    )


def _eligible(array):
    return (
        array["repeat_count"] >= MIN_REPEAT_COUNT
        and len(array["consensus"]) >= MIN_CONSENSUS_LENGTH
        and _ACGT_PATTERN.fullmatch(array["consensus"]) is not None
        and array["median_period"] is not None
    )


def _pair_is_coherent(first, second):
    if not _eligible(first) or not _eligible(second):
        return False

    shorter_length = min(len(first["consensus"]), len(second["consensus"]))
    longer_length = max(len(first["consensus"]), len(second["consensus"]))
    if shorter_length / float(longer_length) < MIN_CONSENSUS_LENGTH_RATIO:
        return False

    identity, _offset = _ungapped_containment_identity(
        first["consensus"],
        second["consensus"],
    )
    if identity < MIN_CONSENSUS_IDENTITY:
        return False

    if _interval_overlap(
        first["source_interval"],
        second["source_interval"],
    ) < shorter_length * MIN_INTERVAL_OVERLAP_CONSENSUS_MULTIPLIER:
        return False

    period_limit = max(
        MAX_PERIOD_DIFFERENCE_BP,
        (
            MAX_PERIOD_DIFFERENCE_FRACTION
            * min(first["median_period"], second["median_period"])
        ),
    )
    if abs(first["median_period"] - second["median_period"]) > period_limit:
        return False

    first_strand = first["strand"]
    second_strand = second["strand"]
    if first_strand is None or second_strand is None:
        return True
    if first_strand == second_strand:
        return True

    return (
        identity >= MIN_OPPOSITE_STRAND_IDENTITY
        and bool(first["repeat_start_set"] & second["repeat_start_set"])
    )


def _locus_id(source_sha256, member_ids):
    identity = {
        "source_sha256": source_sha256,
        "member_ids": sorted(member_ids),
    }
    encoded = json.dumps(
        identity,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return "CRISPR-LOCUS-{}".format(hashlib.sha256(encoded).hexdigest()[:16])


def _locus_record(source_sha256, members):
    members = sorted(
        members,
        key=lambda item: (
            item["source_interval"]["start"],
            item["source_interval"]["end"],
            item["id"],
        ),
    )
    member_ids = [member["id"] for member in members]
    start = min(member["source_interval"]["start"] for member in members)
    end = max(member["source_interval"]["end"] for member in members)
    return {
        "id": _locus_id(source_sha256, member_ids),
        "source_interval": {
            "start": start,
            "end": end,
            "length": end - start + 1,
        },
        "member_count": len(member_ids),
        "member_ids": member_ids,
    }


def build_biological_loci_extension(arrays, source):
    """Build a self-versioned locus extension from canonical report data.

    Only ``Bona-fide`` and ``Possible`` arrays become locus members.  Accepted
    arrays that are not eligible for pairwise merging remain singleton loci.
    Non-accepted arrays are absent from the extension and the input list is
    never mutated.
    """

    source_sha256, source_length = _source_metadata(source)
    if not isinstance(arrays, list):
        raise BiologicalLocusError("arrays must be a list")

    validated = [
        _validated_array(raw, source_length, position)
        for position, raw in enumerate(arrays)
    ]
    ids = [array["id"] for array in validated]
    if len(ids) != len(set(ids)):
        raise BiologicalLocusError("array ids must be unique")

    accepted = [
        array
        for array in validated
        if array["category"] in ACCEPTED_CATEGORIES
    ]
    accepted.sort(key=lambda item: (
        item["source_interval"]["start"],
        item["source_interval"]["end"],
        item["id"],
    ))

    groups = []
    for array in accepted:
        target = None
        if _eligible(array):
            for group in groups:
                if all(_pair_is_coherent(array, member) for member in group):
                    target = group
                    break
        if target is None:
            groups.append([array])
        else:
            target.append(array)

    loci = [
        _locus_record(source_sha256, group)
        for group in groups
    ]
    loci.sort(key=lambda locus: (
        locus["source_interval"]["start"],
        locus["source_interval"]["end"],
        locus["member_ids"],
        locus["id"],
    ))

    singleton_count = sum(
        locus["member_count"] == 1 for locus in loci
    )
    multi_member_count = len(loci) - singleton_count

    return {
        "schema": {
            "name": EXTENSION_SCHEMA_NAME,
            "version": EXTENSION_SCHEMA_VERSION,
        },
        "status": "computed",
        "source_sha256": source_sha256,
        "method": {
            "name": "source_order_complete_linkage_repeat_lattice_v1",
            "accepted_categories": sorted(ACCEPTED_CATEGORIES),
            "grouping": (
                "deterministic source-order complete-linkage; every "
                "multi-member locus is all-pair coherent"
            ),
            "consensus_alignment": (
                "maximum ungapped full containment identity of the shorter "
                "consensus in the longer; leftmost offset wins ties"
            ),
            "cross_origin_grouping": "disabled",
            "thresholds": {
                "minimum_repeat_count": MIN_REPEAT_COUNT,
                "minimum_acgt_consensus_length": MIN_CONSENSUS_LENGTH,
                "minimum_shorter_to_longer_consensus_length_ratio": (
                    MIN_CONSENSUS_LENGTH_RATIO
                ),
                "minimum_ungapped_containment_identity": (
                    MIN_CONSENSUS_IDENTITY
                ),
                "minimum_source_overlap_bp": (
                    "shorter_consensus_length"
                ),
                "maximum_median_repeat_period_difference_bp": (
                    "max(3, 0.10 * minimum_median_repeat_period)"
                ),
                "opposite_strand_minimum_identity": (
                    MIN_OPPOSITE_STRAND_IDENTITY
                ),
                "opposite_strand_requires_exact_shared_repeat_start": True,
            },
        },
        "counts": {
            "raw_accepted_array_count": len(accepted),
            "biological_locus_count": len(loci),
            "singleton_locus_count": singleton_count,
            "multi_member_locus_count": multi_member_count,
        },
        "loci": loci,
    }


__all__ = [
    "BiologicalLocusError",
    "build_biological_loci_extension",
]
