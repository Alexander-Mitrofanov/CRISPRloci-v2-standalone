"""Copy-on-write demotion of accepted low-diversity spacer arrays.

The historical pipeline has five internal category slots and many downstream
contracts depend on that shape.  Tandem-like rows therefore move from an
accepted slot into the existing public ``Low score`` slot.  Their metadata
records the public tag ``tandem`` and the distinct scientific disposition
``likely_tandem``.
"""

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from numbers import Real

from components.likely_tandem import (
    CATEGORY_NAME,
    LIKELY_TANDEM_POLICY,
    LikelyTandemEvidence,
    LikelyTandemPolicy,
    evaluate_likely_tandem,
)


CLASSIFICATION_BASIS_METHOD = "low_spacer_diversity_tandem_demotion"
TRANSACTION_RULE_NAME = "likely-tandem-category-quarantine"
TRANSACTION_RULE_VERSION = "1.1.0"
TANDEM_TAG = "tandem"
_CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
_ACCEPTED_CATEGORY_INDEXES = (0, 2)
_TARGET_CATEGORY_INDEX = 4


@dataclass(frozen=True)
class LikelyTandemQuarantineRecord:
    """One immutable accepted-row quarantine decision."""

    action: str
    disposition: object
    tag: str
    source_category_index: int
    source_category_name: str
    target_category_index: object
    target_category_name: object
    locus_key_json: str
    row_ordinal: int
    score: float
    source_start_0based: int
    source_end_0based_exclusive: int
    repeat_count: int
    spacers_source_forward: tuple
    repeats_source_forward: tuple
    repeat_intervals_0based_half_open: tuple
    spacer_intervals_0based_half_open: tuple
    reconstructed_source_span_sha256: str
    evidence: LikelyTandemEvidence

    def as_dict(self):
        repeats = [
            {
                "index": index,
                "source_interval": {
                    "start": interval[0],
                    "end": interval[1],
                },
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(
                    sequence.encode("utf-8")
                ).hexdigest(),
            }
            for index, (sequence, interval) in enumerate(zip(
                self.repeats_source_forward,
                self.repeat_intervals_0based_half_open,
            ), 1)
        ]
        spacers = [
            {
                "index": index,
                "source_interval": {
                    "start": interval[0],
                    "end": interval[1],
                },
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(
                    sequence.encode("utf-8")
                ).hexdigest(),
            }
            for index, (sequence, interval) in enumerate(zip(
                self.spacers_source_forward,
                self.spacer_intervals_0based_half_open,
            ), 1)
        ]
        return {
            "action": self.action,
            "disposition": self.disposition,
            "tag": self.tag,
            "source_category_index": self.source_category_index,
            "source_category_name": self.source_category_name,
            "target_category_index": self.target_category_index,
            "target_category_name": self.target_category_name,
            "locus_key": json.loads(self.locus_key_json),
            "row_ordinal": self.row_ordinal,
            "score": self.score,
            "source_interval": {
                "coordinate_system": "one_based_inclusive",
                "start": self.source_start_0based + 1,
                "end": self.source_end_0based_exclusive,
                "length": (
                    self.source_end_0based_exclusive
                    - self.source_start_0based
                ),
            },
            "repeat_count": self.repeat_count,
            "spacers_source_forward": list(
                self.spacers_source_forward
            ),
            "source_architecture": {
                "coordinate_system": "zero_based_half_open",
                "reference": "input_forward_strand",
                "repeats": repeats,
                "spacers": spacers,
                "reconstructed_source_span_sha256": (
                    self.reconstructed_source_span_sha256
                ),
            },
            "evidence": self.evidence.as_dict(),
        }


@dataclass(frozen=True)
class LikelyTandemQuarantineStats:
    """Cardinality and conservation audit for one transaction."""

    total_rows_before: int
    total_rows_after: int
    accepted_rows_reviewed: int
    accepted_rows_before: int
    accepted_rows_after: int
    likely_tandem_count: int
    retained_count: int
    bona_fide_quarantined: int
    possible_quarantined: int
    diagnostic_rows_before: int
    diagnostic_rows_after: int

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class _AcceptedLocation:
    category_index: int
    locus_key: object
    locus_key_json: str
    row_ordinal: int
    row: object


def _json_copy(value, label):
    try:
        return json.loads(json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ))
    except (TypeError, ValueError) as error:
        raise ValueError("{} must be JSON-safe".format(label)) from error


def _locus_key_json(locus_key):
    try:
        return json.dumps(
            locus_key,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError) as error:
        raise ValueError(
            "likely-tandem locus keys must be JSON-safe"
        ) from error


def _validate_row(row):
    if not isinstance(row, (list, tuple)) or len(row) not in (2, 3, 4):
        raise ValueError(
            "category rows must contain score, candidate, optional features, "
            "and optional metadata"
        )
    score = row[0]
    if (
            isinstance(score, bool)
            or not isinstance(score, Real)
            or not math.isfinite(float(score))):
        raise ValueError("category row score must be finite")
    candidate = row[1]
    spacers = getattr(candidate, "list_spacers", None)
    if (
            isinstance(spacers, (str, bytes))
            or not isinstance(spacers, (list, tuple))
            or any(not isinstance(spacer, str) for spacer in spacers)):
        raise ValueError("category row candidate has invalid spacers")
    if len(row) == 4 and not isinstance(row[3], dict):
        raise ValueError("category row metadata must be a mapping")


def _candidate_audit_geometry(candidate):
    repeats = getattr(candidate, "list_repeats", None)
    starts = getattr(candidate, "list_repeat_starts", None)
    spacers = tuple(candidate.list_spacers)
    if (
            not isinstance(repeats, (list, tuple))
            or not repeats
            or any(not isinstance(repeat, str) or not repeat
                   for repeat in repeats)
            or not isinstance(starts, (list, tuple))
            or len(starts) != len(repeats)
            or any(isinstance(start, bool) or not isinstance(start, int)
                   or start < 0 for start in starts)
            or len(spacers) != len(repeats) - 1):
        raise ValueError(
            "likely-tandem candidate has invalid source geometry"
        )
    if any(
            starts[index] >= starts[index + 1]
            for index in range(len(starts) - 1)):
        raise ValueError(
            "likely-tandem candidate must use source-forward repeat order"
        )
    source_start = starts[0]
    source_end = starts[-1] + len(repeats[-1])
    if source_end <= source_start:
        raise ValueError(
            "likely-tandem candidate has an empty source interval"
        )
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    spacer_intervals = tuple(
        (repeat_intervals[index][1], starts[index + 1])
        for index in range(len(spacers))
    )
    if any(
            end < start or end - start != len(spacer)
            for spacer, (start, end) in zip(spacers, spacer_intervals)):
        raise ValueError(
            "likely-tandem candidate spacers do not reconstruct source geometry"
        )
    reconstructed_parts = []
    for index, repeat in enumerate(repeats):
        reconstructed_parts.append(repeat)
        if index < len(spacers):
            reconstructed_parts.append(spacers[index])
    reconstructed = "".join(reconstructed_parts)
    if len(reconstructed) != source_end - source_start:
        raise ValueError(
            "likely-tandem candidate components do not fill the source span"
        )
    return (
        source_start,
        source_end,
        len(repeats),
        spacers,
        tuple(repeats),
        repeat_intervals,
        spacer_intervals,
        hashlib.sha256(reconstructed.encode("utf-8")).hexdigest(),
    )


def _index_categories(categories):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("likely-tandem quarantine requires five categories")
    locations = []
    row_ids = set()
    total = 0
    for category_index, category in enumerate(categories):
        if not isinstance(category, dict):
            raise ValueError("every category must be a mapping")
        for locus_key, rows in category.items():
            locus_json = _locus_key_json(locus_key)
            if (
                    isinstance(rows, (str, bytes))
                    or not isinstance(rows, (list, tuple))
                    or not rows):
                raise ValueError("category loci must contain non-empty rows")
            for row_ordinal, row in enumerate(rows):
                _validate_row(row)
                if id(row) in row_ids:
                    raise ValueError(
                        "one scored row occurs more than once in categories"
                    )
                row_ids.add(id(row))
                total += 1
                if category_index in _ACCEPTED_CATEGORY_INDEXES:
                    locations.append(_AcceptedLocation(
                        category_index=category_index,
                        locus_key=locus_key,
                        locus_key_json=locus_json,
                        row_ordinal=row_ordinal,
                        row=row,
                    ))
    return tuple(locations), total


def _quarantined_row(row, evidence, location):
    metadata = dict(row[3]) if len(row) == 4 else {}
    prior_tags = validated_row_tags(row)
    metadata["tags"] = sorted(set(prior_tags + (TANDEM_TAG,)))
    prior_basis = metadata.get("classification_basis")
    if prior_basis is not None:
        prior_basis = _json_copy(prior_basis, "prior classification_basis")
    metadata["classification_basis"] = {
        "method": CLASSIFICATION_BASIS_METHOD,
        "disposition": CATEGORY_NAME,
        "tag": TANDEM_TAG,
        "prior_classification_basis": prior_basis,
        "likely_tandem": evidence.as_dict(),
        "transition": {
            "rule_name": TRANSACTION_RULE_NAME,
            "rule_version": TRANSACTION_RULE_VERSION,
            "source_category_index": location.category_index,
            "source_category_name": _CATEGORY_NAMES[location.category_index],
            "target_category_index": _TARGET_CATEGORY_INDEX,
            "target_category_name": _CATEGORY_NAMES[
                _TARGET_CATEGORY_INDEX
            ],
            "reported_disposition": CATEGORY_NAME,
            "reported_tag": TANDEM_TAG,
            "locus_key": json.loads(location.locus_key_json),
            "row_ordinal": location.row_ordinal,
        },
    }
    prefix = list(row[:3])
    if len(prefix) == 2:
        prefix.append([])
    if isinstance(row, tuple):
        return tuple(prefix) + (metadata,)
    return prefix + [metadata]


def _count_rows(categories, indexes=None):
    selected = range(len(categories)) if indexes is None else indexes
    return sum(
        len(rows)
        for category_index in selected
        for rows in categories[category_index].values()
    )


def validated_row_tags(row):
    """Return deterministic, validated optional output tags for one row."""

    if not isinstance(row, (list, tuple)) or len(row) < 4:
        return ()
    metadata = row[3]
    if not isinstance(metadata, dict):
        raise ValueError("candidate metadata must be a mapping")
    raw_tags = metadata.get("tags")
    if raw_tags is None:
        return ()
    if isinstance(raw_tags, (str, bytes)) or not isinstance(
            raw_tags, (list, tuple)):
        raise ValueError("candidate tags must be a sequence")
    tags = tuple(raw_tags)
    if any(
            not isinstance(tag, str)
            or not tag
            or len(tag) > 64
            or tag[0] not in "abcdefghijklmnopqrstuvwxyz"
            or any(
                character not in (
                    "abcdefghijklmnopqrstuvwxyz0123456789_-"
                )
                for character in tag
            )
            for tag in tags):
        raise ValueError(
            "candidate tags must be lowercase identifier strings"
        )
    if len(tags) != len(set(tags)):
        raise ValueError("candidate tags must be unique")
    return tuple(sorted(tags))


def is_likely_tandem_row(row):
    """Return whether a scored row carries this transaction's disposition."""

    if not isinstance(row, (list, tuple)) or len(row) != 4:
        return False
    metadata = row[3]
    if not isinstance(metadata, dict):
        return False
    basis = metadata.get("classification_basis")
    return (
        isinstance(basis, dict)
        and basis.get("method") == CLASSIFICATION_BASIS_METHOD
        and basis.get("disposition") == CATEGORY_NAME
        and basis.get("tag") == TANDEM_TAG
        and TANDEM_TAG in validated_row_tags(row)
    )


def apply_likely_tandem_quarantine(
        categories,
        policy=LIKELY_TANDEM_POLICY):
    """Move qualifying accepted rows to the diagnostic slot atomically."""

    if not isinstance(policy, LikelyTandemPolicy):
        raise TypeError("policy must be a LikelyTandemPolicy")
    locations, total_before = _index_categories(categories)
    evaluated = tuple(
        (
            location,
            evaluate_likely_tandem(
                location.row[1].list_spacers,
                policy,
            ),
        )
        for location in locations
    )

    new_categories = [
        {locus_key: list(rows) for locus_key, rows in category.items()}
        for category in categories
    ]
    records = []
    quarantined_by_category = {0: 0, 2: 0}
    for location, evidence in evaluated:
        if evidence.likely_tandem:
            (
                source_start,
                source_end,
                repeat_count,
                spacers,
                repeats,
                repeat_intervals,
                spacer_intervals,
                reconstructed_sha256,
            ) = (
                _candidate_audit_geometry(location.row[1])
            )
            source_rows = new_categories[location.category_index][
                location.locus_key
            ]
            matches = [
                index for index, row in enumerate(source_rows)
                if row is location.row
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    "accepted row identity changed during likely-tandem "
                    "quarantine"
                )
            source_rows.pop(matches[0])
            if not source_rows:
                del new_categories[location.category_index][
                    location.locus_key
                ]
            new_categories[_TARGET_CATEGORY_INDEX].setdefault(
                location.locus_key, []
            ).append(_quarantined_row(location.row, evidence, location))
            action = "quarantined"
            disposition = CATEGORY_NAME
            target_index = _TARGET_CATEGORY_INDEX
            target_name = _CATEGORY_NAMES[_TARGET_CATEGORY_INDEX]
            quarantined_by_category[location.category_index] += 1
            records.append(LikelyTandemQuarantineRecord(
                action=action,
                disposition=disposition,
                tag=TANDEM_TAG,
                source_category_index=location.category_index,
                source_category_name=_CATEGORY_NAMES[
                    location.category_index
                ],
                target_category_index=target_index,
                target_category_name=target_name,
                locus_key_json=location.locus_key_json,
                row_ordinal=location.row_ordinal,
                score=float(location.row[0]),
                source_start_0based=source_start,
                source_end_0based_exclusive=source_end,
                repeat_count=repeat_count,
                spacers_source_forward=spacers,
                repeats_source_forward=repeats,
                repeat_intervals_0based_half_open=repeat_intervals,
                spacer_intervals_0based_half_open=spacer_intervals,
                reconstructed_source_span_sha256=reconstructed_sha256,
                evidence=evidence,
            ))
        else:
            action = "retained"
            disposition = None
            target_index = None
            target_name = None

    total_after = _count_rows(new_categories)
    accepted_before = len(locations)
    quarantined = sum(quarantined_by_category.values())
    accepted_after = _count_rows(
        new_categories, _ACCEPTED_CATEGORY_INDEXES
    )
    diagnostic_before = _count_rows(categories, (_TARGET_CATEGORY_INDEX,))
    diagnostic_after = _count_rows(
        new_categories, (_TARGET_CATEGORY_INDEX,)
    )
    if total_after != total_before:
        raise RuntimeError("likely-tandem quarantine changed row cardinality")
    if accepted_after != accepted_before - quarantined:
        raise RuntimeError(
            "likely-tandem quarantine changed accepted cardinality incorrectly"
        )
    if diagnostic_after != diagnostic_before + quarantined:
        raise RuntimeError(
            "likely-tandem quarantine changed diagnostic cardinality incorrectly"
        )
    stats = LikelyTandemQuarantineStats(
        total_rows_before=total_before,
        total_rows_after=total_after,
        accepted_rows_reviewed=accepted_before,
        accepted_rows_before=accepted_before,
        accepted_rows_after=accepted_after,
        likely_tandem_count=quarantined,
        retained_count=accepted_before - quarantined,
        bona_fide_quarantined=quarantined_by_category[0],
        possible_quarantined=quarantined_by_category[2],
        diagnostic_rows_before=diagnostic_before,
        diagnostic_rows_after=diagnostic_after,
    )
    return new_categories, tuple(records), stats


def resolve_likely_tandem_audit_inputs(
        categories,
        records=None,
        stats=None,
        source_categories=None,
        accepted_only=False):
    """Make public writers ledger-complete when audit inputs are omitted.

    Pipeline callers supply all three values from the pre-quarantine
    transaction. Direct writer callers that omit all three are evaluated from
    their provided categories: retained-only inputs receive a mandatory empty
    ledger, while an unquarantined qualifying accepted row is detected by the
    report validation fence and publication fails.
    """

    supplied = tuple(
        value is not None
        for value in (records, stats, source_categories)
    )
    if any(supplied):
        if not all(supplied):
            raise ValueError(
                "likely-tandem records, stats, and source categories must "
                "be supplied together"
            )
        return records, stats, source_categories
    if accepted_only:
        raise ValueError(
            "accepted-only writers require explicit likely-tandem records, "
            "stats, and pre-quarantine source categories"
        )
    _projected, records, stats = apply_likely_tandem_quarantine(categories)
    return records, stats, categories


__all__ = [
    "CLASSIFICATION_BASIS_METHOD",
    "LikelyTandemQuarantineRecord",
    "LikelyTandemQuarantineStats",
    "TANDEM_TAG",
    "TRANSACTION_RULE_NAME",
    "TRANSACTION_RULE_VERSION",
    "apply_likely_tandem_quarantine",
    "is_likely_tandem_row",
    "resolve_likely_tandem_audit_inputs",
    "validated_row_tags",
]
