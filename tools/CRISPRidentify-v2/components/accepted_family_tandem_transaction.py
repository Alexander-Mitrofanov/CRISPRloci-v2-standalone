"""Transactional demotion of accepted coding-tandem array families.

The source-backed tandem validator is intentionally kept separate from
category mutation.  This module supplies the narrow copy-on-write bridge:
every accepted row must have exactly one explicitly complete family context,
every such row is independently revalidated against the original source, and
only a supported row is moved from Bona-fide/Possible to Possible Discarded.

Rows are selected by object identity, never by value equality.  A demoted row
is rebuilt only to attach classification provenance; its score, candidate, and
feature objects are retained unchanged.  All other row objects are retained
by identity in newly allocated category dictionaries and row lists.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass
import json
import math
from numbers import Real
from typing import Optional, Tuple

from components.accepted_family_tandem_validator import (
    AcceptedFamilyTandemDecision,
    validate_accepted_family_tandem,
)
from components.architecture_rescue_coding_context import SourceCDS


RULE_NAME = "accepted-family-tandem-category-demotion"
RULE_VERSION = "1.0.0"
COORDINATE_SYSTEM = "zero_based_half_open"

_ACCEPTED_CATEGORY_INDEXES = (0, 2)
_TARGET_CATEGORY_INDEX = 3
_CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)


@dataclass(frozen=True)
class AcceptedRowFamilyContext:
    """Identity-exact mapping from one accepted row to its complete family.

    ``complete_family_candidates`` is an exhaustive caller contract, not a
    search request.  It must be a duplicate-free tuple and must include the
    accepted row's candidate object itself.  Additional family candidates may
    be unscored source-generated paths, which is how a short accepted fragment
    can be checked against a longer family architecture.
    """

    accepted_row: object
    complete_family_candidates: Tuple[object, ...]


@dataclass(frozen=True)
class AcceptedFamilyTandemDemotionRecord:
    """One immutable, runtime-object-free accepted-row decision record."""

    action: str
    source_category_index: int
    source_category_name: str
    target_category_index: Optional[int]
    target_category_name: Optional[str]
    locus_key_json: str
    row_ordinal: int
    score: float
    repeat_count: int
    accepted_interval: Tuple[int, int]
    decision: AcceptedFamilyTandemDecision

    def as_dict(self):
        """Return a detached JSON-safe representation."""

        return {
            "action": self.action,
            "source_category_index": self.source_category_index,
            "source_category_name": self.source_category_name,
            "target_category_index": self.target_category_index,
            "target_category_name": self.target_category_name,
            "locus_key": json.loads(self.locus_key_json),
            "row_ordinal": self.row_ordinal,
            "score": self.score,
            "repeat_count": self.repeat_count,
            "accepted_interval": list(self.accepted_interval),
            "decision": _json_copy(self.decision.as_dict(), "tandem decision"),
        }


@dataclass(frozen=True)
class AcceptedFamilyTandemDemotionStats:
    """Immutable transaction cardinalities and conservation checks."""

    total_rows_before: int
    total_rows_after: int
    accepted_rows_reviewed: int
    tandem_supported: int
    tandem_not_supported: int
    bona_fide_demoted: int
    possible_demoted: int
    accepted_rows_before: int
    accepted_rows_after: int
    possible_discarded_rows_before: int
    possible_discarded_rows_after: int

    def as_dict(self):
        """Return a JSON-safe statistics mapping."""

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
        serialized = json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        return json.loads(serialized)
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
        raise ValueError("accepted-row locus keys must be JSON-safe") from error


def _validate_source_context(source, source_cds):
    if not isinstance(source, str):
        raise TypeError("source must be a string")
    if not source:
        raise ValueError("source must be non-empty")
    try:
        source.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("source must contain ASCII DNA symbols") from error
    if not isinstance(source_cds, tuple):
        raise TypeError("source_cds must be a tuple of SourceCDS objects")
    for index, cds in enumerate(source_cds):
        if not isinstance(cds, SourceCDS):
            raise TypeError("source_cds[{}] is not a SourceCDS".format(index))
        if cds.end > len(source):
            raise ValueError("source_cds[{}] lies outside the source".format(index))


def _validate_scored_row(row):
    if not isinstance(row, (list, tuple)) or len(row) not in (3, 4):
        raise ValueError(
            "category rows must contain score, candidate, features, "
            "and optional metadata"
        )
    score = row[0]
    if (
        isinstance(score, bool)
        or not isinstance(score, Real)
        or not math.isfinite(float(score))
    ):
        raise ValueError("category row scores must be finite numbers")
    if row[1] is None:
        raise ValueError("category rows must contain a candidate")
    if not isinstance(row[2], (list, tuple)):
        raise ValueError("category row features must be a sequence")
    if len(row) == 4:
        if not isinstance(row[3], Mapping):
            raise ValueError("category row metadata must be a mapping")
        prior = row[3].get("classification_basis")
        if prior is not None:
            if not isinstance(prior, Mapping):
                raise ValueError("prior classification_basis must be a mapping")
            _json_copy(prior, "prior classification_basis")


def _index_categories(categories):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ValueError("tandem demotion requires exactly five categories")

    row_ids = set()
    candidate_ids = set()
    accepted_locations = []
    total_rows = 0
    for category_index, category in enumerate(categories):
        if not isinstance(category, Mapping):
            raise ValueError("each category must be a mapping")
        for locus_key, rows in category.items():
            if not isinstance(rows, (list, tuple)) or not rows:
                raise ValueError(
                    "category locus values must be non-empty row sequences"
                )
            locus_key_json = (
                _locus_key_json(locus_key)
                if category_index in _ACCEPTED_CATEGORY_INDEXES
                else None
            )
            for row_ordinal, row in enumerate(rows):
                _validate_scored_row(row)
                row_id = id(row)
                if row_id in row_ids:
                    raise ValueError(
                        "one scored row occurs more than once in categories"
                    )
                row_ids.add(row_id)
                candidate_id = id(row[1])
                if candidate_id in candidate_ids:
                    raise ValueError("one candidate object has multiple category rows")
                candidate_ids.add(candidate_id)
                total_rows += 1
                if category_index in _ACCEPTED_CATEGORY_INDEXES:
                    accepted_locations.append(
                        _AcceptedLocation(
                            category_index=category_index,
                            locus_key=locus_key,
                            locus_key_json=locus_key_json,
                            row_ordinal=row_ordinal,
                            row=row,
                        )
                    )
    return tuple(accepted_locations), total_rows


def _index_complete_family_contexts(contexts, accepted_locations):
    if isinstance(contexts, (str, bytes)) or contexts is None:
        raise TypeError("accepted_row_families must be an iterable of contexts")
    try:
        contexts = tuple(contexts)
    except TypeError as error:
        raise TypeError(
            "accepted_row_families must be an iterable of contexts"
        ) from error

    accepted_by_id = {id(location.row): location for location in accepted_locations}
    indexed = {}
    for index, context in enumerate(contexts):
        if not isinstance(context, AcceptedRowFamilyContext):
            raise TypeError(
                "accepted_row_families[{}] is not an AcceptedRowFamilyContext".format(
                    index
                )
            )
        row_id = id(context.accepted_row)
        if row_id in indexed:
            raise ValueError("an accepted row has duplicate family contexts")
        location = accepted_by_id.get(row_id)
        if location is None or location.row is not context.accepted_row:
            raise ValueError("family context refers to a non-accepted or foreign row")
        family = context.complete_family_candidates
        if not isinstance(family, tuple) or not family:
            raise ValueError("complete_family_candidates must be a non-empty tuple")
        family_ids = tuple(id(candidate) for candidate in family)
        if len(set(family_ids)) != len(family_ids):
            raise ValueError("complete_family_candidates must not contain duplicates")
        if sum(candidate is context.accepted_row[1] for candidate in family) != 1:
            raise ValueError(
                "a complete family must include its accepted candidate exactly once"
            )
        indexed[row_id] = context

    missing = tuple(
        location for location in accepted_locations if id(location.row) not in indexed
    )
    if missing:
        raise ValueError(
            "{} accepted row(s) are missing complete family context".format(
                len(missing)
            )
        )
    if len(indexed) != len(accepted_locations):
        raise RuntimeError("accepted family context cardinality changed")
    return indexed


def _demoted_row(row, decision, location):
    metadata = dict(row[3]) if len(row) == 4 else {}
    prior_basis = metadata.get("classification_basis")
    if prior_basis is not None:
        prior_basis = _json_copy(prior_basis, "prior classification_basis")
    tandem_basis = _json_copy(decision.as_dict(), "accepted-family tandem basis")
    metadata["classification_basis"] = {
        "method": "source_backed_accepted_family_tandem_demotion",
        "prior_classification_basis": prior_basis,
        "accepted_family_tandem": tandem_basis,
        "transition": {
            "rule_name": RULE_NAME,
            "rule_version": RULE_VERSION,
            "coordinate_system": COORDINATE_SYSTEM,
            "source_category_index": location.category_index,
            "source_category_name": _CATEGORY_NAMES[location.category_index],
            "target_category_index": _TARGET_CATEGORY_INDEX,
            "target_category_name": _CATEGORY_NAMES[_TARGET_CATEGORY_INDEX],
            "locus_key": json.loads(location.locus_key_json),
        },
    }
    if isinstance(row, tuple):
        return tuple(row[:3]) + (metadata,)
    return list(row[:3]) + [metadata]


def _count_rows(categories, category_indexes=None):
    indexes = range(len(categories)) if category_indexes is None else category_indexes
    return sum(
        len(rows)
        for category_index in indexes
        for rows in categories[category_index].values()
    )


def apply_accepted_family_tandem_demotions(
    source, categories, accepted_row_families, source_cds
):
    """Apply independently validated accepted-family tandem demotions.

    Returns ``(new_categories, records, stats)``. ``new_categories`` has the
    same five-category shape as the input, but every category dictionary and
    every per-locus row list is newly allocated. ``records`` and ``stats`` are
    frozen dataclasses whose ``as_dict`` methods are JSON-safe.

    No partial result is exposed: all category, family, metadata, source, CDS,
    and reconstruction checks finish before the copied categories are edited.
    """

    _validate_source_context(source, source_cds)
    accepted_locations, total_rows_before = _index_categories(categories)
    contexts = _index_complete_family_contexts(
        accepted_row_families, accepted_locations
    )

    evaluated = []
    for location in accepted_locations:
        context = contexts[id(location.row)]
        decision = validate_accepted_family_tandem(
            source,
            location.row[1],
            context.complete_family_candidates,
            source_cds,
        )
        evaluated.append((location, decision))

    new_categories = [
        {locus_key: list(rows) for locus_key, rows in category.items()}
        for category in categories
    ]
    records = []
    supported_count = 0
    by_source_category = {0: 0, 2: 0}
    for location, decision in evaluated:
        if decision.supported:
            source_rows = new_categories[location.category_index][location.locus_key]
            matches = [
                index for index, row in enumerate(source_rows) if row is location.row
            ]
            if len(matches) != 1:
                raise RuntimeError("accepted row identity changed during transaction")
            source_rows.pop(matches[0])
            if not source_rows:
                del new_categories[location.category_index][location.locus_key]
            demoted_row = _demoted_row(location.row, decision, location)
            new_categories[_TARGET_CATEGORY_INDEX].setdefault(
                location.locus_key, []
            ).append(demoted_row)
            action = "demoted"
            target_category_index = _TARGET_CATEGORY_INDEX
            target_category_name = _CATEGORY_NAMES[_TARGET_CATEGORY_INDEX]
            supported_count += 1
            by_source_category[location.category_index] += 1
        else:
            action = "retained"
            target_category_index = None
            target_category_name = None

        records.append(
            AcceptedFamilyTandemDemotionRecord(
                action=action,
                source_category_index=location.category_index,
                source_category_name=_CATEGORY_NAMES[location.category_index],
                target_category_index=target_category_index,
                target_category_name=target_category_name,
                locus_key_json=location.locus_key_json,
                row_ordinal=location.row_ordinal,
                score=float(location.row[0]),
                repeat_count=len(location.row[1].list_repeats),
                accepted_interval=decision.accepted_interval,
                decision=decision,
            )
        )

    total_rows_after = _count_rows(new_categories)
    accepted_before = len(accepted_locations)
    accepted_after = _count_rows(new_categories, _ACCEPTED_CATEGORY_INDEXES)
    discarded_before = _count_rows(categories, (_TARGET_CATEGORY_INDEX,))
    discarded_after = _count_rows(new_categories, (_TARGET_CATEGORY_INDEX,))
    if total_rows_after != total_rows_before:
        raise RuntimeError("tandem demotion changed total row cardinality")
    if accepted_after != accepted_before - supported_count:
        raise RuntimeError("tandem demotion changed accepted cardinality incorrectly")
    if discarded_after != discarded_before + supported_count:
        raise RuntimeError("tandem demotion changed discarded cardinality incorrectly")

    stats = AcceptedFamilyTandemDemotionStats(
        total_rows_before=total_rows_before,
        total_rows_after=total_rows_after,
        accepted_rows_reviewed=accepted_before,
        tandem_supported=supported_count,
        tandem_not_supported=accepted_before - supported_count,
        bona_fide_demoted=by_source_category[0],
        possible_demoted=by_source_category[2],
        accepted_rows_before=accepted_before,
        accepted_rows_after=accepted_after,
        possible_discarded_rows_before=discarded_before,
        possible_discarded_rows_after=discarded_after,
    )
    return new_categories, tuple(records), stats


__all__ = (
    "AcceptedFamilyTandemDemotionRecord",
    "AcceptedFamilyTandemDemotionStats",
    "AcceptedRowFamilyContext",
    "apply_accepted_family_tandem_demotions",
)
