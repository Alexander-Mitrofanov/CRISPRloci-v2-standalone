"""Source-first boundary assessment records for CRISPR array output.

This module is intentionally independent of detection, classification, and
reporting.  It normalizes untrusted candidate geometry claims by slicing every
array, repeat, and spacer interval from the complete original input-forward
DNA sequence.

Boundary uncertainty has two independent dimensions:

``occurrence_path``
    Which repeat occurrences belong to the biological array.

``repeat_phase``
    Which source bases delimit each repeat occurrence.

Each dimension is ``unassessed``, ``policy_selected``, or ``ambiguous``.
``policy_selected`` records a reproducible choice among explicit alternatives;
it is not an exact-boundary claim.  The reported geometry therefore always has
the role ``representative``.

All coordinates accepted internally are zero-based and end-exclusive.
Serialized non-empty intervals additionally expose one-based inclusive public
coordinates.  An empty spacer is serialized as a zero-based interbase
boundary because it has no valid one-based inclusive interval.

Alternatives must be supplied as explicit, labeled relationships between
independently source-valid geometries.  This module deliberately has no fuzzy
family identifier input: sharing an upstream family identifier is not evidence
that two candidates are biological boundary alternatives.
"""

from dataclasses import dataclass
import hashlib
import json
import re
from typing import Optional, Tuple

from components.reconstruction_validator import IUPAC_DNA_ALPHABET


SCHEMA_NAME = "crispridentify.boundary_assessment"
SCHEMA_VERSION = "1.0.0"

DIMENSION_OCCURRENCE_PATH = "occurrence_path"
DIMENSION_REPEAT_PHASE = "repeat_phase"
BOUNDARY_DIMENSIONS = (
    DIMENSION_OCCURRENCE_PATH,
    DIMENSION_REPEAT_PHASE,
)

STATUS_UNASSESSED = "unassessed"
STATUS_POLICY_SELECTED = "policy_selected"
STATUS_AMBIGUOUS = "ambiguous"
BOUNDARY_STATUSES = (
    STATUS_UNASSESSED,
    STATUS_POLICY_SELECTED,
    STATUS_AMBIGUOUS,
)

ROLE_REPRESENTATIVE = "representative"
ROLE_ALTERNATIVE = "alternative"

INTERNAL_COORDINATE_SYSTEM = "zero_based_half_open"
PUBLIC_COORDINATE_SYSTEM = "one_based_inclusive"
ZERO_WIDTH_COORDINATE_SYSTEM = "zero_based_interbase"
SOURCE_REFERENCE = "original_input_forward_strand"
SOURCE_HASH_SCOPE = "uppercase_sequence"

_TOKEN_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:/-]*$")
_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_UPPER_IUPAC = frozenset(character.upper() for character in IUPAC_DNA_ALPHABET)


def _json_sha256(payload):
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def _sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _normalize_source(source):
    if not isinstance(source, str):
        raise TypeError("source_sequence must be a string")
    if not source:
        raise ValueError("source_sequence must be non-empty")
    try:
        source.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError("source_sequence must contain ASCII DNA symbols") from error
    normalized = source.upper()
    invalid = sorted(set(normalized).difference(_UPPER_IUPAC))
    if invalid:
        raise ValueError(
            "source_sequence contains non-IUPAC DNA symbols: {}".format(
                ", ".join(repr(symbol) for symbol in invalid)
            )
        )
    return normalized


def _normalize_sha256(value, label):
    if not isinstance(value, str) or _SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError("{} must be a lowercase hexadecimal SHA-256".format(label))
    return value


def source_sha256(source_sequence):
    """Return the canonical hash used by boundary-assessment records."""

    return _sequence_sha256(_normalize_source(source_sequence))


def _as_tuple(value, label):
    if value is None or isinstance(value, (str, bytes)):
        raise TypeError("{} must be an iterable, not a string".format(label))
    try:
        return tuple(value)
    except TypeError as error:
        raise TypeError("{} must be an iterable".format(label)) from error


def _normalize_interval(value, label, source_length, allow_empty):
    coordinates = _as_tuple(value, label)
    if len(coordinates) != 2:
        raise ValueError("{} must contain exactly start and end".format(label))
    start, end = coordinates
    if (
        not isinstance(start, int)
        or isinstance(start, bool)
        or not isinstance(end, int)
        or isinstance(end, bool)
    ):
        raise ValueError("{} coordinates must be integers".format(label))
    if start < 0 or end < 0:
        raise ValueError("{} coordinates must be non-negative".format(label))
    if end < start:
        raise ValueError("{} has reversed coordinates".format(label))
    if not allow_empty and end == start:
        raise ValueError("{} must be non-empty".format(label))
    if end > source_length:
        raise ValueError("{} is out of source bounds".format(label))
    return start, end


def _normalize_sequence(value, label, allow_empty):
    if not isinstance(value, str):
        raise TypeError("{} must be a string".format(label))
    normalized = value.upper()
    if not allow_empty and not normalized:
        raise ValueError("{} must be non-empty".format(label))
    invalid = sorted(set(normalized).difference(_UPPER_IUPAC))
    if invalid:
        raise ValueError(
            "{} contains non-IUPAC DNA symbols: {}".format(
                label,
                ", ".join(repr(symbol) for symbol in invalid),
            )
        )
    return normalized


def _normalize_token(value, label):
    if (
        not isinstance(value, str)
        or _TOKEN_PATTERN.fullmatch(value) is None
    ):
        raise ValueError(
            "{} must be a non-empty stable token using letters, numbers, "
            "'.', '_', ':', '/', or '-'".format(label)
        )
    return value


def _normalize_reason_codes(values, label):
    values = _as_tuple(values, label)
    normalized = tuple(
        _normalize_token(value, "{}[{}]".format(label, index))
        for index, value in enumerate(values)
    )
    if len(normalized) != len(set(normalized)):
        raise ValueError("{} contains duplicate reason codes".format(label))
    return tuple(sorted(normalized))


def _interval_as_dict(interval):
    start, end = interval
    serialized = {
        "internal_0based_half_open": {
            "start": start,
            "end": end,
        },
        "length": end - start,
    }
    if start == end:
        serialized["public_1based_inclusive"] = None
        serialized["zero_width_boundary"] = {
            "position": start,
            "coordinate_system": ZERO_WIDTH_COORDINATE_SYSTEM,
            "reference": SOURCE_REFERENCE,
        }
    else:
        serialized["public_1based_inclusive"] = {
            "start": start + 1,
            "end": end,
        }
    return serialized


def coordinate_system_as_dict():
    """Return the complete internal/public coordinate serialization contract."""

    return {
        "internal": {
            "name": INTERNAL_COORDINATE_SYSTEM,
            "indexing": 0,
            "end_exclusive": True,
            "reference": SOURCE_REFERENCE,
        },
        "public_nonempty_intervals": {
            "name": PUBLIC_COORDINATE_SYSTEM,
            "indexing": 1,
            "end_inclusive": True,
            "reference": SOURCE_REFERENCE,
            "conversion": (
                "public_start=internal_start+1; "
                "public_end=internal_end"
            ),
        },
        "empty_spacers": {
            "name": ZERO_WIDTH_COORDINATE_SYSTEM,
            "indexing": 0,
            "representation": "interbase_boundary",
            "reference": SOURCE_REFERENCE,
            "one_based_inclusive_interval": None,
        },
    }


@dataclass(frozen=True)
class CandidateGeometryInput:
    """An untrusted explicit geometry claim to validate against the source."""

    source_sha256: str
    source_interval: tuple
    repeat_intervals: tuple
    repeat_sequences: tuple
    spacer_intervals: tuple
    spacer_sequences: tuple


@dataclass(frozen=True)
class AlternativeRelationshipInput:
    """An explicitly labeled biological relationship to one alternative."""

    geometry: CandidateGeometryInput
    label: str
    dimensions: tuple
    reason_codes: tuple = ()


@dataclass(frozen=True)
class BoundaryDimensionInput:
    """Caller-supplied assessment for one boundary dimension."""

    status: str = STATUS_UNASSESSED
    selection_policy_id: Optional[str] = None
    reason_codes: tuple = ()


@dataclass(frozen=True)
class ValidatedSourceSlice:
    """One immutable sequence slice independently checked against the source."""

    ordinal: int
    interval: Tuple[int, int]
    sequence: str
    sha256: str

    def as_dict(self):
        return {
            "ordinal": self.ordinal,
            "interval": _interval_as_dict(self.interval),
            "sequence": self.sequence,
            "sha256": self.sha256,
            "sha256_scope": "uppercase_source_slice",
        }


@dataclass(frozen=True)
class ValidatedBoundaryGeometry:
    """One immutable, source-reconstructed repeat/spacer geometry."""

    geometry_id: str
    role: str
    source_sha256: str
    source_interval: Tuple[int, int]
    source_slice_sha256: str
    repeats: Tuple[ValidatedSourceSlice, ...]
    spacers: Tuple[ValidatedSourceSlice, ...]

    @property
    def repeat_intervals(self):
        return tuple(component.interval for component in self.repeats)

    @property
    def spacer_intervals(self):
        return tuple(component.interval for component in self.spacers)

    def as_dict(self):
        return {
            "geometry_id": self.geometry_id,
            "role": self.role,
            "exact_boundary_claim": False,
            "source_sha256": self.source_sha256,
            "source_interval": _interval_as_dict(self.source_interval),
            "source_slice_sha256": self.source_slice_sha256,
            "source_slice_sha256_scope": "uppercase_source_slice",
            "repeat_count": len(self.repeats),
            "spacer_count": len(self.spacers),
            "repeats": [component.as_dict() for component in self.repeats],
            "spacers": [component.as_dict() for component in self.spacers],
        }


@dataclass(frozen=True)
class BoundaryAlternativeRelationship:
    """One explicit relationship between representative and alternative."""

    relationship_id: str
    label: str
    dimensions: Tuple[str, ...]
    reason_codes: Tuple[str, ...]
    geometry: ValidatedBoundaryGeometry

    def as_dict(self):
        return {
            "relationship_id": self.relationship_id,
            "label": self.label,
            "dimensions": list(self.dimensions),
            "reason_codes": list(self.reason_codes),
            "geometry": self.geometry.as_dict(),
        }


@dataclass(frozen=True)
class BoundaryDimensionAssessment:
    """Normalized status and explicitly related alternatives for one dimension."""

    dimension: str
    status: str
    representative_geometry_id: str
    alternative_geometry_ids: Tuple[str, ...]
    selection_policy_id: Optional[str]
    reason_codes: Tuple[str, ...]

    def as_dict(self):
        return {
            "dimension": self.dimension,
            "status": self.status,
            "reported_geometry_role": ROLE_REPRESENTATIVE,
            "representative_geometry_id": self.representative_geometry_id,
            "alternative_geometry_ids": list(self.alternative_geometry_ids),
            "selection_policy_id": self.selection_policy_id,
            "reason_codes": list(self.reason_codes),
        }


@dataclass(frozen=True)
class BoundaryAssessmentIntegrationContract:
    """Machine-readable requirements for future pipeline/report integration."""

    source_requirement: str
    source_hash_scope: str
    internal_coordinates: str
    public_coordinates: str
    reported_geometry_semantics: str
    alternative_relationship_semantics: str
    fuzzy_family_semantics: str
    mutation_semantics: str
    status_semantics: Tuple[str, ...]

    def as_dict(self):
        return {
            "source_requirement": self.source_requirement,
            "source_hash_scope": self.source_hash_scope,
            "internal_coordinates": self.internal_coordinates,
            "public_coordinates": self.public_coordinates,
            "reported_geometry_semantics": self.reported_geometry_semantics,
            "alternative_relationship_semantics": (
                self.alternative_relationship_semantics
            ),
            "fuzzy_family_semantics": self.fuzzy_family_semantics,
            "mutation_semantics": self.mutation_semantics,
            "status_semantics": list(self.status_semantics),
        }


BOUNDARY_ASSESSMENT_INTEGRATION_CONTRACT = (
    BoundaryAssessmentIntegrationContract(
        source_requirement=(
            "complete_original_input_forward_strand_sequence_is_required"
        ),
        source_hash_scope=SOURCE_HASH_SCOPE,
        internal_coordinates=INTERNAL_COORDINATE_SYSTEM,
        public_coordinates=(
            "one_based_inclusive_for_nonempty_intervals_and_"
            "zero_based_interbase_for_empty_spacers"
        ),
        reported_geometry_semantics=(
            "representative_never_implicitly_exact"
        ),
        alternative_relationship_semantics=(
            "explicit_labeled_source_valid_geometries_only"
        ),
        fuzzy_family_semantics=(
            "shared_fuzzy_family_id_is_not_alternative_evidence"
        ),
        mutation_semantics=(
            "pure_normalization_no_candidate_category_or_pipeline_mutation"
        ),
        status_semantics=(
            "unassessed_means_no_explicit_alternative_was_adjudicated",
            "policy_selected_means_policy_chose_the_representative_not_exact",
            "ambiguous_means_one_or_more_explicit_alternatives_remain",
        ),
    )
)


@dataclass(frozen=True)
class BoundaryAssessment:
    """Complete immutable and JSON-safe boundary assessment."""

    assessment_id: str
    source_sha256: str
    source_length: int
    reported_geometry: ValidatedBoundaryGeometry
    alternative_relationships: Tuple[BoundaryAlternativeRelationship, ...]
    occurrence_path: BoundaryDimensionAssessment
    repeat_phase: BoundaryDimensionAssessment
    integration_contract: BoundaryAssessmentIntegrationContract

    def as_dict(self):
        return {
            "schema": {
                "name": SCHEMA_NAME,
                "version": SCHEMA_VERSION,
            },
            "assessment_id": self.assessment_id,
            "source": {
                "length": self.source_length,
                "sha256": self.source_sha256,
                "sha256_scope": SOURCE_HASH_SCOPE,
                "reference": SOURCE_REFERENCE,
            },
            "coordinate_system": coordinate_system_as_dict(),
            "reported_geometry_role": ROLE_REPRESENTATIVE,
            "reported_geometry": self.reported_geometry.as_dict(),
            "dimensions": {
                DIMENSION_OCCURRENCE_PATH: self.occurrence_path.as_dict(),
                DIMENSION_REPEAT_PHASE: self.repeat_phase.as_dict(),
            },
            "alternative_relationships": [
                relationship.as_dict()
                for relationship in self.alternative_relationships
            ],
            "integration_contract": self.integration_contract.as_dict(),
        }


def geometry_input_from_candidate(candidate, expected_source_sha256):
    """Copy a candidate's absolute geometry into an immutable untrusted input.

    The returned object contains no runtime identity.  Validation still occurs
    only when :class:`BoundaryAssessmentNormalizer` receives the original
    source sequence.
    """

    expected_source_sha256 = _normalize_sha256(
        expected_source_sha256,
        "expected_source_sha256",
    )
    required = (
        "list_repeat_starts",
        "list_repeats",
        "list_spacers",
    )
    for field_name in required:
        if not hasattr(candidate, field_name):
            raise TypeError(
                "candidate is missing {!r}".format(field_name)
            )

    starts = _as_tuple(candidate.list_repeat_starts, "list_repeat_starts")
    repeats = _as_tuple(candidate.list_repeats, "list_repeats")
    spacers = _as_tuple(candidate.list_spacers, "list_spacers")
    if not repeats or len(starts) != len(repeats):
        raise ValueError(
            "candidate repeat starts and repeat sequences must be non-empty "
            "and have equal cardinality"
        )
    if any(not isinstance(repeat, str) for repeat in repeats):
        raise TypeError("candidate repeat sequences must be strings")

    repeat_intervals = tuple(
        (start, start + len(repeat))
        if isinstance(start, int) and not isinstance(start, bool)
        else (start, start)
        for start, repeat in zip(starts, repeats)
    )
    spacer_intervals = tuple(
        (left_end, right_start)
        for (_left_start, left_end), (right_start, _right_end)
        in zip(repeat_intervals, repeat_intervals[1:])
    )
    return CandidateGeometryInput(
        source_sha256=expected_source_sha256,
        source_interval=(
            repeat_intervals[0][0],
            repeat_intervals[-1][1],
        ),
        repeat_intervals=repeat_intervals,
        repeat_sequences=tuple(repeats),
        spacer_intervals=spacer_intervals,
        spacer_sequences=tuple(spacers),
    )


def _geometry_id(source_hash, repeat_intervals):
    return "boundary-geometry-sha256:" + _json_sha256(
        {
            "source_sha256": source_hash,
            "repeat_intervals_0based_half_open": [
                list(interval) for interval in repeat_intervals
            ],
        }
    )


def _normalize_geometry(source, source_hash, geometry, role, label):
    if not isinstance(geometry, CandidateGeometryInput):
        raise TypeError("{} must be a CandidateGeometryInput".format(label))
    claim_hash = _normalize_sha256(
        geometry.source_sha256,
        "{}.source_sha256".format(label),
    )
    if claim_hash != source_hash:
        raise ValueError("{} carries a stale source SHA-256".format(label))

    repeat_interval_values = _as_tuple(
        geometry.repeat_intervals,
        "{}.repeat_intervals".format(label),
    )
    repeat_sequence_values = _as_tuple(
        geometry.repeat_sequences,
        "{}.repeat_sequences".format(label),
    )
    spacer_interval_values = _as_tuple(
        geometry.spacer_intervals,
        "{}.spacer_intervals".format(label),
    )
    spacer_sequence_values = _as_tuple(
        geometry.spacer_sequences,
        "{}.spacer_sequences".format(label),
    )
    if not repeat_interval_values:
        raise ValueError("{} must contain at least one repeat".format(label))
    if len(repeat_interval_values) != len(repeat_sequence_values):
        raise ValueError(
            "{} repeat interval/sequence cardinality differs".format(label)
        )
    expected_spacer_count = len(repeat_interval_values) - 1
    if (
        len(spacer_interval_values) != expected_spacer_count
        or len(spacer_sequence_values) != expected_spacer_count
    ):
        raise ValueError(
            "{} must contain exactly one spacer interval and sequence "
            "between adjacent repeats".format(label)
        )

    repeat_intervals = tuple(
        _normalize_interval(
            value,
            "{}.repeat_intervals[{}]".format(label, index),
            len(source),
            allow_empty=False,
        )
        for index, value in enumerate(repeat_interval_values)
    )
    spacer_intervals = tuple(
        _normalize_interval(
            value,
            "{}.spacer_intervals[{}]".format(label, index),
            len(source),
            allow_empty=True,
        )
        for index, value in enumerate(spacer_interval_values)
    )
    source_interval = _normalize_interval(
        geometry.source_interval,
        "{}.source_interval".format(label),
        len(source),
        allow_empty=False,
    )

    for index, (left, right) in enumerate(
        zip(repeat_intervals, repeat_intervals[1:])
    ):
        if right[0] < left[1]:
            raise ValueError(
                "{} repeat intervals overlap at adjacency {}".format(
                    label,
                    index,
                )
            )
        expected = left[1], right[0]
        if spacer_intervals[index] != expected:
            raise ValueError(
                "{} spacer interval {} does not equal the exact source gap "
                "{}".format(label, index, expected)
            )

    expected_source_interval = (
        repeat_intervals[0][0],
        repeat_intervals[-1][1],
    )
    if source_interval != expected_source_interval:
        raise ValueError(
            "{} source interval does not equal its repeat envelope {}".format(
                label,
                expected_source_interval,
            )
        )

    repeat_slices = []
    for index, (interval, claimed_sequence) in enumerate(
        zip(repeat_intervals, repeat_sequence_values)
    ):
        claimed = _normalize_sequence(
            claimed_sequence,
            "{}.repeat_sequences[{}]".format(label, index),
            allow_empty=False,
        )
        observed = source[interval[0]:interval[1]]
        if claimed != observed:
            raise ValueError(
                "{} repeat {} does not match its original source slice "
                "{}".format(label, index, interval)
            )
        repeat_slices.append(
            ValidatedSourceSlice(
                ordinal=index + 1,
                interval=interval,
                sequence=observed,
                sha256=_sequence_sha256(observed),
            )
        )

    spacer_slices = []
    for index, (interval, claimed_sequence) in enumerate(
        zip(spacer_intervals, spacer_sequence_values)
    ):
        claimed = _normalize_sequence(
            claimed_sequence,
            "{}.spacer_sequences[{}]".format(label, index),
            allow_empty=True,
        )
        observed = source[interval[0]:interval[1]]
        if claimed != observed:
            raise ValueError(
                "{} spacer {} does not match its original source slice "
                "{}".format(label, index, interval)
            )
        spacer_slices.append(
            ValidatedSourceSlice(
                ordinal=index + 1,
                interval=interval,
                sequence=observed,
                sha256=_sequence_sha256(observed),
            )
        )

    represented = []
    for index, repeat_slice in enumerate(repeat_slices):
        represented.append(repeat_slice.sequence)
        if index < len(spacer_slices):
            represented.append(spacer_slices[index].sequence)
    observed_source_slice = source[source_interval[0]:source_interval[1]]
    if "".join(represented) != observed_source_slice:
        raise ValueError(
            "{} components do not reconstruct the complete source interval".format(
                label
            )
        )

    return ValidatedBoundaryGeometry(
        geometry_id=_geometry_id(source_hash, repeat_intervals),
        role=role,
        source_sha256=source_hash,
        source_interval=source_interval,
        source_slice_sha256=_sequence_sha256(observed_source_slice),
        repeats=tuple(repeat_slices),
        spacers=tuple(spacer_slices),
    )


def _normalize_dimensions(values, label):
    values = _as_tuple(values, "{}.dimensions".format(label))
    if not values:
        raise ValueError("{} must label at least one dimension".format(label))
    if len(values) != len(set(values)):
        raise ValueError("{} contains duplicate dimensions".format(label))
    invalid = tuple(value for value in values if value not in BOUNDARY_DIMENSIONS)
    if invalid:
        raise ValueError(
            "{} contains unsupported dimensions: {}".format(
                label,
                ", ".join(repr(value) for value in invalid),
            )
        )
    return tuple(
        dimension for dimension in BOUNDARY_DIMENSIONS
        if dimension in values
    )


def _relationship_id(
    representative_geometry_id,
    alternative_geometry_id,
    relationship_label,
    dimensions,
):
    return "boundary-relationship-sha256:" + _json_sha256(
        {
            "representative_geometry_id": representative_geometry_id,
            "alternative_geometry_id": alternative_geometry_id,
            "label": relationship_label,
            "dimensions": list(dimensions),
        }
    )


def _normalize_relationships(
    source,
    source_hash,
    representative,
    alternative_relationships,
):
    inputs = _as_tuple(
        alternative_relationships,
        "alternative_relationships",
    )
    normalized = []
    seen_geometry_ids = set()
    for index, relationship in enumerate(inputs):
        label = "alternative_relationships[{}]".format(index)
        if not isinstance(relationship, AlternativeRelationshipInput):
            raise TypeError(
                "{} must be an AlternativeRelationshipInput".format(label)
            )
        relationship_label = _normalize_token(
            relationship.label,
            "{}.label".format(label),
        )
        dimensions = _normalize_dimensions(relationship.dimensions, label)
        reason_codes = _normalize_reason_codes(
            relationship.reason_codes,
            "{}.reason_codes".format(label),
        )
        geometry = _normalize_geometry(
            source,
            source_hash,
            relationship.geometry,
            ROLE_ALTERNATIVE,
            "{}.geometry".format(label),
        )
        if geometry.geometry_id == representative.geometry_id:
            raise ValueError(
                "{} repeats the representative geometry".format(label)
            )
        if geometry.geometry_id in seen_geometry_ids:
            raise ValueError(
                "{} repeats an already supplied alternative geometry".format(
                    label
                )
            )
        seen_geometry_ids.add(geometry.geometry_id)
        normalized.append(
            BoundaryAlternativeRelationship(
                relationship_id=_relationship_id(
                    representative.geometry_id,
                    geometry.geometry_id,
                    relationship_label,
                    dimensions,
                ),
                label=relationship_label,
                dimensions=dimensions,
                reason_codes=reason_codes,
                geometry=geometry,
            )
        )
    return tuple(
        sorted(
            normalized,
            key=lambda relationship: (
                relationship.geometry.geometry_id,
                relationship.dimensions,
                relationship.label,
            ),
        )
    )


def _normalize_dimension(
    dimension,
    value,
    representative,
    relationships,
):
    if not isinstance(value, BoundaryDimensionInput):
        raise TypeError(
            "{} must be a BoundaryDimensionInput".format(dimension)
        )
    if value.status not in BOUNDARY_STATUSES:
        raise ValueError(
            "{} has unsupported status {!r}".format(dimension, value.status)
        )
    policy_id = (
        None
        if value.selection_policy_id is None
        else _normalize_token(
            value.selection_policy_id,
            "{}.selection_policy_id".format(dimension),
        )
    )
    reason_codes = _normalize_reason_codes(
        value.reason_codes,
        "{}.reason_codes".format(dimension),
    )
    alternative_geometry_ids = tuple(
        relationship.geometry.geometry_id
        for relationship in relationships
        if dimension in relationship.dimensions
    )

    if value.status == STATUS_UNASSESSED:
        if alternative_geometry_ids:
            raise ValueError(
                "{} cannot be unassessed when explicit alternatives are "
                "labeled for that dimension".format(dimension)
            )
        if policy_id is not None:
            raise ValueError(
                "{} unassessed status cannot name a selection policy".format(
                    dimension
                )
            )
    elif value.status == STATUS_AMBIGUOUS:
        if not alternative_geometry_ids:
            raise ValueError(
                "{} ambiguous status requires an explicit alternative".format(
                    dimension
                )
            )
        if policy_id is not None:
            raise ValueError(
                "{} ambiguous status cannot claim policy selection".format(
                    dimension
                )
            )
    elif value.status == STATUS_POLICY_SELECTED:
        if not alternative_geometry_ids:
            raise ValueError(
                "{} policy-selected status requires an explicit alternative".format(
                    dimension
                )
            )
        if policy_id is None:
            raise ValueError(
                "{} policy-selected status requires selection_policy_id".format(
                    dimension
                )
            )
        if not reason_codes:
            raise ValueError(
                "{} policy-selected status requires reason_codes".format(
                    dimension
                )
            )

    return BoundaryDimensionAssessment(
        dimension=dimension,
        status=value.status,
        representative_geometry_id=representative.geometry_id,
        alternative_geometry_ids=alternative_geometry_ids,
        selection_policy_id=policy_id,
        reason_codes=reason_codes,
    )


def _assessment_id(
    source_hash,
    representative,
    relationships,
    occurrence_path,
    repeat_phase,
):
    return "boundary-assessment-sha256:" + _json_sha256(
        {
            "source_sha256": source_hash,
            "representative_geometry_id": representative.geometry_id,
            "relationship_ids": [
                relationship.relationship_id
                for relationship in relationships
            ],
            "dimensions": {
                occurrence_path.dimension: occurrence_path.as_dict(),
                repeat_phase.dimension: repeat_phase.as_dict(),
            },
        }
    )


class BoundaryAssessmentNormalizer:
    """Normalize many geometries against one immutable source sequence.

    Construction performs the only whole-source normalization and SHA-256
    calculation.  :meth:`normalize` then validates only the representative
    and explicitly supplied alternative intervals.  Keeping this context
    explicit avoids both process-global cache state and O(candidate_count *
    source_length) report generation.
    """

    __slots__ = ("_source", "_source_hash", "_source_input")

    def __init__(self, source_sequence, expected_source_sha256):
        self._bind_source(source_sequence)
        expected_hash = _normalize_sha256(
            expected_source_sha256,
            "expected_source_sha256",
        )
        if self._source_hash != expected_hash:
            raise ValueError(
                "expected_source_sha256 is stale for the supplied source "
                "sequence"
            )

    def _bind_source(self, source_sequence):
        source = _normalize_source(source_sequence)
        actual_hash = _sequence_sha256(source)
        self._source = source
        self._source_hash = actual_hash
        # Retain the exact caller object so internal consumers can establish
        # an O(1), collision-independent binding without rescanning the DNA.
        self._source_input = source_sequence

    @classmethod
    def from_source(cls, source_sequence):
        """Create a context and derive its canonical source hash once."""

        instance = cls.__new__(cls)
        instance._bind_source(source_sequence)
        return instance

    @property
    def source_sha256(self):
        return self._source_hash

    @property
    def source_length(self):
        return len(self._source)

    def assert_bound_to(self, source_sequence, expected_source_sha256):
        """Fail unless inputs name the exact source object for this context."""

        expected_hash = _normalize_sha256(
            expected_source_sha256,
            "expected_source_sha256",
        )
        if expected_hash != self._source_hash:
            raise ValueError(
                "boundary normalizer carries a stale source SHA-256"
            )
        if (
            source_sequence is not self._source_input
            and source_sequence is not self._source
        ):
            raise ValueError(
                "boundary normalizer is not bound to this source sequence "
                "object"
            )

    def normalize(
        self,
        reported_geometry,
        alternative_relationships=(),
        occurrence_path=BoundaryDimensionInput(),
        repeat_phase=BoundaryDimensionInput(),
    ):
        """Validate one complete assessment without rescanning the source."""

        representative = _normalize_geometry(
            self._source,
            self._source_hash,
            reported_geometry,
            ROLE_REPRESENTATIVE,
            "reported_geometry",
        )
        relationships = _normalize_relationships(
            self._source,
            self._source_hash,
            representative,
            alternative_relationships,
        )
        normalized_occurrence_path = _normalize_dimension(
            DIMENSION_OCCURRENCE_PATH,
            occurrence_path,
            representative,
            relationships,
        )
        normalized_repeat_phase = _normalize_dimension(
            DIMENSION_REPEAT_PHASE,
            repeat_phase,
            representative,
            relationships,
        )

        return BoundaryAssessment(
            assessment_id=_assessment_id(
                self._source_hash,
                representative,
                relationships,
                normalized_occurrence_path,
                normalized_repeat_phase,
            ),
            source_sha256=self._source_hash,
            source_length=len(self._source),
            reported_geometry=representative,
            alternative_relationships=relationships,
            occurrence_path=normalized_occurrence_path,
            repeat_phase=normalized_repeat_phase,
            integration_contract=BOUNDARY_ASSESSMENT_INTEGRATION_CONTRACT,
        )




__all__ = [
    "AlternativeRelationshipInput",
    "BOUNDARY_ASSESSMENT_INTEGRATION_CONTRACT",
    "BOUNDARY_DIMENSIONS",
    "BOUNDARY_STATUSES",
    "BoundaryAssessment",
    "BoundaryAssessmentIntegrationContract",
    "BoundaryAssessmentNormalizer",
    "BoundaryAlternativeRelationship",
    "BoundaryDimensionAssessment",
    "BoundaryDimensionInput",
    "CandidateGeometryInput",
    "DIMENSION_OCCURRENCE_PATH",
    "DIMENSION_REPEAT_PHASE",
    "INTERNAL_COORDINATE_SYSTEM",
    "PUBLIC_COORDINATE_SYSTEM",
    "ROLE_ALTERNATIVE",
    "ROLE_REPRESENTATIVE",
    "SCHEMA_NAME",
    "SCHEMA_VERSION",
    "SOURCE_HASH_SCOPE",
    "STATUS_AMBIGUOUS",
    "STATUS_POLICY_SELECTED",
    "STATUS_UNASSESSED",
    "ValidatedBoundaryGeometry",
    "ValidatedSourceSlice",
    "ZERO_WIDTH_COORDINATE_SYSTEM",
    "coordinate_system_as_dict",
    "geometry_input_from_candidate",
    "source_sha256",
]
