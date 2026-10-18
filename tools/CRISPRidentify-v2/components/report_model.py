"""Canonical, source-validated data model for CRISPRidentify reports.

Detector indexes are treated as claims, never as authoritative sequence data.
Every repeat, spacer, and complete array in the returned model is independently
re-sliced from the original input FASTA after reconstruction validation.

Public intervals are one-based, inclusive, and always refer to the input
forward strand.  A complete spacer deletion has no interval; it is represented
at an explicitly labelled zero-based interbase boundary.  Predicted biological
orientation is represented separately.
"""

import hashlib
import json
import math
import numbers
import os

from components.initial_competition import competition_evidence
from components.candidate_origins import candidate_origin_evidence
from components.reconstruction_boundary_family import (
    METHOD as BOUNDARY_FAMILY_METHOD, boundary_family_summary,
)
from components.possible_reconstruction import possible_reconstruction_summary
from components.contextual_spacer_rescue import (
    contextual_spacer_summary, contextual_spacer_member_intervals,
)
from components.boundary_assessment import (
    AlternativeRelationshipInput,
    BoundaryAssessment,
    BoundaryAssessmentNormalizer,
    BoundaryDimensionInput,
    CandidateGeometryInput,
    DIMENSION_OCCURRENCE_PATH,
    DIMENSION_REPEAT_PHASE,
    STATUS_AMBIGUOUS,
    geometry_input_from_candidate,
)
from components.cas_output import (
    CasReportError,
    normalize_cas_data,
)
from components.report_provenance_validation import (
    ReportProvenanceError,
    ReportProvenanceValidationContext,
    validate_classification_basis,
    validate_source_intact_recovery_alternatives,
)
from components.reconstruction_validator import (
    read_single_fasta,
    validate_candidate_reconstruction,
)
from components.likely_tandem_transaction import (
    CLASSIFICATION_BASIS_METHOD as LIKELY_TANDEM_BASIS_METHOD,
    TANDEM_TAG,
    validated_row_tags,
)


SCHEMA_NAME = "CRISPRidentify-report"
SCHEMA_VERSION = "1.2.0"

INPUT_CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)

INPUT_CATEGORY_KEYS = (
    "bona-fide",
    "alternative",
    "possible",
    "possible-discarded",
    "low-score",
)

REPORT_CATEGORY_INDEXES = (0, 2, 4)
CATEGORY_NAMES = tuple(
    INPUT_CATEGORY_NAMES[index] for index in REPORT_CATEGORY_INDEXES
)
CATEGORY_KEYS = tuple(
    INPUT_CATEGORY_KEYS[index] for index in REPORT_CATEGORY_INDEXES
)
ARRAY_ID_PREFIXES = {
    "Bona-fide": "bf",
    "Possible": "p",
    "Low score": "ls",
}

CLASSIFICATION_BASIS_METHODS = frozenset((
    LIKELY_TANDEM_BASIS_METHOD,
    "model_score_band",
    "source_backed_accepted_family_tandem_demotion",
    "source_backed_accepted_core_boundary_completion",
    "source_backed_architecture_rescue",
    "source_backed_intact_existence_rescue",
    "source_backed_intact_terminal_conflict_retirement",
    "source_backed_regular_subpath_recovery",
    "source_backed_accepted_family_rescue",
    "source_backed_terminal_repair_completion",
    "source_backed_internal_repeat_completion",
    "source_periodic_lattice_repair",
    "source_backed_candidate_reconstruction",
    "source_backed_possible_reconstruction",
    "source_backed_localized_spacer_context",
    BOUNDARY_FAMILY_METHOD,
))

_BOUNDARY_ASSESSMENT_ABSENT = object()


class ReportModelError(ValueError):
    """The evaluated results cannot be represented without ambiguity."""


class ReportModelValidationError(ReportModelError):
    """A candidate failed reconstruction against the original sequence."""

    def __init__(self, category, category_position, result):
        self.category = category
        self.category_position = category_position
        self.result = result
        issue_codes = ", ".join(issue.code for issue in result.issues)
        super().__init__(
            "cannot report invalid candidate: category={!r}, "
            "category_position={}, issues={}".format(
                category,
                category_position,
                issue_codes,
            )
        )


def _json_safe(value):
    """Return a deterministic JSON-compatible representation."""

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ReportModelError("report values must be finite")
        return value

    item = getattr(value, "item", None)
    if callable(item):
        try:
            scalar = item()
        except (TypeError, ValueError):
            pass
        else:
            if scalar is not value:
                return _json_safe(scalar)

    tolist = getattr(value, "tolist", None)
    if callable(tolist):
        converted = tolist()
        if converted is not value:
            return _json_safe(converted)

    if isinstance(value, dict):
        return {
            str(key): _json_safe(value[key])
            for key in sorted(value, key=lambda item: str(item))
        }
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]

    raise ReportModelError(
        "unsupported report value type: {}".format(type(value).__name__)
    )


def _canonical_json(value):
    """Serialize JSON-safe values so type-distinct claims stay distinct."""

    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def _score_value(value):
    if value is None:
        return None
    item = getattr(value, "item", None)
    if callable(item):
        try:
            value = item()
        except (TypeError, ValueError):
            pass
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ReportModelError("certainty score must be a finite number or None")
    value = float(value)
    if not math.isfinite(value):
        raise ReportModelError("certainty score must be finite")
    return value


def _flatten_feature_vector(value):
    value = _json_safe(value)
    while (
        isinstance(value, list)
        and len(value) == 1
        and isinstance(value[0], list)
    ):
        value = value[0]
    if not isinstance(value, list):
        raise ReportModelError("a model feature vector must be a sequence")
    if any(isinstance(item, (list, dict)) for item in value):
        raise ReportModelError("a model feature vector must be one-dimensional")
    return value


def _feature_sets(record, feature_name_sets):
    if not feature_name_sets:
        return []
    if len(record) < 3:
        raise ReportModelError(
            "candidate record is missing its evaluated feature vectors"
        )

    raw_sets = record[2]
    if isinstance(raw_sets, (str, bytes)) or raw_sets is None:
        raise ReportModelError("candidate feature vectors must be a sequence")
    try:
        raw_sets = list(raw_sets)
    except TypeError:
        raise ReportModelError("candidate feature vectors must be iterable")

    if len(raw_sets) != len(feature_name_sets):
        raise ReportModelError(
            "feature-name sets and evaluated feature vectors have different counts"
        )

    feature_sets = []
    for set_index, (names, raw_vector) in enumerate(
        zip(feature_name_sets, raw_sets),
        1,
    ):
        if isinstance(names, (str, bytes)):
            raise ReportModelError("feature-name sets must be sequences")
        names = list(names)
        values = _flatten_feature_vector(raw_vector)
        if len(names) != len(values):
            raise ReportModelError(
                "feature names and values have different counts in model {}".format(
                    set_index
                )
            )
        feature_sets.append({
            "model_index": set_index,
            "features": [
                {"name": str(name), "value": value}
                for name, value in zip(names, values)
            ],
        })
    return feature_sets


def _report_metadata(record):
    """Return the optional row metadata mapping after structural validation."""

    if len(record) < 4:
        return None
    metadata = record[3]
    if not isinstance(metadata, dict):
        raise ReportModelError(
            "candidate report metadata must be a mapping"
        )
    return metadata


def _classification_basis(
    record,
    source_sequence,
    expected_source_sha256,
    candidate,
    validation_context,
    category_index,
):
    """Return validated optional category-assignment provenance.

    Evaluated rows historically contain only ``score``, ``candidate``, and
    optionally the model feature vectors.  A fourth item may now carry report
    metadata, but only its explicitly named ``classification_basis`` mapping
    enters the canonical report.  The score remains sourced exclusively from
    ``record[0]``.
    """

    metadata = _report_metadata(record)
    if metadata is None:
        return None
    if "classification_basis" not in metadata:
        # Boundary assessment is an independent report concern.  A row may
        # carry only this assessment without inventing classification
        # provenance that the pipeline did not supply.
        if (
            "boundary_assessment" in metadata
            or "source_intact_recovery_alternatives" in metadata
            or "initial_competition" in metadata
        ):
            return None
        raise ReportModelError(
            "candidate report metadata is missing classification_basis"
        )
    raw_basis = metadata["classification_basis"]
    if not isinstance(raw_basis, dict):
        raise ReportModelError("classification_basis must be a mapping")

    method = raw_basis.get("method")
    if (
        not isinstance(method, str)
        or method not in CLASSIFICATION_BASIS_METHODS
    ):
        raise ReportModelError(
            "classification_basis method must be one of: {}".format(
                ", ".join(sorted(CLASSIFICATION_BASIS_METHODS))
            )
        )
    try:
        validated = validate_classification_basis(
            source_sequence,
            expected_source_sha256,
            candidate,
            raw_basis,
            row_score=record[0],
            allowed_methods=CLASSIFICATION_BASIS_METHODS,
            validation_context=validation_context,
            category_index=category_index,
            boundary_assessment=metadata.get("boundary_assessment"),
            current_row=record,
        )
    except ReportProvenanceError as error:
        raise ReportModelError(
            "invalid classification_basis metadata: {}".format(error)
        ) from error
    return _json_safe(validated)


def _classification_tags(record, classification_basis, category_index):
    """Validate optional row tags against their classification provenance."""

    try:
        tags = validated_row_tags(record)
    except ValueError as error:
        raise ReportModelError("invalid candidate tags: {}".format(error)) from error
    tandem_basis = (
        isinstance(classification_basis, dict)
        and classification_basis.get("method") == LIKELY_TANDEM_BASIS_METHOD
    )
    if tandem_basis != (TANDEM_TAG in tags):
        raise ReportModelError(
            "the tandem tag and tandem-demotion basis must occur together"
        )
    if tandem_basis and category_index != 4:
        raise ReportModelError(
            "the tandem tag is only valid on Low score arrays"
        )
    return tags


def _required_mapping(mapping, key, label):
    if not isinstance(mapping, dict):
        raise TypeError("{} must be a mapping".format(label))
    if key not in mapping:
        raise ValueError("{} is missing {!r}".format(label, key))
    value = mapping[key]
    if not isinstance(value, dict):
        raise TypeError("{}.{} must be a mapping".format(label, key))
    return value


def _required_sequence(mapping, key, label):
    if not isinstance(mapping, dict):
        raise TypeError("{} must be a mapping".format(label))
    if key not in mapping:
        raise ValueError("{} is missing {!r}".format(label, key))
    value = mapping[key]
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise TypeError("{}.{} must be a sequence".format(label, key))
    return value


def _serialized_internal_interval(mapping, label):
    internal = _required_mapping(
        mapping,
        "internal_0based_half_open",
        label,
    )
    if "start" not in internal or "end" not in internal:
        raise ValueError(
            "{}.internal_0based_half_open must contain start and end".format(
                label
            )
        )
    return internal["start"], internal["end"]


def _geometry_input_from_boundary_mapping(mapping, label):
    """Rehydrate only explicit geometry claims from a serialized assessment."""

    if not isinstance(mapping, dict):
        raise TypeError("{} must be a mapping".format(label))
    repeats = _required_sequence(mapping, "repeats", label)
    spacers = _required_sequence(mapping, "spacers", label)
    return CandidateGeometryInput(
        source_sha256=mapping.get("source_sha256"),
        source_interval=_serialized_internal_interval(
            _required_mapping(mapping, "source_interval", label),
            "{}.source_interval".format(label),
        ),
        repeat_intervals=tuple(
            _serialized_internal_interval(
                _required_mapping(component, "interval", component_label),
                "{}.interval".format(component_label),
            )
            for index, component in enumerate(repeats)
            for component_label in (
                "{}.repeats[{}]".format(label, index),
            )
        ),
        repeat_sequences=tuple(
            component.get("sequence")
            if isinstance(component, dict) else None
            for component in repeats
        ),
        spacer_intervals=tuple(
            _serialized_internal_interval(
                _required_mapping(component, "interval", component_label),
                "{}.interval".format(component_label),
            )
            for index, component in enumerate(spacers)
            for component_label in (
                "{}.spacers[{}]".format(label, index),
            )
        ),
        spacer_sequences=tuple(
            component.get("sequence")
            if isinstance(component, dict) else None
            for component in spacers
        ),
    )


def _dimension_input_from_boundary_mapping(mapping, dimension):
    dimensions = _required_mapping(mapping, "dimensions", "boundary_assessment")
    value = _required_mapping(
        dimensions,
        dimension,
        "boundary_assessment.dimensions",
    )
    return BoundaryDimensionInput(
        status=value.get("status"),
        selection_policy_id=value.get("selection_policy_id"),
        reason_codes=tuple(
            _required_sequence(
                value,
                "reason_codes",
                "boundary_assessment.dimensions.{}".format(dimension),
            )
        ),
    )


def _boundary_inputs_from_canonical_mapping(mapping):
    """Convert an untrusted canonical mapping back to normalization inputs."""

    reported = _required_mapping(
        mapping,
        "reported_geometry",
        "boundary_assessment",
    )
    relationships = _required_sequence(
        mapping,
        "alternative_relationships",
        "boundary_assessment",
    )
    alternative_inputs = []
    for index, relationship in enumerate(relationships):
        label = "boundary_assessment.alternative_relationships[{}]".format(
            index
        )
        if not isinstance(relationship, dict):
            raise TypeError("{} must be a mapping".format(label))
        alternative_inputs.append(AlternativeRelationshipInput(
            geometry=_geometry_input_from_boundary_mapping(
                _required_mapping(relationship, "geometry", label),
                "{}.geometry".format(label),
            ),
            label=relationship.get("label"),
            dimensions=tuple(
                _required_sequence(relationship, "dimensions", label)
            ),
            reason_codes=tuple(
                _required_sequence(relationship, "reason_codes", label)
            ),
        ))
    return {
        "reported_geometry": _geometry_input_from_boundary_mapping(
            reported,
            "boundary_assessment.reported_geometry",
        ),
        "alternative_relationships": tuple(alternative_inputs),
        "occurrence_path": _dimension_input_from_boundary_mapping(
            mapping,
            DIMENSION_OCCURRENCE_PATH,
        ),
        "repeat_phase": _dimension_input_from_boundary_mapping(
            mapping,
            DIMENSION_REPEAT_PHASE,
        ),
    }


def revalidate_boundary_assessment_for_report(
    source_sequence,
    expected_source_sha256,
    candidate,
    assessment=_BOUNDARY_ASSESSMENT_ABSENT,
    source_normalizer=None,
    additional_alternative_relationships=(),
    additional_reason_codes=("source_intact_recovery_alternative_revalidated",),
):
    """Return a source-reconstructed canonical assessment for one report row.

    ``assessment`` may be either a :class:`BoundaryAssessment` or the exact
    mapping returned by :meth:`BoundaryAssessment.as_dict`.  It is treated as
    untrusted in both cases: every representative and alternative component
    is sliced again from ``source_sequence``, the whole serialized payload is
    compared with the freshly normalized result, and the representative must
    equal the current row candidate.  Independently revalidated transaction
    alternatives may be supplied separately; genuinely new geometries or
    dimensions downgrade the affected dimension to ``ambiguous`` instead of
    leaving a stale ``unassessed`` or ``policy_selected`` claim.
    """

    try:
        if source_normalizer is None:
            source_normalizer = BoundaryAssessmentNormalizer(
                source_sequence,
                expected_source_sha256,
            )
        elif not isinstance(
            source_normalizer,
            BoundaryAssessmentNormalizer,
        ):
            raise TypeError(
                "source_normalizer must be a BoundaryAssessmentNormalizer"
            )
        else:
            source_normalizer.assert_bound_to(
                source_sequence,
                expected_source_sha256,
            )

        candidate_geometry = geometry_input_from_candidate(
            candidate, expected_source_sha256
        )
        candidate_assessment = source_normalizer.normalize(candidate_geometry)
        additional_relationships = tuple(
            additional_alternative_relationships
        )
        if assessment is _BOUNDARY_ASSESSMENT_ABSENT:
            if not additional_relationships:
                return candidate_assessment.as_dict()
            inputs = {
                "reported_geometry": candidate_geometry,
                "alternative_relationships": (),
                "occurrence_path": BoundaryDimensionInput(),
                "repeat_phase": BoundaryDimensionInput(),
            }
        else:
            if isinstance(assessment, BoundaryAssessment):
                claimed_mapping = assessment.as_dict()
            elif isinstance(assessment, dict):
                claimed_mapping = _json_safe(assessment)
            else:
                raise TypeError(
                    "boundary_assessment must be a BoundaryAssessment or its "
                    "canonical mapping"
                )

            source_claim = _required_mapping(
                claimed_mapping,
                "source",
                "boundary_assessment",
            )
            if source_claim.get("sha256") != expected_source_sha256:
                raise ValueError(
                    "boundary_assessment source SHA-256 is stale"
                )
            if source_claim.get("length") != source_normalizer.source_length:
                raise ValueError(
                    "boundary_assessment source length is stale"
                )

            inputs = _boundary_inputs_from_canonical_mapping(claimed_mapping)
            normalized = source_normalizer.normalize(**inputs)
            normalized_mapping = normalized.as_dict()
            if _canonical_json(claimed_mapping) != _canonical_json(
                normalized_mapping
            ):
                raise ValueError(
                    "boundary_assessment is not the exact canonical "
                    "serialization reconstructed from original source"
                )
            if (
                normalized.reported_geometry.as_dict()
                != candidate_assessment.reported_geometry.as_dict()
            ):
                raise ValueError(
                    "boundary_assessment representative geometry does not "
                    "match the current candidate geometry"
                )
            if not additional_relationships:
                return normalized_mapping

        relationships = list(inputs["alternative_relationships"])
        relationship_indexes = {
            tuple(relationship.geometry.repeat_intervals): index
            for index, relationship in enumerate(relationships)
        }
        affected_dimensions = set()
        dimension_order = (
            DIMENSION_OCCURRENCE_PATH,
            DIMENSION_REPEAT_PHASE,
        )
        for relationship in additional_relationships:
            if not isinstance(relationship, AlternativeRelationshipInput):
                raise TypeError(
                    "additional alternatives must be "
                    "AlternativeRelationshipInput objects"
                )
            geometry_key = tuple(relationship.geometry.repeat_intervals)
            existing_index = relationship_indexes.get(geometry_key)
            if existing_index is None:
                relationship_indexes[geometry_key] = len(relationships)
                relationships.append(relationship)
                affected_dimensions.update(relationship.dimensions)
                continue
            existing = relationships[existing_index]
            added_dimensions = set(relationship.dimensions).difference(
                existing.dimensions
            )
            if not added_dimensions:
                continue
            combined_dimensions = tuple(
                dimension for dimension in dimension_order
                if (
                    dimension in existing.dimensions
                    or dimension in relationship.dimensions
                )
            )
            relationships[existing_index] = AlternativeRelationshipInput(
                geometry=existing.geometry,
                label=existing.label,
                dimensions=combined_dimensions,
                reason_codes=tuple(sorted(set(
                    existing.reason_codes + relationship.reason_codes
                ))),
            )
            affected_dimensions.update(added_dimensions)

        def reconciled_dimension(dimension, current):
            if dimension not in affected_dimensions:
                return current
            return BoundaryDimensionInput(
                status=STATUS_AMBIGUOUS,
                reason_codes=tuple(sorted(set(
                    current.reason_codes
                    + tuple(additional_reason_codes)
                ))),
            )

        return source_normalizer.normalize(
            inputs["reported_geometry"],
            alternative_relationships=tuple(relationships),
            occurrence_path=reconciled_dimension(
                DIMENSION_OCCURRENCE_PATH, inputs["occurrence_path"]
            ),
            repeat_phase=reconciled_dimension(
                DIMENSION_REPEAT_PHASE, inputs["repeat_phase"]
            ),
        ).as_dict()
    except ReportModelError:
        raise
    except (
        AttributeError,
        IndexError,
        KeyError,
        TypeError,
        ValueError,
    ) as error:
        raise ReportModelError(
            "invalid boundary_assessment metadata: {}".format(error)
        ) from error


def _indexed_value(mapping, index):
    """Defensively read an integer-positioned legacy annotation mapping."""

    if not isinstance(mapping, dict):
        return None
    if index in mapping:
        return mapping[index]
    string_index = str(index)
    if string_index in mapping:
        return mapping[string_index]
    return None


def _orientation(non_array_data, category_name, category_position, requested):
    if not requested:
        return {
            "status": "not_computed",
            "label": "Orientation was not computed",
            "strand": None,
        }

    strand_root = non_array_data.get("Strand", {})
    category_strands = (
        strand_root.get(category_name, {})
        if isinstance(strand_root, dict)
        else {}
    )
    label = _indexed_value(category_strands, category_position)

    if label == "Forward":
        return {
            "status": "computed",
            "label": "Forward",
            "strand": "+",
        }
    if label in ("Reverse", "Reversed"):
        return {
            "status": "computed",
            "label": str(label),
            "strand": "-",
        }
    if label in (None, "Forward (Orientation was not computed)"):
        return {
            "status": "not_computed",
            "label": "Orientation was not computed",
            "strand": None,
        }
    raise ReportModelError(
        "unknown orientation label for {} candidate {}: {!r}".format(
            category_name,
            category_position,
            label,
        )
    )


def _checked_source_slice(genome, start, end, claimed, label):
    """Validate a zero-based half-open claim, then return source DNA."""

    if (
        not isinstance(start, int)
        or isinstance(start, bool)
        or not isinstance(end, int)
        or isinstance(end, bool)
        or start < 0
        or end <= start
        or end > len(genome)
    ):
        raise ReportModelError(
            "{} has invalid source interval [{}, {})".format(label, start, end)
        )
    if not isinstance(claimed, str):
        raise ReportModelError("{} sequence must be a string".format(label))
    source_sequence = genome[start:end]
    if source_sequence != claimed.upper():
        raise ReportModelError(
            "{} does not match original source DNA at [{}, {})".format(
                label,
                start,
                end,
            )
        )
    return source_sequence


def _public_interval(start_zero, end_exclusive):
    return {
        "start": start_zero + 1,
        "end": end_exclusive,
        "length": end_exclusive - start_zero,
    }


def _candidate_identity(source_sha256, provisional):
    identity = {
        "source_sha256": source_sha256,
        "array": provisional["source_interval"],
        "repeat_intervals": [
            repeat["source_interval"] for repeat in provisional["repeats"]
        ],
        "spacer_intervals": [
            spacer["source_interval"] for spacer in provisional["spacers"]
        ],
        "spacer_deletion_boundaries": [
            event["source_boundary"]
            for event in provisional["spacer_deletions"]
        ],
    }
    encoded = json.dumps(
        identity,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:12]
    interval = provisional["source_interval"]
    return "CRISPR-{}-{}-{}".format(
        interval["start"],
        interval["end"],
        digest,
    )


class ReportModelBuilder:
    """Build one canonical report directly from evaluated candidates."""

    def __init__(
        self,
        file_path,
        categories,
        non_array_data,
        list_feature_names,
        parameters=None,
        flags=None,
        provenance_categories=None,
    ):
        self.file_path = os.fspath(file_path)
        self.categories = categories
        self.non_array_data = non_array_data or {}
        self.list_feature_names = list_feature_names or []
        self.parameters = parameters or {}
        self.flags = flags or {}
        self.provenance_categories = provenance_categories

    def _validate_container(self):
        if isinstance(self.categories, (str, bytes)) or self.categories is None:
            raise ReportModelError("categories must contain five category mappings")
        try:
            category_count = len(self.categories)
        except TypeError:
            raise ReportModelError("categories must contain five category mappings")
        if category_count != len(INPUT_CATEGORY_NAMES):
            raise ReportModelError(
                "expected {} categories, observed {}".format(
                    len(INPUT_CATEGORY_NAMES),
                    category_count,
                )
            )
        for category_name, category in zip(INPUT_CATEGORY_NAMES, self.categories):
            if not isinstance(category, dict):
                raise ReportModelError(
                    "{} category must be an interval mapping".format(category_name)
                )

    def _build_candidate(
        self,
        genome,
        source_sha256,
        source_normalizer,
        provenance_validation_context,
        category_index,
        category_position,
        interval_key,
        record,
    ):
        category_name = INPUT_CATEGORY_NAMES[category_index]
        category_key = INPUT_CATEGORY_KEYS[category_index]
        if not isinstance(record, (list, tuple)) or len(record) < 2:
            raise ReportModelError(
                "malformed {} candidate record at position {}".format(
                    category_name,
                    category_position,
                )
            )

        score, candidate = record[:2]
        validation = validate_candidate_reconstruction(genome, candidate)
        if not validation.valid:
            raise ReportModelValidationError(
                category_name,
                category_position,
                validation,
            )

        repeat_starts = validation.genomic_repeat_starts
        repeats = []
        spacers = []
        spacer_deletions = []
        architecture = []

        for index, (start, claimed_repeat) in enumerate(
            zip(repeat_starts, candidate.list_repeats),
            1,
        ):
            end = start + len(claimed_repeat)
            sequence = _checked_source_slice(
                genome,
                start,
                end,
                claimed_repeat,
                "repeat {}".format(index),
            )
            repeat_record = {
                "index": index,
                "source_interval": _public_interval(start, end),
                "sequence": sequence,
            }
            repeats.append(repeat_record)
            architecture.append({
                "type": "repeat",
                "index": index,
                "source_interval": dict(repeat_record["source_interval"]),
                "sequence": sequence,
            })

            spacer_offset = index - 1
            if spacer_offset < len(candidate.list_spacers):
                spacer_start = end
                spacer_end = repeat_starts[index]
                claimed_spacer = candidate.list_spacers[spacer_offset]
                if claimed_spacer == "":
                    if spacer_start != spacer_end:
                        raise ReportModelError(
                            "empty spacer {} does not have a zero-width "
                            "source gap [{}, {})".format(
                                index,
                                spacer_start,
                                spacer_end,
                            )
                        )
                    source_boundary = {
                        "position": spacer_start,
                        "coordinate_system": "0-based interbase",
                        "reference": "input_forward_strand",
                    }
                    deletion_event = {
                        "index": index,
                        "source_boundary": source_boundary,
                        "flanking_repeat_indices": [index, index + 1],
                        "length": 0,
                        "sequence": "",
                    }
                    spacer_deletions.append(deletion_event)
                    architecture.append({
                        "type": "spacer-deletion",
                        "index": index,
                        "source_boundary": dict(source_boundary),
                        "flanking_repeat_indices": [index, index + 1],
                        "length": 0,
                        "sequence": "",
                    })
                else:
                    spacer_sequence = _checked_source_slice(
                        genome,
                        spacer_start,
                        spacer_end,
                        claimed_spacer,
                        "spacer {}".format(index),
                    )
                    spacer_record = {
                        "index": index,
                        "source_interval": _public_interval(
                            spacer_start,
                            spacer_end,
                        ),
                        "sequence": spacer_sequence,
                    }
                    spacers.append(spacer_record)
                    architecture.append({
                        "type": "spacer",
                        "index": index,
                        "source_interval": dict(
                            spacer_record["source_interval"]
                        ),
                        "sequence": spacer_sequence,
                    })

        array_start = validation.array_start
        array_end = validation.array_end
        claimed_array = "".join(
            element["sequence"] for element in architecture
        )
        array_sequence = _checked_source_slice(
            genome,
            array_start,
            array_end,
            claimed_array,
            "complete array",
        )

        consensus = candidate.consensus
        if not isinstance(consensus, str):
            raise ReportModelError("candidate consensus must be a string")

        metadata = _report_metadata(record)
        classification_basis = _classification_basis(
            record,
            genome,
            source_sha256,
            candidate,
            provenance_validation_context,
            category_index,
        )
        classification_tags = _classification_tags(
            record,
            classification_basis,
            category_index,
        )
        source_intact_alternatives = None
        derived_boundary_relationships = ()
        family_summary = boundary_family_summary(record)
        if (
            metadata is not None
            and "source_intact_recovery_alternatives" in metadata
        ):
            try:
                source_intact_alternatives = (
                    validate_source_intact_recovery_alternatives(
                        genome,
                        source_sha256,
                        candidate,
                        metadata["source_intact_recovery_alternatives"],
                        validation_context=(
                            provenance_validation_context
                        ),
                        category_index=category_index,
                    )
                )
            except ReportProvenanceError as error:
                raise ReportModelError(
                    "invalid source_intact_recovery_alternatives metadata: "
                    "{}".format(error)
                ) from error
            derived_boundary_relationships = (
                source_intact_alternatives.boundary_relationships
            )
        if family_summary is not None:
            family_relationships = []
            for member in classification_basis["boundary_family"]["members"]:
                if member["location"]["geometry_id"] == family_summary["display_geometry_id"]:
                    continue
                intervals = tuple(map(tuple, member["fresh"]["repeat_intervals"]))
                gaps = tuple((a[1], b[0]) for a, b in zip(intervals, intervals[1:]))
                geometry = CandidateGeometryInput(
                    source_sha256=source_sha256,
                    source_interval=(intervals[0][0], intervals[-1][1]),
                    repeat_intervals=intervals,
                    repeat_sequences=tuple(genome[a:b] for a, b in intervals),
                    spacer_intervals=gaps,
                    spacer_sequences=tuple(genome[a:b] for a, b in gaps),
                )
                family_relationships.append(AlternativeRelationshipInput(
                    geometry=geometry, label="coupled_boundary_family_member",
                    dimensions=(DIMENSION_REPEAT_PHASE,),
                    reason_codes=("shared_ordered_starts_uncertain_repeat_ends",),
                ))
            derived_boundary_relationships += tuple(family_relationships)
        displayed_intervals = tuple(
            (start, start + len(repeat)) for start, repeat in
            zip(candidate.list_repeat_starts, candidate.list_repeats)
        )
        contextual_relationships = []
        for member_intervals in contextual_spacer_member_intervals(record):
            intervals = tuple(map(tuple, member_intervals))
            if intervals == displayed_intervals:
                continue
            gaps = tuple((a[1], b[0]) for a, b in zip(intervals, intervals[1:]))
            contextual_relationships.append(AlternativeRelationshipInput(
                geometry=CandidateGeometryInput(
                    source_sha256=source_sha256,
                    source_interval=(intervals[0][0], intervals[-1][1]),
                    repeat_intervals=intervals,
                    repeat_sequences=tuple(genome[a:b] for a, b in intervals),
                    spacer_intervals=gaps,
                    spacer_sequences=tuple(genome[a:b] for a, b in gaps),
                ),
                label="coupled_contextual_spacer_member",
                dimensions=(DIMENSION_REPEAT_PHASE,),
                reason_codes=("shared_ordered_starts_uncertain_repeat_ends",),
            ))
        derived_boundary_relationships += tuple(contextual_relationships)
        boundary_assessment = revalidate_boundary_assessment_for_report(
            genome,
            source_sha256,
            candidate,
            (
                metadata["boundary_assessment"]
                if (
                    metadata is not None
                    and "boundary_assessment" in metadata
                )
                else _BOUNDARY_ASSESSMENT_ABSENT
            ),
            source_normalizer=source_normalizer,
            additional_alternative_relationships=(
                derived_boundary_relationships
            ),
            additional_reason_codes=(("shared_ordered_starts_uncertain_repeat_ends",)
                if family_summary is not None or contextual_relationships
                else ("source_intact_recovery_alternative_revalidated",)),
        )
        provisional = {
            "id": None,
            "category": category_name,
            "category_key": category_key,
            "category_position": category_position,
            "legacy_interval_key": _json_safe(interval_key),
            "certainty_score": _score_value(score),
            "source_interval": _public_interval(array_start, array_end),
            "source_sequence": array_sequence,
            "orientation": _orientation(
                self.non_array_data,
                category_name,
                category_position,
                bool(self.flags.get("flag_strand", False)),
            ),
            "consensus_repeat": consensus,
            "repeat_count": len(repeats),
            "spacer_count": len(candidate.list_spacers),
            "non_deleted_spacer_count": len(spacers),
            "spacer_deletion_count": len(spacer_deletions),
            "repeats": repeats,
            "spacers": spacers,
            "spacer_deletions": spacer_deletions,
            "architecture": architecture,
            "feature_sets": _feature_sets(record, self.list_feature_names),
            "boundary_assessment": boundary_assessment,
            "validation": {
                "status": "valid",
                "source_reconstructed": True,
                "method": "independent_original_sequence_slicing",
                "issues": [],
            },
            "_category_index": category_index,
        }
        if classification_basis is not None:
            provisional["classification_basis"] = classification_basis
        if family_summary is not None:
            provisional["boundary_family"] = family_summary
        possible_summary = possible_reconstruction_summary(record)
        if possible_summary is not None:
            provisional["possible_reconstruction"] = possible_summary
        contextual_summary = contextual_spacer_summary(record)
        if contextual_summary is not None:
            provisional["contextual_spacer_context"] = contextual_summary
        # An initial equal-score representative is an observed ranking choice,
        # not new classification authority or a uniquely resolved repeat phase.
        try:
            competition = competition_evidence(record, len(genome), source=genome)
            if self.parameters.get("param_reconstruction_policy") == "conservative":
                candidate_origin_evidence(candidate, source_sha256)
        except ValueError as error:
            raise ReportModelError(
                "invalid initial_competition metadata: {}".format(error)
            ) from error
        if competition is not None:
            provisional["initial_competition"] = _json_safe(competition)
        if classification_tags:
            provisional["tags"] = list(classification_tags)
        if source_intact_alternatives is not None:
            provisional["source_intact_recovery_alternatives"] = _json_safe(
                source_intact_alternatives.records
            )
        provisional["id"] = _candidate_identity(source_sha256, provisional)
        return provisional

    def build(self):
        self._validate_container()
        source_header, genome = read_single_fasta(self.file_path)
        source_normalizer = BoundaryAssessmentNormalizer.from_source(genome)
        source_sha256 = source_normalizer.source_sha256
        provenance_categories = (
            self.categories
            if self.provenance_categories is None
            else self.provenance_categories
        )
        provenance_validation_context = (
            ReportProvenanceValidationContext.from_source(
                genome,
                source_sha256,
                categories=provenance_categories,
                internal_completion_evidence={
                    "boundary_records": self.parameters.get("param_repeat_boundary_repair", {}).get("records", []),
                    "internal_records": self.parameters.get("param_internal_repeat_completion", {}).get("records", []),
                    "current_authority_inventory": self.parameters.get("param_internal_repeat_completion", {}).get("current_authority_inventory"),
                },
                accepted_only_projection=(
                    self.provenance_categories is None
                    and
                    self.parameters.get("param_accepted_only") is True
                ),
            )
        )
        try:
            provenance_validation_context.ordinary_occurrence_witnesses()
        except ValueError as error:
            raise ReportModelError("invalid initial_competition metadata: {}".format(error)) from error
        source_id = source_header.split()[0]

        provisional_arrays = []
        for category_index in REPORT_CATEGORY_INDEXES:
            category = self.categories[category_index]
            category_position = 0
            for interval_key in sorted(category):
                records = category[interval_key]
                if isinstance(records, (str, bytes)) or records is None:
                    raise ReportModelError(
                        "{} interval records must be a sequence".format(
                            INPUT_CATEGORY_NAMES[category_index]
                        )
                    )
                for record in records:
                    provisional_arrays.append(
                        self._build_candidate(
                            genome,
                            source_sha256,
                            source_normalizer,
                            provenance_validation_context,
                            category_index,
                            category_position,
                            interval_key,
                            record,
                        )
                    )
                    category_position += 1

        provisional_arrays.sort(key=lambda array: (
            array["source_interval"]["start"],
            array["source_interval"]["end"],
            array["_category_index"],
            array["category_position"],
            array["id"],
        ))

        # Public IDs retain source coordinates and use readable per-category
        # numbering instead of exposing an opaque content-hash suffix.
        category_numbers = {name: 0 for name in CATEGORY_NAMES}
        for array in provisional_arrays:
            category = array["category"]
            category_numbers[category] += 1
            interval = array["source_interval"]
            array["id"] = "CRISPR-{}-{}-{}{}".format(
                interval["start"],
                interval["end"],
                ARRAY_ID_PREFIXES[category],
                category_numbers[category],
            )
            del array["_category_index"]

        category_counts = {
            name: sum(
                1 for array in provisional_arrays
                if array["category"] == name
            )
            for name in CATEGORY_NAMES
        }
        array_count = len(provisional_arrays)
        try:
            cas_report = normalize_cas_data(
                self.non_array_data.get(
                    "Cas", {"status": "not_requested", "systems": []}
                ),
                source_id=source_id,
                source_length=len(genome),
            )
        except CasReportError as error:
            raise ReportModelError(
                "invalid Cas report data: {}".format(error)
            ) from error

        run_metadata = {
            "parameters": _json_safe(self.parameters),
            "analyses": {
                "orientation": (
                    "requested"
                    if self.flags.get("flag_strand", False)
                    else "not_requested"
                ),
                "is_element": (
                    "requested"
                    if self.flags.get("flag_is", False)
                    else "not_requested"
                ),
                "cas_prediction": (
                    "requested"
                    if self.flags.get("flag_cas", False)
                    else "not_requested"
                ),
            },
        }
        runtime_preset = self.parameters.get("param_runtime_preset")
        if runtime_preset is not None:
            run_metadata["configuration"] = _json_safe({
                "schema_version": 1,
                "release_version": self.parameters.get(
                    "param_release_version", "unknown"
                ),
                "preset": runtime_preset,
                "preset_provenance": self.parameters.get(
                    "param_runtime_preset_provenance", {}
                ),
                "effective_config": self.parameters.get(
                    "param_effective_config", {}
                ),
            })

        return {
            "schema": {
                "name": SCHEMA_NAME,
                "version": SCHEMA_VERSION,
            },
            "coordinate_system": {
                "indexing": 1,
                "end_inclusive": True,
                "reference": "input_forward_strand",
                "orientation_is_separate": True,
                "zero_width_boundaries": {
                    "indexing": 0,
                    "representation": "interbase",
                    "reference": "input_forward_strand",
                },
            },
            "source": {
                "id": source_id,
                "header": source_header,
                "file_name": os.path.basename(self.file_path),
                "length": len(genome),
                "sha256": source_sha256,
                "sha256_scope": "uppercase_sequence",
            },
            "run": run_metadata,
            "summary": {
                "array_count": array_count,
                "accepted_array_count": (
                    category_counts["Bona-fide"]
                    + category_counts["Possible"]
                ),
                "category_counts": category_counts,
            },
            "validation": {
                "status": "valid",
                "method": "independent_original_sequence_slicing",
                "validated_array_count": array_count,
                "invalid_array_count": 0,
            },
            "arrays": provisional_arrays,
            "cas": _json_safe(cas_report),
        }
