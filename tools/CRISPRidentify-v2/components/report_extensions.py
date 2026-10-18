"""Pure, backward-compatible extensions for canonical CRISPRidentify reports.

The canonical ``CRISPRidentify-report`` 1.2.0 payload remains the authoritative
record of detector output.  Functions in this module return new dictionaries:
they never mutate the input report, its arrays, classifications, scores, or
summary.  Additional interpretation is attached only below ``extensions``.

The geometry extension deliberately uses only canonical, source-reconstructed
repeat intervals and the matching architecture.  It does not use another
detector, benchmark labels, or comparison outcomes.
"""

import copy
import hashlib
import json
import math
import numbers
import re
import statistics

from components.biological_loci import build_biological_loci_extension
from components.circular_array_join_integration import (
    CircularAcceptedJoinResult,
)
from components.likely_tandem import (
    CATEGORY_NAME as LIKELY_TANDEM_CATEGORY,
    LIKELY_TANDEM_POLICY,
    RULE_NAME as LIKELY_TANDEM_RULE_NAME,
    RULE_VERSION as LIKELY_TANDEM_RULE_VERSION,
    evaluate_likely_tandem,
)
from components.likely_tandem_transaction import (
    CLASSIFICATION_BASIS_METHOD,
    LikelyTandemQuarantineRecord,
    LikelyTandemQuarantineStats,
    TANDEM_TAG,
    TRANSACTION_RULE_NAME,
    TRANSACTION_RULE_VERSION,
    apply_likely_tandem_quarantine,
)
from components.report_model import CATEGORY_NAMES
from components.report_model import SCHEMA_NAME
from components.report_model import SCHEMA_VERSION
from components.sequence_topology import SequenceTopologyEvidence
from components.competing_terminal_extents import build_competing_terminal_extents


ACCEPTED_CATEGORIES = frozenset(("Bona-fide", "Possible"))

CLASSIFICATION_SCHEMA_NAME = (
    "CRISPRidentify-classification-semantics"
)
CLASSIFICATION_SCHEMA_VERSION = "1.3.0"
SOURCE_GEOMETRY_SCHEMA_NAME = "CRISPRidentify-source-geometry"
SOURCE_GEOMETRY_SCHEMA_VERSION = "1.0.0"
SEQUENCE_TOPOLOGY_SCHEMA_NAME = (
    "CRISPRidentify-sequence-topology-evidence"
)
CIRCULAR_ORIGIN_LOCI_SCHEMA_NAME = (
    "CRISPRidentify-circular-origin-loci"
)
CIRCULAR_ORIGIN_LOCI_SCHEMA_VERSION = "1.0.0"
LIKELY_TANDEM_SCHEMA_NAME = "CRISPRidentify-likely-tandem-quarantine"
LIKELY_TANDEM_SCHEMA_VERSION = "1.1.0"
SPACER_GAP_CV_ADVISORY_THRESHOLD = 0.75

_SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class ReportExtensionError(ValueError):
    """A canonical report cannot safely support the requested extensions."""


def _is_plain_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _strict_json_equal(left, right):
    """Compare JSON values without Python's bool/int numeric coercion."""

    try:
        options = {
            "allow_nan": False,
            "ensure_ascii": True,
            "separators": (",", ":"),
            "sort_keys": True,
        }
        return json.dumps(left, **options) == json.dumps(right, **options)
    except (TypeError, ValueError):
        return False


def _validated_array_tags(array, label):
    tags = array.get("tags", [])
    if (
        not isinstance(tags, list)
        or any(
            not isinstance(tag, str)
            or re.fullmatch(r"[a-z][a-z0-9_-]{0,63}", tag) is None
            for tag in tags
        )
        or tags != sorted(set(tags))
    ):
        raise ReportExtensionError(
            "{} tags must be unique sorted lowercase identifiers".format(
                label
            )
        )
    return tuple(tags)


def _source_metadata(report):
    source = report.get("source")
    if not isinstance(source, dict):
        raise ReportExtensionError("report source must be a mapping")

    source_sha256 = source.get("sha256")
    if (
        not isinstance(source_sha256, str)
        or _SHA256_PATTERN.fullmatch(source_sha256) is None
    ):
        raise ReportExtensionError(
            "source.sha256 must be a lowercase hexadecimal SHA-256 digest"
        )

    source_length = source.get("length")
    if not _is_plain_int(source_length) or source_length < 0:
        raise ReportExtensionError(
            "source.length must be a non-negative integer"
        )
    return source, source_sha256, source_length


def _validate_report_envelope(report):
    if not isinstance(report, dict):
        raise ReportExtensionError("report must be a mapping")

    schema = report.get("schema")
    if (
        not isinstance(schema, dict)
        or schema.get("name") != SCHEMA_NAME
        or schema.get("version") != SCHEMA_VERSION
    ):
        raise ReportExtensionError(
            "extensions require {} {}".format(SCHEMA_NAME, SCHEMA_VERSION)
        )

    validation = report.get("validation")
    if (
        not isinstance(validation, dict)
        or validation.get("status") != "valid"
        or validation.get("method")
        != "independent_original_sequence_slicing"
    ):
        raise ReportExtensionError(
            "report must have valid independent source reconstruction"
        )

    arrays = report.get("arrays")
    if not isinstance(arrays, list):
        raise ReportExtensionError("report arrays must be a list")

    observed_counts = {name: 0 for name in CATEGORY_NAMES}
    array_ids = []
    for position, array in enumerate(arrays):
        if not isinstance(array, dict):
            raise ReportExtensionError(
                "array {} must be a mapping".format(position)
            )
        array_id = array.get("id")
        if not isinstance(array_id, str) or not array_id:
            raise ReportExtensionError(
                "array {} must have a non-empty id".format(position)
            )
        category = array.get("category")
        if category not in observed_counts:
            raise ReportExtensionError(
                "array {!r} has an unknown category".format(array_id)
            )
        array_ids.append(array_id)
        _validated_array_tags(array, "array {!r}".format(array_id))
        observed_counts[category] += 1
    if len(array_ids) != len(set(array_ids)):
        raise ReportExtensionError("report array ids must be unique")

    summary = report.get("summary")
    expected_accepted = (
        observed_counts["Bona-fide"] + observed_counts["Possible"]
    )
    if (
        not isinstance(summary, dict)
        or summary.get("array_count") != len(arrays)
        or summary.get("accepted_array_count") != expected_accepted
        or summary.get("category_counts") != observed_counts
    ):
        raise ReportExtensionError(
            "report summary counts disagree with arrays"
        )
    if (
        validation.get("validated_array_count") != len(arrays)
        or validation.get("invalid_array_count") != 0
    ):
        raise ReportExtensionError(
            "report validation counts disagree with arrays"
        )

    source, source_sha256, source_length = _source_metadata(report)
    return arrays, source, source_sha256, source_length


def _interval(value, label, source_length, containing=None):
    if not isinstance(value, dict):
        raise ReportExtensionError("{} must be a mapping".format(label))

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
        raise ReportExtensionError(
            "{} must be a valid one-based inclusive source interval".format(
                label
            )
        )

    interval = {"start": start, "end": end, "length": length}
    if containing is not None and (
        start < containing["start"] or end > containing["end"]
    ):
        raise ReportExtensionError(
            "{} must be contained by its array interval".format(label)
        )
    return interval


def _finite_parameter(parameters, name):
    value = parameters.get(name)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ReportExtensionError(
            "{} must be a finite number or null".format(name)
        )
    if not math.isfinite(float(value)):
        raise ReportExtensionError(
            "{} must be a finite number or null".format(name)
        )
    # Keep configured integer values as integers in the public representation.
    return value if _is_plain_int(value) else float(value)


def _configured_spacer_limits(report):
    run = report.get("run")
    if not isinstance(run, dict):
        raise ReportExtensionError("report run must be a mapping")
    parameters = run.get("parameters")
    if not isinstance(parameters, dict):
        raise ReportExtensionError("run.parameters must be a mapping")

    minimum_name = "param_min_avg_spacer_length"
    maximum_name = "param_max_avg_spacer_length"
    minimum = _finite_parameter(parameters, minimum_name)
    maximum = _finite_parameter(parameters, maximum_name)
    if (
        minimum is not None
        and maximum is not None
        and float(minimum) > float(maximum)
    ):
        raise ReportExtensionError(
            "configured minimum mean spacer length exceeds the maximum"
        )

    return {
        "minimum": minimum,
        "maximum": maximum,
        "unit": "bp",
        "minimum_parameter": minimum_name,
        "maximum_parameter": maximum_name,
        "role": (
            "configuration context only; these values and the diagnostics "
            "below do not reclassify an array"
        ),
    }


def _descriptive_statistics(values):
    """Return deterministic population statistics for a numeric sequence."""

    values = list(values)
    if not values:
        return {
            "count": 0,
            "mean": None,
            "median": None,
            "median_absolute_deviation": None,
            "minimum": None,
            "maximum": None,
            "range": None,
            "population_standard_deviation": None,
            "coefficient_of_variation": None,
        }

    mean = float(statistics.fmean(values))
    median = float(statistics.median(values))
    mad = float(statistics.median(
        [abs(value - median) for value in values]
    ))
    population_sd = float(statistics.pstdev(values))
    return {
        "count": len(values),
        "mean": mean,
        "median": median,
        "median_absolute_deviation": mad,
        "minimum": min(values),
        "maximum": max(values),
        "range": max(values) - min(values),
        "population_standard_deviation": population_sd,
        "coefficient_of_variation": (
            population_sd / mean if mean != 0.0 else None
        ),
    }


def _repeat_geometry(array, source_length):
    array_id = array.get("id")
    if not isinstance(array_id, str) or not array_id:
        raise ReportExtensionError(
            "accepted array must have a non-empty id"
        )

    validation = array.get("validation")
    if (
        not isinstance(validation, dict)
        or validation.get("status") != "valid"
        or validation.get("source_reconstructed") is not True
        or validation.get("method")
        != "independent_original_sequence_slicing"
    ):
        raise ReportExtensionError(
            "accepted array {!r} lacks valid source reconstruction".format(
                array_id
            )
        )

    array_interval = _interval(
        array.get("source_interval"),
        "array {!r} source_interval".format(array_id),
        source_length,
    )

    repeats = array.get("repeats")
    repeat_count = array.get("repeat_count")
    if (
        not isinstance(repeats, list)
        or not _is_plain_int(repeat_count)
        or repeat_count < 1
        or repeat_count != len(repeats)
    ):
        raise ReportExtensionError(
            "accepted array {!r} repeat_count disagrees with repeats".format(
                array_id
            )
        )

    repeat_intervals = []
    repeat_sequences = []
    previous_end = None
    for position, repeat in enumerate(repeats, 1):
        if not isinstance(repeat, dict):
            raise ReportExtensionError(
                "array {!r} repeat {} must be a mapping".format(
                    array_id,
                    position,
                )
            )
        if repeat.get("index") != position:
            raise ReportExtensionError(
                "array {!r} repeat indexes must be consecutive".format(
                    array_id
                )
            )
        interval = _interval(
            repeat.get("source_interval"),
            "array {!r} repeat {} source_interval".format(
                array_id,
                position,
            ),
            source_length,
            containing=array_interval,
        )
        sequence = repeat.get("sequence")
        if (
            not isinstance(sequence, str)
            or len(sequence) != interval["length"]
        ):
            raise ReportExtensionError(
                "array {!r} repeat {} sequence length disagrees with its "
                "source interval".format(array_id, position)
            )
        if previous_end is not None and interval["start"] <= previous_end:
            raise ReportExtensionError(
                "array {!r} repeats must be ordered and non-overlapping".format(
                    array_id
                )
            )
        repeat_intervals.append(interval)
        repeat_sequences.append(sequence)
        previous_end = interval["end"]

    if (
        repeat_intervals[0]["start"] != array_interval["start"]
        or repeat_intervals[-1]["end"] != array_interval["end"]
    ):
        raise ReportExtensionError(
            "array {!r} interval must span its first through last repeat".format(
                array_id
            )
        )

    architecture = array.get("architecture")
    expected_architecture_length = 2 * repeat_count - 1
    if (
        not isinstance(architecture, list)
        or len(architecture) != expected_architecture_length
    ):
        raise ReportExtensionError(
            "array {!r} architecture does not cover every repeat and gap".format(
                array_id
            )
        )

    spacer_gap_lengths = []
    deletion_slot_indices = []
    expected_spacers = []
    expected_deletions = []
    for repeat_position, repeat_interval in enumerate(repeat_intervals, 1):
        architecture_repeat = architecture[2 * (repeat_position - 1)]
        if (
            not isinstance(architecture_repeat, dict)
            or architecture_repeat.get("type") != "repeat"
            or architecture_repeat.get("index") != repeat_position
            or architecture_repeat.get("source_interval") != repeat_interval
            or architecture_repeat.get("sequence")
            != repeat_sequences[repeat_position - 1]
        ):
            raise ReportExtensionError(
                "array {!r} architecture repeat {} disagrees with repeats".format(
                    array_id,
                    repeat_position,
                )
            )

        if repeat_position == repeat_count:
            continue

        right = repeat_intervals[repeat_position]
        gap_length = right["start"] - repeat_interval["end"] - 1
        # Overlap was rejected above, so this is always non-negative.
        spacer_gap_lengths.append(gap_length)
        gap = architecture[2 * repeat_position - 1]
        if not isinstance(gap, dict) or gap.get("index") != repeat_position:
            raise ReportExtensionError(
                "array {!r} architecture gap {} is malformed".format(
                    array_id,
                    repeat_position,
                )
            )

        if gap_length == 0:
            expected_boundary = {
                "position": repeat_interval["end"],
                "coordinate_system": "0-based interbase",
                "reference": "input_forward_strand",
            }
            if (
                gap.get("type") != "spacer-deletion"
                or gap.get("source_boundary") != expected_boundary
                or gap.get("flanking_repeat_indices")
                != [repeat_position, repeat_position + 1]
                or gap.get("length") != 0
                or gap.get("sequence") != ""
            ):
                raise ReportExtensionError(
                    "array {!r} zero-width gap {} is not a canonical "
                    "spacer deletion".format(array_id, repeat_position)
                )
            deletion_slot_indices.append(repeat_position)
            expected_deletions.append({
                "index": repeat_position,
                "source_boundary": expected_boundary,
                "flanking_repeat_indices": [
                    repeat_position,
                    repeat_position + 1,
                ],
                "length": 0,
                "sequence": "",
            })
        else:
            expected_interval = {
                "start": repeat_interval["end"] + 1,
                "end": right["start"] - 1,
                "length": gap_length,
            }
            sequence = gap.get("sequence")
            if (
                gap.get("type") != "spacer"
                or gap.get("source_interval") != expected_interval
                or not isinstance(sequence, str)
                or len(sequence) != gap_length
            ):
                raise ReportExtensionError(
                    "array {!r} source gap {} disagrees with its "
                    "architecture spacer".format(array_id, repeat_position)
                )
            expected_spacers.append({
                "index": repeat_position,
                "source_interval": expected_interval,
                "sequence": sequence,
            })

    spacer_count = array.get("spacer_count")
    non_deleted_count = array.get("non_deleted_spacer_count")
    deletion_count = array.get("spacer_deletion_count")
    if (
        not _is_plain_int(spacer_count)
        or spacer_count != repeat_count - 1
        or non_deleted_count != len(expected_spacers)
        or deletion_count != len(expected_deletions)
        or non_deleted_count + deletion_count != spacer_count
        or array.get("spacers") != expected_spacers
        or array.get("spacer_deletions") != expected_deletions
    ):
        raise ReportExtensionError(
            "array {!r} spacer/deletion records or counts disagree with "
            "its architecture".format(array_id)
        )

    source_sequence = array.get("source_sequence")
    reconstructed_sequence = "".join(
        element.get("sequence", "") for element in architecture
    )
    if (
        not isinstance(source_sequence, str)
        or len(source_sequence) != array_interval["length"]
        or source_sequence != reconstructed_sequence
    ):
        raise ReportExtensionError(
            "array {!r} source_sequence disagrees with its validated "
            "architecture".format(array_id)
        )

    repeat_starts = [
        interval["start"] for interval in repeat_intervals
    ]
    repeat_ends = [
        interval["end"] for interval in repeat_intervals
    ]
    repeat_lengths = [
        interval["length"] for interval in repeat_intervals
    ]
    periods = [
        right - left
        for left, right in zip(repeat_starts, repeat_starts[1:])
    ]
    non_deleted_spacer_gap_lengths = [
        length for length in spacer_gap_lengths if length > 0
    ]

    return {
        "array_id": array_id,
        "category": array["category"],
        "source_interval": array_interval,
        "repeat_intervals": repeat_intervals,
        "repeat_starts_1_based": repeat_starts,
        "repeat_inclusive_ends_1_based": repeat_ends,
        "repeat_lengths_bp": repeat_lengths,
        "repeat_start_to_start_periods_bp": periods,
        "spacer_gap_lengths_bp": spacer_gap_lengths,
        "non_deleted_spacer_gap_lengths_bp": (
            non_deleted_spacer_gap_lengths
        ),
        "spacer_deletion_slot_indices": deletion_slot_indices,
        "statistics": {
            "repeat_length_bp": _descriptive_statistics(repeat_lengths),
            "repeat_start_to_start_period_bp": (
                _descriptive_statistics(periods)
            ),
            "spacer_gap_length_bp": (
                _descriptive_statistics(spacer_gap_lengths)
            ),
            "non_deleted_spacer_gap_length_bp": (
                _descriptive_statistics(non_deleted_spacer_gap_lengths)
            ),
        },
    }


def _geometry_advisories(geometry, configured_limits):
    advisories = []
    minimum = configured_limits["minimum"]
    # Benchmark-derived robust diagnostics use present spacers only.
    # Zero-width deletion slots remain in the separate all-gap statistics,
    # but must not manufacture short medians or high variability here.
    spacer_stats = geometry["statistics"][
        "non_deleted_spacer_gap_length_bp"
    ]
    if (
        minimum is not None
        and spacer_stats["count"] > 0
        and spacer_stats["mean"] >= float(minimum)
        and spacer_stats["median"] < float(minimum)
    ):
        advisories.append({
            "code": (
                "median_spacer_gap_below_configured_mean_minimum"
            ),
            "severity": "informational",
            "decision_effect": "none",
            "configured_minimum_mean_spacer_length_bp": minimum,
            "observed_mean_spacer_gap_length_bp": spacer_stats["mean"],
            "observed_median_spacer_gap_length_bp": spacer_stats["median"],
            "message": (
                "Median non-deleted source spacer length is below the "
                "configured minimum for mean spacer length while the "
                "non-deleted arithmetic mean meets it. This is a "
                "non-decision diagnostic, not a replay of the legacy "
                "all-slot detector predicate, and it does not change the "
                "array category or certainty score."
            ),
        })
    spacer_cv = spacer_stats["coefficient_of_variation"]
    if (
        spacer_cv is not None
        and spacer_cv > SPACER_GAP_CV_ADVISORY_THRESHOLD
    ):
        advisories.append({
            "code": "spacer_gap_coefficient_of_variation_above_threshold",
            "severity": "informational",
            "decision_effect": "none",
            "threshold": SPACER_GAP_CV_ADVISORY_THRESHOLD,
            "comparison": "greater_than",
            "observed_coefficient_of_variation": spacer_cv,
            "message": (
                "Non-deleted source spacer-length coefficient of variation "
                "is above 0.75. This is a non-decision diagnostic and does "
                "not change the array category or certainty score."
            ),
        })
    return advisories


def build_classification_semantics_extension(
        likely_tandem_audit_required=False):
    """Return stable wording for fields that are commonly over-interpreted."""

    return {
        "schema": {
            "name": CLASSIFICATION_SCHEMA_NAME,
            "version": CLASSIFICATION_SCHEMA_VERSION,
        },
        "status": "declarative",
        "likely_tandem_audit": {
            "required": bool(likely_tandem_audit_required),
            "extension": "likely_tandem",
        },
        "applies_to_report_schema": {
            "name": SCHEMA_NAME,
            "version": SCHEMA_VERSION,
        },
        "certainty_score": {
            "is_probability": False,
            "aggregation": {
                "method": "arithmetic_mean",
                "inputs": "selected_models_class_1_outputs",
            },
            "initial_score_bands": [
                {"condition": "score < 0.5"},
                {"condition": "0.5 <= score < 0.75"},
                {"condition": "score >= 0.75"},
            ],
            "category_assignment": {
                "bands_are_initial_partition_only": True,
                "additional_steps": [
                    "structural filtering",
                    "within-region winner competition",
                    "final low-spacer-diversity tandem demotion",
                ],
                "interpretation": (
                    "Final category is not determined by score band alone; "
                    "structural filtering and within-region winner "
                    "competition, followed by the final bounded tandem "
                    "screen, can affect category."
                ),
            },
        },
        "feature_sets": {
            "role": "classifier_inputs",
        },
        "classification_basis": {
            "status": "optional_per_array",
            "absent_meaning": (
                "Legacy row; category-assignment provenance was not "
                "explicitly recorded."
            ),
            "methods": {
                "model_score_band": (
                    "Category assignment followed the model-score band and "
                    "the existing structural and within-region cascade. "
                    "The per-array structural-filter record distinguishes a "
                    "true pass from terminal-dominance evidence that was "
                    "insufficient at small sample size. Source-backed geometry "
                    "derivation may be recorded separately and does not by "
                    "itself imply rescue."
                ),
                "source_backed_internal_repeat_completion": (
                    "A complete source-graph analysis uniquely recovered internal "
                    "repeat copies. A checked single-core replacement retained "
                    "ordinary filters and the incumbent score threshold; retained "
                    "cores remain lower evidence with their original scores."
                ),
                "source_backed_architecture_rescue": (
                    "The raw model score was retained, while acceptance was "
                    "authorized separately by the named conservative "
                    "source-backed architecture policy."
                ),
                "source_backed_possible_reconstruction": (
                    "An accepted Possible core authorized a checked longer "
                    "retained path. The result stays Possible. Its full-path "
                    "score and original core confidence are recorded separately; "
                    "this does not promote the array to Bona-fide."
                ),
                "source_backed_localized_spacer_context": (
                    "One internal block of identical spacers was reviewed in "
                    "the context of the complete path and both maximal flanks. "
                    "Each geometry retains its own score; flank checks are "
                    "correlated diagnostic context. Acceptance is capped at "
                    "Possible, with tied complete boundary alternatives retained."
                ),
                "source_backed_accepted_core_boundary_completion": (
                    "An already accepted, strictly adjudicated core authorized "
                    "locus existence. Source-only terminal evidence completed "
                    "its boundaries, and the final geometry was freshly scored "
                    "and rechecked without letting that score change the "
                    "accepted category."
                ),
                "low_spacer_diversity_tandem_demotion": (
                    "A source-replayed accepted row qualified under the "
                    "complete-link low-spacer-diversity rule. Its unmodified "
                    "score is retained, while the public category becomes "
                    "Low score with tag tandem and scientific disposition "
                    "likely_tandem."
                ),
            },
            "certainty_score_relationship": (
                "classification_basis never replaces, calibrates, or "
                "rewrites certainty_score."
            ),
        },
        "array_tags": {
            "status": "optional_sorted_identifiers",
            "tandem": (
                "Occurs only on a Low score array whose classification_basis "
                "replays the likely-tandem demotion transaction."
            ),
        },
        "validation": {
            "scope": "source_reconstruction_only",
            "interpretation": (
                "Validation confirms reconstruction from the original source "
                "sequence only; it does not establish biological truth, "
                "classification correctness, or score calibration."
            ),
        },
    }


def build_source_geometry_extension(report):
    """Build accepted-array geometry from canonical source intervals."""

    arrays, _source, source_sha256, source_length = (
        _validate_report_envelope(report)
    )
    configured_limits = _configured_spacer_limits(report)

    accepted = []
    for position, array in enumerate(arrays):
        if not isinstance(array, dict):
            raise ReportExtensionError(
                "array {} must be a mapping".format(position)
            )
        category = array.get("category")
        if not isinstance(category, str):
            raise ReportExtensionError(
                "array {} category must be a string".format(position)
            )
        if category in ACCEPTED_CATEGORIES:
            accepted.append(_repeat_geometry(array, source_length))

    array_ids = [array["array_id"] for array in accepted]
    if len(array_ids) != len(set(array_ids)):
        raise ReportExtensionError("accepted array ids must be unique")

    accepted.sort(key=lambda array: (
        array["source_interval"]["start"],
        array["source_interval"]["end"],
        array["array_id"],
    ))
    for geometry in accepted:
        geometry["source_sha256"] = source_sha256
        geometry["advisories"] = _geometry_advisories(
            geometry,
            configured_limits,
        )

    return {
        "schema": {
            "name": SOURCE_GEOMETRY_SCHEMA_NAME,
            "version": SOURCE_GEOMETRY_SCHEMA_VERSION,
        },
        "status": "computed",
        "source_sha256": source_sha256,
        "scope": {
            "accepted_categories": sorted(ACCEPTED_CATEGORIES),
            "array_count": len(accepted),
        },
        "coordinate_system": {
            "repeat_starts": "1-based input-forward-strand coordinates",
            "repeat_ends": (
                "1-based inclusive input-forward-strand coordinates"
            ),
            "spacer_deletions": (
                "zero-length slots between immediately adjacent repeats"
            ),
        },
        "derivation": {
            "inputs": (
                "source-reconstructed canonical architecture and repeat "
                "intervals only"
            ),
            "decision_effect": "none",
        },
        "diagnostic_thresholds": {
            "spacer_gap_coefficient_of_variation": {
                "threshold": SPACER_GAP_CV_ADVISORY_THRESHOLD,
                "comparison": "greater_than",
                "applies_to": (
                    "non-deleted source spacer slots only; zero-length "
                    "deletion slots remain in the separate all-gap "
                    "statistics"
                ),
                "decision_effect": "none",
            },
        },
        "configured_mean_spacer_length_limits": configured_limits,
        "arrays": accepted,
    }


def _validated_sequence_topology_extension(report, evidence):
    if not isinstance(evidence, SequenceTopologyEvidence):
        raise ReportExtensionError(
            "sequence_topology_evidence has the wrong type"
        )
    source, source_sha256, source_length = _source_metadata(report)
    if (
            evidence.source_id != source.get("id")
            or evidence.source_length != source_length
            or evidence.source_sha256 != source_sha256):
        raise ReportExtensionError(
            "sequence topology evidence does not match report source"
        )
    extension = evidence.as_dict()
    if extension.get("schema_name") != SEQUENCE_TOPOLOGY_SCHEMA_NAME:
        raise ReportExtensionError(
            "sequence topology evidence schema is unsupported"
        )
    return extension


def _zero_based_report_interval(record, label):
    interval = record.get("source_interval")
    if not isinstance(interval, dict):
        raise ReportExtensionError("{} source_interval is missing".format(label))
    start = interval.get("start")
    end = interval.get("end")
    length = interval.get("length")
    if (
            not _is_plain_int(start)
            or not _is_plain_int(end)
            or not _is_plain_int(length)
            or start < 1
            or end < start
            or length != end - start + 1):
        raise ReportExtensionError(
            "{} source_interval is invalid".format(label)
        )
    return start - 1, end


def _report_array_geometry(array):
    interval = _zero_based_report_interval(array, "array")
    repeats = array.get("repeats")
    if not isinstance(repeats, list) or not repeats:
        raise ReportExtensionError("report array repeats are invalid")
    repeat_intervals = tuple(
        _zero_based_report_interval(repeat, "repeat")
        for repeat in repeats
    )
    return interval, repeat_intervals


def _child_array_ids(report, decision):
    accepted = [
        array for array in report["arrays"]
        if array.get("category") in ACCEPTED_CATEGORIES
    ]
    by_geometry = {}
    for array in accepted:
        key = _report_array_geometry(array)
        by_geometry.setdefault(key, []).append(array.get("id"))

    child_ids_by_hash = {}
    for child in decision.child_geometries:
        key = (child.source_interval, child.repeat_intervals)
        matches = by_geometry.get(key, ())
        if (
                len(matches) != 1
                or not isinstance(matches[0], str)
                or not matches[0]):
            raise ReportExtensionError(
                "circular joined child does not identify exactly one "
                "accepted core report array"
            )
        child_ids_by_hash[child.geometry_sha256] = matches[0]

    ordered_hashes = (
        decision.proposal.linearized_view.child_geometry_sha256s_in_order
    )
    try:
        return [child_ids_by_hash[value] for value in ordered_hashes]
    except KeyError as error:
        raise ReportExtensionError(
            "circular joined view refers to an unknown child geometry"
        ) from error


def _joined_locus_record(report, decision):
    proposal = decision.proposal
    score = float(decision.score_evidence.score)
    head, tail = sorted(
        decision.child_geometries,
        key=lambda child: child.source_interval,
    )
    source_length = report["source"]["length"]
    if proposal.orientation == "forward":
        expected_locus_intervals = (
            (tail.source_interval[0], source_length),
            (0, head.source_interval[1]),
        )
    else:
        expected_locus_intervals = (
            (0, head.source_interval[1]),
            (tail.source_interval[0], source_length),
        )
    locus_segments = proposal.locus_source_segments
    if tuple(
            segment.interval for segment in locus_segments
    ) != expected_locus_intervals:
        raise ReportExtensionError(
            "circular locus segments disagree with child geometry"
        )
    if any(
            start < 0 or end < start or end > source_length
            for start, end in expected_locus_intervals):
        raise ReportExtensionError(
            "circular locus segment lies outside report source"
        )
    view = proposal.linearized_view
    return {
        "id": "CRISPR-CIRCULAR-{}".format(
            proposal.proposal_sha256[:16]
        ),
        "category": "Bona-fide" if score >= 0.75 else "Possible",
        "certainty_score": score,
        "source_traversal": {
            "direction": (
                "input_forward"
                if proposal.orientation == "forward"
                else "input_reverse_complement"
            ),
            "crosses_origin": True,
            "internal_orientation_token": proposal.orientation,
        },
        "biological_orientation": {
            "status": "not_assessed_for_joined_locus",
            "strand": None,
            "child_predictions_are_not_inherited": True,
        },
        "coordinate_system": "zero_based_half_open_segmented",
        "ordered_source_segments": [
            segment.as_dict() for segment in locus_segments
        ],
        "ordered_child_source_segments": [
            segment.as_dict()
            for segment in proposal.ordered_source_segments
        ],
        "crosses_origin": True,
        "repeat_count": len(view.repeats),
        "spacer_count": len(view.spacers),
        "repeat_source_intervals": [
            list(interval) for interval in view.repeat_source_intervals
        ],
        "junction_spacer": proposal.junction_spacer.as_dict(),
        "linearized_view_sha256": view.view_sha256,
        "linearized_content_record_sha256": view.content_sha256,
        "proposal_sha256": proposal.proposal_sha256,
        "child_array_ids_in_traversal_order": _child_array_ids(
            report, decision
        ),
        "child_arrays_retained_in_core_report": True,
        "fresh_structural_gate": decision.structural_decision.as_dict(),
        "fresh_tandem_gate": decision.tandem_decision.as_dict(),
        "fresh_score_evidence": decision.score_evidence.as_dict(),
        "validation": {
            "status": "valid",
            "source_reconstructed": True,
            "method": (
                "independent_original_sequence_slicing_with_modulo_origin_"
                "junction"
            ),
            "invalid_start_greater_than_end_interval_used": False,
        },
    }


def _build_circular_origin_loci_extension(
        report, topology_evidence, result):
    if not isinstance(result, CircularAcceptedJoinResult):
        raise ReportExtensionError(
            "circular_array_join_result has the wrong type"
        )
    _source, source_sha256, _source_length = _source_metadata(report)
    if result.source_sha256 != source_sha256:
        raise ReportExtensionError(
            "circular join result does not match report source"
        )
    if not isinstance(topology_evidence, SequenceTopologyEvidence):
        raise ReportExtensionError(
            "circular join reporting requires sequence topology evidence"
        )
    decisions = tuple(result.decisions)
    for decision in decisions:
        if (
                decision.source_sha256 != source_sha256
                or decision.topology_evidence.source_sha256 != source_sha256
                or decision.topology_evidence.topology
                != topology_evidence.topology):
            raise ReportExtensionError(
                "circular join decision has stale topology/source evidence"
            )
    joined = [
        _joined_locus_record(report, decision)
        for decision in result.joined
    ]
    return {
        "schema": {
            "name": CIRCULAR_ORIGIN_LOCI_SCHEMA_NAME,
            "version": CIRCULAR_ORIGIN_LOCI_SCHEMA_VERSION,
        },
        "source_sha256": source_sha256,
        "coordinate_system": {
            "indexing": 0,
            "end_exclusive": True,
            "reference": "input_forward_strand",
            "origin_crossing_representation": (
                "ordered_non_wrapping_source_segments"
            ),
        },
        "decision_effect": (
            "adds canonical circular-locus annotations; core linear child "
            "arrays and their classifications remain unchanged"
        ),
        "attempt_count": len(decisions),
        "joined_locus_count": len(joined),
        "joined_loci": joined,
        "attempts": [decision.as_dict() for decision in decisions],
    }


def _validated_likely_tandem_architecture(value, outer_interval):
    architecture = value.get("source_architecture")
    if (
            not isinstance(architecture, dict)
            or architecture.get("coordinate_system")
            != "zero_based_half_open"
            or architecture.get("reference") != "input_forward_strand"):
        raise ReportExtensionError(
            "likely-tandem source architecture contract is invalid"
        )
    repeats = architecture.get("repeats")
    spacers = architecture.get("spacers")
    expected_spacers = value.get("spacers_source_forward")
    if (
            not isinstance(repeats, list)
            or not isinstance(spacers, list)
            or not isinstance(expected_spacers, list)
            or len(repeats) != len(spacers) + 1
            or len(spacers) != len(expected_spacers)):
        raise ReportExtensionError(
            "likely-tandem source architecture cardinality is invalid"
        )

    outer_start = outer_interval["start"] - 1
    outer_end = outer_interval["end"]

    def validated_elements(elements, expected_sequences, label):
        result = []
        for position, element in enumerate(elements, 1):
            if (
                    not isinstance(element, dict)
                    or not _is_plain_int(element.get("index"))
                    or element.get("index") != position):
                raise ReportExtensionError(
                    "likely-tandem {} indexes are invalid".format(label)
                )
            sequence = element.get("sequence")
            interval = element.get("source_interval")
            if (
                    not isinstance(sequence, str)
                    or not sequence
                    or not isinstance(interval, dict)
                    or set(interval) != {"start", "end"}
                    or not _is_plain_int(interval.get("start"))
                    or not _is_plain_int(interval.get("end"))
                    or interval["start"] < outer_start
                    or interval["end"] <= interval["start"]
                    or interval["end"] > outer_end
                    or interval["end"] - interval["start"] != len(sequence)
                    or element.get("sequence_sha256") != hashlib.sha256(
                        sequence.encode("utf-8")
                    ).hexdigest()):
                raise ReportExtensionError(
                    "likely-tandem {} source element is invalid".format(label)
                )
            if (
                    expected_sequences is not None
                    and sequence != expected_sequences[position - 1]):
                raise ReportExtensionError(
                    "likely-tandem spacer ledger copies disagree"
                )
            result.append((interval["start"], interval["end"], sequence))
        return result

    repeat_values = validated_elements(repeats, None, "repeat")
    spacer_values = validated_elements(
        spacers, expected_spacers, "spacer"
    )
    if repeat_values[0][0] != outer_start or repeat_values[-1][1] != outer_end:
        raise ReportExtensionError(
            "likely-tandem architecture does not span its source interval"
        )
    reconstructed = []
    for index, repeat in enumerate(repeat_values):
        reconstructed.append(repeat[2])
        if index < len(spacer_values):
            spacer = spacer_values[index]
            if repeat[1] != spacer[0] or spacer[1] != repeat_values[index + 1][0]:
                raise ReportExtensionError(
                    "likely-tandem architecture is not source-contiguous"
                )
            reconstructed.append(spacer[2])
    reconstructed = "".join(reconstructed)
    if (
            len(reconstructed) != outer_interval["length"]
            or architecture.get("reconstructed_source_span_sha256")
            != hashlib.sha256(reconstructed.encode("utf-8")).hexdigest()):
        raise ReportExtensionError(
            "likely-tandem reconstructed source-span hash is invalid"
        )
    return reconstructed


def _likely_tandem_core_architecture_matches(array, ledger_record):
    """Return whether one tagged core row exactly mirrors its ledger geometry."""

    ledger_architecture = ledger_record.get("source_architecture")
    if not isinstance(array, dict) or not isinstance(
            ledger_architecture, dict):
        return False
    core_repeats = array.get("repeats")
    core_spacers = array.get("spacers")
    ledger_repeats = ledger_architecture.get("repeats")
    ledger_spacers = ledger_architecture.get("spacers")
    if not all(isinstance(items, list) for items in (
            core_repeats, core_spacers, ledger_repeats, ledger_spacers)):
        return False
    if (
            len(core_repeats) != len(ledger_repeats)
            or len(core_spacers) != len(ledger_spacers)
            or array.get("repeat_count") != len(core_repeats)
            or array.get("spacer_count") != len(core_spacers)
            or array.get("non_deleted_spacer_count") != len(core_spacers)
            or array.get("spacer_deletion_count") != 0
            or array.get("spacer_deletions") != []):
        return False

    def elements_match(core_elements, ledger_elements):
        for position, (core, ledger) in enumerate(
                zip(core_elements, ledger_elements), 1):
            core_interval = (
                core.get("source_interval")
                if isinstance(core, dict) else None
            )
            ledger_interval = (
                ledger.get("source_interval")
                if isinstance(ledger, dict) else None
            )
            if not isinstance(core_interval, dict) or not isinstance(
                    ledger_interval, dict):
                return False
            core_start = core_interval.get("start")
            core_end = core_interval.get("end")
            core_length = core_interval.get("length")
            ledger_start = ledger_interval.get("start")
            ledger_end = ledger_interval.get("end")
            sequence = core.get("sequence")
            if (
                    any(
                        not _is_plain_int(value)
                        for value in (
                            core_start, core_end, core_length,
                            ledger_start, ledger_end,
                        )
                    )
                    or not isinstance(sequence, str)
                    or not _is_plain_int(core.get("index"))
                    or not _is_plain_int(ledger.get("index"))
                    or core.get("index") != position
                    or ledger.get("index") != position
                    or sequence != ledger.get("sequence")
                    or core_start - 1 != ledger_start
                    or core_end != ledger_end
                    or core_length != len(sequence)):
                return False
        return True

    if (
            not elements_match(core_repeats, ledger_repeats)
            or not elements_match(core_spacers, ledger_spacers)):
        return False
    expected_architecture = []
    for position, repeat in enumerate(core_repeats):
        expected_architecture.append(dict(repeat, type="repeat"))
        if position < len(core_spacers):
            expected_architecture.append(dict(
                core_spacers[position], type="spacer"
            ))
    reconstructed = []
    for position, repeat in enumerate(ledger_repeats):
        reconstructed.append(repeat["sequence"])
        if position < len(ledger_spacers):
            reconstructed.append(ledger_spacers[position]["sequence"])
    reconstructed = "".join(reconstructed)
    return (
        _strict_json_equal(array.get("architecture"), expected_architecture)
        and array.get("source_sequence") == reconstructed
        and ledger_architecture.get("reconstructed_source_span_sha256")
        == hashlib.sha256(reconstructed.encode("utf-8")).hexdigest()
    )


def _validated_likely_tandem_source(report, source_sequence, record_count):
    """Bind a non-empty accepted-only ledger to the original source bytes."""

    if record_count == 0 and source_sequence is None:
        return None
    if not isinstance(source_sequence, str):
        raise ReportExtensionError(
            "non-empty likely-tandem ledgers require the original source sequence"
        )
    try:
        encoded = source_sequence.encode("ascii")
    except UnicodeEncodeError as error:
        raise ReportExtensionError(
            "likely-tandem source sequence must be ASCII"
        ) from error
    source = report.get("source", {})
    if (
            source_sequence != source_sequence.upper()
            or len(source_sequence) != source.get("length")
            or hashlib.sha256(encoded).hexdigest() != source.get("sha256")):
        raise ReportExtensionError(
            "likely-tandem source sequence does not match report source"
        )
    return source_sequence


def _validate_no_likely_tandem_accepted_arrays(arrays):
    """Fail if a qualifying row survived in an accepted report category."""

    for array in arrays:
        if array.get("category") not in ACCEPTED_CATEGORIES:
            continue
        spacers = tuple(
            element.get("sequence")
            for element in array.get("architecture", ())
            if element.get("type") in ("spacer", "spacer-deletion")
        )
        evidence = evaluate_likely_tandem(spacers)
        if evidence.likely_tandem:
            raise ReportExtensionError(
                "accepted report array {!r} still qualifies as likely "
                "tandem; quarantine is incomplete".format(array.get("id"))
            )


def _validate_exhaustive_likely_tandem_replay(
        source_categories, records, stats):
    """Bind positive and retained decisions to all pre-quarantine rows."""

    if source_categories is None:
        raise ReportExtensionError(
            "likely-tandem ledgers require pre-quarantine source categories"
        )
    try:
        _categories, expected_records, expected_stats = (
            apply_likely_tandem_quarantine(source_categories)
        )
        expected_payloads = sorted(
            json.dumps(
                record.as_dict(),
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
            for record in expected_records
        )
        observed_payloads = sorted(
            json.dumps(
                record.as_dict(),
                allow_nan=False,
                ensure_ascii=True,
                separators=(",", ":"),
                sort_keys=True,
            )
            for record in records
        )
    except (AttributeError, TypeError, ValueError, RuntimeError) as error:
        raise ReportExtensionError(
            "likely-tandem source-category replay failed"
        ) from error
    if expected_stats != stats or expected_payloads != observed_payloads:
        raise ReportExtensionError(
            "likely-tandem transaction is incomplete or internally "
            "inconsistent under exhaustive source-category replay"
        )


def _build_likely_tandem_extension(
        report, records, stats, source_sequence=None,
        source_categories=None):
    """Build a source-bound ledger for tagged Low-score demotions."""

    if not isinstance(stats, LikelyTandemQuarantineStats):
        raise ReportExtensionError("likely_tandem_stats has the wrong type")
    stats_values = stats.as_dict()
    if any(
            not _is_plain_int(value) or value < 0
            for value in stats_values.values()):
        raise ReportExtensionError(
            "likely-tandem transaction stats must be non-negative integers"
        )
    if isinstance(records, (str, bytes)):
        raise ReportExtensionError(
            "likely_tandem_records must be quarantine records"
        )
    try:
        records = tuple(records)
    except TypeError as error:
        raise ReportExtensionError(
            "likely_tandem_records must be quarantine records"
        ) from error
    if any(
            not isinstance(record, LikelyTandemQuarantineRecord)
            for record in records):
        raise ReportExtensionError(
            "likely_tandem_records contains the wrong type"
        )
    if len(records) != stats.likely_tandem_count:
        raise ReportExtensionError(
            "likely-tandem ledger count disagrees with transaction stats"
        )

    source_sequence = _validated_likely_tandem_source(
        report, source_sequence, len(records)
    )

    _source, source_sha256, source_length = _source_metadata(report)
    ordered = sorted(
        records,
        key=lambda record: (
            record.source_start_0based,
            record.source_end_0based_exclusive,
            record.source_category_index,
            record.locus_key_json,
            record.row_ordinal,
        ),
    )
    ledger = []
    observed_keys = set()
    for position, record in enumerate(ordered, 1):
        value = record.as_dict()
        evidence = value.get("evidence")
        interval = _interval(
            value.get("source_interval"),
            "likely-tandem record {} source_interval".format(position),
            source_length,
        )
        spacers = value.get("spacers_source_forward")
        if (
                value.get("action") != "quarantined"
                or value.get("disposition") != LIKELY_TANDEM_CATEGORY
                or value.get("tag") != TANDEM_TAG
                or not _is_plain_int(value.get("source_category_index"))
                or value.get("source_category_index") not in (0, 2)
                or value.get("source_category_name") not in (
                    "Bona-fide", "Possible"
                )
                or value.get("source_category_name") != {
                    0: "Bona-fide", 2: "Possible"
                }.get(value.get("source_category_index"))
                or value.get("target_category_index") != 4
                or value.get("target_category_name") != "Low score"
                or not _is_plain_int(value.get("row_ordinal"))
                or value.get("row_ordinal") < 0
                or isinstance(value.get("score"), bool)
                or not isinstance(value.get("score"), numbers.Real)
                or not math.isfinite(float(value["score"]))
                or not isinstance(evidence, dict)
                or evidence.get("decision") != LIKELY_TANDEM_CATEGORY
                or evidence.get("complete_link_threshold_satisfied") is not True
                or not isinstance(spacers, list)
                or any(not isinstance(spacer, str) or not spacer
                       for spacer in spacers)
                or value.get("repeat_count") != len(spacers) + 1):
            raise ReportExtensionError(
                "likely-tandem ledger record is internally inconsistent"
            )
        reconstructed = _validated_likely_tandem_architecture(
            value, interval
        )
        if (
                source_sequence is not None
                and source_sequence[interval["start"] - 1:interval["end"]]
                != reconstructed):
            raise ReportExtensionError(
                "likely-tandem ledger does not match its original source slice"
            )
        replayed = evaluate_likely_tandem(spacers).as_dict()
        if not _strict_json_equal(replayed, evidence):
            raise ReportExtensionError(
                "likely-tandem evidence does not replay from ledger spacers"
            )
        record_key = (
            interval["start"],
            interval["end"],
            value["source_category_index"],
            record.locus_key_json,
            value["row_ordinal"],
        )
        if record_key in observed_keys:
            raise ReportExtensionError(
                "likely-tandem ledger contains duplicate records"
            )
        observed_keys.add(record_key)
        value["id"] = "LT{}".format(position)
        value["source_interval"] = interval
        ledger.append(value)

    accepted_only = (
        report.get("run", {}).get("parameters", {}).get(
            "param_accepted_only"
        ) is True
    )
    tagged_arrays = [
        array for array in report.get("arrays", ())
        if TANDEM_TAG in array.get("tags", ())
    ]
    if accepted_only:
        if tagged_arrays:
            raise ReportExtensionError(
                "accepted-only report contains tandem-tagged Low score arrays"
            )
        for value in ledger:
            value["reported_array_id"] = None
            value["core_projection"] = "omitted_by_accepted_only"
    else:
        unmatched = list(tagged_arrays)
        for value in ledger:
            interval = value["source_interval"]
            matches = []
            for index, array in enumerate(unmatched):
                basis = array.get("classification_basis")
                transition = (
                    basis.get("transition")
                    if isinstance(basis, dict) else None
                )
                if (
                    array.get("category") == "Low score"
                    and isinstance(transition, dict)
                    and basis.get("method") == CLASSIFICATION_BASIS_METHOD
                    and basis.get("disposition") == LIKELY_TANDEM_CATEGORY
                    and basis.get("tag") == TANDEM_TAG
                    and _strict_json_equal(
                        basis.get("likely_tandem"), value["evidence"]
                    )
                    and transition.get("rule_name")
                    == TRANSACTION_RULE_NAME
                    and transition.get("rule_version")
                    == TRANSACTION_RULE_VERSION
                    and _strict_json_equal(
                        transition.get("source_category_index"),
                        value["source_category_index"],
                    )
                    and transition.get("source_category_name")
                    == value["source_category_name"]
                    and transition.get("target_category_index") == 4
                    and transition.get("target_category_name")
                    == "Low score"
                    and _strict_json_equal(
                        transition.get("row_ordinal"), value["row_ordinal"]
                    )
                    and _strict_json_equal(
                        transition.get("locus_key"), value.get("locus_key")
                    )
                    and transition.get("reported_disposition")
                    == LIKELY_TANDEM_CATEGORY
                    and transition.get("reported_tag") == TANDEM_TAG
                    and _strict_json_equal(
                        array.get("legacy_interval_key"),
                        value.get("locus_key"),
                    )
                    and _strict_json_equal(
                        array.get("source_interval"), interval
                    )
                    and _likely_tandem_core_architecture_matches(
                        array, value
                    )
                    and not isinstance(
                        array.get("certainty_score"), bool
                    )
                    and isinstance(
                        array.get("certainty_score"), numbers.Real
                    )
                    and math.isfinite(float(array["certainty_score"]))
                    and _strict_json_equal(
                        array.get("certainty_score"), value["score"]
                    )
                ):
                    matches.append(index)
            if len(matches) != 1:
                raise ReportExtensionError(
                    "likely-tandem ledger does not map one-to-one to a "
                    "tagged Low score array"
                )
            array = unmatched.pop(matches[0])
            value["reported_array_id"] = array["id"]
            value["core_projection"] = "reported_as_low_score"
        if unmatched:
            raise ReportExtensionError(
                "report has tandem-tagged Low score arrays without ledger rows"
            )

    _validate_exhaustive_likely_tandem_replay(
        source_categories, records, stats
    )

    source_category_counts = {
        "Bona-fide": sum(
            value["source_category_name"] == "Bona-fide"
            for value in ledger
        ),
        "Possible": sum(
            value["source_category_name"] == "Possible"
            for value in ledger
        ),
    }
    category_counts = report.get("summary", {}).get("category_counts", {})
    expected_reported_low_score = (
        0 if accepted_only else stats.diagnostic_rows_after
    )
    if (
            stats.accepted_rows_after
            != stats.accepted_rows_before - len(ledger)
            or stats.accepted_rows_reviewed != stats.accepted_rows_before
            or stats.retained_count != stats.accepted_rows_after
            or stats.bona_fide_quarantined
            + stats.possible_quarantined != len(ledger)
            or stats.bona_fide_quarantined
            != source_category_counts["Bona-fide"]
            or stats.possible_quarantined
            != source_category_counts["Possible"]
            or stats.total_rows_before != stats.total_rows_after
            or stats.diagnostic_rows_after
            != stats.diagnostic_rows_before + len(ledger)
            or report["summary"]["accepted_array_count"]
            != stats.accepted_rows_after
            or category_counts.get("Low score")
            != expected_reported_low_score):
        raise ReportExtensionError(
            "likely-tandem transaction stats are internally inconsistent"
        )
    return {
        "schema": {
            "name": LIKELY_TANDEM_SCHEMA_NAME,
            "version": LIKELY_TANDEM_SCHEMA_VERSION,
        },
        "source_sha256": source_sha256,
        "disposition": LIKELY_TANDEM_CATEGORY,
        "tag": TANDEM_TAG,
        "reported_category": "Low score",
        "status": "applied",
        "rule": {
            "name": LIKELY_TANDEM_RULE_NAME,
            "version": LIKELY_TANDEM_RULE_VERSION,
            "linkage": "complete_link_all_spacer_pairs",
            "identity": "1 - global_levenshtein/max(sequence_lengths)",
            "policy": LIKELY_TANDEM_POLICY.as_dict(),
        },
        "decision_effect": (
            "demotes qualifying accepted arrays to tagged Low score output; "
            "accepted-only projection omits them while preserving this ledger"
        ),
        "external_tool_support_used_for_decision": False,
        "transaction_stats": stats.as_dict(),
        "array_count": len(ledger),
        "arrays": ledger,
    }


def build_report_extensions(
        report,
        include_biological_loci=False,
        sequence_topology_evidence=None,
        circular_array_join_result=None,
        likely_tandem_records=None,
        likely_tandem_stats=None,
        likely_tandem_source_categories=None,
        source_sequence=None,
        competing_terminal_extent_context=None):
    """Return an ``extensions`` mapping without changing ``report``."""

    arrays, source, _source_sha256, _source_length = (
        _validate_report_envelope(report)
    )
    tandem_inputs = (
        likely_tandem_records,
        likely_tandem_stats,
        likely_tandem_source_categories,
    )
    extensions = {
        "classification_semantics": (
            build_classification_semantics_extension(
                likely_tandem_audit_required=all(
                    value is not None for value in tandem_inputs
                )
            )
        ),
        "source_geometry": build_source_geometry_extension(report),
    }
    if any(value is not None for value in tandem_inputs):
        _validate_no_likely_tandem_accepted_arrays(arrays)
    if include_biological_loci:
        extensions["biological_loci"] = (
            build_biological_loci_extension(arrays, source)
        )
    if sequence_topology_evidence is not None:
        extensions["sequence_topology"] = (
            _validated_sequence_topology_extension(
                report, sequence_topology_evidence
            )
        )
    if circular_array_join_result is not None:
        extensions["circular_origin_loci"] = (
            _build_circular_origin_loci_extension(
                report,
                sequence_topology_evidence,
                circular_array_join_result,
            )
        )
    if ((likely_tandem_records is None)
            != (likely_tandem_stats is None)):
        raise ReportExtensionError(
            "likely-tandem records and stats must be supplied together"
        )
    if (
            likely_tandem_records is None
            and report.get("run", {}).get("parameters", {}).get(
                "param_accepted_only"
            ) is True):
        raise ReportExtensionError(
            "accepted-only reports require the likely-tandem audit ledger"
        )
    if (
        likely_tandem_records is None
        and any(
            TANDEM_TAG in _validated_array_tags(array, "array")
            for array in arrays
        )
    ):
        raise ReportExtensionError(
            "tandem-tagged arrays require the likely-tandem ledger"
        )
    if likely_tandem_records is not None:
        extensions["likely_tandem"] = _build_likely_tandem_extension(
            report,
            likely_tandem_records,
            likely_tandem_stats,
            source_sequence=source_sequence,
            source_categories=likely_tandem_source_categories,
        )
    if competing_terminal_extent_context is not None:
        competing = build_competing_terminal_extents(
            report, source_sequence=source_sequence,
            **competing_terminal_extent_context)
        if competing is not None:
            extensions["competing_terminal_extents"] = competing
    return extensions


def attach_report_extensions(
        report,
        include_biological_loci=False,
        sequence_topology_evidence=None,
        circular_array_join_result=None,
        likely_tandem_records=None,
        likely_tandem_stats=None,
        likely_tandem_source_categories=None,
        source_sequence=None,
        competing_terminal_extent_context=None):
    """Return a deep-copied report with one new top-level extension object."""

    if not isinstance(report, dict):
        raise ReportExtensionError("report must be a mapping")
    if "extensions" in report:
        raise ReportExtensionError(
            "report already contains an extensions object"
        )

    extensions = build_report_extensions(
        report,
        include_biological_loci=include_biological_loci,
        sequence_topology_evidence=sequence_topology_evidence,
        circular_array_join_result=circular_array_join_result,
        likely_tandem_records=likely_tandem_records,
        likely_tandem_stats=likely_tandem_stats,
        likely_tandem_source_categories=likely_tandem_source_categories,
        source_sequence=source_sequence,
        competing_terminal_extent_context=competing_terminal_extent_context,
    )
    result = copy.deepcopy(report)
    result["extensions"] = extensions
    return result


__all__ = [
    "ReportExtensionError",
    "attach_report_extensions",
    "build_classification_semantics_extension",
    "build_report_extensions",
    "build_source_geometry_extension",
]
