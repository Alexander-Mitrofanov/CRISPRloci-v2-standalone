"""Fail-closed validation of source-backed report provenance.

The recovery transaction deliberately serializes runtime-object-free evidence.
This module treats that evidence as untrusted when a canonical report is
built.  Every serialized source interval and digest is replayed against the
current input sequence.  Strict-existence decisions and conflict topology are
also independently recomputed from the serialized geometry.

The direct-tandem record does not serialize the complete family candidates or
all source-CDS annotations used by its original decision.  Report validation
therefore verifies its fixed rule, source binding, every interval/digest, and
its internal interval relationships; it does not claim to rerun peptide
selection from inputs that are no longer present.
"""

from dataclasses import asdict, dataclass, field
import hashlib
import json
import math
from numbers import Real
from types import SimpleNamespace

import components.accepted_family_tandem_validator as tandem_module
from components.architecture_rescue_coding_context import SourceCDS
from components.architecture_rescue_policy import (
    AUDITED_ARCHITECTURE_RESCUE_POLICY,
)
from components.accepted_source_intact_conflict import (
    ACCEPTED_SOURCE_CONFLICT_POLICY,
    _relationship as recompute_conflict_relationship,
)
from components.boundary_assessment import (
    AlternativeRelationshipInput,
    BoundaryAssessmentNormalizer,
    BoundaryDimensionInput,
    CandidateGeometryInput,
    DIMENSION_OCCURRENCE_PATH,
    DIMENSION_REPEAT_PHASE,
    STATUS_AMBIGUOUS,
    STATUS_POLICY_SELECTED,
    geometry_input_from_candidate,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
)
from components.reconstruction_validator import validate_candidate_reconstruction
from components.source_intact_existence_rescue import (
    SOURCE_INTACT_EXISTENCE_POLICY,
    SourceIntactEvaluationContext,
    evaluate_source_intact_existence_in_context,
)
from components.source_intact_terminal_authorization import (
    GUARDED_FOUR_SPACER,
    REJECTED as TERMINAL_AUTHORIZATION_REJECTED,
    SOURCE_INTACT_TERMINAL_AUTHORIZATION_RULE,
    STATISTICALLY_PASSED,
)
from components.source_periodic_lattice_repair import (
    SourcePeriodicLatticeRule,
    SourcePeriodicLatticeSafetyEvidence,
    propose_source_periodic_lattice_repair,
    replay_source_periodic_lattice_proposal,
    source_periodic_lattice_geometry_prefilter,
)
from components.source_regular_subpath_recovery import (
    DIAGNOSTIC_PARENT,
    SourceRegularSubpathContext,
    propose_source_regular_subpaths_in_context,
)
from components.source_regular_subpath_recovery_transaction import (
    RegularSubpathProposalTandemEvidence,
    SourceRegularSubpathDiagnosticRowContext,
    canonicalize_regular_subpath_row_value,
    plan_source_regular_subpath_recovery_batch,
)
from components.likely_tandem import (
    CATEGORY_NAME as LIKELY_TANDEM_DISPOSITION,
    evaluate_likely_tandem,
)
from components.likely_tandem_transaction import (
    CLASSIFICATION_BASIS_METHOD as LIKELY_TANDEM_METHOD,
    TANDEM_TAG,
    TRANSACTION_RULE_NAME as LIKELY_TANDEM_TRANSACTION_RULE_NAME,
    TRANSACTION_RULE_VERSION as LIKELY_TANDEM_TRANSACTION_RULE_VERSION,
)
from components.fast_accepted_family_rescue import (
    CLASSIFICATION_METHOD as SOURCE_ACCEPTED_FAMILY_METHOD,
    replay_classification_evidence as replay_accepted_family_evidence,
)
from components.source_superlocus_augmentation import (
    replay_source_superlocus_path_local_convergence_proof,
)
from components.internal_repeat_completion import InternalSearchContext
from components.internal_local_core_bridge import validate_declared_family_against_run
from components.internal_repeat_transaction import (
    METHOD as INTERNAL_REPEAT_COMPLETION_METHOD, replay_internal_completion_basis,
    replay_retained_internal_core,
)
from components.terminal_repair_completion import (
    METHOD as TERMINAL_REPAIR_COMPLETION_METHOD,
    replay_terminal_completion_basis,
)


SOURCE_RESCUE_METHOD = "source_backed_intact_existence_rescue"
SOURCE_RETIREMENT_METHOD = (
    "source_backed_intact_terminal_conflict_retirement"
)
SOURCE_PERIODIC_LATTICE_METHOD = "source_periodic_lattice_repair"
SOURCE_ARCHITECTURE_RESCUE_METHOD = "source_backed_architecture_rescue"
SOURCE_REGULAR_SUBPATH_METHOD = "source_backed_regular_subpath_recovery"
SOURCE_METHODS = frozenset((
    SOURCE_RESCUE_METHOD,
    SOURCE_RETIREMENT_METHOD,
    SOURCE_PERIODIC_LATTICE_METHOD,
    SOURCE_ARCHITECTURE_RESCUE_METHOD,
    SOURCE_REGULAR_SUBPATH_METHOD,
    SOURCE_ACCEPTED_FAMILY_METHOD,
    LIKELY_TANDEM_METHOD,
))
MAX_CLASSIFICATION_BASIS_DEPTH = 8

_SOURCE_RECORD_RULE = {
    "name": "source-intact-recovery-transaction",
    "version": "1.1.0",
}
_COORDINATE_SYSTEM = "zero_based_half_open"
_SHA256_LENGTH = 64
_CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
_ACCEPTED_CATEGORY_INDEXES = (0, 2)
_RECOVERY_SOURCE_CATEGORY_INDEXES = (3, 4)

_PERIODIC_TRANSACTION_RULE = {
    "name": "source-periodic-lattice-category-transaction",
    "version": "1.0.0",
}
_PERIODIC_CLASSIFICATION_CEILING = "Possible"
_PERIODIC_ACTION = "repaired_possible"
_REGULAR_TRANSACTION_RULE = {
    "name": "source-regular-subpath-recovery-transaction",
    "version": "1.2.0",
}
_REGULAR_SNAPSHOT_SCHEMA = {
    "name": "crispridentify.source_regular_subpath_report_replay_snapshot",
    "version": "1.0.0",
}
_REGULAR_ROW_NORMALIZER_VERSION = "1.0.0"
_REGULAR_BASIS_KEYS = frozenset((
    "method",
    "rule",
    "source_sha256",
    "classification_cap",
    "bona_fide_admission_allowed",
    "exact_boundary_certified",
    "boundary_assessment_id",
    "report_replay_snapshot",
    "parent",
    "proposal",
    "tandem_validation",
    "accepted_overlap_protection",
    "fresh_classifier",
))
_REGULAR_SNAPSHOT_KEYS = frozenset((
    "schema",
    "transaction_rule",
    "coordinate_system",
    "source_sha256",
    "source_length",
    "source_cds",
    "accepted_geometries",
    "reviewed_diagnostic_parents",
    "pretransaction_row_manifest",
    "transaction_plan",
    "bindings",
    "invariants",
    "snapshot_sha256",
))
_REGULAR_GEOMETRY_KEYS = frozenset((
    "coordinate_system",
    "source_interval",
    "repeat_count",
    "repeat_intervals",
    "repeat_starts",
    "repeat_sequences",
    "spacer_sequences",
    "consensus",
    "source_slice_sha256",
    "geometry_sha256",
))
_REGULAR_FAMILY_KEYS = frozenset((
    "caller_declared_complete",
    "candidate_count",
    "candidate_geometries",
))
_REGULAR_MANIFEST_KEYS = frozenset((
    "category_index",
    "category_name",
    "locus_key",
    "row_ordinal",
    "row_field_count",
    "score",
    "geometry",
    "features",
    "features_sha256",
    "metadata",
    "metadata_sha256",
    "row_binding_sha256",
))
_ARCHITECTURE_RESCUE_KEYS = frozenset((
    "method",
    "source_sha256",
    "geometry_derivation",
    "raw_score_band",
    "policy",
    "decision",
    "source_superlocus_id",
    "source_superlocus_convergence",
    "lineage_path_id",
    "parent_path_ids",
    "family_path_ids",
    "seed_roles",
    "transition_action",
    "transition_reason",
    "dominated_paths",
    "target_locus_key",
    "migrated_locus_keys",
    "coding_context",
    "source_phase_arbitration",
))
_ARCHITECTURE_CONVERGENCE_KEYS = frozenset((
    "status",
    "reason",
    "path_local_status",
    "path_local_reason",
    "path_local_proof",
    "architecture_commit_eligible",
))
_ARCHITECTURE_PHASE_KEYS = frozenset((
    "policy",
    "coordinate_system",
    "superlocus_id",
    "status",
    "reason",
    "winner_path_id",
    "representative_path_id",
    "boundary_status",
    "reported_geometry_role",
    "retained_path_ids",
    "removed_path_ids",
    "competitor_path_ids",
    "source_evidence",
))
_PERIODIC_RECORD_KEYS = frozenset((
    "rule",
    "action",
    "reason_codes",
    "coordinate_system",
    "classification_ceiling",
    "source_category_index",
    "source_category_name",
    "source_locus_key",
    "source_row_ordinal",
    "target_category_index",
    "target_category_name",
    "target_locus_key",
    "parent_score",
    "fresh_score",
    "parent_repeat_intervals",
    "representative_repeat_intervals",
    "geometry_prefilter",
    "parent_direct_tandem",
    "proposal",
    "complete_family_tandem_reviews",
    "fresh_scoring_metadata",
    "boundary_assessment_id",
    "boundary_status",
    "exact_boundary_certified",
    "collision_with",
))
_PERIODIC_BASIS_KEYS = frozenset((
    "method",
    "prior_classification_basis",
    "prior_boundary_assessment",
    "source_periodic_lattice_repair",
    "transition",
    "classification_ceiling",
    "boundary_status",
    "exact_boundary_certified",
))
_PERIODIC_REVIEW_KEYS = frozenset((
    "role",
    "alternative_index",
    "repeat_intervals",
    "decision",
))

_GEOMETRY_KEYS = frozenset((
    "geometry_id",
    "interval",
    "repeat_intervals",
    "spacer_intervals",
    "repeat_count",
    "spacer_count",
    "source_slice_sha256",
    "repeat_source_slice_sha256s",
    "spacer_source_slice_sha256s",
))
_EXISTENCE_KEYS = frozenset((
    "decision",
    "eligible",
    "rejection_codes",
    "reasons",
    "metrics",
    "boundary_status",
    "policy",
))
_EVIDENCE_KEYS = frozenset((
    "source_sha256",
    "geometry",
    "family_geometries",
    "existence",
    "structural_filter",
    "terminal_authorization",
    "supplied_tandem",
    "recomputed_tandem",
    "invalid_evidence_codes",
    "authorized",
    "terminal_status_required",
    "boundary_status",
    "exact_boundary_certified",
    "classifier_score_used",
    "candidate_length_preference_used",
))
_TERMINAL_AUTHORIZATION_KEYS = frozenset((
    "authorized",
    "mode",
    "reason_codes",
    "observed_terminal_status",
    "nonempty_spacer_count",
    "begin_dominance",
    "end_dominance",
    "source_intact_existence_eligible",
    "structural_filter_passed",
    "rule",
))
_TARGET_KEYS = frozenset((
    "category_index",
    "category_name",
    "locus_key",
    "row_ordinal",
    "geometry",
))
_RECORD_KEYS = frozenset((
    "rule",
    "action",
    "reason_codes",
    "context_key",
    "source_category_index",
    "source_category_name",
    "source_locus_key",
    "source_row_ordinal",
    "target_category_index",
    "target_category_name",
    "target_locus_key",
    "score",
    "evidence",
    "accepted_targets",
    "conflict_decisions",
    "alternative_relationships",
    "boundary_status",
    "exact_boundary_certified",
    "classifier_score_used",
    "candidate_length_preference_used",
))
_ALTERNATIVE_KEYS = frozenset((
    "source_sha256",
    "relationship",
    "overlap_interval",
    "first_geometry",
    "second_geometry",
    "boundary_status",
    "representative_selected",
    "classifier_score_compared",
    "candidate_length_preference_used",
))
_CONFLICT_KEYS = frozenset((
    "decision",
    "reason",
    "source_authorized",
    "accepted_action",
    "source_intact_action",
    "boundary_status",
    "accepted_category",
    "source_sha256",
    "accepted_snapshot",
    "source_intact_snapshot",
    "accepted_architecture",
    "source_intact_architecture",
    "relationship",
    "accepted_tandem",
    "source_intact_tandem",
    "invalid_evidence_codes",
    "policy",
    "coordinate_system",
    "classifier_score_compared",
    "candidate_length_preference_used",
    "exact_boundary_certified",
))
_ARCHITECTURE_KEYS = frozenset((
    "strict_existence_eligible",
    "strict_existence_rejection_codes",
    "weakness_domains",
    "structural_filter_passed",
    "structural_filter_terminal_status",
    "strict_existence",
    "structural_filter",
))
_RELATIONSHIP_KEYS = frozenset((
    "topology",
    "overlap_interval",
    "overlap_bp",
    "accepted_overlap_fraction",
    "source_intact_overlap_fraction",
    "exact_shared_repeat_intervals",
    "repeat_overlap_pairs",
    "accepted_exact_repeat_subpath",
    "source_intact_exact_repeat_subpath",
    "single_facing_terminal_repeat_collision",
))
_SNAPSHOT_KEYS = _GEOMETRY_KEYS.difference(("geometry_id",))


class ReportProvenanceError(ValueError):
    """A serialized classification claim cannot be independently replayed."""


@dataclass(frozen=True)
class _Geometry:
    mapping: dict
    interval: tuple
    repeat_intervals: tuple
    spacer_intervals: tuple
    geometry_id: str


@dataclass(frozen=True)
class ValidatedSourceIntactAlternatives:
    """Canonical records plus boundary relationships proven by those records."""

    records: tuple
    boundary_relationships: tuple


def _canonical_json(value):
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )


def _sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _mapping(value, label):
    if not isinstance(value, dict):
        raise ReportProvenanceError("{} must be a mapping".format(label))
    return value


def _sequence(value, label):
    if isinstance(value, (str, bytes)) or not isinstance(value, (list, tuple)):
        raise ReportProvenanceError("{} must be a sequence".format(label))
    return tuple(value)


def _keys(mapping, expected, label):
    _mapping(mapping, label)
    observed = frozenset(mapping)
    if observed != expected:
        missing = sorted(expected.difference(observed))
        extra = sorted(observed.difference(expected))
        raise ReportProvenanceError(
            "{} has noncanonical keys (missing={}, extra={})".format(
                label, missing, extra
            )
        )


def _bool(value, label):
    if not isinstance(value, bool):
        raise ReportProvenanceError("{} must be boolean".format(label))
    return value


def _integer(value, label, minimum=0):
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < minimum
    ):
        raise ReportProvenanceError(
            "{} must be an integer >= {}".format(label, minimum)
        )
    return value


def _finite(value, label):
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not math.isfinite(float(value))
    ):
        raise ReportProvenanceError("{} must be finite".format(label))
    return float(value)


def _hash(value, label):
    if (
        not isinstance(value, str)
        or len(value) != _SHA256_LENGTH
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ReportProvenanceError(
            "{} must be a lowercase hexadecimal SHA-256".format(label)
        )
    return value


@dataclass(frozen=True)
class ReportProvenanceValidationContext:
    """One source binding and existence context reused for a whole report."""

    source_sequence: str
    source_sha256: str
    existence_context: SourceIntactEvaluationContext
    _source_input: object = field(repr=False, compare=False)
    _categories: object = field(default=None, repr=False, compare=False)
    _accepted_only_projection: bool = field(
        default=False,
        repr=False,
        compare=False,
    )
    _internal_completion_evidence: object = field(default=None, repr=False, compare=False)
    _occurrence_replay_context: dict = field(default_factory=dict, repr=False, compare=False)
    _internal_replay_context: dict = field(default_factory=dict, repr=False, compare=False)
    _regular_replay_cache: dict = field(
        default_factory=dict,
        repr=False,
        compare=False,
    )

    @classmethod
    def from_source(
        cls,
        source_sequence,
        expected_source_sha256,
        *,
        categories=None,
        accepted_only_projection=False,
        internal_completion_evidence=None,
    ):
        if not isinstance(accepted_only_projection, bool):
            raise ReportProvenanceError(
                "accepted_only_projection must be an explicit boolean"
            )
        expected = _hash(
            expected_source_sha256, "expected_source_sha256"
        )
        existence_context = SourceIntactEvaluationContext(source_sequence)
        if existence_context.source_sha256 != expected:
            raise ReportProvenanceError(
                "report provenance source binding is stale"
            )
        return cls(
            source_sequence=existence_context.source_sequence,
            source_sha256=existence_context.source_sha256,
            existence_context=existence_context,
            _source_input=source_sequence,
            _categories=categories,
            _accepted_only_projection=accepted_only_projection,
            _internal_completion_evidence=internal_completion_evidence,
        )

    def ordinary_occurrence_witnesses(self):
        # Report-scoped value keys, not another owner of candidate objects.
        from components.initial_competition import validate_ordinary_witness_inventory
        if "keys" not in self._occurrence_replay_context:
            self._occurrence_replay_context["keys"] = validate_ordinary_witness_inventory(
                self.categories, self.source_sequence)
        return self._occurrence_replay_context["keys"]

    def assert_bound_to(self, source_sequence, expected_source_sha256):
        if expected_source_sha256 != self.source_sha256:
            raise ReportProvenanceError(
                "report provenance context carries a stale source hash"
            )
        if (
            source_sequence is self._source_input
            or source_sequence is self.source_sequence
        ):
            return self.source_sequence
        if (
            not isinstance(source_sequence, str)
            or source_sequence.upper() != self.source_sequence
        ):
            raise ReportProvenanceError(
                "report provenance context is bound to another source"
            )
        return self.source_sequence

    @property
    def categories(self):
        if self._categories is None:
            raise ReportProvenanceError(
                "report provenance context has no category graph"
            )
        return self._categories

    @property
    def accepted_only_projection(self):
        return self._accepted_only_projection


def _validation_context(
    source_sequence,
    expected_source_sha256,
    validation_context,
):
    if validation_context is None:
        return ReportProvenanceValidationContext.from_source(
            source_sequence, expected_source_sha256
        )
    if not isinstance(
        validation_context, ReportProvenanceValidationContext
    ):
        raise ReportProvenanceError(
            "validation_context must be a ReportProvenanceValidationContext"
        )
    validation_context.assert_bound_to(
        source_sequence, expected_source_sha256
    )
    return validation_context


def _interval(value, source_length, label, allow_none=False):
    if value is None and allow_none:
        return None
    values = _sequence(value, label)
    if len(values) != 2:
        raise ReportProvenanceError(
            "{} must contain exactly start and end".format(label)
        )
    start = _integer(values[0], "{}.start".format(label))
    end = _integer(values[1], "{}.end".format(label))
    if end <= start or end > source_length:
        raise ReportProvenanceError(
            "{} is not a non-empty in-bounds half-open interval".format(label)
        )
    return start, end


def _intervals(value, source_length, label):
    return tuple(
        _interval(item, source_length, "{}[{}]".format(label, index))
        for index, item in enumerate(_sequence(value, label))
    )


def _expected_geometry_id(repeats, spacers, repeat_hashes, spacer_hashes):
    payload = {
        "coordinate_system": _COORDINATE_SYSTEM,
        "repeat_intervals": repeats,
        "spacer_intervals": spacers,
        "repeat_hashes": repeat_hashes,
        "spacer_hashes": spacer_hashes,
    }
    return _sha256(_canonical_json(payload))


def _validate_geometry(source, mapping, label, *, snapshot=False):
    expected_keys = _SNAPSHOT_KEYS if snapshot else _GEOMETRY_KEYS
    _keys(mapping, expected_keys, label)
    interval = _interval(mapping["interval"], len(source), label + ".interval")
    repeats = _intervals(
        mapping["repeat_intervals"],
        len(source),
        label + ".repeat_intervals",
    )
    spacers = _intervals(
        mapping["spacer_intervals"],
        len(source),
        label + ".spacer_intervals",
    )
    if not repeats:
        raise ReportProvenanceError(label + " must contain repeats")
    if any(left[1] > right[0] for left, right in zip(repeats, repeats[1:])):
        raise ReportProvenanceError(label + " repeats overlap or are unsorted")
    expected_spacers = tuple(
        (left[1], right[0]) for left, right in zip(repeats, repeats[1:])
    )
    if spacers != expected_spacers:
        raise ReportProvenanceError(
            label + " spacers are not the exact inter-repeat source gaps"
        )
    if interval != (repeats[0][0], repeats[-1][1]):
        raise ReportProvenanceError(
            label + " envelope does not match its terminal repeats"
        )
    if mapping["repeat_count"] != len(repeats):
        raise ReportProvenanceError(label + " repeat_count is stale")
    if mapping["spacer_count"] != len(spacers):
        raise ReportProvenanceError(label + " spacer_count is stale")

    repeat_hashes = tuple(
        _hash(item, label + ".repeat_source_slice_sha256s")
        for item in _sequence(
            mapping["repeat_source_slice_sha256s"],
            label + ".repeat_source_slice_sha256s",
        )
    )
    spacer_hashes = tuple(
        _hash(item, label + ".spacer_source_slice_sha256s")
        for item in _sequence(
            mapping["spacer_source_slice_sha256s"],
            label + ".spacer_source_slice_sha256s",
        )
    )
    expected_repeat_hashes = tuple(
        _sha256(source[start:end]) for start, end in repeats
    )
    expected_spacer_hashes = tuple(
        _sha256(source[start:end]) for start, end in spacers
    )
    if repeat_hashes != expected_repeat_hashes:
        raise ReportProvenanceError(label + " repeat source hashes are stale")
    if spacer_hashes != expected_spacer_hashes:
        raise ReportProvenanceError(label + " spacer source hashes are stale")
    if (
        _hash(mapping["source_slice_sha256"], label + ".source_slice_sha256")
        != _sha256(source[interval[0]:interval[1]])
    ):
        raise ReportProvenanceError(label + " source slice hash is stale")

    geometry_id = (
        _expected_geometry_id(
            repeats, spacers, repeat_hashes, spacer_hashes
        )
        if snapshot
        else _hash(mapping["geometry_id"], label + ".geometry_id")
    )
    if (
        not snapshot
        and geometry_id
        != _expected_geometry_id(
            repeats, spacers, repeat_hashes, spacer_hashes
        )
    ):
        raise ReportProvenanceError(label + " geometry_id is stale")
    return _Geometry(mapping, interval, repeats, spacers, geometry_id)


def _geometry_from_candidate(source, candidate):
    result = validate_candidate_reconstruction(source, candidate)
    if not result.valid:
        raise ReportProvenanceError(
            "reported candidate no longer reconstructs from source"
        )
    repeats = tuple(candidate.list_repeats)
    starts = tuple(result.genomic_repeat_starts)
    repeat_intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    spacer_intervals = tuple(
        (left[1], right[0])
        for left, right in zip(repeat_intervals, repeat_intervals[1:])
    )
    repeat_hashes = tuple(
        _sha256(source[start:end]) for start, end in repeat_intervals
    )
    spacer_hashes = tuple(
        _sha256(source[start:end]) for start, end in spacer_intervals
    )
    interval = (repeat_intervals[0][0], repeat_intervals[-1][1])
    mapping = {
        "geometry_id": _expected_geometry_id(
            repeat_intervals,
            spacer_intervals,
            repeat_hashes,
            spacer_hashes,
        ),
        "interval": list(interval),
        "repeat_intervals": [list(item) for item in repeat_intervals],
        "spacer_intervals": [list(item) for item in spacer_intervals],
        "repeat_count": len(repeat_intervals),
        "spacer_count": len(spacer_intervals),
        "source_slice_sha256": _sha256(source[interval[0]:interval[1]]),
        "repeat_source_slice_sha256s": list(repeat_hashes),
        "spacer_source_slice_sha256s": list(spacer_hashes),
    }
    return _Geometry(
        mapping,
        interval,
        repeat_intervals,
        spacer_intervals,
        mapping["geometry_id"],
    )


def _same_geometry(first, second):
    return (
        first.interval == second.interval
        and first.repeat_intervals == second.repeat_intervals
        and first.spacer_intervals == second.spacer_intervals
        and first.geometry_id == second.geometry_id
    )


def _dummy_candidate(source, geometry, consensus_length):
    repeats = tuple(source[start:end] for start, end in geometry.repeat_intervals)
    spacers = tuple(source[start:end] for start, end in geometry.spacer_intervals)
    return SimpleNamespace(
        list_repeat_starts=tuple(start for start, _end in geometry.repeat_intervals),
        list_repeats=repeats,
        list_spacers=spacers,
        consensus="A" * consensus_length,
    )


def _validate_existence(
    source,
    source_hash,
    mapping,
    geometry,
    label,
    validation_context,
):
    _keys(mapping, _EXISTENCE_KEYS, label)
    policy = _mapping(mapping["policy"], label + ".policy")
    if policy != SOURCE_INTACT_EXISTENCE_POLICY.as_dict():
        raise ReportProvenanceError(label + " uses a stale or altered policy")
    metrics = _mapping(mapping["metrics"], label + ".metrics")
    consensus_length = _integer(
        metrics.get("consensus_length"),
        label + ".metrics.consensus_length",
        1,
    )
    candidate = _dummy_candidate(source, geometry, consensus_length)
    expected = evaluate_source_intact_existence_in_context(
        validation_context.existence_context,
        candidate,
        policy=SOURCE_INTACT_EXISTENCE_POLICY,
    ).as_dict()
    if _canonical_json(mapping) != _canonical_json(expected):
        raise ReportProvenanceError(
            label + " does not equal independently recomputed source evidence"
        )
    if metrics.get("source_sha256") != source_hash:
        raise ReportProvenanceError(label + " source SHA-256 is stale")
    return mapping


def _validate_terminal(source, geometry, terminal, label):
    expected_keys = frozenset((
        "status",
        "reason",
        "nonempty_spacer_count",
        "minimum_sample_size",
        "maximum_allowed_dominance",
        "begin_dominance",
        "end_dominance",
    ))
    _keys(terminal, expected_keys, label)
    spacers = [
        source[start:end] for start, end in geometry.spacer_intervals
        if end > start
    ]
    if terminal["nonempty_spacer_count"] != len(spacers):
        raise ReportProvenanceError(label + " spacer count is stale")
    if terminal["minimum_sample_size"] != 5:
        raise ReportProvenanceError(label + " sample threshold is stale")
    if terminal["maximum_allowed_dominance"] != 0.8:
        raise ReportProvenanceError(label + " dominance threshold is stale")
    if not spacers:
        begin = end = None
        status, reason = "failed", "no_nonempty_spacers"
    else:
        def dominance(values):
            return max(values.count(item) for item in set(values)) / float(
                len(values)
            )

        begin = dominance([item[0] for item in spacers])
        end = dominance([item[-1] for item in spacers])
        if len(spacers) < 5:
            status = "insufficient_sample"
            reason = "terminal_dominance_not_authorized_at_small_n"
        elif begin > 0.8 or end > 0.8:
            status = "failed"
            reason = "terminal_base_dominance_above_threshold"
        else:
            status = "passed"
            reason = "terminal_base_dominance_within_threshold"
    expected = {
        "status": status,
        "reason": reason,
        "nonempty_spacer_count": len(spacers),
        "minimum_sample_size": 5,
        "maximum_allowed_dominance": 0.8,
        "begin_dominance": begin,
        "end_dominance": end,
    }
    if _canonical_json(terminal) != _canonical_json(expected):
        raise ReportProvenanceError(label + " is not source-derived")


def _validate_structural(source, geometry, mapping, label):
    expected_keys = frozenset((
        "decision",
        "structural_valid",
        "predicate_results",
        "terminal_dominance",
    ))
    _keys(mapping, expected_keys, label)
    _bool(mapping["structural_valid"], label + ".structural_valid")
    predicates = _sequence(
        mapping["predicate_results"], label + ".predicate_results"
    )
    names = []
    for index, predicate in enumerate(predicates):
        item_label = "{}.predicate_results[{}]".format(label, index)
        _keys(predicate, frozenset(("predicate", "status")), item_label)
        if not isinstance(predicate["predicate"], str):
            raise ReportProvenanceError(item_label + ".predicate must be text")
        if predicate["status"] not in (
            "passed", "failed", "insufficient_sample"
        ):
            raise ReportProvenanceError(item_label + ".status is invalid")
        names.append(predicate["predicate"])
    if len(names) != len(set(names)):
        raise ReportProvenanceError(label + " contains duplicate predicates")
    terminal = mapping["terminal_dominance"]
    if terminal is None:
        if mapping["structural_valid"]:
            raise ReportProvenanceError(
                label + " omits terminal evidence for a valid structure"
            )
    else:
        _validate_terminal(
            source,
            geometry,
            _mapping(terminal, label + ".terminal_dominance"),
            label + ".terminal_dominance",
        )
        terminal_status = terminal["status"]
        named_terminal = [
            item["status"] for item in predicates
            if item["predicate"] == "spacer_terminal_dominance"
        ]
        if named_terminal != [terminal_status]:
            raise ReportProvenanceError(
                label + " terminal predicate contradicts its evidence"
            )
    expected_pass = (
        mapping["structural_valid"]
        and terminal is not None
        and terminal["status"] in ("passed", "insufficient_sample")
        and all(
            item["status"] == "passed"
            for item in predicates
            if item["predicate"] != "spacer_terminal_dominance"
        )
    )
    if mapping["decision"] != ("passed" if expected_pass else "rejected"):
        raise ReportProvenanceError(label + " decision is internally stale")
    return mapping


def _authorization_fraction(value):
    return (
        isinstance(value, Real)
        and not isinstance(value, bool)
        and math.isfinite(float(value))
        and 0.0 <= float(value) <= 1.0
    )


def _validate_terminal_authorization(
    mapping,
    existence,
    structural,
    label,
):
    """Replay the explicit ordinary-or-four-spacer authorization policy."""

    _keys(mapping, _TERMINAL_AUTHORIZATION_KEYS, label)
    rule = SOURCE_INTACT_TERMINAL_AUTHORIZATION_RULE
    if mapping["rule"] != rule.as_dict():
        raise ReportProvenanceError(label + " uses a stale or altered rule")

    terminal = structural["terminal_dominance"]
    if not isinstance(terminal, dict):
        raise ReportProvenanceError(label + " lacks terminal evidence")
    metrics = existence["metrics"]
    reasons = []
    existence_eligible = existence["eligible"] is True
    structural_passed = structural["decision"] == "passed"
    if not existence_eligible:
        reasons.append("source_intact_existence_not_eligible")
    if not structural_passed:
        reasons.append("structural_filter_not_passed")

    status = terminal["status"]
    if status == "passed":
        authorized = not reasons
        mode = STATISTICALLY_PASSED if authorized else TERMINAL_AUTHORIZATION_REJECTED
    else:
        if status != "insufficient_sample":
            reasons.append("terminal_status_not_authorized")
        if terminal["nonempty_spacer_count"] != rule.small_sample_spacer_count:
            reasons.append("small_sample_terminal_spacer_count_differs")
        if (
            terminal["minimum_sample_size"]
            != rule.structural_minimum_terminal_sample_size
        ):
            reasons.append("terminal_statistical_minimum_differs")
        if metrics["repeat_count"] != rule.small_sample_repeat_count:
            reasons.append("small_sample_repeat_count_differs")
        if (
            metrics["supplied_spacer_count"]
            != rule.small_sample_spacer_count
        ):
            reasons.append("small_sample_supplied_spacer_count_differs")
        if (
            metrics["nonempty_supplied_spacer_count"]
            != rule.small_sample_spacer_count
        ):
            reasons.append("small_sample_nonempty_spacer_count_differs")
        if metrics["distinct_spacer_count"] != rule.small_sample_spacer_count:
            reasons.append("small_sample_spacers_not_all_distinct")
        if metrics["exact_duplicate_spacer_pair_count"] != 0:
            reasons.append("small_sample_duplicate_spacer_pair_present")
        for side in ("begin", "end"):
            value = terminal[side + "_dominance"]
            if not _authorization_fraction(value):
                reasons.append(
                    "small_sample_{}_dominance_invalid".format(side)
                )
            elif value > rule.maximum_small_sample_edge_dominance:
                reasons.append(
                    "small_sample_{}_dominance_above_maximum".format(side)
                )
        authorized = not reasons
        mode = GUARDED_FOUR_SPACER if authorized else TERMINAL_AUTHORIZATION_REJECTED

    expected = {
        "authorized": authorized,
        "mode": mode,
        "reason_codes": sorted(set(reasons)),
        "observed_terminal_status": status,
        "nonempty_spacer_count": terminal["nonempty_spacer_count"],
        "begin_dominance": terminal["begin_dominance"],
        "end_dominance": terminal["end_dominance"],
        "source_intact_existence_eligible": existence_eligible,
        "structural_filter_passed": structural_passed,
        "rule": rule.as_dict(),
    }
    if _canonical_json(mapping) != _canonical_json(expected):
        raise ReportProvenanceError(
            label + " does not equal independently replayed authorization"
        )
    return mapping


def _tandem_rule():
    return {
        "name": tandem_module.RULE_NAME,
        "version": tandem_module.RULE_VERSION,
        "minimum_family_envelope_coverage": (
            tandem_module.MINIMUM_FAMILY_ENVELOPE_COVERAGE
        ),
        "minimum_cds_envelope_coverage": (
            tandem_module.MINIMUM_CDS_ENVELOPE_COVERAGE
        ),
        "minimum_cds_confidence": tandem_module.MINIMUM_CDS_CONFIDENCE,
        "minimum_adjacent_period_count": (
            tandem_module.MINIMUM_ADJACENT_PERIOD_COUNT
        ),
        "minimum_period": tandem_module.MINIMUM_PERIOD,
        "maximum_period": tandem_module.MAXIMUM_PERIOD,
        "period_must_be_divisible_by_three": True,
        "minimum_adjacent_peptide_identity": (
            tandem_module.MINIMUM_ADJACENT_PEPTIDE_IDENTITY
        ),
        "minimum_mean_adjacent_peptide_identity": (
            tandem_module.MINIMUM_MEAN_ADJACENT_PEPTIDE_IDENTITY
        ),
        "minimum_mean_all_pair_peptide_identity": (
            tandem_module.MINIMUM_MEAN_ALL_PAIR_PEPTIDE_IDENTITY
        ),
        "minimum_mean_adjacent_dna_identity": (
            tandem_module.MINIMUM_MEAN_ADJACENT_DNA_IDENTITY
        ),
        "translation_scope": "selected_source_cds_strand_and_frame_only",
        "translation_must_be_stop_and_x_free": True,
        "all_full_period_units_must_be_inside_selected_cds": True,
        "classifier_features_used": False,
        "six_frame_maximization_used": False,
    }


_TANDEM_KEYS = frozenset((
    "supported",
    "decision",
    "reason",
    "rule_name",
    "rule_version",
    "coordinate_system",
    "accepted_interval",
    "family_evidence_interval",
    "family_envelope_coverage",
    "selected_cds_interval",
    "selected_cds_strand",
    "selected_cds_confidence",
    "selected_cds_envelope_coverage",
    "period",
    "period_unit_intervals",
    "translation_frame_offset",
    "peptide_lengths",
    "adjacent_peptide_identities",
    "minimum_adjacent_peptide_identity",
    "mean_adjacent_peptide_identity",
    "mean_all_pair_peptide_identity",
    "adjacent_dna_identities",
    "mean_adjacent_dna_identity",
    "family_unique_path_count",
    "source_sha256",
    "accepted_source_slice_sha256",
    "family_source_slice_sha256",
    "selected_cds_source_slice_sha256",
    "period_unit_source_slice_sha256s",
    "rule",
))


def _optional_slice_hash(source, interval, digest, label):
    if interval is None:
        if digest is not None:
            raise ReportProvenanceError(label + " hash exists without interval")
        return
    if _hash(digest, label + ".sha256") != _sha256(
        source[interval[0]:interval[1]]
    ):
        raise ReportProvenanceError(label + " source hash is stale")


def _validate_tandem(
    source,
    source_hash,
    mapping,
    expected_geometry,
    label,
    family_geometries=(),
):
    _keys(mapping, _TANDEM_KEYS, label)
    if mapping["rule"] != _tandem_rule():
        raise ReportProvenanceError(label + " fixed rule is stale")
    if (
        mapping["rule_name"] != tandem_module.RULE_NAME
        or mapping["rule_version"] != tandem_module.RULE_VERSION
        or mapping["coordinate_system"] != _COORDINATE_SYSTEM
        or mapping["source_sha256"] != source_hash
    ):
        raise ReportProvenanceError(label + " rule/source binding is stale")
    accepted = _interval(
        mapping["accepted_interval"], len(source), label + ".accepted_interval"
    )
    if accepted != expected_geometry.interval:
        raise ReportProvenanceError(label + " accepted interval is stale")
    if (
        _hash(
            mapping["accepted_source_slice_sha256"],
            label + ".accepted_source_slice_sha256",
        )
        != _sha256(source[accepted[0]:accepted[1]])
    ):
        raise ReportProvenanceError(label + " accepted slice hash is stale")

    family = _interval(
        mapping["family_evidence_interval"],
        len(source),
        label + ".family_evidence_interval",
        allow_none=True,
    )
    _optional_slice_hash(
        source,
        family,
        mapping["family_source_slice_sha256"],
        label + ".family_evidence_interval",
    )
    if family is not None and family_geometries and all(
        family != geometry.interval for geometry in family_geometries
    ):
        raise ReportProvenanceError(
            label + " family interval is absent from complete family geometry"
        )
    cds = _interval(
        mapping["selected_cds_interval"],
        len(source),
        label + ".selected_cds_interval",
        allow_none=True,
    )
    _optional_slice_hash(
        source,
        cds,
        mapping["selected_cds_source_slice_sha256"],
        label + ".selected_cds_interval",
    )
    units = _intervals(
        mapping["period_unit_intervals"],
        len(source),
        label + ".period_unit_intervals",
    )
    unit_hashes = tuple(
        _hash(item, label + ".period_unit_source_slice_sha256s")
        for item in _sequence(
            mapping["period_unit_source_slice_sha256s"],
            label + ".period_unit_source_slice_sha256s",
        )
    )
    if unit_hashes != tuple(_sha256(source[a:b]) for a, b in units):
        raise ReportProvenanceError(label + " period-unit hashes are stale")
    if cds is None and units:
        raise ReportProvenanceError(label + " period units lack a source CDS")
    if cds is not None and any(
        not (cds[0] <= start < end <= cds[1]) for start, end in units
    ):
        raise ReportProvenanceError(label + " period unit escapes selected CDS")
    period = mapping["period"]
    if period is not None:
        period = _integer(period, label + ".period", 1)
        if period % 3 or any(end - start != period for start, end in units):
            raise ReportProvenanceError(label + " period/unit lengths disagree")
    elif units:
        raise ReportProvenanceError(label + " period units lack a period")
    family_unique_path_count = _integer(
        mapping["family_unique_path_count"],
        label + ".family_unique_path_count",
        1,
    )
    if family_geometries and family_unique_path_count != len({
        geometry.geometry_id for geometry in family_geometries
    }):
        raise ReportProvenanceError(
            label + " family path count is stale"
        )
    _bool(mapping["supported"], label + ".supported")
    return mapping


def _validate_snapshot(source, mapping, label):
    return _validate_geometry(source, mapping, label, snapshot=True)


def _validate_architecture(
    source,
    source_hash,
    mapping,
    geometry,
    label,
    validation_context,
):
    _keys(mapping, _ARCHITECTURE_KEYS, label)
    existence = _validate_existence(
        source,
        source_hash,
        _mapping(mapping["strict_existence"], label + ".strict_existence"),
        geometry,
        label + ".strict_existence",
        validation_context,
    )
    structural = _validate_structural(
        source,
        geometry,
        _mapping(mapping["structural_filter"], label + ".structural_filter"),
        label + ".structural_filter",
    )
    if mapping["strict_existence_eligible"] != existence["eligible"]:
        raise ReportProvenanceError(label + " existence eligibility is stale")
    if mapping["strict_existence_rejection_codes"] != existence[
        "rejection_codes"
    ]:
        raise ReportProvenanceError(label + " rejection codes are stale")
    if mapping["structural_filter_passed"] != (
        structural["decision"] == "passed"
    ):
        raise ReportProvenanceError(label + " structural decision is stale")
    terminal = structural["terminal_dominance"]
    expected_status = None if terminal is None else terminal["status"]
    if mapping["structural_filter_terminal_status"] != expected_status:
        raise ReportProvenanceError(label + " terminal status is stale")
    _sequence(mapping["weakness_domains"], label + ".weakness_domains")
    return mapping


def _validate_conflict(
    source,
    source_hash,
    mapping,
    recovery_geometry,
    targets,
    label,
    validation_context,
):
    _keys(mapping, _CONFLICT_KEYS, label)
    if (
        mapping["source_sha256"] != source_hash
        or mapping["coordinate_system"] != _COORDINATE_SYSTEM
        or _canonical_json(mapping["policy"])
        != _canonical_json(ACCEPTED_SOURCE_CONFLICT_POLICY.as_dict())
    ):
        raise ReportProvenanceError(label + " source/rule binding is stale")
    accepted = _validate_snapshot(
        source,
        _mapping(mapping["accepted_snapshot"], label + ".accepted_snapshot"),
        label + ".accepted_snapshot",
    )
    recovery = _validate_snapshot(
        source,
        _mapping(
            mapping["source_intact_snapshot"],
            label + ".source_intact_snapshot",
        ),
        label + ".source_intact_snapshot",
    )
    if (
        recovery.interval != recovery_geometry.interval
        or recovery.repeat_intervals != recovery_geometry.repeat_intervals
    ):
        raise ReportProvenanceError(label + " recovery snapshot is stale")
    if not any(
        accepted.interval == target.interval
        and accepted.repeat_intervals == target.repeat_intervals
        for target in targets
    ):
        raise ReportProvenanceError(
            label + " accepted snapshot has no accepted target"
        )
    _validate_architecture(
        source,
        source_hash,
        _mapping(
            mapping["accepted_architecture"],
            label + ".accepted_architecture",
        ),
        accepted,
        label + ".accepted_architecture",
        validation_context,
    )
    _validate_architecture(
        source,
        source_hash,
        _mapping(
            mapping["source_intact_architecture"],
            label + ".source_intact_architecture",
        ),
        recovery,
        label + ".source_intact_architecture",
        validation_context,
    )
    claimed_relationship = _mapping(
        mapping["relationship"], label + ".relationship"
    )
    _keys(claimed_relationship, _RELATIONSHIP_KEYS, label + ".relationship")
    expected_relationship = recompute_conflict_relationship(
        SimpleNamespace(
            interval=accepted.interval,
            repeat_intervals=accepted.repeat_intervals,
        ),
        SimpleNamespace(
            interval=recovery.interval,
            repeat_intervals=recovery.repeat_intervals,
        ),
    ).as_dict()
    if _canonical_json(claimed_relationship) != _canonical_json(
        expected_relationship
    ):
        raise ReportProvenanceError(label + " conflict topology is stale")
    for field_name, geometry in (
        ("accepted_tandem", accepted),
        ("source_intact_tandem", recovery),
    ):
        value = mapping[field_name]
        if value is not None:
            _validate_tandem(
                source,
                source_hash,
                _mapping(value, label + "." + field_name),
                geometry,
                label + "." + field_name,
            )
    _sequence(mapping["invalid_evidence_codes"], label + ".invalid_evidence_codes")
    for key in (
        "source_authorized",
        "classifier_score_compared",
        "candidate_length_preference_used",
        "exact_boundary_certified",
    ):
        _bool(mapping[key], label + "." + key)
    return accepted, claimed_relationship


def _relationship_kind(first, second):
    if first.repeat_intervals == second.repeat_intervals:
        return "identical_geometry"
    first_repeats = set(first.repeat_intervals)
    second_repeats = set(second.repeat_intervals)
    if first_repeats.issubset(second_repeats):
        return "first_exact_repeat_subpath"
    if second_repeats.issubset(first_repeats):
        return "second_exact_repeat_subpath"
    return "overlapping_boundary_or_phase_alternatives"


def _overlap(first, second):
    start = max(first[0], second[0])
    end = min(first[1], second[1])
    return None if end <= start else (start, end)


def _validate_alternative(source, source_hash, mapping, label):
    _keys(mapping, _ALTERNATIVE_KEYS, label)
    if mapping["source_sha256"] != source_hash:
        raise ReportProvenanceError(label + " source SHA-256 is stale")
    first = _validate_geometry(
        source,
        _mapping(mapping["first_geometry"], label + ".first_geometry"),
        label + ".first_geometry",
    )
    second = _validate_geometry(
        source,
        _mapping(mapping["second_geometry"], label + ".second_geometry"),
        label + ".second_geometry",
    )
    if first.geometry_id >= second.geometry_id:
        raise ReportProvenanceError(label + " geometries are not canonical")
    overlap = _overlap(first.interval, second.interval)
    if overlap is None or tuple(mapping["overlap_interval"]) != overlap:
        raise ReportProvenanceError(label + " overlap is stale or disjoint")
    if mapping["relationship"] != _relationship_kind(first, second):
        raise ReportProvenanceError(label + " relationship is stale")
    if (
        mapping["boundary_status"] != "AMBIGUOUS"
        or mapping["representative_selected"] is not False
        or mapping["classifier_score_compared"] is not False
        or mapping["candidate_length_preference_used"] is not False
    ):
        raise ReportProvenanceError(label + " overstates alternative selection")
    return first, second


def _validate_recovery_record(
    source,
    source_hash,
    mapping,
    label,
    validation_context,
):
    _keys(mapping, _RECORD_KEYS, label)
    if mapping["rule"] != _SOURCE_RECORD_RULE:
        raise ReportProvenanceError(label + " rule is stale")
    for key in (
        "exact_boundary_certified",
        "classifier_score_used",
        "candidate_length_preference_used",
    ):
        if mapping[key] is not False:
            raise ReportProvenanceError(label + "." + key + " must be false")
    if mapping["boundary_status"] != "AMBIGUOUS":
        raise ReportProvenanceError(label + " must retain ambiguous boundaries")
    action = mapping["action"]
    if action not in (
        "admit_possible",
        "replace_conflicting_possible",
        "attach_ambiguous_alternative",
    ):
        raise ReportProvenanceError(label + " action is not reportable")
    source_category_index = _integer(
        mapping["source_category_index"],
        label + ".source_category_index",
    )
    if (
        source_category_index not in _RECOVERY_SOURCE_CATEGORY_INDEXES
        or mapping["source_category_name"]
        != _CATEGORY_NAMES[source_category_index]
    ):
        raise ReportProvenanceError(
            label + " source category index/name is inconsistent"
        )
    if action in ("admit_possible", "replace_conflicting_possible"):
        expected_target = (2, _CATEGORY_NAMES[2])
    else:
        expected_target = (None, None)
    if (
        mapping["target_category_index"],
        mapping["target_category_name"],
    ) != expected_target:
        raise ReportProvenanceError(
            label + " action and target category are inconsistent"
        )
    _finite(mapping["score"], label + ".score")
    _integer(mapping["source_row_ordinal"], label + ".source_row_ordinal")
    _sequence(mapping["reason_codes"], label + ".reason_codes")
    _sequence(mapping["accepted_targets"], label + ".accepted_targets")
    _sequence(mapping["conflict_decisions"], label + ".conflict_decisions")
    _sequence(
        mapping["alternative_relationships"],
        label + ".alternative_relationships",
    )

    evidence = _mapping(mapping["evidence"], label + ".evidence")
    _keys(evidence, _EVIDENCE_KEYS, label + ".evidence")
    if evidence["source_sha256"] != source_hash:
        raise ReportProvenanceError(label + ".evidence source hash is stale")
    geometry = _validate_geometry(
        source,
        _mapping(evidence["geometry"], label + ".evidence.geometry"),
        label + ".evidence.geometry",
    )
    family = tuple(
        _validate_geometry(
            source,
            _mapping(item, "{}.evidence.family_geometries[{}]".format(
                label, index
            )),
            "{}.evidence.family_geometries[{}]".format(label, index),
        )
        for index, item in enumerate(_sequence(
            evidence["family_geometries"],
            label + ".evidence.family_geometries",
        ))
    )
    if not family or not any(_same_geometry(geometry, item) for item in family):
        raise ReportProvenanceError(
            label + ".evidence complete family omits the candidate geometry"
        )
    family_order = tuple((item.interval, item.geometry_id) for item in family)
    if family_order != tuple(sorted(family_order)):
        raise ReportProvenanceError(
            label + ".evidence family geometries are not canonical"
        )
    family_ids = tuple(item.geometry_id for item in family)
    if (
        len(family_ids) != len(set(family_ids))
        or sum(_same_geometry(geometry, item) for item in family) != 1
    ):
        raise ReportProvenanceError(
            label + ".evidence family geometries are not unique"
        )
    existence = _validate_existence(
        source,
        source_hash,
        _mapping(evidence["existence"], label + ".evidence.existence"),
        geometry,
        label + ".evidence.existence",
        validation_context,
    )
    structural = _validate_structural(
        source,
        geometry,
        _mapping(
            evidence["structural_filter"],
            label + ".evidence.structural_filter",
        ),
        label + ".evidence.structural_filter",
    )
    terminal_authorization = _validate_terminal_authorization(
        _mapping(
            evidence["terminal_authorization"],
            label + ".evidence.terminal_authorization",
        ),
        existence,
        structural,
        label + ".evidence.terminal_authorization",
    )
    supplied = _validate_tandem(
        source,
        source_hash,
        _mapping(evidence["supplied_tandem"], label + ".evidence.supplied_tandem"),
        geometry,
        label + ".evidence.supplied_tandem",
        family,
    )
    recomputed = _validate_tandem(
        source,
        source_hash,
        _mapping(
            evidence["recomputed_tandem"],
            label + ".evidence.recomputed_tandem",
        ),
        geometry,
        label + ".evidence.recomputed_tandem",
        family,
    )
    if _canonical_json(supplied) != _canonical_json(recomputed):
        raise ReportProvenanceError(label + " tandem evidence was not recomputed")
    if (
        evidence["authorized"] is not True
        or evidence["invalid_evidence_codes"] != []
        or evidence["terminal_status_required"]
        != "explicit_policy_authorization"
        or evidence["boundary_status"] != "AMBIGUOUS"
        or evidence["exact_boundary_certified"] is not False
        or evidence["classifier_score_used"] is not False
        or evidence["candidate_length_preference_used"] is not False
        or evidence["existence"]["eligible"] is not True
        or structural["decision"] != "passed"
        or terminal_authorization["authorized"] is not True
        or recomputed["supported"] is not False
    ):
        raise ReportProvenanceError(
            label + " is not a complete authorized recovery record"
        )

    targets = tuple(
        _validate_target(source, item, "{}.accepted_targets[{}]".format(
            label, index
        ))
        for index, item in enumerate(mapping["accepted_targets"])
    )
    target_order = tuple(
        (
            item["category_index"],
            _canonical_json(item["locus_key"]),
            item["row_ordinal"],
        )
        for item in mapping["accepted_targets"]
    )
    if target_order != tuple(sorted(target_order)):
        raise ReportProvenanceError(
            label + " accepted targets are not canonical"
        )
    if len(target_order) != len(set(target_order)):
        raise ReportProvenanceError(
            label + " accepted targets are duplicated"
        )
    conflicts = tuple(
        _validate_conflict(
            source,
            source_hash,
            _mapping(item, "{}.conflict_decisions[{}]".format(label, index)),
            geometry,
            targets,
            "{}.conflict_decisions[{}]".format(label, index),
            validation_context,
        )
        for index, item in enumerate(mapping["conflict_decisions"])
    )
    alternatives = tuple(
        _validate_alternative(
            source,
            source_hash,
            _mapping(item, "{}.alternative_relationships[{}]".format(
                label, index
            )),
            "{}.alternative_relationships[{}]".format(label, index),
        )
        for index, item in enumerate(mapping["alternative_relationships"])
    )
    return geometry, targets, conflicts, alternatives


def _validate_target(source, mapping, label):
    _keys(mapping, _TARGET_KEYS, label)
    category_index = _integer(
        mapping["category_index"], label + ".category_index"
    )
    if (
        category_index not in _ACCEPTED_CATEGORY_INDEXES
        or mapping["category_name"] != _CATEGORY_NAMES[category_index]
    ):
        raise ReportProvenanceError(label + " category index/name disagree")
    _integer(mapping["row_ordinal"], label + ".row_ordinal")
    return _validate_geometry(
        source,
        _mapping(mapping["geometry"], label + ".geometry"),
        label + ".geometry",
    )


def _boundary_required_mapping(mapping, key, label):
    mapping = _mapping(mapping, label)
    if key not in mapping:
        raise ReportProvenanceError(
            "{} is missing {!r}".format(label, key)
        )
    return _mapping(mapping[key], "{}.{}".format(label, key))


def _boundary_required_sequence(mapping, key, label):
    mapping = _mapping(mapping, label)
    if key not in mapping:
        raise ReportProvenanceError(
            "{} is missing {!r}".format(label, key)
        )
    return _sequence(mapping[key], "{}.{}".format(label, key))


def _boundary_internal_interval(mapping, label):
    internal = _boundary_required_mapping(
        mapping,
        "internal_0based_half_open",
        label,
    )
    if "start" not in internal or "end" not in internal:
        raise ReportProvenanceError(
            label + ".internal_0based_half_open lacks start or end"
        )
    return internal["start"], internal["end"]


def _boundary_geometry_input(mapping, label):
    mapping = _mapping(mapping, label)
    repeats = _boundary_required_sequence(mapping, "repeats", label)
    spacers = _boundary_required_sequence(mapping, "spacers", label)

    def component_interval(component, component_label):
        component = _mapping(component, component_label)
        return _boundary_internal_interval(
            _boundary_required_mapping(
                component,
                "interval",
                component_label,
            ),
            component_label + ".interval",
        )

    return CandidateGeometryInput(
        source_sha256=mapping.get("source_sha256"),
        source_interval=_boundary_internal_interval(
            _boundary_required_mapping(mapping, "source_interval", label),
            label + ".source_interval",
        ),
        repeat_intervals=tuple(
            component_interval(
                component,
                "{}.repeats[{}]".format(label, index),
            )
            for index, component in enumerate(repeats)
        ),
        repeat_sequences=tuple(
            _mapping(
                component,
                "{}.repeats[{}]".format(label, index),
            ).get("sequence")
            for index, component in enumerate(repeats)
        ),
        spacer_intervals=tuple(
            component_interval(
                component,
                "{}.spacers[{}]".format(label, index),
            )
            for index, component in enumerate(spacers)
        ),
        spacer_sequences=tuple(
            _mapping(
                component,
                "{}.spacers[{}]".format(label, index),
            ).get("sequence")
            for index, component in enumerate(spacers)
        ),
    )


def _boundary_dimension_input(mapping, dimension, label):
    dimensions = _boundary_required_mapping(mapping, "dimensions", label)
    value = _boundary_required_mapping(
        dimensions,
        dimension,
        label + ".dimensions",
    )
    return BoundaryDimensionInput(
        status=value.get("status"),
        selection_policy_id=value.get("selection_policy_id"),
        reason_codes=_boundary_required_sequence(
            value,
            "reason_codes",
            "{}.dimensions.{}".format(label, dimension),
        ),
    )


def _revalidate_boundary_claim(
    source,
    source_hash,
    candidate,
    raw_mapping,
    label,
):
    """Rehydrate and independently normalize one canonical boundary claim."""

    try:
        mapping = _mapping(raw_mapping, label)
        source_claim = _boundary_required_mapping(mapping, "source", label)
        if (
            source_claim.get("sha256") != source_hash
            or source_claim.get("length") != len(source)
        ):
            raise ReportProvenanceError(label + " source binding is stale")
        reported = _boundary_required_mapping(
            mapping,
            "reported_geometry",
            label,
        )
        relationships = []
        for index, raw_relationship in enumerate(
            _boundary_required_sequence(
                mapping,
                "alternative_relationships",
                label,
            )
        ):
            relationship_label = (
                "{}.alternative_relationships[{}]".format(label, index)
            )
            raw_relationship = _mapping(
                raw_relationship,
                relationship_label,
            )
            relationships.append(AlternativeRelationshipInput(
                geometry=_boundary_geometry_input(
                    _boundary_required_mapping(
                        raw_relationship,
                        "geometry",
                        relationship_label,
                    ),
                    relationship_label + ".geometry",
                ),
                label=raw_relationship.get("label"),
                dimensions=_boundary_required_sequence(
                    raw_relationship,
                    "dimensions",
                    relationship_label,
                ),
                reason_codes=_boundary_required_sequence(
                    raw_relationship,
                    "reason_codes",
                    relationship_label,
                ),
            ))
        normalizer = BoundaryAssessmentNormalizer(source, source_hash)
        normalized = normalizer.normalize(
            _boundary_geometry_input(reported, label + ".reported_geometry"),
            alternative_relationships=tuple(relationships),
            occurrence_path=_boundary_dimension_input(
                mapping,
                DIMENSION_OCCURRENCE_PATH,
                label,
            ),
            repeat_phase=_boundary_dimension_input(
                mapping,
                DIMENSION_REPEAT_PHASE,
                label,
            ),
        )
        expected_mapping = normalized.as_dict()
        if _canonical_json(mapping) != _canonical_json(expected_mapping):
            raise ReportProvenanceError(
                label + " is not its canonical source reconstruction"
            )
        candidate_geometry = normalizer.normalize(
            geometry_input_from_candidate(candidate, source_hash)
        ).reported_geometry.as_dict()
        if expected_mapping["reported_geometry"] != candidate_geometry:
            raise ReportProvenanceError(
                label + " representative differs from its candidate"
            )
        return expected_mapping
    except ReportProvenanceError:
        raise
    except (AttributeError, IndexError, KeyError, TypeError, ValueError) as error:
        raise ReportProvenanceError(
            "{} cannot be source-revalidated: {}".format(label, error)
        ) from error


def _architecture_nonempty_string(value, label):
    if not isinstance(value, str) or not value:
        raise ReportProvenanceError(
            "{} must be a non-empty string".format(label)
        )
    return value


def _architecture_path_id_sequence(value, label):
    values = _sequence(value, label)
    if (
            any(
                not isinstance(item, str) or not item
                for item in values
            )
            or tuple(sorted(set(values))) != values):
        raise ReportProvenanceError(
            "{} must be sorted unique path IDs".format(label)
        )
    return values


def _architecture_geometry_input_from_intervals(
        source, source_hash, raw_intervals, label):
    intervals = _intervals(raw_intervals, len(source), label)
    if not intervals:
        raise ReportProvenanceError(label + " must contain repeats")
    if any(
            right[0] < left[1]
            for left, right in zip(intervals, intervals[1:])):
        raise ReportProvenanceError(
            label + " repeat intervals overlap or are nonmonotonic"
        )
    spacer_intervals = tuple(
        (left[1], right[0])
        for left, right in zip(intervals, intervals[1:])
    )
    return CandidateGeometryInput(
        source_sha256=source_hash,
        source_interval=(intervals[0][0], intervals[-1][1]),
        repeat_intervals=intervals,
        repeat_sequences=tuple(source[start:end] for start, end in intervals),
        spacer_intervals=spacer_intervals,
        spacer_sequences=tuple(
            source[start:end] for start, end in spacer_intervals
        ),
    )


def _architecture_exact_contiguous_subpath(shorter, longer):
    shorter = tuple(shorter)
    longer = tuple(longer)
    if not shorter or len(shorter) > len(longer):
        return False
    return any(
        longer[offset:offset + len(shorter)] == shorter
        for offset in range(len(longer) - len(shorter) + 1)
    )


def _architecture_phase_alternative_geometries(
        source, source_hash, candidate_geometry, phase_arbitration,
        expected_path_id, expected_superlocus_id):
    if phase_arbitration is None:
        return (), None
    label = "classification_basis.source_phase_arbitration"
    phase = _mapping(phase_arbitration, label)
    _keys(phase, _ARCHITECTURE_PHASE_KEYS, label)
    if (
            phase["coordinate_system"] != _COORDINATE_SYSTEM
            or phase["reported_geometry_role"] != "representative"
            or phase["superlocus_id"] != expected_superlocus_id):
        raise ReportProvenanceError(
            label + " source/superlocus geometry contract is stale"
        )
    _mapping(phase["policy"], label + ".policy")
    _architecture_nonempty_string(phase["reason"], label + ".reason")
    retained = _architecture_path_id_sequence(
        phase["retained_path_ids"], label + ".retained_path_ids"
    )
    removed = _architecture_path_id_sequence(
        phase["removed_path_ids"], label + ".removed_path_ids"
    )
    competitors = _architecture_path_id_sequence(
        phase["competitor_path_ids"], label + ".competitor_path_ids"
    )
    if (
            not retained
            or expected_path_id not in retained
            or set(retained).intersection(removed)
            or set(retained).union(removed) != set(competitors)
            or phase["representative_path_id"] != min(retained)):
        raise ReportProvenanceError(
            label + " retained/removed path partition is stale"
        )
    status = phase["status"]
    if status == "selected":
        if (
                phase["boundary_status"] != STATUS_POLICY_SELECTED
                or phase["winner_path_id"] not in retained):
            raise ReportProvenanceError(
                label + " selected phase fields are inconsistent"
            )
    elif status == "ambiguous_representative":
        if (
                phase["boundary_status"] != STATUS_AMBIGUOUS
                or phase["winner_path_id"] is not None):
            raise ReportProvenanceError(
                label + " ambiguous phase fields are inconsistent"
            )
    else:
        raise ReportProvenanceError(label + " status is unsupported")

    source_evidence = _mapping(
        phase["source_evidence"], label + ".source_evidence"
    )
    _keys(source_evidence, frozenset((
        "classifier_used",
        "winner_margin",
        "background_window",
        "background_probabilities",
        "rule",
        "candidates",
    )), label + ".source_evidence")
    _mapping(
        source_evidence["background_probabilities"],
        label + ".source_evidence.background_probabilities",
    )
    _mapping(source_evidence["rule"], label + ".source_evidence.rule")
    geometries = []
    for index, raw_candidate in enumerate(_sequence(
            source_evidence["candidates"],
            label + ".source_evidence.candidates")):
        candidate_label = "{}.source_evidence.candidates[{}]".format(
            label, index
        )
        record = _mapping(raw_candidate, candidate_label)
        _keys(record, frozenset((
            "geometry", "source_score", "boundary_evidence",
        )), candidate_label)
        if record["source_score"] is None:
            if not (
                    status == "ambiguous_representative"
                    and phase["reason"]
                    == "geometry_outside_source_phase_rule"):
                raise ReportProvenanceError(
                    candidate_label + ".source_score is unexpectedly null"
                )
        else:
            _finite(
                record["source_score"], candidate_label + ".source_score"
            )
        _mapping(
            record["boundary_evidence"],
            candidate_label + ".boundary_evidence",
        )
        geometry_input = _architecture_geometry_input_from_intervals(
            source,
            source_hash,
            record["geometry"],
            candidate_label + ".geometry",
        )
        geometries.append(tuple(geometry_input.repeat_intervals))
    if (
            len(geometries) < 2
            or len(set(geometries)) != len(geometries)
            or candidate_geometry not in geometries):
        raise ReportProvenanceError(
            label + " lacks distinct explicit source phase alternatives"
        )
    return tuple(
        geometry for geometry in geometries
        if geometry != candidate_geometry
    ), status


def canonical_source_architecture_rescue_boundary_assessment(
        source_sequence, expected_source_sha256, candidate,
        replayed_path_local_proof, source_phase_arbitration=None, *,
        expected_path_id, expected_superlocus_id, source_normalizer=None):
    """Build the one canonical boundary record for a local envelope.

    The caller must first replay ``replayed_path_local_proof``.  Every bounded
    and confirmation geometry below is nevertheless reconstructed from the
    original source by :class:`BoundaryAssessmentNormalizer`.  Occurrence
    membership is policy-selected among the complete distinct witnessed
    subpaths.  Repeat phase remains unassessed unless the existing cross-family
    arbitration record supplies at least one explicit source-valid phase.
    """

    try:
        source_hash = _hash(
            expected_source_sha256, "expected_source_sha256"
        )
        proof = _mapping(
            replayed_path_local_proof, "path_local_convergence_proof"
        )
        if (
                proof.get("status") != "converged"
                or proof.get("reason")
                != "finite_contained_terminal_closure"
                or proof.get("superlocus_id") != expected_superlocus_id
                or proof.get("envelope_path_id") != expected_path_id):
            raise ReportProvenanceError(
                "path-local proof does not certify the reported envelope"
            )
        if source_normalizer is None:
            source_normalizer = BoundaryAssessmentNormalizer(
                source_sequence, source_hash
            )
        else:
            source_normalizer.assert_bound_to(
                source_sequence, source_hash
            )
        representative_input = geometry_input_from_candidate(
            candidate, source_hash
        )
        representative_geometry = tuple(
            representative_input.repeat_intervals
        )
        envelope_geometry = tuple(
            _intervals(
                proof.get("envelope_geometry"),
                len(source_sequence),
                "path_local_convergence_proof.envelope_geometry",
            )
        )
        if representative_geometry != envelope_geometry:
            raise ReportProvenanceError(
                "path-local envelope differs from its reported candidate"
            )

        descriptors = {}
        witness_memberships = (
            ("bounded_paths", "path_local_bounded_witness"),
            ("confirmation_paths", "path_local_confirmation_witness"),
        )
        for witness_key, reason_code in witness_memberships:
            for index, raw_record in enumerate(_sequence(
                    proof.get(witness_key),
                    "path_local_convergence_proof." + witness_key)):
                record_label = (
                    "path_local_convergence_proof.{}[{}]".format(
                        witness_key, index
                    )
                )
                record = _mapping(raw_record, record_label)
                geometry = tuple(_intervals(
                    record.get("source_geometry"),
                    len(source_sequence),
                    record_label + ".source_geometry",
                ))
                if geometry == representative_geometry:
                    continue
                if not _architecture_exact_contiguous_subpath(
                        geometry, representative_geometry):
                    raise ReportProvenanceError(
                        record_label + " is not dominated by the envelope"
                    )
                descriptor = descriptors.setdefault(geometry, {
                    "dimensions": set(),
                    "reason_codes": set(),
                    "occurrence": False,
                    "phase": False,
                })
                descriptor["dimensions"].add(DIMENSION_OCCURRENCE_PATH)
                descriptor["reason_codes"].update((
                    reason_code,
                    "exact_contiguous_source_subpath",
                ))
                descriptor["occurrence"] = True
        if not any(
                descriptor["occurrence"]
                for descriptor in descriptors.values()):
            raise ReportProvenanceError(
                "path-local policy selection lacks occurrence alternatives"
            )

        phase_geometries, phase_status = (
            _architecture_phase_alternative_geometries(
                source_sequence,
                source_hash,
                representative_geometry,
                source_phase_arbitration,
                expected_path_id,
                expected_superlocus_id,
            )
        )
        for geometry in phase_geometries:
            descriptor = descriptors.setdefault(geometry, {
                "dimensions": set(),
                "reason_codes": set(),
                "occurrence": False,
                "phase": False,
            })
            descriptor["dimensions"].add(DIMENSION_REPEAT_PHASE)
            descriptor["reason_codes"].add(
                "source_phase_arbitration_explicit_alternative"
            )
            descriptor["phase"] = True

        relationships = []
        for ordinal, geometry in enumerate(sorted(descriptors), 1):
            descriptor = descriptors[geometry]
            if descriptor["occurrence"] and descriptor["phase"]:
                relationship_label = (
                    "path-local-occurrence-and-phase-alternative-{:04d}"
                ).format(ordinal)
            elif descriptor["occurrence"]:
                relationship_label = (
                    "path-local-occurrence-subpath-{:04d}"
                ).format(ordinal)
            else:
                relationship_label = (
                    "source-phase-alternative-{:04d}"
                ).format(ordinal)
            relationships.append(AlternativeRelationshipInput(
                geometry=_architecture_geometry_input_from_intervals(
                    source_sequence,
                    source_hash,
                    geometry,
                    "architecture_boundary_alternative[{}]".format(
                        ordinal - 1
                    ),
                ),
                label=relationship_label,
                dimensions=tuple(
                    dimension for dimension in (
                        DIMENSION_OCCURRENCE_PATH,
                        DIMENSION_REPEAT_PHASE,
                    )
                    if dimension in descriptor["dimensions"]
                ),
                reason_codes=tuple(sorted(descriptor["reason_codes"])),
            ))

        occurrence_path = BoundaryDimensionInput(
            status=STATUS_POLICY_SELECTED,
            selection_policy_id=(
                "source-superlocus-path-local-envelope:1.0.0"
            ),
            reason_codes=(
                "finite_contained_terminal_closure",
                "unique_exact_contiguous_envelope_selected",
            ),
        )
        if phase_status == "selected":
            repeat_phase = BoundaryDimensionInput(
                status=STATUS_POLICY_SELECTED,
                selection_policy_id=(
                    "source-cross-family-phase-arbitration"
                ),
                reason_codes=(
                    "reported_phase_selected_by_source_policy",
                    "reported_geometry_is_representative_not_exact",
                ),
            )
        elif phase_status == "ambiguous_representative":
            repeat_phase = BoundaryDimensionInput(
                status=STATUS_AMBIGUOUS,
                reason_codes=(
                    "plausible_source_valid_repeat_phases_retained",
                    "reported_geometry_is_representative_not_exact",
                ),
            )
        else:
            repeat_phase = BoundaryDimensionInput()

        return source_normalizer.normalize(
            representative_input,
            alternative_relationships=tuple(relationships),
            occurrence_path=occurrence_path,
            repeat_phase=repeat_phase,
        ).as_dict()
    except ReportProvenanceError:
        raise
    except (AttributeError, IndexError, KeyError, TypeError, ValueError) as error:
        raise ReportProvenanceError(
            "cannot build canonical architecture boundary: {}".format(error)
        ) from error


def _periodic_candidate_from_intervals(source, value, label):
    intervals = _intervals(value, len(source), label)
    if not intervals:
        raise ReportProvenanceError(label + " must contain repeats")
    try:
        candidate = build_source_backed_candidate_from_boundaries(
            source,
            None,
            tuple(start for start, _end in intervals),
            tuple(end for _start, end in intervals),
        )
    except (TypeError, ValueError) as error:
        raise ReportProvenanceError(
            label + " cannot reconstruct a source-backed parent"
        ) from error
    return candidate, intervals


def _periodic_boundary_assessment(
    source,
    source_hash,
    parent,
    representative,
    proposal,
    replayed,
):
    """Independently rebuild the transaction's exact boundary semantics."""

    relationships = []
    selected = proposal.representative
    parent_dimensions = []
    if (
        len(proposal.parent_repeat_intervals)
        != len(selected.repeat_intervals)
        or selected.terminal_added_slots
        or selected.internal_missing_slots
    ):
        parent_dimensions.append(DIMENSION_OCCURRENCE_PATH)
    if any(
        mapping.parent_interval != mapping.proposed_interval
        for mapping in selected.parent_mappings
    ):
        parent_dimensions.append(DIMENSION_REPEAT_PHASE)
    if parent_dimensions:
        relationships.append(AlternativeRelationshipInput(
            geometry=geometry_input_from_candidate(parent, source_hash),
            label="accepted-parent-path",
            dimensions=tuple(parent_dimensions),
            reason_codes=(
                "accepted_parent_source_replayed",
                "periodic_lattice_policy_selected_representative",
            ),
        ))
    for index, candidate in enumerate(replayed):
        if index == proposal.representative_index:
            continue
        relationships.append(AlternativeRelationshipInput(
            geometry=geometry_input_from_candidate(candidate, source_hash),
            label="periodic-lattice-phase-alternative-{}".format(index),
            dimensions=(DIMENSION_REPEAT_PHASE,),
            reason_codes=(
                "source_valid_plausible_repeat_phase",
                "repeat_starts_and_count_invariant",
            ),
        ))

    has_occurrence = any(
        DIMENSION_OCCURRENCE_PATH in relationship.dimensions
        for relationship in relationships
    )
    has_phase = any(
        DIMENSION_REPEAT_PHASE in relationship.dimensions
        for relationship in relationships
    )
    occurrence = (
        BoundaryDimensionInput(
            status=STATUS_POLICY_SELECTED,
            selection_policy_id=(
                "source-periodic-lattice-phase-then-terminal:"
                + proposal.rule_version
            ),
            reason_codes=(
                "occurrence_count_source_lattice_supported",
                "reported_path_selected_by_periodic_lattice_policy",
            ),
        )
        if has_occurrence
        else BoundaryDimensionInput()
    )
    repeat_phase = (
        BoundaryDimensionInput(
            status=STATUS_AMBIGUOUS,
            reason_codes=(
                "plausible_source_valid_repeat_phases_retained",
                "reported_geometry_is_representative_not_exact",
            ),
        )
        if has_phase
        else BoundaryDimensionInput()
    )
    return BoundaryAssessmentNormalizer(source, source_hash).normalize(
        geometry_input_from_candidate(representative, source_hash),
        alternative_relationships=tuple(relationships),
        occurrence_path=occurrence,
        repeat_phase=repeat_phase,
    ).as_dict()


def _validate_periodic_tandem_records(
    source,
    source_hash,
    parent_geometry,
    alternative_geometries,
    parent_direct,
    raw_reviews,
    safety_evidence,
    label,
):
    parent_direct = _validate_tandem(
        source,
        source_hash,
        _mapping(parent_direct, label + ".parent_direct_tandem"),
        parent_geometry,
        label + ".parent_direct_tandem",
        (parent_geometry,),
    )
    if (
        parent_direct["supported"] is not False
        or parent_direct["decision"] != "not_supported"
        or parent_direct["reason"] != safety_evidence.tandem_decision_reason
    ):
        raise ReportProvenanceError(
            label + " parent tandem authorization is inconsistent"
        )

    unique_family = []
    seen = set()
    for geometry in (parent_geometry,) + tuple(alternative_geometries):
        if geometry.geometry_id in seen:
            continue
        seen.add(geometry.geometry_id)
        unique_family.append(geometry)
    expected_entries = (("parent", None, parent_geometry),) + tuple(
        ("alternative", index, geometry)
        for index, geometry in enumerate(alternative_geometries)
    )
    reviews = _sequence(raw_reviews, label + ".complete_family_tandem_reviews")
    if len(reviews) != len(expected_entries):
        raise ReportProvenanceError(
            label + " complete-family tandem review count is stale"
        )
    for index, (raw_review, expected) in enumerate(
        zip(reviews, expected_entries)
    ):
        review_label = "{}.complete_family_tandem_reviews[{}]".format(
            label,
            index,
        )
        review = _mapping(raw_review, review_label)
        _keys(review, _PERIODIC_REVIEW_KEYS, review_label)
        role, alternative_index, geometry = expected
        if (
            review["role"] != role
            or review["alternative_index"] != alternative_index
            or _intervals(
                review["repeat_intervals"],
                len(source),
                review_label + ".repeat_intervals",
            ) != geometry.repeat_intervals
        ):
            raise ReportProvenanceError(
                review_label + " does not identify the expected family member"
            )
        decision = _validate_tandem(
            source,
            source_hash,
            _mapping(review["decision"], review_label + ".decision"),
            geometry,
            review_label + ".decision",
            tuple(unique_family),
        )
        if (
            decision["supported"] is not False
            or decision["decision"] != "not_supported"
        ):
            raise ReportProvenanceError(
                review_label + " does not authorize a CRISPR repair"
            )


def _validate_periodic_classification_basis(
    source,
    source_hash,
    candidate,
    basis,
    *,
    row_score,
    allowed_methods,
    depth,
    validation_context,
    category_index,
    boundary_assessment,
):
    """Regenerate an accepted periodic-lattice repair from original DNA."""

    label = "classification_basis"
    _keys(basis, _PERIODIC_BASIS_KEYS, label)
    if category_index != 2:
        raise ReportProvenanceError(
            "periodic-lattice repair may only classify Possible rows"
        )
    if (
        basis["classification_ceiling"] != _PERIODIC_CLASSIFICATION_CEILING
        or basis["boundary_status"] != "AMBIGUOUS"
        or basis["exact_boundary_certified"] is not False
    ):
        raise ReportProvenanceError(
            "periodic-lattice basis overstates its classification or boundary"
        )

    record_label = label + ".source_periodic_lattice_repair"
    record = _mapping(basis["source_periodic_lattice_repair"], record_label)
    _keys(record, _PERIODIC_RECORD_KEYS, record_label)
    source_category = _integer(
        record["source_category_index"],
        record_label + ".source_category_index",
    )
    if (
        source_category not in _ACCEPTED_CATEGORY_INDEXES
        or record["source_category_name"] != _CATEGORY_NAMES[source_category]
    ):
        raise ReportProvenanceError(record_label + " source category is stale")
    _integer(record["source_row_ordinal"], record_label + ".source_row_ordinal")
    try:
        _canonical_json(record["source_locus_key"])
        _canonical_json(record["fresh_scoring_metadata"])
    except (TypeError, ValueError) as error:
        raise ReportProvenanceError(
            record_label + " contains non-JSON transaction metadata"
        ) from error
    expected_transition = {
        "source_category_index": source_category,
        "source_category_name": _CATEGORY_NAMES[source_category],
        "target_category_index": 2,
        "target_category_name": _CATEGORY_NAMES[2],
    }
    if _mapping(basis["transition"], label + ".transition") != expected_transition:
        raise ReportProvenanceError(label + " periodic transition is stale")
    if (
        record["rule"] != _PERIODIC_TRANSACTION_RULE
        or record["action"] != _PERIODIC_ACTION
        or record["reason_codes"]
        != ["source_periodic_lattice_repaired_possible"]
        or record["coordinate_system"] != _COORDINATE_SYSTEM
        or record["classification_ceiling"]
        != _PERIODIC_CLASSIFICATION_CEILING
        or record["target_category_index"] != 2
        or record["target_category_name"] != _CATEGORY_NAMES[2]
        or record["boundary_status"] != "AMBIGUOUS"
        or record["exact_boundary_certified"] is not False
        or record["collision_with"] != []
    ):
        raise ReportProvenanceError(
            record_label + " is not a collision-free Possible repair"
        )

    parent_score = _finite(record["parent_score"], record_label + ".parent_score")
    fresh_score = _finite(record["fresh_score"], record_label + ".fresh_score")
    if row_score is None or fresh_score != _finite(row_score, "row_score"):
        raise ReportProvenanceError(
            record_label + " fresh score differs from the reported row"
        )
    parent, parent_intervals = _periodic_candidate_from_intervals(
        source,
        record["parent_repeat_intervals"],
        record_label + ".parent_repeat_intervals",
    )
    parent_geometry = _geometry_from_candidate(source, parent)

    proposal_mapping = _mapping(
        record["proposal"],
        record_label + ".proposal",
    )
    rule_mapping = _boundary_required_mapping(
        proposal_mapping,
        "rule",
        record_label + ".proposal",
    )
    _keys(
        rule_mapping,
        frozenset(SourcePeriodicLatticeRule.__dataclass_fields__),
        record_label + ".proposal.rule",
    )
    safety_mapping = _boundary_required_mapping(
        proposal_mapping,
        "safety_evidence",
        record_label + ".proposal",
    )
    _keys(
        safety_mapping,
        frozenset(SourcePeriodicLatticeSafetyEvidence.__dataclass_fields__),
        record_label + ".proposal.safety_evidence",
    )
    try:
        rule = SourcePeriodicLatticeRule(**rule_mapping)
        safety = SourcePeriodicLatticeSafetyEvidence(**safety_mapping)
        expected_screen = source_periodic_lattice_geometry_prefilter(
            source,
            parent,
            rule=rule,
        )
        expected_proposal = propose_source_periodic_lattice_repair(
            source,
            parent,
            safety,
            rule=rule,
        )
    except (TypeError, ValueError) as error:
        raise ReportProvenanceError(
            record_label + " cannot regenerate its periodic-lattice proposal"
        ) from error
    if _canonical_json(record["geometry_prefilter"]) != _canonical_json(
        asdict(expected_screen)
    ):
        raise ReportProvenanceError(
            record_label + " geometry prefilter is not independently reproducible"
        )
    if _canonical_json(proposal_mapping) != _canonical_json(
        asdict(expected_proposal)
    ):
        raise ReportProvenanceError(
            record_label + " proposal differs from exhaustive source replay"
        )
    if (
        expected_proposal.status != "proposed"
        or expected_proposal.classification_ceiling
        != _PERIODIC_CLASSIFICATION_CEILING
        or expected_proposal.exact_boundary_status != "ambiguous"
        or expected_proposal.representative is None
        or not expected_proposal.repeat_count_invariant
        or not expected_proposal.repeat_start_locations_invariant
        or expected_proposal.source_sha256 != source_hash
        or expected_proposal.parent_repeat_intervals != parent_intervals
        or safety.source_sha256 != source_hash
        or safety.coding_context_assessed is not True
        or safety.tandem_evidence_assessed is not True
        or safety.coding_tandem_supported is not False
    ):
        raise ReportProvenanceError(
            record_label + " proposal is not an authorized ambiguous repair"
        )
    try:
        replayed = replay_source_periodic_lattice_proposal(
            source,
            expected_proposal,
        )
    except (TypeError, ValueError) as error:
        raise ReportProvenanceError(
            record_label + " alternatives do not replay from source"
        ) from error
    alternative_geometries = tuple(
        _geometry_from_candidate(source, alternative)
        for alternative in replayed
    )
    current = _geometry_from_candidate(source, candidate)
    representative = expected_proposal.representative
    recorded_representative = _intervals(
        record["representative_repeat_intervals"],
        len(source),
        record_label + ".representative_repeat_intervals",
    )
    if (
        current.repeat_intervals != representative.repeat_intervals
        or recorded_representative != current.repeat_intervals
    ):
        raise ReportProvenanceError(
            record_label + " representative differs from the reported candidate"
        )
    if _canonical_json(record["target_locus_key"]) != _canonical_json(
        list(current.interval)
    ):
        raise ReportProvenanceError(
            record_label + " target locus is not the representative envelope"
        )

    _validate_periodic_tandem_records(
        source,
        source_hash,
        parent_geometry,
        alternative_geometries,
        record["parent_direct_tandem"],
        record["complete_family_tandem_reviews"],
        safety,
        record_label,
    )

    prior_boundary = basis["prior_boundary_assessment"]
    if prior_boundary is not None:
        _revalidate_boundary_claim(
            source,
            source_hash,
            parent,
            prior_boundary,
            label + ".prior_boundary_assessment",
        )
    prior = basis["prior_classification_basis"]
    if prior is not None:
        validate_classification_basis(
            source,
            source_hash,
            parent,
            _mapping(prior, label + ".prior_classification_basis"),
            row_score=parent_score,
            allowed_methods=allowed_methods,
            depth=depth + 1,
            validation_context=validation_context,
            category_index=source_category,
            boundary_assessment=prior_boundary,
        )

    expected_boundary = _periodic_boundary_assessment(
        source,
        source_hash,
        parent,
        candidate,
        expected_proposal,
        replayed,
    )
    claimed_boundary = _revalidate_boundary_claim(
        source,
        source_hash,
        candidate,
        boundary_assessment,
        "boundary_assessment",
    )
    if _canonical_json(claimed_boundary) != _canonical_json(expected_boundary):
        raise ReportProvenanceError(
            "boundary_assessment differs from periodic-lattice replay"
        )
    if record["boundary_assessment_id"] != expected_boundary["assessment_id"]:
        raise ReportProvenanceError(
            record_label + " boundary assessment ID is stale"
        )
    return basis


def _architecture_expected_score_band(score):
    if score < 0.5:
        return "score_below_0.5"
    if score >= 0.75:
        return "score_at_least_0.75"
    return "score_0.5_to_below_0.75"


def _validate_architecture_coding_context(
        source, candidate_geometry, raw_context, label):
    context = _mapping(raw_context, label)
    _keys(context, frozenset((
        "coordinate_system",
        "candidate_interval",
        "coding_overlap_fraction",
        "context_cds_score",
        "selected_cds",
    )), label)
    expected_interval = (
        candidate_geometry.repeat_intervals[0][0],
        candidate_geometry.repeat_intervals[-1][1],
    )
    if (
            context["coordinate_system"] != _COORDINATE_SYSTEM
            or _interval(
                context["candidate_interval"],
                len(source),
                label + ".candidate_interval",
            ) != expected_interval):
        raise ReportProvenanceError(
            label + " candidate geometry binding is stale"
        )
    overlap = _finite(
        context["coding_overlap_fraction"],
        label + ".coding_overlap_fraction",
    )
    if not 0.0 <= overlap <= 1.0:
        raise ReportProvenanceError(
            label + ".coding_overlap_fraction must be in [0, 1]"
        )
    _finite(context["context_cds_score"], label + ".context_cds_score")
    selected = context["selected_cds"]
    if selected is None:
        if overlap != 0.0:
            raise ReportProvenanceError(
                label + " nonzero coding overlap lacks a selected CDS"
            )
        return
    selected = _mapping(selected, label + ".selected_cds")
    _keys(selected, frozenset((
        "start", "end", "strand", "confidence",
    )), label + ".selected_cds")
    selected_interval = _interval(
        (selected["start"], selected["end"]),
        len(source),
        label + ".selected_cds.interval",
    )
    if selected["strand"] not in ("+", "-"):
        raise ReportProvenanceError(label + " selected CDS strand is invalid")
    _finite(selected["confidence"], label + ".selected_cds.confidence")
    overlap_bp = max(
        0,
        min(expected_interval[1], selected_interval[1])
        - max(expected_interval[0], selected_interval[0]),
    )
    expected_overlap = overlap_bp / float(
        expected_interval[1] - expected_interval[0]
    )
    if not math.isclose(overlap, expected_overlap, rel_tol=0.0, abs_tol=1e-12):
        raise ReportProvenanceError(label + " coding overlap fraction is stale")


def _validate_architecture_decision(
        basis, candidate_geometry, row_score, label):
    policy = _mapping(basis["policy"], label + ".policy")
    expected_policy = AUDITED_ARCHITECTURE_RESCUE_POLICY.as_dict()
    if _canonical_json(policy) != _canonical_json(expected_policy):
        raise ReportProvenanceError(label + " policy is stale")
    decision = _mapping(basis["decision"], label + ".decision")
    if (
            decision.get("decision") != "eligible"
            or decision.get("eligible") is not True
            or decision.get("requires_review") is not False
            or decision.get("rejection_codes") != []
            or _canonical_json(decision.get("policy"))
            != _canonical_json(policy)):
        raise ReportProvenanceError(label + " decision is not eligible")
    stats = _mapping(decision.get("stats"), label + ".decision.stats")
    expected_start = candidate_geometry.repeat_intervals[0][0]
    expected_end = candidate_geometry.repeat_intervals[-1][1]
    if (
            stats.get("array_start") != expected_start
            or stats.get("array_end") != expected_end
            or stats.get("repeat_count")
            != len(candidate_geometry.repeat_intervals)
            or stats.get("spacer_count")
            != len(candidate_geometry.spacer_intervals)):
        raise ReportProvenanceError(
            label + " decision geometry statistics are stale"
        )
    classifier_score = _finite(
        stats.get("classifier_score"),
        label + ".decision.stats.classifier_score",
    )
    if not math.isclose(
            classifier_score, row_score, rel_tol=0.0, abs_tol=1e-12):
        raise ReportProvenanceError(
            label + " decision classifier score differs from its row"
        )


def _validate_architecture_rescue_classification_basis(
        source, source_hash, candidate, basis, *, row_score,
        category_index, boundary_assessment):
    label = "classification_basis"
    _keys(basis, _ARCHITECTURE_RESCUE_KEYS, label)
    if (
            basis["method"] != SOURCE_ARCHITECTURE_RESCUE_METHOD
            or basis["source_sha256"] != source_hash
            or basis["geometry_derivation"]
            != "source_backed_augmentation"):
        raise ReportProvenanceError(
            "architecture rescue source/method binding is stale"
        )
    score = _finite(row_score, "reported architecture rescue score")
    if basis["raw_score_band"] != _architecture_expected_score_band(score):
        raise ReportProvenanceError(
            "architecture rescue raw score band is stale"
        )
    candidate_geometry = geometry_input_from_candidate(candidate, source_hash)
    _validate_architecture_decision(
        basis, candidate_geometry, score, label
    )
    _validate_architecture_coding_context(
        source,
        candidate_geometry,
        basis["coding_context"],
        label + ".coding_context",
    )

    superlocus_id = _architecture_nonempty_string(
        basis["source_superlocus_id"],
        label + ".source_superlocus_id",
    )
    lineage_path_id = _architecture_nonempty_string(
        basis["lineage_path_id"], label + ".lineage_path_id"
    )
    parent_path_ids = _architecture_path_id_sequence(
        basis["parent_path_ids"], label + ".parent_path_ids"
    )
    family_path_ids = _architecture_path_id_sequence(
        basis["family_path_ids"], label + ".family_path_ids"
    )
    if (
            lineage_path_id not in family_path_ids
            or not set(parent_path_ids).issubset(family_path_ids)):
        raise ReportProvenanceError(
            "architecture rescue lineage family is incomplete"
        )
    seed_roles = _mapping(basis["seed_roles"], label + ".seed_roles")
    if (
            any(
                not isinstance(path_id, str)
                or path_id not in family_path_ids
                or role not in ("primary", "satellite")
                for path_id, role in seed_roles.items()
            )
            or list(seed_roles) != sorted(seed_roles)):
        raise ReportProvenanceError(
            "architecture rescue seed roles are noncanonical"
        )
    if (
            basis["transition_action"] not in ("add", "replace")
            or basis["transition_reason"] != "source_architecture_eligible"):
        raise ReportProvenanceError(
            "architecture rescue transition is not committable"
        )
    _sequence(basis["migrated_locus_keys"], label + ".migrated_locus_keys")
    dominated_geometries = []
    for index, raw_geometry in enumerate(_sequence(
            basis["dominated_paths"], label + ".dominated_paths")):
        geometry = tuple(_intervals(
            raw_geometry,
            len(source),
            "{}.dominated_paths[{}]".format(label, index),
        ))
        if not _architecture_exact_contiguous_subpath(
                geometry, candidate_geometry.repeat_intervals):
            raise ReportProvenanceError(
                "architecture rescue dominated path is not a subpath"
            )
        dominated_geometries.append(geometry)
    if len(set(dominated_geometries)) != len(dominated_geometries):
        raise ReportProvenanceError(
            "architecture rescue dominated paths are duplicated"
        )

    convergence = _mapping(
        basis["source_superlocus_convergence"],
        label + ".source_superlocus_convergence",
    )
    _keys(
        convergence,
        _ARCHITECTURE_CONVERGENCE_KEYS,
        label + ".source_superlocus_convergence",
    )
    _architecture_nonempty_string(
        convergence["reason"],
        label + ".source_superlocus_convergence.reason",
    )
    _architecture_nonempty_string(
        convergence["path_local_reason"],
        label + ".source_superlocus_convergence.path_local_reason",
    )
    if convergence["architecture_commit_eligible"] is not True:
        raise ReportProvenanceError(
            "architecture rescue row was not commit eligible"
        )

    global_status = convergence["status"]
    local_status = convergence["path_local_status"]
    proof = convergence["path_local_proof"]
    phase = basis["source_phase_arbitration"]
    phase_status = None
    if phase is not None:
        _phase_geometries, phase_status = (
            _architecture_phase_alternative_geometries(
                source,
                source_hash,
                tuple(candidate_geometry.repeat_intervals),
                phase,
                lineage_path_id,
                superlocus_id,
            )
        )
    expected_category = (
        2
        if phase_status == "ambiguous_representative"
        else (0 if score >= 0.75 else 2)
    )
    if category_index != expected_category:
        raise ReportProvenanceError(
            "architecture rescue row is in the wrong report category"
        )

    if global_status == "nonconverged":
        if local_status != "converged" or not isinstance(proof, dict):
            raise ReportProvenanceError(
                "nonconverged architecture rescue lacks a local proof"
            )
        try:
            replayed = (
                replay_source_superlocus_path_local_convergence_proof(
                    source,
                    proof,
                    expected_candidate=candidate,
                    expected_path_id=lineage_path_id,
                    expected_superlocus_id=superlocus_id,
                )
            )
        except ValueError as error:
            raise ReportProvenanceError(
                "architecture rescue path-local replay failed: {}".format(
                    error
                )
            ) from error
        bounded_path_ids = tuple(sorted(
            record["path_id"] for record in replayed["bounded_paths"]
        ))
        envelope_witness = next(
            record for record in replayed["bounded_paths"]
            if record["path_id"] == lineage_path_id
        )
        if (
                convergence["path_local_reason"] != replayed["reason"]
                or family_path_ids != bounded_path_ids
                or parent_path_ids
                != tuple(envelope_witness["parent_path_ids"])):
            raise ReportProvenanceError(
                "architecture rescue lineage differs from path-local replay"
            )
        expected_boundary = (
            canonical_source_architecture_rescue_boundary_assessment(
                source,
                source_hash,
                candidate,
                replayed,
                phase,
                expected_path_id=lineage_path_id,
                expected_superlocus_id=superlocus_id,
            )
        )
        claimed_boundary = _revalidate_boundary_claim(
            source,
            source_hash,
            candidate,
            boundary_assessment,
            "boundary_assessment",
        )
        if _canonical_json(claimed_boundary) != _canonical_json(
                expected_boundary):
            raise ReportProvenanceError(
                "boundary_assessment differs from path-local replay"
            )
    elif global_status == "converged":
        if (
                local_status != "not_required"
                or convergence["path_local_reason"]
                != "superlocus_converged"
                or proof is not None):
            raise ReportProvenanceError(
                "globally converged rescue falsely claims local replay"
            )
        if boundary_assessment is not None:
            _revalidate_boundary_claim(
                source,
                source_hash,
                candidate,
                boundary_assessment,
                "boundary_assessment",
            )
    elif global_status == "legacy_unspecified":
        if (
                local_status != "legacy_unspecified"
                or convergence["reason"] != "legacy_unspecified"
                or convergence["path_local_reason"] != "legacy_unspecified"
                or proof is not None):
            raise ReportProvenanceError(
                "legacy architecture convergence fields are inconsistent"
            )
        if boundary_assessment is not None:
            _revalidate_boundary_claim(
                source,
                source_hash,
                candidate,
                boundary_assessment,
                "boundary_assessment",
            )
    else:
        raise ReportProvenanceError(
            "architecture rescue convergence status is unsupported"
        )
    return basis


@dataclass(frozen=True)
class _RegularGeometry:
    mapping: dict
    candidate: object
    interval: tuple
    repeat_intervals: tuple
    geometry_sha256: str


@dataclass(frozen=True)
class _RegularReplayState:
    snapshot: dict
    regenerated_batch: object
    pretransaction_categories: tuple
    additions: tuple
    requests_by_id: dict
    drafts_by_request_id: dict
    parents_by_locator: dict


class _HashableJsonDict(dict):
    def __hash__(self):
        return hash(_canonical_json(self))


def _regular_hash(value, label):
    try:
        encoded = _canonical_json(value).encode("ascii")
    except (TypeError, ValueError, UnicodeEncodeError) as error:
        raise ReportProvenanceError(
            label + " is not canonical finite JSON"
        ) from error
    return hashlib.sha256(encoded).hexdigest()


def _regular_hashable_json(value):
    if isinstance(value, list):
        return tuple(_regular_hashable_json(item) for item in value)
    if isinstance(value, dict):
        return _HashableJsonDict({
            key: _regular_hashable_json(item)
            for key, item in value.items()
        })
    return value


def _regular_geometry_mapping_from_candidate(source, candidate, label):
    try:
        result = validate_candidate_reconstruction(source, candidate)
        starts = tuple(candidate.list_repeat_starts)
        claimed_repeats = tuple(candidate.list_repeats)
        tuple(candidate.list_spacers)
        consensus = candidate.consensus
    except (AttributeError, TypeError, ValueError) as error:
        raise ReportProvenanceError(
            label + " is not a standard source-reconstructable candidate"
        ) from error
    if not result.valid:
        raise ReportProvenanceError(label + " does not reconstruct from source")
    if not isinstance(consensus, str) or not consensus:
        raise ReportProvenanceError(label + " consensus must be non-empty text")
    repeats = tuple(
        source[start:start + len(repeat)]
        for start, repeat in zip(starts, claimed_repeats)
    )
    intervals = tuple(
        (start, start + len(repeat))
        for start, repeat in zip(starts, repeats)
    )
    spacers = tuple(
        source[left[1]:right[0]]
        for left, right in zip(intervals, intervals[1:])
    )
    interval = intervals[0][0], intervals[-1][1]
    geometry_record = {
        "coordinate_system": _COORDINATE_SYSTEM,
        "repeat_starts": starts,
        "repeats": repeats,
        "spacers": spacers,
        "consensus": consensus,
    }
    geometry_sha256 = _regular_hash(geometry_record, label + ".geometry")
    return {
        "coordinate_system": _COORDINATE_SYSTEM,
        "source_interval": list(interval),
        "repeat_count": len(intervals),
        "repeat_intervals": [list(item) for item in intervals],
        "repeat_starts": list(starts),
        "repeat_sequences": list(repeats),
        "spacer_sequences": list(spacers),
        "consensus": consensus,
        "source_slice_sha256": _sha256(source[interval[0]:interval[1]]),
        "geometry_sha256": geometry_sha256,
    }


def _regular_geometry_from_mapping(source, raw_mapping, label):
    mapping = _mapping(raw_mapping, label)
    _keys(mapping, _REGULAR_GEOMETRY_KEYS, label)
    repeat_intervals = _intervals(
        mapping["repeat_intervals"],
        len(source),
        label + ".repeat_intervals",
    )
    if not repeat_intervals:
        raise ReportProvenanceError(label + " contains no repeats")
    starts = tuple(
        _integer(value, "{}.repeat_starts[{}]".format(label, index))
        for index, value in enumerate(
            _sequence(mapping["repeat_starts"], label + ".repeat_starts")
        )
    )
    if starts != tuple(start for start, _end in repeat_intervals):
        raise ReportProvenanceError(label + " repeat starts are stale")
    repeats = _sequence(mapping["repeat_sequences"], label + ".repeat_sequences")
    spacers = _sequence(mapping["spacer_sequences"], label + ".spacer_sequences")
    if any(not isinstance(item, str) for item in repeats + spacers):
        raise ReportProvenanceError(label + " sequences must be text")
    if len(repeats) != len(repeat_intervals) or len(spacers) != len(repeats) - 1:
        raise ReportProvenanceError(label + " component cardinality is stale")
    if tuple(
        source[start:end] for start, end in repeat_intervals
    ) != repeats:
        raise ReportProvenanceError(label + " repeat sequence is not source-derived")
    expected_spacers = tuple(
        source[left[1]:right[0]]
        for left, right in zip(repeat_intervals, repeat_intervals[1:])
    )
    if spacers != expected_spacers:
        raise ReportProvenanceError(label + " spacer sequence is not source-derived")
    consensus = mapping["consensus"]
    if not isinstance(consensus, str) or not consensus:
        raise ReportProvenanceError(label + " consensus must be non-empty text")
    candidate = SimpleNamespace(
        list_repeat_starts=starts,
        list_repeats=repeats,
        list_spacers=spacers,
        consensus=consensus,
    )
    expected = _regular_geometry_mapping_from_candidate(source, candidate, label)
    if _canonical_json(mapping) != _canonical_json(expected):
        raise ReportProvenanceError(
            label + " differs from independent original-source replay"
        )
    return _RegularGeometry(
        mapping=mapping,
        candidate=candidate,
        interval=tuple(expected["source_interval"]),
        repeat_intervals=repeat_intervals,
        geometry_sha256=expected["geometry_sha256"],
    )


def _regular_complete_family(source, raw_family, parent_candidate, label):
    family = _mapping(raw_family, label)
    _keys(family, _REGULAR_FAMILY_KEYS, label)
    if family["caller_declared_complete"] is not True:
        raise ReportProvenanceError(label + " is not caller-declared complete")
    raw_geometries = _sequence(
        family["candidate_geometries"],
        label + ".candidate_geometries",
    )
    if family["candidate_count"] != len(raw_geometries) or not raw_geometries:
        raise ReportProvenanceError(label + " candidate count is stale")
    parent_mapping = _regular_geometry_mapping_from_candidate(
        source,
        parent_candidate,
        label + ".parent",
    )
    parent_hash = parent_mapping["geometry_sha256"]
    candidates = []
    geometries = []
    parent_count = 0
    seen = set()
    for index, raw_geometry in enumerate(raw_geometries):
        geometry = _regular_geometry_from_mapping(
            source,
            raw_geometry,
            "{}.candidate_geometries[{}]".format(label, index),
        )
        if geometry.geometry_sha256 in seen:
            raise ReportProvenanceError(label + " contains duplicate geometry")
        seen.add(geometry.geometry_sha256)
        if geometry.geometry_sha256 == parent_hash:
            parent_count += 1
            candidates.append(parent_candidate)
        else:
            candidates.append(geometry.candidate)
        geometries.append(geometry)
    if parent_count != 1:
        raise ReportProvenanceError(label + " must contain the parent exactly once")
    return tuple(candidates), tuple(geometries)


def _regular_source_cds(source, raw_records):
    records = _sequence(raw_records, "report_replay_snapshot.source_cds")
    parsed = []
    expected_keys = frozenset(SourceCDS.__dataclass_fields__)
    for index, raw_record in enumerate(records):
        label = "report_replay_snapshot.source_cds[{}]".format(index)
        record = _mapping(raw_record, label)
        _keys(record, expected_keys, label)
        try:
            cds = SourceCDS(**record)
        except (TypeError, ValueError) as error:
            raise ReportProvenanceError(label + " is invalid") from error
        if cds.end > len(source):
            raise ReportProvenanceError(label + " lies outside the source")
        parsed.append(cds)
    ordered = tuple(sorted(set(parsed), key=lambda cds: (
        cds.start,
        cds.end,
        cds.is_complement,
        cds.confidence,
    )))
    if tuple(parsed) != ordered:
        raise ReportProvenanceError(
            "report_replay_snapshot.source_cds is duplicated or noncanonical"
        )
    return ordered


def _regular_row_basis(row):
    if not isinstance(row, (list, tuple)) or len(row) not in (2, 3, 4):
        raise ReportProvenanceError("category row shape is invalid")
    if len(row) in (2, 3):
        return None
    metadata = _mapping(row[3], "category row metadata")
    basis = metadata.get("classification_basis")
    return None if basis is None else _mapping(basis, "classification_basis")


def _regular_partition_categories(categories):
    if not isinstance(categories, (list, tuple)) or len(categories) != 5:
        raise ReportProvenanceError("regular-subpath replay requires five categories")
    originals = []
    additions = []
    for category_index, category in enumerate(categories):
        if not isinstance(category, dict):
            raise ReportProvenanceError("regular-subpath category must be a mapping")
        copied = {}
        for locus_key, rows in category.items():
            if isinstance(rows, (str, bytes)) or not isinstance(rows, (list, tuple)):
                raise ReportProvenanceError("category locus rows must be a sequence")
            retained = []
            for row_ordinal, row in enumerate(rows):
                basis = _regular_row_basis(row)
                if basis is not None and basis.get("method") == SOURCE_REGULAR_SUBPATH_METHOD:
                    if category_index != 2:
                        raise ReportProvenanceError(
                            "regular-subpath admission escaped the Possible category"
                        )
                    additions.append({
                        "category_index": category_index,
                        "locus_key": locus_key,
                        "row_ordinal": row_ordinal,
                        "row": row,
                        "basis": basis,
                    })
                else:
                    retained.append(row)
            if retained:
                copied[locus_key] = retained
        originals.append(copied)
    return tuple(originals), tuple(additions)


def _regular_restore_projected_diagnostics(source, snapshot, categories):
    if any(categories[index] for index in (1, 3, 4)):
        raise ReportProvenanceError(
            "accepted-only projection must remove every diagnostic category"
        )
    restored = [dict(category) for category in categories]
    manifest = _sequence(
        snapshot["pretransaction_row_manifest"],
        "snapshot.pretransaction_row_manifest",
    )
    for index, raw_record in enumerate(manifest):
        label = "snapshot.pretransaction_row_manifest[{}]".format(index)
        record = _mapping(raw_record, label)
        _keys(record, _REGULAR_MANIFEST_KEYS, label)
        category_index = _integer(
            record["category_index"], label + ".category_index"
        )
        row_ordinal = _integer(
            record["row_ordinal"], label + ".row_ordinal"
        )
        row_field_count = _integer(
            record["row_field_count"], label + ".row_field_count"
        )
        if (
            category_index >= len(_CATEGORY_NAMES)
            or record["category_name"] != _CATEGORY_NAMES[category_index]
            or row_field_count not in (3, 4)
        ):
            raise ReportProvenanceError(
                label + " category or row shape is stale"
            )
        geometry = _regular_geometry_from_mapping(
            source,
            record["geometry"],
            label + ".geometry",
        )
        score = _finite(record["score"], label + ".score")
        try:
            features = canonicalize_regular_subpath_row_value(
                record["features"], label + ".features"
            )
            metadata = canonicalize_regular_subpath_row_value(
                record["metadata"], label + ".metadata"
            )
        except ValueError as error:
            raise ReportProvenanceError(
                label + " row values are not canonical finite data"
            ) from error
        if (
            not isinstance(features, list)
            or (row_field_count == 3 and metadata is not None)
            or (row_field_count == 4 and not isinstance(metadata, dict))
            or record["features_sha256"]
            != _regular_hash(features, label + ".features")
            or record["metadata_sha256"]
            != _regular_hash(metadata, label + ".metadata")
        ):
            raise ReportProvenanceError(label + " row value binding is stale")
        unsigned = dict(record)
        supplied_binding = unsigned.pop("row_binding_sha256")
        if supplied_binding != _regular_hash(unsigned, label):
            raise ReportProvenanceError(label + " row binding is stale")
        if category_index not in (1, 3, 4):
            continue
        locus_key = _regular_hashable_json(record["locus_key"])
        try:
            hash(locus_key)
        except TypeError as error:
            raise ReportProvenanceError(
                label + " locus key is not hashable after replay"
            ) from error
        rows = restored[category_index].setdefault(locus_key, [])
        if row_ordinal != len(rows):
            raise ReportProvenanceError(
                label + " row ordinal is noncontiguous"
            )
        row = [score, geometry.candidate, features]
        if row_field_count == 4:
            row.append(metadata)
        rows.append(row)
    return tuple(restored)


def _regular_resolve_locator(categories, raw_locator, label):
    locator = _mapping(raw_locator, label)
    expected_keys = frozenset((
        "category_index",
        "category_name",
        "locus_key",
        "row_ordinal",
        "context_key",
    ))
    _keys(locator, expected_keys, label)
    category_index = _integer(locator["category_index"], label + ".category_index")
    if (
        category_index >= len(_CATEGORY_NAMES)
        or locator["category_name"] != _CATEGORY_NAMES[category_index]
    ):
        raise ReportProvenanceError(label + " category index/name is stale")
    row_ordinal = _integer(locator["row_ordinal"], label + ".row_ordinal")
    locus_json = _canonical_json(locator["locus_key"])
    matches = tuple(
        rows for locus_key, rows in categories[category_index].items()
        if _canonical_json(locus_key) == locus_json
    )
    if len(matches) != 1 or row_ordinal >= len(matches[0]):
        raise ReportProvenanceError(label + " does not resolve to one current row")
    return matches[0][row_ordinal], category_index


def _regular_snapshot_bindings(snapshot):
    bindings = _mapping(snapshot["bindings"], "report_replay_snapshot.bindings")
    expected_keys = frozenset((
        "row_value_normalizer_version",
        "source_cds_sha256",
        "accepted_geometries_sha256",
        "reviewed_diagnostic_parents_sha256",
        "pretransaction_row_manifest_sha256",
        "transaction_plan_sha256",
    ))
    _keys(bindings, expected_keys, "report_replay_snapshot.bindings")
    if bindings["row_value_normalizer_version"] != _REGULAR_ROW_NORMALIZER_VERSION:
        raise ReportProvenanceError("regular-subpath row normalizer version is stale")
    expected = {
        "source_cds_sha256": _regular_hash(
            snapshot["source_cds"], "snapshot source CDS"
        ),
        "accepted_geometries_sha256": _regular_hash(
            snapshot["accepted_geometries"], "snapshot accepted geometries"
        ),
        "reviewed_diagnostic_parents_sha256": _regular_hash(
            snapshot["reviewed_diagnostic_parents"], "snapshot reviewed parents"
        ),
        "pretransaction_row_manifest_sha256": _regular_hash(
            snapshot["pretransaction_row_manifest"], "snapshot row manifest"
        ),
        "transaction_plan_sha256": _regular_hash(
            snapshot["transaction_plan"], "snapshot transaction plan"
        ),
    }
    for key, value in expected.items():
        if bindings[key] != value:
            raise ReportProvenanceError(
                "report_replay_snapshot binding {!r} is stale".format(key)
            )


def _regular_validate_snapshot_header(source, source_hash, raw_snapshot):
    snapshot = _mapping(raw_snapshot, "report_replay_snapshot")
    _keys(snapshot, _REGULAR_SNAPSHOT_KEYS, "report_replay_snapshot")
    if (
        snapshot["schema"] != _REGULAR_SNAPSHOT_SCHEMA
        or snapshot["transaction_rule"] != _REGULAR_TRANSACTION_RULE
        or snapshot["coordinate_system"] != _COORDINATE_SYSTEM
        or snapshot["source_sha256"] != source_hash
        or snapshot["source_length"] != len(source)
    ):
        raise ReportProvenanceError(
            "report_replay_snapshot schema, rule, or source binding is stale"
        )
    unsigned = dict(snapshot)
    supplied_hash = unsigned.pop("snapshot_sha256")
    if supplied_hash != _regular_hash(unsigned, "report_replay_snapshot"):
        raise ReportProvenanceError("report_replay_snapshot digest is stale")
    _regular_snapshot_bindings(snapshot)
    expected_invariants = {
        "accepted_rows_are_complete_category_0_and_2_snapshot": True,
        "all_pretransaction_rows_source_replayed": True,
        "all_reviewed_parent_locators_canonical": True,
        "features_and_metadata_embedded_and_sha256_bound": True,
        "runtime_object_identity_serialized": False,
    }
    if snapshot["invariants"] != expected_invariants:
        raise ReportProvenanceError("report_replay_snapshot invariants are stale")
    return snapshot


def _regular_expected_boundary(source, source_hash, parent, proposal):
    return BoundaryAssessmentNormalizer(source, source_hash).normalize(
        geometry_input_from_candidate(proposal, source_hash),
        alternative_relationships=(
            AlternativeRelationshipInput(
                geometry=geometry_input_from_candidate(parent, source_hash),
                label="diagnostic-parent-occurrence-path",
                dimensions=(DIMENSION_OCCURRENCE_PATH,),
                reason_codes=(
                    "parent_source_replayed",
                    "proposal_is_policy_selected_parent_subpath",
                ),
            ),
        ),
        occurrence_path=BoundaryDimensionInput(
            status=STATUS_POLICY_SELECTED,
            selection_policy_id="source-regular-subpath-recovery:1.0.0",
            reason_codes=(
                "maximal_regular_segment_selected_between_extreme_breaks",
                "parent_occurrence_path_retained_as_alternative",
            ),
        ),
        repeat_phase=BoundaryDimensionInput(),
    ).as_dict()


def _regular_replay_snapshot(
    source,
    source_hash,
    snapshot,
    validation_context,
):
    snapshot = _regular_validate_snapshot_header(source, source_hash, snapshot)
    cache_key = snapshot["snapshot_sha256"]
    cached = validation_context._regular_replay_cache.get(cache_key)
    if cached is not None:
        if _canonical_json(cached.snapshot) != _canonical_json(snapshot):
            raise ReportProvenanceError(
                "one regular-subpath snapshot digest has two payloads"
            )
        return cached

    categories = validation_context.categories
    pretransaction, additions = _regular_partition_categories(categories)
    if validation_context.accepted_only_projection:
        pretransaction = _regular_restore_projected_diagnostics(
            source,
            snapshot,
            pretransaction,
        )
    source_cds = _regular_source_cds(source, snapshot["source_cds"])

    accepted_rows = tuple(
        row
        for category_index in _ACCEPTED_CATEGORY_INDEXES
        for rows in pretransaction[category_index].values()
        for row in rows
    )
    accepted_candidates = tuple(row[1] for row in accepted_rows)
    actual_accepted = sorted(
        (
            _regular_geometry_mapping_from_candidate(
                source,
                candidate,
                "pretransaction accepted candidate",
            )
            for candidate in accepted_candidates
        ),
        key=lambda geometry: (
            geometry["source_interval"],
            geometry["geometry_sha256"],
        ),
    )
    raw_accepted = _sequence(
        snapshot["accepted_geometries"],
        "snapshot.accepted_geometries",
    )
    for index, geometry in enumerate(raw_accepted):
        _regular_geometry_from_mapping(
            source,
            geometry,
            "snapshot.accepted_geometries[{}]".format(index),
        )
    if _canonical_json(actual_accepted) != _canonical_json(raw_accepted):
        raise ReportProvenanceError(
            "accepted geometry snapshot differs from current pretransaction rows"
        )

    plan = _mapping(snapshot["transaction_plan"], "snapshot.transaction_plan")
    _keys(
        plan,
        frozenset(("drafts", "score_requests")),
        "snapshot.transaction_plan",
    )
    raw_drafts = _sequence(plan["drafts"], "snapshot.transaction_plan.drafts")
    raw_requests = _sequence(
        plan["score_requests"],
        "snapshot.transaction_plan.score_requests",
    )
    for index, raw_draft in enumerate(raw_drafts):
        label = "snapshot.transaction_plan.drafts[{}]".format(index)
        draft = _mapping(raw_draft, label)
        for field_name in ("parent_geometry", "proposal_geometry"):
            geometry = draft.get(field_name)
            if geometry is not None:
                _regular_geometry_from_mapping(
                    source,
                    geometry,
                    label + "." + field_name,
                )
    for index, raw_request in enumerate(raw_requests):
        label = "snapshot.transaction_plan.score_requests[{}]".format(index)
        request = _mapping(raw_request, label)
        for field_name in ("parent_geometry", "proposal_geometry"):
            if field_name not in request:
                raise ReportProvenanceError(label + " omits " + field_name)
            _regular_geometry_from_mapping(
                source,
                request[field_name],
                label + "." + field_name,
            )

    context = SourceRegularSubpathContext(source)
    diagnostic_contexts = []
    parents_by_locator = {}
    reviewed = _sequence(
        snapshot["reviewed_diagnostic_parents"],
        "snapshot.reviewed_diagnostic_parents",
    )
    reviewed_keys = frozenset((
        "context_key",
        "parent_locator",
        "target_locus_key",
        "parent_geometry",
        "complete_family",
    ))
    for index, raw_review in enumerate(reviewed):
        label = "snapshot.reviewed_diagnostic_parents[{}]".format(index)
        review = _mapping(raw_review, label)
        _keys(review, reviewed_keys, label)
        parent_row, parent_category = _regular_resolve_locator(
            pretransaction,
            review["parent_locator"],
            label + ".parent_locator",
        )
        if parent_category not in (1, 3, 4):
            raise ReportProvenanceError(label + " parent is not diagnostic")
        locator = review["parent_locator"]
        if locator["context_key"] != review["context_key"]:
            raise ReportProvenanceError(label + " context key is stale")
        parent_mapping = _regular_geometry_mapping_from_candidate(
            source,
            parent_row[1],
            label + ".current_parent",
        )
        if _canonical_json(parent_mapping) != _canonical_json(
            review["parent_geometry"]
        ):
            raise ReportProvenanceError(label + " parent geometry is stale")
        family_candidates, _family_geometries = _regular_complete_family(
            source,
            review["complete_family"],
            parent_row[1],
            label + ".complete_family",
        )
        parent_identifier = _regular_hash(
            locator,
            label + ".parent_locator",
        )[:24]
        planner = propose_source_regular_subpaths_in_context(
            context,
            parent_row[1],
            parent_identifier=parent_identifier,
            parent_disposition=DIAGNOSTIC_PARENT,
            accepted_candidates=accepted_candidates,
        )
        parent_tandem = tandem_module.validate_accepted_family_tandem(
            source,
            parent_row[1],
            family_candidates,
            source_cds,
        )
        proposal_evidence = tuple(
            RegularSubpathProposalTandemEvidence(
                candidate_geometry_sha256=(
                    proposal.tandem_veto.candidate_geometry_sha256
                ),
                decision=tandem_module.validate_accepted_family_tandem(
                    source,
                    proposal.candidate,
                    family_candidates + (proposal.candidate,),
                    source_cds,
                ),
            )
            for proposal in planner.proposals
        )
        target_key = _regular_hashable_json(review["target_locus_key"])
        diagnostic_contexts.append(SourceRegularSubpathDiagnosticRowContext(
            context_key=review["context_key"],
            diagnostic_row=parent_row,
            target_locus_key=target_key,
            complete_family_candidates=family_candidates,
            supplied_parent_tandem_decision=parent_tandem,
            supplied_proposal_tandem_evidence=proposal_evidence,
            family_complete=True,
        ))
        locator_key = _canonical_json(locator)
        if locator_key in parents_by_locator:
            raise ReportProvenanceError("reviewed parent locator is duplicated")
        parents_by_locator[locator_key] = parent_row

    try:
        regenerated = plan_source_regular_subpath_recovery_batch(
            context,
            pretransaction,
            tuple(diagnostic_contexts),
            accepted_rows,
            source_cds,
        )
    except (RuntimeError, TypeError, ValueError) as error:
        raise ReportProvenanceError(
            "regular-subpath transaction plan cannot be regenerated"
        ) from error
    regenerated_snapshot = regenerated.report_replay_snapshot
    if _canonical_json(regenerated_snapshot) != _canonical_json(snapshot):
        raise ReportProvenanceError(
            "regular-subpath sealed snapshot differs from full batch replay"
        )

    requests_by_id = {}
    for raw_request in raw_requests:
        request = _mapping(raw_request, "snapshot score request")
        request_id = request.get("request_id")
        if not isinstance(request_id, str) or not request_id:
            raise ReportProvenanceError("snapshot score request ID is invalid")
        if request_id in requests_by_id:
            raise ReportProvenanceError("snapshot score request ID is duplicated")
        requests_by_id[request_id] = request
    drafts_by_request_id = {}
    for raw_draft in raw_drafts:
        draft = _mapping(raw_draft, "snapshot draft")
        request_id = draft.get("score_request_id")
        if request_id is None:
            continue
        if request_id in drafts_by_request_id:
            raise ReportProvenanceError("one request has multiple drafts")
        drafts_by_request_id[request_id] = draft
    if set(drafts_by_request_id) != set(requests_by_id):
        raise ReportProvenanceError("snapshot drafts and requests disagree")

    additions_by_id = {}
    for addition in additions:
        basis = addition["basis"]
        if _canonical_json(basis.get("report_replay_snapshot")) != _canonical_json(
            snapshot
        ):
            raise ReportProvenanceError(
                "regular-subpath admissions do not share one sealed snapshot"
            )
        fresh = _mapping(
            basis.get("fresh_classifier"),
            "regular-subpath fresh_classifier",
        )
        request_id = fresh.get("request_id")
        if request_id not in requests_by_id or request_id in additions_by_id:
            raise ReportProvenanceError(
                "regular-subpath admissions and score requests disagree"
            )
        request = requests_by_id[request_id]
        current_geometry = _regular_geometry_mapping_from_candidate(
            source,
            addition["row"][1],
            "regular-subpath admitted candidate",
        )
        if _canonical_json(current_geometry) != _canonical_json(
            request["proposal_geometry"]
        ):
            raise ReportProvenanceError(
                "regular-subpath admission changed proposal geometry"
            )
        if _canonical_json(addition["locus_key"]) != _canonical_json(
            request["target_locus_key"]
        ):
            raise ReportProvenanceError(
                "regular-subpath admission is stored under the wrong target locus"
            )
        additions_by_id[request_id] = addition
    if set(additions_by_id) != set(requests_by_id):
        raise ReportProvenanceError(
            "regular-subpath row conservation or admission cardinality failed"
        )

    state = _RegularReplayState(
        snapshot=snapshot,
        regenerated_batch=regenerated,
        pretransaction_categories=pretransaction,
        additions=tuple(additions_by_id.values()),
        requests_by_id=requests_by_id,
        drafts_by_request_id=drafts_by_request_id,
        parents_by_locator=parents_by_locator,
    )
    validation_context._regular_replay_cache[cache_key] = state
    return state


def _regular_row_metadata(row, label):
    if len(row) == 3:
        return {}
    return _mapping(row[3], label)


def _regular_validate_nested_basis(
    source,
    source_hash,
    candidate,
    row,
    raw_basis,
    raw_boundary,
    *,
    allowed_methods,
    depth,
    validation_context,
    category_index,
    label,
):
    if raw_boundary is not None:
        _revalidate_boundary_claim(
            source,
            source_hash,
            candidate,
            raw_boundary,
            label + ".boundary_assessment",
        )
    if raw_basis is None:
        return
    validate_classification_basis(
        source,
        source_hash,
        candidate,
        _mapping(raw_basis, label + ".classification_basis"),
        row_score=row[0],
        allowed_methods=allowed_methods,
        depth=depth + 1,
        validation_context=validation_context,
        category_index=category_index,
        boundary_assessment=raw_boundary,
        current_row=row,
    )


def _validate_regular_subpath_classification_basis(
    source,
    source_hash,
    candidate,
    basis,
    *,
    row_score,
    allowed_methods,
    depth,
    validation_context,
    category_index,
    boundary_assessment,
    current_row,
):
    """Regenerate one sealed v1.2 regular-subpath admission end to end."""

    label = "classification_basis"
    _keys(basis, _REGULAR_BASIS_KEYS, label)
    if category_index != 2 or current_row is None:
        raise ReportProvenanceError(
            "regular-subpath recovery may only classify a current Possible row"
        )
    if (
        basis["method"] != SOURCE_REGULAR_SUBPATH_METHOD
        or basis["rule"] != _REGULAR_TRANSACTION_RULE
        or basis["source_sha256"] != source_hash
        or basis["classification_cap"] != "Possible"
        or basis["bona_fide_admission_allowed"] is not False
        or basis["exact_boundary_certified"] is not False
    ):
        raise ReportProvenanceError(
            "regular-subpath rule, source, or classification cap is stale"
        )
    score = _finite(row_score, "regular-subpath reported score")
    if not 0.0 <= score <= 1.0:
        raise ReportProvenanceError(
            "regular-subpath reported score must be in [0, 1]"
        )
    if (
        not isinstance(current_row, (list, tuple))
        or len(current_row) != 4
        or current_row[1] is not candidate
    ):
        raise ReportProvenanceError(
            "regular-subpath current row binding is invalid"
        )
    try:
        canonicalize_regular_subpath_row_value(
            current_row[2], "regular-subpath current features"
        )
        canonicalize_regular_subpath_row_value(
            current_row[3], "regular-subpath current metadata"
        )
    except ValueError as error:
        raise ReportProvenanceError(
            "regular-subpath current row is not canonical finite data"
        ) from error

    state = _regular_replay_snapshot(
        source,
        source_hash,
        basis["report_replay_snapshot"],
        validation_context,
    )
    fresh = _mapping(basis["fresh_classifier"], label + ".fresh_classifier")
    _keys(
        fresh,
        frozenset((
            "request_id",
            "score",
            "evidence_method",
            "prior_classification_basis",
            "prior_boundary_assessment",
        )),
        label + ".fresh_classifier",
    )
    request_id = fresh["request_id"]
    if (
        not isinstance(request_id, str)
        or not request_id
        or request_id not in state.requests_by_id
        or not isinstance(fresh["evidence_method"], str)
        or not fresh["evidence_method"]
        or _finite(fresh["score"], "fresh classifier score") != score
    ):
        raise ReportProvenanceError(
            "regular-subpath fresh classifier binding is stale"
        )
    matching_additions = tuple(
        addition
        for addition in state.additions
        if addition["row"] is current_row
    )
    if (
        len(matching_additions) != 1
        or matching_additions[0]["basis"] is not basis
        or matching_additions[0]["category_index"] != 2
    ):
        raise ReportProvenanceError(
            "regular-subpath basis is not bound to its current category row"
        )
    addition = matching_additions[0]
    if addition["basis"]["fresh_classifier"]["request_id"] != request_id:
        raise ReportProvenanceError(
            "regular-subpath current row refers to another score request"
        )

    request = state.requests_by_id[request_id]
    draft = state.drafts_by_request_id[request_id]
    if (
        draft.get("action") != "scoring_required"
        or draft.get("reason_codes") != ["fresh_scoring_required"]
        or draft.get("collision_peer_proposal_ids") != []
        or draft.get("proposal_id") != request_id
        or draft.get("score_request_id") != request_id
    ):
        raise ReportProvenanceError(
            "regular-subpath admission was not an uncollided scoring request"
        )
    current_geometry = _regular_geometry_mapping_from_candidate(
        source,
        candidate,
        "regular-subpath reported proposal",
    )
    if _canonical_json(current_geometry) != _canonical_json(
        request["proposal_geometry"]
    ):
        raise ReportProvenanceError(
            "regular-subpath reported proposal differs from its request"
        )

    locator = request["parent_locator"]
    locator_key = _canonical_json(locator)
    parent_row = state.parents_by_locator.get(locator_key)
    if parent_row is None:
        raise ReportProvenanceError(
            "regular-subpath request does not resolve to a preserved parent"
        )
    parent_row_again, parent_category = _regular_resolve_locator(
        state.pretransaction_categories,
        locator,
        label + ".parent.locator",
    )
    if parent_row_again is not parent_row or parent_category not in (1, 3, 4):
        raise ReportProvenanceError(
            "regular-subpath parent is not preserved in a diagnostic category"
        )
    parent_metadata = _regular_row_metadata(
        parent_row, "regular-subpath parent metadata"
    )
    parent_prior = parent_metadata.get("classification_basis")
    parent_boundary = parent_metadata.get("boundary_assessment")
    _regular_validate_nested_basis(
        source,
        source_hash,
        parent_row[1],
        parent_row,
        parent_prior,
        parent_boundary,
        allowed_methods=allowed_methods,
        depth=depth,
        validation_context=validation_context,
        category_index=parent_category,
        label=label + ".parent.prior",
    )
    _regular_validate_nested_basis(
        source,
        source_hash,
        candidate,
        current_row,
        fresh["prior_classification_basis"],
        fresh["prior_boundary_assessment"],
        allowed_methods=allowed_methods,
        depth=depth,
        validation_context=validation_context,
        category_index=2,
        label=label + ".fresh_classifier.prior",
    )

    expected_boundary = _regular_expected_boundary(
        source,
        source_hash,
        parent_row[1],
        candidate,
    )
    claimed_boundary = _revalidate_boundary_claim(
        source,
        source_hash,
        candidate,
        boundary_assessment,
        "boundary_assessment",
    )
    if _canonical_json(claimed_boundary) != _canonical_json(expected_boundary):
        raise ReportProvenanceError(
            "boundary_assessment differs from regular-subpath replay"
        )
    if basis["boundary_assessment_id"] != expected_boundary["assessment_id"]:
        raise ReportProvenanceError(
            "regular-subpath boundary assessment ID is stale"
        )

    expected_basis = {
        "method": SOURCE_REGULAR_SUBPATH_METHOD,
        "rule": _REGULAR_TRANSACTION_RULE,
        "source_sha256": source_hash,
        "classification_cap": "Possible",
        "bona_fide_admission_allowed": False,
        "exact_boundary_certified": False,
        "boundary_assessment_id": expected_boundary["assessment_id"],
        "report_replay_snapshot": state.snapshot,
        "parent": {
            "locator": draft["parent_locator"],
            "score": draft["parent_score"],
            "prior_classification_basis": parent_prior,
            "preserved_in_diagnostic_category": True,
            "geometry": draft["parent_geometry"],
        },
        "proposal": {
            "proposal_id": draft["proposal_id"],
            "geometry": draft["proposal_geometry"],
            "occurrence_relationship": draft["occurrence_relationship"],
            "boundary_status": (
                "policy_selected_occurrence_path_repeat_phase_unassessed"
            ),
            "planner_boundary_status": "ambiguous_parent_phase",
            "planner_boundary_status_interpretation": (
                "conservative_planner_limitation_not_normalized_phase_evidence"
            ),
            "exact_biological_repeat_phase": False,
            "exact_array_boundaries": False,
            "normalized_boundary_dimensions": {
                "occurrence_path": "policy_selected",
                "repeat_phase": "unassessed",
            },
            "structural_planner": draft["structural_planner"],
        },
        "tandem_validation": {
            "parent": draft["parent_tandem_audit"],
            "proposal": draft["proposal_tandem_audit"],
            "complete_family": draft["complete_family"],
        },
        "accepted_overlap_protection": {
            "accepted_rows_identity_complete": True,
            "accepted_geometry_sha256s": list(
                state.regenerated_batch.accepted_geometry_sha256s
            ),
        },
        "fresh_classifier": fresh,
    }
    if _canonical_json(basis) != _canonical_json(expected_basis):
        raise ReportProvenanceError(
            "regular-subpath classification basis differs from sealed replay"
        )
    return basis


_LIKELY_TANDEM_BASIS_KEYS = frozenset((
    "method",
    "disposition",
    "tag",
    "prior_classification_basis",
    "likely_tandem",
    "transition",
))
_LIKELY_TANDEM_TRANSITION_KEYS = frozenset((
    "rule_name",
    "rule_version",
    "source_category_index",
    "source_category_name",
    "target_category_index",
    "target_category_name",
    "reported_disposition",
    "reported_tag",
    "locus_key",
    "row_ordinal",
))


def _validate_likely_tandem_classification_basis(
    source,
    source_hash,
    candidate,
    basis,
    *,
    row_score,
    allowed_methods,
    depth,
    validation_context,
    category_index,
    boundary_assessment,
):
    """Replay one tagged Low-score demotion against its source category."""

    _keys(basis, _LIKELY_TANDEM_BASIS_KEYS, "classification_basis")
    if category_index != 4:
        raise ReportProvenanceError(
            "tandem demotion may only classify Low score rows"
        )
    if (
        basis["method"] != LIKELY_TANDEM_METHOD
        or basis["disposition"] != LIKELY_TANDEM_DISPOSITION
        or basis["tag"] != TANDEM_TAG
    ):
        raise ReportProvenanceError(
            "tandem demotion method, disposition, or tag is stale"
        )
    transition = _mapping(
        basis["transition"], "classification_basis.transition"
    )
    _keys(
        transition,
        _LIKELY_TANDEM_TRANSITION_KEYS,
        "classification_basis.transition",
    )
    source_category_index = _integer(
        transition["source_category_index"],
        "classification_basis.transition.source_category_index",
    )
    if (
        source_category_index not in _ACCEPTED_CATEGORY_INDEXES
        or transition["source_category_name"]
        != _CATEGORY_NAMES[source_category_index]
        or transition["target_category_index"] != 4
        or transition["target_category_name"] != "Low score"
        or transition["reported_disposition"]
        != LIKELY_TANDEM_DISPOSITION
        or transition["reported_tag"] != TANDEM_TAG
        or transition["rule_name"]
        != LIKELY_TANDEM_TRANSACTION_RULE_NAME
        or transition["rule_version"]
        != LIKELY_TANDEM_TRANSACTION_RULE_VERSION
    ):
        raise ReportProvenanceError(
            "tandem demotion transition is stale"
        )

    categories = validation_context.categories
    locus_json = _canonical_json(transition["locus_key"])
    matches = tuple(
        rows
        for locus_key, rows in categories[source_category_index].items()
        if _canonical_json(locus_key) == locus_json
    )
    row_ordinal = _integer(
        transition["row_ordinal"],
        "classification_basis.transition.row_ordinal",
    )
    if len(matches) != 1 or row_ordinal >= len(matches[0]):
        raise ReportProvenanceError(
            "tandem demotion source row does not resolve uniquely"
        )
    source_row = matches[0][row_ordinal]
    if (
        not isinstance(source_row, (list, tuple))
        or len(source_row) not in (2, 3, 4)
        or row_score is None
        or _finite(source_row[0], "tandem source score")
        != _finite(row_score, "reported tandem score")
        or not _same_geometry(
            _geometry_from_candidate(source, source_row[1]),
            _geometry_from_candidate(source, candidate),
        )
    ):
        raise ReportProvenanceError(
            "tandem demotion changed score or candidate geometry"
        )

    prior = basis["prior_classification_basis"]
    source_prior = _regular_row_basis(source_row)
    if _canonical_json(prior) != _canonical_json(source_prior):
        raise ReportProvenanceError(
            "tandem demotion prior classification basis is stale"
        )
    expected_evidence = evaluate_likely_tandem(
        candidate.list_spacers
    ).as_dict()
    if _canonical_json(basis["likely_tandem"]) != _canonical_json(
        expected_evidence
    ):
        raise ReportProvenanceError(
            "tandem demotion evidence does not replay"
        )
    if (
        expected_evidence.get("decision") != LIKELY_TANDEM_DISPOSITION
        or expected_evidence.get("complete_link_threshold_satisfied") is not True
    ):
        raise ReportProvenanceError(
            "tandem demotion evidence no longer qualifies"
        )
    if prior is not None:
        validate_classification_basis(
            source,
            source_hash,
            source_row[1],
            prior,
            row_score=source_row[0],
            allowed_methods=allowed_methods,
            depth=depth + 1,
            validation_context=validation_context,
            category_index=source_category_index,
            boundary_assessment=boundary_assessment,
            current_row=source_row,
        )
    return basis


def validate_classification_basis(
    source_sequence,
    expected_source_sha256,
    candidate,
    raw_basis,
    *,
    row_score=None,
    allowed_methods,
    depth=0,
    validation_context=None,
    category_index=None,
    boundary_assessment=None,
    current_row=None,
):
    """Validate one classification basis, recursively, and return it unchanged."""

    if depth > MAX_CLASSIFICATION_BASIS_DEPTH:
        raise ReportProvenanceError(
            "classification_basis nesting exceeds depth {}".format(
                MAX_CLASSIFICATION_BASIS_DEPTH
            )
        )
    basis = _mapping(raw_basis, "classification_basis")
    method = basis.get("method")
    if method not in allowed_methods:
        raise ReportProvenanceError("classification_basis method is invalid")

    # Accepted rows can carry the independently computed coding-tandem
    # clearance added after their original classification transaction.  It
    # is orthogonal to the method-specific basis, so validate it against the
    # reported source geometry and then replay the basis without that wrapper.
    tandem_review = basis.get("accepted_family_tandem_review")
    if tandem_review is not None:
        if category_index not in _ACCEPTED_CATEGORY_INDEXES:
            raise ReportProvenanceError(
                "accepted-family tandem review is only valid on accepted rows"
            )
        core_basis = dict(basis)
        del core_basis["accepted_family_tandem_review"]
        validated_core = validate_classification_basis(
            source_sequence,
            expected_source_sha256,
            candidate,
            core_basis,
            row_score=row_score,
            allowed_methods=allowed_methods,
            depth=depth,
            validation_context=validation_context,
            category_index=category_index,
            boundary_assessment=boundary_assessment,
            current_row=current_row,
        )
        current = _geometry_from_candidate(source_sequence, candidate)
        _validate_tandem(
            source_sequence,
            expected_source_sha256,
            _mapping(
                tandem_review,
                "classification_basis.accepted_family_tandem_review",
            ),
            current,
            "classification_basis.accepted_family_tandem_review",
        )
        validated = dict(validated_core)
        validated["accepted_family_tandem_review"] = tandem_review
        return validated

    context = validation_context
    source = source_sequence
    if method == "source_backed_boundary_uncertain_reconstruction":
        from components.reconstruction_transaction import validate_boundary_family_basis
        report_context = _validation_context(
            source_sequence, expected_source_sha256, validation_context,
        )
        try:
            return validate_boundary_family_basis(
                source_sequence, candidate, basis, row_score, category_index,
                current_row, report_context.categories,
                ordinary_witness_keys=report_context.ordinary_occurrence_witnesses(),
            )
        except (ValueError, KeyError, TypeError, AttributeError) as error:
            raise ReportProvenanceError("boundary family source replay failed: {}".format(error)) from error
    if method == "source_backed_localized_spacer_context":
        from components.contextual_spacer_rescue import validate_contextual_spacer_basis
        report_context = _validation_context(
            source_sequence, expected_source_sha256, validation_context,
        )
        try:
            return validate_contextual_spacer_basis(
                source_sequence, candidate, basis, row_score, category_index,
                current_row, report_context.categories,
                ordinary_witness_keys=report_context.ordinary_occurrence_witnesses(),
            )
        except (ValueError, TypeError, KeyError, AttributeError) as error:
            raise ReportProvenanceError(
                "localized spacer context source replay failed: {}".format(error)
            ) from error
    if method == "source_backed_possible_reconstruction":
        from components.possible_reconstruction import validate_possible_reconstruction_basis
        report_context = _validation_context(
            source_sequence, expected_source_sha256, validation_context,
        )
        try:
            return validate_possible_reconstruction_basis(
                source_sequence, candidate, basis, row_score, category_index,
                current_row, report_context.categories,
                ordinary_witness_keys=report_context.ordinary_occurrence_witnesses(),
            )
        except (ValueError, TypeError, KeyError, AttributeError) as error:
            raise ReportProvenanceError(
                "Possible reconstruction source replay failed: {}".format(error)
            ) from error
    if method == "source_backed_candidate_reconstruction":
        from components.reconstruction_transaction import validate_reconstruction_basis
        report_context = _validation_context(
            source_sequence, expected_source_sha256, validation_context,
        )
        try:
            return validate_reconstruction_basis(
                source_sequence, candidate, basis, row_score, category_index,
                current_row, report_context.categories,
                ordinary_witness_keys=report_context.ordinary_occurrence_witnesses(),
            )
        except (ValueError, TypeError, KeyError, AttributeError) as error:
            raise ReportProvenanceError(
                "candidate reconstruction source replay failed: {}".format(error)
            ) from error
    if method == INTERNAL_REPEAT_COMPLETION_METHOD:
        report_context = _validation_context(source_sequence, expected_source_sha256, validation_context)
        if "source" not in report_context._internal_replay_context:
            report_context._internal_replay_context["source"] = InternalSearchContext(source_sequence)
        internal_context = report_context._internal_replay_context["source"]
        retained = basis.get("retained_core", False)
        required = {"method", "internal_repeat_completion"}
        if retained is True:
            required.update(("retained_core", "prior_classification_basis"))
        _keys(basis, frozenset(required), "classification_basis")
        try:
            completion = basis["internal_repeat_completion"]
            if "local_core_composition" in completion:
                validate_declared_family_against_run(completion["local_core_composition"],
                    tuple(map(tuple, completion["completed_geometry"])),
                    report_context._internal_completion_evidence)
            if retained is True:
                replay_retained_internal_core(source_sequence, candidate,
                    basis["internal_repeat_completion"], row_score, category_index, internal_context)
                validate_classification_basis(
                    source_sequence, expected_source_sha256, candidate,
                    basis["prior_classification_basis"], row_score=row_score,
                    allowed_methods=allowed_methods, depth=depth + 1,
                    validation_context=validation_context,
                    category_index=basis["internal_repeat_completion"]["core_category_index"],
                )
                return basis
            replay_internal_completion_basis(
                source_sequence, candidate, basis["internal_repeat_completion"],
                row_score, category_index, internal_context,
            )
        except (ValueError, TypeError, KeyError, AttributeError) as error:
            raise ReportProvenanceError("internal completion source replay failed: {}".format(error)) from error
        return basis
    if method == TERMINAL_REPAIR_COMPLETION_METHOD:
        _keys(basis, frozenset(("method", "terminal_repair_completion")), "classification_basis")
        try:
            replay_terminal_completion_basis(
                source_sequence, candidate,
                basis["terminal_repair_completion"], row_score, category_index,
            )
        except (ValueError, TypeError, KeyError) as error:
            raise ReportProvenanceError("terminal completion source replay failed: {}".format(error)) from error
        return basis
    if method in SOURCE_METHODS:
        context = _validation_context(
            source_sequence,
            expected_source_sha256,
            validation_context,
        )
        source = context.source_sequence

    if method == LIKELY_TANDEM_METHOD:
        return _validate_likely_tandem_classification_basis(
            source,
            expected_source_sha256,
            candidate,
            basis,
            row_score=row_score,
            allowed_methods=allowed_methods,
            depth=depth,
            validation_context=context,
            category_index=category_index,
            boundary_assessment=boundary_assessment,
        )

    if method == SOURCE_PERIODIC_LATTICE_METHOD:
        return _validate_periodic_classification_basis(
            source,
            expected_source_sha256,
            candidate,
            basis,
            row_score=row_score,
            allowed_methods=allowed_methods,
            depth=depth,
            validation_context=context,
            category_index=category_index,
            boundary_assessment=boundary_assessment,
        )

    if method == SOURCE_ARCHITECTURE_RESCUE_METHOD:
        return _validate_architecture_rescue_classification_basis(
            source,
            expected_source_sha256,
            candidate,
            basis,
            row_score=row_score,
            category_index=category_index,
            boundary_assessment=boundary_assessment,
        )

    if method == SOURCE_REGULAR_SUBPATH_METHOD:
        return _validate_regular_subpath_classification_basis(
            source,
            expected_source_sha256,
            candidate,
            basis,
            row_score=row_score,
            allowed_methods=allowed_methods,
            depth=depth,
            validation_context=context,
            category_index=category_index,
            boundary_assessment=boundary_assessment,
            current_row=current_row,
        )

    if method == SOURCE_ACCEPTED_FAMILY_METHOD:
        if category_index != 2:
            raise ReportProvenanceError(
                "accepted-family rescue may only classify Possible rows"
            )
        _keys(
            basis,
            frozenset((
                "method",
                "classification_ceiling",
                "boundary_status",
                "exact_boundary_certified",
                "accepted_family_rescue",
            )),
            "classification_basis",
        )
        if (
                basis["classification_ceiling"] != "Possible"
                or basis["boundary_status"] != "AMBIGUOUS"
                or basis["exact_boundary_certified"] is not False):
            raise ReportProvenanceError(
                "accepted-family rescue overstates classification/boundary"
            )
        if row_score is None:
            raise ReportProvenanceError(
                "accepted-family rescue lacks a reported fresh score"
            )
        try:
            replay_accepted_family_evidence(
                source,
                candidate,
                _mapping(
                    basis["accepted_family_rescue"],
                    "classification_basis.accepted_family_rescue",
                ),
                row_score,
            )
        except (TypeError, ValueError) as error:
            raise ReportProvenanceError(
                "accepted-family rescue source replay failed: {}".format(
                    error
                )
            ) from error
        return basis

    prior = basis.get("prior_classification_basis")
    if prior is not None:
        validate_classification_basis(
            source_sequence,
            expected_source_sha256,
            candidate,
            _mapping(prior, "prior_classification_basis"),
            row_score=row_score,
            allowed_methods=allowed_methods,
            depth=depth + 1,
            validation_context=context,
            category_index=category_index,
        )
    if method not in SOURCE_METHODS:
        return basis

    current = _geometry_from_candidate(source, candidate)
    if method == SOURCE_RESCUE_METHOD:
        if category_index != 2:
            raise ReportProvenanceError(
                "source-intact rescue may only classify Possible rows"
            )
        _keys(
            basis,
            frozenset((
                "method",
                "prior_classification_basis",
                "source_intact_recovery",
                "boundary_status",
                "exact_boundary_certified",
            )),
            "classification_basis",
        )
        if (
            basis["boundary_status"] != "AMBIGUOUS"
            or basis["exact_boundary_certified"] is not False
        ):
            raise ReportProvenanceError(
                "source-intact rescue overstates its boundary status"
            )
    else:
        if category_index != 3:
            raise ReportProvenanceError(
                "source-intact retirement may only classify discarded rows"
            )
        _keys(
            basis,
            frozenset((
                "method",
                "prior_classification_basis",
                "source_intact_recovery",
                "transition",
            )),
            "classification_basis",
        )
        transition = _mapping(basis["transition"], "classification_basis.transition")
        expected_transition = {
            "source_category_index": 2,
            "source_category_name": "Possible",
            "target_category_index": 3,
            "target_category_name": "Possible Discarded",
        }
        if transition != expected_transition:
            raise ReportProvenanceError(
                "terminal-conflict retirement transition is stale"
            )

    record = _mapping(
        basis["source_intact_recovery"],
        "classification_basis.source_intact_recovery",
    )
    recovery, targets, _conflicts, _alternatives = _validate_recovery_record(
        source,
        expected_source_sha256,
        record,
        "classification_basis.source_intact_recovery",
        context,
    )
    if method == SOURCE_RESCUE_METHOD:
        if record["action"] not in (
            "admit_possible", "replace_conflicting_possible"
        ):
            raise ReportProvenanceError(
                "source-intact rescue record has the wrong action"
            )
        if not _same_geometry(current, recovery):
            raise ReportProvenanceError(
                "source-intact rescue geometry differs from reported candidate"
            )
        if row_score is not None and float(record["score"]) != float(row_score):
            raise ReportProvenanceError(
                "source-intact rescue score differs from reported row"
            )
        consensus = getattr(candidate, "consensus", None)
        if (
            not isinstance(consensus, str)
            or len(consensus)
            != record["evidence"]["existence"]["metrics"]["consensus_length"]
        ):
            raise ReportProvenanceError(
                "source-intact rescue consensus claim is stale"
            )
    else:
        if record["action"] != "replace_conflicting_possible":
            raise ReportProvenanceError(
                "terminal-conflict retirement record has the wrong action"
            )
        if not any(_same_geometry(current, target) for target in targets):
            raise ReportProvenanceError(
                "retired candidate is absent from accepted targets"
            )
    return basis


def _boundary_input(source, source_hash, geometry, label, dimensions):
    return AlternativeRelationshipInput(
        geometry=CandidateGeometryInput(
            source_sha256=source_hash,
            source_interval=geometry.interval,
            repeat_intervals=geometry.repeat_intervals,
            repeat_sequences=tuple(
                source[start:end] for start, end in geometry.repeat_intervals
            ),
            spacer_intervals=geometry.spacer_intervals,
            spacer_sequences=tuple(
                source[start:end] for start, end in geometry.spacer_intervals
            ),
        ),
        label=label,
        dimensions=dimensions,
        reason_codes=("source_intact_recovery_relationship_revalidated",),
    )


def validate_source_intact_recovery_alternatives(
    source_sequence,
    expected_source_sha256,
    candidate,
    raw_records,
    *,
    validation_context=None,
    category_index=None,
):
    """Validate attached recovery records and derive only proven alternatives."""

    if category_index not in (0, 2):
        raise ReportProvenanceError(
            "source-intact alternatives require an accepted report row"
        )
    context = _validation_context(
        source_sequence,
        expected_source_sha256,
        validation_context,
    )
    source = context.source_sequence
    records = _sequence(
        raw_records, "source_intact_recovery_alternatives"
    )
    current = _geometry_from_candidate(source, candidate)
    relationships = []
    canonical = []
    for index, value in enumerate(records):
        label = "source_intact_recovery_alternatives[{}]".format(index)
        record = _mapping(value, label)
        recovery, targets, conflicts, explicit = _validate_recovery_record(
            source,
            expected_source_sha256,
            record,
            label,
            context,
        )
        if record["action"] != "attach_ambiguous_alternative":
            raise ReportProvenanceError(
                label + " is not an attached-alternative action"
            )
        matching_target_indexes = tuple(
            target_index
            for target_index, target in enumerate(targets)
            if _same_geometry(current, target)
        )
        if not matching_target_indexes:
            raise ReportProvenanceError(
                label + " does not identify the reported accepted anchor"
            )
        if not any(
            record["accepted_targets"][target_index]["category_index"]
            == category_index
            for target_index in matching_target_indexes
        ):
            raise ReportProvenanceError(
                label + " identifies the accepted anchor in another category"
            )

        relationship_proven = False
        for accepted, relationship in conflicts:
            if (
                accepted.interval == current.interval
                and accepted.repeat_intervals == current.repeat_intervals
                and relationship["topology"] in (
                    "source_intact_contains_accepted",
                    "accepted_contains_source_intact",
                )
                and (
                    relationship["accepted_exact_repeat_subpath"]
                    or relationship["source_intact_exact_repeat_subpath"]
                )
            ):
                relationships.append(_boundary_input(
                    source,
                    expected_source_sha256,
                    recovery,
                    "source_intact_exact_repeat_subpath",
                    (DIMENSION_OCCURRENCE_PATH,),
                ))
                relationship_proven = True
        for first, second in explicit:
            if _same_geometry(first, current):
                alternative = second
            elif _same_geometry(second, current):
                alternative = first
            else:
                continue
            first_set = set(current.repeat_intervals)
            second_set = set(alternative.repeat_intervals)
            if first_set.issubset(second_set) or second_set.issubset(first_set):
                dimensions = (DIMENSION_OCCURRENCE_PATH,)
                relation_label = "source_intact_exact_repeat_subpath"
            elif len(first_set) == len(second_set):
                dimensions = (DIMENSION_REPEAT_PHASE,)
                relation_label = "source_intact_alternate_repeat_phase"
            else:
                dimensions = (
                    DIMENSION_OCCURRENCE_PATH,
                    DIMENSION_REPEAT_PHASE,
                )
                relation_label = "source_intact_overlapping_boundary_alternative"
            relationships.append(_boundary_input(
                source,
                expected_source_sha256,
                alternative,
                relation_label,
                dimensions,
            ))
            relationship_proven = True
        if not relationship_proven:
            raise ReportProvenanceError(
                label + " lacks a source-revalidated same-locus relationship"
            )
        canonical.append(record)

    unique = {}
    for relationship in relationships:
        key = (
            relationship.geometry.repeat_intervals,
            relationship.label,
            relationship.dimensions,
        )
        unique[key] = relationship
    ordered_relationships = tuple(
        unique[key] for key in sorted(unique, key=lambda item: repr(item))
    )
    return ValidatedSourceIntactAlternatives(
        records=tuple(canonical),
        boundary_relationships=ordered_relationships,
    )
