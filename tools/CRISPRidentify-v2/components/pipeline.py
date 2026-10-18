from collections import Counter, deque
import json
import math
import os

from components.accepted_core_boundary_attestation import (
    BoundaryPhaseCompletenessEvidence,
    BoundaryPhaseLineageEntry,
    KNOWN_LINEAGE_DERIVATIONS,
    attest_accepted_core_boundary_phase,
)
from components.accepted_core_boundary_completion import (
    BoundaryGeometryWitness,
    CoreAuthorizationEvidence,
    FinalBoundaryEvidence,
    SourceGeometry,
    WitnessContextEvidence,
)
from components.accepted_core_boundary_transaction import (
    AcceptedCoreBoundaryContext,
    BoundaryCompletionScoreResult,
    apply_accepted_core_boundary_completion_batch,
    plan_accepted_core_boundary_completion_batch,
)
from components.accepted_family_tandem_transaction import (
    AcceptedRowFamilyContext,
    apply_accepted_family_tandem_demotions,
)
from components.accepted_family_tandem_validator import (
    validate_accepted_family_tandem,
)
from components.circular_array_join_integration import (
    evaluate_accepted_circular_edge_joins,
)
from components.likely_tandem_transaction import (
    apply_likely_tandem_quarantine,
)
from components.reconstruction_transaction import reconstruct_accepted_arrays
from components.possible_reconstruction import reconstruct_possible_arrays
from components.contextual_spacer_rescue import rescue_contextual_spacer_arrays
from components.architecture_rescue_coding_context import (
    batch_source_coding_context_evidence,
    predict_source_cds,
)
from components.architecture_rescue_integration import (
    arbitrate_cross_family_phase_proposals,
    candidate_repeat_path,
    evaluate_scored_architecture_candidate,
    plan_architecture_rescue_transitions,
    select_augmentation_seed_ids,
    select_lineage_family_proposals,
)
from components.boundary_assessment import BoundaryAssessmentNormalizer
from components.module_detection import Detection
from components.module_detection_refinement import DetectionRefinement
from components.repeat_boundary_repair import POLICY_VERSION, source_geometry
from components.internal_repeat_completion import InternalSearchContext, POLICY as INTERNAL_POLICY
from components.internal_completion_authority import capture_current_authority_inventory
from components.internal_repeat_transaction import complete_internal_repeats
from components.terminal_repair_completion import (
    plan_terminal_repair_completion,
    apply_terminal_repair_completion,
)
from components.module_evaluation import (
    ArrayEvaluation,
    build_advanced_fuzzy_search_filter,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
    collapse_source_backed_phase_equivalent_candidates,
    source_alignment_cache_info,
)
from components.module_evaluated_arrays_enhancement import EvaluatedArraysEnhancement
from components.module_evaluated_arrays_enhancement import (
    build_source_backed_terminal_variants,
)
from components.module_non_array_computations import NonArrayComputations
from components.module_output_maker import OutputMaker
from components.folder_scheduler import pipeline_result_key
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.fast_run_rescue import candidate_score_fingerprint
from components.fast_run_rescue import candidate_source_interval
from components.fast_run_rescue import compact_feature_values
from components.fast_run_rescue import select_unresolved_components
from components.fast_two_repeat_terminal_repair import (
    MAXIMUM_ADMISSION_REPEAT_DERIVED_SPACER_COVERAGE,
    MINIMUM_ADMISSION_MODEL_SCORE,
    fast_two_repeat_admission_allowed,
    flank_continuation_identities,
    generate_source_winners,
)
from components.fast_accepted_family_rescue import (
    classification_basis as accepted_family_classification_basis,
    completed_score_is_admissible as accepted_family_score_is_admissible,
    generate_accepted_family_rescue_proposals,
    select_exact_template_raw_seed_parents,
)
from components.reconstruction_validator import (
    read_single_fasta,
    validate_candidate_reconstruction,
)
from components.report_provenance_validation import (
    ReportProvenanceError,
    canonical_source_architecture_rescue_boundary_assessment,
)
from components.sequence_topology import (
    to_circular_topology_evidence,
    unknown_sequence_topology,
)
from components.source_superlocus_augmentation import (
    _validated_source_path_family_keys,
    augment_source_superloci_to_fixed_point,
    replay_source_superlocus_path_local_convergence_proof,
)
from components.source_intact_recovery_transaction import (
    SourceBackedRowFamilyEvidence,
    SourceIntactRecoveryCandidateContext,
    SourceIntactRecoveryLocusContext,
    apply_source_intact_recovery_plan,
    plan_source_intact_recoveries,
    validate_source_intact_alternative_relationship,
)
from components.source_periodic_lattice_category_transaction import (
    SourcePeriodicLatticeScoreResult,
    apply_source_periodic_lattice_repairs,
    plan_source_periodic_lattice_repairs,
)
from components.source_regular_subpath_recovery import (
    DIAGNOSTIC_PARENT,
    SourceRegularSubpathContext,
    propose_source_regular_subpaths_in_context,
)
from components.source_regular_subpath_recovery_transaction import (
    RegularSubpathFreshScoreResult,
    RegularSubpathProposalTandemEvidence,
    SourceRegularSubpathDiagnosticRowContext,
    apply_source_regular_subpath_recovery_batch,
    plan_source_regular_subpath_recovery_batch,
)
CATEGORY_NAMES = (
    "Bona-fide",
    "Alternative",
    "Possible",
    "Possible Discarded",
    "Low score",
)
MAX_FAST_TWO_REPEAT_ROUTED_PARENTS_PER_SOURCE = 8
_SOURCE_PHASE_COLLAPSE_COUNT_FIELDS = (
    "groups_considered",
    "groups_collapsed",
    "groups_ambiguous",
    "groups_outside_phase_rule",
    "groups_weak_boundary_evidence",
    "groups_weak_score_margin",
    "candidates_removed",
    "candidates_input",
    "candidates_output",
)


def _empty_source_phase_collapse_stats():
    return dict(
        {field: 0 for field in _SOURCE_PHASE_COLLAPSE_COUNT_FIELDS},
        enabled=False,
    )


class _CandidateScoreSnapshot:
    __slots__ = ("candidate", "entry", "row")

    def __init__(self, entry):
        self.entry = entry
        self.candidate = entry[1]
        self.row = tuple(entry)


class CandidateScoreCache:
    """Reuse scores only while the exact scored entry remains unchanged.

    Enhancement replaces the candidate at index one of an existing score row.
    The row object's identity is therefore a stable provenance token even when
    structural filtering moves the row to another category. Equality with the
    original candidate is an independent fence: provenance alone never permits
    reuse after sequence, alignment, spacer, or coordinate changes.
    """

    def __init__(self, categories):
        self._snapshots = {}
        self.reused_count = 0
        self.recomputed_count = 0
        for _category_index, _interval, _entry_index, entry in self._entries(
                categories):
            self._snapshots.setdefault(
                id(entry), _CandidateScoreSnapshot(entry)
            )

    @staticmethod
    def _entries(categories):
        for category_index, category in enumerate(categories):
            for interval in sorted(category):
                for entry_index, entry in enumerate(category[interval]):
                    if not isinstance(entry, (tuple, list)) or len(entry) < 2:
                        raise RuntimeError(
                            "cannot cache a malformed scored candidate entry"
                        )
                    yield category_index, interval, entry_index, entry

    @staticmethod
    def _candidate_is_unchanged(snapshot, entry):
        return (
            snapshot is not None
            and snapshot.entry is entry
            and snapshot.candidate == entry[1]
        )

    def merge_after_enhancement(self, categories, score_changed_candidates):
        """Return all score rows in legacy traversal order.

        ``score_changed_candidates`` receives only candidates that failed the
        provenance-plus-structural-equality fence. It must return the ordinary
        ``{interval: [score rows...]}`` representation in candidate order.
        """

        slots = []
        changed_candidates = {}
        reused_count = 0
        recomputed_count = 0

        for _category_index, interval, _entry_index, entry in self._entries(
                categories):
            candidate = entry[1]
            snapshot = self._snapshots.get(id(entry))
            if self._candidate_is_unchanged(snapshot, entry):
                cached_row = list(snapshot.row)
                # Preserve the final post-enhancement object even when it is a
                # newly constructed but structurally equal candidate.
                cached_row[1] = candidate
                slots.append((interval, candidate, cached_row))
                reused_count += 1
            else:
                changed_candidates.setdefault(interval, []).append(candidate)
                slots.append((interval, candidate, None))
                recomputed_count += 1

        if changed_candidates:
            fresh_rows = score_changed_candidates(changed_candidates)
        else:
            fresh_rows = {}

        unexpected_intervals = set(fresh_rows).difference(changed_candidates)
        if unexpected_intervals:
            raise RuntimeError(
                "fresh evaluation returned unexpected intervals: {}".format(
                    sorted(unexpected_intervals)
                )
            )
        fresh_queues = {
            interval: deque(fresh_rows.get(interval, ()))
            for interval in changed_candidates
        }

        merged_rows = {}
        for interval, candidate, cached_row in slots:
            if cached_row is None:
                queue = fresh_queues[interval]
                if not queue:
                    raise RuntimeError(
                        "fresh evaluation omitted a changed candidate"
                    )
                row = list(queue.popleft())
                if len(row) < 2 or row[1] is not candidate:
                    raise RuntimeError(
                        "fresh evaluation changed candidate order or identity"
                    )
                row[1] = candidate
            else:
                row = cached_row
            merged_rows.setdefault(interval, []).append(row)

        if any(fresh_queues.values()):
            raise RuntimeError("fresh evaluation returned surplus score rows")

        self.reused_count = reused_count
        self.recomputed_count = recomputed_count
        return merged_rows


def retain_requested_categories(categories, accepted_only):
    """Apply the output-only accepted-category projection in place."""

    if accepted_only:
        for category_index in (1, 3, 4):
            categories[category_index] = {}
    return categories


class CandidateReconstructionError(RuntimeError):
    """A candidate failed a sequence/coordinate correctness fence."""

    def __init__(self, stage, category, interval, candidate_index, result):
        self.stage = stage
        self.result = result
        issue_codes = ", ".join(issue.code for issue in result.issues)
        message = (
            "candidate reconstruction failed at {}: category={!r}, "
            "interval={!r}, candidate_index={}, issues={}"
        ).format(stage, category, interval, candidate_index, issue_codes)
        super().__init__(message)



class Pipeline:
    def __init__(self, result_folder_path, pickle_folder_path, json_folder_path, file_path,
                 parameters, flags, flag_dev_mode,
                 sequence_topology_manifest=None, retain_candidate_evidence=False):
        self.retain_candidate_evidence = retain_candidate_evidence
        self._candidate_evidence_categories = None
        self.result_folder_path = (
            result_folder_path + "/" + pipeline_result_key(file_path)
        )
        try:
            os.unlink(os.path.join(self.result_folder_path, "report.json"))
        except FileNotFoundError:
            pass
        self.pickle_folder_path = pickle_folder_path
        self.json_folder_path = json_folder_path
        self.file_path = file_path
        self.list_features = [list(COMPACT_FEATURE_SCHEMA)]
        self.flags = flags
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode
        self.sequence_topology_manifest = sequence_topology_manifest

        self.header = None
        self.dict_fuzzy_crisprs = {}
        self.genome_sequence = None
        self.dict_crispr_candidates = {}
        self.categories = {}
        self.non_array_data = {}
        self.sequence_topology_evidence = None
        self.circular_array_join_result = None
        self.likely_tandem_records = ()
        self.likely_tandem_stats = None
        self.likely_tandem_source_categories = None
        self.reconstruction_evidence = None
        self.possible_reconstruction_evidence = None
        self.contextual_spacer_evidence = None
        self.pending_contextual_spacers = []
        self._candidate_score_cache = None
        self.reevaluation_stats = {"reused": 0, "recomputed": 0}
        self.source_phase_collapse_provenance = ()
        self.repeat_boundary_repair_records = []
        self.repeat_boundary_repair_stats = {}
        self.source_enhancement_lineage = ()
        self.source_phase_collapse_stats = (
            _empty_source_phase_collapse_stats()
        )
        self.source_phase_collapse_stage_stats = {}
        self.source_terminal_completion_stats = {
            "enabled": False,
            "variants_scored": 0,
            "variants_committed": 0,
            "variants_rejected": 0,
            "committed_repeats": 0,
        }
        self.source_architecture_rescue_lineage = ()
        self.source_architecture_rescue_transitions = ()
        self.source_architecture_phase_arbitration = ()
        self.source_accepted_family_tandem_records = ()
        self._source_accepted_family_tandem_by_candidate_id = {}
        self.source_accepted_family_tandem_stats = {
            "enabled": False,
            "accepted_rows_reviewed": 0,
            "tandem_supported": 0,
        }
        self.source_accepted_core_boundary_records = ()
        self.source_accepted_core_boundary_phase_attestations = ()
        self._source_accepted_core_boundary_family_records = ()
        self.source_accepted_core_boundary_stats = {
            "enabled": False,
            "accepted_rows_planned": 0,
            "scoring_requests": 0,
            "completed_rows": 0,
        }
        self.source_intact_recovery_records = ()
        self.source_intact_recovery_stats = {
            "enabled": False,
            "reviewed": 0,
            "admitted_possible": 0,
            "replaced_possible": 0,
            "alternatives_attached": 0,
            "abstained": 0,
            "rejected": 0,
            "lineage_only_excluded": 0,
        }
        self.source_intact_recovery_cost_audit = None
        self._source_intact_recovery_contexts_snapshot = ()
        self.source_periodic_lattice_repair_records = ()
        self.source_periodic_lattice_repair_stats = {
            "enabled": False,
            "accepted_rows_reviewed": 0,
            "scoring_requests": 0,
            "repaired_rows": 0,
        }
        self.source_periodic_lattice_repair_cost_audit = None
        self.source_regular_subpath_recovery_records = ()
        self.source_regular_subpath_recovery_stats = {
            "enabled": False,
            "parents_reviewed": 0,
            "structural_proposals": 0,
            "scoring_requests": 0,
            "possible_admissions": 0,
            "bona_fide_admissions": 0,
        }
        self.source_regular_subpath_recovery_batch_audit = None
        self.source_regular_subpath_recovery_cost_audit = None
        self._source_intact_lineage_only_excluded = 0
        self._source_intact_context_tandem_evaluations = 0
        self._source_architecture_baseline_rows = {}
        self._source_architecture_augmentation_stats = {}
        self._source_cds = None
        self.source_architecture_rescue_stats = {
            "enabled": False,
            "seed_count": 0,
            "primary_seed_count": 0,
            "satellite_seed_count": 0,
            "generated_candidate_count": 0,
            "generated_candidates_scored": 0,
            "context_candidates_checked": 0,
            "proposals_considered": 0,
            "proposals_after_phase_arbitration": 0,
            "phase_conflicts_selected": 0,
            "phase_conflicts_ambiguous_representative": 0,
            "phase_conflicts_quarantined": 0,
            "phase_proposals_removed": 0,
            "calls_added": 0,
            "calls_replaced": 0,
            "calls_abstained": 0,
        }
        self._post_enhancement_scored_rows = {}
        self._final_scored_rows = {}
        self.fast_run_rescue_stats = {
            "enabled": False,
            "triggered_components": 0,
            "selected_clusters": 0,
            "delta_patterns": 0,
            "fresh_scores": 0,
            "reused_scores": 0,
        }
        self.fast_run_phase_repair_stats = {
            "enabled": False,
            "promotion_parents": 0,
            "admission_parents": 0,
            "variants_scored": 0,
            "promotions_committed": 0,
            "admissions_committed": 0,
        }
        self.fast_two_repeat_terminal_repair_stats = {
            "enabled": False,
            "seeds_retained": 0,
            "source_winners": 0,
            "fresh_scores": 0,
            "ordinary_score_passes": 0,
            "admissions_committed": 0,
        }
        self.fast_accepted_family_rescue_stats = {
            "enabled": False,
            "proposals": 0,
            "fresh_scores": 0,
            "score_rejections": 0,
            "admissions_committed": 0,
        }
        self._fast_accepted_family_anchor_categories = None
        self._raw_detection = None
        self._initial_scored_rows = {}

        self._get_header()
        self._run_detection()
        self._run_detection_refinement()
        self._validate_refined_candidates()
        self._collapse_source_backed_phase_equivalent_candidates()
        self._run_evaluation()
        self._run_fast_routed_start_end_rescue()
        # Rescue has transferred every live candidate into scored categories.
        # Rebind owning references: the raw detector still shares some maps and
        # must retain its rejected witnesses for the later fast repair stages.
        self.dict_fuzzy_crisprs = {}
        self.dict_crispr_candidates = {}
        self._initial_scored_rows = {}
        self._results_enhancement()
        self._validate_candidate_categories("post-enhancement")
        self._normalize_post_enhancement_source_phase()
        self._validate_candidate_categories(
            "post-enhancement-source-phase-normalization"
        )
        self._reevaluate_enhanced_candidates(
            project_requested_categories=False
        )
        self._validate_candidate_categories("post-reevaluation")
        self._candidate_score_cache = None
        if self.flags.get("flag_fast_accepted_family_rescue", False):
            self._fast_accepted_family_anchor_categories = (
                self._copy_categories(self.categories)
            )
        self._run_fast_two_repeat_terminal_repair()
        self._validate_candidate_categories(
            "post-fast-two-repeat-terminal-repair"
        )
        self._run_fast_accepted_family_rescue()
        self._validate_candidate_categories(
            "post-fast-accepted-family-rescue"
        )
        # Both seed consumers are finished. Keep their copied audit records,
        # but release the raw fuzzy inventory and the pre-family snapshot.
        self._raw_detection = None
        self._fast_accepted_family_anchor_categories = None
        self._run_fast_bounded_phase_repair()
        self._validate_candidate_categories("post-fast-phase-repair")
        self._complete_attested_terminal_repairs()
        self._validate_candidate_categories("post-terminal-repair-completion")
        self._complete_pending_internal_repeats()
        self._validate_candidate_categories("post-internal-repeat-completion")
        if self._run_configured_source_backed_augmentation():
            self._validate_candidate_categories("post-architecture-rescue")
        else:
            self._validate_candidate_categories(
                "post-source-augmentation-disabled"
            )
        self._annotate_model_score_classification_bases()
        self._demote_accepted_family_tandems_transactionally()
        self._validate_candidate_categories(
            "post-accepted-family-tandem-demotion"
        )
        self._complete_accepted_core_boundaries_transactionally()
        self._validate_candidate_categories(
            "post-accepted-core-boundary-completion"
        )
        self._recover_source_intact_candidates_transactionally()
        self._validate_candidate_categories(
            "post-source-intact-recovery"
        )
        self._repair_source_periodic_lattices_transactionally()
        self._validate_candidate_categories(
            "post-source-periodic-lattice-repair"
        )
        self._recover_source_regular_subpaths_transactionally()
        self._validate_candidate_categories(
            "post-source-regular-subpath-recovery"
        )
        self._reconstruct_accepted_arrays()
        self._validate_candidate_categories("post-candidate-reconstruction")
        # Subsequent assessment/publication consumes categories and immutable
        # transaction evidence, not these historical scoring workspaces.
        self._candidate_score_cache = None
        self._post_enhancement_scored_rows = {}
        self._final_scored_rows = {}
        self._quarantine_likely_tandems_transactionally()
        self._validate_candidate_categories(
            "post-likely-tandem-quarantine"
        )
        self._evaluate_circular_origin_joins()
        if self.retain_candidate_evidence:
            # Snapshot containers before accepted-only projection clears them.
            # Candidate geometry/scores remain unmodified by this observation.
            self._candidate_evidence_categories = tuple(
                {key: tuple(rows) for key, rows in category.items()}
                for category in self.categories
            )
        self._drop_nonaccepted_categories()
        self._validate_candidate_categories("post-output-projection")
        self._run_non_crispr_computation()
        self._validate_and_write_output()

    def _get_header(self):
        header, self.genome_sequence = read_single_fasta(self.file_path)
        self.header = ">{}\n".format(header)
        source_id = header.split()[0]
        topology_manifest = getattr(
            self, "sequence_topology_manifest", None
        )
        if topology_manifest is None:
            self.sequence_topology_evidence = unknown_sequence_topology(
                source_id,
                self.genome_sequence,
            )
        else:
            self.sequence_topology_evidence = (
                topology_manifest.evidence_for(
                    source_id,
                    self.genome_sequence,
                )
            )

    def _reconstruct_accepted_arrays(self):
        """Keep core-supported reconstruction outside ordinary score competition."""
        policy = self.parameters.get("param_reconstruction_policy", "legacy")
        if policy == "legacy":
            return
        if policy != "conservative":
            raise ValueError("unknown reconstruction policy: {}".format(policy))

        def source_cds():
            if self._source_cds is None:
                self._source_cds = predict_source_cds(self.genome_sequence)
            return self._source_cds

        self.categories, self.reconstruction_evidence = reconstruct_accepted_arrays(
            self.genome_sequence, self.categories, self.parameters,
            self._score_candidates_complete, source_cds,
        )
        # The lower-confidence lane cannot change the completed Bona-fide
        # transaction. Its proposals must respect those accepted extents.
        self.categories, self.possible_reconstruction_evidence = reconstruct_possible_arrays(
            self.genome_sequence, self.categories, self.parameters,
            self._score_candidates_complete, source_cds,
        )
        self.categories, self.contextual_spacer_evidence = rescue_contextual_spacer_arrays(
            self.genome_sequence, self.categories,
            getattr(self, "pending_contextual_spacers", ()), self.parameters,
            self._score_candidates_complete, source_cds,
        )
        self.pending_contextual_spacers = []
        self._final_scored_rows = self._scored_rows_from_categories(self.categories)

    def _run_detection(self):
        print("1. Run initial array detection")
        detection = Detection(file_path=self.file_path,
                              flags=self.flags,
                              parameters=self.parameters,
                              flag_dev_mode=self.flag_dev_mode)
        if (
                self.flags.get("flag_fast_run", False)
                and (
                    self.flags.get("flag_fast_run_rescue", False)
                    or self.flags.get(
                        "flag_fast_two_repeat_terminal_repair", False
                    )
                    or self.flags.get(
                        "flag_fast_accepted_family_rescue", False
                    )
                )):
            self._raw_detection = detection
        self.dict_fuzzy_crisprs = detection.output()

    def _run_detection_refinement(self):
        print("2. Refine detected arrays")
        self.internal_search_context = InternalSearchContext(self.genome_sequence)
        det_ref = DetectionRefinement(dict_fuzzy_crisprs=self.dict_fuzzy_crisprs,
                                      parameters=self.parameters,
                                      flag_dev_mode=self.flag_dev_mode,
                                      source_dna=self.genome_sequence,
                                      capture_contextual_spacers=(
                                          self.parameters.get("param_reconstruction_policy")
                                          == "conservative"
                                      ),
                                      internal_search_context=self.internal_search_context)
        self.dict_crispr_candidates = det_ref.output()
        self._record_repeat_boundary_repair("initial", det_ref)

    def _record_repeat_boundary_repair(self, stage, refinement):
        contextual = getattr(self, "pending_contextual_spacers", [])
        contextual.extend(getattr(refinement, "pending_contextual_spacers", ()))
        self.pending_contextual_spacers = contextual
        pending = getattr(self, "pending_internal_repairs", [])
        pending.extend(getattr(refinement, "pending_internal_repairs", ()))
        self.pending_internal_repairs = pending
        internal_records = getattr(self, "internal_repair_records", [])
        internal_records.extend(getattr(refinement, "internal_repair_records", ()))
        self.internal_repair_records = internal_records
        records = getattr(self, "repeat_boundary_repair_records", [])
        records.extend(dict(record, stage=stage)
                       for record in refinement.boundary_repair_records)
        self.repeat_boundary_repair_records = records
        stats = getattr(self, "repeat_boundary_repair_stats", {})
        for name, count in refinement.boundary_repair_stats.items():
            stats[name] = stats.get(name, 0) + count
        self.repeat_boundary_repair_stats = stats
        snapshots = getattr(self, "repeat_boundary_repair_cache_stage_snapshots", [])
        snapshots.append({"stage": stage, "caches": source_alignment_cache_info()})
        self.repeat_boundary_repair_cache_stage_snapshots = snapshots

    def _repeat_boundary_repair_audit(self):
        """Retain proposal lineage independently of later object rebuilds."""
        outcomes = {}
        for category_index, category in enumerate(self.categories):
            for rows in category.values():
                for row in rows:
                    outcomes.setdefault(source_geometry(row[1]), []).append({
                        "category": CATEGORY_NAMES[category_index],
                        "certainty_score": float(row[0]),
                    })
        records = []
        for original in getattr(self, "repeat_boundary_repair_records", ()):
            record = dict(original)
            geometry = record.get("child_geometry")
            if geometry is not None:
                record["final_geometry_outcomes"] = outcomes.get(geometry, [])
            records.append(record)
        return {"policy": POLICY_VERSION,
                "stats": getattr(self, "repeat_boundary_repair_stats", {}),
                "records": records,
                "source_alignment_caches": source_alignment_cache_info(),
                "cache_stage_snapshots": getattr(self, "repeat_boundary_repair_cache_stage_snapshots", []),
                "completion_transactions": getattr(self, "terminal_repair_completion_records", [])}

    def _complete_pending_internal_repeats(self):
        structural = build_advanced_fuzzy_search_filter(self.parameters)
        # Retain an immutable snapshot before the transaction can replace rows.
        # This is independent of the bridge proof and its filtered family.
        self.internal_current_authority_inventory_json = json.dumps(
            capture_current_authority_inventory(
                self.genome_sequence, self.categories,
                getattr(self, "pending_internal_repairs", ())), sort_keys=True)
        self.categories, self.internal_completion_transactions = complete_internal_repeats(
            self.genome_sequence, self.categories,
            getattr(self, "pending_internal_repairs", ()), structural,
            self._score_candidates_complete,
            search_context=getattr(self, "internal_search_context", None),
            boundary_records=getattr(self, "repeat_boundary_repair_records", None),
        )
        self._final_scored_rows = self._scored_rows_from_categories(self.categories)
        self._post_enhancement_scored_rows = self._final_scored_rows

    def _internal_repeat_audit(self):
        context = getattr(self, "internal_search_context", None)
        return {"policy": INTERNAL_POLICY,
                "stats": dict(context.stats) if context is not None else {},
                "records": getattr(self, "internal_repair_records", []),
                "current_authority_inventory": json.loads(getattr(
                    self, "internal_current_authority_inventory_json", "null")),
                "completion_transactions": getattr(self, "internal_completion_transactions", [])}

    def _complete_attested_terminal_repairs(self):
        structural = build_advanced_fuzzy_search_filter(self.parameters)
        plan = plan_terminal_repair_completion(
            self.genome_sequence, self.categories,
            getattr(self, "repeat_boundary_repair_records", ()), structural,
        )
        self.categories, self.terminal_repair_completion_records = (
            apply_terminal_repair_completion(
                self.genome_sequence, self.categories, plan, structural,
                self._score_candidates_complete,
            )
        )
        self._final_scored_rows = self._scored_rows_from_categories(self.categories)
        self._post_enhancement_scored_rows = self._final_scored_rows

    def _collapse_source_backed_phase_equivalent_candidates(self):
        """Choose phase coordinates from source evidence before scoring."""

        print("2b. Collapse source-backed phase-equivalent candidates")
        collapsed, provenance, stats = (
            collapse_source_backed_phase_equivalent_candidates(
                self.genome_sequence,
                self.dict_crispr_candidates,
                skip_source_scores_when_boundaries_reject=(
                    getattr(self, "flags", {}).get("output_profile")
                    == "machine"
                ),
            )
        )
        self.dict_crispr_candidates = collapsed
        self._record_source_phase_collapse(
            "pre_evaluation", provenance, stats
        )

    def _record_source_phase_collapse(self, stage, provenance, stats):
        """Append one named collapse stage and recompute aggregate counts."""

        stage_stats = getattr(
            self, "source_phase_collapse_stage_stats", None
        )
        if stage_stats is None:
            stage_stats = {}
            self.source_phase_collapse_stage_stats = stage_stats
        if stage in stage_stats:
            raise RuntimeError(
                "source phase-collapse stage recorded twice: {}".format(
                    stage
                )
            )

        labelled_stats = dict(stats)
        labelled_stats["stage"] = stage
        stage_stats[stage] = labelled_stats

        aggregate = _empty_source_phase_collapse_stats()
        for recorded_stats in stage_stats.values():
            aggregate["enabled"] = (
                aggregate["enabled"] or recorded_stats["enabled"]
            )
            for field in _SOURCE_PHASE_COLLAPSE_COUNT_FIELDS:
                aggregate[field] += recorded_stats[field]
        self.source_phase_collapse_stats = aggregate

        prior_provenance = tuple(getattr(
            self, "source_phase_collapse_provenance", ()
        ))
        labelled_provenance = tuple(
            dict(record, stage=stage) for record in provenance
        )
        self.source_phase_collapse_provenance = (
            prior_provenance + labelled_provenance
        )

    def _category_candidates_for_source_phase(self):
        """Return phase input plus exact category-row locations by locus."""

        if not isinstance(self.categories, (tuple, list)) or len(
                self.categories) != len(CATEGORY_NAMES):
            raise RuntimeError(
                "source phase normalization requires five categories"
            )

        candidates_by_locus = {}
        slots_by_locus = {}
        for category_index, category in enumerate(self.categories):
            if not isinstance(category, dict):
                raise RuntimeError(
                    "source phase normalization requires category mappings"
                )
            for locus in sorted(category):
                for entry_index, entry in enumerate(category[locus]):
                    if not isinstance(entry, (tuple, list)) or len(entry) < 2:
                        raise RuntimeError(
                            "cannot normalize a malformed scored candidate row"
                        )
                    candidate = entry[1]
                    candidates_by_locus.setdefault(locus, []).append(
                        candidate
                    )
                    slots_by_locus.setdefault(locus, []).append((
                        category_index,
                        locus,
                        entry_index,
                        entry,
                        candidate,
                    ))
        return candidates_by_locus, slots_by_locus

    def _normalize_post_enhancement_source_phase(self):
        """Collapse phase variants added by enhancement without stale scores."""

        print("4a. Normalize post-enhancement source phase")
        candidates_by_locus, slots_by_locus = (
            self._category_candidates_for_source_phase()
        )
        collapsed_by_locus, provenance, stats = (
            collapse_source_backed_phase_equivalent_candidates(
                self.genome_sequence,
                candidates_by_locus,
                skip_source_scores_when_boundaries_reject=(
                    getattr(self, "flags", {}).get("output_profile")
                    == "machine"
                ),
            )
        )

        actions_by_locus = {}
        enriched_provenance = []
        for raw_record in provenance:
            record = dict(raw_record)
            if record["status"] != "collapsed":
                enriched_provenance.append(record)
                continue

            locus = record["locus"]
            member_indexes = tuple(record["member_candidate_indexes"])
            output_index = record["output_candidate_index"]
            slots = slots_by_locus[locus]
            if (
                    output_index not in member_indexes
                    or len(member_indexes) != record["member_count"]
                    or len(set(member_indexes)) != len(member_indexes)
                    or any(index < 0 or index >= len(slots)
                           for index in member_indexes)):
                raise RuntimeError(
                    "invalid source phase-collapse row lineage"
                )

            actions = actions_by_locus.setdefault(locus, {})
            if any(index in actions for index in member_indexes):
                raise RuntimeError(
                    "overlapping source phase-collapse row lineage"
                )
            for index in member_indexes:
                actions[index] = None

            retained_entry = slots[output_index][3]
            canonical_entry = list(retained_entry)
            canonical_entry[1] = record["canonical_candidate"]
            actions[output_index] = canonical_entry

            record.update({
                "retained_category_index": slots[output_index][0],
                "retained_locus_key": slots[output_index][1],
                "retained_entry_index": slots[output_index][2],
                "retained_entry_id": id(retained_entry),
                "canonical_entry_id": id(canonical_entry),
                "removed_entry_ids": tuple(
                    id(slots[index][3])
                    for index in member_indexes
                    if index != output_index
                ),
            })
            enriched_provenance.append(record)

        category_updates = {}
        for locus, actions in actions_by_locus.items():
            updated_entries = {}
            expected_candidates = []
            for slot_index, slot in enumerate(slots_by_locus[locus]):
                category_index, locus_key, _entry_index, entry, _candidate = slot
                if slot_index in actions:
                    entry = actions[slot_index]
                    if entry is None:
                        continue
                updated_entries.setdefault(
                    (category_index, locus_key), []
                ).append(entry)
                expected_candidates.append(entry[1])

            collapsed_candidates = collapsed_by_locus[locus]
            if (
                    len(expected_candidates) != len(collapsed_candidates)
                    or any(expected is not collapsed
                           for expected, collapsed in zip(
                               expected_candidates, collapsed_candidates))):
                raise RuntimeError(
                    "source phase-collapse row remapping diverged from output"
                )

            affected_locations = {
                (slot[0], slot[1]) for slot in slots_by_locus[locus]
            }
            for category_index, locus_key in affected_locations:
                category_updates[(category_index, locus_key)] = (
                    updated_entries.get((category_index, locus_key), ())
                )

        for (category_index, locus_key), entries in category_updates.items():
            category = self.categories[category_index]
            if entries:
                current_entries = category[locus_key]
                if isinstance(current_entries, list):
                    current_entries[:] = entries
                else:
                    category[locus_key] = list(entries)
            else:
                del category[locus_key]

        self._record_source_phase_collapse(
            "post_enhancement", enriched_provenance, stats
        )

    def _run_evaluation(self):
        print("3. Evaluate candidates")
        if not any(self.dict_crispr_candidates.values()):
            self._initial_scored_rows = {}
            self.categories = [{}, {}, {}, {}, {}]
        else:
            evaluation = ArrayEvaluation(
                dict_crispr_array_candidates=self.dict_crispr_candidates,
                parameters=self.parameters,
                flag_dev_mode=self.flag_dev_mode,
            )
            self._initial_scored_rows = evaluation.scored_rows()
            self.categories = evaluation.output()
        self._candidate_score_cache = CandidateScoreCache(self.categories)

    @staticmethod
    def _score_rows_with_structural_cache(
            candidates, cached_rows, score_fresh_candidates):
        """Score only geometries absent from an exact first-pass cache."""

        cache = {}
        for rows in cached_rows.values():
            for row in rows:
                if not isinstance(row, (tuple, list)) or len(row) < 3:
                    raise RuntimeError("fast rescue score cache row is malformed")
                fingerprint = candidate_score_fingerprint(row[1])
                cache.setdefault(fingerprint, []).append(row)

        slots = []
        fresh_candidates = {}
        reused = 0
        recomputed = 0
        for interval, interval_candidates in candidates.items():
            for candidate in interval_candidates:
                fingerprint = candidate_score_fingerprint(candidate)
                cached = next((
                    row for row in cache.get(fingerprint, ())
                    if row[1] == candidate
                ), None)
                if cached is None:
                    fresh_candidates.setdefault(interval, []).append(candidate)
                    slots.append((interval, candidate, None))
                    recomputed += 1
                else:
                    reused_row = list(cached)
                    reused_row[1] = candidate
                    slots.append((interval, candidate, reused_row))
                    reused += 1

        fresh_rows = (
            score_fresh_candidates(fresh_candidates)
            if fresh_candidates else {}
        )
        unexpected = set(fresh_rows).difference(fresh_candidates)
        if unexpected:
            raise RuntimeError(
                "fast rescue scoring returned unexpected intervals: {}".format(
                    sorted(unexpected)
                )
            )
        queues = {
            interval: deque(fresh_rows.get(interval, ()))
            for interval in fresh_candidates
        }
        merged = {}
        for interval, candidate, cached in slots:
            if cached is None:
                queue = queues[interval]
                if not queue:
                    raise RuntimeError(
                        "fast rescue scoring omitted a fresh candidate"
                    )
                row = list(queue.popleft())
                if len(row) < 2 or row[1] is not candidate:
                    raise RuntimeError(
                        "fast rescue scoring changed candidate order"
                    )
            else:
                row = cached
            merged.setdefault(interval, []).append(row)
        if any(queues.values()):
            raise RuntimeError("fast rescue scoring returned surplus rows")
        return merged, {"reused": reused, "recomputed": recomputed}

    def _run_fast_routed_start_end_rescue(self):
        """Expand only classifier-routed unresolved fast-run components."""

        enabled = bool(
            self.flags.get("flag_fast_run", False)
            and self.flags.get("flag_fast_run_rescue", False)
        )
        if not enabled:
            return
        print("3b. Rescue routed fast-run components")
        detection = self._raw_detection
        if detection is None:
            raise RuntimeError("fast rescue requires retained raw detection")
        baseline_accepted_signatures = (
            self._accepted_structural_signatures(self.categories)
        )
        baseline_accepted_intervals = (
            self._accepted_intervals_by_identity(self.categories)
        )

        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        routing_records, routing_stats = select_unresolved_components(
            self._initial_scored_rows, structural_filter
        )
        routing_telemetry = [
            {
                "source_interval": list(record["source_interval"]),
                "trigger_score": record["trigger_score"],
                "features": dict(record["features"]),
            }
            for record in routing_records
        ]
        source_intervals = tuple(
            record["source_interval"] for record in routing_records
        )
        selected_clusters, attribution_stats = (
            detection.cluster_indices_for_intervals(source_intervals)
        )
        if not selected_clusters:
            self.fast_run_rescue_stats = dict(
                routing_stats,
                enabled=True,
                selected_clusters=0,
                delta_patterns=0,
                fresh_scores=0,
                routed_components=routing_telemetry,
                accepted_overlap_filtered=0,
                accepted_arrays_preserved=True,
                transaction_aborted=False,
                reused_scores=sum(
                    len(rows) for rows in self._initial_scored_rows.values()
                ),
                **attribution_stats
            )
            print(
                "FAST_RUN_RESCUE_STATS {}".format(json.dumps(
                    self.fast_run_rescue_stats, sort_keys=True
                ))
            )
            return

        merged_fuzzy, detection_stats = detection.run_start_end_delta(
            selected_clusters,
            repeat_seeds_by_cluster=(
                self._fast_rescue_seed_scope(
                    routing_records,
                    selected_clusters,
                    self.flags.get(
                        "fast_run_rescue_seed_scope", "cluster"
                    ),
                )
                if self.flags.get(
                    "fast_run_rescue_seed_scope", "cluster"
                ) != "cluster"
                else None
            ),
        )
        refinement = DetectionRefinement(
            dict_fuzzy_crisprs=merged_fuzzy,
            parameters=self.parameters,
            flag_dev_mode=self.flag_dev_mode,
            source_dna=self.genome_sequence,
            internal_search_context=getattr(self, "internal_search_context", None),
            capture_contextual_spacers=(
                self.parameters.get("param_reconstruction_policy") == "conservative"
            ),
            consume_input=True,
        )
        del merged_fuzzy
        self._record_repeat_boundary_repair("fast_rescue", refinement)
        combined_candidates = refinement.output()
        # Output candidates and pending/audit records have been transferred.
        # The refinement owner otherwise pins all intermediate fuzzy maps
        # through feature scoring and the final rescue transaction.
        del refinement
        combined_candidates, collapse_provenance, collapse_stats = (
            collapse_source_backed_phase_equivalent_candidates(
                self.genome_sequence,
                combined_candidates,
                skip_source_scores_when_boundaries_reject=(
                    getattr(self, "flags", {}).get("output_profile")
                    == "machine"
                ),
            )
        )
        raw_fingerprints = {
            candidate_score_fingerprint(row[1])
            for rows in self._initial_scored_rows.values()
            for row in rows
        }
        refined_before_locality = sum(
            len(candidates) for candidates in combined_candidates.values()
        )
        locality_filtered = {}
        new_candidates_in_routed_components = 0
        accepted_overlap_filtered = 0
        for interval, candidates in combined_candidates.items():
            retained = []
            for candidate in candidates:
                fingerprint = candidate_score_fingerprint(candidate)
                if fingerprint in raw_fingerprints:
                    retained.append(candidate)
                    continue
                source_interval = candidate_source_interval(candidate)
                if self._strictly_overlaps_accepted(
                        candidate, baseline_accepted_intervals):
                    accepted_overlap_filtered += 1
                    continue
                if any(
                        source_interval[0] < routed[1]
                        and routed[0] < source_interval[1]
                        for routed in source_intervals):
                    retained.append(candidate)
                    new_candidates_in_routed_components += 1
            if retained:
                locality_filtered[interval] = retained
        combined_candidates = locality_filtered
        del raw_fingerprints
        merged_rows, score_stats = self._score_rows_with_structural_cache(
            combined_candidates,
            self._initial_scored_rows,
            self._score_candidates_complete,
        )
        categories = ArrayEvaluation.categorize_scored_rows(
            merged_rows,
            parameters=self.parameters,
            flag_dev_mode=self.flag_dev_mode,
        )
        rescued_accepted_signatures = (
            self._accepted_structural_signatures(categories)
        )
        missing_accepted = baseline_accepted_signatures - (
            rescued_accepted_signatures
        )
        if missing_accepted:
            self.fast_run_rescue_stats = dict(
                routing_stats,
                enabled=True,
                selected_clusters=len(selected_clusters),
                fresh_scores=score_stats["recomputed"],
                reused_scores=score_stats["reused"],
                routed_components=routing_telemetry,
                accepted_overlap_filtered=accepted_overlap_filtered,
                accepted_arrays_preserved=False,
                transaction_aborted=True,
                missing_accepted_arrays=sum(missing_accepted.values()),
                **attribution_stats
            )
            self.fast_run_rescue_stats.update(detection_stats)
            print(
                "FAST_RUN_RESCUE_STATS {}".format(json.dumps(
                    self.fast_run_rescue_stats, sort_keys=True
                ))
            )
            return

        self.dict_fuzzy_crisprs = {}
        self.dict_crispr_candidates = combined_candidates
        self._initial_scored_rows = merged_rows
        self.categories = categories
        self._candidate_score_cache = CandidateScoreCache(self.categories)
        self._record_source_phase_collapse(
            "fast_rescue", collapse_provenance, collapse_stats
        )
        self._validate_refined_candidates()
        self.fast_run_rescue_stats = dict(
            routing_stats,
            enabled=True,
            selected_clusters=len(selected_clusters),
            fresh_scores=score_stats["recomputed"],
            reused_scores=score_stats["reused"],
            refined_components=len(combined_candidates),
            routed_components=routing_telemetry,
            refined_candidates_before_locality=refined_before_locality,
            refined_candidates_after_locality=sum(
                len(candidates)
                for candidates in combined_candidates.values()
            ),
            new_candidates_in_routed_components=(
                new_candidates_in_routed_components
            ),
            accepted_overlap_filtered=accepted_overlap_filtered,
            accepted_arrays_preserved=True,
            transaction_aborted=False,
            missing_accepted_arrays=0,
            **attribution_stats
        )
        self.fast_run_rescue_stats.update(detection_stats)
        print(
            "FAST_RUN_RESCUE_STATS {}".format(json.dumps(
                self.fast_run_rescue_stats, sort_keys=True
            ))
        )

    @staticmethod
    def _fast_rescue_edit_distance(first, second):
        """Return exact Levenshtein distance for short repeat sequences."""

        if len(first) > len(second):
            first, second = second, first
        previous = list(range(len(first) + 1))
        for second_index, second_base in enumerate(second, 1):
            current = [second_index]
            for first_index, first_base in enumerate(first, 1):
                current.append(min(
                    current[-1] + 1,
                    previous[first_index] + 1,
                    previous[first_index - 1] + (
                        first_base != second_base
                    ),
                ))
            previous = current
        return previous[-1]

    def _fast_rescue_seed_scope(
            self, routing_records, selected_clusters, scope_name):
        """Select routed lineages or their three-seed local families."""

        if scope_name not in ("lineage", "locus", "family"):
            raise ValueError("fast rescue seed scope is invalid")
        records_by_cluster = {index: [] for index in selected_clusters}
        for record in routing_records:
            attributed, _stats = (
                self._raw_detection.cluster_indices_for_intervals((
                    record["source_interval"],
                ))
            )
            for cluster_index in attributed:
                if cluster_index in records_by_cluster:
                    records_by_cluster[cluster_index].append(record)

        scope = {}
        for cluster_index in selected_clusters:
            records = records_by_cluster[cluster_index]
            trigger_seeds = set()
            reference_repeats = []
            for record in records:
                candidate = record["trigger_candidate"]
                reference_repeats.extend(candidate.list_repeats)
                consensus = getattr(candidate, "consensus", None)
                if isinstance(consensus, str) and consensus:
                    reference_repeats.append(consensus.replace("-", ""))
                provenance = getattr(candidate, "source_path_provenance", ())
                if provenance is None:
                    provenance = ()
                for lineage in provenance:
                    repeat = getattr(lineage, "repeat_candidate", None)
                    if not isinstance(repeat, str) or not repeat:
                        raise RuntimeError(
                            "fast rescue lineage has no valid repeat seed"
                        )
                    trigger_seeds.add(repeat)

            cluster = self._raw_detection.cluster_sequences[cluster_index]
            direct = [
                repeat for repeat in cluster.tuple_repeats
                if repeat in trigger_seeds
            ]
            if scope_name == "lineage":
                scope[cluster_index] = tuple(direct) or cluster.tuple_repeats
                continue

            if scope_name == "locus":
                routed_intervals = tuple(
                    record["source_interval"] for record in records
                )
                local_sequences = {
                    occurrence.sequence
                    for occurrence in self._raw_detection.clusters[
                        cluster_index
                    ].list_cluster_occurrences
                    if any(
                        occurrence.begin < interval[1]
                        and interval[0] < occurrence.end
                        for interval in routed_intervals
                    )
                }
                selected = [
                    repeat for repeat in cluster.tuple_repeats
                    if repeat in trigger_seeds or repeat in local_sequences
                ]
                scope[cluster_index] = (
                    tuple(selected) or cluster.tuple_repeats
                )
                continue

            ranked = sorted(
                enumerate(cluster.tuple_repeats),
                key=lambda item: (
                    min(
                        self._fast_rescue_edit_distance(
                            item[1], reference
                        )
                        for reference in reference_repeats
                    ) if reference_repeats else len(item[1]),
                    item[0],
                ),
            )
            selected = list(direct)
            selected_set = set(selected)
            family_limit = max(3, len(selected))
            for _seed_index, repeat in ranked:
                if repeat in selected_set:
                    continue
                selected.append(repeat)
                selected_set.add(repeat)
                if len(selected) >= family_limit:
                    break
            scope[cluster_index] = tuple(selected) or cluster.tuple_repeats
        return scope

    @staticmethod
    def _accepted_locations(categories):
        locations = {}
        for category_index in (0, 2):
            for locus_key, rows in categories[category_index].items():
                for row in rows:
                    candidate_id = id(row[1])
                    if candidate_id in locations:
                        raise RuntimeError(
                            "accepted candidate identity is not unique"
                        )
                    locations[candidate_id] = (category_index, locus_key)
        return locations

    @staticmethod
    def _accepted_intervals_by_identity(categories):
        intervals = {}
        for category_index in (0, 2):
            for rows in categories[category_index].values():
                for row in rows:
                    candidate_id = id(row[1])
                    if candidate_id in intervals:
                        raise RuntimeError(
                            "accepted candidate identity is not unique"
                        )
                    intervals[candidate_id] = candidate_source_interval(row[1])
        return intervals

    @staticmethod
    def _accepted_structural_signatures(categories):
        signatures = Counter()
        for category_index in (0, 2):
            for rows in categories[category_index].values():
                for row in rows:
                    signatures[(
                        category_index,
                        candidate_score_fingerprint(row[1]),
                    )] += 1
        return signatures

    @staticmethod
    def _strictly_overlaps_accepted(
            candidate, accepted_intervals, allowed_identity=None):
        interval = candidate_source_interval(candidate)
        return any(
            candidate_id != allowed_identity
            and interval[0] < accepted_interval[1]
            and accepted_interval[0] < interval[1]
            for candidate_id, accepted_interval in accepted_intervals.items()
        )

    def _build_fast_phase_variants(
            self, parent_records, structural_filter, baseline_fingerprints):
        records = []
        seen = set(baseline_fingerprints)
        stats = {
            "hypotheses_attempted": 0,
            "invalid_geometry_skipped": 0,
            "duplicate_geometry_skipped": 0,
            "structural_filter_skipped": 0,
        }
        for parent in parent_records:
            candidate = parent["candidate"]
            starts = tuple(candidate.list_repeat_starts)
            ends = tuple(
                start + len(repeat)
                for start, repeat in zip(starts, candidate.list_repeats)
            )
            maximum = 6 if parent["mode"] == "promotion" else 3
            for left in range(maximum + 1):
                for right in range(maximum + 1):
                    if not left and not right:
                        continue
                    stats["hypotheses_attempted"] += 1
                    if parent["mode"] == "promotion":
                        variant_starts = tuple(start - left for start in starts)
                        variant_ends = tuple(end + right for end in ends)
                    else:
                        variant_starts = tuple(start + left for start in starts)
                        variant_ends = tuple(end - right for end in ends)
                    valid_geometry = (
                        all(
                            0 <= start < end <= len(self.genome_sequence)
                            for start, end in zip(
                                variant_starts, variant_ends
                            )
                        )
                        and all(
                            end <= next_start
                            for end, next_start in zip(
                                variant_ends, variant_starts[1:]
                            )
                        )
                    )
                    if not valid_geometry:
                        stats["invalid_geometry_skipped"] += 1
                        continue
                    variant = build_source_backed_candidate_from_boundaries(
                        self.genome_sequence,
                        None,
                        variant_starts,
                        variant_ends,
                    )
                    fingerprint = candidate_score_fingerprint(variant)
                    if fingerprint in seen:
                        stats["duplicate_geometry_skipped"] += 1
                        continue
                    seen.add(fingerprint)
                    if not structural_filter(variant):
                        stats["structural_filter_skipped"] += 1
                        continue
                    records.append({
                        **parent,
                        "variant": variant,
                        "left": left,
                        "right": right,
                    })
        return records, stats

    def _run_fast_bounded_phase_repair(self):
        """Classifier-arbitrate a hard-capped source-backed phase grid."""

        enabled = bool(
            self.flags.get("flag_fast_run", False)
            and self.flags.get("flag_fast_run_rescue", False)
            and self.flags.get("flag_fast_run_phase_repair", False)
        )
        if not enabled:
            return
        if not self.parameters.get(
                "param_split_competitive_components", False):
            raise RuntimeError(
                "fast phase repair requires independent source components"
            )
        print("4c. Repair bounded fast-run repeat phase")
        baseline_rows = self._post_enhancement_scored_rows
        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )

        promotion_parents = []
        for locus_key in sorted(self.categories[2]):
            for row in self.categories[2][locus_key]:
                features = compact_feature_values(row)
                candidate = row[1]
                eligible = getattr(
                    candidate, "source_path_classification_eligible", True
                )
                if not isinstance(eligible, bool):
                    raise RuntimeError(
                        "fast phase parent eligibility is not boolean"
                    )
                if (
                        eligible
                        and 0.5 <= float(row[0]) < 0.60
                        and features["repeat_similarity"] >= 0.80
                        and features["number_mismatches"] <= 3.0
                        # A uniform phase truncation transfers repeat-edge
                        # bases into every spacer. Restrict the expensive
                        # extension grid to that observable short-repeat /
                        # long-spacer signature.
                        and features["repeat_len"] <= 30.0
                        and 50.0 <= features["avg_spacer_len"] <= 60.0
                        and structural_filter(candidate)):
                    source_interval = candidate_source_interval(candidate)
                    promotion_parents.append({
                        "mode": "promotion",
                        "locus_key": locus_key,
                        "base_row": row,
                        "candidate": candidate,
                        "base_score": float(row[0]),
                        "source_interval": source_interval,
                    })
        promotion_parents.sort(key=lambda record: (
            record["base_score"],
            record["source_interval"][0],
            record["source_interval"][1],
            repr(record["locus_key"]),
        ))
        promotion_parents = promotion_parents[:8]

        routed, _routing_stats = select_unresolved_components(
            baseline_rows, structural_filter
        )
        row_by_candidate_id = {
            id(row[1]): row
            for rows in baseline_rows.values()
            for row in rows
        }
        admission_parents = []
        for record in routed[:8]:
            candidate = record["trigger_candidate"]
            base_row = row_by_candidate_id.get(id(candidate))
            if base_row is None or base_row[1] is not candidate:
                raise RuntimeError(
                    "routed phase parent lost score-row identity"
                )
            admission_parents.append({
                "mode": "admission",
                "locus_key": record["locus_key"],
                "base_row": base_row,
                "candidate": candidate,
                "base_score": float(base_row[0]),
                "source_interval": record["source_interval"],
            })

        parents = promotion_parents + admission_parents
        baseline_fingerprints = {
            candidate_score_fingerprint(row[1])
            for rows in baseline_rows.values()
            for row in rows
        }
        variant_records, construction_stats = (
            self._build_fast_phase_variants(
                parents, structural_filter, baseline_fingerprints
            )
        )
        if variant_records:
            score_key = "fast-bounded-phase-repair"
            scored = self._score_candidates_complete({
                score_key: [
                    record["variant"] for record in variant_records
                ]
            })
            if set(scored) != {score_key}:
                raise RuntimeError(
                    "fast phase scoring returned an unexpected batch"
                )
            scored_rows = list(scored[score_key])
            if len(scored_rows) != len(variant_records):
                raise RuntimeError(
                    "fast phase scoring changed candidate cardinality"
                )
            for record, row in zip(variant_records, scored_rows):
                if len(row) < 2 or row[1] is not record["variant"]:
                    raise RuntimeError(
                        "fast phase scoring changed candidate identity"
                    )
                record["scored_row"] = row

        proposals_by_parent = {}
        for record in variant_records:
            score = float(record["scored_row"][0])
            threshold = 0.75 if record["mode"] == "promotion" else 0.5
            if score < threshold:
                continue
            parent_id = id(record["candidate"])
            current = proposals_by_parent.get(parent_id)
            rank = (
                score,
                -(record["left"] + record["right"]),
                -record["left"],
                -record["right"],
            )
            if current is None or rank > current[0]:
                proposals_by_parent[parent_id] = (rank, record)

        committed_rows = self._copy_scored_rows(baseline_rows)
        committed_categories = self.categories
        committed_records = []
        rejected_proposals = 0
        overlap_rejected_proposals = 0
        proposals = [
            record for _rank, record in proposals_by_parent.values()
        ]
        proposals.sort(key=lambda record: (
            0 if record["mode"] == "promotion" else 1,
            record["source_interval"][0],
            record["source_interval"][1],
            repr(record["locus_key"]),
        ))
        for record in proposals:
            prior_locations = self._accepted_locations(committed_categories)
            prior_intervals = self._accepted_intervals_by_identity(
                committed_categories
            )
            prior_count = len(prior_locations)
            allowed_overlap_id = (
                id(record["candidate"])
                if record["mode"] == "promotion" else None
            )
            if self._strictly_overlaps_accepted(
                    record["variant"], prior_intervals, allowed_overlap_id):
                rejected_proposals += 1
                overlap_rejected_proposals += 1
                continue
            trial_rows = self._copy_scored_rows(committed_rows)
            trial_rows.setdefault(record["locus_key"], []).append(
                record["scored_row"]
            )
            trial_categories = ArrayEvaluation.categorize_scored_rows(
                trial_rows,
                parameters=self.parameters,
                flag_dev_mode=self.flag_dev_mode,
            )
            trial_locations = self._accepted_locations(trial_categories)
            variant_id = id(record["variant"])
            base_id = id(record["candidate"])
            preserved_ids = set(prior_locations)
            expected_count = prior_count + 1
            valid = True
            if record["mode"] == "promotion":
                preserved_ids.discard(base_id)
                expected_count = prior_count
                valid = (
                    trial_locations.get(variant_id, (None,))[0] == 0
                    and base_id not in trial_locations
                )
            else:
                valid = (
                    trial_locations.get(variant_id, (None,))[0] in (0, 2)
                )
            valid = (
                valid
                and len(trial_locations) == expected_count
                and all(
                    trial_locations.get(candidate_id) == location
                    for candidate_id, location in prior_locations.items()
                    if candidate_id in preserved_ids
                )
            )
            if not valid:
                rejected_proposals += 1
                continue
            committed_rows = trial_rows
            committed_categories = trial_categories
            committed_records.append({
                "mode": record["mode"],
                "source_interval_before": list(record["source_interval"]),
                "source_interval_after": list(candidate_source_interval(
                    record["variant"]
                )),
                "base_score": record["base_score"],
                "variant_score": float(record["scored_row"][0]),
                "left": record["left"],
                "right": record["right"],
            })

        self.categories = committed_categories
        self._post_enhancement_scored_rows = committed_rows
        self._final_scored_rows = committed_rows
        self._candidate_score_cache = CandidateScoreCache(self.categories)
        self.fast_run_phase_repair_stats = {
            "enabled": True,
            "promotion_parents": len(promotion_parents),
            "admission_parents": len(admission_parents),
            **construction_stats,
            "variants_scored": len(variant_records),
            "proposals": len(proposals),
            "proposals_rejected": rejected_proposals,
            "overlap_rejected_proposals": overlap_rejected_proposals,
            "promotions_committed": sum(
                record["mode"] == "promotion"
                for record in committed_records
            ),
            "admissions_committed": sum(
                record["mode"] == "admission"
                for record in committed_records
            ),
            "committed": committed_records,
        }
        print(
            "FAST_RUN_PHASE_REPAIR_STATS {}".format(json.dumps(
                self.fast_run_phase_repair_stats, sort_keys=True
            ))
        )

    def _run_fast_two_repeat_terminal_repair(self):
        """Add only ordinary-model-accepted source winners transactionally."""

        enabled = bool(
            self.flags.get("flag_fast_run", False)
            and self.flags.get(
                "flag_fast_two_repeat_terminal_repair", False
            )
        )
        if not enabled:
            return
        if not self.parameters.get("param_accepted_only", False):
            raise RuntimeError(
                "fast two-repeat terminal repair requires accepted-only mode"
            )
        if self._raw_detection is None:
            raise RuntimeError(
                "fast two-repeat terminal repair lost raw detection"
            )

        print("4c. Repair narrow two-repeat terminal candidates")
        seeds, seed_stats = (
            self._raw_detection.two_repeat_terminal_repair_seeds()
        )
        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        seeds, routing_stats = self._route_fast_two_repeat_terminal_seeds(
            seeds
        )
        winner_records, generation_stats = generate_source_winners(
            self.genome_sequence, seeds, structural_filter
        )

        candidate_batches = {}
        for record in winner_records:
            candidate = record["candidate"]
            candidate_batches.setdefault(
                candidate_source_interval(candidate), []
            ).append(candidate)
        scored = self._score_candidates_complete(candidate_batches)
        unexpected = set(scored).difference(candidate_batches)
        if unexpected:
            raise RuntimeError(
                "two-repeat repair scoring returned unexpected intervals"
            )
        rows_by_candidate_id = {}
        for interval, candidates in candidate_batches.items():
            rows = scored.get(interval, ())
            if len(rows) != len(candidates):
                raise RuntimeError(
                    "two-repeat repair scoring changed cardinality"
                )
            for candidate, row in zip(candidates, rows):
                if (
                        not isinstance(row, (tuple, list))
                        or len(row) < 2
                        or row[1] is not candidate):
                    raise RuntimeError(
                        "two-repeat repair scoring changed identity or order"
                    )
                rows_by_candidate_id[id(candidate)] = row

        ordinary_model_passes = 0
        stage_model_score_rejections = 0
        repeat_derived_spacer_rejections = 0
        flank_continuation_rejections = 0
        proposals = []
        for record in winner_records:
            row = rows_by_candidate_id[id(record["candidate"])]
            score = float(row[0])
            if score < 0.5:
                continue
            ordinary_model_passes += 1
            if score < MINIMUM_ADMISSION_MODEL_SCORE:
                stage_model_score_rejections += 1
                continue
            repeat_derived_fraction = record[
                "maximum_repeat_derived_spacer_coverage"
            ]
            if repeat_derived_fraction >= (
                    MAXIMUM_ADMISSION_REPEAT_DERIVED_SPACER_COVERAGE):
                repeat_derived_spacer_rejections += 1
                continue
            flank_evidence = flank_continuation_identities(
                self.genome_sequence, record["candidate"]
            )
            record.update(flank_evidence)
            if not fast_two_repeat_admission_allowed(
                    score,
                    repeat_derived_fraction,
                    flank_evidence["maximum_repeat_flank_identity"],
                    flank_evidence["maximum_spacer_flank_identity"]):
                flank_continuation_rejections += 1
                continue
            record["scored_row"] = row
            proposals.append(record)

        committed_rows = self._copy_scored_rows(
            self._post_enhancement_scored_rows
        )
        committed_categories = self.categories
        committed = []
        overlap_rejections = 0
        transaction_rejections = 0
        for record in proposals:
            candidate = record["candidate"]
            prior_locations = self._accepted_locations(
                committed_categories
            )
            prior_intervals = self._accepted_intervals_by_identity(
                committed_categories
            )
            if self._strictly_overlaps_accepted(
                    candidate, prior_intervals):
                overlap_rejections += 1
                continue
            trial_rows = self._copy_scored_rows(committed_rows)
            locus_key = candidate_source_interval(candidate)
            trial_rows.setdefault(locus_key, []).append(
                record["scored_row"]
            )
            trial_categories = ArrayEvaluation.categorize_scored_rows(
                trial_rows,
                parameters=self.parameters,
                flag_dev_mode=self.flag_dev_mode,
            )
            trial_locations = self._accepted_locations(trial_categories)
            candidate_id = id(candidate)
            valid = (
                trial_locations.get(candidate_id, (None,))[0] in (0, 2)
                and len(trial_locations) == len(prior_locations) + 1
                and all(
                    trial_locations.get(identity) == location
                    for identity, location in prior_locations.items()
                )
            )
            if not valid:
                transaction_rejections += 1
                continue
            committed_rows = trial_rows
            committed_categories = trial_categories
            committed.append({
                "source_interval": list(locus_key),
                "seed_interval": list(record["seed_interval"]),
                "score": float(record["scored_row"][0]),
                "source_phase_score": record["source_phase_score"],
                "maximum_repeat_derived_spacer_coverage": record[
                    "maximum_repeat_derived_spacer_coverage"
                ],
                "repeat_dust_triplet_score": record[
                    "repeat_dust_triplet_score"
                ],
                "maximum_repeat_flank_identity": record[
                    "maximum_repeat_flank_identity"
                ],
                "maximum_spacer_flank_identity": record[
                    "maximum_spacer_flank_identity"
                ],
            })

        self._post_enhancement_scored_rows = committed_rows
        self._final_scored_rows = committed_rows
        self.categories = committed_categories
        self._candidate_score_cache = CandidateScoreCache(self.categories)
        self.fast_two_repeat_terminal_repair_stats = {
            "enabled": True,
            **seed_stats,
            **routing_stats,
            **generation_stats,
            "seeds_retained": len(seeds),
            "source_winners": len(winner_records),
            "fresh_scores": sum(map(len, scored.values())),
            "ordinary_score_passes": ordinary_model_passes,
            "stage_admission_policy": (
                "score>=0.57_repeat_derived_spacer_coverage<0.75_"
                "reject_repeat_flank>=0.90_and_spacer_flank>=0.55_"
                "within_separate_250bp_source_flanks"
            ),
            "stage_model_score_rejections": (
                stage_model_score_rejections
            ),
            "repeat_derived_spacer_rejections": (
                repeat_derived_spacer_rejections
            ),
            "flank_continuation_rejections": (
                flank_continuation_rejections
            ),
            "stage_admission_policy_passes": len(proposals),
            "overlap_rejections": overlap_rejections,
            "transaction_rejections": transaction_rejections,
            "admissions_committed": len(committed),
            "committed": committed,
        }
        print(
            "FAST_TWO_REPEAT_TERMINAL_REPAIR_STATS {}".format(
                json.dumps(
                    self.fast_two_repeat_terminal_repair_stats,
                    sort_keys=True,
                )
            )
        )

    def _run_fast_accepted_family_rescue(self):
        """Admit only source-replayed, accepted-family-supported satellites."""

        enabled = bool(
            self.flags.get("flag_fast_run", False)
            and self.flags.get("flag_fast_accepted_family_rescue", False)
        )
        if not enabled:
            return
        if not self.parameters.get("param_accepted_only", False):
            raise RuntimeError(
                "fast accepted-family rescue requires accepted-only mode"
            )
        if self._raw_detection is None:
            raise RuntimeError(
                "fast accepted-family rescue lost raw detection"
            )
        if self._fast_accepted_family_anchor_categories is None:
            raise RuntimeError(
                "fast accepted-family rescue lost its pre-rescue anchors"
            )

        print("4d. Rescue source-backed accepted-family satellites")
        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        raw_seeds, raw_capture_stats = (
            self._raw_detection.two_repeat_terminal_repair_seeds()
        )
        raw_parents, raw_selection_stats = (
            select_exact_template_raw_seed_parents(
                self.genome_sequence,
                self.categories,
                raw_seeds,
                structural_filter,
                accepted_categories=(
                    self._fast_accepted_family_anchor_categories
                ),
            )
        )
        raw_parent_batches = {}
        for parent_record in raw_parents:
            parent = parent_record.candidate
            raw_parent_batches.setdefault(
                candidate_source_interval(parent), []
            ).append(parent)
        raw_parent_scored = self._score_candidates_complete(
            raw_parent_batches
        )
        raw_parent_rows_by_candidate_id = {}
        for interval, parents in raw_parent_batches.items():
            rows = raw_parent_scored.get(interval, ())
            if len(rows) != len(parents):
                raise RuntimeError(
                    "accepted-family raw-parent scoring changed cardinality"
                )
            for parent, row in zip(parents, rows):
                if (
                        not isinstance(row, (tuple, list))
                        or len(row) < 2
                        or row[1] is not parent):
                    raise RuntimeError(
                        "accepted-family raw-parent scoring changed identity/order"
                    )
                raw_parent_rows_by_candidate_id[id(parent)] = row
        raw_seed_parent_rows = tuple(
            (
                raw_parent_rows_by_candidate_id[id(record.candidate)],
                record,
            )
            for record in raw_parents
        )
        proposals, generation_stats = (
            generate_accepted_family_rescue_proposals(
                self.genome_sequence,
                self.categories,
                structural_filter,
                raw_seed_parent_rows=raw_seed_parent_rows,
                accepted_categories=(
                    self._fast_accepted_family_anchor_categories
                ),
            )
        )

        batches = {}
        for proposal in proposals:
            candidate = proposal.candidate
            batches.setdefault(
                candidate_source_interval(candidate), []
            ).append(candidate)
        scored = self._score_candidates_complete(batches)
        rows_by_candidate_id = {}
        for interval, candidates in batches.items():
            rows = scored.get(interval, ())
            if len(rows) != len(candidates):
                raise RuntimeError(
                    "accepted-family rescue scoring changed cardinality"
                )
            for candidate, row in zip(candidates, rows):
                if (
                        not isinstance(row, (tuple, list))
                        or len(row) < 3
                        or row[1] is not candidate):
                    raise RuntimeError(
                        "accepted-family rescue scoring changed identity/order"
                    )
                rows_by_candidate_id[id(candidate)] = row

        committed_categories = self._copy_categories(self.categories)
        committed_rows = self._copy_scored_rows(
            self._post_enhancement_scored_rows
        )
        committed = []
        score_rejections = 0
        overlap_rejections = 0
        transaction_rejections = 0
        for proposal in proposals:
            candidate = proposal.candidate
            raw_row = rows_by_candidate_id[id(candidate)]
            score = float(raw_row[0])
            if not accepted_family_score_is_admissible(score):
                score_rejections += 1
                continue
            prior_locations = self._accepted_locations(
                committed_categories
            )
            prior_intervals = self._accepted_intervals_by_identity(
                committed_categories
            )
            if self._strictly_overlaps_accepted(candidate, prior_intervals):
                overlap_rejections += 1
                continue

            metadata = {
                "classification_basis": (
                    accepted_family_classification_basis(
                        self.genome_sequence, proposal, score
                    )
                )
            }
            admitted_row = list(raw_row[:3]) + [metadata]
            trial_categories = self._copy_categories(committed_categories)
            locus_key = candidate_source_interval(candidate)
            trial_categories[2].setdefault(locus_key, []).append(
                admitted_row
            )
            trial_locations = self._accepted_locations(trial_categories)
            candidate_id = id(candidate)
            valid = (
                trial_locations.get(candidate_id) == (2, locus_key)
                and len(trial_locations) == len(prior_locations) + 1
                and all(
                    trial_locations.get(identity) == location
                    for identity, location in prior_locations.items()
                )
            )
            if not valid:
                transaction_rejections += 1
                continue

            trial_rows = self._copy_scored_rows(committed_rows)
            trial_rows.setdefault(locus_key, []).append(admitted_row)
            committed_categories = trial_categories
            committed_rows = trial_rows
            committed.append({
                "arm": proposal.arm,
                "source_interval": list(locus_key),
                "repeat_intervals": [
                    list(interval) for interval in proposal.repeat_intervals
                ],
                "score": score,
                "parent_score": proposal.parent_score,
                "family_anchor_count": len(proposal.anchors),
                "family_template": proposal.family_template,
            })

        self.categories = committed_categories
        self._post_enhancement_scored_rows = committed_rows
        self._final_scored_rows = committed_rows
        self._candidate_score_cache = CandidateScoreCache(self.categories)
        self.fast_accepted_family_rescue_stats = {
            "enabled": True,
            **{
                "raw_seed_capture_{}".format(key): value
                for key, value in raw_capture_stats.items()
            },
            **raw_selection_stats,
            **generation_stats,
            "raw_seed_parent_fresh_scores": sum(
                map(len, raw_parent_scored.values())
            ),
            "fresh_scores": sum(map(len, scored.values())),
            "total_fresh_scores": (
                sum(map(len, raw_parent_scored.values()))
                + sum(map(len, scored.values()))
            ),
            "score_rejections": score_rejections,
            "overlap_rejections": overlap_rejections,
            "transaction_rejections": transaction_rejections,
            "admissions_committed": len(committed),
            "committed": committed,
        }
        print(
            "FAST_ACCEPTED_FAMILY_RESCUE_STATS {}".format(
                json.dumps(
                    self.fast_accepted_family_rescue_stats,
                    sort_keys=True,
                )
            )
        )

    def _route_fast_two_repeat_terminal_seeds(self, seeds):
        """Score cheap parents and cap the expensive source lattice at eight."""

        accepted_intervals = self._accepted_intervals_by_identity(
            self.categories
        )
        records = []
        candidates = {}
        candidate_to_seed = {}
        overlapping_accepted = 0
        represented_seed_clusters = 0
        for seed in seeds:
            cluster_index = seed.cluster_index
            if (
                    not isinstance(cluster_index, int)
                    or isinstance(cluster_index, bool)
                    or not 0 <= cluster_index < len(
                        self._raw_detection.cluster_sequences
                    )):
                raise RuntimeError(
                    "two-repeat repair seed lost its source cluster"
                )
            cluster = self._raw_detection.cluster_sequences[cluster_index]
            cluster_interval = (cluster.start, cluster.end)
            if any(
                    cluster_interval[0] < accepted_interval[1]
                    and accepted_interval[0] < cluster_interval[1]
                    for accepted_interval in accepted_intervals.values()):
                represented_seed_clusters += 1
                continue
            parent = build_source_backed_candidate_from_boundaries(
                self.genome_sequence,
                None,
                seed.repeat_starts,
                seed.repeat_ends,
            )
            if self._strictly_overlaps_accepted(
                    parent, accepted_intervals):
                overlapping_accepted += 1
                continue
            interval = candidate_source_interval(parent)
            candidates.setdefault(interval, []).append(parent)
            candidate_to_seed[id(parent)] = seed

        scored = self._score_candidates_complete(candidates)
        for interval, parents in candidates.items():
            rows = scored.get(interval, ())
            if len(rows) != len(parents):
                raise RuntimeError(
                    "two-repeat parent routing changed cardinality"
                )
            for parent, row in zip(parents, rows):
                if (
                        not isinstance(row, (tuple, list))
                        or len(row) < 2
                        or row[1] is not parent):
                    raise RuntimeError(
                        "two-repeat parent routing changed identity or order"
                    )
                features = compact_feature_values(row)
                score = float(row[0])
                if not (
                        score >= 0.35
                        and features["repeat_similarity"] >= 0.85
                        and features["number_mismatches"] <= 1.0
                        and features["blast_score_1"] >= 15.0
                        and 21.0 <= features["repeat_len"] <= 40.0
                        and 24.0 <= features["avg_spacer_len"] <= 50.0):
                    continue
                records.append({
                    "seed": candidate_to_seed[id(parent)],
                    "interval": interval,
                    "score": score,
                    "features": features,
                })

        records.sort(key=lambda record: (
            -record["score"],
            -record["features"]["repeat_similarity"],
            -record["features"]["blast_score_1"],
            record["features"]["number_mismatches"],
            record["interval"],
        ))
        eligible_before_cap = len(records)
        records = records[
            :MAX_FAST_TWO_REPEAT_ROUTED_PARENTS_PER_SOURCE
        ]
        return tuple(record["seed"] for record in records), {
            "parent_routing_policy": (
                "score>=0.35_similarity>=0.85_mismatches<=1_"
                "blast>=15_top8_per_source"
            ),
            "parents_scored_for_routing": sum(map(len, scored.values())),
            "parents_overlapping_accepted_before_routing": (
                overlapping_accepted
            ),
            "parents_in_clusters_with_accepted_calls": (
                represented_seed_clusters
            ),
            "parents_eligible_before_cap": eligible_before_cap,
            "parents_routed": len(records),
            "parent_route_cap": (
                MAX_FAST_TWO_REPEAT_ROUTED_PARENTS_PER_SOURCE
            ),
            "parent_route_cap_skipped": (
                eligible_before_cap - len(records)
            ),
            "routed_parents": [
                {
                    "source_interval": list(record["interval"]),
                    "score": record["score"],
                    "repeat_similarity": record["features"][
                        "repeat_similarity"
                    ],
                    "number_mismatches": record["features"][
                        "number_mismatches"
                    ],
                    "blast_score_1": record["features"]["blast_score_1"],
                }
                for record in records
            ],
        }

    def _evaluate_candidates(self, candidates):
        if not any(candidates.values()):
            return [{}, {}, {}, {}, {}]
        ae = ArrayEvaluation(dict_crispr_array_candidates=candidates,
                             parameters=self.parameters,
                             flag_dev_mode=self.flag_dev_mode)
        return ae.output()

    def _score_candidates_complete(self, candidates):
        """Score a candidate subset with the compact model."""

        if not any(candidates.values()):
            return {}
        evaluation = ArrayEvaluation(
            dict_crispr_array_candidates=candidates,
            parameters=self.parameters,
            flag_dev_mode=self.flag_dev_mode,
        )
        return evaluation.scored_rows()

    def _candidate_mapping_from_categories(self):
        """Recover every candidate so enhanced arrays can be scored afresh."""

        candidates = {}
        for category in self.categories:
            for interval in sorted(category):
                interval_candidates = candidates.setdefault(interval, [])
                for entry in category[interval]:
                    if not isinstance(entry, (tuple, list)) or len(entry) < 2:
                        raise RuntimeError(
                            "cannot re-evaluate a malformed candidate entry"
                        )
                    interval_candidates.append(entry[1])
        return candidates

    @staticmethod
    def _scored_rows_from_categories(categories):
        """Flatten categorized rows for compatibility-stage provenance."""

        scored_rows = {}
        for category in categories:
            for interval in sorted(category):
                scored_rows.setdefault(interval, []).extend(
                    category[interval]
                )
        return scored_rows

    def _reevaluate_enhanced_candidates(
            self, project_requested_categories=True):
        """Refresh features, scores, ranking, and categories after enhancement."""

        print("4b. Re-evaluate enhanced arrays")
        score_cache = getattr(self, "_candidate_score_cache", None)
        if score_cache is None:
            # Compatibility path for callers that enter the pipeline between
            # stages without the first-pass provenance snapshot.
            candidates = self._candidate_mapping_from_categories()
            self.categories = self._evaluate_candidates(candidates)
            merged_rows = self._scored_rows_from_categories(self.categories)
            self.reevaluation_stats = {
                "reused": 0,
                "recomputed": sum(len(rows) for rows in candidates.values()),
            }
            self._post_enhancement_scored_rows = merged_rows
            self._final_scored_rows = merged_rows
            if project_requested_categories:
                self._drop_nonaccepted_categories()
            return
        else:
            merged_rows = score_cache.merge_after_enhancement(
                self.categories, self._score_candidates_complete
            )
            self.reevaluation_stats = {
                "reused": score_cache.reused_count,
                "recomputed": score_cache.recomputed_count,
            }
        self._post_enhancement_scored_rows = merged_rows
        self._final_scored_rows = merged_rows
        self.categories = ArrayEvaluation.categorize_scored_rows(
            merged_rows,
            parameters=self.parameters,
            flag_dev_mode=self.flag_dev_mode,
        )
        if project_requested_categories:
            self._drop_nonaccepted_categories()

    @staticmethod
    def _accepted_entry_count(categories):
        return sum(
            len(entries)
            for category_index in (0, 2)
            for entries in categories[category_index].values()
        )

    @staticmethod
    def _accepted_candidate_ids(categories):
        return {
            id(entry[1])
            for category_index in (0, 2)
            for entries in categories[category_index].values()
            for entry in entries
        }

    @staticmethod
    def _copy_scored_rows(scored_rows):
        return {
            interval: list(rows) for interval, rows in scored_rows.items()
        }

    @staticmethod
    def _copy_categories(categories):
        return [
            {
                interval: list(rows)
                for interval, rows in category.items()
            }
            for category in categories
        ]

    @staticmethod
    def _rows_by_candidate_id(scored_rows):
        indexed = {}
        for rows in scored_rows.values():
            for row in rows:
                if not isinstance(row, (list, tuple)) or len(row) < 3:
                    raise RuntimeError(
                        "architecture rescue received a malformed score row"
                    )
                candidate_id = id(row[1])
                previous = indexed.get(candidate_id)
                if previous is not None and previous is not row:
                    raise RuntimeError(
                        "one candidate object has multiple scored rows"
                    )
                indexed[candidate_id] = row
        return indexed

    @staticmethod
    def _remove_category_rows_by_identity(categories, row_ids):
        for category in categories:
            for interval in tuple(category):
                retained = [
                    row for row in category[interval]
                    if id(row) not in row_ids
                ]
                if retained:
                    category[interval] = retained
                else:
                    del category[interval]

    @staticmethod
    def _migrate_category_locus_keys(
            categories, target_locus_key, migrated_locus_keys):
        """Coalesce all representations of one rescued region under one key."""

        for category in categories:
            target_rows = list(category.get(target_locus_key, ()))
            target_row_ids = {id(row) for row in target_rows}
            for locus_key in migrated_locus_keys:
                if locus_key == target_locus_key:
                    continue
                for row in category.pop(locus_key, ()):
                    if id(row) not in target_row_ids:
                        target_rows.append(row)
                        target_row_ids.add(id(row))
            if target_rows:
                category[target_locus_key] = target_rows

    def _source_architecture_rescue_enabled(self):
        """Return the effective source-augmentation switch.

        Architecture rescue can add terminal repeats, so the older
        ``--degenerated False`` contract is an unconditional outer gate.
        """

        return (
            self.flags.get("flag_degenerated", True)
            and self.parameters.get(
                "param_source_backed_architecture_rescue", False
            )
        )

    def _run_configured_source_backed_augmentation(self):
        """Run the selected rescue or preserve post-enhancement rows exactly."""

        if not self._source_architecture_rescue_enabled():
            self._final_scored_rows = self._post_enhancement_scored_rows
            return False
        self._run_source_backed_architecture_rescue_transactionally()
        return True

    @staticmethod
    def _validated_coding_evidence_by_candidate_id(
            source_sequence, candidates, coding_batch, decisions):
        """Validate source-wide context cardinality, order, and coordinates."""

        source_length = getattr(coding_batch, "source_length", None)
        if (
            not isinstance(source_length, int)
            or isinstance(source_length, bool)
            or source_length != len(source_sequence)
        ):
            raise RuntimeError(
                "coding-context batch has the wrong source length"
            )
        try:
            evidence_records = tuple(coding_batch.candidates)
        except (AttributeError, TypeError) as error:
            raise RuntimeError(
                "coding-context batch has no iterable candidate evidence"
            ) from error
        if len(evidence_records) != len(candidates):
            raise RuntimeError(
                "coding-context evidence changed candidate cardinality"
            )
        if len({id(candidate) for candidate in candidates}) != len(candidates):
            raise RuntimeError(
                "coding-context candidates must have unique object identities"
            )

        evidence_by_candidate_id = {}
        for candidate_index, (candidate, evidence) in enumerate(zip(
                candidates, evidence_records)):
            observed_index = getattr(evidence, "candidate_index", None)
            if (
                not isinstance(observed_index, int)
                or isinstance(observed_index, bool)
                or observed_index != candidate_index
            ):
                raise RuntimeError(
                    "coding-context evidence changed candidate order/index"
                )
            decision = decisions[id(candidate)]
            expected_interval = (
                decision.stats.array_start,
                decision.stats.array_end,
            )
            observed_interval = (
                getattr(evidence, "candidate_start", None),
                getattr(evidence, "candidate_end", None),
            )
            if (
                any(
                    not isinstance(value, int) or isinstance(value, bool)
                    for value in observed_interval
                )
                or observed_interval != expected_interval
            ):
                raise RuntimeError(
                    "coding-context evidence changed candidate interval"
                )

            overlap_fraction = getattr(
                evidence, "coding_overlap_fraction", None
            )
            context_score = getattr(evidence, "context_cds_score", None)
            if (
                isinstance(overlap_fraction, bool)
                or not isinstance(overlap_fraction, (int, float))
                or not math.isfinite(float(overlap_fraction))
                or not 0.0 <= float(overlap_fraction) <= 1.0
            ):
                raise RuntimeError(
                    "coding-context overlap fraction is invalid"
                )
            if (
                isinstance(context_score, bool)
                or not isinstance(context_score, (int, float))
                or not math.isfinite(float(context_score))
                or not 0.0 <= float(context_score) <= 100.0
            ):
                raise RuntimeError("coding-context CDS score is invalid")

            selected_cds = getattr(evidence, "selected_cds", None)
            if selected_cds is None:
                if (
                    float(overlap_fraction) != 0.0
                    or float(context_score) != 0.0
                ):
                    raise RuntimeError(
                        "coding-context zero evidence has no selected CDS"
                    )
            else:
                cds_start = getattr(selected_cds, "start", None)
                cds_end = getattr(selected_cds, "end", None)
                is_complement = getattr(
                    selected_cds, "is_complement", None
                )
                confidence = getattr(selected_cds, "confidence", None)
                strand = getattr(selected_cds, "strand", None)
                expected_strand = "-" if is_complement is True else "+"
                if (
                    not isinstance(cds_start, int)
                    or isinstance(cds_start, bool)
                    or not isinstance(cds_end, int)
                    or isinstance(cds_end, bool)
                    or cds_start < 0
                    or cds_end <= cds_start
                    or cds_end > source_length
                    or not isinstance(is_complement, bool)
                    or strand != expected_strand
                    or isinstance(confidence, bool)
                    or not isinstance(confidence, (int, float))
                    or not math.isfinite(float(confidence))
                    or not 0.0 <= float(confidence) <= 100.0
                ):
                    raise RuntimeError(
                        "coding-context selected CDS is invalid"
                    )
                expected_overlap = max(
                    0,
                    min(observed_interval[1], cds_end)
                    - max(observed_interval[0], cds_start),
                ) / float(observed_interval[1] - observed_interval[0])
                if expected_overlap <= 0.0:
                    raise RuntimeError(
                        "coding-context selected CDS does not overlap candidate"
                    )
                if not math.isclose(
                    float(overlap_fraction),
                    expected_overlap,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise RuntimeError(
                        "coding-context overlap does not match selected CDS"
                    )
                if not math.isclose(
                    float(context_score),
                    float(confidence),
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ):
                    raise RuntimeError(
                        "coding-context score does not match selected CDS"
                    )
            evidence_by_candidate_id[id(candidate)] = evidence
        return evidence_by_candidate_id

    @staticmethod
    def _coding_evidence_classification_basis(evidence):
        selected_cds = evidence.selected_cds
        return {
            "coordinate_system": "zero_based_half_open",
            "candidate_interval": [
                evidence.candidate_start,
                evidence.candidate_end,
            ],
            "coding_overlap_fraction": float(
                evidence.coding_overlap_fraction
            ),
            "context_cds_score": float(evidence.context_cds_score),
            "selected_cds": (
                None
                if selected_cds is None
                else {
                    "start": selected_cds.start,
                    "end": selected_cds.end,
                    "strand": selected_cds.strand,
                    "confidence": float(selected_cds.confidence),
                }
            ),
        }

    def _rescued_row_with_basis(
            self, transition, source_normalizer=None):
        row = transition["row"]
        if len(row) > 4:
            raise RuntimeError("scored row has unsupported trailing fields")
        metadata = {}
        if len(row) == 4:
            if not isinstance(row[3], dict):
                raise RuntimeError("scored row metadata must be a mapping")
            metadata.update(row[3])
        score = float(row[0])
        if source_normalizer is None:
            source_normalizer = BoundaryAssessmentNormalizer.from_source(
                self.genome_sequence
            )
        else:
            source_normalizer.assert_bound_to(
                self.genome_sequence, source_normalizer.source_sha256
            )
        lineage = transition["lineage"]
        global_status = lineage.get(
            "superlocus_convergence_status", "legacy_unspecified"
        )
        global_reason = lineage.get(
            "superlocus_convergence_reason", "legacy_unspecified"
        )
        if global_status == "converged":
            default_path_local_status = "not_required"
            default_path_local_reason = "superlocus_converged"
        else:
            default_path_local_status = "legacy_unspecified"
            default_path_local_reason = "legacy_unspecified"
        path_local_status = lineage.get(
            "path_local_convergence_status", default_path_local_status
        )
        path_local_reason = lineage.get(
            "path_local_convergence_reason", default_path_local_reason
        )
        path_local_proof = lineage.get("path_local_convergence_proof")
        architecture_commit_eligible = lineage.get(
            "architecture_commit_eligible", True
        )
        if architecture_commit_eligible is not True:
            raise RuntimeError(
                "architecture rescue attempted to commit an ineligible path"
            )

        replayed_path_local_proof = None
        if global_status == "nonconverged":
            if (
                    path_local_status != "converged"
                    or not isinstance(path_local_proof, dict)):
                raise RuntimeError(
                    "nonconverged architecture rescue lacks a local proof"
                )
            try:
                replayed_path_local_proof = (
                    replay_source_superlocus_path_local_convergence_proof(
                        self.genome_sequence,
                        path_local_proof,
                        expected_candidate=row[1],
                        expected_path_id=lineage["path_id"],
                        expected_superlocus_id=lineage["superlocus_id"],
                    )
                )
            except ValueError as error:
                raise RuntimeError(
                    "path-local architecture proof replay failed"
                ) from error
            bounded_path_ids = tuple(sorted(
                record["path_id"]
                for record in replayed_path_local_proof["bounded_paths"]
            ))
            envelope_witness = next(
                record
                for record in replayed_path_local_proof["bounded_paths"]
                if record["path_id"] == lineage["path_id"]
            )
            if (
                    tuple(transition["family_path_ids"])
                    != bounded_path_ids
                    or tuple(lineage.get("parent_path_ids", ()))
                    != tuple(envelope_witness["parent_path_ids"])):
                raise RuntimeError(
                    "path-local architecture lineage differs from replay"
                )
        elif global_status == "converged":
            if (
                    path_local_status != "not_required"
                    or path_local_reason != "superlocus_converged"
                    or path_local_proof is not None):
                raise RuntimeError(
                    "globally converged architecture path claims local replay"
                )
        elif global_status == "legacy_unspecified":
            if (
                    path_local_status != "legacy_unspecified"
                    or path_local_reason != "legacy_unspecified"
                    or path_local_proof is not None):
                raise RuntimeError(
                    "legacy architecture convergence state is inconsistent"
                )
        else:
            raise RuntimeError(
                "architecture rescue convergence status is unsupported"
            )

        phase_arbitration = transition.get("source_phase_arbitration")
        classification_basis = {
            "method": "source_backed_architecture_rescue",
            "source_sha256": source_normalizer.source_sha256,
            "geometry_derivation": "source_backed_augmentation",
            "raw_score_band": (
                "score_below_0.5"
                if score < 0.5
                else (
                    "score_at_least_0.75"
                    if score >= 0.75 else "score_0.5_to_below_0.75"
                )
            ),
            "policy": transition["decision"].policy.as_dict(),
            "decision": transition["decision"].as_dict(),
            "source_superlocus_id": lineage["superlocus_id"],
            "source_superlocus_convergence": {
                "status": global_status,
                "reason": global_reason,
                "path_local_status": path_local_status,
                "path_local_reason": path_local_reason,
                "path_local_proof": path_local_proof,
                "architecture_commit_eligible": (
                    architecture_commit_eligible
                ),
            },
            "lineage_path_id": lineage["path_id"],
            "parent_path_ids": list(
                lineage.get("parent_path_ids", ())
            ),
            "family_path_ids": list(transition["family_path_ids"]),
            "seed_roles": dict(sorted(
                transition["family_seed_roles"].items()
            )),
            "transition_action": transition["action"],
            "transition_reason": transition["reason"],
            "dominated_paths": [
                [list(repeat_interval) for repeat_interval in
                 candidate_repeat_path(dominated_row[1])]
                for dominated_row in transition["dominated_rows"]
            ],
            "target_locus_key": transition["target_locus_key"],
            "migrated_locus_keys": list(
                transition["migrated_locus_keys"]
            ),
            "coding_context": (
                Pipeline._coding_evidence_classification_basis(
                    transition["coding_evidence"]
                )
            ),
            "source_phase_arbitration": phase_arbitration,
        }
        metadata["classification_basis"] = classification_basis
        if replayed_path_local_proof is not None:
            try:
                metadata["boundary_assessment"] = (
                    canonical_source_architecture_rescue_boundary_assessment(
                        self.genome_sequence,
                        source_normalizer.source_sha256,
                        row[1],
                        replayed_path_local_proof,
                        phase_arbitration,
                        expected_path_id=lineage["path_id"],
                        expected_superlocus_id=lineage["superlocus_id"],
                        source_normalizer=source_normalizer,
                    )
                )
            except ReportProvenanceError as error:
                raise RuntimeError(
                    "path-local architecture boundary construction failed"
                ) from error
        return list(row[:3]) + [metadata]

    def _run_source_backed_architecture_rescue_transactionally(self):
        """Add only source-valid, freshly scored, policy-eligible arrays."""

        print("4c. Assemble and adjudicate source-backed array architectures")
        baseline_categories = self.categories
        baseline_rows = self._post_enhancement_scored_rows
        self._source_architecture_baseline_rows = baseline_rows
        stats = dict(self.source_architecture_rescue_stats)
        stats["enabled"] = True
        if not baseline_rows:
            self._final_scored_rows = baseline_rows
            self.source_architecture_rescue_stats = stats
            return

        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        seed_ids, _seed_decisions, seed_roles = (
            select_augmentation_seed_ids(
                self.genome_sequence,
                baseline_rows,
                self.list_features,
                structural_filter,
            )
        )
        stats.update({
            "seed_count": len(seed_ids),
            "primary_seed_count": sum(
                role == "primary" for role in seed_roles.values()
            ),
            "satellite_seed_count": sum(
                role == "satellite" for role in seed_roles.values()
            ),
        })
        if not seed_ids:
            self._final_scored_rows = baseline_rows
            self.source_architecture_rescue_stats = stats
            return

        baseline_candidates = {
            interval: [row[1] for row in rows]
            for interval, rows in baseline_rows.items()
        }
        outputs, lineage, augmentation_stats = (
            augment_source_superloci_to_fixed_point(
                self.genome_sequence,
                baseline_candidates,
                seed_ids,
            )
        )
        self._source_architecture_augmentation_stats = augmentation_stats
        baseline_candidate_ids = {
            id(row[1])
            for rows in baseline_rows.values()
            for row in rows
        }
        generated = {
            superlocus_id: [
                candidate for candidate in candidates
                if id(candidate) not in baseline_candidate_ids
            ]
            for superlocus_id, candidates in outputs.items()
        }
        generated = {
            superlocus_id: candidates
            for superlocus_id, candidates in generated.items()
            if candidates
        }
        stats["generated_candidate_count"] = sum(
            len(candidates) for candidates in generated.values()
        )

        fresh_rows = self._score_candidates_complete(generated)
        if set(fresh_rows) != set(generated):
            raise RuntimeError(
                "architecture scoring returned unexpected superloci"
            )
        for superlocus_id, candidates in generated.items():
            rows = fresh_rows[superlocus_id]
            if len(rows) != len(candidates):
                raise RuntimeError(
                    "architecture scoring changed candidate cardinality"
                )
            for candidate, row in zip(candidates, rows):
                if len(row) < 2 or row[1] is not candidate:
                    raise RuntimeError(
                        "architecture scoring changed candidate order or identity"
                    )
        stats["generated_candidates_scored"] = sum(
            len(rows) for rows in fresh_rows.values()
        )

        adjudication_lineage = tuple(
            entry for entry in lineage
            if entry.get("architecture_commit_eligible", True) is True
        )
        stats["lineage_candidates_quarantined"] = (
            len(lineage) - len(adjudication_lineage)
        )

        rows_by_candidate_id = self._rows_by_candidate_id(baseline_rows)
        for candidate_id, row in self._rows_by_candidate_id(
                fresh_rows).items():
            if candidate_id in rows_by_candidate_id:
                raise RuntimeError(
                    "generated architecture reused a baseline candidate object"
                )
            rows_by_candidate_id[candidate_id] = row

        preliminary_decisions = {}
        context_candidates = []
        for entry in adjudication_lineage:
            candidate = entry["candidate"]
            row = rows_by_candidate_id.get(id(candidate))
            if row is None:
                raise RuntimeError("lineage candidate has no scored row")
            decision = evaluate_scored_architecture_candidate(
                self.genome_sequence,
                row,
                self.list_features,
                structural_filter,
            )
            preliminary_decisions[id(candidate)] = decision
            if decision.eligible:
                context_candidates.append(candidate)

        coding_by_candidate_id = {}
        if context_candidates:
            coding_batch = batch_source_coding_context_evidence(
                self.genome_sequence,
                context_candidates,
            )
            source_cds = getattr(coding_batch, "cds", None)
            if isinstance(source_cds, tuple):
                self._source_cds = source_cds
            coding_by_candidate_id = (
                self._validated_coding_evidence_by_candidate_id(
                    self.genome_sequence,
                    context_candidates,
                    coding_batch,
                    preliminary_decisions,
                )
            )
        stats["context_candidates_checked"] = len(context_candidates)

        final_decisions = {}
        for entry in adjudication_lineage:
            candidate = entry["candidate"]
            candidate_id = id(candidate)
            coding_evidence = coding_by_candidate_id.get(candidate_id)
            if coding_evidence is None:
                final_decisions[candidate_id] = preliminary_decisions[
                    candidate_id
                ]
                continue
            final_decisions[candidate_id] = (
                evaluate_scored_architecture_candidate(
                    self.genome_sequence,
                    rows_by_candidate_id[candidate_id],
                    self.list_features,
                    structural_filter,
                    coding_evidence=(
                        coding_evidence.coding_overlap_fraction,
                        coding_evidence.context_cds_score,
                    ),
                )
            )

        proposals = select_lineage_family_proposals(
            lineage,
            rows_by_candidate_id,
            final_decisions,
            seed_roles,
        )
        proposals = tuple({
            **proposal,
            "coding_evidence": coding_by_candidate_id[
                id(proposal["row"][1])
            ],
        } for proposal in proposals)
        stats["proposals_considered"] = len(proposals)
        proposals, phase_arbitration = (
            arbitrate_cross_family_phase_proposals(
                self.genome_sequence,
                proposals,
            )
        )
        self.source_architecture_phase_arbitration = phase_arbitration
        stats.update({
            "proposals_after_phase_arbitration": len(proposals),
            "phase_conflicts_selected": sum(
                record["status"] == "selected"
                for record in phase_arbitration
            ),
            "phase_conflicts_ambiguous_representative": sum(
                record["status"] == "ambiguous_representative"
                for record in phase_arbitration
            ),
            "phase_conflicts_quarantined": sum(
                record["status"] == "quarantined"
                for record in phase_arbitration
            ),
            "phase_proposals_removed": sum(
                len(record["removed_path_ids"])
                for record in phase_arbitration
            ),
        })
        transitions = plan_architecture_rescue_transitions(
            baseline_categories,
            proposals,
        )

        final_categories = self._copy_categories(baseline_categories)
        committed_rows = []
        removed_accepted_row_ids = set()
        source_normalizer = BoundaryAssessmentNormalizer.from_source(
            self.genome_sequence
        )
        for transition in transitions:
            action = transition["action"]
            if action not in ("add", "replace"):
                continue
            row_ids_to_remove = {id(transition["row"])}
            row_ids_to_remove.update(
                id(row) for row in transition["dominated_rows"]
            )
            self._migrate_category_locus_keys(
                final_categories,
                transition["target_locus_key"],
                transition["migrated_locus_keys"],
            )
            self._remove_category_rows_by_identity(
                final_categories, row_ids_to_remove
            )
            committed_row = self._rescued_row_with_basis(
                transition, source_normalizer=source_normalizer
            )
            final_categories[transition["target_category_index"]].setdefault(
                transition["target_locus_key"], []
            ).append(committed_row)
            committed_rows.append((transition, committed_row))
            removed_accepted_row_ids.update(
                id(row) for row in transition["dominated_rows"]
            )

        baseline_accepted = self._accepted_candidate_ids(
            baseline_categories
        )
        expected_accepted = {
            candidate_id
            for category_index in (0, 2)
            for rows in baseline_categories[category_index].values()
            for row in rows
            for candidate_id in (id(row[1]),)
            if id(row) not in removed_accepted_row_ids
        }
        expected_accepted.update(
            id(committed_row[1])
            for _transition, committed_row in committed_rows
        )
        final_accepted = self._accepted_candidate_ids(final_categories)
        if final_accepted != expected_accepted:
            raise RuntimeError(
                "architecture rescue violated the accepted-call transaction"
            )
        if not baseline_accepted.issubset(
                final_accepted.union({
                    id(row[1])
                    for transition, _row in committed_rows
                    for row in transition["dominated_rows"]
                })):
            raise RuntimeError(
                "architecture rescue removed an unrelated accepted call"
            )

        stats.update({
            "calls_added": sum(
                transition["action"] == "add"
                for transition in transitions
            ),
            "calls_replaced": sum(
                transition["action"] == "replace"
                for transition in transitions
            ),
            "calls_abstained": sum(
                transition["action"] == "abstain"
                for transition in transitions
            ),
            "augmentation": augmentation_stats,
        })
        final_rows = self._copy_scored_rows(baseline_rows)
        for superlocus_id, rows in fresh_rows.items():
            if superlocus_id in final_rows:
                raise RuntimeError(
                    "architecture superlocus collides with a legacy score key"
                )
            final_rows[superlocus_id] = list(rows)
        self.source_architecture_rescue_lineage = tuple(lineage)
        self.source_architecture_rescue_transitions = tuple(transitions)
        self.source_architecture_rescue_stats = stats
        self._final_scored_rows = final_rows
        self.categories = final_categories

    @staticmethod
    def _candidate_family_sort_key(candidate):
        """Return a source-geometry-first deterministic family order."""

        try:
            starts = tuple(candidate.list_repeat_starts)
            repeats = tuple(candidate.list_repeats)
            spacers = tuple(candidate.list_spacers)
        except (AttributeError, TypeError) as error:
            raise RuntimeError(
                "accepted-family candidate is malformed"
            ) from error
        if not starts or len(starts) != len(repeats):
            raise RuntimeError(
                "accepted-family candidate has inconsistent repeat geometry"
            )
        return (
            tuple(
                (start, start + len(repeat))
                for start, repeat in zip(starts, repeats)
            ),
            repeats,
            spacers,
        )

    @staticmethod
    def _accepted_row_family_contexts(
            categories, baseline_rows, lineage, augmentation_stats,
            completed_family_records=()):
        """Attest complete runtime families for every accepted row.

        Baseline candidates from each legacy locus are always family peers.
        Source-superlocus records then join every legacy locus that the
        augmentation graph connected. Finally, every eligible, generated,
        and quarantined lineage candidate is included. This deliberately
        uses the complete candidate universe, not only rescue proposals.
        """

        if not isinstance(baseline_rows, dict):
            raise RuntimeError(
                "accepted-family baseline rows must be a mapping"
            )
        try:
            lineage = tuple(lineage)
        except TypeError as error:
            raise RuntimeError(
                "accepted-family lineage must be iterable"
            ) from error
        if not isinstance(augmentation_stats, dict):
            raise RuntimeError(
                "accepted-family augmentation stats must be a mapping"
            )

        objects = {}
        parent = {}
        order = []

        def add(candidate):
            candidate_id = id(candidate)
            previous = objects.get(candidate_id)
            if previous is not None and previous is not candidate:
                raise RuntimeError(
                    "accepted-family candidate identity collision"
                )
            if candidate_id not in objects:
                objects[candidate_id] = candidate
                parent[candidate_id] = candidate_id
                order.append(candidate_id)
            return candidate_id

        def find(candidate_id):
            root = candidate_id
            while parent[root] != root:
                root = parent[root]
            while parent[candidate_id] != candidate_id:
                next_id = parent[candidate_id]
                parent[candidate_id] = root
                candidate_id = next_id
            return root

        def join(candidates):
            candidate_ids = [add(candidate) for candidate in candidates]
            if not candidate_ids:
                return
            root = find(candidate_ids[0])
            for candidate_id in candidate_ids[1:]:
                other_root = find(candidate_id)
                if other_root != root:
                    parent[other_root] = root

        candidates_by_locus = {}
        for locus_key, rows in baseline_rows.items():
            if not isinstance(rows, (list, tuple)) or not rows:
                raise RuntimeError(
                    "accepted-family baseline loci need non-empty rows"
                )
            candidates = []
            for row in rows:
                if (
                    not isinstance(row, (list, tuple))
                    or len(row) < 2
                    or row[1] is None
                ):
                    raise RuntimeError(
                        "accepted-family baseline row is malformed"
                    )
                candidates.append(row[1])
            candidates_by_locus[locus_key] = tuple(candidates)
            join(candidates)

        lineage_by_superlocus = {}
        lineage_loci_by_superlocus = {}
        for entry in lineage:
            if not isinstance(entry, dict):
                raise RuntimeError(
                    "accepted-family lineage entries must be mappings"
                )
            candidate = entry.get("candidate")
            superlocus_id = entry.get("superlocus_id")
            if candidate is None or not isinstance(superlocus_id, str):
                raise RuntimeError(
                    "accepted-family lineage entry is incomplete"
                )
            lineage_by_superlocus.setdefault(
                superlocus_id, []
            ).append(candidate)
            details = entry.get("details", {})
            if details is None:
                details = {}
            if not isinstance(details, dict):
                raise RuntimeError(
                    "accepted-family lineage details must be a mapping"
                )
            legacy_locus_keys = details.get("legacy_locus_keys", ())
            if legacy_locus_keys is None:
                legacy_locus_keys = ()
            try:
                legacy_locus_keys = tuple(legacy_locus_keys)
            except TypeError as error:
                raise RuntimeError(
                    "accepted-family legacy locus keys must be iterable"
                ) from error
            lineage_loci_by_superlocus.setdefault(
                superlocus_id, []
            ).extend(legacy_locus_keys)

        superlocus_records = augmentation_stats.get("superloci", ())
        if superlocus_records is None:
            superlocus_records = ()
        try:
            superlocus_records = tuple(superlocus_records)
        except TypeError as error:
            raise RuntimeError(
                "accepted-family superlocus stats must be iterable"
            ) from error

        observed_superloci = set()
        for record in superlocus_records:
            if not isinstance(record, dict):
                raise RuntimeError(
                    "accepted-family superlocus records must be mappings"
                )
            superlocus_id = record.get("superlocus_id")
            if not isinstance(superlocus_id, str):
                raise RuntimeError(
                    "accepted-family superlocus ID is missing"
                )
            if superlocus_id in observed_superloci:
                raise RuntimeError(
                    "accepted-family superlocus ID is duplicated"
                )
            observed_superloci.add(superlocus_id)
            legacy_locus_keys = record.get("legacy_locus_keys", ())
            try:
                legacy_locus_keys = tuple(legacy_locus_keys)
            except TypeError as error:
                raise RuntimeError(
                    "accepted-family legacy locus keys must be iterable"
                ) from error
            members = list(
                lineage_by_superlocus.get(superlocus_id, ())
            )
            for locus_key in legacy_locus_keys:
                if locus_key not in candidates_by_locus:
                    raise RuntimeError(
                        "accepted-family superlocus refers to an "
                        "unknown baseline locus"
                    )
                members.extend(candidates_by_locus[locus_key])
            join(members)

        for superlocus_id, candidates in lineage_by_superlocus.items():
            members = list(candidates)
            for locus_key in lineage_loci_by_superlocus.get(
                    superlocus_id, ()):
                if locus_key not in candidates_by_locus:
                    raise RuntimeError(
                        "accepted-family lineage refers to an "
                        "unknown baseline locus"
                    )
                members.extend(candidates_by_locus[locus_key])
            join(members)

        try:
            completed_family_records = tuple(completed_family_records)
        except TypeError as error:
            raise RuntimeError(
                "accepted-core completed family records must be iterable"
            ) from error
        completed_candidate_ids = set()
        for record in completed_family_records:
            if not isinstance(record, dict):
                raise RuntimeError(
                    "accepted-core completed family records must be mappings"
                )
            completed_candidate = record.get("candidate")
            authorizing_core = record.get("authorizing_core")
            members = record.get("complete_family_candidates")
            if not isinstance(members, tuple):
                raise RuntimeError(
                    "accepted-core completed family must be an explicit tuple"
                )
            member_ids = tuple(id(member) for member in members)
            if (
                    completed_candidate is None
                    or authorizing_core is None
                    or len(set(member_ids)) != len(member_ids)
                    or not any(
                        member is completed_candidate for member in members
                    )
                    or not any(
                        member is authorizing_core for member in members
                    )
                    or id(completed_candidate) in completed_candidate_ids):
                raise RuntimeError(
                    "accepted-core completed family record is inconsistent"
                )
            if id(authorizing_core) not in objects:
                raise RuntimeError(
                    "accepted-core completion refers to an unknown core family"
                )
            if any(
                    member is not completed_candidate
                    and id(member) not in objects
                    for member in members):
                raise RuntimeError(
                    "accepted-core completion added an unknown family member"
                )
            completed_candidate_ids.add(id(completed_candidate))
            join(members)

        component_members = {}
        for candidate_id in order:
            component_members.setdefault(
                find(candidate_id), []
            ).append(objects[candidate_id])
        for members in component_members.values():
            members.sort(key=Pipeline._candidate_family_sort_key)

        accepted_rows = []
        if not isinstance(categories, (list, tuple)) or len(categories) != 5:
            raise RuntimeError(
                "accepted-family categories must contain five mappings"
            )
        for category_index in (0, 2):
            category = categories[category_index]
            if not isinstance(category, dict):
                raise RuntimeError(
                    "accepted-family accepted categories must be mappings"
                )
            for rows in category.values():
                if not isinstance(rows, (list, tuple)) or not rows:
                    raise RuntimeError(
                        "accepted-family category rows must be non-empty"
                    )
                accepted_rows.extend(rows)

        contexts = []
        for row in accepted_rows:
            if (
                not isinstance(row, (list, tuple))
                or len(row) < 2
                or row[1] is None
            ):
                raise RuntimeError(
                    "accepted-family accepted row is malformed"
                )
            candidate_id = id(row[1])
            if candidate_id not in objects:
                raise RuntimeError(
                    "accepted candidate is absent from the complete "
                    "baseline/lineage family universe"
                )
            family = tuple(component_members[find(candidate_id)])
            contexts.append(AcceptedRowFamilyContext(row, family))
        return tuple(contexts)

    def _demote_accepted_family_tandems_transactionally(self):
        """Demote only source-CDS-supported coding tandem families."""

        if not self._source_architecture_rescue_enabled():
            return False

        print("4d. Validate accepted families against source coding tandems")
        contexts = self._accepted_row_family_contexts(
            self.categories,
            self._source_architecture_baseline_rows,
            self.source_architecture_rescue_lineage,
            self._source_architecture_augmentation_stats,
            getattr(
                self,
                "_source_accepted_core_boundary_family_records",
                (),
            ),
        )
        source_cds = self._source_cds
        if source_cds is None:
            source_cds = (
                predict_source_cds(self.genome_sequence)
                if contexts
                else ()
            )
            self._source_cds = source_cds
        categories, records, stats = (
            apply_accepted_family_tandem_demotions(
                self.genome_sequence,
                self.categories,
                contexts,
                source_cds,
            )
        )
        if len(contexts) != len(records):
            raise RuntimeError(
                "accepted-family tandem decision cardinality changed"
            )
        self._source_accepted_family_tandem_by_candidate_id = {
            id(context.accepted_row[1]): record.decision
            for context, record in zip(contexts, records)
        }
        self.categories = self._attach_retained_tandem_reviews(
            categories,
            contexts,
            records,
        )
        self.source_accepted_family_tandem_records = records
        self.source_accepted_family_tandem_stats = {
            "enabled": True,
            **stats.as_dict(),
        }
        return stats.tandem_supported > 0

    def _accepted_core_boundary_contexts(self, structural_filter):
        """Build complete source-family witnesses for every retained core."""

        family_contexts = self._accepted_row_family_contexts(
            self.categories,
            self._source_architecture_baseline_rows,
            self.source_architecture_rescue_lineage,
            self._source_architecture_augmentation_stats,
            getattr(
                self,
                "_source_accepted_core_boundary_family_records",
                (),
            ),
        )
        accepted_candidates = {
            id(row[1]): row[1]
            for category_index in (0, 2)
            for rows in self.categories[category_index].values()
            for row in rows
        }
        contexts = []
        families_by_row_id = {}
        for family_context in family_contexts:
            row = family_context.accepted_row
            candidate = row[1]
            tandem_decision = (
                self._source_accepted_family_tandem_by_candidate_id.get(
                    id(candidate)
                )
            )
            if tandem_decision is None:
                raise RuntimeError(
                    "accepted core is missing complete-family tandem evidence"
                )
            if tandem_decision.supported:
                raise RuntimeError(
                    "a coding-tandem-supported row remained accepted"
                )

            core_decision = evaluate_scored_architecture_candidate(
                self.genome_sequence,
                row,
                self.list_features,
                structural_filter,

            )
            family_by_id = {
                id(member): member
                for member in family_context.complete_family_candidates
            }
            family_by_id.pop(id(candidate), None)
            ordered_witnesses = sorted(
                family_by_id.values(),
                key=self._candidate_family_sort_key,
            )
            witnesses = tuple(
                BoundaryGeometryWitness(
                    witness,
                    "boundary-family-{:04d}-{}-{}-r{}".format(
                        witness_index,
                        witness.list_repeat_starts[0],
                        witness.list_repeat_starts[-1]
                        + len(witness.list_repeats[-1]),
                        len(witness.list_repeats),
                    ),
                    accepted=id(witness) in accepted_candidates,
                )
                for witness_index, witness in enumerate(ordered_witnesses)
            )
            contexts.append(AcceptedCoreBoundaryContext(
                accepted_row=row,
                core_decision=core_decision,
                core_authorization=CoreAuthorizationEvidence(
                    coding_context_assessed=True,
                    tandem_evidence_assessed=True,
                    coding_tandem_clear=True,
                ),
                witnesses=witnesses,
                witness_context=WitnessContextEvidence(
                    source_family_complete=True,
                    accepted_statuses_complete=True,
                ),
            ))
            families_by_row_id[id(row)] = tuple(
                family_context.complete_family_candidates
            )
        return tuple(contexts), families_by_row_id

    @staticmethod
    def _boundary_phase_completeness_dict(completeness):
        return {
            "source_family_complete": completeness.source_family_complete,
            "lineage_complete": completeness.lineage_complete,
            "same_occurrence_alternatives_complete": (
                completeness.same_occurrence_alternatives_complete
            ),
            "augmentation_converged": completeness.augmentation_converged,
        }

    @staticmethod
    def _boundary_phase_lineage_dict(lineage_entries):
        records = []
        for entry in lineage_entries:
            geometry = entry.source_geometry
            records.append({
                "path_id": entry.path_id,
                "family_id": entry.family_id,
                "derivation": entry.derivation,
                "parent_path_ids": list(entry.parent_path_ids),
                "coordinate_system": "zero_based_half_open",
                "source_interval": list(geometry.source_interval),
                "repeat_count": geometry.repeat_count,
                "repeat_intervals": [
                    list(interval) for interval in geometry.repeat_intervals
                ],
            })
        return records

    def _accepted_core_boundary_phase_inputs(
            self, request, draft, complete_family_candidates):
        """Build source-replayed lineage and explicit completeness evidence.

        Completeness is derived only from the augmentation's bounded runtime
        universe and convergence records. Missing, inconsistent, or
        non-converged state produces false flags; the pure attestor then fails
        closed without granting exact boundaries.
        """

        try:
            family = tuple(complete_family_candidates)
        except TypeError:
            family = ()
        candidates = (
            family
            if any(member is request.candidate for member in family)
            else family + (request.candidate,)
        )
        candidate_ids = {id(candidate) for candidate in candidates}
        family_candidate_ids = {
            id(candidate) for candidate in family
        }

        reconstruction_complete = True
        observed_geometries = {}
        for candidate in candidates:
            try:
                observed_geometries[id(candidate)] = (
                    SourceGeometry.from_candidate(
                        self.genome_sequence, candidate
                    )
                )
            except Exception:
                reconstruction_complete = False

        try:
            runtime_lineage = tuple(
                self.source_architecture_rescue_lineage
            )
        except (AttributeError, TypeError):
            runtime_lineage = ()
        augmentation_stats = getattr(
            self, "_source_architecture_augmentation_stats", None
        )
        if not isinstance(augmentation_stats, dict):
            augmentation_stats = {}

        runtime_lineage_valid = True
        runtime_by_candidate_id = {}
        runtime_path_ids = set()
        for entry in runtime_lineage:
            if not isinstance(entry, dict):
                runtime_lineage_valid = False
                continue
            candidate = entry.get("candidate")
            path_id = entry.get("path_id")
            superlocus_id = entry.get("superlocus_id")
            if (
                    candidate is None
                    or not isinstance(path_id, str)
                    or not path_id
                    or not isinstance(superlocus_id, str)
                    or not superlocus_id
                    or path_id in runtime_path_ids
                    or id(candidate) in runtime_by_candidate_id):
                runtime_lineage_valid = False
                continue
            runtime_path_ids.add(path_id)
            runtime_by_candidate_id[id(candidate)] = entry

        try:
            completed_family_records = tuple(getattr(
                self,
                "_source_accepted_core_boundary_family_records",
                (),
            ))
        except TypeError:
            completed_family_records = ()
            runtime_lineage_valid = False
        completed_by_candidate_id = {}
        for record in completed_family_records:
            if not isinstance(record, dict):
                runtime_lineage_valid = False
                continue
            completed_candidate = record.get("candidate")
            authorizing_core = record.get("authorizing_core")
            request_id = record.get("request_id")
            members = record.get("complete_family_candidates")
            if (
                    completed_candidate is None
                    or authorizing_core is None
                    or not isinstance(request_id, str)
                    or not request_id
                    or not isinstance(members, tuple)
                    or not any(
                        member is completed_candidate for member in members
                    )
                    or not any(
                        member is authorizing_core for member in members
                    )
                    or id(completed_candidate)
                    in completed_by_candidate_id):
                runtime_lineage_valid = False
                continue
            completed_by_candidate_id[id(completed_candidate)] = record

        family_runtime_entries = tuple(
            runtime_by_candidate_id[candidate_id]
            for candidate_id in family_candidate_ids
            if candidate_id in runtime_by_candidate_id
        )
        relevant_superlocus_ids = {
            entry["superlocus_id"] for entry in family_runtime_entries
        }

        superlocus_records = augmentation_stats.get("superloci", ())
        if not isinstance(superlocus_records, (list, tuple)):
            superlocus_records = ()
            runtime_lineage_valid = False
        records_by_superlocus = {}
        for record in superlocus_records:
            if (
                    not isinstance(record, dict)
                    or not isinstance(record.get("superlocus_id"), str)
                    or record["superlocus_id"] in records_by_superlocus):
                runtime_lineage_valid = False
                continue
            records_by_superlocus[record["superlocus_id"]] = record

        relevant_superlocus_id = (
            next(iter(relevant_superlocus_ids))
            if len(relevant_superlocus_ids) == 1
            else None
        )
        relevant_record = records_by_superlocus.get(
            relevant_superlocus_id
        )
        relevant_record_values = (
            relevant_record if relevant_record is not None else {}
        )

        baseline_rows = getattr(
            self, "_source_architecture_baseline_rows", None
        )
        baseline_valid = isinstance(baseline_rows, dict)
        baseline_count = 0
        baseline_loci_by_candidate_id = {}
        if baseline_valid:
            for locus_key, rows in baseline_rows.items():
                if not isinstance(rows, (list, tuple)) or not rows:
                    baseline_valid = False
                    continue
                baseline_count += len(rows)
                for row in rows:
                    if (
                            not isinstance(row, (list, tuple))
                            or len(row) < 2
                            or row[1] is None):
                        baseline_valid = False
                        continue
                    baseline_loci_by_candidate_id.setdefault(
                        id(row[1]), set()
                    ).add(locus_key)

        stats_counts_valid = (
            baseline_valid
            and augmentation_stats.get("source_validated_candidate_count")
            == baseline_count
            and augmentation_stats.get("lineage_candidate_count")
            == len(runtime_lineage)
            and augmentation_stats.get("superlocus_count")
            == len(superlocus_records)
        )

        relevant_runtime_entries = tuple(
            entry for entry in runtime_by_candidate_id.values()
            if entry["superlocus_id"] == relevant_superlocus_id
        )
        relevant_runtime_ids = {
            id(entry["candidate"]) for entry in relevant_runtime_entries
        }
        bounded_count_matches = (
            relevant_record is not None
            and relevant_record.get("bounded_candidate_count")
            == len(relevant_runtime_entries)
        )
        family_contains_runtime_universe = (
            relevant_superlocus_id is not None
            and relevant_runtime_ids.issubset(family_candidate_ids)
        )
        record_locus_keys = (
            set(relevant_record.get("legacy_locus_keys", ()))
            if relevant_record is not None
            and isinstance(
                relevant_record.get("legacy_locus_keys", ()),
                (list, tuple),
            )
            else set()
        )
        family_has_runtime_or_baseline_origin = all(
            candidate_id in runtime_by_candidate_id
            or candidate_id in baseline_loci_by_candidate_id
            or candidate_id in completed_by_candidate_id
            for candidate_id in family_candidate_ids
        )
        family_baseline_loci_match = all(
            baseline_loci_by_candidate_id.get(candidate_id, set())
            .issubset(record_locus_keys)
            for candidate_id in family_candidate_ids
        )

        source_family_complete = all((
            reconstruction_complete,
            stats_counts_valid,
            runtime_lineage_valid,
            relevant_record is not None,
            family_contains_runtime_universe,
            bounded_count_matches,
            family_has_runtime_or_baseline_origin,
            family_baseline_loci_match,
        ))

        convergence_counts = tuple(
            augmentation_stats.get(field_name)
            for field_name in (
                "confirmed_converged_superlocus_count",
                "nonconverged_superlocus_count",
                "commit_eligible_superlocus_count",
            )
        )
        convergence_counts_valid = all(
            isinstance(value, int)
            and not isinstance(value, bool)
            and value >= 0
            for value in convergence_counts
        )
        global_convergence_counts_consistent = (
            stats_counts_valid
            and convergence_counts_valid
            and convergence_counts[0] + convergence_counts[1]
            == len(superlocus_records)
            and convergence_counts[2] == convergence_counts[0]
        )
        augmentation_converged = all((
            global_convergence_counts_consistent,
            relevant_record is not None,
            relevant_record_values.get("convergence_status") == "converged",
            relevant_record_values.get(
                "architecture_commit_eligible"
            ) is True,
            relevant_record_values.get(
                "confirmation_semantic_state_changed"
            ) is False,
            relevant_record_values.get(
                "confirmation_terminal_geometries_added"
            ) == 0,
            relevant_record_values.get(
                "confirmation_shared_geometries_added"
            ) == 0,
            all(
                entry.get("superlocus_convergence_status") == "converged"
                and entry.get("architecture_commit_eligible") is True
                for entry in relevant_runtime_entries
            ),
        ))

        family_id = "runtime-source-family:{}".format(
            relevant_superlocus_id or "unresolved"
        )
        lineage_entries = []
        built_path_ids = set()
        lineage_mapping_complete = (
            reconstruction_complete and runtime_lineage_valid
        )
        ordered_candidates = sorted(
            candidates,
            key=self._candidate_family_sort_key,
        )
        path_id_by_candidate_id = {}
        for candidate_index, candidate in enumerate(ordered_candidates):
            if candidate is request.candidate:
                continue
            runtime_entry = runtime_by_candidate_id.get(id(candidate))
            completion_record = completed_by_candidate_id.get(id(candidate))
            if completion_record is not None:
                path_id = "{}:proposed".format(
                    completion_record["request_id"]
                )
            elif runtime_entry is None:
                geometry = observed_geometries.get(id(candidate))
                if geometry is None:
                    lineage_mapping_complete = False
                    continue
                path_id = "legacy-observed-{:04d}-{}-{}-r{}".format(
                    candidate_index,
                    geometry.source_interval[0],
                    geometry.source_interval[1],
                    geometry.repeat_count,
                )
            else:
                path_id = runtime_entry["path_id"]
            if path_id in built_path_ids:
                lineage_mapping_complete = False
            built_path_ids.add(path_id)
            path_id_by_candidate_id[id(candidate)] = path_id

        for candidate in ordered_candidates:
            if candidate is request.candidate:
                continue
            geometry = observed_geometries.get(id(candidate))
            path_id = path_id_by_candidate_id.get(id(candidate))
            if geometry is None or path_id is None:
                lineage_mapping_complete = False
                continue
            runtime_entry = runtime_by_candidate_id.get(id(candidate))
            completion_record = completed_by_candidate_id.get(id(candidate))
            if completion_record is not None:
                derivation = "accepted_core_boundary_completion"
                parent_path_id = path_id_by_candidate_id.get(
                    id(completion_record["authorizing_core"])
                )
                if parent_path_id is None:
                    parent_path_ids = ()
                    lineage_mapping_complete = False
                else:
                    parent_path_ids = (parent_path_id,)
            elif runtime_entry is None:
                derivation = "legacy_observed"
                parent_path_ids = ()
            else:
                derivation = runtime_entry.get("stage")
                raw_parent_ids = runtime_entry.get("parent_path_ids")
                if not isinstance(raw_parent_ids, (list, tuple)):
                    raw_parent_ids = ()
                    lineage_mapping_complete = False
                parent_path_ids = tuple(raw_parent_ids)
                if (
                        derivation not in KNOWN_LINEAGE_DERIVATIONS
                        or any(
                            not isinstance(parent_id, str)
                            or not parent_id
                            for parent_id in parent_path_ids
                        )):
                    lineage_mapping_complete = False
                try:
                    runtime_geometry = tuple(
                        tuple(interval)
                        for interval in runtime_entry["source_geometry"]
                    )
                except (KeyError, TypeError, ValueError):
                    runtime_geometry = None
                if runtime_geometry != geometry.repeat_intervals:
                    lineage_mapping_complete = False
            lineage_entries.append(BoundaryPhaseLineageEntry(
                candidate=candidate,
                source_geometry=geometry,
                path_id=path_id,
                family_id=family_id,
                derivation=derivation,
                parent_path_ids=parent_path_ids,
            ))

        proposed_geometry = observed_geometries.get(id(request.candidate))
        parent_path_id = path_id_by_candidate_id.get(id(draft.core_row[1]))
        proposed_path_id = "{}:proposed".format(request.request_id)
        if (
                proposed_geometry is None
                or parent_path_id is None
                or proposed_path_id in built_path_ids):
            lineage_mapping_complete = False
        else:
            lineage_entries.append(BoundaryPhaseLineageEntry(
                candidate=request.candidate,
                source_geometry=proposed_geometry,
                path_id=proposed_path_id,
                family_id=family_id,
                derivation="accepted_core_boundary_completion",
                parent_path_ids=(parent_path_id,),
            ))
            built_path_ids.add(proposed_path_id)

        all_parent_ids_known = all(
            set(entry.parent_path_ids).issubset(built_path_ids)
            for entry in lineage_entries
        )
        lineage_complete = all((
            source_family_complete,
            lineage_mapping_complete,
            len(lineage_entries) == len(candidate_ids),
            len(built_path_ids) == len(lineage_entries),
            all_parent_ids_known,
        ))
        completeness = BoundaryPhaseCompletenessEvidence(
            source_family_complete=source_family_complete,
            lineage_complete=lineage_complete,
            same_occurrence_alternatives_complete=all((
                source_family_complete,
                lineage_complete,
                augmentation_converged,
                bounded_count_matches,
            )),
            augmentation_converged=augmentation_converged,
        )
        return tuple(lineage_entries), completeness

    @classmethod
    def _boundary_phase_provenance(
            cls, attestation, completeness, lineage_entries):
        value = attestation.as_dict()
        value["completeness_evidence"] = (
            cls._boundary_phase_completeness_dict(completeness)
        )
        value["lineage"] = cls._boundary_phase_lineage_dict(
            lineage_entries
        )
        return value

    @staticmethod
    def _scored_row_with_boundary_phase_provenance(row, provenance):
        if not isinstance(row, (list, tuple)) or len(row) not in (3, 4):
            raise RuntimeError("boundary scorer returned a malformed row")
        metadata = dict(row[3]) if len(row) == 4 else {}
        metadata["accepted_core_boundary_phase_attestation"] = provenance
        fields = list(row[:3]) + [metadata]
        return tuple(fields) if isinstance(row, tuple) else fields

    @staticmethod
    def _attach_boundary_phase_classification_provenance(
            categories, batch, provenance_by_request_id):
        annotated = [
            {locus_key: list(rows) for locus_key, rows in category.items()}
            for category in categories
        ]
        for draft in batch.drafts:
            if draft.score_request is None:
                continue
            provenance = provenance_by_request_id[
                draft.score_request.request_id
            ]
            rows = annotated[draft.category_index][
                draft._runtime_locus_key
            ]
            row = rows[draft.row_ordinal]
            metadata = dict(row[3]) if len(row) == 4 else {}
            basis = metadata.get("classification_basis")
            if not isinstance(basis, dict):
                raise RuntimeError(
                    "boundary phase provenance needs a classification basis"
                )
            basis = dict(basis)
            basis["accepted_core_boundary_phase_attestation"] = provenance
            metadata["classification_basis"] = basis
            fields = list(row[:3]) + [metadata]
            rows[draft.row_ordinal] = (
                tuple(fields) if isinstance(row, tuple) else fields
            )
        return annotated

    def _complete_accepted_core_boundaries_transactionally(self):
        """Extend strict accepted cores without re-deciding their existence."""

        if not self._source_architecture_rescue_enabled():
            return False

        print("4e. Complete accepted-core boundaries from source evidence")
        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        contexts, families_by_row_id = (
            self._accepted_core_boundary_contexts(structural_filter)
        )
        source_cds = self._source_cds
        if source_cds is None:
            source_cds = (
                predict_source_cds(self.genome_sequence)
                if contexts
                else ()
            )
            self._source_cds = source_cds
        batch = plan_accepted_core_boundary_completion_batch(
            self.genome_sequence,
            self.categories,
            contexts,
        )

        fresh_rows = []
        if batch.score_requests:
            score_key = "accepted-core-boundary-batch"
            scored = self._score_candidates_complete({
                score_key: list(batch.proposed_candidates),
            })
            if set(scored) != {score_key}:
                raise RuntimeError(
                    "boundary scoring returned an unexpected batch"
                )
            fresh_rows = list(scored[score_key])
            if len(fresh_rows) != len(batch.score_requests):
                raise RuntimeError(
                    "boundary scoring changed request cardinality"
                )

        draft_by_request_id = {
            draft.score_request.request_id: draft
            for draft in batch.drafts
            if draft.score_request is not None
        }
        results = []
        phase_attestations = []
        phase_provenance_by_request_id = {}
        for request, row in zip(batch.score_requests, fresh_rows):
            if len(row) < 2 or row[1] is not request.candidate:
                raise RuntimeError(
                    "boundary scoring changed candidate order or identity"
                )
            draft = draft_by_request_id[request.request_id]
            complete_family = tuple(
                families_by_row_id[id(draft.core_row)]
            )
            lineage_entries, completeness = (
                self._accepted_core_boundary_phase_inputs(
                    request,
                    draft,
                    complete_family,
                )
            )
            phase_attestation = attest_accepted_core_boundary_phase(
                self.genome_sequence,
                request.candidate,
                complete_family,
                lineage_entries,
                completeness,
            )
            phase_provenance = self._boundary_phase_provenance(
                phase_attestation,
                completeness,
                lineage_entries,
            )
            phase_provenance_by_request_id[
                request.request_id
            ] = phase_provenance
            phase_attestations.append({
                "request_id": request.request_id,
                "category_index": draft.category_index,
                "category_name": draft.category_name,
                "locus_key": draft.locus_key,
                "row_ordinal": draft.row_ordinal,
                "phase_attestation": phase_provenance,
            })

            family = list(complete_family)
            if all(member is not request.candidate for member in family):
                family.append(request.candidate)
            tandem_decision = validate_accepted_family_tandem(
                self.genome_sequence,
                request.candidate,
                tuple(family),
                source_cds,
            )
            structural_evidence = structural_filter.evidence(
                request.candidate
            )
            reason_codes = []
            if not structural_evidence.passes:
                reason_codes.append("completed_structural_filter_rejected")
            if tandem_decision.supported:
                reason_codes.append("completed_family_coding_tandem_supported")
            if not phase_attestation.phase_unambiguous:
                reason_codes.extend(
                    "boundary_phase_attestation_{}".format(reason)
                    for reason in phase_attestation.reason_codes
                )
            scored_row = self._scored_row_with_boundary_phase_provenance(
                row, phase_provenance
            )
            results.append(BoundaryCompletionScoreResult(
                request_id=request.request_id,
                scored_row=scored_row,
                final_evidence=FinalBoundaryEvidence(
                    fresh_score=float(row[0]),
                    source_reconstruction_valid=True,
                    structural_valid=structural_evidence.passes,
                    phase_unambiguous=(
                        phase_attestation.phase_unambiguous
                    ),
                    coding_context_assessed=True,
                    tandem_evidence_assessed=True,
                    coding_tandem_clear=not tandem_decision.supported,
                    reason_codes=tuple(reason_codes),
                ),
            ))

        categories, records, stats = (
            apply_accepted_core_boundary_completion_batch(
                self.genome_sequence,
                self.categories,
                batch,
                tuple(results),
            )
        )
        completed_family_records = []
        for draft, record in zip(batch.drafts, records):
            if record.action != "completed":
                continue
            request = draft.score_request
            family = list(families_by_row_id[id(draft.core_row)])
            if all(member is not request.candidate for member in family):
                family.append(request.candidate)
            completed_family_records.append({
                "request_id": request.request_id,
                "source_sha256": request.source_sha256,
                "candidate": request.candidate,
                "authorizing_core": draft.core_row[1],
                "complete_family_candidates": tuple(family),
            })
        self.categories = (
            self._attach_boundary_phase_classification_provenance(
                categories,
                batch,
                phase_provenance_by_request_id,
            )
        )
        self.source_accepted_core_boundary_records = records
        self.source_accepted_core_boundary_phase_attestations = tuple(
            phase_attestations
        )
        self._source_accepted_core_boundary_family_records = (
            tuple(getattr(
                self,
                "_source_accepted_core_boundary_family_records",
                (),
            ))
            + tuple(completed_family_records)
        )
        self.source_accepted_core_boundary_stats = {
            "enabled": True,
            **stats.as_dict(),
            "phase_exact": sum(
                record["phase_attestation"]["status"] == "EXACT"
                for record in phase_attestations
            ),
            "phase_ambiguous": sum(
                record["phase_attestation"]["status"] == "AMBIGUOUS"
                for record in phase_attestations
            ),
            "phase_unassessed": sum(
                record["phase_attestation"]["status"] == "UNASSESSED"
                for record in phase_attestations
            ),
        }
        return stats.completed_rows > 0

    def _source_intact_recovery_contexts(self, source_cds):
        """Build exhaustive recovery contexts from source geometry.

        The historical fuzzy-search locus key is retained only as the target
        output key.  It never defines a recovery family: source-replayed,
        overlapping candidate envelopes do.  Consequently two disjoint
        arrays remain independent even when an old detector key or one
        source-superlocus happens to contain both.
        """

        source = self.genome_sequence
        self._source_intact_lineage_only_excluded = 0
        self._source_intact_context_tandem_evaluations = 0
        baseline_rows = getattr(
            self, "_source_architecture_baseline_rows", None
        )
        runtime_lineage = getattr(
            self, "source_architecture_rescue_lineage", ()
        )
        augmentation_stats = getattr(
            self, "_source_architecture_augmentation_stats", None
        )
        completed_records = getattr(
            self, "_source_accepted_core_boundary_family_records", ()
        )
        if not isinstance(baseline_rows, dict):
            raise RuntimeError(
                "source-intact recovery baseline rows must be a mapping"
            )
        if not isinstance(augmentation_stats, dict):
            augmentation_stats = {}
        try:
            runtime_lineage = tuple(runtime_lineage)
            completed_records = tuple(completed_records)
        except TypeError as error:
            raise RuntimeError(
                "source-intact recovery lineage must be iterable"
            ) from error

        candidates = {}
        intervals = {}
        baseline_loci_by_id = {}
        superloci_by_id = {}

        def add_candidate(candidate):
            candidate_id = id(candidate)
            previous = candidates.get(candidate_id)
            if previous is not None and previous is not candidate:
                raise RuntimeError(
                    "source-intact recovery candidate identity collision"
                )
            if candidate_id in candidates:
                return candidate_id
            reconstruction = validate_candidate_reconstruction(
                source, candidate
            )
            if not reconstruction.valid:
                raise RuntimeError(
                    "source-intact recovery family failed source replay"
                )
            repeat_intervals = tuple(
                (start, start + len(repeat))
                for start, repeat in zip(
                    reconstruction.genomic_repeat_starts,
                    candidate.list_repeats,
                )
            )
            if not repeat_intervals:
                raise RuntimeError(
                    "source-intact recovery candidate has no repeats"
                )
            candidates[candidate_id] = candidate
            intervals[candidate_id] = (
                repeat_intervals[0][0], repeat_intervals[-1][1]
            )
            return candidate_id

        baseline_valid = True
        baseline_count = 0
        for locus_key, rows in baseline_rows.items():
            if not isinstance(rows, (list, tuple)) or not rows:
                baseline_valid = False
                continue
            baseline_count += len(rows)
            for row in rows:
                if (
                        not isinstance(row, (list, tuple))
                        or len(row) < 2
                        or row[1] is None):
                    baseline_valid = False
                    continue
                candidate_id = add_candidate(row[1])
                baseline_loci_by_id.setdefault(
                    candidate_id, set()
                ).add(locus_key)

        lineage_valid = True
        lineage_candidate_ids = set()
        lineage_path_ids = set()
        runtime_entries_by_superlocus = {}
        for entry in runtime_lineage:
            if not isinstance(entry, dict):
                lineage_valid = False
                continue
            candidate = entry.get("candidate")
            path_id = entry.get("path_id")
            superlocus_id = entry.get("superlocus_id")
            if (
                    candidate is None
                    or not isinstance(path_id, str)
                    or not path_id
                    or not isinstance(superlocus_id, str)
                    or not superlocus_id
                    or path_id in lineage_path_ids
                    or id(candidate) in lineage_candidate_ids):
                lineage_valid = False
                continue
            candidate_id = add_candidate(candidate)
            lineage_path_ids.add(path_id)
            lineage_candidate_ids.add(candidate_id)
            superloci_by_id.setdefault(candidate_id, set()).add(
                superlocus_id
            )
            runtime_entries_by_superlocus.setdefault(
                superlocus_id, []
            ).append(entry)
            claimed_geometry = entry.get("source_geometry")
            if claimed_geometry is not None:
                observed = tuple(
                    (start, start + len(repeat))
                    for start, repeat in zip(
                        candidate.list_repeat_starts,
                        candidate.list_repeats,
                    )
                )
                try:
                    claimed_geometry = tuple(
                        tuple(item) for item in claimed_geometry
                    )
                except TypeError:
                    lineage_valid = False
                else:
                    if claimed_geometry != observed:
                        lineage_valid = False

        superlocus_records = augmentation_stats.get("superloci", ())
        if not isinstance(superlocus_records, (list, tuple)):
            superlocus_records = ()
            lineage_valid = False
        records_by_superlocus = {}
        for record in superlocus_records:
            if (
                    not isinstance(record, dict)
                    or not isinstance(record.get("superlocus_id"), str)
                    or record["superlocus_id"] in records_by_superlocus):
                lineage_valid = False
                continue
            superlocus_id = record["superlocus_id"]
            records_by_superlocus[superlocus_id] = record
            legacy_loci = record.get("legacy_locus_keys", ())
            if not isinstance(legacy_loci, (list, tuple)):
                lineage_valid = False
                continue
            legacy_loci = set(legacy_loci)
            for candidate_id, candidate_loci in baseline_loci_by_id.items():
                if candidate_loci.intersection(legacy_loci):
                    superloci_by_id.setdefault(
                        candidate_id, set()
                    ).add(superlocus_id)

        completed_valid = True
        for record in completed_records:
            if not isinstance(record, dict):
                completed_valid = False
                continue
            completed = record.get("candidate")
            core = record.get("authorizing_core")
            family = record.get("complete_family_candidates")
            if (
                    completed is None
                    or core is None
                    or not isinstance(family, tuple)
                    or not any(member is completed for member in family)
                    or not any(member is core for member in family)):
                completed_valid = False
                continue
            for member in family:
                add_candidate(member)
            completed_id = id(completed)
            core_superloci = superloci_by_id.get(id(core), set())
            if not core_superloci:
                completed_valid = False
            superloci_by_id.setdefault(completed_id, set()).update(
                core_superloci
            )

        current_rows = []
        current_candidate_ids = set()
        for category_index, category in enumerate(self.categories):
            if not isinstance(category, dict):
                raise RuntimeError(
                    "source-intact recovery categories must be mappings"
                )
            for locus_key, rows in category.items():
                if not isinstance(rows, (list, tuple)) or not rows:
                    raise RuntimeError(
                        "source-intact recovery loci need non-empty rows"
                    )
                for row_ordinal, row in enumerate(rows):
                    if (
                            not isinstance(row, (list, tuple))
                            or len(row) not in (3, 4)
                            or row[1] is None):
                        raise RuntimeError(
                            "source-intact recovery received malformed row"
                        )
                    candidate_id = id(row[1])
                    if candidate_id not in candidates:
                        raise RuntimeError(
                            "category candidate is absent from the complete "
                            "baseline/lineage/completion universe"
                        )
                    if candidate_id in current_candidate_ids:
                        raise RuntimeError(
                            "one source-intact candidate has multiple rows"
                        )
                    current_candidate_ids.add(candidate_id)
                    current_rows.append((
                        category_index,
                        locus_key,
                        row_ordinal,
                        row,
                        candidate_id,
                    ))

        lineage_only_candidate_ids = set()
        for candidate_id, candidate in candidates.items():
            provenance = getattr(candidate, "source_path_provenance", None)
            if provenance is not None:
                _validated_source_path_family_keys(source, candidate)
            explicit_lineage_only = getattr(
                candidate, "source_path_lineage_only", False
            )
            if not isinstance(explicit_lineage_only, bool):
                raise RuntimeError(
                    "source path lineage-only marker must be boolean"
                )
            if explicit_lineage_only or (
                    provenance is not None
                    and any(
                        record.derivation == "local_chain_family_parent"
                        for record in provenance
                    )):
                lineage_only_candidate_ids.add(candidate_id)

        parent = {candidate_id: candidate_id for candidate_id in candidates}

        def find(candidate_id):
            root = candidate_id
            while parent[root] != root:
                root = parent[root]
            while parent[candidate_id] != candidate_id:
                next_id = parent[candidate_id]
                parent[candidate_id] = root
                candidate_id = next_id
            return root

        def join(first_id, second_id):
            first_root = find(first_id)
            second_root = find(second_id)
            if first_root != second_root:
                parent[second_root] = first_root

        active = []
        for candidate_id in sorted(
                candidates,
                key=lambda item: (
                    intervals[item],
                    self._candidate_family_sort_key(candidates[item]),
                )):
            start, end = intervals[candidate_id]
            active = [
                other_id for other_id in active
                if intervals[other_id][1] > start
            ]
            for other_id in active:
                other_start, other_end = intervals[other_id]
                if (
                        candidate_id not in lineage_only_candidate_ids
                        and other_id not in lineage_only_candidate_ids
                        and max(start, other_start) < min(end, other_end)):
                    join(candidate_id, other_id)
            active.append(candidate_id)

        members_by_root = {}
        for candidate_id in candidates:
            members_by_root.setdefault(find(candidate_id), []).append(
                candidates[candidate_id]
            )
        for family in members_by_root.values():
            family.sort(key=self._candidate_family_sort_key)
        family_by_root = {}
        for root, members in members_by_root.items():
            member_ids = {id(member) for member in members}
            member_superloci = set()
            for member_id in member_ids:
                member_superloci.update(superloci_by_id.get(member_id, ()))
            evidence = list(members)
            if not member_ids.issubset(lineage_only_candidate_ids):
                for parent_id in lineage_only_candidate_ids:
                    if parent_id in member_ids:
                        continue
                    if not member_superloci.intersection(
                            superloci_by_id.get(parent_id, ())):
                        continue
                    parent_interval = intervals[parent_id]
                    if any(
                            max(parent_interval[0], intervals[member_id][0])
                            < min(parent_interval[1], intervals[member_id][1])
                            for member_id in member_ids):
                        evidence.append(candidates[parent_id])
            evidence.sort(key=self._candidate_family_sort_key)
            family_by_root[root] = tuple(evidence)

        stats_counts_valid = (
            baseline_valid
            and augmentation_stats.get("source_validated_candidate_count")
            == baseline_count
            and augmentation_stats.get("lineage_candidate_count")
            == len(runtime_lineage)
            and augmentation_stats.get("superlocus_count")
            == len(superlocus_records)
        )
        convergence_counts = tuple(
            augmentation_stats.get(field_name)
            for field_name in (
                "confirmed_converged_superlocus_count",
                "nonconverged_superlocus_count",
                "commit_eligible_superlocus_count",
            )
        )
        convergence_counts_valid = all(
            isinstance(value, int)
            and not isinstance(value, bool)
            and value >= 0
            for value in convergence_counts
        )
        convergence_globally_consistent = (
            stats_counts_valid
            and convergence_counts_valid
            and convergence_counts[0] + convergence_counts[1]
            == len(superlocus_records)
            and convergence_counts[2] == convergence_counts[0]
        )

        recovery_by_root = {}
        accepted_rows = []
        for item in current_rows:
            category_index, _locus, _ordinal, _row, candidate_id = item
            root = find(candidate_id)
            if category_index in (3, 4):
                if candidate_id in lineage_only_candidate_ids:
                    self._source_intact_lineage_only_excluded += 1
                    continue
                recovery_by_root.setdefault(root, []).append(item)
            elif category_index in (0, 2):
                accepted_rows.append(item)

        tandem_by_candidate_id = {}
        complete_family_by_candidate_id = {}

        def exact_source_geometry_key(candidate):
            try:
                consensus = candidate.consensus
            except (AttributeError, TypeError) as error:
                raise RuntimeError(
                    "source-intact family candidate has no consensus"
                ) from error
            if not isinstance(consensus, str) or not consensus:
                raise RuntimeError(
                    "source-intact family candidate consensus is invalid"
                )
            return self._candidate_family_sort_key(candidate) + (consensus,)

        def complete_family_for(candidate_id):
            """Return one source geometry per member, retaining the subject.

            Architecture completion can materialize two distinct candidate
            objects for the same exact source path. Downstream recovery
            attestations require geometry-unique families and require their
            subject to occur exactly once, so collapse only exact geometry
            duplicates here while deterministically preferring the subject.
            """

            family = complete_family_by_candidate_id.get(candidate_id)
            if family is None:
                subject = candidates[candidate_id]
                subject_key = exact_source_geometry_key(subject)
                unique_by_geometry = {subject_key: subject}
                for member in family_by_root[find(candidate_id)]:
                    unique_by_geometry.setdefault(
                        exact_source_geometry_key(member), member
                    )
                family = tuple(
                    unique_by_geometry[key]
                    for key in sorted(unique_by_geometry)
                )
                complete_family_by_candidate_id[candidate_id] = family
            return family

        def tandem_for(candidate_id):
            decision = tandem_by_candidate_id.get(candidate_id)
            if decision is None:
                family = complete_family_for(candidate_id)
                decision = validate_accepted_family_tandem(
                    source,
                    candidates[candidate_id],
                    family,
                    source_cds,
                )
                self._source_intact_context_tandem_evaluations += 1
                tandem_by_candidate_id[candidate_id] = decision
            return decision

        contexts = []
        for root, recovery_items in sorted(
                recovery_by_root.items(),
                key=lambda item: min(
                    intervals[row_item[4]][0]
                    for row_item in item[1]
                )):
            family = family_by_root[root]
            recovery_items = tuple(sorted(
                recovery_items,
                key=lambda item: (
                    intervals[item[4]], item[0], repr(item[1]), item[2]
                ),
            ))
            recovery_intervals = tuple(
                intervals[item[4]] for item in recovery_items
            )
            directly_overlapping_accepted = tuple(
                item for item in accepted_rows
                if any(
                    max(intervals[item[4]][0], recovery_interval[0])
                    < min(intervals[item[4]][1], recovery_interval[1])
                    for recovery_interval in recovery_intervals
                )
            )
            recovery_contexts = tuple(
                SourceIntactRecoveryCandidateContext(
                    scored_row=item[3],
                    target_locus_key=item[1],
                    complete_family_candidates=complete_family_for(item[4]),
                    tandem_decision=tandem_for(item[4]),
                )
                for item in recovery_items
            )
            accepted_contexts = tuple(
                SourceBackedRowFamilyEvidence(
                    scored_row=item[3],
                    complete_family_candidates=complete_family_for(item[4]),
                    tandem_decision=tandem_for(item[4]),
                )
                for item in directly_overlapping_accepted
            )
            relationships = []
            for left_index, left in enumerate(recovery_items):
                for right in recovery_items[left_index + 1:]:
                    left_interval = intervals[left[4]]
                    right_interval = intervals[right[4]]
                    if max(left_interval[0], right_interval[0]) >= min(
                            left_interval[1], right_interval[1]):
                        continue
                    relationships.append(
                        validate_source_intact_alternative_relationship(
                            source, left[3], right[3]
                        )
                    )

            relevant_superloci = set()
            for member in family:
                relevant_superloci.update(
                    superloci_by_id.get(id(member), ())
                )
            relevant_records = tuple(
                records_by_superlocus.get(superlocus_id)
                for superlocus_id in sorted(relevant_superloci)
            )
            local_augmentation_complete = all((
                convergence_globally_consistent,
                bool(relevant_superloci),
                all(record is not None for record in relevant_records),
                all(
                    record is not None
                    and record.get("convergence_status") == "converged"
                    and record.get("architecture_commit_eligible") is True
                    and record.get("bounded_candidate_count") == len(
                        runtime_entries_by_superlocus.get(
                            record["superlocus_id"], ()
                        )
                    )
                    and record.get(
                        "confirmation_semantic_state_changed"
                    ) is False
                    and record.get(
                        "confirmation_terminal_geometries_added"
                    ) == 0
                    and record.get(
                        "confirmation_shared_geometries_added"
                    ) == 0
                    for record in relevant_records
                ),
                all(
                    entry.get("superlocus_convergence_status")
                    == "converged"
                    and entry.get("architecture_commit_eligible") is True
                    for superlocus_id in relevant_superloci
                    for entry in runtime_entries_by_superlocus.get(
                        superlocus_id, ()
                    )
                ),
            ))
            contexts.append(SourceIntactRecoveryLocusContext(
                context_key={
                    "coordinate_system": "zero_based_half_open",
                    "source_geometry_component": [
                        min(interval[0] for interval in recovery_intervals),
                        max(interval[1] for interval in recovery_intervals),
                    ],
                },
                recovery_candidates=recovery_contexts,
                overlapping_accepted=accepted_contexts,
                alternative_relationships=tuple(relationships),
                family_complete=(
                    baseline_valid
                    and completed_valid
                    and current_candidate_ids.issubset(candidates)
                ),
                lineage_complete=(
                    lineage_valid and stats_counts_valid
                ),
                augmentation_complete=local_augmentation_complete,
                accepted_overlap_complete=True,
            ))
        return tuple(contexts)

    def _recover_source_intact_candidates_transactionally(self):
        """Recover only source-attested classifier discards as Possible."""

        if not self._source_architecture_rescue_enabled():
            return False
        print("4f. Review source-intact classifier-discarded arrays")
        source_cds = self._source_cds
        if source_cds is None:
            source_cds = predict_source_cds(self.genome_sequence)
            self._source_cds = source_cds
        contexts = self._source_intact_recovery_contexts(source_cds)
        # Preserve the exact, completeness-attested diagnostic family graph.
        # Later source transactions may replace accepted representatives, so
        # rebuilding this graph from the mutated categories would be neither
        # complete nor identity-equivalent to the reviewed source lineage.
        self._source_intact_recovery_contexts_snapshot = contexts
        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        plan = plan_source_intact_recoveries(
            self.genome_sequence,
            self.categories,
            contexts,
            structural_filter,
            source_cds,
        )
        categories, records, stats = apply_source_intact_recovery_plan(
            self.genome_sequence, self.categories, plan
        )
        self.categories = categories
        self.source_intact_recovery_records = records
        self.source_intact_recovery_stats = {
            "enabled": True,
            **stats.as_dict(),
            "lineage_only_excluded": getattr(
                self, "_source_intact_lineage_only_excluded", 0
            ),
        }
        cost_audit = plan.cost_audit.as_dict()
        context_tandem_count = getattr(
            self, "_source_intact_context_tandem_evaluations", 0
        )
        cost_audit.update({
            "context_builder_direct_tandem_evaluations": (
                context_tandem_count
            ),
            "estimated_total_full_source_hashes": (
                cost_audit[
                    "estimated_residual_full_source_hashes_during_planning"
                ]
                + cost_audit["source_context_full_hashes"]
                + context_tandem_count
            ),
            "pipeline_cost_note": (
                "the context builder computes one direct tandem decision per "
                "reviewed recovery or overlapping accepted candidate; the "
                "planner independently recomputes it to reject stale evidence"
            ),
        })
        self.source_intact_recovery_cost_audit = cost_audit
        return stats.admitted_possible > 0 or stats.replaced_possible > 0

    def _repair_source_periodic_lattices_transactionally(self):
        """Repair accepted periodic occurrence paths under a Possible ceiling."""

        if not self._source_architecture_rescue_enabled():
            return False
        print("4g. Repair source-periodic repeat lattices")
        source_cds_populated = self._source_cds is None
        if source_cds_populated:
            self._source_cds = predict_source_cds(self.genome_sequence)
        source_cds = self._source_cds
        batch = plan_source_periodic_lattice_repairs(
            self.genome_sequence,
            self.categories,
            source_cds,
        )

        fresh_rows = ()
        if batch.score_requests:
            score_key = "source-periodic-lattice-repair-batch"
            scored = self._score_candidates_complete({
                score_key: list(batch.proposed_candidates),
            })
            if set(scored) != {score_key}:
                raise RuntimeError(
                    "periodic-lattice scoring returned an unexpected batch"
                )
            fresh_rows = tuple(scored[score_key])
            if len(fresh_rows) != len(batch.score_requests):
                raise RuntimeError(
                    "periodic-lattice scoring changed request cardinality"
                )

        score_results = []
        for request, row in zip(batch.score_requests, fresh_rows):
            if (
                not isinstance(row, (tuple, list))
                or len(row) < 2
                or row[1] is not request.candidate
            ):
                raise RuntimeError(
                    "periodic-lattice scoring changed candidate order or identity"
                )
            score_results.append(SourcePeriodicLatticeScoreResult(
                request_id=request.request_id,
                scored_row=row,
            ))

        categories, records, stats = apply_source_periodic_lattice_repairs(
            self.genome_sequence,
            self.categories,
            source_cds,
            batch,
            tuple(score_results),
        )
        self.categories = categories
        self.source_periodic_lattice_repair_records = records
        self.source_periodic_lattice_repair_stats = {
            "enabled": True,
            **stats.as_dict(),
        }
        self.source_periodic_lattice_repair_cost_audit = {
            **batch.cost_audit.as_dict(),
            "pipeline_populated_source_cds": source_cds_populated,
            "apply_proposal_replays": stats.apply_proposal_replays,
            "apply_direct_tandem_revalidations": (
                stats.apply_direct_tandem_revalidations
            ),
        }
        return stats.repaired_rows > 0

    def _recover_source_regular_subpaths_transactionally(self):
        """Recover source-regular diagnostic subpaths as Possible only."""

        if not self._source_architecture_rescue_enabled():
            return False
        print("4h. Recover regular subpaths of irregular diagnostics")
        if self._source_cds is None:
            self._source_cds = predict_source_cds(self.genome_sequence)
        source_cds = tuple(self._source_cds)
        source_context = SourceRegularSubpathContext(self.genome_sequence)

        accepted_rows = tuple(
            row
            for category_index in (0, 2)
            for rows in self.categories[category_index].values()
            for row in rows
        )
        accepted_candidates = tuple(row[1] for row in accepted_rows)

        diagnostic_locations = {}
        for category_index in (1, 3, 4):
            for locus_key, rows in self.categories[category_index].items():
                for row_ordinal, row in enumerate(rows):
                    row_id = id(row)
                    if row_id in diagnostic_locations:
                        raise RuntimeError(
                            "one diagnostic row has multiple category locations"
                        )
                    diagnostic_locations[row_id] = (
                        category_index,
                        locus_key,
                        row_ordinal,
                        row,
                    )

        diagnostic_contexts = []
        preflight_parent_count = 0
        preflight_proposal_count = 0
        preflight_tandem_count = 0
        source_intact_contexts = tuple(getattr(
            self, "_source_intact_recovery_contexts_snapshot", ()
        ))
        for locus_ordinal, locus_context in enumerate(
                source_intact_contexts):
            if not all((
                    locus_context.family_complete,
                    locus_context.lineage_complete,
                    locus_context.augmentation_complete,
                    locus_context.accepted_overlap_complete,
            )):
                continue
            for candidate_ordinal, candidate_context in enumerate(
                    locus_context.recovery_candidates):
                location = diagnostic_locations.get(
                    id(candidate_context.scored_row)
                )
                if (
                        location is None
                        or location[3] is not candidate_context.scored_row):
                    # Source-intact recovery may already have admitted or
                    # retired this exact row.  Only rows still preserved as
                    # diagnostics are eligible for subpath review.
                    continue
                (
                    category_index,
                    locus_key,
                    row_ordinal,
                    diagnostic_row,
                ) = location
                preflight_parent_count += 1
                preflight = propose_source_regular_subpaths_in_context(
                    source_context,
                    diagnostic_row[1],
                    parent_identifier=(
                        "pipeline-preflight-{}-{}-{}-{}".format(
                            locus_ordinal,
                            candidate_ordinal,
                            category_index,
                            row_ordinal,
                        )
                    ),
                    parent_disposition=DIAGNOSTIC_PARENT,
                    accepted_candidates=accepted_candidates,
                )
                proposal_evidence = []
                proposal_hashes = set()
                for proposal in preflight.proposals:
                    geometry_sha256 = (
                        proposal.tandem_veto.candidate_geometry_sha256
                    )
                    if geometry_sha256 in proposal_hashes:
                        raise RuntimeError(
                            "regular-subpath preflight emitted duplicate geometry"
                        )
                    proposal_hashes.add(geometry_sha256)
                    decision = validate_accepted_family_tandem(
                        self.genome_sequence,
                        proposal.candidate,
                        candidate_context.complete_family_candidates
                        + (proposal.candidate,),
                        source_cds,
                    )
                    proposal_evidence.append(
                        RegularSubpathProposalTandemEvidence(
                            candidate_geometry_sha256=geometry_sha256,
                            decision=decision,
                        )
                    )
                    preflight_tandem_count += 1
                preflight_proposal_count += len(preflight.proposals)
                diagnostic_contexts.append(
                    SourceRegularSubpathDiagnosticRowContext(
                        context_key={
                            "stage": "source_regular_subpath_recovery",
                            "source_intact_context": locus_context.context_key,
                            "source_intact_locus_ordinal": locus_ordinal,
                            "source_intact_candidate_ordinal": (
                                candidate_ordinal
                            ),
                            "diagnostic_category_index": category_index,
                            "diagnostic_locus_key": locus_key,
                            "diagnostic_row_ordinal": row_ordinal,
                        },
                        diagnostic_row=diagnostic_row,
                        target_locus_key=(
                            candidate_context.target_locus_key
                        ),
                        complete_family_candidates=(
                            candidate_context.complete_family_candidates
                        ),
                        supplied_parent_tandem_decision=(
                            candidate_context.tandem_decision
                        ),
                        supplied_proposal_tandem_evidence=tuple(
                            proposal_evidence
                        ),
                        family_complete=True,
                    )
                )

        batch = plan_source_regular_subpath_recovery_batch(
            source_context,
            self.categories,
            tuple(diagnostic_contexts),
            accepted_rows,
            source_cds,
        )

        score_results = {}
        if batch.score_requests:
            score_key = "source-regular-subpath-recovery-batch"
            scored = self._score_candidates_complete({
                score_key: list(batch.proposed_candidates),
            })
            if set(scored) != {score_key}:
                raise RuntimeError(
                    "regular-subpath scoring returned an unexpected batch"
                )
            fresh_rows = tuple(scored[score_key])
            if len(fresh_rows) != len(batch.score_requests):
                raise RuntimeError(
                    "regular-subpath scoring changed request cardinality"
                )
            for request, row in zip(batch.score_requests, fresh_rows):
                if (
                        not isinstance(row, (tuple, list))
                        or len(row) < 2
                        or row[1] is not request.candidate):
                    raise RuntimeError(
                        "regular-subpath scoring changed candidate order or identity"
                    )
                score_results[request.request_id] = (
                    request,
                    RegularSubpathFreshScoreResult(
                        request_id=request.request_id,
                        scored_row=row,
                        evidence_method=(
                            "pipeline_compact_model_rescore_batch_v2"
                        ),
                    ),
                )

        def fresh_scorer(request):
            try:
                expected_request, result = score_results[request.request_id]
            except KeyError as error:
                raise RuntimeError(
                    "regular-subpath scorer requested an unknown candidate"
                ) from error
            if expected_request is not request:
                raise RuntimeError(
                    "regular-subpath score request identity changed"
                )
            return result

        categories, records, stats = (
            apply_source_regular_subpath_recovery_batch(
                source_context,
                self.categories,
                batch,
                fresh_scorer,
            )
        )
        self.categories = categories
        self.source_regular_subpath_recovery_batch_audit = batch.as_dict()
        self.source_regular_subpath_recovery_records = records
        self.source_regular_subpath_recovery_stats = {
            "enabled": True,
            **stats.as_dict(),
        }
        self.source_regular_subpath_recovery_cost_audit = {
            "source_context_full_hashes": 1,
            "preflight_parents": preflight_parent_count,
            "preflight_structural_proposals": preflight_proposal_count,
            "preflight_direct_tandem_validations": preflight_tandem_count,
            "transaction_parent_tandem_revalidations": (
                stats.parents_reviewed
            ),
            "transaction_proposal_tandem_revalidations": (
                stats.structural_proposals
            ),
            "fresh_scoring_batches": int(bool(batch.score_requests)),
            "fresh_scoring_requests": len(batch.score_requests),
        }
        return stats.possible_admissions > 0

    def _quarantine_likely_tandems_transactionally(self):
        """Soft-filter final accepted rows by source-backed spacer diversity."""

        if self.likely_tandem_stats is not None:
            raise RuntimeError(
                "likely-tandem quarantine may run exactly once per pipeline"
            )
        print("4i. Quarantine likely tandem-repeat representations")
        self.likely_tandem_source_categories = self._copy_categories(
            self.categories
        )
        categories, records, stats = apply_likely_tandem_quarantine(
            self.categories
        )
        self.categories = categories
        self.likely_tandem_records = records
        self.likely_tandem_stats = stats

    def _score_circular_join_candidate(self, candidate):
        """Return one fresh normal-ensemble score for a linearized view."""

        repeat_starts = tuple(candidate.list_repeat_starts)
        repeats = tuple(candidate.list_repeats)
        if not repeat_starts or len(repeat_starts) != len(repeats):
            raise RuntimeError("circular scoring candidate is malformed")
        interval = (0, repeat_starts[-1] + len(repeats[-1]))
        rows_by_interval = self._score_candidates_complete({
            interval: [candidate],
        })
        if set(rows_by_interval) != {interval}:
            raise RuntimeError(
                "circular view scoring returned an unexpected interval"
            )
        rows = rows_by_interval[interval]
        if (
                len(rows) != 1
                or not isinstance(rows[0], (tuple, list))
                or len(rows[0]) < 2
                or rows[0][1] is not candidate):
            raise RuntimeError(
                "circular view scoring changed candidate identity/cardinality"
            )
        return rows[0][0]

    def _evaluate_circular_origin_joins(self):
        """Assess accepted edge children under explicit source topology."""

        print("4j. Assess explicit-topology circular origin joins")
        evidence = self.sequence_topology_evidence
        if evidence is None:
            raise RuntimeError("sequence topology evidence was not resolved")
        circular_evidence = to_circular_topology_evidence(
            evidence,
            evidence.source_id,
            self.genome_sequence,
        )
        self.circular_array_join_result = (
            evaluate_accepted_circular_edge_joins(
                self.genome_sequence,
                self.categories,
                circular_evidence,
                self.parameters,
                self._score_circular_join_candidate,
            )
        )

    @staticmethod
    def _raw_score_band(score):
        score = float(score)
        if score < 0.5:
            return "score_below_0.5"
        if score >= 0.75:
            return "score_at_least_0.75"
        return "score_0.5_to_below_0.75"

    def _annotate_model_score_classification_bases(self):
        """Make every ordinary final classification explicitly auditable."""

        structural_filter = build_advanced_fuzzy_search_filter(
            self.parameters
        )
        annotated = []
        for category_index, category in enumerate(self.categories):
            copied_category = {}
            for locus_key, rows in category.items():
                copied_rows = []
                for row in rows:
                    if (
                            not isinstance(row, (list, tuple))
                            or len(row) not in (3, 4)):
                        raise RuntimeError(
                            "classification provenance received a malformed row"
                        )
                    metadata = dict(row[3]) if len(row) == 4 else {}
                    if "classification_basis" not in metadata:
                        evidence = structural_filter.evidence(row[1])
                        metadata["classification_basis"] = {
                            "method": "model_score_band",
                            "raw_score_band": self._raw_score_band(row[0]),
                            "final_category_before_source_transactions": (
                                CATEGORY_NAMES[category_index]
                            ),
                            "legacy_locus_key": locus_key,
                            "structural_filter": evidence.as_dict(),
                            "coordinate_system": "zero_based_half_open",
                        }
                        copied_row = list(row[:3]) + [metadata]
                    else:
                        copied_row = row
                    copied_rows.append(copied_row)
                copied_category[locus_key] = copied_rows
            annotated.append(copied_category)
        self.categories = annotated

    @staticmethod
    def _attach_retained_tandem_reviews(categories, contexts, records):
        """Attach negative tandem findings without changing their category."""

        contexts = tuple(contexts)
        records = tuple(records)
        if len(contexts) != len(records):
            raise RuntimeError(
                "accepted-family tandem contexts and records diverged"
            )
        review_by_row_id = {}
        for context, record in zip(contexts, records):
            if record.action == "retained":
                review_by_row_id[id(context.accepted_row)] = (
                    context.accepted_row,
                    record.decision.as_dict(),
                )
            elif record.action != "demoted":
                raise RuntimeError("unknown accepted-family tandem action")

        reviewed = []
        matched = set()
        for category in categories:
            copied_category = {}
            for locus_key, rows in category.items():
                copied_rows = []
                for row in rows:
                    review = review_by_row_id.get(id(row))
                    if review is None:
                        copied_rows.append(row)
                        continue
                    expected_row, decision = review
                    if row is not expected_row:
                        raise RuntimeError(
                            "accepted-family retained row identity changed"
                        )
                    if len(row) != 4 or not isinstance(row[3], dict):
                        raise RuntimeError(
                            "retained tandem review needs classification metadata"
                        )
                    metadata = dict(row[3])
                    basis = metadata.get("classification_basis")
                    if not isinstance(basis, dict):
                        raise RuntimeError(
                            "retained tandem review needs a classification basis"
                        )
                    basis = dict(basis)
                    basis["accepted_family_tandem_review"] = decision
                    metadata["classification_basis"] = basis
                    copied_rows.append(list(row[:3]) + [metadata])
                    matched.add(id(row))
                copied_category[locus_key] = copied_rows
            reviewed.append(copied_category)
        if matched != set(review_by_row_id):
            raise RuntimeError(
                "a retained accepted row was absent after tandem review"
            )
        return reviewed

    def _complete_source_backed_terminal_variants_transactionally(self):
        """Commit only completed variants that remain accepted winners."""

        if not (
                self.flags.get("flag_degenerated", True)
                and self.parameters.get(
                    "param_source_backed_terminal_completion", False
                )):
            self._final_scored_rows = self._post_enhancement_scored_rows
            return

        print("4c. Transactionally complete terminal repeats")
        baseline_categories = self.categories
        baseline_rows = self._post_enhancement_scored_rows
        baseline_accepted_count = self._accepted_entry_count(
            baseline_categories
        )
        variant_records, generation_stats = (
            build_source_backed_terminal_variants(
                self.genome_sequence, baseline_categories
            )
        )
        stats = dict(generation_stats)
        stats.update({
            "enabled": True,
            "variants_scored": 0,
            "variants_committed": 0,
            "variants_rejected": 0,
            "committed_repeats": 0,
        })
        if not variant_records:
            self.source_terminal_completion_stats = stats
            self._final_scored_rows = baseline_rows
            return

        candidates = {}
        records_by_interval = {}
        for record in variant_records:
            interval = record["interval"]
            candidates.setdefault(interval, []).append(record["variant"])
            records_by_interval.setdefault(interval, []).append(record)

        fresh_rows = self._score_candidates_complete(candidates)
        if set(fresh_rows) != set(candidates):
            raise RuntimeError(
                "terminal variant scoring returned unexpected intervals"
            )
        for interval, records in records_by_interval.items():
            rows = fresh_rows.get(interval, ())
            if len(rows) != len(records):
                raise RuntimeError(
                    "terminal variant scoring changed candidate cardinality"
                )
            for record, row in zip(records, rows):
                if (
                        not isinstance(row, (tuple, list))
                        or len(row) < 2
                        or row[1] is not record["variant"]):
                    raise RuntimeError(
                        "terminal variant scoring changed order or identity"
                    )
                record["variant_row"] = row

        row_locations = {}
        for interval, rows in baseline_rows.items():
            for row_index, row in enumerate(rows):
                token = id(row)
                if token in row_locations:
                    raise RuntimeError(
                        "baseline score row identity is not unique"
                    )
                row_locations[token] = (interval, row_index, row)

        trial_rows = self._copy_scored_rows(baseline_rows)
        for record in variant_records:
            location = row_locations.get(id(record["base_entry"]))
            if (
                    location is None
                    or location[2] is not record["base_entry"]
                    or location[0] != record["interval"]):
                raise RuntimeError(
                    "terminal variant base row lost score provenance"
                )
            trial_rows[location[0]][location[1]] = record["variant_row"]

        trial_categories = ArrayEvaluation.categorize_scored_rows(
            trial_rows,
            parameters=self.parameters,
            flag_dev_mode=self.flag_dev_mode,
        )
        accepted_variant_ids = self._accepted_candidate_ids(trial_categories)
        committed = [
            record for record in variant_records
            if id(record["variant"]) in accepted_variant_ids
        ]

        final_rows = self._copy_scored_rows(baseline_rows)
        for record in committed:
            interval, row_index, _row = row_locations[id(
                record["base_entry"]
            )]
            final_rows[interval][row_index] = record["variant_row"]
        final_categories = ArrayEvaluation.categorize_scored_rows(
            final_rows,
            parameters=self.parameters,
            flag_dev_mode=self.flag_dev_mode,
        )

        final_accepted_ids = self._accepted_candidate_ids(final_categories)
        if any(
                id(record["variant"]) not in final_accepted_ids
                for record in committed):
            raise RuntimeError(
                "a committed terminal variant lost final interval acceptance"
            )
        if self._accepted_entry_count(final_categories) != baseline_accepted_count:
            raise RuntimeError(
                "transactional terminal completion changed accepted-call count"
            )

        stats.update({
            "variants_scored": len(variant_records),
            "variants_committed": len(committed),
            "variants_rejected": len(variant_records) - len(committed),
            "committed_repeats": sum(
                len(record["safe_proposals"]) for record in committed
            ),
        })
        self.source_terminal_completion_stats = stats
        self._final_scored_rows = final_rows
        self.categories = final_categories

    def _drop_nonaccepted_categories(self):
        retain_requested_categories(
            self.categories,
            self.parameters.get("param_accepted_only", False),
        )

    def _results_enhancement(self):
        print("4. Enhance evaluated arrays")
        a_enh = EvaluatedArraysEnhancement(file_path=self.file_path,
                                           categories=self.categories,
                                           parameters=self.parameters,
                                           flag_dev_mode=self.flag_dev_mode,
                                           flag_degenerated=self.flags.get(
                                               "flag_degenerated", True
                                           ))
        self.categories = a_enh.output()
        self.source_enhancement_lineage = tuple(
            a_enh.source_geometry_lineage
        )

    def _run_non_crispr_computation(self):
        print("5. Complement arrays with additional info")
        nac = NonArrayComputations(
            file_path=self.file_path,
            categories=self.categories,
            flags_non_arrays_computations=self.flags,
        )
        self.non_array_data = nac.output()

    def _validate_refined_candidates(self):
        """Fail before evaluation when refinement produced invalid indexes."""

        for interval, candidates in self.dict_crispr_candidates.items():
            for candidate_index, candidate in enumerate(candidates):
                result = validate_candidate_reconstruction(
                    self.genome_sequence,
                    candidate,
                    coordinate_offset=0,
                )
                if not result.valid:
                    raise CandidateReconstructionError(
                        "post-refinement",
                        "Unclassified",
                        interval,
                        candidate_index,
                        result,
                    )

    def _validate_candidate_categories(self, stage):
        """Validate every scored candidate at a named pipeline boundary."""

        for category_index, category in enumerate(self.categories):
            category_name = (
                CATEGORY_NAMES[category_index]
                if category_index < len(CATEGORY_NAMES)
                else "category {}".format(category_index)
            )
            for interval, entries in category.items():
                for candidate_index, entry in enumerate(entries):
                    if not isinstance(entry, (tuple, list)) or len(entry) < 2:
                        raise RuntimeError(
                            "malformed candidate entry at {}: ".format(stage) +
                            "category={!r}, interval={!r}, "
                            "candidate_index={}".format(
                                category_name, interval, candidate_index
                            )
                        )
                    candidate = entry[1]
                    result = validate_candidate_reconstruction(
                        self.genome_sequence,
                        candidate,
                        coordinate_offset=0,
                    )
                    if not result.valid:
                        raise CandidateReconstructionError(
                            stage,
                            category_name,
                            interval,
                            candidate_index,
                            result,
                        )

    def _validate_and_write_output(self):
        self._validate_candidate_categories("pre-output")
        self._write_output()

    def _write_output(self):
        print("6. Write down the results")
        OutputMaker(file_path=self.file_path,
                    parameters=dict(self.parameters, param_internal_repeat_completion=self._internal_repeat_audit(), param_repeat_boundary_repair=(
                        self._repeat_boundary_repair_audit()
                    )),
                    flags=self.flags,
                    result_path=self.result_folder_path,
                    pickle_result_path=self.pickle_folder_path,
                    json_result_path=self.json_folder_path,
                    categories=self.categories,
                    non_array_data=self.non_array_data,
                    list_features=self.list_features,
                    header=self.header,
                    sequence_topology_evidence=(
                        self.sequence_topology_evidence
                    ),
                    circular_array_join_result=(
                        self.circular_array_join_result
                    ),
                    likely_tandem_records=self.likely_tandem_records,
                    likely_tandem_stats=self.likely_tandem_stats,
                    likely_tandem_source_categories=(
                        self.likely_tandem_source_categories
                    ),
                    candidate_evidence_categories=(
                        getattr(self, "_candidate_evidence_categories", None)
                    ),
                    reconstruction_evidence=(
                        getattr(self, "reconstruction_evidence", None)
                    ),
                    possible_reconstruction_evidence=(
                        getattr(self, "possible_reconstruction_evidence", None)
                    ),
                    contextual_spacer_evidence=(
                        getattr(self, "contextual_spacer_evidence", None)
                    ))
