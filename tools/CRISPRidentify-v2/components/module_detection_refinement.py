from components.components_detection_refinement import SameStartEndFilter
from components.components_detection_refinement import AdvancedFuzzySearchFilter
from components.components_detection_refinement import CrisprCandidate
from components.components_detection_refinement import merged_source_path_provenance
from components.components_detection import RejectedFuzzyCandidateWitness
from components.components_detection import RejectedFuzzyCandidateTopology
from components.repeat_boundary_repair import propose_terminal_repairs
from components.internal_repeat_completion import InternalSearchContext, propose_internal_repairs
import hashlib
from components.candidate_origins import inherit_candidate_origins
from components.localized_spacer_context import capture_pending_context


class _RejectedFuzzyIntervalWitness:
    """One expanded topology run used only for locus connectivity."""

    __slots__ = ("start", "end", "source_path_family_parent")

    def __init__(self, start, end, source_path_family_parent):
        self.start = start
        self.end = end
        self.source_path_family_parent = source_path_family_parent


class _SourceLocusKey:
    """Identity key carrying the exact envelope of one candidate component.

    Separate upstream clusters can legitimately resolve to the same source
    interval, so equality intentionally remains object identity.  The public
    refinement output still aggregates candidates by the ordinary
    ``(start, end)`` tuple.
    """

    def __init__(self, start, end):
        self.start = start
        self.end = end


class DetectionRefinement:
    def __init__(self, dict_fuzzy_crisprs, parameters, flag_dev_mode,
                 source_dna=None, internal_search_context=None, *,
                 consume_input=False, capture_contextual_spacers=False):
        # Only an explicitly transferred map may be drained. The detector's
        # original inventory still supplies later seed and attribution stages.
        if type(consume_input) is not bool:
            raise TypeError("refinement input consumption flag must be boolean")
        self.consume_input = consume_input
        if type(capture_contextual_spacers) is not bool:
            raise TypeError("contextual spacer capture flag must be boolean")
        self.capture_contextual_spacers = capture_contextual_spacers
        self.pending_contextual_spacers = []
        self.dict_fuzzy_crisprs = dict_fuzzy_crisprs
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode
        self.dict_fuzzy_crisprs_refined_st_end = {}
        self.dict_fuzzy_crisprs_fully_refined = {}
        self.accepted_only = bool(parameters.get("param_accepted_only", False))
        self.fallback_keys = set()
        self.compact_witness_empty_keys = set()
        self.source_dna = source_dna
        if source_dna is not None and not isinstance(source_dna, str):
            raise ValueError("source DNA must be a string")
        self.source_sha256 = (
            hashlib.sha256(source_dna.encode("ascii")).hexdigest()
            if source_dna is not None else None
        )
        self.boundary_repair_records = []
        self.boundary_repair_stats = {}
        self.internal_parent_candidates = []
        self.internal_repair_records = []
        self.pending_internal_repairs = []
        self.internal_search_context = internal_search_context or (
            InternalSearchContext(source_dna) if source_dna is not None else None
        )

        if self.internal_search_context is not None and self.internal_search_context.source != source_dna:
            raise ValueError("internal search context differs from refinement source")

        self._filter_out_same_start_end_cases()
        self._split_disjoint_source_loci()
        self._filter_out_non_crispr_cases()
        self._collect_pending_internal_repairs()
        self._reformat_ac_crispr_candidates()

    def _filter_out_same_start_end_cases(self):
        ssef = SameStartEndFilter(self.dict_fuzzy_crisprs)
        self.dict_fuzzy_crisprs_refined_st_end = ssef.output()
        if self.consume_input:
            self.dict_fuzzy_crisprs.clear()

    def _normalization_components(self):
        mapping = self.dict_fuzzy_crisprs_refined_st_end
        if self.consume_input:
            # Preserve global phase and component order while transferring
            # each independently owned list to the normalization loop.
            for key in tuple(mapping):
                yield key, mapping.pop(key)
        else:
            yield from mapping.items()

    @staticmethod
    def _fuzzy_source_interval(candidate):
        """Return a validated zero-based, half-open fuzzy-candidate span."""

        if isinstance(candidate, (
                RejectedFuzzyCandidateWitness,
                _RejectedFuzzyIntervalWitness)):
            return candidate.start, candidate.end

        starts = candidate.list_absolute_start
        repeats = candidate.list_repeats
        if not starts or len(starts) != len(repeats):
            raise ValueError(
                "every fuzzy repeat must have one source start coordinate"
            )

        previous_end = None
        for start, repeat in zip(starts, repeats):
            if (
                    not isinstance(start, int)
                    or isinstance(start, bool)
                    or start < 0
                    or not isinstance(repeat, str)
                    or not repeat):
                raise ValueError("fuzzy candidate has invalid source geometry")
            end = start + len(repeat)
            if previous_end is not None and start < previous_end:
                raise ValueError(
                    "fuzzy candidate repeats overlap or are not monotonic"
                )
            previous_end = end
        return starts[0], previous_end

    @classmethod
    def _source_locus_components(cls, candidates):
        """Partition only candidates whose source envelopes are disjoint.

        Overlap is deliberately transitive: alternate repeat phases and
        partial representations of the same biological locus must continue to
        compete.  A strict source gap, however, proves that two representations
        cannot describe the same bases and therefore must not suppress each
        other during fallback or classifier winner selection.
        """

        indexed = sorted(
            (
                cls._fuzzy_source_interval(candidate),
                candidate_index,
                candidate,
            )
            for candidate_index, candidate in enumerate(candidates)
        )
        if not indexed:
            return []

        components = []
        current = [indexed[0]]
        current_end = indexed[0][0][1]
        for record in indexed[1:]:
            (start, end), _candidate_index, _candidate = record
            if start < current_end:
                current.append(record)
                current_end = max(current_end, end)
            else:
                components.append(current)
                current = [record]
                current_end = end
        components.append(current)

        # Preserve the historical candidate traversal order inside each locus.
        return [
            [record[2] for record in sorted(component, key=lambda item: item[1])]
            for component in components
        ]

    @staticmethod
    def _is_local_chain_family_parent(candidate):
        if isinstance(candidate, (
                RejectedFuzzyCandidateWitness,
                _RejectedFuzzyIntervalWitness)):
            return candidate.source_path_family_parent
        provenance = merged_source_path_provenance(candidate)
        if provenance is None:
            return False
        return any(
            getattr(record, "derivation", None)
            == "local_chain_family_parent"
            for record in provenance
        )

    def _split_disjoint_source_loci(self):
        grouped = {}

        def retain_materialized(key, component):
            materialized = [
                candidate for candidate in component
                if not isinstance(candidate, (
                    RejectedFuzzyCandidateWitness,
                    _RejectedFuzzyIntervalWitness,
                    RejectedFuzzyCandidateTopology,
                ))
            ]
            grouped[key] = materialized
            if component and not materialized:
                self.compact_witness_empty_keys.add(key)

        for cluster_key, candidates in (
                self.dict_fuzzy_crisprs_refined_st_end.items()):
            expanded_candidates = []
            for candidate in candidates:
                if not isinstance(
                        candidate, RejectedFuzzyCandidateTopology):
                    expanded_candidates.append(candidate)
                    continue
                expanded_candidates.extend(
                    _RejectedFuzzyIntervalWitness(start, end, False)
                    for start, end in candidate.local_intervals
                )
                expanded_candidates.extend(
                    _RejectedFuzzyIntervalWitness(start, end, True)
                    for start, end in candidate.family_parent_intervals
                )
            candidates = expanded_candidates
            family_parents = []
            local_paths = []
            for candidate in candidates:
                destination = (
                    family_parents
                    if self._is_local_chain_family_parent(candidate)
                    else local_paths
                )
                destination.append(candidate)
            components = self._source_locus_components(local_paths)
            components.extend(self._source_locus_components(family_parents))
            if len(components) <= 1:
                retain_materialized(cluster_key, candidates)
                continue
            for component in components:
                intervals = [
                    self._fuzzy_source_interval(candidate)
                    for candidate in component
                ]
                locus_key = _SourceLocusKey(
                    min(interval[0] for interval in intervals),
                    max(interval[1] for interval in intervals),
                )
                retain_materialized(locus_key, component)
        self.dict_fuzzy_crisprs_refined_st_end = grouped

    def _filter_out_non_crispr_cases(self):
        self.param_min_avg_repeat_length = self.parameters["param_min_avg_repeat_length"]
        self.param_max_avg_repeat_length = self.parameters["param_max_avg_repeat_length"]
        self.param_min_avg_spacer_length = self.parameters["param_min_avg_spacer_length"]
        self.param_max_avg_spacer_length = self.parameters["param_max_avg_spacer_length"]
        self.param_min_repeats = self.parameters["param_min_repeats"]
        self.param_max_identical_spacers = self.parameters["param_max_identical_spacers"]
        self.param_max_identical_cluster_spacers = self.parameters["param_max_identical_cluster_spacers"]

        self.afsf = AdvancedFuzzySearchFilter(min_column_dominance_repeat=0.6,
                                              max_spacer_length=140, max_column_dominance_spacer=0.8,
                                              max_allowed_consecutive_spacers=self.param_max_identical_cluster_spacers,
                                              max_allowed_same_spacers=self.param_max_identical_spacers,
                                              max_inconsistent_columns=5,
                                              min_avg_repeat_length=self.param_min_avg_repeat_length,
                                              max_avg_repeat_length=self.param_max_avg_repeat_length,
                                              min_avg_spacer_length=self.param_min_avg_spacer_length,
                                              max_avg_spacer_length=self.param_max_avg_spacer_length,
                                              min_repeats=self.param_min_repeats)

        for key, values in self._normalization_components():
            if not values and key in self.compact_witness_empty_keys:
                self.dict_fuzzy_crisprs_fully_refined[key] = []
                continue
            # Structural validation deliberately permits an empty spacer
            # claim.  It can represent a complete spacer-deletion event, but
            # only the later source-reconstruction fence can prove that the
            # flanking repeats meet at exactly one interbase boundary.
            #
            # Normalize every well-formed fuzzy alignment before applying any
            # biological predicate.  ``CrisprCandidate`` removes columns that
            # contain only alignment gaps.  Filtering the raw alignment first
            # can therefore reject a source-valid path that passes after this
            # lossless normalization, then let a lower-error failing path win
            # the fallback competition.
            normalized = []
            for index, value in enumerate(values):
                if self.afsf._candidate_structure_is_valid(value):
                    normalized.append((
                        value.number_errors if self.consume_input else value,
                        self._normalize_fuzzy_candidate(value),
                    ))
                if self.consume_input:
                    # Normalized geometry and validated lineage have been
                    # transferred. Only the original error scalar is needed
                    # for the unchanged stable fallback competition.
                    values[index] = None
                del value
            if self.consume_input:
                values.clear()

            repairs, repair_records, repair_stats = propose_terminal_repairs(
                self.source_dna, [candidate for _, candidate in normalized],
                self.afsf, source_sha256=self.source_sha256,
            )
            self.boundary_repair_records.extend(
                dict(record, locus=(key.start, key.end))
                for record in repair_records
            )
            for name, count in repair_stats.items():
                self.boundary_repair_stats[name] = (
                    self.boundary_repair_stats.get(name, 0) + count
                )
            self.internal_parent_candidates.extend(candidate for _, candidate in normalized)
            self.internal_parent_candidates.extend(repairs)
            if self.capture_contextual_spacers and self.source_dna is not None:
                for _, candidate in normalized:
                    pending = capture_pending_context(
                        self.source_dna, self.source_sha256, (key.start, key.end),
                        candidate, self.afsf)
                    if pending is not None:
                        self.pending_contextual_spacers.append(pending)
            normalized_passes = [
                candidate
                for value, candidate in normalized
                if self.afsf(candidate)
            ] + repairs
            if not normalized_passes:
                sorted_by_num_errors = sorted(
                    normalized,
                    key=(lambda item: item[0]) if self.consume_input
                    else (lambda item: item[0].number_errors),
                )
                if sorted_by_num_errors:
                    candidate_fewer_mismatches = sorted_by_num_errors[0][1]
                    self.dict_fuzzy_crisprs_fully_refined[key] = [candidate_fewer_mismatches]
                    self.fallback_keys.add(key)
            else:
                self.dict_fuzzy_crisprs_fully_refined[key] = normalized_passes

    def _collect_pending_internal_repairs(self):
        if self.internal_search_context is None:
            return
        ordered = sorted(self.internal_parent_candidates, key=lambda candidate: (
            candidate.list_repeat_starts[0],
            candidate.list_repeat_starts[-1] + len(candidate.list_repeats[-1]),
        ))
        components = []
        maximum_end = -1
        for candidate in ordered:
            start = candidate.list_repeat_starts[0]
            end = candidate.list_repeat_starts[-1] + len(candidate.list_repeats[-1])
            if start >= maximum_end:
                components.append([])
            components[-1].append(candidate)
            maximum_end = max(maximum_end, end)
        for component in components:
            children, records = propose_internal_repairs(
                self.internal_search_context, component, self.afsf,
            )
            self.pending_internal_repairs.extend(children)
            self.internal_repair_records.extend(records)

    @staticmethod
    def _normalize_fuzzy_candidate(fuzzy):
        """Return the losslessly normalized candidate used for filtering."""

        candidate = CrisprCandidate(
            fuzzy.list_repeats,
            fuzzy.list_gaped_repeats,
            fuzzy.list_spacers,
            fuzzy.list_absolute_start,
        )
        if hasattr(fuzzy, "source_path_restriction_reasons"):
            candidate.source_path_restriction_reasons = fuzzy.source_path_restriction_reasons
        if hasattr(fuzzy, "source_path_classification_reason"):
            reasons = list(getattr(candidate, "source_path_restriction_reasons", ()))
            reasons.append(fuzzy.source_path_classification_reason)
            candidate.source_path_restriction_reasons = tuple(reasons)
        source_path_provenance = merged_source_path_provenance(fuzzy)
        if source_path_provenance is not None:
            DetectionRefinement._attach_source_path_role(
                candidate, source_path_provenance
            )
        # Explicit restrictions remain restrictive even without fuzzy lineage.
        if hasattr(fuzzy, "source_path_classification_eligible"):
            eligible = fuzzy.source_path_classification_eligible
            if type(eligible) is not bool or not eligible:
                candidate.source_path_classification_eligible = eligible
                candidate.source_path_classification_reason = getattr(
                    fuzzy, "source_path_classification_reason",
                    "inherited_source_path_restriction",
                )
        if getattr(fuzzy, "boundary_repair_records", ()):
            candidate.boundary_repair_records = fuzzy.boundary_repair_records
        if hasattr(fuzzy, "competition_origins"):
            inherit_candidate_origins(candidate, (fuzzy,))
        return candidate

    @staticmethod
    def _attach_source_path_role(candidate, source_path_provenance):
        """Bind validated lineage and its output-eligibility contract."""

        candidate.source_path_provenance = source_path_provenance
        lineage_only = any(
            record.derivation == "local_chain_family_parent"
            for record in source_path_provenance
        )
        prior_eligibility = getattr(candidate, "source_path_classification_eligible", True)
        candidate.source_path_classification_eligible = (
            False if lineage_only else prior_eligibility
        )
        if lineage_only:
            candidate.source_path_classification_reason = (
                "broken_local_edge_lineage_parent"
            )

    def _reformat_ac_crispr_candidates(self):
        self.dict_crispr_candidates = {}
        candidates_by_value = {}
        nonindexed_value_keys = set()
        for key, list_fuzzy in self.dict_fuzzy_crisprs_fully_refined.items():
            new_key = (key.start, key.end)
            list_crispr_candidates = self.dict_crispr_candidates.setdefault(
                new_key, []
            )
            by_value = candidates_by_value.setdefault(new_key, {})
            for fuzzy in list_fuzzy:
                candidate = (
                    fuzzy
                    if isinstance(fuzzy, CrisprCandidate)
                    else self._normalize_fuzzy_candidate(fuzzy)
                )
                if (self.accepted_only and key in self.fallback_keys
                        and not self.afsf(candidate)):
                    continue
                equality_fields = (
                    candidate.list_repeats,
                    candidate.list_repeats_gaped,
                    candidate.list_spacers,
                    candidate.list_repeat_starts,
                )
                canonical_value = (
                    type(candidate) is CrisprCandidate
                    and all(
                        type(field) is list for field in equality_fields
                    )
                    and all(
                        type(value) is str
                        for field in equality_fields[:3]
                        for value in field
                    )
                    and all(
                        type(value) is int
                        for value in candidate.list_repeat_starts
                    )
                )
                if not canonical_value:
                    nonindexed_value_keys.add(new_key)
                use_value_index = (
                    canonical_value
                    and new_key not in nonindexed_value_keys
                )
                if use_value_index:
                    candidate_key = tuple(
                        tuple(field) for field in equality_fields
                    )
                    existing = by_value.get(candidate_key)
                else:
                    # Preserve the historical equality contract for public
                    # callers that supply tuple/list subclasses, custom
                    # equality, or unhashable coordinate values. Production
                    # candidates use the canonical indexed representation.
                    candidate_key = None
                    existing = next((
                        prior
                        for prior in list_crispr_candidates
                        if candidate == prior
                    ), None)
                if existing is None:
                    list_crispr_candidates.append(candidate)
                    if use_value_index:
                        by_value[candidate_key] = candidate
                else:
                    inherit_candidate_origins(existing, (existing, candidate))
                    merged_provenance = merged_source_path_provenance(
                        existing, candidate
                    )
                    if merged_provenance is not None:
                        self._attach_source_path_role(
                            existing, merged_provenance
                        )
                    incoming_eligibility = getattr(
                        candidate, "source_path_classification_eligible", True,
                    )
                    if type(incoming_eligibility) is not bool or not incoming_eligibility:
                        existing.source_path_classification_eligible = incoming_eligibility
                    repair_records = list(getattr(existing, "boundary_repair_records", ()))
                    for record in getattr(candidate, "boundary_repair_records", ()):
                        if record not in repair_records:
                            repair_records.append(record)
                    if repair_records:
                        existing.boundary_repair_records = tuple(repair_records)

    def output(self):
        return self.dict_crispr_candidates
