import multiprocessing
from collections import defaultdict, deque
from multiprocessing import Pool

from components.components_detection import VmatchRun
from components.components_detection import ClusterMaker
from components.components_detection import FilterApproximationClusters
from components.components_detection import StartEndEnhancementClusters
from components.components_detection import IntermediateEnhancementClusters
from components.components_detection import ClusterSequence
from components.components_detection import FuzzySearch
from components.components_detection import RejectedFuzzyCandidateWitness
from components.components_detection import RejectedFuzzyCandidateTopology
from components.components_detection import TwoRepeatTerminalRepairSeed
from components.components_detection_refinement import AdvancedFuzzySearchFilter
from components.components_detection_refinement import merged_source_path_provenance
from components import literal_fuzzy
from components.fuzzy_lattice_storage import LatticeDecodeCache
from components.fuzzy_boundary_lattice import bind_fuzzy_boundary_source
from components.fuzzy_search_utils import aligned_regions_from_spans
from components.fuzzy_search_utils import anchored_subfamilies
from components.fuzzy_search_utils import substring_families
from components.reconstruction_validator import read_single_fasta
from components.fast_run_rescue import MAX_DELTA_PATTERNS
from components.fast_run_rescue import round_robin_start_end_repeat_variants
from components.runtime_resources import configure_native_threads


FAMILY_MEMBER_CHUNK_SIZE = 8
# Keep enough independent cluster batches for the small per-genome worker
# pool, while amortizing source encoding, IPC, and native scratch allocation.
# Results are flattened in the original pattern order, so this changes task
# geometry only—not matcher semantics or candidate order.
NATIVE_PATTERN_BATCH_SIZE = 128
FUZZY_TASK_CHUNK_SIZE = 8
MAX_TWO_REPEAT_TERMINAL_REPAIR_SEEDS_PER_CLUSTER = 128
_FUZZY_SOURCE_DNA = None
_FUZZY_BOUNDARY_SOURCE = None
_FUZZY_BUILD_BOUNDARY_PATHS = None


def _native_worker_pool(*args, **kwargs):
    """Start clean matcher workers without inheriting retained parent arrays."""
    return multiprocessing.get_context("spawn").Pool(*args, **kwargs)


def _ordered_bounded_results(pool, worker, tasks, window):
    """Consume in input order with at most ``window`` outstanding batches.

    Unlike an ordered imap, a slow first task cannot accumulate the rest of
    the campaign in the result cache. The window bounds scheduling only;
    every task and every candidate is still consumed without truncation.
    """
    if window < 1:
        raise ValueError("worker window must be positive")
    tasks = iter(tasks)
    pending = deque()
    for _ in range(window):
        try:
            task = next(tasks)
        except StopIteration:
            break
        pending.append(pool.apply_async(worker, (task,)))
    while pending:
        yield pending.popleft().get()
        try:
            task = next(tasks)
        except StopIteration:
            continue
        pending.append(pool.apply_async(worker, (task,)))


def _accepted_repeat_length_impossible(repeat, parameters):
    """Prove that no <=3 insertion/deletion match can pass repeat length."""

    if parameters is None:
        return False
    minimum = parameters["param_min_avg_repeat_length"]
    maximum = parameters["param_max_avg_repeat_length"]
    return len(repeat) + 3 < minimum or len(repeat) - 3 > maximum


def _initialize_fuzzy_source(source_dna, build_boundary_paths=True):
    """Install one immutable whole-genome coordinate authority per worker."""

    if not isinstance(source_dna, str):
        raise ValueError("fuzzy search source DNA must be a string")
    if not isinstance(build_boundary_paths, bool):
        raise TypeError("fuzzy boundary path worker flag must be boolean")
    global _FUZZY_SOURCE_DNA, _FUZZY_BOUNDARY_SOURCE
    global _FUZZY_BUILD_BOUNDARY_PATHS
    _FUZZY_SOURCE_DNA = source_dna
    _FUZZY_BOUNDARY_SOURCE = bind_fuzzy_boundary_source(source_dna)
    _FUZZY_BUILD_BOUNDARY_PATHS = build_boundary_paths


def _fuzzy_source_dna():
    if _FUZZY_SOURCE_DNA is None:
        raise RuntimeError("fuzzy search worker has no original source DNA")
    return _FUZZY_SOURCE_DNA


def _initialize_fuzzy_worker(source_dna, build_boundary_paths=True):
    """Each detector process uses one native thread within the shared budget."""
    configure_native_threads(1)
    _initialize_fuzzy_source(source_dna, build_boundary_paths)


def _fuzzy_boundary_source():
    if _FUZZY_BOUNDARY_SOURCE is None:
        raise RuntimeError(
            "fuzzy search worker has no bound source identity"
        )
    return _FUZZY_BOUNDARY_SOURCE


def _fuzzy_build_boundary_paths():
    if not isinstance(_FUZZY_BUILD_BOUNDARY_PATHS, bool):
        raise RuntimeError(
            "fuzzy search worker has no boundary-path policy"
        )
    return _FUZZY_BUILD_BOUNDARY_PATHS


class Detection:
    def __init__(self, file_path, flags, parameters, flag_dev_mode):
        self.file_path = file_path
        self.flags = flags
        self.parameters = parameters
        self.flag_parallel = flags["flag_parallel"]
        self.flag_cpu = flags["flag_cpu"]
        self.flag_fast_run = flags["flag_fast_run"]
        self.flag_fast_run_rescue = flags.get("flag_fast_run_rescue", False)
        if not isinstance(self.flag_fast_run_rescue, bool):
            raise ValueError("flag_fast_run_rescue must be boolean")
        self.flag_fast_two_repeat_terminal_repair = flags.get(
            "flag_fast_two_repeat_terminal_repair", False
        )
        if not isinstance(
                self.flag_fast_two_repeat_terminal_repair, bool):
            raise ValueError(
                "flag_fast_two_repeat_terminal_repair must be boolean"
            )
        self.flag_fast_accepted_family_rescue = flags.get(
            "flag_fast_accepted_family_rescue", False
        )
        if not isinstance(self.flag_fast_accepted_family_rescue, bool):
            raise ValueError(
                "flag_fast_accepted_family_rescue must be boolean"
            )
        self.fast_run_rescue_seed_scope = flags.get(
            "fast_run_rescue_seed_scope", "cluster"
        )
        if self.fast_run_rescue_seed_scope not in (
                "lineage", "locus", "family", "cluster"):
            raise ValueError(
                "fast_run_rescue_seed_scope must be lineage, locus, family, "
                "or cluster"
            )
        self.fast_run_rescue_trim_profile = flags.get(
            "fast_run_rescue_trim_profile", "full"
        )
        if self.fast_run_rescue_trim_profile not in (
                "full", "total3", "one_sided"):
            raise ValueError(
                "fast_run_rescue_trim_profile must be full, total3, or "
                "one_sided"
            )
        self.fast_run_seed_profile = flags.get(
            "fast_run_seed_profile", "balanced"
        )
        if self.fast_run_seed_profile not in (
                "balanced", "legacy", "sensitive"):
            raise ValueError(
                "fast_run_seed_profile must be balanced, legacy, or sensitive"
            )
        self.flag_enhancement_min_max = flags["flag_enhancement_min_max"]
        self.flag_enhancement_start_end = flags["flag_enhancement_start_end"]
        self.flag_enhancement_intermediate = flags.get(
            "flag_enhancement_intermediate", True
        )
        self.flag_fuzzy_boundary_paths = flags.get(
            "flag_fuzzy_boundary_paths", True
        )
        if not isinstance(self.flag_fuzzy_boundary_paths, bool):
            raise ValueError("flag_fuzzy_boundary_paths must be boolean")
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode
        self.accepted_only = bool(
            parameters.get("param_accepted_only", False)
        )
        self.compact_fuzzy_filter_parameters = (
            {
                name: parameters[name]
                for name in (
                    "param_min_avg_repeat_length",
                    "param_max_avg_repeat_length",
                    "param_min_avg_spacer_length",
                    "param_max_avg_spacer_length",
                    "param_min_repeats",
                    "param_max_identical_spacers",
                    "param_max_identical_cluster_spacers",
                )
            }
            if self.accepted_only else None
        )
        self.fuzzy_search_mode = parameters.get("param_fuzzy_search_mode", "native")
        self.fuzzy_search_flank = parameters.get("param_fuzzy_search_flank", 100)
        self.fuzzy_family_padding = parameters.get("param_fuzzy_family_padding", 15)
        if self.fuzzy_search_mode not in ("full", "family", "native"):
            raise ValueError(
                "param_fuzzy_search_mode must be 'full', 'family', or 'native'"
            )
        if self.fuzzy_search_flank < 0 or self.fuzzy_family_padding < 0:
            raise ValueError("fuzzy search flank and padding must be non-negative")

        self.clusters = []
        self.cluster_sequences = []
        self.dict_fuzzy_crisprs = {}
        self.compact_fuzzy_topology_telemetry = {
            "raw_rejected_candidates": 0,
            "reduced_rejected_candidates": 0,
            "interval_runs": 0,
        }
        self.two_repeat_terminal_repair_seed_stats = {
            "enabled": bool(
                self.flag_fast_run
                and self.accepted_only
                and (
                    self.flag_fast_two_repeat_terminal_repair
                    or self.flag_fast_accepted_family_rescue
                )
            ),
            "captured": 0,
            "retained": 0,
            "clusters_over_cap": 0,
            "seeds_dropped_by_cluster_cap": 0,
            "seed_cap_per_cluster": (
                MAX_TWO_REPEAT_TERMINAL_REPAIR_SEEDS_PER_CLUSTER
            ),
        }
        self.fast_run_rescue_stats = {
            "enabled": bool(
                self.flag_fast_run and self.flag_fast_run_rescue
            ),
            "selected_clusters": 0,
            "delta_patterns": 0,
            "delta_pattern_cap": MAX_DELTA_PATTERNS,
            "delta_pattern_cap_reached": False,
            "delta_candidates": 0,
        }

        self._get_complete_dna()
        self._run_cluster_detection()
        self._extract_cluster_sequences()
        self._run_array_detection()

    def _get_complete_dna(self):
        header, self.dna = read_single_fasta(self.file_path)
        self.input_header = ">{}\n".format(header)
        self.dna_length = len(self.dna)

    def _run_cluster_detection(self):
        vr = VmatchRun(
            self.file_path,
            self.flag_fast_run,
            getattr(self, "fast_run_seed_profile", "balanced"),
        )
        list_repeats_from_vmatch = vr.output()
        #print("list vmatch repeats", list_repeats_from_vmatch)

        cm = ClusterMaker(list_repeats_from_vmatch, self.dna)
        self.clusters = cm.output()

        if self.flag_fast_run:
            return

        if self.flag_enhancement_min_max:
            fa = FilterApproximationClusters(self.clusters)
            self.clusters = fa.output()

        if self.flag_enhancement_start_end:
            st = StartEndEnhancementClusters(self.clusters)
            self.clusters = st.output()

        if (
                self.flag_enhancement_intermediate
                and (
                    self.flag_enhancement_min_max
                    or self.flag_enhancement_start_end
                )):
            ie = IntermediateEnhancementClusters(self.clusters)
            self.clusters = ie.output()

    def _extract_cluster_sequences(self):
        flank = getattr(self, "fuzzy_search_flank", 100)
        if (
                not isinstance(flank, int)
                or isinstance(flank, bool)
                or flank < 0):
            raise ValueError("fuzzy search flank must be a non-negative integer")
        for cluster in self.clusters:
            if (
                    not isinstance(cluster.begin, int)
                    or isinstance(cluster.begin, bool)
                    or not isinstance(cluster.end, int)
                    or isinstance(cluster.end, bool)
                    or cluster.begin < 0
                    or cluster.end <= cluster.begin
                    or cluster.end > len(self.dna)):
                raise ValueError(
                    "cluster has invalid coordinates in the original sequence"
                )
            seq_start = max(0, cluster.begin - flank)
            seq_end = min(len(self.dna), cluster.end + flank)
            cluster_seq = self.dna[seq_start:seq_end]
            if len(cluster_seq) != seq_end - seq_start:
                raise ValueError("cluster sequence does not reconstruct its span")
            tup_cluster_dif_rep = tuple(cluster.list_clust_dif_rep_seq)
            if not tup_cluster_dif_rep or any(
                    not isinstance(repeat, str) or not repeat
                    for repeat in tup_cluster_dif_rep):
                raise ValueError("cluster repeat seeds must be non-empty strings")

            self.cluster_sequences.append(ClusterSequence(cluster_seq, seq_start, seq_end, tup_cluster_dif_rep))

    @staticmethod
    def _parallel_run_fuzzy_run(input_tuple):
        repeat, sequence, start, weighted_error = input_tuple
        return FuzzySearch(
            sequence,
            start,
            repeat,
            weighted_error,
            source_dna=_fuzzy_source_dna(),
            boundary_source_binding=_fuzzy_boundary_source(),
            _build_boundary_paths=_fuzzy_build_boundary_paths(),
        )

    @staticmethod
    def _parallel_run_fuzzy_native(input_tuple):
        repeat, sequence, start, weighted_error = input_tuple
        if literal_fuzzy.supports(repeat):
            matches = literal_fuzzy.finditer(repeat, sequence)
        else:
            pattern = FuzzySearch.compile_pattern(repeat, weighted_error)
            matches = list(pattern.finditer(sequence))

        return FuzzySearch(
            sequence, start, repeat, weighted_error,
            source_dna=_fuzzy_source_dna(),
            precomputed_matches=matches,
            boundary_source_binding=_fuzzy_boundary_source(),
            _build_boundary_paths=_fuzzy_build_boundary_paths(),
        )

    @staticmethod
    def _parallel_run_fuzzy_native_batch(input_tuple):
        """Search one contiguous repeat batch and materialize its candidates."""

        if len(input_tuple) == 4:
            repeats, sequence, start, weighted_error = input_tuple
            compact_parameters = None
            length_impossible = (False,) * len(repeats)
            compact_topology = False
            capture_two_repeat_terminal_repair = False
        elif len(input_tuple) == 6:
            (
                repeats,
                sequence,
                start,
                weighted_error,
                compact_parameters,
                length_impossible,
            ) = input_tuple
            # Six fields are the pre-quotient public worker contract.  It has
            # no cluster-wide proof that repeat patterns are unique, so keep
            # historical witnesses unless an explicit seventh flag is given.
            compact_topology = False
            capture_two_repeat_terminal_repair = False
        elif len(input_tuple) == 7:
            (
                repeats,
                sequence,
                start,
                weighted_error,
                compact_parameters,
                length_impossible,
                compact_topology,
            ) = input_tuple
            capture_two_repeat_terminal_repair = False
        elif len(input_tuple) == 8:
            (
                repeats,
                sequence,
                start,
                weighted_error,
                compact_parameters,
                length_impossible,
                compact_topology,
                capture_two_repeat_terminal_repair,
            ) = input_tuple
        else:
            raise ValueError("native fuzzy worker tuple has invalid arity")
        if not isinstance(compact_topology, bool):
            raise TypeError("compact fuzzy topology flag must be boolean")
        if not isinstance(capture_two_repeat_terminal_repair, bool):
            raise TypeError(
                "two-repeat terminal repair capture flag must be boolean"
            )
        if len(length_impossible) != len(repeats):
            raise ValueError(
                "every fuzzy repeat needs one length-feasibility decision"
            )
        compact_filter = (
            AdvancedFuzzySearchFilter.from_parameters(compact_parameters)
            if compact_parameters is not None else None
        )
        native_support = tuple(
            literal_fuzzy.supports(repeat) for repeat in repeats
        )
        native_repeats = tuple(
            repeat
            for repeat, supported in zip(repeats, native_support)
            if supported
        )
        native_matches = iter(
            literal_fuzzy.iter_match_batches(native_repeats, sequence)
            if native_repeats else ()
        )

        candidates_by_repeat = []
        for repeat, supported, impossible in zip(
                repeats, native_support, length_impossible):
            if supported:
                matches = next(native_matches)
            else:
                pattern = FuzzySearch.compile_pattern(
                    repeat, weighted_error
                )
                matches = list(pattern.finditer(sequence))

            fuzzy_arguments = {
                "source_dna": _fuzzy_source_dna(),
                "precomputed_matches": matches,
                "boundary_source_binding": _fuzzy_boundary_source(),
                "_build_boundary_paths": _fuzzy_build_boundary_paths(),
            }
            if compact_filter is not None:
                fuzzy_arguments["_defer_alignment"] = True
            fuzzy_search = FuzzySearch(
                sequence,
                start,
                repeat,
                weighted_error,
                **fuzzy_arguments
            )
            candidates_by_repeat.append(
                Detection._array_candidates_from_search(
                    fuzzy_search,
                    compact_filter=compact_filter,
                    compact_acceptance_impossible=impossible,
                    compact_topology=compact_topology,
                    capture_two_repeat_terminal_repair=(
                        capture_two_repeat_terminal_repair
                    ),
                )
            )
            # These raw-match owners are obsolete after candidate conversion.
            del matches, fuzzy_arguments, fuzzy_search

        # This worker publishes a flat, fully enumerated candidate inventory.
        # Keep every lattice, path identity and lineage record, but consume the
        # redundant sibling-object enumeration caches before IPC. The direct
        # FuzzySearch family API remains repeatable for its other callers.
        for candidates in candidates_by_repeat:
            for candidate in candidates:
                detach = getattr(candidate, "detach_materialized_path_alternatives", None)
                if callable(detach):
                    detach()
        return candidates_by_repeat

    @staticmethod
    def _parallel_find_family_anchor(input_tuple):
        (group_index, cluster_index, sequence, sequence_start, anchor,
         weighted_error) = input_tuple

        anchor_pattern = FuzzySearch.compile_pattern(anchor, weighted_error)
        anchor_matches = list(anchor_pattern.finditer(sequence))
        fuzzy_search = FuzzySearch(
            sequence, sequence_start, anchor, weighted_error,
            source_dna=_fuzzy_source_dna(),
            precomputed_matches=anchor_matches,
            boundary_source_binding=_fuzzy_boundary_source(),
            _build_boundary_paths=_fuzzy_build_boundary_paths(),
        )
        occurrence_paths = (
            fuzzy_search.complete_occurrence_path_family()
            if fuzzy_search.match_hit else ()
        )
        anchor_spans = sorted({
            (
                start - sequence_start,
                start - sequence_start + len(repeat),
            )
            for path in occurrence_paths
            for start, repeat in zip(
                path.list_absolute_start, path.list_repeats
            )
        })
        return (group_index, cluster_index, anchor,
                fuzzy_search, anchor_spans)

    @staticmethod
    def _parallel_run_fuzzy_regions(input_tuple):
        (cluster_index, sequence, sequence_start, anchor, anchor_spans,
         members, weighted_error, padding) = input_tuple
        results = []

        for repeat in members:
            regions = aligned_regions_from_spans(
                repeat, anchor, anchor_spans, len(sequence), padding
            )
            repeat_pattern = FuzzySearch.compile_pattern(repeat, weighted_error)
            matches = []
            for region_start, region_end in regions:
                matches.extend(repeat_pattern.finditer(
                    sequence, region_start, region_end
                ))

            fuzzy_search = FuzzySearch(
                sequence, sequence_start, repeat, weighted_error,
                source_dna=_fuzzy_source_dna(),
                precomputed_matches=matches,
                boundary_source_binding=_fuzzy_boundary_source(),
                _build_boundary_paths=_fuzzy_build_boundary_paths(),
            )
            results.append((repeat, fuzzy_search))

        return cluster_index, results

    @staticmethod
    def _is_array_candidate(fuzzy_search):
        return (
            fuzzy_search is not None
            and fuzzy_search.match_hit
            and len(fuzzy_search.list_repeats) > 1
        )

    @staticmethod
    def _array_candidates_from_search(
            fuzzy_search, compact_filter=None,
            compact_acceptance_impossible=False,
            compact_topology=False,
            capture_two_repeat_terminal_repair=False):
        """Return every occurrence path and its maximal local chains.

        The ``getattr`` fallback keeps lightweight test doubles and external
        callers that implement the historical FuzzySearch interface working.
        """

        if fuzzy_search is None or not fuzzy_search.match_hit:
            return []
        complete_occurrence_path_family = getattr(
            fuzzy_search, "complete_occurrence_path_family", None
        )
        searches = (
            complete_occurrence_path_family()
            if complete_occurrence_path_family is not None
            else [fuzzy_search]
        )
        candidates = []
        for search in searches:
            complete_local_chain_family = getattr(
                search, "complete_local_chain_family", None
            )
            split_into_local_chains = getattr(
                search, "split_into_local_chains", None
            )
            candidates.extend(
                complete_local_chain_family()
                if complete_local_chain_family is not None
                else (
                    split_into_local_chains()
                    if split_into_local_chains is not None
                    else [search]
                )
            )
        candidates = [
            candidate for candidate in candidates
            if Detection._is_array_candidate(candidate)
        ]
        unique = []
        by_geometry = {}
        for candidate in candidates:
            if not all(hasattr(candidate, attribute) for attribute in (
                    "list_absolute_start", "list_spacers",
                    "list_fuzzy_counts")):
                unique.append(candidate)
                continue
            geometry = (
                candidate.repeat_candidate,
                tuple(candidate.list_absolute_start),
                tuple(candidate.list_repeats),
                tuple(candidate.list_spacers),
                tuple(tuple(counts)
                      for counts in candidate.list_fuzzy_counts),
            )
            existing = by_geometry.get(geometry)
            if existing is None:
                by_geometry[geometry] = candidate
                unique.append(candidate)
                continue

            existing_records = getattr(
                existing, "source_path_provenance", ()
            )
            candidate_records = getattr(
                candidate, "source_path_provenance", ()
            )
            if existing_records or candidate_records:
                merged_records = list(existing_records)
                merged_records.extend(
                    record for record in sorted(set(candidate_records))
                    if record not in merged_records
                )
                existing.source_path_provenance = tuple(merged_records)
                validator = getattr(
                    existing, "_validate_source_path_provenance", None
                )
                if validator is not None:
                    validator()
        if compact_filter is None:
            return unique

        if compact_topology:
            return Detection._compact_topology_candidates(
                unique,
                compact_filter,
                compact_acceptance_impossible,
                capture_two_repeat_terminal_repair,
            )

        compact = []
        for candidate in unique:
            candidate.detach_materialized_path_alternatives()
            lineage_parent = any(
                record.derivation == "local_chain_family_parent"
                for record in candidate.source_path_provenance
            )
            if (
                    lineage_parent
                    or compact_acceptance_impossible
                    or not compact_filter._passes_non_gapped_predicates(
                        candidate
                    )):
                seed = Detection._two_repeat_terminal_repair_seed(
                    candidate,
                    compact_filter,
                    compact_acceptance_impossible,
                    lineage_parent,
                    capture_two_repeat_terminal_repair,
                )
                compact.append(candidate.rejected_candidate_witness(seed))
            else:
                compact.append(candidate._materialize_alignment())
        return compact

    @staticmethod
    def _compact_topology_candidates(
            candidates, compact_filter, compact_acceptance_impossible,
            capture_two_repeat_terminal_repair=False):
        """Resolve exact local classes before discarding rejected key data."""

        best_by_class = {}
        raw_rejected = 0
        for candidate in candidates:
            candidate.detach_materialized_path_alternatives()
            candidate._validate_source_path_provenance()
            class_key = candidate.compact_same_start_end_token_key()
            family_parent = any(
                record.derivation == "local_chain_family_parent"
                for record in candidate.source_path_provenance
            )
            rejected = (
                family_parent
                or compact_acceptance_impossible
                or not compact_filter._passes_non_gapped_predicates(
                    candidate
                )
            )
            raw_rejected += int(rejected)
            current = best_by_class.get(class_key)
            if current is None:
                best_by_class[class_key] = {
                    "candidate": candidate,
                    "rejected": rejected,
                    "family_parent": family_parent,
                }
                continue
            if current["rejected"] != rejected:
                raise ValueError(
                    "compact structural decision changed inside one "
                    "SameStartEnd equivalence class"
                )
            if rejected:
                current["family_parent"] = (
                    current["family_parent"] or family_parent
                )
                if (
                        candidate.number_errors
                        < current["candidate"].number_errors):
                    current["candidate"] = candidate
                continue

            merged_provenance = merged_source_path_provenance(
                current["candidate"], candidate
            )
            if (
                    candidate.number_errors
                    < current["candidate"].number_errors):
                current["candidate"] = candidate
            if merged_provenance is not None:
                current["candidate"].source_path_provenance = (
                    merged_provenance
                )

        survivors = []
        local_intervals = []
        family_parent_intervals = []
        two_repeat_terminal_repair_seeds = []
        reduced_rejected = 0
        for record in best_by_class.values():
            candidate = record["candidate"]
            if not record["rejected"]:
                survivors.append(candidate._materialize_alignment())
                continue
            reduced_rejected += 1
            interval = (
                candidate.list_absolute_start[0],
                candidate.list_absolute_start[-1]
                + len(candidate.list_repeats[-1]),
            )
            destination = (
                family_parent_intervals
                if record["family_parent"] else local_intervals
            )
            destination.append(interval)
            seed = Detection._two_repeat_terminal_repair_seed(
                candidate,
                compact_filter,
                compact_acceptance_impossible,
                record["family_parent"],
                capture_two_repeat_terminal_repair,
            )
            if seed is not None:
                two_repeat_terminal_repair_seeds.append(seed)

        if reduced_rejected:
            survivors.append(RejectedFuzzyCandidateTopology(
                tuple(local_intervals),
                tuple(family_parent_intervals),
                raw_rejected,
                reduced_rejected,
                tuple(two_repeat_terminal_repair_seeds),
            ))
        return survivors

    @staticmethod
    def _two_repeat_terminal_repair_seed(
            candidate, compact_filter, compact_acceptance_impossible,
            family_parent, enabled):
        """Keep only a narrow, sole-repeat-count-failure source geometry."""

        if (
                not enabled
                or family_parent
                or compact_acceptance_impossible
                or len(candidate.list_repeats) != 2
                or len(candidate.list_spacers) != 1):
            return None
        repeat_lengths = tuple(map(len, candidate.list_repeats))
        spacer_length = len(candidate.list_spacers[0])
        if (
                min(repeat_lengths) < 21
                or max(repeat_lengths) > 40
                or not 24 <= spacer_length <= 50):
            return None

        materialized = candidate._materialize_alignment()
        evidence = compact_filter.evidence(materialized)
        failed = tuple(
            name for name, status in evidence.predicate_results
            if status == "failed"
        )
        if (
                not evidence.structural_valid
                or failed != ("minimum_repeat_count",)):
            return None
        starts = tuple(materialized.list_absolute_start)
        ends = tuple(
            start + len(repeat)
            for start, repeat in zip(starts, materialized.list_repeats)
        )
        return TwoRepeatTerminalRepairSeed(starts, ends)

    def _run_full_array_detection(self, weighted_error, worker=None):
        if worker is None:
            worker = Detection._parallel_run_fuzzy_run

        if self.flag_parallel:
            cluster_task_counts = [
                len(cluster_sequence.tuple_repeats)
                for cluster_sequence in self.cluster_sequences
            ]
            total_tasks = sum(cluster_task_counts)
            if total_tasks == 0:
                for cluster_sequence in self.cluster_sequences:
                    self.dict_fuzzy_crisprs[cluster_sequence] = []
                return

            max_possible = multiprocessing.cpu_count()
            num_workers_suggested = (
                max_possible if self.flag_cpu == "ALL" else int(self.flag_cpu)
            )
            num_workers = max(
                1, min(num_workers_suggested, max_possible, total_tasks)
            )
            chunksize = min(FUZZY_TASK_CHUNK_SIZE, total_tasks)
            input_tuples = [
                (repeat, cluster_sequence.sequence,
                 cluster_sequence.start, weighted_error)
                for cluster_sequence in self.cluster_sequences
                for repeat in cluster_sequence.tuple_repeats
            ]

            # One ordered map removes a synchronization barrier per cluster.
            # Chunk pickles also memoize the shared sequence object, avoiding
            # one full sequence serialization for every repeat candidate.
            with Pool(
                    num_workers,
                    initializer=_initialize_fuzzy_worker,
                    initargs=(
                        self.dna, self._fuzzy_boundary_paths_enabled()
                    )) as pool:
                fuzzy_results = pool.map(
                    worker, input_tuples, chunksize=chunksize
                )

            result_offset = 0
            for cluster_sequence, task_count in zip(
                    self.cluster_sequences, cluster_task_counts):
                cluster_results = fuzzy_results[
                    result_offset:result_offset + task_count
                ]
                result_offset += task_count
                self.dict_fuzzy_crisprs[cluster_sequence] = [
                    candidate
                    for result in cluster_results
                    for candidate in self._array_candidates_from_search(result)
                ]
        else:
            for cluster_sequence in self.cluster_sequences:
                fuzzy_results = []
                for repeat in cluster_sequence.tuple_repeats:
                    result = worker((
                        repeat, cluster_sequence.sequence,
                        cluster_sequence.start, weighted_error
                    ))
                    fuzzy_results.extend(
                        self._array_candidates_from_search(result)
                    )
                self.dict_fuzzy_crisprs[cluster_sequence] = fuzzy_results

    def _run_native_array_detection(self, weighted_error):
        cluster_batch_counts = []
        input_tuples = []
        compact_parameters = getattr(
            self, "compact_fuzzy_filter_parameters", None
        )

        for cluster_sequence in self.cluster_sequences:
            batch_count = 0
            compact_topology = (
                compact_parameters is not None
                and len(set(cluster_sequence.tuple_repeats))
                == len(cluster_sequence.tuple_repeats)
            )
            for offset in range(
                    0,
                    len(cluster_sequence.tuple_repeats),
                    NATIVE_PATTERN_BATCH_SIZE):
                repeats = tuple(cluster_sequence.tuple_repeats[
                    offset:offset + NATIVE_PATTERN_BATCH_SIZE
                ])
                length_impossible = tuple(
                    _accepted_repeat_length_impossible(
                        repeat, compact_parameters
                    )
                    for repeat in repeats
                )
                input_tuple = (
                    repeats,
                    cluster_sequence.sequence,
                    cluster_sequence.start,
                    weighted_error,
                )
                if compact_parameters is not None:
                    input_tuple += (
                        compact_parameters,
                        length_impossible,
                        compact_topology,
                    )
                    if (
                            self.flag_fast_two_repeat_terminal_repair
                            or self.flag_fast_accepted_family_rescue):
                        input_tuple += (True,)
                input_tuples.append(input_tuple)
                batch_count += 1
            cluster_batch_counts.append(batch_count)

        total_tasks = len(input_tuples)
        if total_tasks == 0:
            for cluster_sequence in self.cluster_sequences:
                self.dict_fuzzy_crisprs[cluster_sequence] = []
            return

        telemetry = {
            "raw_rejected_candidates": 0,
            "reduced_rejected_candidates": 0,
            "interval_runs": 0,
        }

        def collect(batch_results):
            batches = iter(batch_results)
            lattice_cache = LatticeDecodeCache()
            for cluster_sequence, batch_count in zip(
                    self.cluster_sequences, cluster_batch_counts):
                candidates = []
                for _ in range(batch_count):
                    batch = next(batches)
                    lattice_memo = {}
                    for repeat_candidates in batch:
                        for candidate in repeat_candidates:
                            if isinstance(candidate, FuzzySearch):
                                candidate._compact_lattice_storage(lattice_cache, lattice_memo)
                            candidates.append(candidate)
                            if isinstance(candidate, RejectedFuzzyCandidateTopology):
                                for name in telemetry:
                                    telemetry[name] += getattr(candidate, name)
                    del batch, lattice_memo
                self.dict_fuzzy_crisprs[cluster_sequence] = candidates

        if self.flag_parallel:
            max_possible = multiprocessing.cpu_count()
            num_workers_suggested = (
                max_possible if self.flag_cpu == "ALL" else int(self.flag_cpu)
            )
            num_workers = max(
                1, min(num_workers_suggested, max_possible, total_tasks)
            )
            with _native_worker_pool(
                    num_workers,
                    initializer=_initialize_fuzzy_worker,
                    initargs=(
                        self.dna, self._fuzzy_boundary_paths_enabled()
                    )) as pool:
                collect(_ordered_bounded_results(
                    pool,
                    Detection._parallel_run_fuzzy_native_batch,
                    input_tuples,
                    window=num_workers,
                ))
        else:
            collect(
                Detection._parallel_run_fuzzy_native_batch(input_tuple)
                for input_tuple in input_tuples
            )
        self.compact_fuzzy_topology_telemetry = telemetry

    def _family_tasks(self, weighted_error):
        padding = getattr(self, "fuzzy_family_padding", 15)
        tasks = []
        for cluster_index, cluster_sequence in enumerate(self.cluster_sequences):
            families = substring_families(cluster_sequence.tuple_repeats)
            for family in families:
                for anchor, members in anchored_subfamilies(family):
                    tasks.append((
                        cluster_index, cluster_sequence.sequence,
                        cluster_sequence.start, anchor, tuple(members),
                        weighted_error, padding
                    ))
        return tasks

    @staticmethod
    def _build_regional_tasks(tasks, anchor_results):
        anchor_spans_by_group = {
            group_index: anchor_spans
            for (group_index, _cluster_index, _anchor,
                 _fuzzy_search, anchor_spans) in anchor_results
        }
        regional_tasks = []
        for group_index, task in enumerate(tasks):
            (cluster_index, sequence, sequence_start, anchor, members,
             weighted_error, padding) = task
            derived_members = [repeat for repeat in members if repeat != anchor]
            for offset in range(0, len(derived_members),
                                FAMILY_MEMBER_CHUNK_SIZE):
                member_chunk = tuple(derived_members[
                    offset:offset + FAMILY_MEMBER_CHUNK_SIZE
                ])
                regional_tasks.append((
                    cluster_index, sequence, sequence_start, anchor,
                    anchor_spans_by_group[group_index], member_chunk,
                    weighted_error, padding
                ))
        return regional_tasks

    @staticmethod
    def _family_worker_capacity(tasks):
        regional_chunks = 0
        for task in tasks:
            members = task[4]
            derived_count = max(0, len(members) - 1)
            regional_chunks += (
                derived_count + FAMILY_MEMBER_CHUNK_SIZE - 1
            ) // FAMILY_MEMBER_CHUNK_SIZE
        return max(len(tasks), regional_chunks, 1)

    def _run_family_array_detection(self, weighted_error):
        tasks = self._family_tasks(weighted_error)
        if not tasks:
            for cluster_sequence in self.cluster_sequences:
                self.dict_fuzzy_crisprs[cluster_sequence] = []
            return

        anchor_tasks = [
            (group_index, task[0], task[1], task[2], task[3],
             weighted_error)
            for group_index, task in enumerate(tasks)
        ]

        if self.flag_parallel:
            max_possible = multiprocessing.cpu_count()
            num_workers_suggested = (
                max_possible if self.flag_cpu == "ALL" else int(self.flag_cpu)
            )
            num_workers = max(
                1, min(num_workers_suggested, max_possible,
                       self._family_worker_capacity(tasks))
            )
            with Pool(
                    num_workers,
                    initializer=_initialize_fuzzy_worker,
                    initargs=(
                        self.dna, self._fuzzy_boundary_paths_enabled()
                    )) as pool:
                anchor_results = pool.map(
                    Detection._parallel_find_family_anchor,
                    anchor_tasks, chunksize=1
                )
                regional_tasks = self._build_regional_tasks(
                    tasks, anchor_results
                )
                regional_results = pool.map(
                    Detection._parallel_run_fuzzy_regions,
                    regional_tasks, chunksize=1
                ) if regional_tasks else []
        else:
            anchor_results = [
                Detection._parallel_find_family_anchor(task)
                for task in anchor_tasks
            ]
            regional_tasks = self._build_regional_tasks(tasks, anchor_results)
            regional_results = [
                Detection._parallel_run_fuzzy_regions(task)
                for task in regional_tasks
            ]

        results_by_cluster = [defaultdict(deque)
                              for _ in self.cluster_sequences]
        for (_group_index, cluster_index, anchor,
             fuzzy_search, _anchor_spans) in anchor_results:
            results_by_cluster[cluster_index][anchor].append(fuzzy_search)
        for cluster_index, repeat_results in regional_results:
            for repeat, fuzzy_search in repeat_results:
                results_by_cluster[cluster_index][repeat].append(fuzzy_search)

        # Family processing is deliberately longest-first. Restore the exact
        # candidate order used by the full-domain implementation before the
        # refinement stages see these objects.
        for cluster_index, cluster_sequence in enumerate(self.cluster_sequences):
            ordered_results = []
            pending = results_by_cluster[cluster_index]
            for repeat in cluster_sequence.tuple_repeats:
                if not pending[repeat]:
                    continue
                fuzzy_search = pending[repeat].popleft()
                ordered_results.extend(
                    self._array_candidates_from_search(fuzzy_search)
                )
            self.dict_fuzzy_crisprs[cluster_sequence] = ordered_results

    def _run_array_detection(self):
        weighted_error = "{i<=3,d<=3,s<=3,i+d+s<=6}"
        _initialize_fuzzy_source(
            self.dna, self._fuzzy_boundary_paths_enabled()
        )
        mode = getattr(self, "fuzzy_search_mode", "native")
        if mode == "family":
            self._run_family_array_detection(weighted_error)
        elif mode == "native":
            self._run_native_array_detection(weighted_error)
        elif mode == "full":
            self._run_full_array_detection(weighted_error)
        else:
            raise ValueError(
                "fuzzy search mode must be 'full', 'family', or 'native'"
            )

    def _fuzzy_boundary_paths_enabled(self):
        enabled = getattr(self, "flag_fuzzy_boundary_paths", True)
        if not isinstance(enabled, bool):
            raise ValueError("flag_fuzzy_boundary_paths must be boolean")
        return enabled

    @staticmethod
    def _raw_candidate_source_intervals(candidate):
        """Return source envelopes represented by one raw detection value."""

        if isinstance(candidate, RejectedFuzzyCandidateTopology):
            return tuple(candidate.local_intervals) + tuple(
                candidate.family_parent_intervals
            )
        starts = getattr(candidate, "list_absolute_start", None)
        repeats = getattr(candidate, "list_repeats", None)
        if starts is not None or repeats is not None:
            if (
                    not isinstance(starts, (tuple, list))
                    or not isinstance(repeats, (tuple, list))
                    or not starts
                    or len(starts) != len(repeats)):
                raise ValueError("raw fuzzy candidate has invalid geometry")
            return ((starts[0], starts[-1] + len(repeats[-1])),)
        start = getattr(candidate, "start", None)
        end = getattr(candidate, "end", None)
        if (
                isinstance(start, int)
                and not isinstance(start, bool)
                and isinstance(end, int)
                and not isinstance(end, bool)
                and 0 <= start < end):
            return ((start, end),)
        raise ValueError("raw detection value has no source interval")

    @staticmethod
    def _strictly_overlaps(first, second):
        return first[0] < second[1] and second[0] < first[1]

    def cluster_indices_for_intervals(self, source_intervals):
        """Attribute routed source components to their raw seed clusters."""

        intervals = tuple(source_intervals)
        for interval in intervals:
            if (
                    not isinstance(interval, (tuple, list))
                    or len(interval) != 2
                    or not all(
                        isinstance(value, int) and not isinstance(value, bool)
                        for value in interval
                    )
                    or not 0 <= interval[0] < interval[1] <= len(self.dna)):
                raise ValueError("routed source interval is invalid")

        raw_intervals = []
        for cluster_index, cluster_sequence in enumerate(
                self.cluster_sequences):
            represented = []
            for candidate in self.dict_fuzzy_crisprs.get(
                    cluster_sequence, ()):
                represented.extend(
                    self._raw_candidate_source_intervals(candidate)
                )
            raw_intervals.append(tuple(represented))

        selected = []
        seen = set()
        exact_attributions = 0
        fallback_attributions = 0
        unmatched_intervals = 0
        for interval in intervals:
            matches = [
                cluster_index
                for cluster_index, represented in enumerate(raw_intervals)
                if any(
                    self._strictly_overlaps(interval, raw_interval)
                    for raw_interval in represented
                )
            ]
            if matches:
                exact_attributions += 1
            else:
                matches = [
                    cluster_index
                    for cluster_index, cluster in enumerate(self.clusters)
                    if self._strictly_overlaps(
                        interval, (cluster.begin, cluster.end)
                    )
                ]
                if matches:
                    fallback_attributions += 1
                else:
                    unmatched_intervals += 1
            for cluster_index in matches:
                if cluster_index not in seen:
                    seen.add(cluster_index)
                    selected.append(cluster_index)
        return tuple(selected), {
            "exact_interval_attributions": exact_attributions,
            "fallback_interval_attributions": fallback_attributions,
            "unmatched_intervals": unmatched_intervals,
        }

    def run_start_end_delta(
            self, selected_cluster_indices, repeat_seeds_by_cluster=None):
        """Search only new 4x4-trim seeds and return a merged raw mapping.

        The initial clusters, result mapping, and compact-topology telemetry
        remain immutable. Delta results are copied under the corresponding
        original ``ClusterSequence`` key so ordinary refinement performs the
        same raw/delta competition as a monolithic run.
        """

        try:
            selected = tuple(selected_cluster_indices)
        except TypeError as error:
            raise ValueError("selected cluster IDs must be iterable") from error
        if any(
                not isinstance(index, int)
                or isinstance(index, bool)
                or index < 0
                or index >= len(self.cluster_sequences)
                for index in selected):
            raise ValueError("selected cluster ID is out of range")
        if len(set(selected)) != len(selected):
            raise ValueError("selected cluster IDs must be unique")
        if repeat_seeds_by_cluster is not None:
            if not isinstance(repeat_seeds_by_cluster, dict):
                raise ValueError("repeat seed scope must be a mapping")
            unexpected = set(repeat_seeds_by_cluster).difference(selected)
            if unexpected:
                raise ValueError(
                    "repeat seed scope contains an unselected cluster"
                )

        merged = {
            cluster_sequence: list(values)
            for cluster_sequence, values in self.dict_fuzzy_crisprs.items()
        }
        selected_records = []
        raw_seed_patterns = 0
        selected_seed_patterns = 0
        for cluster_index in selected:
            raw_cluster = self.cluster_sequences[cluster_index]
            raw_seed_patterns += len(raw_cluster.tuple_repeats)
            scoped_repeats = (
                tuple(repeat_seeds_by_cluster.get(
                    cluster_index, raw_cluster.tuple_repeats
                ))
                if repeat_seeds_by_cluster is not None
                else raw_cluster.tuple_repeats
            )
            if (
                    not scoped_repeats
                    or any(seed not in raw_cluster.tuple_repeats
                           for seed in scoped_repeats)):
                raise ValueError(
                    "repeat seed scope must be a non-empty cluster subset"
                )
            selected_seed_patterns += len(scoped_repeats)
            selected_records.append((
                cluster_index, raw_cluster, scoped_repeats
            ))

        allocated_variants, cap_reached = (
            round_robin_start_end_repeat_variants(
                tuple(record[2] for record in selected_records),
                maximum=MAX_DELTA_PATTERNS,
                trim_profile=self.fast_run_rescue_trim_profile,
            )
        )
        delta_records = []
        delta_pattern_counts = []
        for record, variants in zip(selected_records, allocated_variants):
            cluster_index, raw_cluster, _scoped_repeats = record
            delta_pattern_counts.append({
                "cluster_index": cluster_index,
                "patterns": len(variants),
            })
            if not variants:
                continue
            delta_cluster = ClusterSequence(
                raw_cluster.sequence,
                raw_cluster.start,
                raw_cluster.end,
                variants,
            )
            delta_records.append((raw_cluster, delta_cluster, variants))

        if not delta_records:
            stats = {
                "enabled": True,
                "selected_clusters": len(selected),
                "searched_clusters": 0,
                "delta_patterns": 0,
                "delta_pattern_cap": MAX_DELTA_PATTERNS,
                "delta_pattern_cap_reached": cap_reached,
                "delta_pattern_counts": delta_pattern_counts,
                "trim_profile": self.fast_run_rescue_trim_profile,
                "delta_candidates": 0,
                "raw_seed_patterns": raw_seed_patterns,
                "selected_seed_patterns": selected_seed_patterns,
                "delta_topology": {
                    "raw_rejected_candidates": 0,
                    "reduced_rejected_candidates": 0,
                    "interval_runs": 0,
                },
            }
            self.fast_run_rescue_stats = stats
            return merged, stats

        raw_cluster_sequences = self.cluster_sequences
        raw_mapping = self.dict_fuzzy_crisprs
        raw_telemetry = self.compact_fuzzy_topology_telemetry
        try:
            self.cluster_sequences = [
                record[1] for record in delta_records
            ]
            self.dict_fuzzy_crisprs = {}
            self.compact_fuzzy_topology_telemetry = {
                "raw_rejected_candidates": 0,
                "reduced_rejected_candidates": 0,
                "interval_runs": 0,
            }
            self._run_array_detection()
            delta_mapping = self.dict_fuzzy_crisprs
            delta_telemetry = dict(self.compact_fuzzy_topology_telemetry)
        finally:
            self.cluster_sequences = raw_cluster_sequences
            self.dict_fuzzy_crisprs = raw_mapping
            self.compact_fuzzy_topology_telemetry = raw_telemetry

        delta_candidates = 0
        for raw_cluster, delta_cluster, _variants in delta_records:
            candidates = list(delta_mapping.get(delta_cluster, ()))
            merged.setdefault(raw_cluster, []).extend(candidates)
            delta_candidates += len(candidates)
        stats = {
            "enabled": True,
            "selected_clusters": len(selected),
            "searched_clusters": len(delta_records),
            "delta_patterns": sum(
                len(record[2]) for record in delta_records
            ),
            "delta_pattern_cap": MAX_DELTA_PATTERNS,
            "delta_pattern_cap_reached": cap_reached,
            "delta_pattern_counts": delta_pattern_counts,
            "trim_profile": self.fast_run_rescue_trim_profile,
            "delta_candidates": delta_candidates,
            "raw_seed_patterns": raw_seed_patterns,
            "selected_seed_patterns": selected_seed_patterns,
            "delta_topology": delta_telemetry,
        }
        self.fast_run_rescue_stats = stats
        return merged, stats

    def output(self):
        return self.dict_fuzzy_crisprs

    def two_repeat_terminal_repair_seeds(self):
        """Return deterministic opt-in seeds with a per-cluster fail-closed cap."""

        enabled = bool(
            self.flag_fast_run
            and self.accepted_only
            and (
                self.flag_fast_two_repeat_terminal_repair
                or self.flag_fast_accepted_family_rescue
            )
        )
        if not enabled:
            return (), dict(self.two_repeat_terminal_repair_seed_stats)

        retained = []
        captured = 0
        clusters_over_cap = 0
        dropped = 0
        for cluster_index, cluster_sequence in enumerate(
                self.cluster_sequences):
            cluster_seeds = []
            for candidate in self.dict_fuzzy_crisprs.get(
                    cluster_sequence, ()):
                if isinstance(candidate, RejectedFuzzyCandidateTopology):
                    cluster_seeds.extend(
                        candidate.two_repeat_terminal_repair_seeds
                    )
                elif isinstance(candidate, RejectedFuzzyCandidateWitness):
                    seed = candidate.two_repeat_terminal_repair_seed
                    if seed is not None:
                        cluster_seeds.append(seed)
            cluster_seeds = tuple(sorted(set(cluster_seeds)))
            captured += len(cluster_seeds)
            if (
                    len(cluster_seeds)
                    > MAX_TWO_REPEAT_TERMINAL_REPAIR_SEEDS_PER_CLUSTER):
                clusters_over_cap += 1
                dropped += len(cluster_seeds)
                continue
            retained.extend(
                TwoRepeatTerminalRepairSeed(
                    seed.repeat_starts,
                    seed.repeat_ends,
                    cluster_index,
                )
                for seed in cluster_seeds
            )
        retained = tuple(sorted(set(retained)))
        stats = {
            "enabled": True,
            "captured": captured,
            "retained": len(retained),
            "clusters_over_cap": clusters_over_cap,
            "seeds_dropped_by_cluster_cap": dropped,
            "seed_cap_per_cluster": (
                MAX_TWO_REPEAT_TERMINAL_REPAIR_SEEDS_PER_CLUSTER
            ),
        }
        self.two_repeat_terminal_repair_seed_stats = stats
        return retained, dict(stats)
