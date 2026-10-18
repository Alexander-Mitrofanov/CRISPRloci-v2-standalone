from os.path import basename
from components.components_evaluated_arrays_enhancement import IterativeDegeneratedSearch
from components.candidate_origins import inherit_candidate_origins
from components.components_evaluated_arrays_enhancement import create_boundaries_for_intervals
from components.components_evaluated_arrays_enhancement import ArrayRefinerInsertionsDeletions
from components.components_evaluated_arrays_enhancement import apply_source_backed_weak_terminal_repeat_proposal
from components.components_evaluated_arrays_enhancement import apply_source_backed_terminal_repeat_proposals
from components.components_evaluated_arrays_enhancement import source_backed_weak_terminal_repeat_proposals
from components.components_evaluated_arrays_enhancement import source_backed_terminal_repeat_proposals
from components.components_evaluated_arrays_enhancement import validate_candidate_source_reconstruction
from components.components_detection_refinement import AdvancedFuzzySearchFilter
from components.reconstruction_validator import read_single_fasta


def _candidate_source_interval(candidate):
    return (
        candidate.list_repeat_starts[0],
        candidate.list_repeat_starts[-1]
        + len(candidate.list_repeats[-1]),
    )


def _intervals_overlap(first, second):
    return first[0] < second[1] and second[0] < first[1]


def _proposal_added_span(source_interval, proposal):
    """Return the complete source span annexed by one proposal."""

    repeat_start, repeat_end = proposal["repeat_interval"]
    if proposal["side"] == "left":
        if repeat_end > source_interval[0]:
            raise ValueError(
                "left terminal proposal overlaps its input array"
            )
        return repeat_start, source_interval[0]
    if proposal["side"] == "right":
        if repeat_start < source_interval[1]:
            raise ValueError(
                "right terminal proposal overlaps its input array"
            )
        return source_interval[1], repeat_end
    raise ValueError("terminal proposal side must be left or right")


def build_source_backed_terminal_variants(full_dna, categories):
    """Build source-validated terminal variants without mutating scored rows.

    Only the currently accepted Bona-fide and Possible rows are eligible.
    The complete added span (repeat plus bridge spacer) may not cross another
    accepted array. Proposals whose added spans conflict with each other are
    rejected symmetrically, so the result does not depend on traversal order.

    The returned variants are deliberately uncommitted. The pipeline must
    score them together with the untouched post-legacy rows and commit only
    variants that remain accepted after interval-level categorization.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    if not isinstance(categories, (tuple, list)) or len(categories) != 5:
        raise ValueError(
            "terminal completion requires the five evaluation categories"
        )

    records = []
    for category_index, category_name in ((0, "Bona-fide"), (2, "Possible")):
        category = categories[category_index]
        for interval in sorted(category):
            for array_index, entry in enumerate(category[interval]):
                if not isinstance(entry, (tuple, list)) or len(entry) < 2:
                    raise ValueError(
                        "terminal completion received a malformed score row"
                    )
                candidate = entry[1]
                source_interval = _candidate_source_interval(candidate)
                records.append({
                    "category": category_name,
                    "interval": interval,
                    "array_index": array_index,
                    "base_entry": entry,
                    "base_candidate": candidate,
                    "source_interval": source_interval,
                    "proposals": source_backed_terminal_repeat_proposals(
                        full_dna, candidate
                    ),
                })

    proposed = sum(len(record["proposals"]) for record in records)
    proposal_records = []
    for record in records:
        owner_proposals = []
        for proposal in record["proposals"]:
            proposal_record = {
                "index": len(proposal_records),
                "owner": record,
                "proposal": proposal,
                "added_span": _proposal_added_span(
                    record["source_interval"], proposal
                ),
            }
            proposal_records.append(proposal_record)
            owner_proposals.append(proposal_record)
        record["proposal_records"] = tuple(owner_proposals)

    collision_reasons = {}
    for proposal_record in proposal_records:
        if any(
                other is not proposal_record["owner"]
                and _intervals_overlap(
                    proposal_record["added_span"],
                    other["source_interval"],
                )
                for other in records):
            collision_reasons.setdefault(
                proposal_record["index"], set()
            ).add("accepted_peer")

    peer_safe_proposals = [
        proposal_record for proposal_record in proposal_records
        if proposal_record["index"] not in collision_reasons
    ]
    for left_index, left in enumerate(peer_safe_proposals):
        for right in peer_safe_proposals[left_index + 1:]:
            if (
                    left["owner"] is not right["owner"]
                    and _intervals_overlap(
                        left["added_span"], right["added_span"]
                    )):
                collision_reasons.setdefault(
                    left["index"], set()
                ).add("proposal_peer")
                collision_reasons.setdefault(
                    right["index"], set()
                ).add("proposal_peer")

    variants = []
    variant_repeats = 0
    for record in records:
        safe_proposals = tuple(
            proposal_record["proposal"]
            for proposal_record in record["proposal_records"]
            if proposal_record["index"] not in collision_reasons
        )
        if not safe_proposals:
            continue
        variant = apply_source_backed_terminal_repeat_proposals(
            full_dna, record["base_candidate"], safe_proposals
        )
        validate_candidate_source_reconstruction(
            full_dna,
            variant.list_repeats,
            variant.list_repeat_starts,
            variant.list_spacers,
        )
        variants.append({
            "category": record["category"],
            "interval": record["interval"],
            "array_index": record["array_index"],
            "base_entry": record["base_entry"],
            "base_candidate": record["base_candidate"],
            "variant": variant,
            "source_interval": record["source_interval"],
            "safe_proposals": safe_proposals,
        })
        variant_repeats += len(safe_proposals)

    stats = {
        "candidate_count": len(records),
        "proposed_repeats": proposed,
        "collision_skipped_repeats": len(collision_reasons),
        "accepted_peer_collision_skipped_repeats": sum(
            "accepted_peer" in reasons
            for reasons in collision_reasons.values()
        ),
        "proposal_collision_skipped_repeats": sum(
            "proposal_peer" in reasons
            for reasons in collision_reasons.values()
        ),
        "variants_built": len(variants),
        "variant_repeats": variant_repeats,
    }
    return tuple(variants), stats


def _candidate_structure_key(candidate):
    """Return the complete immutable representation used for deduplication."""

    return (
        tuple(candidate.list_repeat_starts),
        tuple(candidate.list_repeats),
        tuple(candidate.list_repeats_gaped),
        tuple(candidate.list_spacers),
    )


def build_source_backed_terminal_candidate_variant_records(
        full_dna, candidates_by_locus, eligible_candidate_ids=None,
        include_weak=False):
    """Build ordered, collision-safe variants for selective rescoring.

    This is a policy-free construction primitive, not an active pipeline
    stage.  Eligibility is supplied exclusively as candidate object IDs.  All
    candidates are source-validated and retained as collision obstacles, but
    the terminal proposal rule is called only for eligible originals.  Passing
    ``None`` preserves the original all-candidate behavior.

    For each eligible source-valid original, the existing one-step terminal
    rule may contribute at most one combined left/right variant. When
    ``include_weak`` is true and the ordinary rule finds nothing, an isolated
    source-only weak rule may contribute exactly one terminal. Weak variants
    are labelled for fail-closed coding-context authorization downstream.
    Returned records follow locus insertion order and candidate order, and
    contain the base, variant, locus, and exact safe proposals needed for
    selective scoring.

    The complete span annexed by a proposal (bridge spacer plus repeat) may
    not overlap a candidate owned by another locus.  Proposals from different
    loci that overlap one another are rejected on both sides.  Alternative
    representations within one locus deliberately do not block each other;
    they still compete during ordinary evaluation.  Exact duplicate candidate
    structures are not returned.
    """

    if not isinstance(full_dna, str):
        raise ValueError("source DNA must be a string")
    if not isinstance(candidates_by_locus, dict):
        raise ValueError("terminal augmentation requires a candidate mapping")
    if not isinstance(include_weak, bool):
        raise ValueError("include_weak must be a boolean")
    eligibility_limited = eligible_candidate_ids is not None
    if eligibility_limited:
        try:
            eligible_candidate_ids = frozenset(eligible_candidate_ids)
        except TypeError as error:
            raise ValueError(
                "eligible candidate IDs must be an iterable of object IDs"
            ) from error
        if any(
                not isinstance(candidate_id, int)
                or isinstance(candidate_id, bool)
                or candidate_id < 0
                for candidate_id in eligible_candidate_ids):
            raise ValueError(
                "eligible candidate IDs must be non-negative integers"
            )

    records = []
    for locus_key, candidates in candidates_by_locus.items():
        if not isinstance(candidates, (tuple, list)):
            raise ValueError(
                "each terminal-augmentation locus must contain a candidate list"
            )
        for candidate_index, candidate in enumerate(candidates):
            validate_candidate_source_reconstruction(
                full_dna,
                candidate.list_repeats,
                candidate.list_repeat_starts,
                candidate.list_spacers,
            )
            source_interval = _candidate_source_interval(candidate)
            records.append({
                "locus_key": locus_key,
                "candidate_index": candidate_index,
                "base_candidate": candidate,
                "source_interval": source_interval,
            })

    observed_candidate_ids = {
        id(record["base_candidate"]) for record in records
    }
    for record in records:
        candidate = record["base_candidate"]
        eligible = (
            not eligibility_limited
            or id(candidate) in eligible_candidate_ids
        )
        record["eligible"] = eligible
        proposals = ()
        terminal_tier = None
        if eligible:
            proposals = source_backed_terminal_repeat_proposals(
                full_dna, candidate
            )
            if proposals:
                terminal_tier = "strong"
            elif include_weak:
                proposals = source_backed_weak_terminal_repeat_proposals(
                    full_dna, candidate
                )
                if proposals:
                    terminal_tier = "weak"
        record["proposals"] = proposals
        record["terminal_tier"] = terminal_tier

    proposal_records = []
    for record in records:
        owner_proposals = []
        for proposal in record["proposals"]:
            proposal_record = {
                "index": len(proposal_records),
                "owner": record,
                "proposal": proposal,
                "added_span": _proposal_added_span(
                    record["source_interval"], proposal
                ),
            }
            proposal_records.append(proposal_record)
            owner_proposals.append(proposal_record)
        record["proposal_records"] = tuple(owner_proposals)

    collision_reasons = {}
    for proposal_record in proposal_records:
        owner = proposal_record["owner"]
        if any(
                other["locus_key"] != owner["locus_key"]
                and _intervals_overlap(
                    proposal_record["added_span"],
                    other["source_interval"],
                )
                for other in records):
            collision_reasons.setdefault(
                proposal_record["index"], set()
            ).add("candidate_peer")

    peer_safe_proposals = [
        proposal_record for proposal_record in proposal_records
        if proposal_record["index"] not in collision_reasons
    ]
    for left_index, left in enumerate(peer_safe_proposals):
        for right in peer_safe_proposals[left_index + 1:]:
            if (
                    left["owner"]["locus_key"]
                    != right["owner"]["locus_key"]
                    and _intervals_overlap(
                        left["added_span"], right["added_span"]
                    )):
                collision_reasons.setdefault(
                    left["index"], set()
                ).add("proposal_peer")
                collision_reasons.setdefault(
                    right["index"], set()
                ).add("proposal_peer")

    seen_structures = {
        locus_key: {
            _candidate_structure_key(candidate) for candidate in candidates
        }
        for locus_key, candidates in candidates_by_locus.items()
    }
    variant_records = []
    variants_built = 0
    weak_variants_built = 0
    duplicate_variants_skipped = 0
    weak_duplicate_variants_skipped = 0
    variant_repeats = 0
    for record in records:
        safe_proposals = tuple(
            proposal_record["proposal"]
            for proposal_record in record["proposal_records"]
            if proposal_record["index"] not in collision_reasons
        )
        if not safe_proposals:
            continue
        if record["terminal_tier"] == "weak":
            if len(safe_proposals) != 1:
                raise RuntimeError(
                    "a weak terminal variant must contain exactly one proposal"
                )
            variant = apply_source_backed_weak_terminal_repeat_proposal(
                full_dna,
                record["base_candidate"],
                safe_proposals[0],
            )
            weak_variants_built += 1
        else:
            variant = apply_source_backed_terminal_repeat_proposals(
                full_dna, record["base_candidate"], safe_proposals
            )
        validate_candidate_source_reconstruction(
            full_dna,
            variant.list_repeats,
            variant.list_repeat_starts,
            variant.list_spacers,
        )
        variants_built += 1
        structure_key = _candidate_structure_key(variant)
        locus_structures = seen_structures[record["locus_key"]]
        if structure_key in locus_structures:
            duplicate_variants_skipped += 1
            if record["terminal_tier"] == "weak":
                weak_duplicate_variants_skipped += 1
            continue
        locus_structures.add(structure_key)
        variant_records.append({
            "locus_key": record["locus_key"],
            "candidate_index": record["candidate_index"],
            "base_candidate": record["base_candidate"],
            "variant": variant,
            "source_interval": record["source_interval"],
            "safe_proposals": safe_proposals,
            "terminal_tier": record["terminal_tier"],
        })
        variant_repeats += len(safe_proposals)

    eligible_candidate_count = sum(
        record["eligible"] for record in records
    )
    stats = {
        "candidate_count": len(records),
        "eligibility_limited": eligibility_limited,
        "eligible_candidate_count": eligible_candidate_count,
        "ineligible_candidate_count": (
            len(records) - eligible_candidate_count
        ),
        "unmatched_eligible_candidate_ids": (
            len(eligible_candidate_ids.difference(observed_candidate_ids))
            if eligibility_limited else 0
        ),
        "proposed_repeats": len(proposal_records),
        "collision_skipped_repeats": len(collision_reasons),
        "candidate_peer_collision_skipped_repeats": sum(
            "candidate_peer" in reasons
            for reasons in collision_reasons.values()
        ),
        "proposal_collision_skipped_repeats": sum(
            "proposal_peer" in reasons
            for reasons in collision_reasons.values()
        ),
        "variants_built": variants_built,
        "variants_added": len(variant_records),
        "variant_repeats": variant_repeats,
        "duplicate_variants_skipped": duplicate_variants_skipped,
    }
    if include_weak:
        stats.update({
            "weak_proposed_repeats": sum(
                record["owner"]["terminal_tier"] == "weak"
                for record in proposal_records
            ),
            "weak_variants_built": weak_variants_built,
            "weak_variants_added": sum(
                record["terminal_tier"] == "weak"
                for record in variant_records
            ),
            "weak_duplicate_variants_skipped": (
                weak_duplicate_variants_skipped
            ),
        })
    return tuple(variant_records), stats




class EvaluatedArraysEnhancement:
    def __init__(self, file_path, categories, parameters, flag_dev_mode,
                 flag_degenerated=True):
        self.file_path = file_path
        self.categories = categories
        self.parameters = parameters
        self.flag_dev_mode = flag_dev_mode
        self.flag_degenerated = flag_degenerated
        self.bona_fide_arrays = categories[0]
        self.alternative_arrays = categories[1]
        self.possible_arrays = categories[2]

        self.dict_arrays_into_categories_enhanced = {}
        self.source_geometry_lineage = []

        self._get_complete_dna()
        if self.flag_degenerated:
            self._search_missed_or_degenerated_repeats()
        self._refine_nucleotides_repeat_spacer()
        self._filter_enhanced()

    def _get_complete_dna(self):
        header, self.dna = read_single_fasta(self.file_path)
        self.input_header = ">{}\n".format(header)
        self.dna_length = len(self.dna)

    @staticmethod
    def _candidate_geometry(candidate):
        return tuple(
            (start, start + len(repeat))
            for start, repeat in zip(
                candidate.list_repeat_starts,
                candidate.list_repeats,
            )
        )

    @staticmethod
    def _derived_row(base_row, candidate):
        if not isinstance(base_row, (list, tuple)) or len(base_row) not in (3, 4):
            raise ValueError("enhancement rows must contain three or four fields")
        derived = list(base_row[:3])
        derived[1] = candidate
        if len(base_row) == 4:
            if not isinstance(base_row[3], dict):
                raise ValueError("enhancement row metadata must be a mapping")
            derived.append(dict(base_row[3]))
        return derived

    def _retain_derived_geometry(
            self, category_index, category, locus_key, base_row,
            derived_candidate, stage):
        """Append a unique source-valid variant without replacing its parent."""

        validate_candidate_source_reconstruction(
            self.dna,
            derived_candidate.list_repeats,
            derived_candidate.list_repeat_starts,
            derived_candidate.list_spacers,
        )
        base_candidate = base_row[1]
        if derived_candidate is not base_candidate:
            inherit_candidate_origins(derived_candidate, (base_candidate,))
        base_key = _candidate_structure_key(base_candidate)
        derived_key = _candidate_structure_key(derived_candidate)
        if derived_key == base_key:
            return base_row

        rows = category[locus_key]
        matching_rows = [
            row for row in rows
            if _candidate_structure_key(row[1]) == derived_key
        ]
        if matching_rows:
            derived_row = matching_rows[0]
            inherit_candidate_origins(derived_row[1], (derived_row[1], derived_candidate))
            action = "deduplicated_existing_geometry"
        else:
            derived_row = self._derived_row(base_row, derived_candidate)
            rows.append(derived_row)
            action = "retained_derived_geometry"

        self.source_geometry_lineage.append({
            "stage": stage,
            "action": action,
            "category_index": category_index,
            "locus_key": locus_key,
            "base_row": base_row,
            "base_candidate": base_candidate,
            "base_geometry": self._candidate_geometry(base_candidate),
            "derived_row": derived_row,
            "derived_candidate": derived_row[1],
            "derived_geometry": self._candidate_geometry(derived_row[1]),
            "coordinate_system": "zero_based_half_open",
        })
        return derived_row

    def _search_missed_or_degenerated_repeats(self):
        for category_index, category in enumerate([
                self.bona_fide_arrays,
                self.alternative_arrays,
                self.possible_arrays,
        ]):
            intervals = []
            arrays_for_intervals = []

            for interval, list_data in category.items():
                intervals.append(interval)
                arrays_for_intervals.append(list(list_data))

            boundaries = create_boundaries_for_intervals(
                intervals, 500, source_length=self.dna_length
            )

            for interval, rows_in_interval, boundary in zip(
                    intervals, arrays_for_intervals, boundaries):
                for base_row in rows_in_interval:
                    array = base_row[1]
                    consensus = array.consensus
                    list_repeats = array.list_repeats
                    list_repeats_starts = array.list_repeat_starts
                    list_spacers = array.list_spacers


                    ids = IterativeDegeneratedSearch(full_dna=self.dna,
                                                     repeat_seq_candidate=consensus,
                                                     spacer_margin=self.parameters["param_spacer_margin_degenerated_search"],
                                                     repeat_seq_candidate_gaped=None,
                                                     list_repeats_starts=list_repeats_starts,
                                                     list_repeats=list_repeats,
                                                     list_spacers=list_spacers,
                                                     start_flanking_region_left=boundary[0],
                                                     end_flanking_region_right=boundary[1],
                                                     allowed_max_editing_distance=self.parameters["param_max_edit_distance"],
                                                     iterative_size_flanking_region=150,
                                                     prevent_long_spacers=True,
                                                     attempt_to_improve_initial_array=True)

                    new_crispr_candidate = ids.output()
                    if self.flag_dev_mode:
                        if array != new_crispr_candidate:
                            with open("log.txt", "a") as f:
                                acc_num = basename(self.file_path).split(".")[0]
                                f.write(f"Iteractive degenerated search {acc_num}\n")
                                f.write(array.dot_repr())
                                f.write("\n\n")
                                f.write(new_crispr_candidate.dot_repr())
                                f.write("\n\n")

                    """except Exception:
                        new_crispr_candidate = array

                        if self.flag_dev_mode:
                            with open("log_error.txt", "a") as f:
                                acc_num = basename(self.file_path).split(".")[0]
                                f.write(f"Iteractive degenerated search error {acc_num}\n")
                                f.write(array.dot_repr())
                                f.write("\n\n")"""

                    self._retain_derived_geometry(
                        category_index,
                        category,
                        interval,
                        base_row,
                        new_crispr_candidate,
                        "iterative_degenerated_search",
                    )

    def _refine_nucleotides_repeat_spacer(self):
        for category_index, category in enumerate([
                self.bona_fide_arrays,
                self.alternative_arrays,
                self.possible_arrays,
        ]):
            for interval, list_data in category.items():
                base_rows = list(list_data)
                for base_row in base_rows:
                    array = base_row[1]
                    arid = ArrayRefinerInsertionsDeletions(array)
                    new_crispr_candidate = arid.output()

                    if self.flag_dev_mode:
                        if array != new_crispr_candidate:
                            with open("log.txt", "a") as f:
                                acc_num = basename(self.file_path).split(".")[0]
                                f.write(f"Array refinement {acc_num}\n")
                                f.write(array.dot_repr())
                                f.write("\n\n")
                                f.write(new_crispr_candidate.dot_repr())
                                f.write("\n\n")

                    self._retain_derived_geometry(
                        category_index,
                        category,
                        interval,
                        base_row,
                        new_crispr_candidate,
                        "insertion_deletion_refinement",
                    )

    def _filter_enhanced(self):
        self.param_min_avg_repeat_length = self.parameters["param_min_avg_repeat_length"]
        self.param_max_avg_repeat_length = self.parameters["param_max_avg_repeat_length"]
        self.param_min_avg_spacer_length = self.parameters["param_min_avg_spacer_length"]
        self.param_max_avg_spacer_length = self.parameters["param_max_avg_spacer_length"]
        self.param_min_repeats = self.parameters["param_min_repeats"]
        self.param_max_identical_spacers = self.parameters["param_max_identical_spacers"]
        self.param_max_identical_cluster_spacers = self.parameters["param_max_identical_cluster_spacers"]

        afsf = AdvancedFuzzySearchFilter(min_column_dominance_repeat=0.6,
                                         max_spacer_length=140, max_column_dominance_spacer=0.8,
                                         max_allowed_consecutive_spacers=self.param_max_identical_cluster_spacers,
                                         max_allowed_same_spacers=self.param_max_identical_spacers,
                                         max_inconsistent_columns=5,
                                         min_avg_repeat_length=self.param_min_avg_repeat_length,
                                         max_avg_repeat_length=self.param_max_avg_repeat_length,
                                         min_avg_spacer_length=self.param_min_avg_spacer_length,
                                         max_avg_spacer_length=self.param_max_avg_spacer_length,
                                         min_repeats=self.param_min_repeats)

        bona_fide_not_filtered = self.categories[0]
        alternative_not_filtered = self.categories[1]
        possible_not_filtered = self.categories[2]
        low_score = self.categories[4]

        bona_fide_filtered = {}
        alternative_filtered = {}
        possible_filtered = {}

        for not_filtered_category, filtered_category in zip([bona_fide_not_filtered, alternative_not_filtered, possible_not_filtered],
                                                            [bona_fide_filtered, alternative_filtered, possible_filtered]):
            for key, value in not_filtered_category.items():
                for crispr_tuple in value:
                    crispr = crispr_tuple[1]
                    if not afsf(crispr):
                        if key in low_score:
                            low_score[key].append(crispr_tuple)
                        else:
                            low_score[key] = [crispr_tuple]
                    else:
                        if key not in filtered_category:
                            filtered_category[key] = [crispr_tuple]
                        else:
                            filtered_category[key].append(crispr_tuple)

        self.categories[0] = bona_fide_filtered
        self.categories[1] = alternative_filtered
        self.categories[2] = possible_filtered
        self.categories[4] = low_score

    def output(self):
        return self.categories
