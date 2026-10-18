import collections
from dataclasses import asdict, dataclass
import json
from functools import wraps
from itertools import groupby
from components.candidate_origins import mark_ordinary_birth

from components.components_detection import FuzzySourcePathProvenance
from components.components_detection import RejectedFuzzyCandidateWitness
from components.components_detection import RejectedFuzzyCandidateTopology


def merged_source_path_provenance(*candidates):
    """Return the deterministic union of source-path records, if present."""

    if len(candidates) == 1:
        candidate = candidates[0]
        provenance = getattr(candidate, "source_path_provenance", None)
        if (
                isinstance(provenance, tuple)
                and len(provenance) == 1
                and isinstance(provenance[0], FuzzySourcePathProvenance)):
            validator = getattr(
                candidate, "_validate_source_path_provenance", None
            )
            if validator is not None:
                if not callable(validator):
                    raise ValueError("source path validator is not callable")
                validator()
            return provenance

    records = []
    provenance_present = False
    for candidate in candidates:
        provenance = getattr(candidate, "source_path_provenance", None)
        if provenance is None:
            continue
        validator = getattr(
            candidate, "_validate_source_path_provenance", None
        )
        if validator is not None:
            if not callable(validator):
                raise ValueError("source path validator is not callable")
            validator()
        provenance_present = True
        if (
                not isinstance(provenance, tuple)
                or any(
                    not isinstance(record, FuzzySourcePathProvenance)
                    for record in provenance
                )):
            raise ValueError("source path provenance is malformed")
        records.extend(provenance)
    if not provenance_present:
        return None
    try:
        return tuple(sorted(set(records)))
    except TypeError as error:
        raise ValueError("source path provenance is not immutable") from error


class SameStartEndFilter:
    def __init__(self, dict_crispr_candidates):
        self.dict_crispr_candidates = dict_crispr_candidates
        self.dict_filtered_start_end_crispr_candidates = {}

        self._filter_fuzzy_searches_same_start_end()

    def _filter_fuzzy_searches_same_start_end(self):
        for cluster_seq, list_fuzzy_s in self.dict_crispr_candidates.items():
            # Only collapse candidates whose complete source geometry and
            # alignment representation are identical. Sharing the first/last
            # repeat start or the seed repeat is not enough: those shortcuts
            # can discard candidates with different internal or terminal
            # coordinates.
            best_by_geometry = {}
            topology_markers = []
            for fuzzy_s in list_fuzzy_s:
                if isinstance(
                        fuzzy_s, RejectedFuzzyCandidateTopology):
                    topology_markers.append(fuzzy_s)
                    continue
                if isinstance(fuzzy_s, RejectedFuzzyCandidateWitness):
                    key = (
                        "compact_fuzzy_alignment",
                        fuzzy_s.same_start_end_key,
                    )
                else:
                    compact_key = getattr(
                        fuzzy_s, "compact_same_start_end_key", None
                    )
                    key = (
                        (
                            "compact_fuzzy_alignment",
                            compact_key(),
                        )
                        if callable(compact_key)
                        else (
                            "legacy_dense_alignment",
                            tuple(fuzzy_s.list_absolute_start),
                            tuple(fuzzy_s.list_repeats),
                            tuple(fuzzy_s.list_spacers),
                            tuple(fuzzy_s.list_gaped_repeats),
                            getattr(
                                fuzzy_s, "gaped_matching_pattern", None
                            ),
                        )
                    )
                current_best = best_by_geometry.get(key)
                if current_best is None:
                    best_by_geometry[key] = fuzzy_s
                    continue

                current_is_witness = isinstance(
                    current_best, RejectedFuzzyCandidateWitness
                )
                incoming_is_witness = isinstance(
                    fuzzy_s, RejectedFuzzyCandidateWitness
                )
                if current_is_witness or incoming_is_witness:
                    if current_is_witness != incoming_is_witness:
                        raise ValueError(
                            "compact structural decision changed inside one "
                            "SameStartEnd equivalence class"
                        )
                    family_parent = (
                        current_best.source_path_family_parent
                        or fuzzy_s.source_path_family_parent
                    )
                    if fuzzy_s.number_errors < current_best.number_errors:
                        current_best = fuzzy_s
                    best_by_geometry[key] = current_best._replace(
                        source_path_family_parent=family_parent
                    )
                    continue

                merged_provenance = merged_source_path_provenance(
                    current_best, fuzzy_s
                )
                all_reasons = []
                for member in (current_best, fuzzy_s):
                    prior = getattr(member, "source_path_restriction_reasons", ())
                    if not isinstance(prior, (list, tuple)):
                        prior = (None,)
                    all_reasons.extend(prior)
                    if hasattr(member, "source_path_classification_reason"):
                        all_reasons.append(member.source_path_classification_reason)
                restricted = [
                    candidate for candidate in (current_best, fuzzy_s)
                    if hasattr(candidate, "source_path_classification_eligible")
                    and (type(candidate.source_path_classification_eligible) is not bool
                         or not candidate.source_path_classification_eligible)
                ]
                if fuzzy_s.number_errors < current_best.number_errors:
                    current_best = fuzzy_s
                    best_by_geometry[key] = current_best
                if all_reasons:
                    current_best.source_path_restriction_reasons = tuple(all_reasons)
                if merged_provenance is not None:
                    current_best.source_path_provenance = merged_provenance
                if restricted:
                    restriction = next((candidate for candidate in restricted
                                        if type(candidate.source_path_classification_eligible)
                                        is not bool), restricted[0])
                    current_best.source_path_classification_eligible = (
                        restriction.source_path_classification_eligible
                    )
                    current_best.source_path_classification_reason = getattr(
                        restriction, "source_path_classification_reason",
                        "inherited_source_path_restriction",
                    )

            filtered = list(best_by_geometry.values())
            if topology_markers:
                filtered.append(RejectedFuzzyCandidateTopology.combine(
                    topology_markers
                ))
            self.dict_filtered_start_end_crispr_candidates[cluster_seq] = (
                filtered
            )

    def output(self):
        return self.dict_filtered_start_end_crispr_candidates


#          For filtering out non CRISPR cases
#############################################################
#############################################################
DEBUG_MODE = False


@dataclass(frozen=True)
class SpacerTerminalDominanceEvidence:
    """Finite-sample evidence for spacer terminal-base dominance."""

    status: str
    reason: str
    nonempty_spacer_count: int
    minimum_sample_size: int
    maximum_allowed_dominance: float
    begin_dominance: object
    end_dominance: object

    @property
    def passes_filter(self):
        return self.status in ("passed", "insufficient_sample")

    def as_dict(self):
        return asdict(self)


@dataclass(frozen=True)
class AdvancedFuzzySearchFilterEvidence:
    """Exhaustive structural-filter evidence for one candidate."""

    decision: str
    structural_valid: bool
    predicate_results: tuple
    terminal_dominance: object

    @property
    def passes(self):
        return self.decision == "passed"

    def as_dict(self):
        return {
            "decision": self.decision,
            "structural_valid": self.structural_valid,
            "predicate_results": [
                {"predicate": name, "status": status}
                for name, status in self.predicate_results
            ],
            "terminal_dominance": (
                None
                if self.terminal_dominance is None
                else self.terminal_dominance.as_dict()
            ),
        }


def printing_if_filtered(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        result = function(*args, **kwargs)
        if DEBUG_MODE:
            if not result:
                with open("filtered_results.txt", "a") as f:
                    f.write("\n\n")
                    f.write("\n".join([str(arg) for arg in args]))
                    f.write("\n\n")
                    f.write(function.__name__)
                    f.write("\n\n")

        return result
    return wrapper


class AdvancedFuzzySearchFilter:
    # With two spacers, two independent random terminal bases agree 25% of
    # the time. Applying the 0.8 dominance gate at that sample size therefore
    # rejects about 43.75% of otherwise unrelated begin/end pairs by chance.
    # Five observations are the first size at which ``> 0.8`` requires all
    # bases to agree and the two-ended chance rejection rate falls below 1%.
    MIN_SPACERS_FOR_TERMINAL_DOMINANCE = 5

    @classmethod
    def from_parameters(cls, parameters):
        """Build the canonical detector/evaluator structural conjunction."""

        required = (
            "param_min_avg_repeat_length",
            "param_max_avg_repeat_length",
            "param_min_avg_spacer_length",
            "param_max_avg_spacer_length",
            "param_min_repeats",
            "param_max_identical_spacers",
            "param_max_identical_cluster_spacers",
        )
        missing = tuple(name for name in required if name not in parameters)
        if missing:
            raise ValueError(
                "structural filter parameters are incomplete: {}".format(
                    ", ".join(missing)
                )
            )
        return cls(
            min_column_dominance_repeat=0.6,
            max_spacer_length=140,
            max_column_dominance_spacer=0.8,
            max_allowed_consecutive_spacers=parameters[
                "param_max_identical_cluster_spacers"
            ],
            max_allowed_same_spacers=parameters[
                "param_max_identical_spacers"
            ],
            max_inconsistent_columns=5,
            min_avg_repeat_length=parameters[
                "param_min_avg_repeat_length"
            ],
            max_avg_repeat_length=parameters[
                "param_max_avg_repeat_length"
            ],
            min_avg_spacer_length=parameters[
                "param_min_avg_spacer_length"
            ],
            max_avg_spacer_length=parameters[
                "param_max_avg_spacer_length"
            ],
            min_repeats=parameters["param_min_repeats"],
        )

    def __init__(self, min_column_dominance_repeat, min_avg_spacer_length,
                 max_spacer_length, max_column_dominance_spacer, max_allowed_consecutive_spacers,
                 max_allowed_same_spacers, max_inconsistent_columns, min_avg_repeat_length,
                 max_avg_repeat_length, max_avg_spacer_length, min_repeats):

        self.column_dominance = min_column_dominance_repeat
        self.min_avg_spacer_length = min_avg_spacer_length
        self.max_spacer_length = max_spacer_length
        self.max_column_dominance_spacer = max_column_dominance_spacer
        self.max_allowed_consecutive_spacers = max_allowed_consecutive_spacers
        self.max_allowed_same_spacers = max_allowed_same_spacers
        self.max_inconsistent_columns = max_inconsistent_columns
        self.min_avg_repeat_length = min_avg_repeat_length
        self.max_avg_repeat_length = max_avg_repeat_length
        self.max_avg_spacer_length = max_avg_spacer_length
        self.min_number_repeats = min_repeats

    @printing_if_filtered
    def _filter_by_column(self, candidate):
        def find_first_three_columns():
            list_three_columns = []
            list_gaped_repeats = candidate.list_gaped_repeats
            for index in range(len(list_gaped_repeats[0])):
                column_vec = [repeat[index] for repeat in list_gaped_repeats]
                column_gaps = sum([1 for x in column_vec if (x == " ")])
                percentage_gaps = column_gaps / len(column_vec)
                if percentage_gaps < 0.5:
                    list_three_columns.append(column_vec)
                if len(list_three_columns) == 3:
                    return list_three_columns
            return list_three_columns

        def find_last_three_columns():
            list_three_columns = []
            list_gaped_repeats = candidate.list_gaped_repeats
            for index in range(len(list_gaped_repeats[0]) - 1, -1, -1):
                column_vec = [repeat[index] for repeat in list_gaped_repeats]
                column_gaps = sum([1 for x in column_vec if (x == " ")])
                percentage_gaps = column_gaps/len(column_vec)
                if percentage_gaps < 0.5:
                    list_three_columns.append(column_vec)
                if len(list_three_columns) == 3:
                    return list_three_columns
            return list_three_columns

        first_columns = find_first_three_columns()
        last_columns = find_last_three_columns()
        if len(first_columns) != 3 or len(last_columns) != 3:
            return False

        for column in first_columns:
            column_characters = [x for x in column if (x not in (" ", "-"))]
            if column_characters:
                most_freq_char = max(column_characters, key=column_characters.count)
                most_freq_char_freq = column_characters.count(most_freq_char)
                freq = most_freq_char_freq/len(column_characters)
                if len(column) <= 4:
                    if freq < 0.49:
                        return False
                else:
                    if freq < self.column_dominance:
                        return False
            else:
                return False

        for column in last_columns:
            column_characters = [x for x in column if (x not in (" ", "-"))]
            if column_characters:
                most_freq_char = max(column_characters, key=column_characters.count)
                most_freq_char_freq = column_characters.count(most_freq_char)
                freq = most_freq_char_freq/len(column_characters)
                if len(column) <= 4:
                    if freq < 0.49:
                        return False
                else:
                    if freq < self.column_dominance:
                        return False
            else:
                return False
        return True

    @printing_if_filtered
    def _filter_by_min_avg_spacer(self, candidate):
        list_spacers = candidate.list_spacers
        avg_len = sum(len(x) for x in list_spacers) / len(list_spacers)
        if avg_len >= self.min_avg_spacer_length:
            return True
        return False

    @printing_if_filtered
    def _filter_by_max_spacer(self, candidate):
        list_spacers = candidate.list_spacers
        long_spacers = [spacer for spacer in list_spacers if len(spacer) > self.max_spacer_length]
        if len(long_spacers) / len(list_spacers) > 0.3:
            return False
        if len(long_spacers) > 3:
            return False
        return True

    @printing_if_filtered
    def _filter_by_spacer_begin_end_similarity(self, candidate):
        return self.spacer_terminal_dominance_evidence(
            candidate
        ).passes_filter

    def spacer_terminal_dominance_evidence(self, candidate):
        """Return pass/fail/insufficient-sample evidence without hiding N."""

        nonempty_spacers = [
            spacer for spacer in candidate.list_spacers if spacer
        ]
        if not nonempty_spacers:
            return SpacerTerminalDominanceEvidence(
                status="failed",
                reason="no_nonempty_spacers",
                nonempty_spacer_count=0,
                minimum_sample_size=self.MIN_SPACERS_FOR_TERMINAL_DOMINANCE,
                maximum_allowed_dominance=self.max_column_dominance_spacer,
                begin_dominance=None,
                end_dominance=None,
            )

        def dominance(values):
            most_frequent = max(values, key=values.count)
            return values.count(most_frequent) / float(len(values))

        begin_dominance = dominance([
            spacer[0] for spacer in nonempty_spacers
        ])
        end_dominance = dominance([
            spacer[-1] for spacer in nonempty_spacers
        ])
        if (
                len(nonempty_spacers)
                < self.MIN_SPACERS_FOR_TERMINAL_DOMINANCE):
            status = "insufficient_sample"
            reason = "terminal_dominance_not_authorized_at_small_n"
        elif (
                begin_dominance > self.max_column_dominance_spacer
                or end_dominance > self.max_column_dominance_spacer):
            status = "failed"
            reason = "terminal_base_dominance_above_threshold"
        else:
            status = "passed"
            reason = "terminal_base_dominance_within_threshold"
        return SpacerTerminalDominanceEvidence(
            status=status,
            reason=reason,
            nonempty_spacer_count=len(nonempty_spacers),
            minimum_sample_size=self.MIN_SPACERS_FOR_TERMINAL_DOMINANCE,
            maximum_allowed_dominance=self.max_column_dominance_spacer,
            begin_dominance=begin_dominance,
            end_dominance=end_dominance,
        )

    @printing_if_filtered
    def _filter_by_the_same_spacer(self, candidate):
        list_spacers = candidate.list_spacers
        list_spacers = [s for s in list_spacers if s]
        if not list_spacers:
            return False
        groups = [len(list(group)) for key, group in groupby(list_spacers)]
        if self.max_allowed_consecutive_spacers:
            if max(groups) > self.max_allowed_consecutive_spacers:
                return False

        list_sorted_spacers = sorted(list_spacers)
        groups_sorted = [len(list(group)) for key, group in groupby(list_sorted_spacers)]
        if self.max_allowed_same_spacers:
            if max(groups_sorted) > self.max_allowed_same_spacers:
                return False
        return True

    @printing_if_filtered
    def _filter_by_overall_repeat_consistency(self, candidate):
        list_column_consistency = []
        list_repeats_gaped = candidate.list_gaped_repeats
        for index, _ in enumerate(list_repeats_gaped[0]):
            column = [repeat[index] for repeat in list_repeats_gaped]
            column_characters = [x for x in column if (x not in (" ", "-"))]
            try:
                most_freq_char = max(column_characters, key=column_characters.count)
                most_freq_char_freq = column_characters.count(most_freq_char)
                freq = most_freq_char_freq / len(column_characters)
                list_column_consistency.append(freq)
            except ValueError:
                pass

        number_inconsistent = sum(1 for x in list_column_consistency if x < 0.66)
        if number_inconsistent > self.max_inconsistent_columns:
            return False
        return True

    @printing_if_filtered
    def _filter_min_number_repeats(self, candidate):
        list_repeats = candidate.list_repeats
        if len(list_repeats) >= self.min_number_repeats:
            return True
        return False

    @printing_if_filtered
    def _filter_min_avg_repeat_length(self, candidate):
        list_repeats = candidate.list_repeats
        avg_len = sum(len(x) for x in list_repeats) / len(list_repeats)
        if avg_len >= self.min_avg_repeat_length:
            return True
        return False

    @printing_if_filtered
    def _filter_max_avg_repeat_length(self, candidate):
        list_repeats = candidate.list_repeats
        avg_len = sum(len(x) for x in list_repeats) / len(list_repeats)
        if avg_len <= self.max_avg_repeat_length:
            return True
        return False

    @printing_if_filtered
    def _filter_max_avg_spacer_length(self, candidate):
        list_spacers = candidate.list_spacers
        spacers_to_average = (
            list_spacers[1:-1] if len(list_spacers) > 4 else list_spacers
        )
        if not spacers_to_average:
            return False
        avg_len = (
            sum(len(spacer) for spacer in spacers_to_average)
            / len(spacers_to_average)
        )
        return avg_len <= self.max_avg_spacer_length

    @printing_if_filtered
    def _filter_min_repeat_length(self, candidate):
        list_repeats = candidate.list_repeats
        avg_len = sum(len(x) for x in list_repeats) / len(list_repeats)
        if avg_len >= self.min_avg_repeat_length:
            return True
        return False

    @staticmethod
    def _candidate_structure_is_valid(candidate):
        """Reject known malformed geometry without hiding programming errors."""
        repeats = candidate.list_repeats
        gapped_repeats = candidate.list_gaped_repeats
        spacers = candidate.list_spacers

        if not repeats or not gapped_repeats:
            return False
        if len(gapped_repeats) != len(repeats):
            return False
        if len(spacers) != len(repeats) - 1:
            return False
        if any(not isinstance(sequence, str) for sequence in repeats):
            return False
        if any(not isinstance(sequence, str) for sequence in gapped_repeats):
            return False
        if any(not isinstance(sequence, str) for sequence in spacers):
            return False
        if any(not repeat for repeat in repeats):
            return False
        gapped_lengths = {len(repeat) for repeat in gapped_repeats}
        return len(gapped_lengths) == 1 and next(iter(gapped_lengths)) > 0

    def _passes_non_gapped_predicates(self, candidate):
        """Evaluate mandatory predicates that do not need gapped repeats."""

        predicates = (
            self._filter_min_number_repeats,
            self._filter_min_avg_repeat_length,
            self._filter_max_avg_repeat_length,
            self._filter_max_avg_spacer_length,
            self._filter_by_min_avg_spacer,
            self._filter_by_max_spacer,
            self._filter_by_the_same_spacer,
            self._filter_by_spacer_begin_end_similarity,
        )
        return all(predicate(candidate) for predicate in predicates)

    def evidence(self, candidate):
        """Evaluate every named structural gate for audit/report provenance."""

        structural_valid = self._candidate_structure_is_valid(candidate)
        if not structural_valid:
            return AdvancedFuzzySearchFilterEvidence(
                decision="rejected",
                structural_valid=False,
                predicate_results=(("candidate_structure", "failed"),),
                terminal_dominance=None,
            )

        terminal = self.spacer_terminal_dominance_evidence(candidate)
        predicates = (
            ("minimum_repeat_count", self._filter_min_number_repeats),
            ("minimum_mean_repeat_length", self._filter_min_avg_repeat_length),
            ("maximum_mean_repeat_length", self._filter_max_avg_repeat_length),
            ("maximum_mean_spacer_length", self._filter_max_avg_spacer_length),
            ("minimum_mean_spacer_length", self._filter_by_min_avg_spacer),
            ("maximum_long_spacer_burden", self._filter_by_max_spacer),
            ("identical_spacer_burden", self._filter_by_the_same_spacer),
            ("repeat_terminal_columns", self._filter_by_column),
            ("overall_repeat_consistency", self._filter_by_overall_repeat_consistency),
        )
        results = [
            (name, "passed" if predicate(candidate) else "failed")
            for name, predicate in predicates
        ]
        results.insert(7, (
            "spacer_terminal_dominance",
            terminal.status,
        ))
        passes = (
            terminal.passes_filter
            and all(status == "passed" for name, status in results
                    if name != "spacer_terminal_dominance")
        )
        return AdvancedFuzzySearchFilterEvidence(
            decision="passed" if passes else "rejected",
            structural_valid=True,
            predicate_results=tuple(results),
            terminal_dominance=terminal,
        )

    def __call__(self, candidate):
        # These predicates form a pure conjunction. Put constant-time length
        # and count gates first so candidates that are structurally impossible
        # never pay for gapped-column scans.
        if not self._candidate_structure_is_valid(candidate):
            return
        if not self._passes_non_gapped_predicates(candidate):
            return
        if not self._filter_by_column(candidate):
            return
        if not self._filter_by_overall_repeat_consistency(candidate):
            return
        return candidate

#              CRISPR Candidate
#####################################################
#####################################################
class CrisprConsensus(object):
    def __init__(self, list_repeats_gaped):
        self.list_repeats_gaped = list_repeats_gaped

        self.num_different_repeat_length = None
        self.consensus = None
        self.consensus_no_gap = None
        self.len_consensus = None
        self.number_repeats = None

        self._check_repeat_length()
        self._compute_consensus()

    def _check_repeat_length(self):
        list_lengths = [len(repeat) for repeat in self.list_repeats_gaped]
        self.num_different_repeat_length = len(set(list_lengths))

    def _compute_consensus(self):
        if self.num_different_repeat_length == 0:
            print('Got repeats of 0 length')
        elif self.num_different_repeat_length != 1:
            print('Got a case with different repeat lengths')
            for rep_gapped in self.list_repeats_gaped:
                print(rep_gapped)
        else:
            self.consensus = ''
            for char_ind, _ in enumerate(self.list_repeats_gaped[0]):
                list_char_in_column = [repeat[char_ind] for repeat in self.list_repeats_gaped]
                counter = collections.Counter(list_char_in_column)
                freq = counter.most_common()
                most_common_char = freq[0][0] if freq[0][0] != '-' else freq[1][0]
                self.consensus += most_common_char

        self.consensus_no_gap = self.consensus.replace(' ', '').replace('+', '')
        self.len_consensus = len(self.consensus_no_gap)

    def output(self):
        return self.consensus_no_gap, self.consensus


class CrisprCandidate(object):
    def __init__(self, list_repeats, list_repeats_gaped, list_spacers, list_repeat_starts):
        self.list_repeats = list_repeats
        self.list_repeats_gaped = list_repeats_gaped
        self.list_spacers = list_spacers
        self.list_repeat_starts = list_repeat_starts

        self._validate_alignment_structure()

        self.list_repeat_mismatches = []
        self.list_mismatches_indexes = []

        self.consensus = None
        self.consensus_gaped = None
        self.total_mismatches = None

        self._filter_redundant_insertion_deletions()
        self._compute_consensus()
        self._compute_mismatches()

        self.list_gaped_repeats = self.list_repeats_gaped
        mark_ordinary_birth(self)

    def _validate_alignment_structure(self):
        repeat_count = len(self.list_repeats)
        if repeat_count == 0:
            raise ValueError("a CRISPR candidate must contain repeats")
        if len(self.list_repeats_gaped) != repeat_count:
            raise ValueError("every repeat must have one gapped representation")
        if len(self.list_repeat_starts) != repeat_count:
            raise ValueError("every repeat must have one start coordinate")
        if len(self.list_spacers) != repeat_count - 1:
            raise ValueError("spacer count must be one less than repeat count")
        if any(not isinstance(repeat, str) or not repeat
               for repeat in self.list_repeats):
            raise ValueError("repeat sequences must be non-empty strings")
        if any(not isinstance(spacer, str) for spacer in self.list_spacers):
            raise ValueError("spacer sequences must be strings")
        if any(not isinstance(repeat, str) or not repeat
               for repeat in self.list_repeats_gaped):
            raise ValueError("gapped repeats must be non-empty strings")

        gapped_lengths = {len(repeat) for repeat in self.list_repeats_gaped}
        if len(gapped_lengths) != 1:
            raise ValueError("gapped repeats must have equal alignment lengths")
        for index, (repeat, gapped_repeat) in enumerate(zip(
                self.list_repeats, self.list_repeats_gaped)):
            reconstructed = gapped_repeat.replace(" ", "").replace("-", "")
            if reconstructed != repeat:
                raise ValueError(
                    "gapped repeat {} does not reconstruct its source repeat".format(
                        index
                    )
                )

    def _filter_redundant_insertion_deletions(self):
        def _fix_repeats(list_repeats, list_bad_indexes_to_fix):
            list_repeats_new = []
            for repeat in list_repeats:
                list_repeats_new.append(_fix_repeat(repeat, list_bad_indexes_to_fix))

            return list_repeats_new

        def _fix_repeat(repeat, list_bad_indexes_to_fix):
            new_repeat = ''
            for index, char in enumerate(repeat):
                if index not in list_bad_indexes_to_fix:
                    new_repeat += char

            return new_repeat

        list_bad_indexes = []
        for char_ind, _ in enumerate(self.list_repeats_gaped[0]):
            list_char_in_column = [repeat[char_ind] for repeat in self.list_repeats_gaped]
            chars = set(list_char_in_column)

            if chars == {' '} or chars == {'-'}:
                list_bad_indexes.append(char_ind)

        if list_bad_indexes:
            self.list_repeats_gaped = _fix_repeats(self.list_repeats_gaped, list_bad_indexes)

    def _compute_consensus(self):
        self.consensus, self.consensus_gaped = CrisprConsensus(self.list_repeats_gaped).output()

    def _compute_mismatches(self):
        def _compute_mismatches_repeat(gaped_repeat):
            substitutions = 0
            insertions = 0
            deletions = 0
            list_mismatches_indexes_one_repeat = []
            for index, char_repeat, char_con_repeat in zip(range(len(gaped_repeat)),
                                                           gaped_repeat,
                                                           self.consensus_gaped):

                if char_con_repeat == ' ':
                    if char_repeat != ' ':
                        insertions += 1
                        list_mismatches_indexes_one_repeat.append(index)
                else:
                    if char_repeat == char_con_repeat:
                        pass
                    else:
                        if char_repeat == '-':
                            deletions += 1
                            list_mismatches_indexes_one_repeat.append(index)
                        elif char_repeat == ' ':
                            deletions += 1
                        else:
                            substitutions += 1
                            list_mismatches_indexes_one_repeat.append(index)

            return substitutions, insertions, deletions, list_mismatches_indexes_one_repeat

        for gaped_repeat in self.list_repeats_gaped:
            s, i, d, list_mismatches_indexes_one_repeat = _compute_mismatches_repeat(gaped_repeat)
            total = s + i + d
            repeat_stats = [s, i, d, total]
            self.list_repeat_mismatches.append(repeat_stats)
            self.list_mismatches_indexes.append(list_mismatches_indexes_one_repeat)

        self.total_mismatches = sum([x[3] for x in self.list_repeat_mismatches])

    def dot_repeat(self, gaped_repeat):
        string = ''
        substitutions = 0
        insertions = 0
        deletions = 0
        for char_repeat, char_consensus in zip(gaped_repeat, self.consensus_gaped):
            if char_consensus == ' ':
                string += char_repeat
                if char_repeat != ' ':
                    insertions += 1
            else:
                if char_repeat == char_consensus:
                    string += '.'
                else:
                    string += char_repeat
                    if char_repeat == '-':
                        deletions += 1
                    elif char_repeat == ' ':
                        deletions += 1
                    else:
                        substitutions += 1
        return string, substitutions, insertions, deletions

    def dot_repr(self):
        string = ''
        g_s, g_i, g_d = 0, 0, 0
        max_length_start_index = max(len(str(start)) for start in self.list_repeat_starts) + 3
        max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3

        for index, gaped_repeat in enumerate(self.list_repeats_gaped):
            repeat_start_index = self.list_repeat_starts[index] + 1
            n_gaps_after_start = max_length_start_index - len(str(repeat_start_index))

            if index == len(self.list_spacers):
                spacer = ""
            else:
                spacer = self.list_spacers[index]
            n_gaps_after_spacer = max_length_spacer - len(spacer)

            dotted_repeats, s, i, d = self.dot_repeat(gaped_repeat)
            errors = "   s:{} i:{} d:{}".format(s, i, d)
            g_s += s
            g_i += i
            g_d += d

            string += "{}{}{}  {}{}{}\n".format(repeat_start_index,
                                                " " * n_gaps_after_start,
                                                dotted_repeats, spacer,
                                                " " * n_gaps_after_spacer,
                                                errors)

        string += "_" * 100 + "\n"

        string += " " * max_length_start_index + self.consensus_gaped
        string += " " * (max_length_spacer + 2) + "   s:{} i:{} d:{}".format(g_s, g_i, g_d) + "\n"

        return string

    def dot_repr_web_server(self):
        string = ''
        g_s, g_i, g_d = 0, 0, 0
        max_length_start_index = max(len(str(start)) for start in self.list_repeat_starts) + 3
        max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3

        for index, gaped_repeat in enumerate(self.list_repeats_gaped):
            repeat_start_index = self.list_repeat_starts[index] + 1
            n_gaps_after_start = max_length_start_index - len(str(repeat_start_index))

            if index == len(self.list_spacers):
                spacer = ""
            else:
                spacer = "$" + self.list_spacers[index] + "$"
            n_gaps_after_spacer = max_length_spacer - len(spacer)

            dotted_repeats, s, i, d = self.dot_repeat(gaped_repeat)
            errors = "   s:{} i:{} d:{}".format(s, i, d)
            g_s += s
            g_i += i
            g_d += d

            string += "{}{}{}  {}{}{}\n".format(repeat_start_index,
                                                " " * n_gaps_after_start,
                                                dotted_repeats, spacer,
                                                " " * n_gaps_after_spacer,
                                                errors)

        string += "_" * 100 + "\n"

        string += " " * max_length_start_index + self.consensus_gaped
        string += " " * (max_length_spacer + 2) + "   s:{} i:{} d:{}".format(g_s, g_i, g_d) + "\n"

        string += "_" * 100 + "\n"

        string += "consensus: " + self.consensus + "\n"

        return string

    def write_file(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.dot_repr())

    def write_as_json(self, filename):
        dict_to_write = {"repeat_begins": self.list_repeat_starts,
                         "repeats": self.list_repeats,
                         "repeats_gaped": self.list_repeats_gaped,
                         "spacers": self.list_spacers}

        with open(filename, 'w') as outfile:
            json.dump(dict_to_write, outfile)

    def compute_stats(self):
        start = self.list_repeat_starts[0] + 1
        end = self.list_repeat_starts[-1] + len(self.list_repeats[-1])
        avg_repeat = len(self.consensus)
        avg_spacer = int(sum((len(spacer) for spacer in self.list_spacers)) / len(self.list_spacers))
        number_repeats = len(self.list_repeats)
        return {"start": start, "end": end, "avg_repeat": avg_repeat,
                "avg_spacer": avg_spacer, "number_repeats": number_repeats}

    @classmethod
    def init_from_json(cls, file_name):
        with open(file_name) as json_file:
            dict_data = json.load(json_file)

            list_repeas = dict_data["repeats"]
            list_repeats_starts = dict_data["repeat_begins"]
            list_spacers = dict_data["spacers"]
            list_repeats_gaped = dict_data["repeats_gaped"]

        return cls(list_repeats=list_repeas, list_spacers=list_spacers,
                   list_repeats_gaped=list_repeats_gaped, list_repeat_starts=list_repeats_starts)

    def __repr__(self):
        return self.dot_repr()

    def __eq__(self, other):
        if not isinstance(other, CrisprCandidate):
            return NotImplemented
        return (
            self.list_repeats == other.list_repeats
            and self.list_repeats_gaped == other.list_repeats_gaped
            and self.list_spacers == other.list_spacers
            and self.list_repeat_starts == other.list_repeat_starts
        )

    def __ne__(self, other):
        return not self.__eq__(other)
