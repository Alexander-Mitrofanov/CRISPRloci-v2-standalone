import os
import subprocess
import hashlib
from collections import namedtuple
from os.path import basename
from os import listdir
from os.path import isfile
import regex

from components.fuzzy_alignment import normalize_fuzzy_changes
from components.fuzzy_alignment_merge import build_pairwise_fuzzy_alignment
from components.fuzzy_alignment_merge import merge_fuzzy_alignments
from components.fuzzy_alignment_merge import merge_pairwise_fuzzy_alignments
from components.fuzzy_boundary_lattice import (
    FUZZY_BOUNDARY_LATTICE_SCHEMA,
    FuzzyBoundaryLattice,
    FuzzyBoundarySearchPolicy,
    bind_fuzzy_boundary_source,
    build_fuzzy_boundary_alternative_plan,
)
from components.fuzzy_occurrence_paths import OCCURRENCE_PATH_SCHEMA
from components.fuzzy_occurrence_paths import OccurrencePathLattice
from components.fuzzy_occurrence_paths import build_occurrence_path_plan
from components.fuzzy_lattice_storage import CompressedLattice
from components.reconstruction_validator import read_single_fasta
from components.runtime_resources import external_thread_budget


ClusterSequence = namedtuple("ClusterSequence", "sequence, start, end, tuple_repeats")

# Keep the local-chain contract aligned with the individual-spacer ceiling
# used by ``AdvancedFuzzySearchFilter``. A larger source gap cannot be a
# direct repeat-to-repeat edge in one local array candidate; it starts a new
# maximal chain instead.
MAX_LOCAL_CHAIN_SPACER_LENGTH = 140


FuzzySourceOccurrence = namedtuple(
    "FuzzySourceOccurrence",
    "occurrence_id start end sequence",
)
class FuzzySourcePathProvenance(namedtuple(
        "FuzzySourcePathProvenance",
        (
            "schema family_id repeat_candidate parent_geometry "
            "parent_spacers occurrence_ids derivation source_index_span"
        ))):
    """Legacy fuzzy path with an explicit baseline origin label."""

    __slots__ = ()
    occurrence_path_kind = "baseline_fuzzy"

    def __eq__(self, other):
        if isinstance(other, FuzzySourcePathProvenance):
            return (
                self.occurrence_path_kind == other.occurrence_path_kind
                and tuple.__eq__(self, other)
            )
        return tuple.__eq__(self, other)

    def __hash__(self):
        return hash(tuple(self))


class FuzzyExactOverlaySourcePathProvenance(FuzzySourcePathProvenance):
    """Source-path record derived from a complete exact overlay lattice."""

    __slots__ = ()
    occurrence_path_kind = "source_exact_overlay"


class FuzzyBoundaryAlternativeSourcePathProvenance(
        FuzzySourcePathProvenance):
    """Source path selected from a complete forced-boundary lattice."""

    __slots__ = ()
    occurrence_path_kind = "forced_boundary_alternative"


FUZZY_SOURCE_PATH_SCHEMA = "crispridentify.fuzzy-source-path.v1"


class RejectedFuzzyCandidateWitness(namedtuple(
        "RejectedFuzzyCandidateWitness",
        (
            "start end same_start_end_key "
            "source_path_family_parent number_errors "
            "two_repeat_terminal_repair_seed"
        ))):
    """Compact accepted-only witness for a structurally impossible path.

    Rejected fuzzy paths still participate in the historical transitive
    source-locus partition.  Keeping only their half-open envelope and exact
    SameStartEnd alignment class preserves those grouping semantics without
    retaining their gapped alignment, provenance graph, or match objects.
    """

    __slots__ = ()

    def __new__(cls, start, end, same_start_end_key,
                source_path_family_parent, number_errors,
                two_repeat_terminal_repair_seed=None):
        if (
                not isinstance(start, int)
                or isinstance(start, bool)
                or not isinstance(end, int)
                or isinstance(end, bool)
                or start < 0
                or end <= start):
            raise ValueError("rejected fuzzy witness has an invalid interval")
        if not isinstance(same_start_end_key, tuple):
            raise TypeError(
                "rejected fuzzy witness alignment key must be a tuple"
            )
        try:
            hash(same_start_end_key)
        except TypeError as error:
            raise ValueError(
                "rejected fuzzy witness alignment key must be immutable"
            ) from error
        if not isinstance(source_path_family_parent, bool):
            raise TypeError(
                "rejected fuzzy witness source-path role must be boolean"
            )
        if (
                not isinstance(number_errors, int)
                or isinstance(number_errors, bool)
                or number_errors < 0):
            raise ValueError(
                "rejected fuzzy witness error count must be non-negative"
            )
        if (
                two_repeat_terminal_repair_seed is not None
                and not isinstance(
                    two_repeat_terminal_repair_seed,
                    TwoRepeatTerminalRepairSeed,
                )):
            raise TypeError(
                "rejected fuzzy witness repair seed has an invalid type"
            )
        return super(RejectedFuzzyCandidateWitness, cls).__new__(
            cls,
            start,
            end,
            same_start_end_key,
            source_path_family_parent,
            number_errors,
            two_repeat_terminal_repair_seed,
        )


class TwoRepeatTerminalRepairSeed(namedtuple(
        "TwoRepeatTerminalRepairSeed",
        "repeat_starts repeat_ends cluster_index")):
    """Compact source geometry for one opt-in two-repeat repair attempt."""

    __slots__ = ()

    def __new__(cls, repeat_starts, repeat_ends, cluster_index=None):
        repeat_starts = tuple(repeat_starts)
        repeat_ends = tuple(repeat_ends)
        if len(repeat_starts) != 2 or len(repeat_ends) != 2:
            raise ValueError("two-repeat repair seed must contain two repeats")
        if any(
                not isinstance(value, int) or isinstance(value, bool)
                for value in repeat_starts + repeat_ends):
            raise TypeError("two-repeat repair boundaries must be integers")
        if not (
                0 <= repeat_starts[0] < repeat_ends[0]
                <= repeat_starts[1] < repeat_ends[1]):
            raise ValueError("two-repeat repair seed geometry is invalid")
        if (
                cluster_index is not None
                and (
                    not isinstance(cluster_index, int)
                    or isinstance(cluster_index, bool)
                    or cluster_index < 0
                )):
            raise ValueError(
                "two-repeat repair cluster index must be non-negative"
            )
        return super(TwoRepeatTerminalRepairSeed, cls).__new__(
            cls, repeat_starts, repeat_ends, cluster_index
        )

    @property
    def source_interval(self):
        return self.repeat_starts[0], self.repeat_ends[-1]


def _strictly_merge_intervals(intervals):
    """Return sorted half-open components under strict source overlap."""

    ordered = sorted(intervals)
    if not ordered:
        return ()
    merged = []
    current_start, current_end = ordered[0]
    for start, end in ordered[1:]:
        if start < current_end:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    return tuple(merged)


class RejectedFuzzyCandidateTopology(namedtuple(
        "RejectedFuzzyCandidateTopology",
        (
            "local_intervals family_parent_intervals "
            "raw_rejected_candidates reduced_rejected_candidates "
            "two_repeat_terminal_repair_seeds"
        ))):
    """Exact rejected-path quotient used by accepted-only native detection.

    SameStartEnd equivalence is resolved inside one repeat-pattern worker
    before this marker is constructed.  Downstream rejected paths are observed
    only through strict-overlap locus connectivity and their parent/local role,
    so one interval hull per connected component is sufficient evidence.
    """

    __slots__ = ()

    def __new__(cls, local_intervals=(), family_parent_intervals=(),
                raw_rejected_candidates=0, reduced_rejected_candidates=0,
                two_repeat_terminal_repair_seeds=()):
        normalized_roles = []
        for label, intervals in (
                ("local", local_intervals),
                ("family-parent", family_parent_intervals)):
            normalized = []
            for interval in intervals:
                if (
                        not isinstance(interval, tuple)
                        or len(interval) != 2):
                    raise ValueError(
                        "{} rejected topology interval is malformed".format(
                            label
                        )
                    )
                start, end = interval
                if (
                        not isinstance(start, int)
                        or isinstance(start, bool)
                        or not isinstance(end, int)
                        or isinstance(end, bool)
                        or start < 0
                        or end <= start):
                    raise ValueError(
                        "{} rejected topology interval is invalid".format(
                            label
                        )
                    )
                normalized.append((start, end))
            normalized_roles.append(_strictly_merge_intervals(normalized))
        for label, value in (
                ("raw", raw_rejected_candidates),
                ("reduced", reduced_rejected_candidates)):
            if (
                    not isinstance(value, int)
                    or isinstance(value, bool)
                    or value < 0):
                raise ValueError(
                    "{} rejected topology count must be non-negative".format(
                        label
                    )
                )
        if reduced_rejected_candidates > raw_rejected_candidates:
            raise ValueError(
                "rejected topology reduction cannot increase candidates"
            )
        if reduced_rejected_candidates and not any(normalized_roles):
            raise ValueError("rejected topology lost every source interval")
        normalized_seeds = tuple(sorted(set(
            two_repeat_terminal_repair_seeds
        )))
        if any(
                not isinstance(seed, TwoRepeatTerminalRepairSeed)
                for seed in normalized_seeds):
            raise TypeError("rejected topology repair seed has an invalid type")
        return super(RejectedFuzzyCandidateTopology, cls).__new__(
            cls,
            normalized_roles[0],
            normalized_roles[1],
            raw_rejected_candidates,
            reduced_rejected_candidates,
            normalized_seeds,
        )

    @property
    def interval_runs(self):
        return len(self.local_intervals) + len(self.family_parent_intervals)

    @classmethod
    def combine(cls, markers):
        markers = tuple(markers)
        return cls(
            tuple(
                interval
                for marker in markers
                for interval in marker.local_intervals
            ),
            tuple(
                interval
                for marker in markers
                for interval in marker.family_parent_intervals
            ),
            sum(marker.raw_rejected_candidates for marker in markers),
            sum(marker.reduced_rejected_candidates for marker in markers),
            tuple(
                seed
                for marker in markers
                for seed in marker.two_repeat_terminal_repair_seeds
            ),
        )


_ValidatedFuzzyAlignmentToken = namedtuple(
    "_ValidatedFuzzyAlignmentToken",
    "normalized_changes pairwise_alignment class_id",
)


def _freeze_fuzzy_match_state(value):
    """Take a mutation-sensitive snapshot of nested matcher edit state."""

    if isinstance(value, tuple):
        return ("tuple", tuple(_freeze_fuzzy_match_state(item)
                               for item in value))
    if isinstance(value, list):
        return ("list", tuple(_freeze_fuzzy_match_state(item)
                              for item in value))
    if isinstance(value, (set, frozenset)):
        return (
            "set",
            frozenset(_freeze_fuzzy_match_state(item) for item in value),
        )
    if isinstance(value, dict):
        return (
            "dict",
            frozenset(
                (_freeze_fuzzy_match_state(key),
                 _freeze_fuzzy_match_state(item))
                for key, item in value.items()
            ),
        )
    try:
        hash(value)
    except TypeError:
        return ("opaque", type(value).__name__, repr(value))
    return ("atom", type(value), value)


class _FuzzyAlignmentTokenCache:
    """Worker-local normalization/validation cache shared by sibling paths."""

    def __init__(self):
        self._by_match_identity = {}
        self._class_ids = {}
        self._tokens_by_logical_class = {}

    @staticmethod
    def _snapshot(match):
        snapshot_hook = getattr(
            match, "_fuzzy_alignment_cache_snapshot", None
        )
        direct_state = snapshot_hook() if callable(snapshot_hook) else None
        if direct_state is None:
            counts_state = _freeze_fuzzy_match_state(match.fuzzy_counts)
            changes_state = _freeze_fuzzy_match_state(match.fuzzy_changes)
        else:
            counts_state, changes_state = direct_state
        return (
            match.start(),
            match.end(),
            match.group(),
            counts_state,
            changes_state,
            getattr(match, "deletion_position_shift", None),
        )

    def token_for(self, repeat_candidate, sequence_start, match,
                  expected_absolute_start, expected_repeat):
        snapshot = self._snapshot(match)
        identity = id(match)
        cached = self._by_match_identity.get(identity)
        if (
                cached is not None
                and cached[0] is match
                and cached[1] == snapshot):
            return cached[2]

        absolute_start = sequence_start + match.start()
        if (
                absolute_start != expected_absolute_start
                or match.group() != expected_repeat
                or match.end() - match.start() != len(expected_repeat)):
            raise ValueError(
                "deferred fuzzy match no longer reconstructs its candidate"
            )
        normalized = normalize_fuzzy_changes(match)
        canonical = tuple(
            tuple(sorted(indexes)) for indexes in normalized
        )
        logical_class = (
            repeat_candidate,
            expected_absolute_start,
            expected_repeat,
            canonical,
        )
        validated = self._tokens_by_logical_class.get(logical_class)
        if validated is None:
            pairwise = build_pairwise_fuzzy_alignment(
                repeat_candidate, expected_repeat, canonical
            )
            class_id = self._class_ids.setdefault(
                logical_class, len(self._class_ids)
            )
            self._tokens_by_logical_class[logical_class] = (
                pairwise, class_id
            )
        else:
            pairwise, class_id = validated
        token = _ValidatedFuzzyAlignmentToken(
            canonical, pairwise, class_id
        )
        # Retain the match strongly alongside its id so id reuse cannot alias
        # an unrelated matcher object during a long worker batch.
        self._by_match_identity[identity] = (match, snapshot, token)
        return token


def _within_one_edit(first, second):
    """Return whether two strings are separated by at most one edit.

    Vmatch is invoked with ``-e 1``. Rechecking that contract against the two
    original source slices makes both reported occurrence indexes independently
    testable instead of trusting the external program coordinate columns.
    """

    if abs(len(first) - len(second)) > 1:
        return False
    if len(first) == len(second):
        return sum(left != right for left, right in zip(first, second)) <= 1

    shorter, longer = (
        (first, second) if len(first) < len(second) else (second, first)
    )
    short_index = 0
    long_index = 0
    skipped = False
    while short_index < len(shorter) and long_index < len(longer):
        if shorter[short_index] == longer[long_index]:
            short_index += 1
            long_index += 1
        elif skipped:
            return False
        else:
            skipped = True
            long_index += 1
    return True


class V_Repeat(namedtuple(
        "V_Repeat",
        "begin_first begin_second length sequence length_first length_second")):
    """One Vmatch pair with separate source extents for both occurrences.

    ``length`` remains the length of ``sequence`` for compatibility with the
    fuzzy-search seed.  Vmatch can report different aligned lengths for the
    two genomic occurrences, so cluster bounds must use ``length_first`` and
    ``length_second`` rather than pretending both spans have ``length`` bases.
    """

    __slots__ = ()

    def __new__(cls, begin_first, begin_second, length, sequence,
                length_first=None, length_second=None):
        for label, value in (
                ("first occurrence start", begin_first),
                ("second occurrence start", begin_second),
                ("seed length", length)):
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError("{} must be an integer".format(label))
        if begin_first < 0 or begin_second < 0:
            raise ValueError("Vmatch occurrence starts must be non-negative")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError("Vmatch seed sequence must be a non-empty string")
        if length != len(sequence):
            raise ValueError("Vmatch seed length does not match its sequence")
        if length_first is None:
            length_first = length
        if length_second is None:
            length_second = length
        for label, value in (
                ("first occurrence length", length_first),
                ("second occurrence length", length_second)):
            if (
                    not isinstance(value, int)
                    or isinstance(value, bool)
                    or value <= 0):
                raise ValueError("{} must be a positive integer".format(label))
        return super(V_Repeat, cls).__new__(
            cls,
            begin_first,
            begin_second,
            length,
            sequence,
            length_first,
            length_second,
        )


class V_RepeatOccurrence(namedtuple(
        "V_RepeatOccurrence",
        "begin length sequence source_sequence v_repeat side")):
    """One exact source occurrence carrying its Vmatch seed provenance."""

    __slots__ = ()

    def __new__(cls, begin, length, sequence, source_sequence, v_repeat, side):
        if (
                not isinstance(begin, int)
                or isinstance(begin, bool)
                or begin < 0):
            raise ValueError("Vmatch occurrence start must be non-negative")
        if (
                not isinstance(length, int)
                or isinstance(length, bool)
                or length <= 0):
            raise ValueError("Vmatch occurrence length must be positive")
        if not isinstance(sequence, str) or not sequence:
            raise ValueError("Vmatch occurrence seed must be non-empty")
        if source_sequence is not None and (
                not isinstance(source_sequence, str)
                or len(source_sequence) != length):
            raise ValueError(
                "Vmatch occurrence source slice must match its exact length"
            )
        if not isinstance(v_repeat, V_Repeat):
            raise ValueError("Vmatch occurrence must retain its pair record")
        if side not in (0, 1):
            raise ValueError("Vmatch occurrence side must be 0 or 1")
        return super(V_RepeatOccurrence, cls).__new__(
            cls, begin, length, sequence, source_sequence, v_repeat, side
        )

    @property
    def end(self):
        return self.begin + self.length


def _v_repeat_sort_key(repeat):
    return (
        min(repeat.begin_first, repeat.begin_second),
        max(repeat.begin_first, repeat.begin_second),
        repeat.length,
        repeat.length_first,
        repeat.length_second,
        repeat.sequence,
        repeat.begin_first,
        repeat.begin_second,
    )


def _occurrence_sort_key(occurrence):
    return (
        occurrence.begin,
        occurrence.end,
        occurrence.sequence,
        occurrence.side,
        _v_repeat_sort_key(occurrence.v_repeat),
    )


def _occurrences_for_v_repeat(repeat, source_dna=None):
    occurrences = []
    for side, (begin, length) in enumerate((
            (repeat.begin_first, repeat.length_first),
            (repeat.begin_second, repeat.length_second))):
        end = begin + length
        source_sequence = None
        occurrence_seed = repeat.sequence
        if source_dna is not None:
            if end > len(source_dna):
                raise ValueError(
                    "Vmatch occurrence is outside the original source DNA"
                )
            source_sequence = source_dna[begin:end]
            if len(source_sequence) != length:
                raise ValueError(
                    "Vmatch occurrence does not reconstruct its source span"
                )
            # Vmatch reports one pair-level seed, historically copied onto
            # both occurrences.  The two aligned source spans can differ by
            # the permitted edit, so derive each fuzzy-search seed from its
            # own source side.  Keep the existing pair seed length as a hard
            # ceiling; a deletion-shortened side contributes its full span.
            occurrence_seed = source_sequence[:min(repeat.length, length)]
            if (
                    not occurrence_seed
                    or source_dna[
                        begin:begin + len(occurrence_seed)
                    ] != occurrence_seed):
                raise ValueError(
                    "Vmatch side seed does not reconstruct its source prefix"
                )
        occurrences.append(V_RepeatOccurrence(
            begin,
            length,
            occurrence_seed,
            source_sequence,
            repeat,
            side,
        ))
    return occurrences


def _deduplicate_occurrences(occurrences):
    """Collapse repeated Vmatch edges without discarding distinct seeds."""

    unique = {}
    for occurrence in sorted(occurrences, key=_occurrence_sort_key):
        key = (
            occurrence.begin,
            occurrence.length,
            occurrence.sequence,
        )
        previous = unique.get(key)
        if previous is not None:
            if (
                    previous.source_sequence is not None
                    and occurrence.source_sequence is not None
                    and previous.source_sequence != occurrence.source_sequence):
                raise ValueError(
                    "duplicate Vmatch occurrence has conflicting source slices"
                )
            if (
                    previous.source_sequence is None
                    and occurrence.source_sequence is not None):
                unique[key] = occurrence
        else:
            unique[key] = occurrence
    return sorted(unique.values(), key=_occurrence_sort_key)

#                    INITIAL REPEAT CANDIDATES
#########################################################################
#########################################################################


class VmatchRun:
    # Ask Vmatch for one sentinel record above the accepted domain.  This
    # makes the historical million-record ceiling observable instead of
    # silently importing a truncated candidate set.
    RESULT_LIMIT = 1_000_000
    FAST_SEED_SETTINGS = {
        "balanced": "23 -10 100",
        "legacy": "23 25 60",
        "sensitive": "21 -10 100",
    }

    def __init__(
            self, file_path, flag_fast_run, fast_seed_profile="legacy"):
        self.file_path = file_path
        self.flag_fast_run = flag_fast_run
        self.fast_seed_profile = fast_seed_profile
        self.list_repeat_candidates = []

        self._get_run_settings()
        self._create_folders()
        self._create_input_uppercase()
        self._mkvtree_command()
        self._vmatch_command()
        self._take_vmatch_results_new()
        self._clean_after_mkv()

    def _get_run_settings(self):
        if self.fast_seed_profile not in self.FAST_SEED_SETTINGS:
            raise ValueError(
                "fast seed profile must be balanced, legacy, or sensitive"
            )
        if self.flag_fast_run:
            self.settings = {"vmatch -l": self.FAST_SEED_SETTINGS[
                                 self.fast_seed_profile
                             ],
                             "vmatch -e": "1",
                             "vmatch -evalue": "1",
                             "flanking_region": "100",
                             "max_length_between_repeats": "100"}
        else:
            self.settings = {"vmatch -l": "21 -10 100",
                             "vmatch -e": "1",
                             "vmatch -evalue": "1",
                             "flanking_region": "100",
                             "max_length_between_repeats": "100"}

    def _create_folders(self):
        self.file_base = basename(self.file_path)

    def _create_input_uppercase(self):
        header, self.dna = read_single_fasta(self.file_path)
        self.input_header = ">{}\n".format(header)

        with open("new_input.fa", "w") as f:
            f.write(self.input_header)
            f.write(self.dna)

    def _get_double_dna(self):
        header, self.dna = read_single_fasta(self.file_path)
        self.input_header = ">{}\n".format(header)
        self.dna_length = len(self.dna)
        self.dna += self.dna

    def _mkvtree_command(self):
        command = [
            "mkvtree", "-db", "new_input.fa", "-dna", "-pl", "-lcp",
            # ``vmatch -s leftseq`` requires these five tables.  ``bck`` and
            # ``sti1`` are not consulted for this query shape; omitting them
            # reduced genome-8 temporary index writes by 13.6% with identical
            # ordered Vmatch output and coordinates.
            "-suf", "-tis", "-ois", "-bwt",
        ]
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def _vmatch_command(self):
        l_flag_val = self.settings["vmatch -l"]
        e_flag_val = self.settings["vmatch -e"]
        e_value_flag_val = self.settings["vmatch -evalue"]

        command = (
            ["vmatch", "-l"]
            + l_flag_val.split()
            + [
                "-evalue", e_value_flag_val,
                "-e", e_flag_val,
                "-s", "leftseq",
                "-absolute",
                "-nodist",
                "-noevalue",
                "-noscore",
                "-noidentity",
                "-sort", "ia",
                "-best", str(self.RESULT_LIMIT + 1),
                "new_input.fa",
            ]
        )
        with open("vmatch_result_new.txt", "w") as output_handle:
            subprocess.run(
                command,
                check=True,
                stdout=output_handle,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )


    def _take_vmatch_results_new(self):
        cutoff = 55
        file_to_open = "vmatch_result_new.txt"
        with open(file_to_open, "r") as stream:
            record_count = sum(
                1 for line in stream if line.startswith(">")
            )
        if record_count > self.RESULT_LIMIT:
            raise RuntimeError(
                "Vmatch candidate domain exceeded {:,} records; refusing "
                "silent -best truncation".format(self.RESULT_LIMIT)
            )

        def records(stream):
            header = None
            sequence_parts = []
            for raw_line in stream:
                if raw_line.startswith(">"):
                    if header is not None:
                        yield header, "".join(sequence_parts)
                    header = raw_line
                    sequence_parts = []
                elif header is not None and raw_line.strip():
                    sequence_parts.append(raw_line.strip())
            if header is not None:
                yield header, "".join(sequence_parts)

        with open(file_to_open, "r") as stream:
            for line, reported_sequence in records(stream):
                # Vmatch wraps `-s leftseq` output at 60 columns. Joining the
                # complete FASTA-like record is required before checking the
                # occurrence length or truncating the fuzzy-search seed.
                line = line.replace(">", "> ")
                line = line.replace(" +", " ")
                list_param = line.split()
                if len(list_param) < 6:
                    raise ValueError(
                        "malformed Vmatch coordinate header: {}".format(
                            line.strip()
                        )
                    )
                raw_length_first = int(list_param[1])
                raw_length_second = int(list_param[4])
                if len(reported_sequence) != raw_length_first:
                    raise ValueError(
                        "Vmatch left sequence length does not match its header"
                    )
                begin_first = int(list_param[2])
                begin_second = int(list_param[5])
                for begin, length in (
                        (begin_first, raw_length_first),
                        (begin_second, raw_length_second)):
                    if begin < 0 or begin + length > len(self.dna):
                        raise ValueError(
                            "Vmatch occurrence lies outside the source sequence"
                        )
                first_source = self.dna[
                    begin_first:begin_first + raw_length_first
                ].upper()
                second_source = self.dna[
                    begin_second:begin_second + raw_length_second
                ].upper()
                if first_source != reported_sequence.upper():
                    raise ValueError(
                        "Vmatch first occurrence does not match its complete "
                        "source slice"
                    )
                if not _within_one_edit(first_source, second_source):
                    raise ValueError(
                        "Vmatch occurrence coordinates do not reconstruct a "
                        "pair within the configured edit distance"
                    )

                seed_sequence = reported_sequence[:cutoff]
                seed_length = len(seed_sequence)
                repeat = V_Repeat(
                    begin_first,
                    begin_second,
                    seed_length,
                    seed_sequence,
                    raw_length_first,
                    raw_length_second,
                )
                self.list_repeat_candidates.append(repeat)

    def _clean_after_mkv(self):
        """Cleans the folder moving the intermediate files into
        the corresponding temp folder"""

        onlyfiles = [f for f in listdir('.') if isfile(f)]
        to_move = []

        for el in onlyfiles:
            if "new_input.fa" in el:
                to_move.append(el)

        for el in to_move:
            os.remove(el)
        try:
            os.remove("vmatch_result_new.txt")
        except Exception:
            pass

    def output(self):
        return self.list_repeat_candidates


#                       CLUSTERS
###########################################################################
###########################################################################


class ClusterCandidate:
    def __init__(self, list_cluster_v_repeats, list_cluster_occurrences=None):
        self.list_cluster_repeats = sorted(
            list_cluster_v_repeats,
            key=_v_repeat_sort_key,
        )
        if list_cluster_occurrences is None:
            list_cluster_occurrences = [
                occurrence
                for repeat in self.list_cluster_repeats
                for occurrence in _occurrences_for_v_repeat(repeat)
            ]
        self.list_cluster_occurrences = _deduplicate_occurrences(
            list_cluster_occurrences
        )
        if not self.list_cluster_occurrences:
            raise ValueError("a cluster must contain a Vmatch occurrence")

        self.begin = min(
            occurrence.begin for occurrence in self.list_cluster_occurrences
        )
        self.end = max(
            occurrence.end for occurrence in self.list_cluster_occurrences
        )
        self.goes_through_zero = None

        self.list_clust_dif_rep_seq = sorted({
            occurrence.sequence
            for occurrence in self.list_cluster_occurrences
        })

    @classmethod
    def from_occurrences(cls, occurrences):
        occurrences = _deduplicate_occurrences(occurrences)
        repeat_records = []
        seen = set()
        for occurrence in occurrences:
            if occurrence.v_repeat not in seen:
                seen.add(occurrence.v_repeat)
                repeat_records.append(occurrence.v_repeat)
        return cls(repeat_records, occurrences)

    def __repr__(self):
        string = f"Cluster {self.begin} - {self.end}: \n"
        for s in self.list_clust_dif_rep_seq:
            string += str(s) + "\n"
        return string


class CandidateMerger:
    def __init__(self, dna, list_clusters):
        self.dna = dna
        self.list_clusters = sorted(
            list_clusters,
            key=lambda cluster: (cluster.begin, cluster.end),
        )

        #self._merge_by_distance()
        self._merge_base_on_similarity_to_repeat()

    def _merge_by_distance(self):
        list_merged_by_distance = []
        for cluster in self.list_clusters:
            if not list_merged_by_distance:
                list_merged_by_distance.append(cluster)
            else:
                last_cluster = list_merged_by_distance.pop()
                for el in self.merging_close(last_cluster, cluster):
                    list_merged_by_distance.append(el)

        self.list_clusters = list_merged_by_distance

    @staticmethod
    def merging_close(first_cluster, second_cluster):
        if abs(first_cluster.end - second_cluster.begin) <= 100:
            new_cluster_candidate = CandidateMerger._merged_cluster(
                first_cluster, second_cluster
            )
            return new_cluster_candidate,
        else:
            return first_cluster, second_cluster

    @staticmethod
    def _merged_cluster(first_cluster, second_cluster):
        first_occurrences = getattr(
            first_cluster, "list_cluster_occurrences", None
        )
        second_occurrences = getattr(
            second_cluster, "list_cluster_occurrences", None
        )
        if first_occurrences is not None and second_occurrences is not None:
            return ClusterCandidate.from_occurrences(
                list(first_occurrences) + list(second_occurrences)
            )
        return ClusterCandidate(
            first_cluster.list_cluster_repeats
            + second_cluster.list_cluster_repeats
        )

    def _merge_base_on_similarity_to_repeat(self):
        flag_looping = True
        list_clusters = self.list_clusters
        while flag_looping:
            flag_looping, list_clusters = self.iterative_merging(self.dna, list_clusters)

        self.list_clusters = list_clusters

    @staticmethod
    def mergeable(dna, first_cluster, second_cluster):
        first_cluster_end = first_cluster.end
        second_cluster_begin = second_cluster.begin
        if second_cluster_begin < first_cluster_end:
            return True
        # Cluster ends are zero-based and end-exclusive.  The only sequence
        # between two non-overlapping clusters is therefore [first_end,
        # second_begin); including a base from either cluster changes the
        # 100/500-base gates and can count a seed occurrence outside the gap.
        dna_interval = dna[int(first_cluster_end):int(second_cluster_begin)]
        if len(dna_interval) < 100:
            candidate_first = set(first_cluster.list_clust_dif_rep_seq)
            candidate_second = set(second_cluster.list_clust_dif_rep_seq)
            overlap_candidates = candidate_first.intersection(candidate_second)
            if overlap_candidates:
                return True
        elif len(dna_interval) < 500:
            main_repeat_to_check = (
                second_cluster.list_cluster_repeats[0].sequence
            )
            main_repeat_to_check = main_repeat_to_check[3:-3]
            if main_repeat_to_check:
                occurrences = dna_interval.count(main_repeat_to_check)
                if occurrences >= len(dna_interval) / 160:
                    return True
            if len(second_cluster.list_cluster_repeats) > 3:
                second_repeat_to_check = (
                    second_cluster.list_cluster_repeats[1].sequence
                )
                second_repeat_to_check = second_repeat_to_check[3:-3]
                if second_repeat_to_check:
                    occurrences_second = dna_interval.count(second_repeat_to_check)
                    if occurrences_second >= len(dna_interval) / 160:
                        return True
                third_repeat_to_check = (
                    second_cluster.list_cluster_repeats[2].sequence
                )
                third_repeat_to_check = third_repeat_to_check[3:-3]
                if third_repeat_to_check:
                    occurrences_third = dna_interval.count(third_repeat_to_check)
                    if occurrences_third >= len(dna_interval) / 180:
                        return True
        return False

    def iterative_merging(self, dna, list_clusters):
        temp_list_clusters = []
        for index, _ in enumerate(list_clusters[:-1]):
            first_cluster = list_clusters[index]
            second_cluster = list_clusters[index+1]
            if self.mergeable(dna, first_cluster, second_cluster):
                new_cluster_candidate = self._merged_cluster(
                    first_cluster, second_cluster
                )
                temp_list_clusters.append(new_cluster_candidate)
                temp_list_clusters += list_clusters[index+2:]
                return True, temp_list_clusters
            else:
                temp_list_clusters.append(first_cluster)

        if list_clusters:
            temp_list_clusters.append(list_clusters[-1])
        return False, temp_list_clusters

    def output(self):
        return self.list_clusters


class ClusterMaker:
    def __init__(self, list_repeats_from_vmatch, dna):
        self.list_repeats_from_vmatch = list_repeats_from_vmatch
        self.dna = dna
        self.dna_length = len(dna)

        self.list_clusters = []
        self.list_cluster_sequences = []

        self._read_settings()
        self._find_clusters()
        self._merge_candidates_based_on_similarity()
        self._filter_cluster_duplicates()

    def _read_settings(self):
        self.settings = {
                         "max_length_between_repeats": str(
                             MAX_LOCAL_CHAIN_SPACER_LENGTH
                         ),
                         "flanking_region": "100"}

    def _find_clusters(self):
        max_len_b_rep = int(self.settings["max_length_between_repeats"])
        if max_len_b_rep < 0:
            raise ValueError(
                "maximum length between repeat occurrences must be non-negative"
            )

        source_dna = getattr(self, "dna", None)
        ordered_repeat_candidates = sorted(
            self.list_repeats_from_vmatch,
            key=_v_repeat_sort_key,
        )
        ordered_occurrences = _deduplicate_occurrences([
            occurrence
            for repeat_candidate in ordered_repeat_candidates
            for occurrence in _occurrences_for_v_repeat(
                repeat_candidate, source_dna
            )
        ])

        local_occurrences = []
        local_end = None
        for occurrence in ordered_occurrences:
            if not local_occurrences:
                local_occurrences.append(occurrence)
                local_end = occurrence.end
            elif occurrence.begin - local_end <= max_len_b_rep:
                local_occurrences.append(occurrence)
                local_end = max(local_end, occurrence.end)
            else:
                self.list_clusters.append(
                    ClusterCandidate.from_occurrences(local_occurrences)
                )
                local_occurrences = [occurrence]
                local_end = occurrence.end

        if local_occurrences:
            self.list_clusters.append(
                ClusterCandidate.from_occurrences(local_occurrences)
            )

        #print("Clusters were found: got {} clusters".format(len(self.list_clusters)))

    def _merge_candidates_based_on_similarity(self):
        cm = CandidateMerger(self.dna, self.list_clusters)
        self.list_clusters = cm.output()
        #print("Clusters were merged: now got {} clusters".format(len(self.list_clusters)))

    def _filter_cluster_duplicates(self):
        retained_clusters = []
        for cluster_candidate in self.list_clusters:
            if cluster_candidate.begin < self.dna_length:
                if cluster_candidate.end > self.dna_length:
                    cluster_candidate.goes_through_zero = True
                retained_clusters.append(cluster_candidate)
        self.list_clusters = retained_clusters

    def output(self):
        return self.list_clusters


class FilterApproximationClusters:
    def __init__(self, clusters):
        self.clusters = clusters
        self.clusters_filter_enhanced = []

        self._apply_filter_enhancement()

    def _apply_filter_enhancement(self):
        for cluster in self.clusters:
            list_different_repeats = cluster.list_clust_dif_rep_seq
            fa = FilterApproximation(list_different_repeats)
            missing_repeats_from_filter = fa.output()

            list_different_repeats += missing_repeats_from_filter
            list_different_repeats = sorted(set(list_different_repeats))

            cluster.list_clust_dif_rep_seq = list_different_repeats
            self.clusters_filter_enhanced.append(cluster)

    def output(self):
        return self.clusters_filter_enhanced


class FilterApproximation:
    def __init__(self, list_repeats):
        self.list_repeats = list_repeats
        self.list_repeats = sorted(self.list_repeats)
        self.max_seq = None
        self.min_seq = None
        self.list_missing_candidates = []

        self._relative_path_generation()
        self._make_fasta_with_repeats()

        if len(self.list_repeats) > 2:
            self._obtain_min_max_sequences()
            self._compute_all_the_missing_cases()

        self._clean_up()

    def _relative_path_generation(self):
        full_path = os.path.realpath(__file__)
        self.absolute_path_to_tools = "/" + "/".join(full_path.split("/")[:-2]) + "/tools/"

    def _make_fasta_with_repeats(self):
        with open("clustal_repeats.fa", "w") as f:
            lines = "".join([f">r{index}\n{repeat}\n" for index, repeat in enumerate(self.list_repeats, 1)])
            f.write(lines)

    def _run_clustal_omega_repeats(self):
        subprocess.run(
            ["clustalo", "-i", "clustal_repeats.fa", "-o", "loc_align.txt",
             "--threads", str(external_thread_budget())],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def _find_max_clustal_omega_repeat_results(self):
        if self.align_result_repeats:
            if set([line.strip()[0] for line in self.align_result_repeats[::2]]) == {">"}:
                alignments = [line.strip() for line in self.align_result_repeats[1::2]]
                columns = [[al[index] for al in alignments] for index in range(len(alignments[0]))]
                columns = [[x for x in column if x != "-"] for column in columns]
                columns = [sorted(column) for column in columns]
                # most_common = [max(set(lst).difference(set("-")), key=lst.count) for lst in columns]
                most_common = [max(lst, key=lst.count) for lst in columns]
                self.max_seq = "".join(most_common)

    def _find_min_clustal_omega_repeat_results(self):
        index_start = None
        index_end = None

        if self.align_result_repeats:
            if set([line.strip()[0] for line in self.align_result_repeats[::2]]) == {">"}:
                alignments = [line.strip() for line in self.align_result_repeats[1::2]]
                columns = [[al[index] for al in alignments] for index in range(len(alignments[0]))]
                columns = [[x for x in column if x != "-"] for column in columns]
                columns = [sorted(column) for column in columns]
                most_common = [max(lst, key=lst.count) for lst in columns]
                conserved = "".join(most_common)

                for index in range(len(alignments[0])):
                    column = [al[index] for al in alignments]
                    if "-" not in column:
                        if index_start is None:
                            index_start = index
                        index_end = index
                if index_start is not None and index_end is not None:
                    self.min_seq = conserved[index_start:index_end + 1]

    def _read_clustal_omega_repeat_results(self):
        with open("loc_align.txt") as f:
            self.align_result_repeats = f.readlines()
        os.remove("loc_align.txt")

    def _run_clustal_omega_min_max(self):
        with open("min_max.fa", "w") as f:
            f.write(">max\n")
            f.write(f"{self.max_seq}\n")
            f.write(">min\n")
            f.write(f"{self.min_seq}\n")

        subprocess.run(
            ["clustalo", "-i", "min_max.fa", "-o", "min_max_align.txt",
             "--threads", str(external_thread_budget())],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

    def _obtain_max_sequence(self):
        self._run_clustal_omega_repeats()
        self._read_clustal_omega_repeat_results()
        self._find_max_clustal_omega_repeat_results()

    def _obtain_min_sequence(self):
        self._run_clustal_omega_repeats()
        self._read_clustal_omega_repeat_results()
        self._find_min_clustal_omega_repeat_results()

    def _obtain_min_max_sequences(self):
        self._run_clustal_omega_repeats()
        self._read_clustal_omega_repeat_results()
        self._find_max_clustal_omega_repeat_results()
        self._find_min_clustal_omega_repeat_results()

    def _compute_all_the_missing_cases(self):
        if self.max_seq and self.min_seq:
            if len(self.min_seq) > 4:
                if abs(len(self.max_seq) - len(self.min_seq)) < 16:
                    length_longest = len(self.max_seq)
                    length_shortest = len(self.min_seq)
                    all_possible_substrings = [self.max_seq[i:j + 1] for i in range(length_longest)
                                               for j in range(i, length_longest)
                                               if len(self.max_seq[i:j + 1]) > length_shortest]

                    for candidate in all_possible_substrings:
                        if self.min_seq in candidate:
                            self.list_missing_candidates.append(candidate)

    def _clean_up(self):
        try:
            os.remove("clustal_repeats.fa")
        except FileNotFoundError:
            pass

    def output(self):
        return self.list_missing_candidates


class StartEndEnhancementClusters:
    def __init__(self, clusters):
        self.clusters = clusters
        self.clusters_enhanced_start_end = []

        self._enhance_start_end()

    def _enhance_start_end(self):
        for cluster in self.clusters:
            additional_candidates = []
            list_different_repeats = cluster.list_clust_dif_rep_seq
            for repeat_seq in list_different_repeats:
                for i in range(4):
                    for j in range(4):
                        additional_candidate = repeat_seq[i:len(repeat_seq) - j]
                        additional_candidates.append(additional_candidate)

            enhanced_repeats = sorted(set(list_different_repeats + additional_candidates))
            cluster.list_clust_dif_rep_seq = enhanced_repeats
            self.clusters_enhanced_start_end.append(cluster)

    def output(self):
        return self.clusters_enhanced_start_end


class IntermediateEnhancementClusters:
    def __init__(self, clusters):
        self.clusters = clusters
        self.clusters_enhanced_intermediate = []

        self._enhance_intermediate()

    def _enhance_intermediate(self):
        for cluster in self.clusters:
            list_groups_of_repeats = []
            list_different_repeats = cluster.list_clust_dif_rep_seq
            list_sorted_repeats_by_length = sorted(list_different_repeats, key=lambda x: len(x), reverse=True)
            for repeat_seq in list_sorted_repeats_by_length:
                if not list_groups_of_repeats:
                    list_groups_of_repeats = [[repeat_seq]]
                else:
                    self.add_new(list_groups_of_repeats, repeat_seq)

            list_complete_groups = [self.complete_with_all_intermediate(group_to_complete)
                                    for group_to_complete in list_groups_of_repeats]

            list_all_repeats_candidates = [repeat_seq for group in list_complete_groups for repeat_seq in group]
            list_all_repeats_candidates = sorted(set(list_all_repeats_candidates))

            cluster.list_clust_dif_rep_seq = list_all_repeats_candidates
            self.clusters_enhanced_intermediate.append(cluster)

    @staticmethod
    def add_new(groups, new_seq):
        for group_to_check in groups:
            if new_seq in group_to_check[0]:
                group_to_check.append(new_seq)
                break
        else:
            groups.append([new_seq])

    @staticmethod
    def complete_with_all_intermediate(group_to_complete):
        if len(group_to_complete) < 2:
            return group_to_complete
        else:
            group_to_complete = sorted(group_to_complete, key=lambda x: len(x), reverse=True)
            longest_one = group_to_complete[0]
            length_longest = len(longest_one)
            shortest_one = group_to_complete[-1]
            all_possible_substrings = [longest_one[i:j+1] for i in range(length_longest)
                                       for j in range(i, length_longest)
                                       if len(longest_one[i:j+1]) > len(shortest_one)]

            for candidate in all_possible_substrings:
                if shortest_one in candidate:
                    group_to_complete.append(candidate)

        return group_to_complete

    def output(self):
        return self.clusters_enhanced_intermediate




#                           ARRAYS
###########################################################################
###########################################################################


class _FuzzyOccurrencePathAlternatives(tuple):
    """Tuple with structural equality for cross-engine parity checks."""

    def __new__(cls, values=()):
        return super(_FuzzyOccurrencePathAlternatives, cls).__new__(
            cls, values
        )

    def __eq__(self, other):
        if not isinstance(other, tuple) or len(self) != len(other):
            return False
        return all(
            isinstance(right, FuzzySearch)
            and left.__dict__ == right.__dict__
            for left, right in zip(self, other)
        )


class FuzzySearch:
    _MATCHES_NOT_PROVIDED = object()
    _BOUNDARY_ERROR_CONTRACT = "{i<=3,d<=3,s<=3,i+d+s<=6}"
    _BOUNDARY_SEARCH_POLICY = FuzzyBoundarySearchPolicy()

    @property
    def boundary_path_lattice(self):
        value = self.__dict__.get("boundary_path_lattice")
        return value.resolve() if isinstance(value, CompressedLattice) else value

    @boundary_path_lattice.setter
    def boundary_path_lattice(self, value):
        if isinstance(self.__dict__.get("boundary_path_lattice"), CompressedLattice):
            self._source_provenance_validation_seal = None
            self._occurrence_path_validation_seal = None
            self._boundary_path_validation_seal = None
        self.__dict__["boundary_path_lattice"] = value

    @property
    def occurrence_path_lattice(self):
        value = self.__dict__.get("occurrence_path_lattice")
        return value.resolve() if isinstance(value, CompressedLattice) else value

    @occurrence_path_lattice.setter
    def occurrence_path_lattice(self, value):
        if isinstance(self.__dict__.get("occurrence_path_lattice"), CompressedLattice):
            self._source_provenance_validation_seal = None
            self._occurrence_path_validation_seal = None
            self._boundary_path_validation_seal = None
        self.__dict__["occurrence_path_lattice"] = value

    def _has_compressed_lattices(self):
        return any(isinstance(self.__dict__.get(name), CompressedLattice)
                   for name in ("boundary_path_lattice", "occurrence_path_lattice"))

    def _compact_lattice_storage(self, cache, batch_memo):
        """Consume validated native results without dropping lattice evidence.

        The caller owns the batch memo and releases it after that batch. It
        preserves sibling sharing without keeping original lattices alive
        across all subsequent batches.
        """
        replacements = {}
        for name in ("boundary_path_lattice", "occurrence_path_lattice"):
            value = self.__dict__.get(name)
            if value is None or isinstance(value, CompressedLattice):
                continue
            identity = id(value)
            existing = batch_memo.get(identity)
            if existing is None:
                existing = (value, cache.compact(value))
                batch_memo[identity] = existing
            if existing[0] is not value:
                raise ValueError("fuzzy lattice batch identity changed")
            replacements[name] = existing[1]
        self.__dict__.update(replacements)
        self._source_provenance_validation_seal = None
        self._occurrence_path_validation_seal = None
        self._boundary_path_validation_seal = None

    def __init__(self, sequence, sequence_start, repeat, weighted_error,
                 *, source_dna,
                 precomputed_matches=_MATCHES_NOT_PROVIDED,
                 _defer_alignment=False,
                 _alignment_token_cache=None,
                 _build_occurrence_paths=True,
                 _build_boundary_paths=True,
                 _occurrence_path_lattice=None,
                 _occurrence_path_id=None,
                 _occurrence_path_node_ids=(),
                 _occurrence_path_kind=None,
                 boundary_source_binding=None,
                 _boundary_path_lattice=None,
                 _boundary_path_id=None,
                 _boundary_path_node_ids=(),
                 _boundary_path_kind=None):
        if not isinstance(_defer_alignment, bool):
            raise TypeError("deferred fuzzy alignment flag must be boolean")
        if not isinstance(_build_boundary_paths, bool):
            raise TypeError("fuzzy boundary path flag must be boolean")
        if (
                _alignment_token_cache is not None
                and not isinstance(
                    _alignment_token_cache, _FuzzyAlignmentTokenCache
                )):
            raise TypeError("fuzzy alignment token cache is malformed")
        self.repeat_candidate = repeat
        self.pattern = self.pattern_for(repeat, weighted_error)
        self.sequence = sequence
        self.sequence_start = sequence_start
        self.dot_representation = None
        self.occurrence_path_lattice = _occurrence_path_lattice
        self.occurrence_path_id = _occurrence_path_id
        self.occurrence_path_node_ids = tuple(_occurrence_path_node_ids)
        self.occurrence_path_kind = _occurrence_path_kind
        self.occurrence_path_alternatives = (
            _FuzzyOccurrencePathAlternatives()
        )
        self.boundary_path_lattice = _boundary_path_lattice
        self.boundary_path_id = _boundary_path_id
        self.boundary_path_node_ids = tuple(_boundary_path_node_ids)
        self.boundary_path_kind = _boundary_path_kind
        self.boundary_path_alternatives = (
            _FuzzyOccurrencePathAlternatives()
        )
        self._defer_alignment = _defer_alignment
        self._build_boundary_paths = _build_boundary_paths
        self._alignment_materialized = False
        self._deferred_matching_object = None
        self._alignment_token_cache = (
            _alignment_token_cache or _FuzzyAlignmentTokenCache()
            if _defer_alignment else None
        )
        self._validated_pairwise_alignments = None
        self._alignment_class_tokens = None
        self.list_relative_error_indexes = None
        self._source_provenance_validation_seal = None
        self._occurrence_path_validation_seal = None
        self._boundary_path_validation_seal = None
        self._boundary_source_binding = (
            bind_fuzzy_boundary_source(source_dna)
            if boundary_source_binding is None
            else bind_fuzzy_boundary_source(boundary_source_binding)
        )
        bound_source = self._boundary_source_binding.sequence
        if bound_source is not source_dna and bound_source != source_dna:
            raise ValueError(
                "fuzzy boundary source identity is stale"
            )

        self._match(precomputed_matches, source_dna)
        if self.match_hit is True:
            self._get_repeats()
            self._get_spacers()
            self._get_errors()
            self._get_absolute_starts()
            if _defer_alignment:
                self.list_relative_error_indexes = None
                self._deferred_matching_object = tuple(
                    self.matching_object
                )
            else:
                self._materialize_alignment()
            self._initialize_source_path_provenance(source_dna)
            if self.occurrence_path_lattice is not None:
                self._validate_occurrence_path_link()
            if self.boundary_path_lattice is not None:
                self._validate_boundary_path_link()
            if (
                    self.occurrence_path_lattice is None
                    and self.boundary_path_lattice is None
                    and _build_occurrence_paths):
                self._initialize_occurrence_path_alternatives(
                    weighted_error, source_dna
                )
        self.matching_object = None
        # The cluster sequence can be large and is no longer needed after the
        # derived fields above are populated. Dropping it substantially cuts
        # multiprocessing result traffic.
        self.sequence = None
        self._boundary_source_binding = None

    def _initialize_occurrence_path_alternatives(
            self, weighted_error, source_dna):
        """Attach complete exact and fuzzy-boundary alternative paths."""

        plan = build_occurrence_path_plan(
            self.repeat_candidate,
            self.sequence,
            self.matching_object,
            self.sequence_start,
        )
        if plan.lattice is not None:
            self.occurrence_path_lattice = plan.lattice
            self.occurrence_path_id = plan.path_ids[0]
            self.occurrence_path_node_ids = plan.lattice.paths[0].node_ids
            self.occurrence_path_kind = plan.lattice.paths[0].kind
            self._validate_occurrence_path_link()

            # A bounded lattice is evidence that alternatives exist, not a
            # safe candidate universe.  Never materialize a truncated exact
            # path family.
            if plan.lattice.exact_paths_complete:
                alternatives = []
                for matches, path_id, record in zip(
                        plan.match_paths[1:],
                        plan.path_ids[1:],
                        plan.lattice.paths[1:]):
                    alternatives.append(self.__class__(
                        self.sequence,
                        self.sequence_start,
                        self.repeat_candidate,
                        weighted_error,
                        source_dna=source_dna,
                        precomputed_matches=matches,
                        _defer_alignment=self._defer_alignment,
                        _alignment_token_cache=self._alignment_token_cache,
                        _build_occurrence_paths=False,
                        _build_boundary_paths=self._build_boundary_paths,
                        _occurrence_path_lattice=plan.lattice,
                        _occurrence_path_id=path_id,
                        _occurrence_path_node_ids=record.node_ids,
                        _occurrence_path_kind=record.kind,
                        boundary_source_binding=(
                            self._boundary_source_binding
                        ),
                    ))
                self.occurrence_path_alternatives = (
                    _FuzzyOccurrencePathAlternatives(alternatives)
                )

        if self._build_boundary_paths:
            self._initialize_boundary_path_alternatives(
                weighted_error, source_dna
            )

    def _boundary_search_is_applicable(self, weighted_error):
        """Limit the first integration to proven deletion-shortened paths."""

        return (
            weighted_error == self._BOUNDARY_ERROR_CONTRACT
            and len(self.matching_object) >= 2
            and any(
                match.end() - match.start() < len(self.repeat_candidate)
                and tuple(match.fuzzy_counts)[2] > 0
                for match in self.matching_object
            )
        )

    def _initialize_boundary_path_alternatives(
            self, weighted_error, source_dna):
        """Attach only a complete forced-boundary lattice and its paths."""

        if not self._boundary_search_is_applicable(weighted_error):
            return
        plan = build_fuzzy_boundary_alternative_plan(
            self.repeat_candidate,
            self.sequence,
            self.matching_object,
            self.sequence_start,
            source_dna=self._boundary_source_binding,
            policy=self._BOUNDARY_SEARCH_POLICY,
        )
        lattice = plan.lattice
        self.boundary_path_lattice = lattice
        complete = (
            lattice.alternatives_complete
            and lattice.interval_evaluation_complete
            and lattice.node_universe_complete
            and lattice.path_universe_complete
            and lattice.status == "complete"
        )
        if not complete:
            if plan.alternatives:
                raise ValueError(
                    "incomplete fuzzy boundary lattice emitted candidates"
                )
            return

        paths = {path.path_id: path for path in lattice.paths}
        baseline = paths.get(lattice.baseline_path_id)
        if baseline is None or baseline.kind != "baseline_geometry":
            raise ValueError(
                "complete fuzzy boundary lattice lost its baseline"
            )
        self.boundary_path_id = baseline.path_id
        self.boundary_path_node_ids = baseline.node_ids
        self.boundary_path_kind = baseline.kind
        self._validate_boundary_path_link()

        alternatives = []
        for path in plan.alternatives:
            record = paths.get(path.path_id)
            if (
                    record is None
                    or record.kind != "forced_boundary_alternative"
                    or record.node_ids != path.node_ids):
                raise ValueError(
                    "fuzzy boundary alternative is absent from its lattice"
                )
            alternatives.append(self.__class__(
                self.sequence,
                self.sequence_start,
                self.repeat_candidate,
                weighted_error,
                source_dna=source_dna,
                precomputed_matches=path.matches,
                _defer_alignment=self._defer_alignment,
                _alignment_token_cache=self._alignment_token_cache,
                _build_occurrence_paths=False,
                _build_boundary_paths=self._build_boundary_paths,
                boundary_source_binding=self._boundary_source_binding,
                _boundary_path_lattice=lattice,
                _boundary_path_id=path.path_id,
                _boundary_path_node_ids=path.node_ids,
                _boundary_path_kind=record.kind,
            ))
        self.boundary_path_alternatives = (
            _FuzzyOccurrencePathAlternatives(alternatives)
        )

    def _validate_occurrence_path_link(self):
        """Fail closed unless this object replays its linked lattice path."""

        lattice = self.occurrence_path_lattice
        node_ids_state = self.occurrence_path_node_ids
        starts_state = self.list_absolute_start
        repeats_state = self.list_repeats
        counts_state = self.list_fuzzy_counts
        snapshot = (
            (type(lattice), lattice),
            self.repeat_candidate,
            self.occurrence_path_id,
            self.occurrence_path_kind,
            (type(node_ids_state), tuple(node_ids_state)),
            (type(starts_state), tuple(starts_state)),
            (type(repeats_state), tuple(repeats_state)),
            (type(counts_state), tuple(
                tuple(counts) for counts in counts_state
            )),
        )
        if self._occurrence_path_validation_seal == snapshot:
            return
        if (
                not isinstance(lattice, OccurrencePathLattice)
                or lattice.schema != OCCURRENCE_PATH_SCHEMA
                or lattice.repeat_candidate != self.repeat_candidate):
            raise ValueError("fuzzy occurrence lattice is malformed")
        records = {
            record.path_id: record for record in lattice.paths
        }
        record = records.get(self.occurrence_path_id)
        if record is None:
            raise ValueError("fuzzy occurrence path is absent from its lattice")
        if self.occurrence_path_kind != record.kind:
            raise ValueError("fuzzy occurrence path kind is inconsistent")
        node_ids = tuple(self.occurrence_path_node_ids)
        if not node_ids or len(node_ids) != len(self.list_repeats):
            raise ValueError("fuzzy occurrence path has invalid node geometry")

        record_ids = record.node_ids
        try:
            first_index = record_ids.index(node_ids[0])
        except ValueError:
            raise ValueError(
                "fuzzy occurrence path nodes are absent from its record"
            )
        if record_ids[first_index:first_index + len(node_ids)] != node_ids:
            raise ValueError(
                "local fuzzy occurrence path is not a contiguous source path"
            )

        nodes = {node.node_id: node for node in lattice.nodes}
        for node_id, start, repeat, counts in zip(
                node_ids,
                self.list_absolute_start,
                self.list_repeats,
                self.list_fuzzy_counts):
            node = nodes.get(node_id)
            if (
                    node is None
                    or node.start != start
                    or node.end != start + len(repeat)
                    or node.sequence != repeat
                    or node.fuzzy_counts != tuple(counts)):
                raise ValueError(
                    "fuzzy occurrence path does not reconstruct lattice nodes"
                )
        self._occurrence_path_validation_seal = (
            None if self._has_compressed_lattices() else snapshot
        )

    def _validate_boundary_path_link(self):
        """Fail closed unless this object replays its boundary-lattice path."""

        lattice = self.boundary_path_lattice
        node_ids_state = self.boundary_path_node_ids
        starts_state = self.list_absolute_start
        repeats_state = self.list_repeats
        counts_state = self.list_fuzzy_counts
        snapshot = (
            (type(lattice), lattice),
            self.repeat_candidate,
            self.boundary_path_id,
            self.boundary_path_kind,
            (type(node_ids_state), tuple(node_ids_state)),
            (type(starts_state), tuple(starts_state)),
            (type(repeats_state), tuple(repeats_state)),
            (type(counts_state), tuple(
                tuple(counts) for counts in counts_state
            )),
        )
        if self._boundary_path_validation_seal == snapshot:
            return
        if (
                not isinstance(lattice, FuzzyBoundaryLattice)
                or lattice.schema != FUZZY_BOUNDARY_LATTICE_SCHEMA
                or lattice.repeat_candidate != self.repeat_candidate):
            raise ValueError("fuzzy boundary lattice is malformed")
        records = {record.path_id: record for record in lattice.paths}
        record = records.get(self.boundary_path_id)
        if record is None:
            raise ValueError(
                "fuzzy boundary path is absent from its lattice"
            )
        if self.boundary_path_kind != record.kind:
            raise ValueError("fuzzy boundary path kind is inconsistent")
        node_ids = tuple(self.boundary_path_node_ids)
        if not node_ids or len(node_ids) != len(self.list_repeats):
            raise ValueError("fuzzy boundary path has invalid node geometry")
        try:
            first_index = record.node_ids.index(node_ids[0])
        except ValueError:
            raise ValueError(
                "fuzzy boundary path nodes are absent from its record"
            )
        if (
                record.node_ids[
                    first_index:first_index + len(node_ids)
                ] != node_ids):
            raise ValueError(
                "local fuzzy boundary path is not a contiguous source path"
            )
        nodes = {node.node_id: node for node in lattice.nodes}
        for node_id, start, repeat, counts in zip(
                node_ids,
                self.list_absolute_start,
                self.list_repeats,
                self.list_fuzzy_counts):
            node = nodes.get(node_id)
            if (
                    node is None
                    or node.start != start
                    or node.end != start + len(repeat)
                    or node.sequence != repeat
                    or node.fuzzy_counts != tuple(counts)):
                raise ValueError(
                    "fuzzy boundary path does not reconstruct lattice nodes"
                )
        self._boundary_path_validation_seal = (
            None if self._has_compressed_lattices() else snapshot
        )

    def complete_occurrence_path_family(self):
        """Return baseline plus complete exact and fuzzy-boundary paths."""

        family = [self]
        if self.occurrence_path_lattice is not None:
            self._validate_occurrence_path_link()
            if not self.occurrence_path_lattice.exact_paths_complete:
                if self.occurrence_path_alternatives:
                    raise ValueError(
                        "incomplete fuzzy lattice constructed exact candidates"
                    )
            else:
                exact_alternatives = tuple(
                    self.occurrence_path_alternatives
                )
                expected_ids = tuple(
                    record.path_id
                    for record in self.occurrence_path_lattice.paths[1:]
                )
                if tuple(
                        alternative.occurrence_path_id
                        for alternative in exact_alternatives) != expected_ids:
                    raise ValueError(
                        "fuzzy occurrence alternatives are incomplete"
                    )
                for alternative in exact_alternatives:
                    alternative._validate_occurrence_path_link()
                family.extend(exact_alternatives)

        boundary_lattice = self.boundary_path_lattice
        if boundary_lattice is None:
            return family
        boundary_alternatives = tuple(self.boundary_path_alternatives)
        complete = (
            boundary_lattice.alternatives_complete
            and boundary_lattice.interval_evaluation_complete
            and boundary_lattice.node_universe_complete
            and boundary_lattice.path_universe_complete
            and boundary_lattice.status == "complete"
        )
        if not complete:
            if boundary_alternatives:
                raise ValueError(
                    "incomplete fuzzy boundary lattice constructed candidates"
                )
            return family

        self._validate_boundary_path_link()
        expected_boundary_ids = tuple(
            record.path_id for record in boundary_lattice.paths
            if record.path_id != boundary_lattice.baseline_path_id
        )
        if tuple(
                alternative.boundary_path_id
                for alternative in boundary_alternatives
        ) != expected_boundary_ids:
            raise ValueError(
                "fuzzy boundary alternatives are incomplete"
            )
        for alternative in boundary_alternatives:
            alternative._validate_boundary_path_link()

        def source_geometry(candidate):
            return (
                tuple(candidate.list_absolute_start),
                tuple(candidate.list_repeats),
                tuple(candidate.list_spacers),
                tuple(tuple(counts) for counts in candidate.list_fuzzy_counts),
            )

        seen = {source_geometry(candidate) for candidate in family}
        for alternative in boundary_alternatives:
            geometry = source_geometry(alternative)
            if geometry not in seen:
                family.append(alternative)
                seen.add(geometry)
        return family

    def detach_materialized_path_alternatives(self):
        """Release sibling objects after complete-family enumeration.

        Lattice records and this candidate's path/node identity are the
        downstream provenance contract.  The alternative object tuples are
        only an enumeration cache; retaining them on one accepted baseline
        would otherwise keep every rejected sibling and its deferred match
        payload alive across multiprocessing.
        """

        self.occurrence_path_alternatives = (
            _FuzzyOccurrencePathAlternatives()
        )
        self.boundary_path_alternatives = (
            _FuzzyOccurrencePathAlternatives()
        )
        return self

    def __getstate__(self):
        """Never serialize worker-local caches or validation seals."""

        state = self.__dict__.copy()
        # Keep the public pickle contract and old persisted candidates logical;
        # process-local compression/cache objects never cross that boundary.
        for name in ("boundary_path_lattice", "occurrence_path_lattice"):
            if name in state:
                state[name] = getattr(self, name)
        state["_alignment_token_cache"] = None
        state["_validated_pairwise_alignments"] = None
        state["_alignment_class_tokens"] = None
        state["_source_provenance_validation_seal"] = None
        state["_occurrence_path_validation_seal"] = None
        state["_boundary_path_validation_seal"] = None
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        if (
                self._defer_alignment
                and self._deferred_matching_object is not None):
            self._alignment_token_cache = _FuzzyAlignmentTokenCache()
        self._source_provenance_validation_seal = None
        self._occurrence_path_validation_seal = None
        self._boundary_path_validation_seal = None
        if self.match_hit:
            # Seals are never trusted across IPC. Rebuild them only after the
            # unpickled candidate has passed every authoritative validator.
            self._validate_source_path_provenance()

    @staticmethod
    def pattern_for(repeat, weighted_error):
        return "(?e)" + "({})".format(repeat) + weighted_error

    @classmethod
    def compile_pattern(cls, repeat, weighted_error):
        return regex.compile(cls.pattern_for(repeat, weighted_error))

    def _match(self, precomputed_matches, source_dna):
        if (
                not isinstance(self.sequence, str)
                or not isinstance(self.repeat_candidate, str)
                or not self.repeat_candidate):
            raise ValueError("fuzzy search requires non-empty string sequences")
        if (
                not isinstance(self.sequence_start, int)
                or isinstance(self.sequence_start, bool)
                or self.sequence_start < 0):
            raise ValueError(
                "fuzzy search sequence start must be a non-negative integer"
            )
        if not isinstance(source_dna, str):
            raise ValueError(
                "fuzzy search requires the original source DNA sequence"
            )
        sequence_end = self.sequence_start + len(self.sequence)
        if (
                sequence_end > len(source_dna)
                or source_dna[self.sequence_start:sequence_end] != self.sequence):
            raise ValueError(
                "fuzzy search window does not reconstruct its original "
                "source slice"
            )
        if precomputed_matches is self._MATCHES_NOT_PROVIDED:
            self.matching_object = list(regex.finditer(self.pattern, self.sequence))
        else:
            self.matching_object = list(precomputed_matches)
        previous_end = 0
        for match_index, match in enumerate(self.matching_object):
            start = match.start()
            end = match.end()
            if (
                    not isinstance(start, int)
                    or isinstance(start, bool)
                    or not isinstance(end, int)
                    or isinstance(end, bool)
                    or start < previous_end
                    or end <= start
                    or end > len(self.sequence)):
                raise ValueError(
                    "fuzzy match {} has invalid source coordinates".format(
                        match_index
                    )
                )
            if self.sequence[start:end] != match.group():
                raise ValueError(
                    "fuzzy match {} does not reconstruct its source slice".format(
                        match_index
                    )
                )
            absolute_start = self.sequence_start + start
            absolute_end = self.sequence_start + end
            if source_dna[absolute_start:absolute_end] != match.group():
                raise ValueError(
                    "fuzzy match {} does not reconstruct its original source "
                    "slice".format(match_index)
                )
            previous_end = end
        self.match_hit = bool(self.matching_object)

    def _get_repeats(self):
        self.list_repeats = [match.group() for match in self.matching_object]

    def _get_spacers(self):
        spacer_starts = [match.end() for match in self.matching_object][:-1]
        spacer_end = [match.start() for match in self.matching_object][1:]
        spacer_intervals = [(start, end) for start, end in zip(spacer_starts, spacer_end)]
        self.list_spacers = [self.sequence[interval[0]:interval[1]] for interval in spacer_intervals]

    def _get_errors(self):
        self.number_substitutions = 0
        self.number_insertions = 0
        self.number_deletions = 0
        self.list_fuzzy_counts = []

        for match in self.matching_object:
            s, i, d = match.fuzzy_counts
            self.list_fuzzy_counts.append((s, i, d))
            self.number_substitutions += s
            self.number_insertions += i
            self.number_deletions += d

        self.number_errors = self.number_substitutions + self.number_insertions + self.number_deletions

    def _get_absolute_starts(self):
        self.list_absolute_start = [self.sequence_start + match.start() for match in self.matching_object]
        self.start_end = '{}_{}'.format(self.list_absolute_start[0], self.list_absolute_start[-1])

    def _get_relative_error_indexes(self, matches=None):
        if matches is None:
            matches = self.matching_object
        self.list_relative_error_indexes = [
            normalize_fuzzy_changes(match) for match in matches
        ]

    def _load_validated_alignment_tokens(self, matches):
        if len(matches) != len(self.list_repeats):
            raise ValueError(
                "every fuzzy match must have one source repeat"
            )
        tokens = tuple(
            self._alignment_token_cache.token_for(
                self.repeat_candidate,
                self.sequence_start,
                match,
                absolute_start,
                repeat,
            )
            for match, absolute_start, repeat in zip(
                matches,
                self.list_absolute_start,
                self.list_repeats,
            )
        )
        self.list_relative_error_indexes = [
            [list(indexes) for indexes in token.normalized_changes]
            for token in tokens
        ]
        self._validated_pairwise_alignments = tuple(
            token.pairwise_alignment for token in tokens
        )
        self._alignment_class_tokens = tuple(
            token.class_id for token in tokens
        )

    def _ensure_relative_error_indexes(self):
        matches = self._deferred_matching_object
        if matches is not None and self._alignment_token_cache is not None:
            # Recheck the raw-state snapshots even after an earlier lookup.
            # A mutated sibling match must never inherit a stale valid token.
            self._load_validated_alignment_tokens(matches)
            return
        if self.list_relative_error_indexes is not None:
            return
        if matches is None:
            matches = self.matching_object
        if matches is None:
            raise RuntimeError(
                "deferred fuzzy alignment lost its match coordinates"
            )
        self._get_relative_error_indexes(matches)

    def _materialize_alignment(self):
        """Construct the canonical gapped alignment at most once."""

        if self._alignment_materialized:
            return self
        if not self.match_hit:
            return self
        self._ensure_relative_error_indexes()
        self._get_all_insertions_deletions()
        self._get_gaped_repeats()
        self._alignment_materialized = True
        self._deferred_matching_object = None
        self._alignment_token_cache = None
        self._validated_pairwise_alignments = None
        self._alignment_class_tokens = None
        return self

    def compact_same_start_end_token_key(self):
        """Return the interned exact alignment class used inside one worker."""

        self._ensure_relative_error_indexes()
        if self._alignment_class_tokens is None:
            # The topology quotient is deliberately a deferred-native path.
            # Keep this method exact for defensive/direct callers as well.
            return ("sparse", self.compact_same_start_end_key())
        return self._alignment_class_tokens

    def compact_same_start_end_key(self):
        """Return an exact sparse key for the eventual gapped alignment.

        This is an invertible sparse encoding of the historical dense key.
        Removing spaces from the gapped pattern recovers the seed.  For each
        source repeat, pattern-space columns recover insertions, repeat ``-``
        columns recover deletions, and unequal base columns recover
        substitutions.  Conversely those normalized columns reconstruct the
        same dense alignment uniquely.  Sorting mirrors the canonical merge
        and avoids constructing or retaining its dense strings.
        """

        self._ensure_relative_error_indexes()
        if len(self.list_repeats) != len(self.list_relative_error_indexes):
            raise ValueError(
                "every fuzzy match must have one edit-coordinate set"
            )
        # A deferred token lookup has already built and validated every
        # pairwise alignment.  Eager and legacy direct callers retain the
        # authoritative construction fence here.
        if (
                not self._alignment_materialized
                and self._validated_pairwise_alignments is None):
            for matched_text, normalized_changes in zip(
                    self.list_repeats, self.list_relative_error_indexes):
                build_pairwise_fuzzy_alignment(
                    self.repeat_candidate,
                    matched_text,
                    normalized_changes,
                )
        changes = tuple(
            tuple(tuple(sorted(indexes)) for indexes in error_types)
            for error_types in self.list_relative_error_indexes
        )
        return (
            tuple(self.list_absolute_start),
            tuple(self.list_repeats),
            tuple(self.list_spacers),
            self.repeat_candidate,
            changes,
        )

    def rejected_candidate_witness(
            self, two_repeat_terminal_repair_seed=None):
        """Collapse this source-validated path to its grouping evidence."""

        self._validate_source_path_provenance()
        family_parent = any(
            record.derivation == "local_chain_family_parent"
            for record in self.source_path_provenance
        )
        return RejectedFuzzyCandidateWitness(
            self.list_absolute_start[0],
            self.list_absolute_start[-1] + len(self.list_repeats[-1]),
            self.compact_same_start_end_key(),
            family_parent,
            self.number_errors,
            two_repeat_terminal_repair_seed,
        )

    def _get_all_insertions_deletions(self):
        self.list_all_insertions = []
        self.list_all_deletions = []

        for err_type in self.list_relative_error_indexes:
            self.list_all_insertions += self.unique_gaps_end(err_type[1])
            self.list_all_deletions += self.unique_gaps_end(err_type[2])

        self.list_all_insertions = sorted(set(self.list_all_insertions))
        self.list_all_deletions = sorted(set(self.list_all_deletions))

    @staticmethod
    def gap_introduction(sequence, list_of_gaps, gap_char):
        """Insert gap characters at final-string positions in linear time."""
        if not list_of_gaps:
            return sequence

        pieces = []
        source_index = 0
        for inserted_count, gap_position in enumerate(sorted(list_of_gaps)):
            cut = gap_position - inserted_count
            pieces.append(sequence[source_index:cut])
            pieces.append(gap_char)
            source_index = cut
        pieces.append(sequence[source_index:])
        return "".join(pieces)

    @staticmethod
    def apply_insertions_to_deletions(list_insertions, list_deletions):
        list_deletions_new = []
        for element in list_deletions:
            add = sum(element >= i for i in list_insertions)
            list_deletions_new.append(element + add)
        return list_deletions_new

    @staticmethod
    def unique_gaps_end(list_gaps):
        list_gaps_new = []
        for el in list_gaps:
            while True:
                if el in list_gaps_new:
                    el += 1
                else:
                    list_gaps_new.append(el)
                    break
        return list_gaps_new

    def _get_gaped_repeats(self):
        if self._validated_pairwise_alignments is None:
            alignment = merge_fuzzy_alignments(
                self.repeat_candidate,
                self.list_repeats,
                self.list_relative_error_indexes,
            )
        else:
            alignment = merge_pairwise_fuzzy_alignments(
                self.repeat_candidate,
                self._validated_pairwise_alignments,
            )
        self.gaped_matching_pattern, self.list_gaped_repeats = alignment

    @staticmethod
    def _source_occurrence_id(start, end):
        return "source-occurrence:{}:{}".format(start, end)

    @staticmethod
    def _source_family_id(repeat_candidate, geometry, spacers):
        fields = [repeat_candidate]
        fields.extend(
            "{}:{}:{}".format(
                occurrence.start,
                occurrence.end,
                occurrence.sequence,
            )
            for occurrence in geometry
        )
        fields.extend("spacer:{}".format(spacer) for spacer in spacers)
        payload = "\x1f".join(fields).encode("utf-8")
        return "fuzzy-family:" + hashlib.sha256(payload).hexdigest()

    def _initialize_source_path_provenance(self, source_dna):
        """Record one immutable, source-validated fuzzy-search family root."""

        geometry = []
        for start, repeat in zip(
                self.list_absolute_start, self.list_repeats):
            end = start + len(repeat)
            if source_dna[start:end] != repeat:
                raise ValueError(
                    "fuzzy repeat does not reconstruct its original source "
                    "slice"
                )
            geometry.append(FuzzySourceOccurrence(
                self._source_occurrence_id(start, end),
                start,
                end,
                repeat,
            ))

        for index, spacer in enumerate(self.list_spacers):
            spacer_start = geometry[index].end
            spacer_end = geometry[index + 1].start
            if source_dna[spacer_start:spacer_end] != spacer:
                raise ValueError(
                    "fuzzy spacer does not reconstruct its original source "
                    "slice"
                )

        geometry = tuple(geometry)
        parent_spacers = tuple(self.list_spacers)
        occurrence_ids = tuple(
            occurrence.occurrence_id for occurrence in geometry
        )
        if self.boundary_path_kind == "forced_boundary_alternative":
            provenance_type = (
                FuzzyBoundaryAlternativeSourcePathProvenance
            )
        elif self.occurrence_path_kind == "source_exact_overlay":
            provenance_type = FuzzyExactOverlaySourcePathProvenance
        else:
            provenance_type = FuzzySourcePathProvenance
        provenance = provenance_type(
            FUZZY_SOURCE_PATH_SCHEMA,
            self._source_family_id(
                self.repeat_candidate, geometry, parent_spacers
            ),
            self.repeat_candidate,
            geometry,
            parent_spacers,
            occurrence_ids,
            "complete_fuzzy_search",
            (0, len(geometry)),
        )
        self.source_path_provenance = (provenance,)
        self._validate_source_path_provenance()

    def _source_path_validation_snapshot(self):
        records = getattr(self, "source_path_provenance", None)
        starts_state = self.list_absolute_start
        repeats_state = self.list_repeats
        spacers_state = self.list_spacers
        counts_state = self.list_fuzzy_counts
        occurrence_nodes = getattr(self, "occurrence_path_node_ids", ())
        boundary_nodes = getattr(self, "boundary_path_node_ids", ())
        occurrence_lattice = getattr(self, "occurrence_path_lattice", None)
        boundary_lattice = getattr(self, "boundary_path_lattice", None)
        return (
            (type(records), tuple(
                (type(record), record) for record in records
            ) if isinstance(records, tuple) else None),
            self.repeat_candidate,
            (type(starts_state), tuple(starts_state)),
            (type(repeats_state), tuple(repeats_state)),
            (type(spacers_state), tuple(spacers_state)),
            (type(counts_state), tuple(
                tuple(counts) for counts in counts_state
            )),
            (type(occurrence_lattice), occurrence_lattice),
            getattr(self, "occurrence_path_id", None),
            getattr(self, "occurrence_path_kind", None),
            (type(occurrence_nodes), tuple(occurrence_nodes)),
            (type(boundary_lattice), boundary_lattice),
            getattr(self, "boundary_path_id", None),
            getattr(self, "boundary_path_kind", None),
            (type(boundary_nodes), tuple(boundary_nodes)),
        )

    def _validate_source_path_provenance(self):
        """Fail closed unless every recorded path is an exact parent slice."""

        snapshot = self._source_path_validation_snapshot()
        if self._source_provenance_validation_seal == snapshot:
            return
        records = getattr(self, "source_path_provenance", None)
        if (
                not isinstance(records, tuple)
                or not records
                or any(
                    not isinstance(record, FuzzySourcePathProvenance)
                    for record in records
                )):
            raise ValueError("fuzzy source path provenance is malformed")
        if len(set(records)) != len(records):
            raise ValueError("fuzzy source path provenance is duplicated")
        for record in records:
            self._validate_source_path_record(record)
        if getattr(self, "occurrence_path_lattice", None) is not None:
            self._validate_occurrence_path_link()
        if (
                getattr(self, "boundary_path_lattice", None) is not None
                and getattr(self, "boundary_path_id", None) is not None):
            self._validate_boundary_path_link()
        self._source_provenance_validation_seal = (
            None if self._has_compressed_lattices() else snapshot
        )

    def _validate_source_path_record(self, record):
        """Validate one parent derivation against this exact path geometry."""

        if record.schema != FUZZY_SOURCE_PATH_SCHEMA:
            raise ValueError("fuzzy source path has an unsupported schema")
        if (
                record.repeat_candidate != self.repeat_candidate
                or not isinstance(record.parent_geometry, tuple)
                or not record.parent_geometry
                or len(record.parent_spacers)
                != len(record.parent_geometry) - 1):
            raise ValueError("fuzzy source path has malformed parent geometry")
        expected_family_id = self._source_family_id(
            record.repeat_candidate,
            record.parent_geometry,
            record.parent_spacers,
        )
        if record.family_id != expected_family_id:
            raise ValueError("fuzzy source path family identity is inconsistent")

        for occurrence in record.parent_geometry:
            if (
                    not isinstance(occurrence, FuzzySourceOccurrence)
                    or occurrence.occurrence_id
                    != self._source_occurrence_id(
                        occurrence.start, occurrence.end
                    )
                    or occurrence.end
                    != occurrence.start + len(occurrence.sequence)):
                raise ValueError(
                    "fuzzy source path has malformed source occurrences"
                )
        parent_occurrence_ids = tuple(
            occurrence.occurrence_id
            for occurrence in record.parent_geometry
        )
        if len(set(parent_occurrence_ids)) != len(parent_occurrence_ids):
            raise ValueError("fuzzy source occurrence IDs must be unique")

        span = record.source_index_span
        if (
                not isinstance(span, tuple)
                or len(span) != 2
                or any(
                    not isinstance(index, int) or isinstance(index, bool)
                    for index in span
                )
                or span[0] < 0
                or span[1] > len(record.parent_geometry)
                or span[1] - span[0]
                < (2 if record.derivation == "maximal_local_chain" else 1)):
            raise ValueError("fuzzy source path has an invalid parent span")
        selected_geometry = record.parent_geometry[span[0]:span[1]]
        selected_ids = tuple(
            occurrence.occurrence_id for occurrence in selected_geometry
        )
        if record.occurrence_ids != selected_ids:
            raise ValueError(
                "fuzzy source path occurrence IDs do not match its parent span"
            )
        if record.derivation not in (
                "complete_fuzzy_search", "local_chain_family_parent",
                "maximal_local_chain"):
            raise ValueError("fuzzy source path has an unknown derivation")
        if (
                record.derivation == "complete_fuzzy_search"
                and span != (0, len(record.parent_geometry))):
            raise ValueError(
                "complete fuzzy path must contain every parent occurrence"
            )

        if tuple(self.list_absolute_start) != tuple(
                occurrence.start for occurrence in selected_geometry):
            raise ValueError(
                "fuzzy source path starts do not match its parent geometry"
            )
        if tuple(self.list_repeats) != tuple(
                occurrence.sequence for occurrence in selected_geometry):
            raise ValueError(
                "fuzzy source path repeats do not match its parent geometry"
            )
        expected_spacers = record.parent_spacers[span[0]:span[1] - 1]
        if tuple(self.list_spacers) != expected_spacers:
            raise ValueError(
                "fuzzy source path spacers do not match its parent geometry"
            )

    def _slice_local_chain(self, start_index, end_index):
        """Return an independently aligned view of one repeat-hit chain.

        The parent search has already proven every repeat and spacer against
        the original source window. Rebuilding the fuzzy alignment from the
        selected repeats is nevertheless important: an insertion in a remote
        decoy must not add a gap column to an otherwise unchanged local array.
        """

        self._validate_source_path_provenance()
        repeat_count = len(self.list_repeats)
        if (
                not isinstance(start_index, int)
                or isinstance(start_index, bool)
                or not isinstance(end_index, int)
                or isinstance(end_index, bool)
                or start_index < 0
                or end_index > repeat_count
                or end_index - start_index < 2):
            raise ValueError("local fuzzy chain requires at least two repeats")
        if start_index == 0 and end_index == repeat_count:
            return self

        local = self.__class__.__new__(self.__class__)
        local.repeat_candidate = self.repeat_candidate
        local.pattern = self.pattern
        local.sequence = None
        local.sequence_start = self.sequence_start
        local.dot_representation = None
        local.match_hit = True
        local.matching_object = None
        # Copy the storage reference directly: logical access would decode and
        # pin a complete lattice on every derived local chain.
        local.occurrence_path_lattice = self.__dict__.get("occurrence_path_lattice")
        local.occurrence_path_id = getattr(
            self, "occurrence_path_id", None
        )
        local.occurrence_path_kind = getattr(
            self, "occurrence_path_kind", None
        )
        local.occurrence_path_node_ids = tuple(getattr(
            self, "occurrence_path_node_ids", ()
        )[start_index:end_index])
        local.occurrence_path_alternatives = (
            _FuzzyOccurrencePathAlternatives()
        )
        local.boundary_path_lattice = self.__dict__.get("boundary_path_lattice")
        local.boundary_path_id = getattr(
            self, "boundary_path_id", None
        )
        local.boundary_path_kind = getattr(
            self, "boundary_path_kind", None
        )
        local.boundary_path_node_ids = tuple(getattr(
            self, "boundary_path_node_ids", ()
        )[start_index:end_index])
        local.boundary_path_alternatives = (
            _FuzzyOccurrencePathAlternatives()
        )
        local._boundary_source_binding = None
        local._defer_alignment = self._defer_alignment
        local._alignment_materialized = False
        local._alignment_token_cache = self._alignment_token_cache
        local._validated_pairwise_alignments = None
        local._alignment_class_tokens = None
        local._source_provenance_validation_seal = None
        local._occurrence_path_validation_seal = None
        local._boundary_path_validation_seal = None
        if self._alignment_materialized:
            local._deferred_matching_object = None
        else:
            local._deferred_matching_object = tuple(
                self._deferred_matching_object[start_index:end_index]
            )

        local.list_repeats = list(
            self.list_repeats[start_index:end_index]
        )
        local.list_spacers = list(
            self.list_spacers[start_index:end_index - 1]
        )
        local.list_fuzzy_counts = list(
            self.list_fuzzy_counts[start_index:end_index]
        )
        local.number_substitutions = sum(
            counts[0] for counts in local.list_fuzzy_counts
        )
        local.number_insertions = sum(
            counts[1] for counts in local.list_fuzzy_counts
        )
        local.number_deletions = sum(
            counts[2] for counts in local.list_fuzzy_counts
        )
        local.number_errors = (
            local.number_substitutions
            + local.number_insertions
            + local.number_deletions
        )
        local.list_absolute_start = list(
            self.list_absolute_start[start_index:end_index]
        )
        local.start_end = "{}_{}".format(
            local.list_absolute_start[0], local.list_absolute_start[-1]
        )
        if self.list_relative_error_indexes is None:
            local.list_relative_error_indexes = None
        else:
            local.list_relative_error_indexes = [
                [list(indexes) for indexes in error_types]
                for error_types in self.list_relative_error_indexes[
                    start_index:end_index
                ]
            ]
        if self._alignment_materialized:
            local._materialize_alignment()
        derived_records = []
        for parent_record in self.source_path_provenance:
            parent_start = parent_record.source_index_span[0]
            source_span = (
                parent_start + start_index,
                parent_start + end_index,
            )
            derived_records.append(parent_record.__class__(
                parent_record.schema,
                parent_record.family_id,
                parent_record.repeat_candidate,
                parent_record.parent_geometry,
                parent_record.parent_spacers,
                tuple(
                    occurrence.occurrence_id
                    for occurrence in parent_record.parent_geometry[
                        source_span[0]:source_span[1]
                    ]
                ),
                "maximal_local_chain",
                source_span,
            ))
        local.source_path_provenance = tuple(sorted(set(derived_records)))
        local._validate_source_path_provenance()
        return local

    def split_into_local_chains(
            self, max_spacer_length=MAX_LOCAL_CHAIN_SPACER_LENGTH):
        """Partition ordered matches into maximal local candidate chains.

        Adjacent fuzzy matches are connected only when their exact intervening
        source spacer is no longer than ``max_spacer_length``. Runs are maximal
        under that edge rule and remain in source order. Singleton runs are
        evidence hits, not array candidates, so they are omitted.

        When every match belongs to one valid chain, return ``self`` unchanged.
        This preserves the established representation and avoids rebuilding
        its alignment on the overwhelmingly common path.
        """

        if (
                not isinstance(max_spacer_length, int)
                or isinstance(max_spacer_length, bool)
                or max_spacer_length < 0):
            raise ValueError(
                "maximum local-chain spacer length must be a non-negative "
                "integer"
            )
        if not self.match_hit:
            return []

        repeat_count = len(self.list_repeats)
        if repeat_count < 2:
            return []
        break_indexes = [
            spacer_index + 1
            for spacer_index, spacer in enumerate(self.list_spacers)
            if len(spacer) > max_spacer_length
        ]
        if not break_indexes:
            return [self]

        boundaries = [0] + break_indexes + [repeat_count]
        return [
            self._slice_local_chain(start_index, end_index)
            for start_index, end_index in zip(
                boundaries[:-1], boundaries[1:]
            )
            if end_index - start_index >= 2
        ]

    def complete_local_chain_family(
            self, max_spacer_length=MAX_LOCAL_CHAIN_SPACER_LENGTH):
        """Return the source-valid parent followed by every maximal subchain.

        ``split_into_local_chains`` retains its historical subchain-only API.
        Detection uses this complete-family view so a locality split never
        destroys the original occurrence lattice that produced its children.
        """

        self._validate_source_path_provenance()
        local_chains = self.split_into_local_chains(max_spacer_length)
        if len(local_chains) == 1 and local_chains[0] is self:
            return [self]
        parent_records = []
        for record in self.source_path_provenance:
            parent_records.append(record._replace(
                derivation="local_chain_family_parent"
            ))
        self.source_path_provenance = tuple(sorted(set(parent_records)))
        self._validate_source_path_provenance()
        return [self] + local_chains

    def __repr__(self):
        return self.dot_repr()

    @staticmethod
    def apply_dots(repeat_with_gaps, matching_pattern_with_gaps):
        dotted_repeat = ''
        if len(repeat_with_gaps) != len(matching_pattern_with_gaps):
            raise ValueError
        for char_repeat, char_pattern in zip(repeat_with_gaps, matching_pattern_with_gaps):
            if char_repeat == char_pattern:
                if char_repeat == ' ':
                    dotted_repeat += char_repeat
                else:
                    dotted_repeat += '.'
            else:
                dotted_repeat += char_repeat

        return dotted_repeat

    def _compute_dot_representation(self):
        if self.match_hit:
            self._materialize_alignment()
        string = ""
        if self.match_hit:
            max_length_start_index = max(len(str(start)) for start in self.list_absolute_start) + 3
            try:
                max_length_spacer = max(len(spacer) for spacer in self.list_spacers) + 3
            except ValueError:
                max_length_spacer = 0

            for index, gaped_repeat in enumerate(self.list_gaped_repeats):
                repeat_start_index = self.list_absolute_start[index] + 1
                n_gaps_after_start = max_length_start_index - len(str(self.list_absolute_start[index]))

                if index < len(self.list_spacers):
                    spacer = self.list_spacers[index]
                else:
                    spacer = ""
                n_gaps_after_spacer = max_length_spacer - len(spacer)

                s, i, d = self.list_fuzzy_counts[index]
                errors = "s:{} i:{} d{}".format(s, i, d)
                dotted_repeats = self.apply_dots(gaped_repeat, self.gaped_matching_pattern)

                string += "{}{}{}  {}{}{}{}\n".format(repeat_start_index,
                                                      " " * n_gaps_after_start,
                                                      dotted_repeats,
                                                      spacer,
                                                      " " * n_gaps_after_spacer,
                                                      errors,
                                                      self.list_relative_error_indexes[index])

            string += "_" * 100 + "\n"

            string += " " * max_length_start_index + self.gaped_matching_pattern
            string += " " * (max_length_spacer + 2) + "s:{} i:{} d:{}".format(self.number_substitutions,
                                                                              self.number_insertions,
                                                                              self.number_deletions) + "\n"

            string += "_" * 100 + "\n"
            string += " " * max_length_start_index + self.repeat_candidate

        self.dot_representation = string

    def dot_repr(self):
        if self.dot_representation is None:
            self._compute_dot_representation()
        return self.dot_representation
