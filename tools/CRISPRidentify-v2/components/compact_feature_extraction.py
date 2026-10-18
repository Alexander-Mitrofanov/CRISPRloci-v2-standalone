"""Extract exactly the ten inputs used by the compact CRISPR classifier.

Classification uses only the ten features in this module and never invokes
Prodigal, HMMER, or the second historical BLAST database.
"""

from collections import OrderedDict
import math
import os
import re
import tempfile
import threading

import numpy as np

from components.external_tools import ExternalToolError
from components.external_tools import run_external_tool
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.sequence_path_vectorizer import vectorize_sequence_paths


PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BLAST_DATABASE = os.path.join(
    PROJECT_ROOT, "tools", "blasting", "Verified_repeats_dataset1.fa"
)
BLAST_SUFFIXES = ("", ".nhr", ".nin", ".nsq")
EXTERNAL_SCORE_CACHE_MAXSIZE = 16384
_BLAST_COMMAND = "blastn"
_BLAST_WORD_SIZE = "6"
_BLAST_OUTPUT_FORMAT = "6 qseqid bitscore"
_RNAFOLD_COMMAND = "RNAfold"
_RNAFOLD_ARGUMENTS = ("--noLP", "--noPS")
_CANDIDATE_IDENTIFIER = re.compile(r"^candidate_(\d+)$")
_RNAFOLD_ENERGY = re.compile(
    r"\(\s*(-?(?:\d+(?:\.\d*)?|\.\d+))\s*\)\s*$"
)
_SPACER_BIN_RIGHT_BOUNDARIES = (15, 25, 36, 52, 72)


class _BoundedScoreCache:
    """Small process-local, thread-safe LRU for deterministic tool scores."""

    def __init__(self, maximum):
        if not isinstance(maximum, int) or isinstance(maximum, bool):
            raise TypeError("external score cache maximum must be an integer")
        if maximum < 1:
            raise ValueError("external score cache maximum must be positive")
        self.maximum = maximum
        self._values = OrderedDict()
        self._lock = threading.RLock()

    def get_many(self, keys):
        hits = {}
        with self._lock:
            for key in keys:
                if key not in self._values:
                    continue
                value = self._values[key]
                self._values.move_to_end(key)
                hits[key] = value
        return hits

    def put_many(self, items):
        # Materialize before taking the lock so an invalid/lazy input cannot
        # leave a partially updated cache. Tool callers similarly invoke this
        # method only after the complete external result has parsed.
        values = tuple((key, float(value)) for key, value in items)
        with self._lock:
            for key, value in values:
                self._values[key] = value
                self._values.move_to_end(key)
            while len(self._values) > self.maximum:
                self._values.popitem(last=False)

    def clear(self):
        with self._lock:
            self._values.clear()

    def __len__(self):
        with self._lock:
            return len(self._values)


_BLAST_SCORE_CACHE = _BoundedScoreCache(EXTERNAL_SCORE_CACHE_MAXSIZE)
_MFE_SCORE_CACHE = _BoundedScoreCache(EXTERNAL_SCORE_CACHE_MAXSIZE)


def _clear_external_score_caches_for_tests():
    """Clear process-local compact-feature caches for isolated tests."""

    _BLAST_SCORE_CACHE.clear()
    _MFE_SCORE_CACHE.clear()


def _external_score_cache_sizes_for_tests():
    """Return cache cardinalities without exposing mutable cache state."""

    return {
        "blast": len(_BLAST_SCORE_CACHE),
        "mfe": len(_MFE_SCORE_CACHE),
    }


def _flatten_candidates(candidates_by_key):
    return [
        (key, index, candidate)
        for key, candidates in candidates_by_key.items()
        for index, candidate in enumerate(candidates)
    ]


def _result_template(candidates_by_key):
    return {
        key: [0.0] * len(candidates)
        for key, candidates in candidates_by_key.items()
    }


def _write_candidate_fasta(records, path, sequence_getter):
    with open(path, "w") as output:
        for global_index, (_key, _local_index, candidate) in enumerate(records):
            sequence = sequence_getter(candidate)
            if not sequence:
                raise ValueError(
                    "candidate {} produced an empty compact-feature sequence"
                    .format(global_index)
                )
            output.write(">candidate_{}\n{}\n".format(global_index, sequence))


def _consensus_cache_plan(records, cache, policy):
    """Return ordered cache entries, hits, and first-seen unique misses."""

    entries = []
    first_by_key = {}
    for global_index, record in enumerate(records):
        sequence = record[2].consensus
        if not sequence:
            raise ValueError(
                "candidate {} produced an empty compact-feature sequence"
                .format(global_index)
            )
        cache_key = (sequence, policy)
        entries.append((record, cache_key))
        if cache_key not in first_by_key:
            first_by_key[cache_key] = (record, sequence)

    hits = cache.get_many(first_by_key)
    misses = [
        (cache_key, record, sequence)
        for cache_key, (record, sequence) in first_by_key.items()
        if cache_key not in hits
    ]
    return entries, hits, misses


def _assign_consensus_scores(entries, scores_by_cache_key, destination):
    """Fan cached or freshly scored consensus values back to input order."""

    for (key, candidate_index, _candidate), cache_key in entries:
        destination[key][candidate_index] = float(
            scores_by_cache_key[cache_key]
        )


def _blast_database_identity(database):
    """Return a cheap process-local cache identity for the frozen database."""

    identity = []
    for suffix in BLAST_SUFFIXES:
        path = os.path.abspath(database + suffix)
        metadata = os.stat(path)
        identity.append((
            path,
            metadata.st_dev,
            metadata.st_ino,
            metadata.st_size,
            metadata.st_mtime_ns,
        ))
    return tuple(identity)


def _require_blast_database():
    missing = [
        os.path.abspath(BLAST_DATABASE + suffix)
        for suffix in BLAST_SUFFIXES
        if not os.path.isfile(BLAST_DATABASE + suffix)
    ]
    if missing:
        raise FileNotFoundError(
            "compact-model BLAST database is incomplete; missing {}".format(
                ", ".join(missing)
            )
        )
    return os.path.abspath(BLAST_DATABASE)


def _parse_blast_scores(path, candidate_count):
    scores = [0.0] * candidate_count
    with open(path) as output:
        for line_number, raw_line in enumerate(output, 1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 2:
                raise ExternalToolError(
                    "invalid BLAST row {}: {!r}".format(line_number, line)
                )
            match = _CANDIDATE_IDENTIFIER.fullmatch(fields[0])
            if match is None:
                raise ExternalToolError(
                    "BLAST returned an unknown query identifier: {}".format(
                        fields[0]
                    )
                )
            candidate_index = int(match.group(1))
            if not 0 <= candidate_index < candidate_count:
                raise ExternalToolError(
                    "BLAST candidate index is out of range: {}".format(
                        candidate_index
                    )
                )
            try:
                score = float(fields[1])
            except ValueError as error:
                raise ExternalToolError(
                    "BLAST returned a nonnumeric bit score: {}".format(
                        fields[1]
                    )
                ) from error
            if not math.isfinite(score) or score < 0.0:
                raise ExternalToolError(
                    "BLAST returned an invalid bit score: {}".format(score)
                )
            scores[candidate_index] = max(scores[candidate_index], score)
    return scores


def _blast_scores(candidates_by_key):
    records = _flatten_candidates(candidates_by_key)
    result = _result_template(candidates_by_key)
    if not records:
        return result
    database = _require_blast_database()
    policy = (
        _BLAST_COMMAND,
        ("-word_size", _BLAST_WORD_SIZE),
        ("-outfmt", _BLAST_OUTPUT_FORMAT),
        _blast_database_identity(database),
    )
    entries, resolved, misses = _consensus_cache_plan(
        records, _BLAST_SCORE_CACHE, policy
    )
    if not misses:
        _assign_consensus_scores(entries, resolved, result)
        return result

    miss_records = [record for _key, record, _sequence in misses]
    with tempfile.TemporaryDirectory(
            prefix="crispridentify_compact_blast_") as temporary_directory:
        query_path = os.path.join(temporary_directory, "queries.fa")
        output_path = os.path.join(temporary_directory, "blast.tsv")
        _write_candidate_fasta(
            miss_records, query_path, lambda candidate: candidate.consensus
        )
        run_external_tool(
            _BLAST_COMMAND,
            [
                "-query", query_path,
                "-db", database,
                "-word_size", _BLAST_WORD_SIZE,
                "-outfmt", _BLAST_OUTPUT_FORMAT,
                "-out", output_path,
            ],
        )
        scores = _parse_blast_scores(output_path, len(miss_records))

    fresh = {
        cache_key: float(score)
        for (cache_key, _record, _sequence), score in zip(misses, scores)
    }
    # Cache only after the complete subprocess output has parsed successfully.
    _BLAST_SCORE_CACHE.put_many(fresh.items())
    resolved.update(fresh)
    _assign_consensus_scores(entries, resolved, result)
    return result


def _parse_rnafold_scores(path, expected_sequences):
    scores = [None] * len(expected_sequences)
    with open(path) as output:
        lines = [line.rstrip("\n") for line in output]
    if len(lines) % 3 != 0:
        raise ExternalToolError(
            "RNAfold output is not a sequence of three-line records"
        )
    for offset in range(0, len(lines), 3):
        header, folded_sequence, structure_line = lines[offset:offset + 3]
        if not header.startswith(">"):
            raise ExternalToolError("RNAfold record is missing a FASTA header")
        match = _CANDIDATE_IDENTIFIER.fullmatch(header[1:])
        if match is None:
            raise ExternalToolError(
                "RNAfold returned an unknown candidate identifier: {}".format(
                    header[1:]
                )
            )
        candidate_index = int(match.group(1))
        if not 0 <= candidate_index < len(expected_sequences):
            raise ExternalToolError(
                "RNAfold candidate index is out of range: {}".format(
                    candidate_index
                )
            )
        if scores[candidate_index] is not None:
            raise ExternalToolError(
                "RNAfold returned a duplicate candidate index: {}".format(
                    candidate_index
                )
            )
        normalized = folded_sequence.upper().replace("U", "T")
        if normalized != expected_sequences[candidate_index].upper():
            raise ExternalToolError(
                "RNAfold sequence does not match candidate {}".format(
                    candidate_index
                )
            )
        match = _RNAFOLD_ENERGY.search(structure_line)
        if match is None:
            raise ExternalToolError(
                "RNAfold output is missing an MFE score for candidate {}"
                .format(candidate_index)
            )
        score = float(match.group(1))
        if not math.isfinite(score):
            raise ExternalToolError("RNAfold returned a non-finite MFE score")
        scores[candidate_index] = score
    missing = [index for index, score in enumerate(scores) if score is None]
    if missing:
        raise ExternalToolError(
            "RNAfold did not return every candidate: {}".format(missing)
        )
    return scores


def _mfe_scores(candidates_by_key):
    records = _flatten_candidates(candidates_by_key)
    result = _result_template(candidates_by_key)
    if not records:
        return result
    policy = (_RNAFOLD_COMMAND, tuple(_RNAFOLD_ARGUMENTS))
    entries, resolved, misses = _consensus_cache_plan(
        records, _MFE_SCORE_CACHE, policy
    )
    if not misses:
        _assign_consensus_scores(entries, resolved, result)
        return result

    miss_records = [record for _key, record, _sequence in misses]
    expected_sequences = [sequence for _key, _record, sequence in misses]
    with tempfile.TemporaryDirectory(
            prefix="crispridentify_compact_mfe_") as temporary_directory:
        input_path = os.path.join(temporary_directory, "consensus.fa")
        output_path = os.path.join(temporary_directory, "rnafold.txt")
        _write_candidate_fasta(
            miss_records, input_path, lambda candidate: candidate.consensus
        )
        with open(input_path, "rb") as input_file, open(
                output_path, "wb") as output_file:
            run_external_tool(
                _RNAFOLD_COMMAND,
                list(_RNAFOLD_ARGUMENTS),
                stdin=input_file,
                stdout=output_file,
            )
        scores = _parse_rnafold_scores(output_path, expected_sequences)

    fresh = {
        cache_key: float(score)
        for (cache_key, _record, _sequence), score in zip(misses, scores)
    }
    # Cache only after the complete subprocess output has parsed successfully.
    _MFE_SCORE_CACHE.put_many(fresh.items())
    resolved.update(fresh)
    _assign_consensus_scores(entries, resolved, result)
    return result


class _BulkSimilarity:
    """Exact EDeN-compatible similarity for repeats and spacers."""

    def __init__(self, candidates_by_key):
        self.candidates_by_key = candidates_by_key

    @staticmethod
    def _contains_sequence_content(sequences):
        return any(set(sequence).difference(("", "-")) for sequence in sequences)

    @classmethod
    def _eligible_sequences(cls, sequences):
        if not sequences or not cls._contains_sequence_content(sequences):
            return ()

        nonempty_sequences = tuple(
            sequence for sequence in sequences if sequence
        )
        # A zero/one-row kernel has the fixed similarity 1.0 and _similarity
        # returns before consulting the matrix. Do not vectorize a sequence
        # whose row therefore cannot affect any score.
        if len(nonempty_sequences) <= 1:
            return ()
        return nonempty_sequences

    @staticmethod
    def _find_component(parent, sequence):
        root = sequence
        while parent[root] != root:
            root = parent[root]
        while parent[sequence] != sequence:
            next_sequence = parent[sequence]
            parent[sequence] = root
            sequence = next_sequence
        return root

    @classmethod
    def _union_components(cls, parent, first_seen, first, second):
        first_root = cls._find_component(parent, first)
        second_root = cls._find_component(parent, second)
        if first_root == second_root:
            return
        # Retaining the earliest root gives components a deterministic order
        # without affecting vector or kernel arithmetic.
        if first_seen[first_root] > first_seen[second_root]:
            first_root, second_root = second_root, first_root
        parent[second_root] = first_root

    def _component_plan(self, sequence_getter):
        """Group candidates whose exact sequence vectors can be shared.

        A candidate's sequence list is a hyperedge. Two candidates need the
        same in-memory vector batch exactly when their hyperedges are connected
        by an identical sequence. Processing the connected components one at a
        time therefore retains all cross-candidate vector reuse while bounding
        live feature dictionaries by the largest component rather than the
        whole genome-wide candidate set.
        """

        parent = {}
        first_seen = {}
        records = []
        for key, candidates in self.candidates_by_key.items():
            for candidate_index, candidate in enumerate(candidates):
                sequences = self._eligible_sequences(
                    sequence_getter(candidate)
                )
                if not sequences:
                    continue
                for sequence in sequences:
                    if sequence not in parent:
                        parent[sequence] = sequence
                        first_seen[sequence] = len(first_seen)
                anchor = sequences[0]
                for sequence in sequences[1:]:
                    self._union_components(
                        parent, first_seen, anchor, sequence
                    )
                records.append((key, candidate_index, sequences))

        sequences_by_root = {}
        for sequence in first_seen:
            root = self._find_component(parent, sequence)
            sequences_by_root.setdefault(root, []).append(sequence)
        records_by_root = {}
        for record in records:
            root = self._find_component(parent, record[2][0])
            records_by_root.setdefault(root, []).append(record)

        return tuple(
            (sequences, records_by_root[root])
            for root, sequences in sequences_by_root.items()
        )

    @staticmethod
    def _similarity(matrix, indexes):
        if len(indexes) <= 1:
            return 1.0
        selected = matrix[indexes]
        kernel = selected.dot(selected.T).toarray()
        count = len(indexes)
        return float((np.sum(kernel) - count) / (count ** 2 - count))

    def _score_sequence_kind(self, sequence_getter):
        result = {
            key: [1.0] * len(candidates)
            for key, candidates in self.candidates_by_key.items()
        }
        for unique_sequences, records in self._component_plan(
                sequence_getter):
            matrix = vectorize_sequence_paths(unique_sequences)
            sequence_to_index = {
                sequence: index
                for index, sequence in enumerate(unique_sequences)
            }
            for key, candidate_index, sequences in records:
                indexes = [
                    sequence_to_index[sequence] for sequence in sequences
                ]
                result[key][candidate_index] = self._similarity(
                    matrix, indexes
                )
            # SparseFeatureMatrix owns the component's Python feature
            # dictionaries and dot cache. Drop them before vectorizing the next
            # component so peak memory follows the largest component.
            del matrix
        return result

    @staticmethod
    def _repeat_sequences(candidate):
        return candidate.list_repeats

    @staticmethod
    def _spacer_sequences(candidate):
        spacers = list(candidate.list_spacers)
        if spacers and spacers[-1] == "":
            spacers = spacers[:-1]
        if set(spacers) == {"-"}:
            return ()
        return spacers

    def output(self):
        repeats = self._score_sequence_kind(self._repeat_sequences)
        spacers = self._score_sequence_kind(self._spacer_sequences)
        result = {}
        for key, candidates in self.candidates_by_key.items():
            values = []
            for candidate_index, _candidate in enumerate(candidates):
                repeat_similarity = repeats[key][candidate_index]
                if repeat_similarity > 0.999:
                    repeat_similarity = 1.0
                values.append((
                    float(repeat_similarity),
                    float(spacers[key][candidate_index]),
                ))
            result[key] = values
        return result


def _spacer_evenness(spacers):
    bins = [0] * (len(_SPACER_BIN_RIGHT_BOUNDARIES) + 1)
    for spacer in spacers:
        for index, boundary in enumerate(_SPACER_BIN_RIGHT_BOUNDARIES):
            if len(spacer) <= boundary:
                bins[index] += 1
                break
        else:
            bins[-1] += 1
    count = sum(bins)
    if count == 0:
        return 0.0
    normalized = [value / float(count) for value in bins]
    entropy = sum(
        -value * math.log(value + 0.00000001)
        for value in normalized
    )
    return float(round(entropy / math.log(len(bins)), 6))


def _structural_features(candidate, similarity):
    repeats = candidate.list_repeats
    gapped_repeats = candidate.list_repeats_gaped
    spacers = list(candidate.list_spacers)
    average_spacers = spacers[:-1] if spacers and spacers[-1] == "" else spacers
    consensus_gapped = candidate.consensus_gaped
    mismatch_count = sum(
        repeat_character != consensus_character
        for repeat in gapped_repeats
        for repeat_character, consensus_character in zip(
            repeat, consensus_gapped
        )
    )
    nucleotide_count = sum(
        character in "ATGC"
        for repeat in repeats
        for character in repeat
    )
    at_count = sum(
        character in "AT"
        for repeat in repeats
        for character in repeat
    )
    at_richness = (
        float(at_count) / nucleotide_count if nucleotide_count else 0.0
    )
    average_spacer_length = (
        sum(len(spacer) for spacer in average_spacers)
        / float(len(average_spacers))
        if average_spacers else 0.0
    )
    return (
        float(max(len(repeat) for repeat in repeats)),
        float(len(repeats)),
        float(similarity[0]),
        float(at_richness),
        float(average_spacer_length),
        float(similarity[1]),
        float(mismatch_count),
        _spacer_evenness(spacers),
    )


def extract_compact_feature_rows(candidates_by_key):
    """Return ``(key, candidate, row)`` tuples in stable candidate order."""

    similarities = _BulkSimilarity(candidates_by_key).output()
    mfe_scores = _mfe_scores(candidates_by_key)
    blast_scores = _blast_scores(candidates_by_key)
    rows = []
    for key, candidates in candidates_by_key.items():
        for index, candidate in enumerate(candidates):
            values = _structural_features(
                candidate, similarities[key][index]
            ) + (
                float(mfe_scores[key][index]),
                float(blast_scores[key][index]),
            )
            if len(values) != len(COMPACT_FEATURE_SCHEMA):
                raise RuntimeError("compact feature width is inconsistent")
            row = np.asarray(values, dtype=float)
            if not np.isfinite(row).all():
                raise ValueError("compact features must all be finite")
            rows.append((key, candidate, row))
    return rows
