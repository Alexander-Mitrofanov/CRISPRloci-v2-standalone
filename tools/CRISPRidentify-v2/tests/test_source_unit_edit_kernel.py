from concurrent.futures import ThreadPoolExecutor
import itertools
import json
from pathlib import Path
import random
import unittest
from unittest.mock import patch

from Bio.Align import PairwiseAligner

from components import components_evaluated_arrays_enhancement as enhancement
from components.fast_run_rescue import candidate_score_fingerprint
from tests.test_repeat_boundary_repair import synthetic


def independent_distance(first, second):
    """Two-row unit-cost DP independent of production kernels and caches."""
    previous = list(range(len(second) + 1))
    for first_index, first_base in enumerate(first, 1):
        current = [first_index]
        for second_index, second_base in enumerate(second, 1):
            current.append(min(
                previous[second_index] + 1,
                current[second_index - 1] + 1,
                previous[second_index - 1] + (first_base != second_base),
            ))
        previous = current
    return previous[-1]


def all_strings(alphabet, maximum):
    return ["".join(chars) for length in range(maximum + 1)
            for chars in itertools.product(alphabet, repeat=length)]


class SourceUnitEditKernelTests(unittest.TestCase):
    def setUp(self):
        enhancement.clear_source_alignment_caches()

    def tearDown(self):
        enhancement.clear_source_alignment_caches()

    def test_explicit_global_unit_configuration_and_forward_strand(self):
        aligner = PairwiseAligner()
        # Deliberately wrong initial settings must be overwritten completely.
        aligner.mode = "local"
        aligner.match_score = 9
        aligner.mismatch_score = -3
        aligner.wildcard = "N"
        aligner.open_left_insertion_score = 0
        aligner.extend_right_deletion_score = 0
        with patch.object(enhancement, "PairwiseAligner", return_value=aligner):
            scorer = enhancement._make_source_unit_edit_scorer()
        self.assertEqual("global", aligner.mode)
        self.assertEqual("Needleman-Wunsch", aligner.algorithm)
        self.assertIsNone(aligner.substitution_matrix)
        self.assertIsNone(aligner.wildcard)
        self.assertEqual(0, aligner.match_score)
        self.assertEqual(-1, aligner.mismatch_score)
        names = ["{}_{}_{}_score".format(action, location, kind)
                 for action in ("open", "extend")
                 for location in ("internal", "left", "right")
                 for kind in ("insertion", "deletion")]
        self.assertEqual(12, len(names))
        self.assertTrue(all(getattr(aligner, name) == -1 for name in names))
        self.assertEqual(-2, scorer("ACG", "CGT"))  # reverse complement would score zero
        self.assertEqual(-1, scorer("ANA", "AAA"))  # N is a literal, not wildcard
        self.assertEqual(-1, scorer("a", "A"))

    def test_exhaustive_short_pairs_against_independent_dp(self):
        for alphabet, maximum in (("AC", 5), ("AN-", 3)):
            strings = all_strings(alphabet, maximum)
            for first in strings:
                for second in strings:
                    self.assertEqual(independent_distance(first, second),
                                     enhancement._source_repeat_edit_distance_uncached(first, second),
                                     (first, second))

    def test_two_thousand_seeded_pairs_across_full_fast_range(self):
        rng = random.Random(2026092604)
        alphabet = "ACGTNRYKMSWBDHV-?acgt.!"
        for index in range(2000):
            first = "".join(rng.choice(alphabet) for _ in range(rng.randrange(1, 65)))
            if index % 2:
                # Related strings exercise terminal/internal indels and tied
                # repeated motifs, in addition to unrelated random pairs.
                position = rng.randrange(len(first) + 1)
                second = (first[:position] + rng.choice(alphabet) + first[position:])[:64]
                position = rng.randrange(len(second))
                second = second[:position] + second[position + 1:]
            else:
                second = "".join(rng.choice(alphabet) for _ in range(rng.randrange(1, 65)))
            expected = independent_distance(first, second)
            self.assertEqual(expected, enhancement._source_repeat_edit_distance_uncached(first, second))
            self.assertEqual(expected, enhancement._source_repeat_edit_distance_uncached(second, first))

    def test_ascii_literals_indels_ties_and_fallback_boundaries(self):
        pairs = [(chr(code) + "N", "N" + chr(code)) for code in range(128)]
        pairs += [("AAAA", "AAA"), ("ATATAT", "TATATA"),
                  ("TACGT", "ACGT"), ("ACGTT", "ACGT"),
                  ("ACGGT", "ACGT"), ("ACGT", "TACGT"),
                  ("ACGT", "ACGTT"), ("ACGT", "ACGGT"),
                  ("A" * 64, "A" * 63 + "C")]
        for first, second in pairs:
            self.assertEqual(independent_distance(first, second),
                             enhancement._source_repeat_edit_distance_uncached(first, second))
        with patch.object(enhancement, "_SOURCE_UNIT_EDIT_SCORE",
                          side_effect=AssertionError("fallback reached C scorer")):
            for first, second in (("", "ACG"), ("ACG", ""), ("", ""),
                                  ("ACG", "ACG"), ("ACé", "ACG"),
                                  ("A" * 65, "A" * 64 + "C")):
                self.assertEqual(independent_distance(first, second),
                                 enhancement._source_repeat_edit_distance_uncached(first, second))

    def test_unexpected_scalar_outputs_fail_visibly(self):
        for score in (float("nan"), float("inf"), -.5, 1.0, -100.0, False, "-1", None):
            with patch.object(enhancement, "_SOURCE_UNIT_EDIT_SCORE", return_value=score):
                with self.assertRaises(RuntimeError):
                    enhancement._source_repeat_edit_distance_uncached("ACG", "ATG")

    def test_private_fixed_scorer_is_safe_under_shared_thread_calls(self):
        rng = random.Random(1501)
        pairs = [("".join(rng.choice("ACGTN-") for _ in range(64)),
                  "".join(rng.choice("ACGTN-") for _ in range(64)))
                 for _ in range(64)]
        expected = [independent_distance(*pair) for pair in pairs]
        with ThreadPoolExecutor(max_workers=4) as executor:
            actual = list(executor.map(lambda pair: enhancement._source_repeat_edit_distance_uncached(*pair), pairs))
        self.assertEqual(expected, actual)

    def test_medoid_candidate_alignment_and_score_inputs_match_python_kernel(self):
        fixtures = json.loads((Path(__file__).parent / "data" / "repeat_boundary_repair.json").read_text())
        cases = [(fixture["sequence"], fixture["starts"], fixture["ends"]) for fixture in fixtures]
        for kwargs in ({}, {"terminal_gap": True}, {"internal_indel": True}, {"left": True}):
            source, candidate = synthetic(**kwargs)
            cases.append((source, candidate.list_repeat_starts,
                          [start + len(repeat) for start, repeat in zip(candidate.list_repeat_starts, candidate.list_repeats)]))
        for source, starts, ends in cases:
            enhancement.clear_source_alignment_caches()
            with patch.object(enhancement, "_source_repeat_edit_distance_uncached",
                              side_effect=enhancement._source_repeat_edit_distance_python):
                old = enhancement.build_source_backed_candidate_from_boundaries(source, None, starts, ends)
                old_medoid = enhancement._source_only_alignment_seed(old.list_repeats)
            enhancement.clear_source_alignment_caches()
            new = enhancement.build_source_backed_candidate_from_boundaries(source, None, starts, ends)
            self.assertEqual(old_medoid, enhancement._source_only_alignment_seed(new.list_repeats))
            self.assertEqual(vars(old), vars(new))
            self.assertEqual(candidate_score_fingerprint(old), candidate_score_fingerprint(new))


if __name__ == "__main__":
    unittest.main()
