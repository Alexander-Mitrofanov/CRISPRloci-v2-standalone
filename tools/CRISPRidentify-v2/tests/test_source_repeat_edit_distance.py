import itertools
import random
import unittest
from unittest.mock import patch

from components import components_evaluated_arrays_enhancement as enhancement


def _legacy_distance(first, second):
    return sum(
        len(changes)
        for changes in enhancement._minimal_edit_changes(first, second)
    )


class TestSourceRepeatEditDistance(unittest.TestCase):
    def test_empty_and_unequal_length_inputs(self):
        cases = (
            ("", "", 0),
            ("", "ACG", 3),
            ("ACGT", "", 4),
            ("ACGT", "AG", 2),
            ("AAAA", "TTTTTT", 6),
            ("GATTACA", "GCATGCU", 4),
        )
        for first, second, expected in cases:
            with self.subTest(first=first, second=second):
                self.assertEqual(
                    enhancement._source_repeat_edit_distance(first, second),
                    expected,
                )
                self.assertEqual(
                    enhancement._source_repeat_edit_distance(second, first),
                    expected,
                )

    def test_ambiguous_characters_retain_literal_equality_semantics(self):
        cases = (
            ("ACGTNRY-", "ACGTNRY-"),
            ("ANRY-", "ANYY?"),
            ("NNNN", "N-NN"),
            ("acgtn", "ACGTN"),
        )
        for first, second in cases:
            with self.subTest(first=first, second=second):
                self.assertEqual(
                    enhancement._source_repeat_edit_distance(first, second),
                    _legacy_distance(first, second),
                )

    def test_randomized_equivalence_to_legacy_dynamic_programming(self):
        generator = random.Random(20260824)
        alphabet = "ACGTNRYKMSWBDHV-?acgt"
        for case_index in range(2000):
            first = "".join(
                generator.choice(alphabet)
                for _ in range(generator.randrange(56))
            )
            second = "".join(
                generator.choice(alphabet)
                for _ in range(generator.randrange(56))
            )
            with self.subTest(case_index=case_index):
                self.assertEqual(
                    enhancement._source_repeat_edit_distance(first, second),
                    _legacy_distance(first, second),
                )

    def test_long_inputs_use_exact_legacy_fallback(self):
        length = enhancement._SOURCE_REPEAT_BITPARALLEL_MAX_LENGTH + 1
        first = "A" * length + "NRY"
        second = "A" * (length - 2) + "CGTNRY"
        expected = _legacy_distance(first, second)

        original = enhancement._minimal_edit_changes
        with patch.object(
                enhancement, "_minimal_edit_changes", wraps=original) as legacy:
            self.assertEqual(
                enhancement._source_repeat_edit_distance(first, second),
                expected,
            )
        legacy.assert_called_once_with(first, second)

    def test_medoid_is_deterministic_and_uses_short_fast_path(self):
        repeats = ("ACGT", "ACGT", "ACGA", "TCGA")
        expected = "ACGT"
        for ordering in itertools.permutations(repeats):
            with patch.object(
                    enhancement,
                    "_minimal_edit_changes",
                    side_effect=AssertionError("short medoid used legacy DP"),
            ):
                self.assertEqual(
                    enhancement._source_only_alignment_seed(ordering),
                    expected,
                )


if __name__ == "__main__":
    unittest.main()
