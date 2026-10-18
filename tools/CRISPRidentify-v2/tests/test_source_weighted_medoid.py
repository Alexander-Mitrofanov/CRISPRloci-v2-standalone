from collections import Counter
import itertools
import random
import unittest

from components import components_evaluated_arrays_enhancement as enhancement
from tests.test_source_unit_edit_kernel import independent_distance


def full_weighted_reference(repeats):
    counts = Counter(repeats)
    return min(counts, key=lambda repeat: (
        sum(count * independent_distance(repeat, other) for other, count in counts.items()),
        -counts[repeat], -len(repeat), repeat,
    ))


class SourceWeightedMedoidTests(unittest.TestCase):
    def test_exhaustive_small_multisets_keep_original_weighted_ties(self):
        choices = ("A", "C", "AA", "AC", "CA", "N-")
        for count in range(1, 6):
            for repeats in itertools.combinations_with_replacement(choices, count):
                expected = full_weighted_reference(repeats)
                self.assertEqual(expected, enhancement._source_only_alignment_seed(repeats))
                self.assertEqual(expected, enhancement._source_only_alignment_seed(tuple(reversed(repeats))))

    def test_seeded_weighted_indel_ambiguity_and_order_equivalence(self):
        rng = random.Random(2026092605)
        for _ in range(300):
            choices = ["".join(rng.choice("ACGTN-?") for _ in range(rng.randrange(1, 20)))
                       for _ in range(rng.randrange(1, 9))]
            repeats = [repeat for repeat in choices for _ in range(rng.randrange(1, 8))]
            expected = full_weighted_reference(repeats)
            self.assertEqual(expected, enhancement._source_only_alignment_seed(repeats))
            rng.shuffle(repeats)
            self.assertEqual(expected, enhancement._source_only_alignment_seed(repeats))

    def test_tie_breakers_remain_multiplicity_length_then_lexical(self):
        for repeats, expected in ((("A", "AA"), "AA"), (("AC", "CA"), "AC"),
                                   (("A", "A", "C"), "A"),
                                   (("ACGT", "ACGT", "ACGA", "TCGA"), "ACGT")):
            self.assertEqual(expected, full_weighted_reference(repeats))
            self.assertEqual(expected, enhancement._source_only_alignment_seed(repeats))


if __name__ == "__main__":
    unittest.main()
