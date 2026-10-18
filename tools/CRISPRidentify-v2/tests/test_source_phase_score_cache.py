import math
import random
import struct
import unittest

from components.components_evaluated_arrays_enhancement import (
    _CANONICAL_DNA_BASES,
    _canonical_base_counts,
    _jeffreys_column_log_bayes_factor,
    _jeffreys_column_log_bayes_factor_from_counts,
)


def _legacy_score(observations, background_probabilities):
    counts = {base: 0 for base in _CANONICAL_DNA_BASES}
    for character in observations:
        base = character.upper()
        if base in counts:
            counts[base] += 1
    observation_count = sum(counts.values())
    if observation_count < 2:
        return 0.0
    log_column_marginal = (
        math.lgamma(2.0)
        - math.lgamma(observation_count + 2.0)
        + sum(
            math.lgamma(counts[base] + 0.5) - math.lgamma(0.5)
            for base in _CANONICAL_DNA_BASES
        )
    )
    log_background_likelihood = sum(
        counts[base] * math.log(background_probabilities[base])
        for base in _CANONICAL_DNA_BASES
    )
    return log_column_marginal - log_background_likelihood


class SourcePhaseScoreCacheTests(unittest.TestCase):
    def setUp(self):
        _jeffreys_column_log_bayes_factor_from_counts.cache_clear()

    def test_cached_wrapper_is_float_bit_exact_to_legacy_order(self):
        generator = random.Random(20260823)
        backgrounds = (
            {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
            {"A": 0.11, "C": 0.19, "G": 0.29, "T": 0.41},
        )
        observations = ["", "A", "ac", "ACGT", "NN--acgt"]
        observations.extend(
            "".join(
                generator.choice("ACGTN-acgtn")
                for _ in range(generator.randint(0, 80))
            )
            for _ in range(512)
        )
        for background in backgrounds:
            for sequence in observations:
                expected = _legacy_score(iter(sequence), background)
                observed = _jeffreys_column_log_bayes_factor(
                    iter(sequence), background
                )
                self.assertEqual(
                    struct.pack(">d", expected), struct.pack(">d", observed)
                )

    def test_low_evidence_keeps_legacy_background_short_circuit(self):
        self.assertEqual(
            0.0, _jeffreys_column_log_bayes_factor(iter("N-A"), {})
        )

    def test_canonical_base_count_mapping_keeps_public_shape(self):
        self.assertEqual(
            {"A": 2, "C": 2, "G": 2, "T": 2},
            _canonical_base_counts(iter("AaCcGgTtN-")),
        )

    def test_repeated_exact_key_hits_and_cache_stays_bounded(self):
        background = (0.11, 0.19, 0.29, 0.41)
        key = ((7, 3, 2, 1), background)
        first = _jeffreys_column_log_bayes_factor_from_counts(*key)
        second = _jeffreys_column_log_bayes_factor_from_counts(*key)
        self.assertEqual(struct.pack(">d", first), struct.pack(">d", second))
        self.assertEqual(
            1,
            _jeffreys_column_log_bayes_factor_from_counts.cache_info().hits,
        )

        for index in range(5000):
            _jeffreys_column_log_bayes_factor_from_counts(
                (index + 2, index % 5, index % 7, index % 11), background
            )
        info = _jeffreys_column_log_bayes_factor_from_counts.cache_info()
        self.assertEqual(4096, info.maxsize)
        self.assertEqual(4096, info.currsize)


if __name__ == "__main__":
    unittest.main()
