import random
import unittest
from unittest.mock import patch

from components.components_detection import (
    RejectedFuzzyCandidateTopology,
    TwoRepeatTerminalRepairSeed,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
)
from components.fast_two_repeat_terminal_repair import (
    FLANK_CONTINUATION_RADIUS,
    MINIMUM_TANDEM_REPEAT_DERIVED_SPACER_COVERAGE,
    MINIMUM_TANDEM_REPEAT_DUST_SCORE,
    _maximum_same_length_flank_identity,
    build_source_winner,
    fast_two_repeat_admission_allowed,
    flank_continuation_identities,
    generate_source_winners,
    repeat_derived_spacer_coverage,
    repeat_dust_triplet_score,
)
from components.module_detection import Detection
from components.module_evaluation import build_advanced_fuzzy_search_filter


PARAMETERS = {
    "param_min_avg_repeat_length": 21,
    "param_max_avg_repeat_length": 55,
    "param_min_avg_spacer_length": 18,
    "param_max_avg_spacer_length": 78,
    "param_min_repeats": 3,
    "param_max_identical_spacers": 4,
    "param_max_identical_cluster_spacers": 3,
}
STRUCTURAL_FILTER = build_advanced_fuzzy_search_filter(PARAMETERS)


def target_fixture():
    repeat = "TTAAATCCCACATAGTTCAGATAAAAC"
    terminal = "TTAAATCCCACATAGGGGACTGTACA"
    first_spacer = "CACGCAACCTTCCGGGATGGGTTCGAATGTAAGTTCGGT"
    second_spacer = "TGGGAGCCCTGTCAGAAGTCCCACCTGGCCGCCCAAATGGT"
    flank = (
        "GCGATCGTACGATTCGAGCTAGCATCGATCGTACGATCGTAGCTAGCTAC"
        "GATCGTACGATCGTAGCTACGATCGTACGATCGTAGCTACGAT"
    )
    source = (
        flank + repeat + first_spacer + repeat
        + second_spacer + terminal + flank
    )
    first_start = len(flank)
    second_start = first_start + len(repeat) + len(first_spacer)
    seed = TwoRepeatTerminalRepairSeed(
        (first_start, second_start),
        (first_start + len(repeat), second_start + len(repeat)),
    )
    expected = (
        (first_start, first_start + len(repeat)),
        (second_start, second_start + len(repeat)),
        (
            second_start + len(repeat) + len(second_spacer),
            second_start + len(repeat) + len(second_spacer) + len(terminal),
        ),
    )
    return source, seed, expected


class FakeTwoRepeatCandidate:
    def __init__(self, repeat, spacer, starts):
        self.list_repeats = [repeat, repeat]
        self.list_gaped_repeats = [repeat, repeat]
        self.list_spacers = [spacer]
        self.list_absolute_start = list(starts)

    def _materialize_alignment(self):
        return self


class FastTwoRepeatTerminalRepairTests(unittest.TestCase):
    def test_vectorized_flank_identity_matches_literal_reference(self):
        rng = random.Random(8128)
        alphabet = "ACGTNRY"

        def literal_reference(query, flanks):
            maximum = 0
            for flank in flanks:
                for offset in range(len(flank) - len(query) + 1):
                    maximum = max(maximum, sum(
                        first == second
                        for first, second in zip(
                            query, flank[offset:offset + len(query)]
                        )
                    ))
            return maximum / float(len(query))

        fixtures = [
            ("A", ("", "A")),
            ("ACGT", ("AC", "TTACGTAA")),
            ("NNRY", ("NNNNRY", "RYNN")),
        ]
        fixtures.extend(
            (
                "".join(rng.choice(alphabet) for _ in range(query_length)),
                tuple(
                    "".join(rng.choice(alphabet) for _ in range(flank_length))
                    for flank_length in flank_lengths
                ),
            )
            for query_length, flank_lengths in (
                (rng.randint(1, 80), (rng.randint(0, 250), rng.randint(0, 250)))
                for _ in range(100)
            )
        )
        for query, flanks in fixtures:
            self.assertEqual(
                literal_reference(query, flanks),
                _maximum_same_length_flank_identity(query, flanks),
            )

    def test_compact_seed_validates_and_survives_topology_combination(self):
        seed = TwoRepeatTerminalRepairSeed((10, 70), (37, 97))
        first = RejectedFuzzyCandidateTopology(
            ((10, 97),), (), 2, 1, (seed,)
        )
        combined = RejectedFuzzyCandidateTopology.combine((first, first))

        self.assertEqual((10, 97), seed.source_interval)
        self.assertEqual((seed,), combined.two_repeat_terminal_repair_seeds)

    def test_capture_requires_narrow_sole_minimum_repeat_failure(self):
        repeat = "TTAAATCCCACATAGTTCAGATAAAAC"
        spacer = "CACGCAACCTTCCGGGATGGGTTCGAATGTAAGTTCGGT"
        candidate = FakeTwoRepeatCandidate(
            repeat, spacer, (10, 10 + len(repeat) + len(spacer))
        )

        seed = Detection._two_repeat_terminal_repair_seed(
            candidate,
            STRUCTURAL_FILTER,
            compact_acceptance_impossible=False,
            family_parent=False,
            enabled=True,
        )
        self.assertIsInstance(seed, TwoRepeatTerminalRepairSeed)

        candidate.list_spacers = ["C" * 23]
        self.assertIsNone(Detection._two_repeat_terminal_repair_seed(
            candidate,
            STRUCTURAL_FILTER,
            compact_acceptance_impossible=False,
            family_parent=False,
            enabled=True,
        ))

    def test_source_only_builder_recovers_target_terminal_region(self):
        source, seed, expected = target_fixture()
        record, stats, _background = build_source_winner(
            source, seed, STRUCTURAL_FILTER
        )

        self.assertIsNotNone(record)
        # Region recall is the contract: the source-only boundary winner may
        # represent the terminal repeat up to the frozen +/-2 length radius.
        self.assertEqual(expected[:2], record["repeat_intervals"][:2])
        self.assertEqual(expected[2][0], record["repeat_intervals"][2][0])
        self.assertLessEqual(
            abs(expected[2][1] - record["repeat_intervals"][2][1]), 2
        )
        self.assertEqual(0, stats["cap_exhausted"])
        self.assertLess(
            record["repeat_dust_triplet_score"],
            MINIMUM_TANDEM_REPEAT_DUST_SCORE,
        )
        self.assertLess(
            record["maximum_repeat_derived_spacer_coverage"],
            MINIMUM_TANDEM_REPEAT_DERIVED_SPACER_COVERAGE,
        )

    def test_stage_scoped_tandem_conjunction_flags_minisatellite(self):
        repeats = (
            "AATGGATAATGGATGATGGATAATGGATAATTAATTAATC",
            "AATGGATAATGGATAATTAATTAATCAGGTTGAAGACCTC",
            "AATGGATAATGGATAATTAATTAATCAGGTTGAAGACCTC",
        )
        spacers = (
            "AGGTTGAAGACCTCTCCTGATC",
            "TCCTGATCAATGGATAATGGAT",
        )
        source = (
            repeats[0] + spacers[0] + repeats[1]
            + spacers[1] + repeats[2]
        )
        starts = (
            0,
            len(repeats[0]) + len(spacers[0]),
            len(repeats[0]) + len(spacers[0])
            + len(repeats[1]) + len(spacers[1]),
        )
        ends = tuple(
            start + len(repeat) for start, repeat in zip(starts, repeats)
        )
        candidate = build_source_backed_candidate_from_boundaries(
            source, None, starts, ends
        )

        self.assertGreaterEqual(
            repeat_dust_triplet_score(
                candidate,
                "AATGGATAATGGATAATTAATTAATCAGGTTGAAGACCTC",
            ),
            MINIMUM_TANDEM_REPEAT_DUST_SCORE,
        )
        self.assertGreaterEqual(
            repeat_derived_spacer_coverage(candidate),
            MINIMUM_TANDEM_REPEAT_DERIVED_SPACER_COVERAGE,
        )

    def test_flank_continuation_uses_separate_fixed_source_windows(self):
        repeat = "GTTCACTGCCGTACAGGCAGCTTAGAA"
        spacer = "ACGATCCGTTAGTCGACCTGATTCGATACG"
        prefix = repeat + "T" * 17 + spacer
        left_padding = "A" * (FLANK_CONTINUATION_RADIUS - len(prefix))
        candidate_source = repeat + spacer + repeat + spacer + repeat
        source = left_padding + prefix + candidate_source + "C" * 250
        start = len(left_padding) + len(prefix)
        starts = (
            start,
            start + len(repeat) + len(spacer),
            start + 2 * (len(repeat) + len(spacer)),
        )
        ends = tuple(value + len(repeat) for value in starts)
        candidate = build_source_backed_candidate_from_boundaries(
            source, None, starts, ends
        )

        evidence = flank_continuation_identities(source, candidate)

        self.assertEqual(1.0, evidence["maximum_repeat_flank_identity"])
        self.assertEqual(1.0, evidence["maximum_spacer_flank_identity"])

    def test_admission_policy_is_conjunctive_and_fails_closed(self):
        # Weak continuation preserves the weakest desired sealed-region call.
        self.assertTrue(fast_two_repeat_admission_allowed(
            0.5722424172387197, 0.0, 0.619, 0.512
        ))
        self.assertFalse(fast_two_repeat_admission_allowed(
            0.569999, 0.0, 0.0, 0.0
        ))
        self.assertFalse(fast_two_repeat_admission_allowed(
            0.70, 0.75, 0.0, 0.0
        ))
        self.assertFalse(fast_two_repeat_admission_allowed(
            0.70, 0.0, 0.90, 0.55
        ))
        self.assertTrue(fast_two_repeat_admission_allowed(
            0.70, 0.0, 0.90, 0.549
        ))
        self.assertTrue(fast_two_repeat_admission_allowed(
            0.70, 0.0, 0.899, 0.55
        ))
        self.assertFalse(fast_two_repeat_admission_allowed(
            float("nan"), 0.0, 0.0, 0.0
        ))

    def test_source_seed_cap_fails_closed_before_generation(self):
        source, first, _expected = target_fixture()
        second = TwoRepeatTerminalRepairSeed(
            (first.repeat_starts[0] + 1, first.repeat_starts[1] + 1),
            (first.repeat_ends[0] + 1, first.repeat_ends[1] + 1),
        )
        with patch(
                "components.fast_two_repeat_terminal_repair."
                "MAXIMUM_SEEDS_PER_SOURCE",
                1):
            winners, stats = generate_source_winners(
                source, (first, second), STRUCTURAL_FILTER
            )

        self.assertEqual((), winners)
        self.assertEqual(1, stats["source_seed_cap_exhausted"])


if __name__ == "__main__":
    unittest.main()
