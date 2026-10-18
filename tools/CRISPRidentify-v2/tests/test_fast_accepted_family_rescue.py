import copy
import unittest
from unittest.mock import patch

from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries,
)
from components.fast_accepted_family_rescue import (
    ACCEPTED_FAMILY_RESCUE_POLICY,
    EXACT_TEMPLATE_SATELLITE,
    FAMILY_TERMINAL,
    RAW_TWO_REPEAT_SEED_PARENT,
    classification_basis,
    completed_score_is_admissible,
    generate_accepted_family_rescue_proposals,
    replay_classification_evidence,
    select_exact_template_raw_seed_parents,
)
from components.components_detection import TwoRepeatTerminalRepairSeed
from components.module_evaluation import build_advanced_fuzzy_search_filter


PARAMETERS = {
    "param_min_repeats": 3,
    "param_min_avg_repeat_length": 21,
    "param_max_avg_repeat_length": 55,
    "param_min_avg_spacer_length": 18,
    "param_max_avg_spacer_length": 78,
    "param_max_identical_spacers": 4,
    "param_max_identical_cluster_spacers": 3,
}


def append_array(parts, repeats, spacers):
    starts = []
    ends = []
    for index, repeat in enumerate(repeats):
        starts.append(sum(map(len, parts)))
        parts.append(repeat)
        ends.append(sum(map(len, parts)))
        if index < len(spacers):
            parts.append(spacers[index])
    return tuple(starts), tuple(ends)


def candidate(source, starts, ends):
    return build_source_backed_candidate_from_boundaries(
        source, None, starts, ends
    )


class FastAcceptedFamilyRescueTests(unittest.TestCase):
    def setUp(self):
        self.structural = build_advanced_fuzzy_search_filter(PARAMETERS)

    def test_family_phase_and_terminal_requires_strict_source_intact(self):
        family = "GTTTCAATACGCCCTCTGGCGTGACTCCCGTATCGAC"
        target_repeats = (
            "GTTTCAATAAGCCCTCCGCCGCGACTCCCGCATCGAC",
            "GTTGCAATACACCCTCCGGCGTGACTCCCGTATCGAC",
            "GTTTCAATACGCCCTCCGCCGTGACTCCCGTATCATC",
            "GTTTCAATACGCCCTCCGGCGTGGCTCCCGTATCGAC",
            "ATTTCAATGCGCCCACCGGCGTGACTCCCGTATCGAC",
            "GTTTCAATACGCCCTCTGGCGCGACTCCCATATCGAC",
        )
        target_spacers = (
            "GTGCCGCTAACGCGATCATCACCTTCGCGCAGAACCCCTTCTGAC",
            "TCTTTGCAATACCCGGTGACTATTCCCATTTTGAGTGGT",
            "GTTTAAGCAAAGCTCGCCGTGACCAAATCGCCATCTA",
            "TTCTTGCTCCACTCGATCGACATGCTCTGATCTGTGA",
            "GTTGCTGAAGTGGGTTCAGCATGCTGACGGGTCTTTCGAT",
        )
        parts = ["ACGT" * 25]
        anchor_starts, anchor_ends = append_array(
            parts,
            (family, family, family),
            ("AACCGGTTAACCGGTTAACCGGTTAACCGGTT", "TTGGAACCTTGGAACCTTGGAACCTTGGAACC"),
        )
        parts.append("ACGT" * 375)
        full_starts, full_ends = append_array(parts, target_repeats, target_spacers)
        parts.append("TGCA" * 100)
        source = "".join(parts)
        anchor = candidate(source, anchor_starts, anchor_ends)
        parent_starts = full_starts[:-1]
        parent_ends = tuple(start + 23 for start in parent_starts)
        parent = candidate(source, parent_starts, parent_ends)
        categories = [{}, {}, {}, {}, {}]
        categories[0][(anchor_starts[0], anchor_ends[-1])] = [[0.8, anchor, {}]]
        categories[4][(parent_starts[0], parent_ends[-1])] = [[0.21, parent, {}]]

        proposals, stats = generate_accepted_family_rescue_proposals(
            source, categories, self.structural
        )
        self.assertEqual(1, len(proposals))
        proposal = proposals[0]
        self.assertEqual(FAMILY_TERMINAL, proposal.arm)
        self.assertEqual(full_starts, tuple(proposal.candidate.list_repeat_starts))
        self.assertEqual(
            ("spacer_length_cv_above_maximum",),
            tuple(reason.code for reason in proposal.base_source_intact_decision.reasons),
        )
        self.assertTrue(proposal.completed_source_intact_decision.eligible)
        self.assertEqual(1, stats["family_terminal_proposals"])
        self.assertEqual(1, stats["source_intact_context_constructions"])
        basis = classification_basis(source, proposal, 0.329)
        self.assertEqual("source_backed_accepted_family_rescue", basis["method"])
        self.assertEqual(FAMILY_TERMINAL, basis["accepted_family_rescue"]["arm"])
        replay_classification_evidence(
            source,
            proposal.candidate,
            basis["accepted_family_rescue"],
            0.329,
        )
        stale = copy.deepcopy(basis["accepted_family_rescue"])
        stale["terminal_start"] += 1
        with self.assertRaisesRegex(ValueError, "terminal start"):
            replay_classification_evidence(
                source, proposal.candidate, stale, 0.329
            )

    def test_exact_template_arm_uses_independent_parent_offsets(self):
        template = "ATTGAACCTTAACATAGGATGTATTGAAAT"
        parts = ["TGCA" * 25]
        anchors = []
        anchor_spacers = (
            "ACGTCAGTACGATCGTACCTGATCGTACGTA",
            "TTGACCGATGCTAGTCGATGACCTAGTCGAA",
        )
        for _ in range(4):
            starts, ends = append_array(
                parts,
                (template,) * 8,
                tuple(anchor_spacers[index % 2] for index in range(7)),
            )
            anchors.append((starts, ends))
            parts.append("CGTA" * 30)
        parts.append("ACGT" * 100)
        parent_start = sum(map(len, parts))
        parent_repeat_one = "ATTTATCATTGAACCTTAACATAGGATGTATTGAAAT"
        parent_repeat_two = "ATTTTCATTGAACCTTAACATAGGATGTATTGAAAT"
        parent_spacer = "CTTTGCATCCCTGCCCCTCTTACCTGAGC"
        second_spacer = "TATTTTGCTCTTAAATTTGCAAGGTCTGTATTAT"
        terminal = "ATTGAACCTTAACATAGTATATACTAAACC"
        parts.extend((parent_repeat_one, parent_spacer, parent_repeat_two, second_spacer, terminal))
        parts.append("TGCA" * 100)
        source = "".join(parts)

        categories = [{}, {}, {}, {}, {}]
        for starts, ends in anchors:
            accepted = candidate(source, starts, ends)
            categories[2][(starts[0], ends[-1])] = [[0.7, accepted, {}]]
        parent_starts = (
            parent_start,
            parent_start + len(parent_repeat_one) + len(parent_spacer),
        )
        parent_ends = (
            parent_starts[0] + len(parent_repeat_one),
            parent_starts[1] + len(parent_repeat_two),
        )
        raw_seed = TwoRepeatTerminalRepairSeed(
            parent_starts, parent_ends, 7
        )
        raw_parents, raw_stats = select_exact_template_raw_seed_parents(
            source, categories, (raw_seed,), self.structural
        )
        self.assertEqual(1, len(raw_parents))
        self.assertEqual(1, raw_stats["raw_seed_parent_matches"])
        self.assertEqual((7,), raw_parents[0].cluster_indices)

        proposals, stats = generate_accepted_family_rescue_proposals(
            source,
            categories,
            self.structural,
            raw_seed_parent_rows=((
                [0.371, raw_parents[0].candidate, {}],
                raw_parents[0],
            ),),
        )
        self.assertEqual(1, len(proposals))
        proposal = proposals[0]
        self.assertEqual(EXACT_TEMPLATE_SATELLITE, proposal.arm)
        self.assertEqual((7, 6), proposal.source_metrics["template_offsets"])
        self.assertEqual(4, proposal.source_metrics["accepted_family_array_count"])
        self.assertEqual(3, len(proposal.candidate.list_repeats))
        self.assertEqual(RAW_TWO_REPEAT_SEED_PARENT, proposal.parent_origin)
        self.assertEqual((7,), proposal.raw_seed_cluster_indices)
        self.assertEqual(1, stats["exact_template_proposals"])
        basis = classification_basis(source, proposal, 0.399)
        replay_classification_evidence(
            source,
            proposal.candidate,
            basis["accepted_family_rescue"],
            0.399,
        )

    def test_exact_template_arm_requires_three_accepted_arrays(self):
        template = "ATTGAACCTTAACATAGGATGTATTGAAAT"
        parts = ["ACGT" * 20]
        anchors = []
        for _ in range(2):
            starts, ends = append_array(
                parts,
                (template,) * 8,
                tuple(
                    ("ACGT" * 8) if index % 2 == 0 else ("TGCA" * 8)
                    for index in range(7)
                ),
            )
            anchors.append((starts, ends))
            parts.append("GTCA" * 20)
        parent_start = sum(map(len, parts))
        first = "AAAAAAA" + template
        second = "CCCCCC" + template
        parts.extend((first, "GATC" * 8, second, "CTAG" * 9, template))
        source = "".join(parts)
        categories = [{}, {}, {}, {}, {}]
        for starts, ends in anchors:
            accepted = candidate(source, starts, ends)
            categories[2][(starts[0], ends[-1])] = [[0.8, accepted, {}]]
        starts = (parent_start, parent_start + len(first) + 32)
        ends = (starts[0] + len(first), starts[1] + len(second))
        parent = candidate(source, starts, ends)
        categories[4][(starts[0], ends[-1])] = [[0.4, parent, {}]]
        proposals, _stats = generate_accepted_family_rescue_proposals(
            source, categories, self.structural
        )
        self.assertEqual((), proposals)

    def test_raw_seed_global_family_absence_short_circuits_safely(self):
        repeat = "ATTGAACCTTAACATAGGATGTATTGAAAT"
        spacer = "CTTTGCATCCCTGCCCCTCTTACCTGAGC"
        source = "A" * 100 + repeat + spacer + repeat + "C" * 100
        first = 100
        second = first + len(repeat) + len(spacer)
        seed = TwoRepeatTerminalRepairSeed(
            (first, second),
            (first + len(repeat), second + len(repeat)),
            3,
        )

        parents, stats = select_exact_template_raw_seed_parents(
            source, [{}, {}, {}, {}, {}], (seed,), self.structural
        )

        self.assertEqual((), parents)
        self.assertEqual(1, stats["raw_seed_global_exact_family_impossible"])
        self.assertEqual(1, stats["raw_seed_unique_uncategorized_geometries"])
        self.assertEqual(1, stats["raw_seed_global_family_skips"])
        self.assertEqual(1, stats["raw_seed_skipped_geometry_count"])
        self.assertEqual(64, len(stats["raw_seed_skipped_geometry_sha256"]))
        self.assertEqual(
            64, len(stats["raw_seed_skipped_clustered_geometry_sha256"])
        )
        self.assertEqual(0, stats["raw_seed_prefilter_terminal_windows_tested"])
        self.assertEqual(0, stats["raw_seed_parent_matches"])
        self.assertEqual(
            "complete", stats["raw_seed_parent_match_count_status"]
        )

    def test_raw_seed_parent_cap_accepts_eight_and_stops_on_ninth(self):
        template = "ATTGAACCTTAACATAGGATGTATTGAAAT"
        spacer = "CTTTGCATCCCTGCCCCTCTTACCTGAGC"
        parts = ["A" * 100]
        anchors = []
        for _index in range(3):
            starts, ends = append_array(
                parts,
                (template,) * 8,
                tuple(("ACGT" * 8) for _ in range(7)),
            )
            anchors.append((starts, ends))
            parts.append("C" * 100)
        raw_intervals = []
        for _index in range(10):
            first = sum(map(len, parts))
            parts.extend((template, spacer, template, "G" * 100))
            second = first + len(template) + len(spacer)
            raw_intervals.append((
                (first, second),
                (first + len(template), second + len(template)),
            ))
        source = "".join(parts)
        categories = [{}, {}, {}, {}, {}]
        for starts, ends in anchors:
            anchor = candidate(source, starts, ends)
            categories[2][(starts[0], ends[-1])] = [[0.8, anchor, {}]]
        seeds = tuple(
            TwoRepeatTerminalRepairSeed(starts, ends, index)
            for index, (starts, ends) in enumerate(raw_intervals)
        )
        self.assertEqual(
            8,
            ACCEPTED_FAMILY_RESCUE_POLICY.
            maximum_raw_seed_parent_matches_per_source,
        )

        def select(selected_seeds):
            with patch(
                "components.fast_accepted_family_rescue."
                "_template_satellite_proposals",
                return_value=(object(),),
            ) as proposal_builder, patch(
                "components.fast_accepted_family_rescue."
                "_select_unique_source_winner",
                return_value=object(),
            ):
                parents, stats = select_exact_template_raw_seed_parents(
                    source, categories, selected_seeds, self.structural
                )
            return parents, stats, proposal_builder.call_count

        eight, eight_stats, eight_calls = select(seeds[:8])
        self.assertEqual(8, len(eight))
        self.assertEqual(8, eight_calls)
        self.assertEqual(8, eight_stats["raw_seed_parent_matches"])
        self.assertEqual(
            "complete", eight_stats["raw_seed_parent_match_count_status"]
        )
        self.assertNotIn(
            "raw_seed_parent_match_cap_exhausted", eight_stats
        )

        parents, stats, proposal_calls = select(seeds)
        self.assertEqual((), parents)
        self.assertEqual(9, proposal_calls)
        self.assertEqual(1, stats["raw_seed_parent_match_cap_exhausted"])
        self.assertEqual(9, stats["raw_seed_parent_matches_at_cap"])
        self.assertEqual(1, stats["raw_seed_parent_matches_is_lower_bound"])
        self.assertEqual(
            "lower_bound_cap_exhausted",
            stats["raw_seed_parent_match_count_status"],
        )
        self.assertEqual(1, stats["raw_seed_unexamined_geometry_count"])
        self.assertEqual(64, len(stats["raw_seed_unexamined_geometry_sha256"]))
        self.assertEqual(
            64, len(stats["raw_seed_unexamined_clustered_geometry_sha256"])
        )

        reversed_parents, reversed_stats, reversed_calls = select(
            tuple(reversed(seeds))
        )
        self.assertEqual((), reversed_parents)
        self.assertEqual(9, reversed_calls)
        self.assertEqual(
            stats["raw_seed_unexamined_geometry_sha256"],
            reversed_stats["raw_seed_unexamined_geometry_sha256"],
        )
        self.assertEqual(
            stats["raw_seed_unexamined_clustered_geometry_sha256"],
            reversed_stats["raw_seed_unexamined_clustered_geometry_sha256"],
        )

    def test_completed_score_band_fails_closed(self):
        self.assertFalse(completed_score_is_admissible(0.299999))
        self.assertTrue(completed_score_is_admissible(0.30))
        self.assertTrue(completed_score_is_admissible(0.499999))
        self.assertFalse(completed_score_is_admissible(0.50))
        self.assertFalse(completed_score_is_admissible(float("nan")))


if __name__ == "__main__":
    unittest.main()
