import unittest
from unittest.mock import patch

from components.components_detection_refinement import CrisprCandidate
from components.components_evaluated_arrays_enhancement import (
    collapse_source_backed_phase_equivalent_candidates,
)


class MachinePhaseScoreGateTests(unittest.TestCase):
    @staticmethod
    def _candidate(dna, starts, ends):
        repeats = [dna[start:end] for start, end in zip(starts, ends)]
        spacers = [
            dna[end:next_start]
            for end, next_start in zip(ends, starts[1:])
        ]
        return CrisprCandidate(
            list_repeats=repeats,
            list_repeats_gaped=repeats,
            list_spacers=spacers,
            list_repeat_starts=list(starts),
        )

    @classmethod
    def _weak_boundary_family(cls):
        repeat = "C" + "G" * 27 + "T"
        spacers = (
            "T" + "A" * 31,
            "T" + "C" * 31,
            "T" + "G" * 31,
            "T" * 32,
        )
        dna = "ACGT" * 1250
        true_starts = []
        for repeat_index in range(5):
            true_starts.append(len(dna))
            dna += repeat
            if repeat_index < len(spacers):
                dna += spacers[repeat_index]
        dna += "T" + "ACGT" * 1250

        true_starts = tuple(true_starts)
        true_ends = tuple(start + len(repeat) for start in true_starts)
        wrong = cls._candidate(
            dna,
            tuple(start + 1 for start in true_starts),
            tuple(end - 1 for end in true_ends),
        )
        correct = cls._candidate(dna, true_starts, true_ends)
        return dna, [wrong, correct]

    def test_all_boundary_rejections_skip_scores_and_preserve_candidates(self):
        dna, candidates = self._weak_boundary_family()
        with patch(
                "components.components_evaluated_arrays_enhancement."
                "_source_phase_geometry_score") as source_score:
            collapsed, provenance, stats = (
                collapse_source_backed_phase_equivalent_candidates(
                    dna,
                    {(4900, 5400): candidates},
                    skip_source_scores_when_boundaries_reject=True,
                )
            )

        source_score.assert_not_called()
        self.assertEqual(candidates, collapsed[(4900, 5400)])
        self.assertEqual(0, stats["groups_collapsed"])
        self.assertEqual(1, stats["groups_weak_boundary_evidence"])
        self.assertEqual(0, stats["candidates_removed"])
        self.assertEqual("unchanged", provenance[0]["status"])
        self.assertEqual(
            "weak_source_boundary_evidence", provenance[0]["reason"]
        )
        self.assertTrue(provenance[0]["source_scoring_skipped"])
        self.assertTrue(all(
            not evidence["passes"]
            for _geometry, evidence
            in provenance[0]["boundary_evidence_by_geometry"]
        ))

    def test_gate_requires_boolean(self):
        with self.assertRaisesRegex(ValueError, "score gate must be boolean"):
            collapse_source_backed_phase_equivalent_candidates(
                "A", {}, skip_source_scores_when_boundaries_reject=1
            )


if __name__ == "__main__":
    unittest.main()
