import random
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from components.refinement.enhancement_search import (
    build_source_backed_candidate_from_boundaries,
)
from components.refinement.extent_trials import (
    _convergent_path,
    _terminal_hypotheses,
    _width_hypothesis,
    annotate_trials,
    recover_embedded,
)
from components.sequence.utils import reverse_complement


class ExtentTrialsTests(unittest.TestCase):
    def setUp(self):
        self.motif = "TACCTTGCAGTTGACGCTAACGTTCAGC"
        rng = random.Random(220)
        self.dna = "".join(rng.choices("ACGT", k=90))
        self.path = []
        for i in range(9):
            if i:
                self.dna += (
                    "ACGT"[i % 4] + "".join(rng.choices("ACGT", k=30)) + "TGCA"[i % 4]
                )
            start = len(self.dna)
            self.dna += self.motif
            self.path.append((start, len(self.dna)))
        self.dna += "".join(rng.choices("ACGT", k=90))
        self.path = tuple(self.path)

    def candidate(self, dna, path):
        return build_source_backed_candidate_from_boundaries(dna, None, *zip(*path))

    def test_two_distinct_observed_extents_agree_in_either_orientation(self):
        for reverse in (False, True):
            dna = reverse_complement(self.dna) if reverse else self.dna

            def oriented(path):
                return (
                    tuple((len(dna) - b, len(dna) - a) for a, b in reversed(path))
                    if reverse
                    else path
                )

            wanted = oriented(self.path)
            old = oriented(self.path[3:6])
            seeds = [oriented(self.path[:7]), oriented(self.path[:8])]
            donors = {p: [0.9, self.candidate(dna, p)] for p in seeds}
            self.assertEqual(_convergent_path(dna, old, donors, 6)[0], wanted)
            for p, row in donors.items():
                self.assertIsNone(_convergent_path(dna, old, {p: row}, 6))

    def test_terminal_trial_only_adds_one_and_preserves_the_core(self):
        old = self.path[:-1]
        proposals = dict(_terminal_hypotheses(self.dna, self.candidate(self.dna, old)))
        self.assertIn(self.path, proposals)
        self.assertEqual(proposals[self.path]["edit_distance"], 0)
        self.assertTrue(all(len(p) == len(old) + 1 for p in proposals))
        short = self.path[:3]
        self.assertEqual(
            list(_terminal_hypotheses(self.dna, self.candidate(self.dna, short))), []
        )

    def test_recovered_fragment_labels_only_the_new_divergent_repeat(self):
        dna = list(self.dna)
        for start, end in self.path[-2:]:
            for i in range(start, end, 3):
                dna[i] = "ACGT"[("ACGT".index(dna[i]) + 1) % 4]
        dna = "".join(dna)
        anchor = self.path[3:6]
        witnesses = {
            g: [0.85, self.candidate(dna, g), {}]
            for g in (self.path[:7], self.path[:8])
        }
        categories = [
            {"scope": [[0.95, self.candidate(dna, anchor), {}]]},
            {"scope": list(witnesses.values())},
            {},
            {},
            {},
        ]

        def assessment(dna, candidates, *args):
            return [([0.8, c, {}], SimpleNamespace(qualified=True)) for c in candidates]

        module = "components.refinement.extent_trials"
        with (
            patch(module + "._convergent_path", return_value=(self.path, witnesses)),
            patch(module + "._assess", side_effect=assessment),
        ):
            categories = recover_embedded(
                dna,
                categories,
                {"param_max_edit_distance": 6},
                None,
                None,
                categories,
                None,
            )
        labels = categories[0]["scope"][0][1].potential_degenerate_repeats
        self.assertEqual(set(labels), {self.path[-1]})
        self.assertEqual(labels[self.path[-1]]["status"], "potential_degenerate_repeat")
        self.assertEqual(labels[self.path[-1]]["location"], "right")
        self.assertGreater(labels[self.path[-1]]["edit_distance"], 6)

    def test_same_membership_preserves_ordinary_reconstruction_and_records_competing_end(
        self,
    ):
        old = self.path[3:6]
        seeds = [
            {"scope": [[0.95, self.candidate(self.dna, old), {}]]},
            {
                "scope": [
                    [0.9, self.candidate(self.dna, self.path[:n]), {}] for n in (7, 8)
                ]
            },
            {},
            {},
            {},
        ]
        x, y = self.path[-1]
        clipped = (*self.path[:-1], (x, y - 1))
        established = self.candidate(self.dna, clipped)
        baseline = [{"scope": [[0.9, established, {}]]}, {}, {}, {}, {}]
        config = {"param_max_edit_distance": 6}
        deferred = {}
        self.assertIsNone(
            recover_embedded(self.dna, seeds, config, None, None, baseline, deferred)
        )
        self.assertIn(self.path, deferred[clipped])
        self.assertFalse(hasattr(established, "extent_trials"))

        def assess(dna, candidates, *args):
            check = SimpleNamespace(qualified=True, evidence={"decision": "passed"})
            return [([0.9, c, {}], check) for c in candidates]

        with patch("components.refinement.extent_trials._assess", side_effect=assess):
            annotate_trials(self.dna, baseline, config, None, None, deferred)
        (proposal,) = established.extent_trials["proposals"]
        self.assertEqual(proposal["method"], "convergent_boundary_alternative")
        self.assertEqual(proposal["repeat_intervals"], self.path)
        self.assertEqual(
            established.extent_trials["unchanged_repeat_intervals"], clipped
        )

    def test_width_trial_compares_both_anchored_edges(self):
        old = list(self.path)
        a, b = old[3]
        old[3] = (a + 2, b)
        self.assertEqual(_width_hypothesis(self.dna, old, self.motif), self.path)
