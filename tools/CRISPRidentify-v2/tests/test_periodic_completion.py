import random
import unittest

from components.refinement.boundary_panel import configured_consensuses, oriented_consensuses, resolve_phase
from components.refinement.enhancement_search import (
    build_source_backed_candidate_from_boundaries,
)
from components.refinement.periodic_completion import _placement, plan_path, recurrent_allele_trim
from components.refinement.periodic_transaction import _compatible_anchors, _embedded_fragment_donor
from components.sequence.utils import reverse_complement


class PeriodicCompletionTests(unittest.TestCase):
    motif = "TACCTTGCAGTTGACGCTAACGTTCAGC"

    def sample(self, copies, gaps):
        rng = random.Random(91)
        dna = ""
        path = []
        for i, seq in enumerate(copies):
            if i:
                dna += "".join(rng.choices("ACGT", k=gaps[i - 1]))
            path.append((len(dna), len(dna) + len(seq)))
            dna += seq
        return dna, tuple(path)

    def core(self, dna, path):
        return build_source_backed_candidate_from_boundaries(dna, None, *zip(*path))

    def test_internal_bridges_and_single_weak_terminal_in_both_orientations(self):
        variant = "A" + self.motif[1:]
        weak = "CGATACG" + self.motif[7:]
        dna, full = self.sample(
            [
                self.motif,
                variant,
                self.motif,
                self.motif,
                variant,
                self.motif,
                self.motif,
                weak,
            ],
            [32] * 7,
        )
        retained = [full[i] for i in (0, 2, 3, 5, 6)]
        self.assertEqual(full, plan_path(dna, self.core(dna, retained), 6))

        def reverse(path):
            return tuple((len(dna) - b, len(dna) - a) for a, b in reversed(path))

        rc = reverse_complement(dna)
        self.assertEqual(
            reverse(full), plan_path(rc, self.core(rc, reverse(retained)), 6)
        )

    def test_fixed_width_extension_requires_positional_identity(self):
        # Two opposing indels give a cheap edit alignment but displace most
        # source columns. The fixed-width pass must not assert both ends.
        shifted = self.motif[:10] + self.motif[12:] + "GG"
        dna = "A" * 32 + shifted
        self.assertIsNone(_placement(dna, self.motif, 0, None, 32, 0, 6))

    def test_nested_fragment_requires_one_regular_spanning_path(self):
        dna, full = self.sample([self.motif] * 14, [34] * 13)
        inner = [full[i] for i in (3, 5, 7)]
        extension = "G" * 16
        for _, b in inner:
            dna = dna[:b] + extension + dna[b + len(extension):]
        fragment = self.core(dna, [(a, b + len(extension)) for a, b in inner])
        row = [0.8, self.core(dna, full), []]
        self.assertEqual(_embedded_fragment_donor(dna, fragment, [row]), (full, row))
        suffix = [0.8, self.core(dna, full[1:]), []]
        self.assertEqual(_embedded_fragment_donor(dna, fragment, [row, suffix]), (full, row))
        rival = list(full[1:])
        a, b = rival[-2]
        rival[-2] = (a + 1, b + 1)
        self.assertIsNone(_embedded_fragment_donor(dna, fragment, [row, [0.8, self.core(dna, rival), []]]))
        self.assertIsNone(_embedded_fragment_donor(dna, fragment, [[0.7, row[1], []]]))
        self.assertIsNone(_embedded_fragment_donor(dna, fragment, [[*row, {"tags": ["tandem"]}]]))
        row[1].source_path_classification_eligible = False
        self.assertIsNone(_embedded_fragment_donor(dna, fragment, [row]))

    def test_single_base_trim_requires_independent_whole_alleles_in_both_orientations(self):
        variant = "C" + self.motif[1:9] + "G" + self.motif[10:]
        dna, path = self.sample([self.motif, variant] * 4, [32] * 7)
        old = tuple((a - (i == 3), b) for i, (a, b) in enumerate(path))
        support = recurrent_allele_trim(dna, old, self.core(dna, path), .6)
        self.assertEqual(support["trimmed_side"], "left")
        self.assertEqual(support["unchanged_full_allele_counts"][variant], 3)
        def reverse(p):
            return tuple((len(dna) - b, len(dna) - a) for a, b in reversed(p))

        rc = reverse_complement(dna)
        support = recurrent_allele_trim(rc, reverse(old), self.core(rc, reverse(path)), .6)
        self.assertEqual(support["trimmed_side"], "right")
        for seq in (variant[:12] + "C" + variant[13:], "T" + variant[1:]):
            altered = dna[:path[3][0]] + seq + dna[path[3][1]:]
            self.assertIsNone(recurrent_allele_trim(altered, old, self.core(altered, path), .6))

    def test_release_panel_and_source_only_override(self):
        motifs = configured_consensuses(None)
        self.assertEqual(len(motifs), 56)
        self.assertEqual(set(motifs), {reverse_complement(s) for s in motifs})
        self.assertEqual(configured_consensuses("none"), [])
        self.assertEqual(configured_consensuses(self.motif), oriented_consensuses(self.motif))

    def test_phase_requires_exemplars_and_abstains_on_competing_phases(self):
        dna, full = self.sample([self.motif] * 5, [32] * 4)
        clipped = tuple((a + 1, b) for a, b in full)
        motifs = oriented_consensuses(reverse_complement(self.motif))
        self.assertEqual((full, [self.motif]), resolve_phase(dna, clipped, motifs))
        self.assertIsNone(resolve_phase(dna, clipped, motifs + [self.motif[1:]]))
        self.assertIsNone(resolve_phase(dna, clipped[:1], motifs))
        with self.assertRaises(ValueError):
            oriented_consensuses(self.motif + ",N")

    def test_extending_to_a_separate_core_requires_a_spanning_witness(self):
        dna, path = self.sample([self.motif] * 8, [32] * 7)
        first, second = path[:4], path[4:]
        left, right = self.core(dna, first), self.core(dna, second)
        a, b = (0, 100), (101, 200)
        cats = [{a: [[0.9, left, []]], b: [[0.9, right, []]]}, {}, {}, {}, {}]
        slots = [((0, a, 0), cats[0][a][0], first), ((0, b, 0), cats[0][b][0], second)]
        self.assertFalse(_compatible_anchors(path, first, left, slots, cats))
        cats[1][a] = [[0.8, self.core(dna, path), []]]
        self.assertTrue(_compatible_anchors(path, first, left, slots, cats))
