"""Competing one-copy, conserved-block and retained-allele interpretations."""

import random
import unittest
from types import SimpleNamespace

from components.refinement.reviewed_boundary_selection import (
    rank_alignment,
    select,
    select_alignment,
)
from components.sequence.utils import reverse_complement


def sequence(seed, length):
    random_source = random.Random(seed)
    return "".join(random_source.choice("ACGT") for _ in range(length))


def array(repeats, spacers):
    source = "ATCGACTA"
    path = []
    for i, repeat in enumerate(repeats):
        path.append((len(source), len(source) + len(repeat)))
        source += repeat
        if i < len(spacers):
            source += spacers[i]
    return source + "GACTACAT", tuple(path)


def reflect(source, path):
    return tuple((len(source) - b, len(source) - a) for a, b in reversed(path))


class BoundarySelectionTests(unittest.TestCase):
    def check_orientations(self, source, primary, retained, expected, rule):
        for reverse in (False, True):
            if reverse:
                retained = {reflect(source, p): s for p, s in retained.items()}
                primary = reflect(source, primary)
                expected = reflect(source, expected)
                source = source[::-1].translate(str.maketrans("ATCG", "TAGC"))
            choice = select(source, primary, retained)
            self.assertIsNotNone(choice)
            self.assertEqual(choice[:2], (expected, rule))

    def test_distinct_spacers_and_nearby_source_cluster_support_one_longer_edge(self):
        core = sequence(3, 27)
        motif = "GACTGA" + core
        long = ["GATCAA" + core, motif, motif, motif]
        spacers = [sequence(i, 26) + end for i, end in enumerate("ATC")]
        local, p = array(long, spacers)
        primary = tuple((a + 6, b) for a, b in p)
        self.assertIsNone(select(local, primary, {primary: 0.59}))
        distant, _ = array([motif] * 3, [sequence(100, 26), sequence(101, 27)])
        source = local + "A" * 150 + distant
        self.check_orientations(
            source,
            primary,
            {primary: 0.59},
            p,
            "conserved_flank_with_adjacent_exemplars",
        )
        # Repeating the same inter-copy sequence supplies poor tandem controls.
        tandem, _ = array([motif] * 3, [sequence(100, 26)] * 2)
        self.assertIsNone(select(local + tandem, primary, {primary: 0.59}))

    def test_internal_length_outlier_near_edge_keeps_both_mutation_models(self):
        tip = "TCAG" + sequence(71, 22)
        source, p = array(
            [tip, tip[:2] + "A" + tip[2:], tip], [sequence(300, 25), sequence(400, 30)]
        )
        a, b = p[1]
        expected = (p[0], (a + 1, b), p[2])
        self.assertIsNone(select(source, p, {p: 0.65}))
        outlier = source[a:b]
        aligned = tip[:2] + "-" + tip[2:]
        parent = SimpleNamespace(
            consensus=outlier, list_gaped_repeats=(aligned, outlier, aligned)
        )
        chosen = select_alignment(source, p, parent)
        self.assertEqual(chosen[:2], (p, "modal_insertion_alignment"))
        self.assertEqual(
            chosen[2]["alternative_repeat_intervals"], list(map(list, expected))
        )
        self.assertEqual(chosen[2]["insertions_to_modal_consensus"], 1)
        self.assertEqual(chosen[2]["previous_deletions_to_long_consensus"], 2)
        corrected_parent = SimpleNamespace(
            consensus=tip,
            list_gaped_repeats=(
                aligned.replace("-", " "),
                outlier,
                aligned.replace("-", " "),
            ),
        )
        self.assertIsNone(select_alignment(source, p, corrected_parent))
        # Separate full-length copies may carry the same suffix, making them
        # uninformative about which end is justified.
        near, _ = array([outlier] * 3, [sequence(1, 25), sequence(2, 25)])
        self.assertIsNone(select(source + near, p, {p: 0.65}))
        allele = outlier[1:]
        near, _ = array([allele] * 3, [sequence(1, 25) + "G", sequence(2, 25) + "C"])
        self.assertIsNone(select_alignment(source + near, p, parent))
        source += near
        self.check_orientations(
            source, p, {p: 0.65}, expected, "uniform_edge_substitutions"
        )
        self.assertTrue(select(source, p, {p: 0.65})[2]["alternative_single_insertion"])
        source, p = array(
            [tip, tip[:12] + "C" + tip[12:], tip],
            [sequence(300, 25), sequence(400, 30)],
        )
        self.assertIsNone(select(source, p, {p: 0.65}))
        source, p = array([tip, "A" + tip, tip], [sequence(300, 25), sequence(400, 30)])
        self.assertIsNone(select(source, p, {p: 0.65}))

    def test_polymorphic_edge_needs_two_native_subpaths_and_a_variable_outer_column(
        self,
    ):
        a = sequence(13, 30)
        b = a[:14] + ("A" if a[14] != "A" else "T") + a[15:]
        alleles = ["C" + a] * 5 + ["G" + b] * 5
        spacers = [sequence(i, 26) + "GACT"[i % 4] for i in range(9)]
        source, p = array(alleles, spacers)
        primary = tuple((a + 1, b) for a, b in p)
        retained = {primary: 0.74, p[:5]: 0.81, p[5:]: 0.82}
        self.check_orientations(source, primary, retained, p, "retained_allele_flanks")
        self.assertIsNone(select(source, primary, {primary: 0.74, p[:5]: 0.81}))
        weak = {primary: 0.74, p[:5]: 0.81, p[5:]: 0.42}
        self.assertIsNone(select(source, primary, weak))
        # Additional conserved bases would leave the extent of the edge ambiguous.
        spacers = [sequence(i, 26) + ("A" if i < 4 else "T") for i in range(9)]
        source, _ = array(alleles, spacers)
        self.assertIsNone(select(source, primary, retained))

    def test_trimming_requires_one_terminal_polymorphic_base_and_observed_exact_core(
        self,
    ):
        motif = sequence(208, 28)
        mutant = list(motif)
        mutant[8] = "A" if mutant[8] != "A" else "C"
        copies = [motif] * 5 + ["".join(mutant)]
        spacers = [sequence(i, 28) for i in range(5)]
        source, primary = array([s + x for s, x in zip(copies, "AACAAG")], spacers)
        p = tuple((a, b - 1) for a, b in primary)
        witnesses = {primary: 0.90, p: 0.91, p[:5]: 0.95}
        self.check_orientations(
            source, primary, witnesses, p, "retained_polymorphic_edge_trim"
        )
        self.assertIsNone(select(source, primary, {primary: 0.90, p: 0.91}))
        self.assertIsNone(
            select(source, primary, {primary: 0.92, p: 0.91, p[:5]: 0.95})
        )
        source, _ = array([s + x for s, x in zip(copies, "AAAACG")], [spacers[0]] * 5)
        self.assertIsNone(select(source, primary, witnesses))

    def test_unanchored_slippage_uses_qualified_width_on_a_score_tie(self):
        motif = "GCTA" + sequence(176, 22)
        longer = motif[:2] + motif[2] + motif[2:]
        gapped = motif[:2] + "-" + motif[2:]
        flank = sequence(533, 18)
        head = len(flank)
        tail = sequence(701, 20)
        source = flank + motif + sequence(101, 30) + longer + tail + flank + motif
        b = head + len(motif) + 30
        c = b + len(longer) + len(tail) + len(flank)
        old = ((head, head + len(motif)), (b, b + len(longer)), (c, c + len(motif)))
        trimmed = (old[0], (old[1][0] + 1, old[1][1]), old[2])
        parent = SimpleNamespace(
            consensus=longer, list_gaped_repeats=(gapped, longer, gapped)
        )
        for reverse in (False, True):
            dna, path, expected, donor = source, old, trimmed, parent
            if reverse:
                dna, path, expected = (
                    reverse_complement(source),
                    reflect(source, old),
                    reflect(source, trimmed),
                )
                donor = SimpleNamespace(
                    consensus=reverse_complement(longer),
                    list_gaped_repeats=tuple(
                        reverse_complement(x)
                        for x in reversed(parent.list_gaped_repeats)
                    ),
                )
            proposal = select_alignment(dna, path, donor)
            self.assertEqual(proposal[1], "ambiguous_edge_width_competition")
            support = proposal[2]["uniform_width_preference"]
            self.assertEqual(
                len(support["indistinguishable_insertion_positions_0based"]), 2
            )
            for eligible, score, chosen in (
                (True, 0.61, expected),
                (False, 0.67, path),
                (True, 0.60, path),
            ):
                result = rank_alignment(
                    proposal,
                    uniform_eligible=eligible,
                    compact_score=0.61,
                    uniform_score=score,
                )
                self.assertEqual(result[0], chosen)
                other = path if chosen == expected else expected
                self.assertEqual(
                    result[2]["alternative_repeat_intervals"], list(map(list, other))
                )
        # Either independent edge evidence or a separate allele occurrence
        # removes the only grounds for applying this display convention.
        for extra in (motif, longer, longer[1:], reverse_complement(longer)):
            self.assertEqual(
                select_alignment(source + tail + extra, old, parent)[1],
                "modal_insertion_alignment",
            )
        x = head - 1
        changed = source[:x] + ("C" if source[x] == "A" else "A") + source[x + 1 :]
        self.assertEqual(
            select_alignment(changed, old, parent)[1], "modal_insertion_alignment"
        )


if __name__ == "__main__":
    unittest.main()
