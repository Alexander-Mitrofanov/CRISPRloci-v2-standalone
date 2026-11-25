import random
import unittest

from components.refinement.reviewed_membership_selection import rank_extension, select
from components.sequence.utils import reverse_complement


def _spacers(n, width):
    rng = random.Random(839)
    return ["".join(rng.choices("ACGT", k=width)) for _ in range(n)]


def _array(repeats, spacers):
    source, path = "CACTATGC", []
    for i, repeat in enumerate(repeats):
        start = len(source)
        source += repeat
        path.append((start, len(source)))
        if i < len(spacers):
            source += spacers[i]
    return source + "CTAGCGATAC", tuple(path)


def _mutate(seq, *positions):
    result = list(seq)
    for position in positions:
        result[position] = "CGTA"["ACGT".index(result[position])]
    return "".join(result)


class ReviewedMembershipTests(unittest.TestCase):
    motif = "TGGACGTAGATCGTGCTAACTCTAGCCAT"

    def test_local_spacing_and_independent_sequence_place_internal(self):
        m = self.motif
        dna, path = _array([m, _mutate(m, 7, 13, 19), m, m], _spacers(3, 32))
        old = path[:1] + path[2:]
        full, rule, details = select(dna, old, {old: .65, path: .47})
        self.assertEqual((path, "periodic_internal_insertion"), (full, rule))
        self.assertEqual(details["split_spacer_interval"], [path[0][1], path[2][0]])
        self.assertEqual(details["spacer_lengths"], [32] * 3)
        self.assertIsNone(select(dna, old, {old: .65}))

    def test_real_long_spacer_and_repetitive_spacers(self):
        m = self.motif
        for extras, gaps in [
            ("ACGTGCAAGCT" * 2 + "GTCAGTA", _spacers(3, 32)),
            (m, ["GAT" * 10 + "CG"] * 3),
        ]:
            dna, path = _array([m, extras, m, m], gaps)
            old = path[:1] + path[2:]
            self.assertIsNone(select(dna, old, {old: .65, path: .85}))

    def test_clipped_tail_source_supports_two_full_substitution_alignments(self):
        m = self.motif
        reps = [m, m, _mutate(m, 3, 9, 28), _mutate(m, 3, 15, 21)]
        dna, path = _array(reps, _spacers(3, 34))
        old = path[:2] + ((path[2][0], path[2][1] - 1),)
        donor = old + ((path[3][0], path[3][1] - 3),)
        result = select(dna, old, {old: .61, donor: .48})
        self.assertEqual(result[:2], (path, "retained_substitution_width_extension"))

    def test_recurrent_internal_deletion_survives_new_terminal_membership(self):
        m = self.motif
        clipped = m[:8] + m[10:]
        reps = [m, m, clipped, clipped, _mutate(m, 8), _mutate(m, 14),
                _mutate(m, 11), _mutate(m, 8, 16), _mutate(m, 22)]
        dna, path = _array(reps, _spacers(8, 32))
        old, donor = path[:6], path[:7]
        result = select(dna, old, {old: .61, donor: .55})
        self.assertEqual(result[:2], (path, "retained_indel_allele_extension"))
        alt = result[2]["alternative_repeat_intervals"]
        self.assertEqual([b - a for a, b in alt], [len(m)] * 9)
        self.assertEqual(result[2]["alternative_metrics"]["different_copy_indexes_1based"], [3, 4])

        # The same ordered occurrences can carry a second qualified width
        # interpretation. Both metric sets must follow their actual intervals.
        picked, _, details = rank_extension(
            result, uniform_eligible=True, compact_score=.61, uniform_score=.62,
        )
        self.assertEqual(list(map(list, picked)), alt)
        self.assertEqual(details["alternative_repeat_intervals"], list(map(list, path)))
        self.assertEqual(details["spacer_lengths"], [c - b for (_, b), (c, _) in zip(picked, picked[1:])])
        self.assertEqual(details["repeat_edit_distances"], result[2]["alternative_metrics"]["repeat_edit_distances"])
        self.assertEqual(details["alternative_metrics"]["repeat_edit_distances"], result[2]["repeat_edit_distances"])
        self.assertEqual(details["alternative_metrics"]["different_copy_indexes_1based"], [3, 4])
        self.assertEqual(details["boundary_ranking"]["selected_interpretation"], "uniform_repeat_width")

        # The preference is a representation convention. Orientation-sensitive
        # score changes must not decide which genomic bases belong to a copy.
        def reverse(p):
            return tuple((len(dna) - b, len(dna) - a) for a, b in reversed(p))

        rc = select(reverse_complement(dna), reverse(old), {reverse(old): .61, reverse(donor): .55})
        self.assertEqual(rc[0], reverse(path))
        for score in (.59, .61, .62):
            rc_picked, _, proof = rank_extension(rc, uniform_eligible=True, compact_score=.61, uniform_score=score)
            self.assertEqual(rc_picked, reverse(picked))
            self.assertEqual(tuple(map(tuple, proof["alternative_repeat_intervals"])), reverse(path))
        blocked = rank_extension(result, uniform_eligible=False, compact_score=.61, uniform_score=.99)
        self.assertEqual(blocked, result)

    def test_tandem_duplication_must_cross_repeat_edge_and_fix_spacing(self):
        m = self.motif
        short = m[4:]
        spacers = _spacers(6, 32)
        spacers[0] = "A" + short[:8] + spacers[0][9:]
        spacers[0] = "A" + short + spacers[0]
        dna, path = _array([m, _mutate(m, 12), m, m, m, m, _mutate(m, 9)], spacers)
        expected = [(a + 4, b) for a, b in path]
        b = path[0][1]
        displaced = (b + 1, b + 1 + len(short))
        old = (expected[0], displaced, *expected[2:-1])
        eight = (path[0], (b, displaced[1]), *path[1:])
        inventory = {old: .98, path: .9, eight: .84}
        result = select(dna, old, inventory)
        self.assertEqual(result[:2], (path, "tandem_duplication_resegmentation"))
        self.assertEqual(result[2]["displaced_interval"], list(displaced))
        self.assertEqual(result[2]["duplication_adjusted_spacer_lengths"], [32] * 6)
        self.assertIsNone(select(dna, old, {old: .98, path: .9}))
        modified = dna[:displaced[1]+1] + "GACTGACT" + dna[displaced[1]+9:]
        self.assertIsNone(select(modified, old, inventory))
        # A long spacer unrelated to the length of the duplicated unit fails.
        other = list(path)
        other[-1] = (other[-1][0] - 5, other[-1][1] - 5)
        self.assertIsNone(select(dna, old, {old: .98, tuple(other): .9, eight: .84}))

    def test_tied_resegmentation_endpoints_abstain_in_either_inventory_order(self):
        # Both terminal shifts cost the same number of substitutions and
        # have the same spacer deviation. Only one has an exact eight-copy
        # witness; this must not be allowed to break the seven-copy tie.
        m = "A" * 10 + "C" * 9 + "G" * 10
        short = m[4:]
        gaps = _spacers(6, 32)
        gaps[0] = "T" + short + "T" + short[:8] + gaps[0][9:]
        gaps[-1] = gaps[-1][:-1] + "A"
        dna, path = _array([m, _mutate(m, 13), m, m, m, m, m], gaps)
        a, b = path[-1]
        dna = dna[:b] + "G" + dna[b+1:]
        start = path[0][1]
        lost = start + 1, start + 1 + len(short)
        old = ((path[0][0] + 4, start), lost, *((x + 4, y) for x, y in path[2:-1]))
        first = path[:-1] + ((a - 1, b - 1),)
        second = path[:-1] + ((a + 1, b + 1),)
        eight = (first[0], (start, lost[1]), *first[1:])
        self.assertIsNotNone(select(dna, old, {old: .96, first: .9, eight: .85}))
        candidates = {old: .96, first: .9, second: .89, eight: .85}
        self.assertIsNone(select(dna, old, candidates))
        self.assertIsNone(select(dna, old, dict(reversed(candidates.items()))))
