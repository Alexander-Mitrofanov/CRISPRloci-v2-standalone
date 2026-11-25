"""Retained terminal paths may preserve every base of a consistently shorter core."""
import unittest

from components.reconstruction.possible_reconstruction import PossibleReconstructionPolicy
from components.reconstruction.reconstruction_selection import relationship, select_reconstructions
from tests.test_reconstruction_selection import candidate, categories, valid


CORE = "ACGT" * 5
FULL = CORE + "TACG"
STARTS = [100, 160, 220, 280]


def copies(starts, seqs):
    value = candidate(starts, repeats=seqs)
    value.consensus = seqs[0]
    return value


def core():
    return copies(STARTS, [CORE] * 4)


class UniformExpansionTests(unittest.TestCase):
    def test_fixed_left_and_fixed_right_edges(self):
        left = copies([40] + STARTS, [FULL] * 5)
        right = copies([p - 4 for p in STARTS] + [336], ["GTAC" + CORE] * 5)
        for proposed, edge, offsets, mapping, neighbor in (
            (left, "start", [0, 4], [1, 2, 3, 4], 1),
            (right, "end", [-4, 0], [0, 1, 2, 3], 3),
        ):
            with self.subTest(edge=edge):
                r = relationship(core(), proposed)
                self.assertTrue(r["qualified"])
                self.assertEqual(r["correspondence"]["anchor_to_candidate_indexes"], mapping)
                proof = r["uniform_outward_extension"]
                self.assertEqual(proof["fixed_edge"], edge)
                self.assertEqual(proof["donor_minus_core_offsets"], offsets)
                self.assertEqual(proof["flanking_copy_checks"][0]["mapped_donor_index"], neighbor)
                self.assertEqual(proof["flanking_copy_checks"][0]["edit_distance"], 0)

    def test_rejects_nonuniform_shifts_inward_motion_and_internal_occurrence_changes(self):
        proposals = [
            copies(STARTS + [340], [FULL, FULL, FULL[:-1], FULL, FULL]),
            copies(STARTS + [340], [CORE[:-4]] * 5),
            copies([100, 160, 216, 276, 336], [FULL] * 5),
            copies([100, 160, 185, 220, 280, 340], [FULL] * 6),
            copies([100, 160, 280, 340, 400], [FULL] * 5),
            copies(STARTS + [340], [CORE + "TACGTA"] * 5),
        ]
        for proposed in proposals:
            with self.subTest(geometry=proposed.list_repeat_starts):
                self.assertFalse(relationship(core(), proposed)["qualified"])

    def test_added_copy_needs_local_support(self):
        proposed = copies(STARTS + [340], [FULL] * 4 + ["A" * 24])
        self.assertFalse(relationship(core(), proposed)["qualified"])

    def test_preserves_small_core_and_possible_tier_contracts(self):
        proposed = copies(STARTS + [340], [FULL] * 5)
        short_core = copies(STARTS[:3], [CORE] * 3)
        self.assertFalse(relationship(short_core, proposed)["qualified"])
        self.assertFalse(relationship(core(), proposed, PossibleReconstructionPolicy())["qualified"])

    def test_collision_remains_a_veto(self):
        anchor = [0.90, core(), []]
        donor = [0.88, copies(STARTS + [340], [FULL] * 5), []]
        collision = [0.70, copies([345, 410, 470], [CORE] * 3), []]
        result = select_reconstructions(categories(anchor, [donor], [collision]), valid)
        self.assertNotEqual(result.decisions[0].action, "selected")
        self.assertEqual(result.decisions[0].evidence["geometry_rejections"]["overlaps_other_accepted_core"], 1)

    def test_strict_mapping_proof_keeps_its_previous_shape(self):
        proposed = copies(STARTS + [340], [CORE] * 5)
        r = relationship(core(), proposed)
        self.assertTrue(r["qualified"])
        self.assertNotIn("uniform_outward_extension", r)
