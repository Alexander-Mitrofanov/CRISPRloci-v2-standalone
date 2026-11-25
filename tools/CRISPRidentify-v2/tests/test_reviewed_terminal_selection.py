import random
import unittest

from components.refinement.candidates import CrisprCandidate
from components.refinement.reviewed_terminal_selection import select, terminal_consensus_basis
from components.reporting.legacy_exports import _display_candidate
from components.sequence.utils import reverse_complement


def fixture(*, both=False, stronger_short_copy=False, ambiguous=False):
    rng = random.Random(80206)
    motif = "GTGGTCAATGTCGGATGCAGATTAGGAC"
    seq = list(motif)
    for i in (2, 5, 11, 15, 20, 24):
        seq[i] = next(c for c in "ACGT" if c != seq[i])
    hit = "".join(seq)
    if stronger_short_copy:
        # Five substitutions plus a one-base deletion at the inner side.
        # The equal-width placement would also have to absorb a spacer base.
        hit = hit[:13] + hit[14:]
    if ambiguous:
        hit = hit[:10] + "N" + hit[11:]
    dna = "C" * 80
    primary = []
    for _ in range(10):
        start = len(dna)
        dna += motif
        primary.append((start, len(dna)))
        dna += "".join(rng.choices("ACGT", k=32))
    wanted = len(dna), len(dna) + len(hit)
    dna += hit + "C" * 80
    if both:
        dna = dna[:20] + hit + dna[48:]
    return dna, tuple(primary), wanted


def reverse_path(path, length):
    return tuple((length - b, length - a) for a, b in reversed(path))


class ReviewedTerminalSelectionTests(unittest.TestCase):
    def test_only_the_divergent_end_breaks_a_tie_between_intact_alleles(self):
        first, second, terminal = "CGTAGCTA", "CGTTCCTA", "CATTCATA"
        core = [first, first, second, second, second, second, first, first]
        rows = [terminal, *core]
        classification = {
            "method": "source_supported_competing_geometry",
            "rule": "regular_core_divergent_terminal",
            "evidence": {"side": "left", "core_repeat_count": len(core)},
        }
        candidate = CrisprCandidate(rows, rows, ["AGC"] * len(core), list(range(0, len(rows) * 11, 11)))
        record = (0.78, candidate, [], {"classification_basis": classification})
        self.assertEqual(candidate.consensus, second)
        observed = terminal_consensus_basis(rows, classification)
        self.assertEqual(observed["reported_consensus"], first)
        self.assertEqual(observed["scored_full_alignment_consensus"], second)
        self.assertEqual([x["alignment_column_1based"] for x in observed["tied_columns"]], [4, 5])
        for reverse in (False, True):
            view = _display_candidate(record, reverse=reverse)
            expected = reverse_complement(first) if reverse else first
            self.assertEqual(view.consensus, expected)
            self.assertIn(expected, view.dot_repr())
            self.assertIn(expected, view.dot_repr_web_server())
        self.assertEqual(candidate.consensus, second)
        rc = [reverse_complement(r) for r in reversed(rows)]
        classification["evidence"]["side"] = "right"
        opposite = terminal_consensus_basis(rc, classification)
        self.assertEqual(opposite["reported_consensus"], reverse_complement(first))
        self.assertEqual(opposite["scored_full_alignment_consensus"], reverse_complement(second))
        classification["evidence"]["side"] = "left"
        self.assertIsNone(terminal_consensus_basis([terminal, *core[1:], second], classification))
        untied = [terminal, first, first, second, second, second, first, first, first]
        self.assertIsNone(terminal_consensus_basis(untied, classification))
        self.assertIsNone(terminal_consensus_basis(rows, {"method": "model_score_band"}))

    def test_regular_spacers_identify_a_single_divergent_copy_on_either_strand(self):
        source, path, wanted = fixture()
        for dna, core, interval, side in (
            (source, path, wanted, "right"),
            (
                source.translate(str.maketrans("ACGT", "TGCA"))[::-1],
                reverse_path(path, len(source)),
                reverse_path((wanted,), len(source))[0],
                "left",
            ),
        ):
            with self.subTest(side=side):
                result, name, evidence = select(dna, core, {core: 0.8})
                self.assertEqual(result, tuple(sorted((*core, interval))))
                self.assertEqual(name, "regular_core_divergent_terminal")
                self.assertEqual(evidence["side"], side)
                self.assertEqual(evidence["edit_distance"], 6)
                self.assertEqual(evidence["core_same_length_spacers"], 9)
                self.assertGreater(evidence["next_full_width_edits"], 6)

    def test_every_tied_variable_width_interval_is_retained(self):
        dna, core, (a, b) = fixture()
        replacement = next(c for c in "ACGT" if c != dna[b - 1])
        dna = dna[: b - 1] + replacement + dna[b:]
        _, _, evidence = select(dna, core, {core: 0.8})
        self.assertEqual(evidence["coordinate_system"], "zero_based_half_open")
        self.assertTrue(evidence["equal_edit_alternative_boundary"])
        self.assertEqual(
            evidence["equal_edit_alternative_intervals"], [(a, b - 1), (a, b + 1)]
        )

    def test_abstains_when_a_shorter_alignment_fits_better(self):
        dna, core, _ = fixture(stronger_short_copy=True)
        self.assertIsNone(select(dna, core, {core: 0.9}))

    def test_no_forced_end_when_both_sides_qualify(self):
        dna, core, _ = fixture(both=True)
        self.assertIsNone(select(dna, core, {core: 0.9}))

    def test_short_irregular_or_low_score_core_cannot_authorize_extension(self):
        dna, core, _ = fixture()
        irregular = core[:4] + ((core[4][0] + 7, core[4][1] + 7),) + core[5:]
        for path, score in (
            (core[-3:], 0.99),
            (irregular, 0.99),
            (core, 0.5),
            (core, float("nan")),
        ):
            self.assertIsNone(select(dna, path, {path: score}))

    def test_ambiguous_bases_and_motif_recurring_in_spacer_block_extension(self):
        dna, core, _ = fixture(ambiguous=True)
        self.assertIsNone(select(dna, core, {core: 0.9}))
        dna, core, (a, b) = fixture()
        spacer_start = core[0][1]
        repeat_in_spacer = dna[:spacer_start] + dna[a:b] + dna[spacer_start + b - a :]
        self.assertIsNone(select(repeat_in_spacer, core, {core: 0.9}))
