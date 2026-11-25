"""Source-only review alternatives preserve canonical and retained paths."""

import copy
import unittest

from components.evidence.review_hypotheses import compare_paths


def build(repeats, spacers, flank="TTTTTT"):
    text = flank
    path = []
    for i, rep in enumerate(repeats):
        start = len(text)
        text += rep
        path.append((start, len(text)))
        if i < len(spacers):
            text += spacers[i]
    return text + flank, tuple(path)


def selected(comparison, path):
    return next(row for row in comparison["hypotheses"] if row["repeat_intervals"] == path)


class ReviewHypothesesTests(unittest.TestCase):
    def test_internal_repeat_requires_two_supported_spacers_and_keeps_original_score(self):
        source, extended = build(["ACCGT"] * 4, ["AACCCC", "GGAATT", "CCCCAA"])
        primary = extended[:1] + extended[2:]
        scored = [(extended, "observed", "Low score", 0.42)]
        unchanged = copy.deepcopy(scored)
        result = compare_paths(source, primary, scored)
        self.assertEqual(result["primary_repeat_intervals"], primary)
        evidence = selected(result, extended)["evidence"]["retained_internal_repeat"]
        self.assertEqual(evidence["category"], "Low score")
        self.assertEqual(evidence["observed_score"], 0.42)
        self.assertEqual(scored, unchanged)
        source, extended = build(["ACCGT"] * 4, ["A", "GGAATT" * 3, "CCCCAA"])
        self.assertIsNone(compare_paths(source, extended[:1] + extended[2:],
                                        [(extended, "wrong-gaps", "Possible", 0.8)]))

    def test_polymorphic_edge_with_recurrent_internal_alleles(self):
        repeats = ["CACCGTA", "GACCGTT"] * 4
        source, extended = build(repeats, ["CCCATA"] * 7)
        primary = tuple((a + 1, b) for a, b in extended)
        result = compare_paths(source, primary, [])
        record = selected(result, extended)
        evidence = record["evidence"]["allele_associated_flank"]
        self.assertEqual(evidence["internal_repeat_base_counts"],
                         {"ACCGTA": {"C": 4}, "ACCGTT": {"G": 4}})
        self.assertNotIn("observed_score", evidence)

    def test_conserved_flank_stops_on_local_loss_of_support(self):
        source, longer = build(["GGACGT", "GGACGT", "GTACGT", "GGACGT"],
                               ["CCCCA", "CCCCT", "CCCCC"], flank="TTTTG")
        primary = tuple((a + 2, b) for a, b in longer)
        evidence = selected(compare_paths(source, primary, []), longer)["evidence"]
        self.assertEqual(evidence["conserved_flank"]["columns"][-1],
                         {"A": 1, "C": 1, "G": 1, "T": 1})

    def test_indel_alignment_and_bounds_including_contig_start(self):
        source, primary = build(["ACAC", "ACAAC", "ACAC"], ["TCGA", "ATCG"], flank="")
        result = compare_paths(source, primary, [])
        a, b = primary[1]
        candidate = primary[:1] + ((a + 1, b),) + primary[2:]
        evidence = selected(result, candidate)["evidence"]["local_boundary"]
        aligned = evidence["primary_alignment_to_modal"]
        self.assertEqual([s.replace(" ", "").replace("-", "") for s in aligned],
                         ["ACAC", "ACAAC"])
        self.assertTrue(all(0 <= a < b <= len(source)
                            for h in result["hypotheses"] for a, b in h["repeat_intervals"]))
        self.assertIsNone(compare_paths(source, ((-1, 2), (5, 7), (9, 11)), []))

    def test_whole_path_segmentation_explicitly_displaces_a_repeat(self):
        primary = ((2, 6), (7, 11), (60, 64), (74, 78))
        other = ((1, 6), (35, 40), (59, 64), (73, 78), (87, 92))
        source = "AGCT" * 25
        observed = [(other, "one", "Alternative", 0.83), (other, "two", "Alternative", 0.90)]
        result = compare_paths(source, primary, observed)
        self.assertEqual(selected(result, other)["evidence"]["retained_resegmentation"]
                         ["observation_id"], "two")
        self.assertEqual(result, compare_paths(source, primary, reversed(observed)))
