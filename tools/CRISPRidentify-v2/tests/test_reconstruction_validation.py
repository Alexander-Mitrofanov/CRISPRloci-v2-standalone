"""Proposed geometries must earn their own negative-screen evidence."""

import copy
import json
import random
import unittest

from components.architecture_rescue_coding_context import SourceCDS
from components.components_detection_refinement import (
    AdvancedFuzzySearchFilter,
    CrisprCandidate,
)
from components.reconstruction_validation import (
    needs_coding_context,
    validate_reconstruction_candidate,
)


PARAMETERS = {
    "param_min_avg_repeat_length": 21,
    "param_max_avg_repeat_length": 55,
    "param_min_avg_spacer_length": 18,
    "param_max_avg_spacer_length": 78,
    "param_min_repeats": 3,
    "param_max_identical_spacers": 4,
    "param_max_identical_cluster_spacers": 3,
}
REPEAT = "CTCGGCCCCGCGCTCGCGGGGGTTGCTCC"


def candidate_source(repeats, spacers):
    sequence = repeats[0]
    starts = [0]
    for spacer, repeat in zip(spacers, repeats[1:]):
        starts.append(len(sequence) + len(spacer))
        sequence += spacer + repeat
    return sequence, CrisprCandidate(
        list(repeats), list(repeats), list(spacers), starts,
    )


def diverse_candidate(count=21):
    rng = random.Random(19121)
    spacers = ["".join(rng.choices("ACGT", k=32)) for _ in range(count - 1)]
    repeats = [REPEAT] * count
    if count >= 21:
        # The two audited 19->21 inserted motifs: one and three substitutions.
        repeats[15] = "CTCGGCCCCGCGCTCGCGGGGGTTGCCCC"
        repeats[18] = "CTCGGCCCCGCGCTCGCAGGGATTGCTCG"
    return candidate_source(repeats, spacers)


def coding_tandem():
    repeats = ["GCC" * 8] * 11
    first = ("GCT", "ACT", "CCT", "TCT")
    last = ("GCT", "GCA", "GCC", "GCG")
    # Rotate the last codon separately so the ordinary end-dominance gate
    # does not itself reject this direct coding-tandem adversary.
    spacers = [first[index % 4] + "GCC" * 6 + last[(index // 4 + index) % 4]
               for index in range(10)]
    return candidate_source(repeats, spacers)


class ReconstructionCandidateValidationTests(unittest.TestCase):
    def setUp(self):
        self.structural = AdvancedFuzzySearchFilter.from_parameters(PARAMETERS)

    def validate(self, source, candidate, *, metadata=None, **kwargs):
        original = [0.76, candidate, [[3.0]]]
        if metadata is not None:
            original.append(metadata)
        # A different row models pipeline rescoring, without copying the
        # incumbent's confidence or clearance into the proposal.
        fresh = [0.71, candidate, [[float(len(candidate.list_repeats))]]]
        return validate_reconstruction_candidate(
            source, original, fresh, self.structural, **kwargs,
        )

    def test_long_source_backed_internal_repeat_parse_is_not_unknown_tandem(self):
        source, candidate = diverse_candidate()
        before = copy.deepcopy(candidate.__dict__)
        result = self.validate(source, candidate)
        self.assertTrue(result.qualified)
        self.assertEqual("inapplicable", result.evidence["likely_tandem"]["status"])
        self.assertEqual("spacer_count_above_policy_maximum",
                         result.evidence["likely_tandem"]["reason"])
        self.assertEqual("inapplicable", result.evidence["coding_tandem"]["status"])
        self.assertFalse(result.evidence["biological_correctness_certified"])
        self.assertTrue(result.evidence["policy_checks_passed"])
        self.assertEqual(before, candidate.__dict__)
        json.dumps(result.evidence, allow_nan=False)

    def test_source_mismatch_cannot_be_rescued_by_fresh_score(self):
        source, candidate = diverse_candidate()
        result = self.validate("A" + source[1:], candidate)
        self.assertFalse(result.source_valid)
        self.assertFalse(result.qualified)

    def test_alignment_tampering_fails_even_with_valid_source_intervals(self):
        source, candidate = diverse_candidate()
        candidate.list_gaped_repeats = list(candidate.list_gaped_repeats)
        candidate.list_gaped_repeats[0] = "A" + candidate.list_gaped_repeats[0][1:]
        result = self.validate(source, candidate)
        self.assertTrue(result.source_valid)
        self.assertFalse(result.structural_valid)
        self.assertEqual("failed", result.evidence["alignment_reconstruction"]["status"])

    def test_classification_ineligible_source_path_stays_ineligible(self):
        source, candidate = diverse_candidate()
        candidate.source_path_classification_eligible = False
        candidate.source_path_classification_reason = "terminal_proposal_unattested"
        result = self.validate(source, candidate)
        self.assertFalse(result.qualified)
        self.assertEqual("failed", result.evidence["candidate_eligibility"]["status"])

    def test_canonical_structural_filter_rejects_spacer_duplication_burden(self):
        source, candidate = candidate_source([REPEAT] * 7, ["ACGT" * 8] * 6)
        result = self.validate(source, candidate)
        self.assertFalse(result.structural_valid)
        predicates = {p["predicate"]: p["status"] for p in
                      result.evidence["structural_filter"]["predicate_results"]}
        self.assertEqual("failed", predicates["identical_spacer_burden"])

    def test_small_tandem_screen_runs_on_proposed_spacers(self):
        source, candidate = candidate_source([REPEAT] * 4, ["ACGT" * 8] * 3)
        result = self.validate(source, candidate)
        self.assertTrue(result.structural_valid)
        self.assertFalse(result.negatives_passed)
        self.assertEqual("failed", result.evidence["likely_tandem"]["status"])

    def test_small_diverse_spacers_have_actual_evaluated_clearance(self):
        source, candidate = diverse_candidate(4)
        result = self.validate(source, candidate)
        self.assertTrue(result.qualified)
        self.assertEqual("passed", result.evidence["likely_tandem"]["status"])
        self.assertEqual(3, result.evidence["likely_tandem"]["decision"]["pairwise_comparison_count"])

    def test_noncanonical_spacer_cannot_be_counted_as_clearance(self):
        source, candidate = diverse_candidate(4)
        spacers = list(candidate.list_spacers)
        spacers[0] = spacers[0][:10] + "N" + spacers[0][11:]
        source, candidate = candidate_source(candidate.list_repeats, spacers)
        result = self.validate(source, candidate)
        self.assertTrue(result.source_valid)
        self.assertFalse(result.negatives_passed)
        self.assertEqual("not_run", result.evidence["canonical_sequence"]["status"])
        self.assertEqual("not_run", result.evidence["likely_tandem"]["status"])

    def test_known_tandem_metadata_cannot_be_lost_during_rescoring(self):
        source, candidate = diverse_candidate()
        for metadata in (
                {"tags": ["tandem"]},
                {"classification_basis": {"method": "low_spacer_diversity_tandem_demotion"}},
                {"classification_basis": {"accepted_family_tandem_review": {"supported": True}}}):
            with self.subTest(metadata=metadata):
                result = self.validate(source, candidate, metadata=metadata)
                self.assertFalse(result.negatives_passed)
                self.assertEqual("failed", result.evidence["known_negative_annotations"]["status"])

    def test_coding_context_not_run_is_not_an_empty_assessed_result(self):
        source, candidate = coding_tandem()
        self.assertTrue(needs_coding_context(candidate))
        missing = self.validate(source, candidate)
        self.assertTrue(missing.structural_valid)
        self.assertFalse(missing.negatives_passed)
        self.assertEqual("not_run", missing.evidence["coding_tandem"]["status"])
        empty = self.validate(source, candidate, source_cds=())
        self.assertTrue(empty.qualified)
        self.assertEqual("passed", empty.evidence["coding_tandem"]["status"])
        self.assertEqual("eligible_source_cds_not_found", empty.evidence["coding_tandem"]["reason"])
        self.assertTrue(empty.evidence["coding_tandem"]["source_cds_assessed"])

    def test_disabled_applicable_coding_screen_cannot_authorize_promotion(self):
        source, candidate = coding_tandem()
        result = self.validate(source, candidate, source_cds=(), coding_screen_required=False)
        self.assertFalse(result.qualified)
        self.assertEqual("not_run", result.evidence["coding_tandem"]["status"])

    def test_long_coding_tandem_fails_even_when_small_tandem_rule_inapplicable(self):
        source, candidate = coding_tandem()
        result = self.validate(source, candidate, source_cds=(SourceCDS(0, len(source), False, 99.0),))
        self.assertTrue(result.structural_valid)
        self.assertFalse(result.negatives_passed)
        self.assertEqual("inapplicable", result.evidence["likely_tandem"]["status"])
        coding = result.evidence["coding_tandem"]
        self.assertEqual("failed", coding["status"])
        self.assertTrue(coding["decision"]["supported"])
        self.assertEqual(1, len(coding["considered_source_cds"]))

    def test_coding_evidence_omits_irrelevant_source_wide_cds(self):
        source, candidate = coding_tandem()
        source += "A" * 200
        result = self.validate(source, candidate, source_cds=(
            SourceCDS(0, 504, False, 99.0),
            SourceCDS(len(source) - 150, len(source) - 20, False, 95.0),
        ))
        self.assertEqual(1, len(result.evidence["coding_tandem"]["considered_source_cds"]))

    def test_original_and_fresh_rows_must_reference_the_same_candidate(self):
        source, candidate = diverse_candidate()
        with self.assertRaisesRegex(ValueError, "same candidate"):
            validate_reconstruction_candidate(
                source, [0.7, candidate, [[1]]],
                [0.7, copy.deepcopy(candidate), [[1]]], self.structural,
            )

    def test_malformed_recorded_negative_evidence_is_not_clearance(self):
        source, candidate = diverse_candidate()
        with self.assertRaisesRegex(ValueError, "support must be boolean"):
            self.validate(source, candidate, metadata={
                "classification_basis": {
                    "accepted_family_tandem_review": {"supported": "true"},
                },
            })

    def test_missing_or_nonfinite_fresh_score_features_rejected(self):
        source, candidate = diverse_candidate()
        for row in ([float("nan"), candidate, [[1]]], [0.7, candidate, None]):
            with self.subTest(row=row[0]), self.assertRaises(ValueError):
                validate_reconstruction_candidate(
                    source, [0.7, candidate, [[1]]], row, self.structural,
                )


if __name__ == "__main__":
    unittest.main()
