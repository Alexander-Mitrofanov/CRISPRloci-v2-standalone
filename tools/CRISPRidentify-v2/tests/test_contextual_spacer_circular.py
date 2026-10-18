"""The duplicate-block exception does not authorize circular join promotion.

Synthetic source and topology exercise the real join and structural path.
Scoring/CDS callbacks are mocks, never extraction or biological evidence.
"""
import copy
import hashlib
import random
import unittest
from unittest.mock import Mock

from components.circular_array_join import CircularTopologyEvidence
from components.circular_array_join_integration import (
    candidate_from_circular_view,
    evaluate_accepted_circular_edge_joins,
)
from components.components_detection_refinement import CrisprCandidate
from components.localized_spacer_context import localized_block
from components.module_evaluation import build_advanced_fuzzy_search_filter
from tests.test_reconstruction_validation import PARAMETERS, REPEAT, candidate_source


def edge_fixture(*, duplicated):
    rng = random.Random(511173)
    spacers = ["".join(rng.choices("ACGT", k=32)) for _ in range(10)]
    head_source, head = candidate_source([REPEAT] * 3, spacers[:2])
    tail_spacers = spacers[2:]
    if duplicated:
        tail_spacers[2:6] = [tail_spacers[2]] * 4
    tail_source, tail = candidate_source([REPEAT] * 9, tail_spacers)
    prefix, middle, suffix = "ACGT" * 5, "G" * 400, "TGCA" * 5
    source = prefix + head_source + middle + tail_source + suffix

    def shifted(candidate, offset):
        return CrisprCandidate(
            list(candidate.list_repeats), list(candidate.list_repeats_gaped),
            list(candidate.list_spacers),
            [start + offset for start in candidate.list_repeat_starts],
        )

    head = shifted(head, len(prefix))
    tail = shifted(tail, len(prefix) + len(head_source) + len(middle))
    categories = [{}, {}, {}, {}, {}]
    categories[0][(head.list_repeat_starts[0], head.list_repeat_starts[-1]
                   + len(REPEAT))] = [[0.9, head, []]]
    # Model the source geometry/category consumed downstream, not a fake
    # contextual proof. The circular code receives an accepted Possible row.
    categories[2][(tail.list_repeat_starts[0], tail.list_repeat_starts[-1]
                   + len(REPEAT))] = [[0.82, tail, []]]
    topology = CircularTopologyEvidence(
        topology="circular", accession="SYNTHETIC_CONTEXT_CIRCULAR",
        source_sha256=hashlib.sha256(source.encode()).hexdigest(),
        evidence_method="synthetic_test_fixture", evidence_date="2026-10-05",
        evidence_sha256=hashlib.sha256(b"synthetic topology declaration").hexdigest(),
        evidence_complete=True, evidence_ambiguous=False,
    )
    return source, categories, topology, tail


class ContextualSpacerCircularTests(unittest.TestCase):
    def test_intact_internal_duplicate_run_cannot_inherit_context_exception(self):
        structural = build_advanced_fuzzy_search_filter(PARAMETERS)
        # Positive routing control: this edge layout and legal junction reach
        # the later mocked gates when the spacers are diverse.
        source, categories, topology, _tail = edge_fixture(duplicated=False)
        scorer, cds = Mock(return_value=0.99), Mock(return_value=())
        ordinary = evaluate_accepted_circular_edge_joins(
            source, categories, topology, PARAMETERS, scorer,
            source_cds_predictor=cds,
        )
        self.assertEqual(1, len(ordinary.joined))
        scorer.assert_called_once()
        cds.assert_called_once()

        source, categories, topology, tail = edge_fixture(duplicated=True)
        plan = localized_block(tail, structural)
        self.assertEqual([2, 6], plan["spacer_run"])
        self.assertEqual([0, 3], plan["left_repeat_indexes"])
        self.assertEqual([6, 9], plan["right_repeat_indexes"])
        before = copy.deepcopy([row[1].__dict__ for category in categories
                                for rows in category.values() for row in rows])
        scorer, cds = Mock(return_value=0.99), Mock(return_value=())
        result = evaluate_accepted_circular_edge_joins(
            source, categories, topology, PARAMETERS, scorer,
            source_cds_predictor=cds,
        )
        self.assertEqual(1, len(result.decisions))
        self.assertEqual((), result.joined)
        decision = result.decisions[0]
        self.assertEqual("rejected", decision.status)
        self.assertEqual(("fresh_linearized_structural_filter_rejected",),
                         decision.reason_codes)
        self.assertFalse(decision.structural_decision.passed)
        self.assertIsNone(decision.score_evidence)
        self.assertIsNone(decision.tandem_decision)
        scorer.assert_not_called()
        cds.assert_not_called()
        _content, joined = candidate_from_circular_view(
            decision.proposal.linearized_view)
        evidence = structural.evidence(joined).as_dict()
        self.assertEqual(["identical_spacer_burden"],
                         [p["predicate"] for p in evidence["predicate_results"]
                          if p["status"] == "failed"])
        self.assertEqual(before, [row[1].__dict__ for category in categories
                                  for rows in category.values() for row in rows])
        self.assertEqual(0.82, next(iter(categories[2].values()))[0][0])
        self.assertEqual(1, sum(map(len, categories[2].values())))


if __name__ == "__main__":
    unittest.main()
