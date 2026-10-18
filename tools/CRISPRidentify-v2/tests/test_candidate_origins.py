"""Repair-only ownership survives real normalization and rebuild pathways."""
import copy
import hashlib
import json
import unittest

from components.candidate_origins import (
    candidate_origin_evidence, candidate_origin_kind, inherit_candidate_origins,
    mark_repair_origin,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
    collapse_source_backed_phase_equivalent_candidates,
)
from components.module_detection_refinement import DetectionRefinement
from components.module_evaluated_arrays_enhancement import EvaluatedArraysEnhancement
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.repeat_boundary_repair import propose_terminal_repairs, source_geometry
from tests.test_repeat_boundary_repair import synthetic, fuzzy, PARAMETERS, Locus


class CandidateOriginTests(unittest.TestCase):
    def fixture(self, count=7):
        source, parent = synthetic(count=count)
        children, records, _ = propose_terminal_repairs(
            source, [parent], build_advanced_fuzzy_search_filter(PARAMETERS))
        self.assertEqual(len(children), 1)
        return source, parent, children[0], next(r for r in records if r["status"] == "generated")

    def test_factory_bindings_and_attestation_hash_are_independently_reproducible(self):
        source, parent, child, record = self.fixture()
        self.assertEqual(candidate_origin_kind(parent), "ordinary")
        self.assertEqual(candidate_origin_kind(child), "repair_only")
        evidence = candidate_origin_evidence(child, hashlib.sha256(source.encode()).hexdigest())
        self.assertIsNone(evidence["ordinary_birth_sha256"])
        self.assertEqual(evidence["repair_births"][0]["attestation_sha256"], hashlib.sha256(
            json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest())
        identity = [child.list_repeat_starts, child.list_repeats,
                    child.list_repeats_gaped, child.list_spacers]
        self.assertEqual(evidence["current_representation_sha256"], hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest())

    def test_source_rebuild_normalization_and_changed_geometry_keep_repair_ancestry(self):
        source, _parent, child, _record = self.fixture()
        starts, ends = source_geometry(child)
        rebuilt = build(source, child, starts, ends)
        shorter = build(source, child, starts[:-1], ends[:-1])
        normalized = DetectionRefinement._normalize_fuzzy_candidate(fuzzy(child))
        for derived in (rebuilt, shorter, normalized):
            self.assertEqual(candidate_origin_kind(derived), "repair_only")
            self.assertEqual(candidate_origin_evidence(derived)["repair_births"],
                             candidate_origin_evidence(child)["repair_births"])
        self.assertNotEqual(candidate_origin_evidence(shorter)["current_representation_sha256"],
                            candidate_origin_evidence(child)["current_representation_sha256"])

    def test_real_equality_dedup_preserves_independent_ordinary_witness(self):
        source, _parent, child, _record = self.fixture()
        natural = build(source, None, *source_geometry(child))
        for inputs in ([copy.deepcopy(child), copy.deepcopy(natural)],
                       [copy.deepcopy(natural), copy.deepcopy(child)]):
            refinement = DetectionRefinement.__new__(DetectionRefinement)
            refinement.accepted_only = False
            refinement.dict_fuzzy_crisprs_fully_refined = {Locus(0, len(source)): inputs}
            refinement._reformat_ac_crispr_candidates()
            combined = refinement.output()[(0, len(source))]
            self.assertEqual(len(combined), 1)
            evidence = candidate_origin_evidence(combined[0])
            self.assertEqual(evidence["kind"], "ordinary")
            self.assertIsNotNone(evidence["ordinary_birth_sha256"])
            self.assertEqual(len(evidence["repair_births"]), 1)

    def test_real_phase_canonicalization_does_not_launder_repair_birth(self):
        source, parent = synthetic(count=28)
        starts, ends = source_geometry(parent)
        bases = list(source)
        for index, end in enumerate(ends):
            bases[end - 1] = "ACGT"[index % 4]
        source = "".join(bases)
        parent = build(source, None, starts, ends)
        children, _, _ = propose_terminal_repairs(source, [parent],
            build_advanced_fuzzy_search_filter(PARAMETERS))
        collapsed, _, stats = collapse_source_backed_phase_equivalent_candidates(
            source, {(0, len(source)): [parent, children[0]]})
        self.assertEqual(stats["groups_collapsed"], 1)
        self.assertEqual(candidate_origin_kind(collapsed[(0, len(source))][0]), "repair_only")

    def test_enhancement_rebinds_provisional_ordinary_constructor_to_real_parent(self):
        source, _parent, child, _record = self.fixture()
        starts, ends = source_geometry(child)
        provisional = build(source, None, starts[:-1], ends[:-1])
        self.assertEqual(candidate_origin_kind(provisional), "ordinary")
        stage = EvaluatedArraysEnhancement.__new__(EvaluatedArraysEnhancement)
        stage.dna, stage.source_geometry_lineage = source, []
        key = (0, len(source))
        base = [.8, child, []]
        category = {key: [base]}
        stage._retain_derived_geometry(0, category, key, base, provisional, "test_enhancement")
        self.assertEqual(candidate_origin_kind(category[key][1][1]), "repair_only")

    def test_stale_parent_and_wrong_source_fail_closed_without_default_legacy_exceptions(self):
        source, _parent, child, record = self.fixture()
        wrong = copy.deepcopy(record)
        wrong["source_sha256"] = "f" * 64
        counterfeit = copy.deepcopy(child)
        mark_repair_origin(counterfeit, wrong)
        with self.assertRaisesRegex(ValueError, "current source"):
            candidate_origin_evidence(counterfeit, hashlib.sha256(source.encode()).hexdigest())
        mutated = copy.deepcopy(child)
        mutated.list_spacers[0] += "A"
        target = copy.deepcopy(child)
        inherit_candidate_origins(target, (mutated,))  # Legacy transformation records, not raises.
        with self.assertRaisesRegex(ValueError, "unresolved"):
            candidate_origin_kind(target)


if __name__ == "__main__":
    unittest.main()
