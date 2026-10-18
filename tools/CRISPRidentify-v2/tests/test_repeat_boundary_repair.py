import copy
import hashlib
import json
from collections import namedtuple
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch

from components.components_detection_refinement import CrisprCandidate
from components.components_evaluated_arrays_enhancement import (
    _minimal_edit_changes,
    _source_only_alignment_seed,
    build_source_backed_candidate_from_boundaries as build,
    collapse_source_backed_phase_equivalent_candidates,
    validate_candidate_source_reconstruction,
)
from components.fuzzy_alignment_merge import merge_fuzzy_alignments
from components.fast_run_rescue import candidate_score_fingerprint
from components.module_detection_refinement import DetectionRefinement
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.pipeline import CandidateScoreCache, Pipeline
from components.report_model import ReportModelBuilder
from components.repeat_boundary_repair import (
    MAX_PARENTS_PER_LOCUS, propose_terminal_repairs, source_geometry,
)


PARAMETERS = {
    "param_min_repeats": 3,
    "param_min_avg_repeat_length": 21,
    "param_max_avg_repeat_length": 55,
    "param_min_avg_spacer_length": 18,
    "param_max_avg_spacer_length": 78,
    "param_max_identical_spacers": 4,
    "param_max_identical_cluster_spacers": 3,
}
CORE = "GCCGTGGTTTCCCTACCGATTTGCCAGTGGTAGGTT"
Locus = namedtuple("Locus", "start end")


def synthetic(left=False, right=True, terminal_gap=False, core=CORE,
              count=7, identical_spacers=False, internal_indel=False):
    rng = random.Random(617)
    noise = "GGGCTAA"
    dna = "ACGT" * 10
    starts, ends = [], []
    for index in range(count):
        starts.append(len(dna))
        if left:
            dna += noise[index % 7]
        repeat = core
        if internal_indel and index == 0:
            repeat = core[:15] + core[16:]
        dna += repeat
        if right and not (terminal_gap and index == 6):
            dna += noise[index % 7]
        ends.append(len(dna))
        dna += ("AGCT" * 8 if identical_spacers else
                "".join(rng.choice("ACGT") for _ in range(32)))
    return dna, build(dna, None, starts, ends)


def fuzzy(candidate):
    result = copy.deepcopy(candidate)
    result.list_absolute_start = result.list_repeat_starts
    result.number_errors = result.total_mismatches
    return result


class RepeatBoundaryRepairTests(unittest.TestCase):
    def setUp(self):
        self.structural = build_advanced_fuzzy_search_filter(PARAMETERS)

    def repair(self, dna, candidates):
        return propose_terminal_repairs(dna, candidates, self.structural)

    def assert_reconstructs(self, dna, candidate):
        validate_candidate_source_reconstruction(
            dna, candidate.list_repeats, candidate.list_repeat_starts,
            candidate.list_spacers,
        )
        self.assertTrue(self.structural.evidence(candidate).passes)

    def test_saved_source_geometry_fixtures(self):
        fixtures = json.loads((Path(__file__).parent / "data" /
                               "repeat_boundary_repair.json").read_text())
        self.assertEqual(3, len(fixtures))
        for fixture in fixtures:
            with self.subTest(accession=fixture["accession"]):
                dna = fixture["sequence"]
                self.assertEqual(fixture["sequence_sha256"],
                                 hashlib.sha256(dna.encode()).hexdigest())
                parent = build(dna, None, fixture["starts"], fixture["ends"])
                children, records, stats = self.repair(dna, [parent])
                self.assertEqual(1, len(children))
                child = children[0]
                self.assertEqual(len(parent.list_repeats), len(child.list_repeats))
                self.assertEqual(fixture["expected_repeat_length"], len(child.consensus))
                self.assertEqual(parent.list_repeat_starts, child.list_repeat_starts)
                self.assertEqual(tuple(end - delta[1] for end, delta in zip(
                    fixture["ends"], records[0]["consumed_bases"])),
                                 source_geometry(child)[1])
                self.assertEqual(15, stats["hypotheses_considered"])
                self.assertEqual(1, records[0]["right_trim_columns"])
                self.assert_reconstructs(dna, child)

    def test_rebuild_memoization_preserves_exact_heterogeneous_and_indel_state(self):
        for kwargs in ({"right": False}, {}, {"terminal_gap": True},
                       {"internal_indel": True}, {"left": True, "right": True}):
            dna, rebuilt = synthetic(**kwargs)
            repeats = rebuilt.list_repeats
            seed = _source_only_alignment_seed(repeats)
            # Independent per-occurrence replay is the pre-optimization path.
            _, gapped = merge_fuzzy_alignments(seed, repeats, [
                _minimal_edit_changes(seed, repeat) for repeat in repeats
            ])
            uncached = CrisprCandidate(repeats, gapped, rebuilt.list_spacers,
                                      rebuilt.list_repeat_starts)
            self.assertEqual(uncached, rebuilt)
            self.assertEqual(vars(uncached), vars(rebuilt))

    def test_left_right_both_and_internal_indel_keep_all_source_bases(self):
        for left, right, indel in ((True, False, False), (False, True, False),
                                    (True, True, False), (False, True, True)):
            with self.subTest(left=left, right=right, indel=indel):
                dna, parent = synthetic(left, right, internal_indel=indel)
                original = copy.deepcopy(parent)
                children, records, _ = self.repair(dna, [parent])
                self.assertEqual(1, len(children))
                self.assertEqual(original, parent)
                child = children[0]
                self.assertEqual(7, len(child.list_repeats))
                self.assertEqual(CORE, child.consensus)
                self.assertEqual(int(left), records[0]["left_trim_columns"])
                self.assertEqual(int(right), records[0]["right_trim_columns"])
                for index, spacer in enumerate(child.list_spacers):
                    expected = ((parent.list_repeats[index][-1] if right else "") +
                                parent.list_spacers[index] +
                                (parent.list_repeats[index + 1][0] if left else ""))
                    self.assertEqual(expected, spacer)
                self.assert_reconstructs(dna, child)

    def test_terminal_gap_consumes_zero_nucleotides(self):
        dna, parent = synthetic(terminal_gap=True)
        children, records, _ = self.repair(dna, [parent])
        self.assertEqual(1, len(children))
        self.assertEqual(0, records[0]["consumed_bases"][-1][1])
        self.assertEqual(source_geometry(parent)[1][-1],
                         source_geometry(children[0])[1][-1])
        self.assertEqual([CORE] * 7, children[0].list_repeats)
        self.assert_reconstructs(dna, children[0])

    def test_legacy_alignment_and_input_order_cannot_change_geometry(self):
        dna, parent = synthetic()
        # The same source nucleotides can be distributed over sparse tail
        # columns. The old raw alignment gate ignored those columns and saw
        # this representation as terminal-passing, unlike the canonical one.
        gapped = [repeat[:-1] + "".join(
            repeat[-1] if repeat[-1] == base else " " for base in "ACGT"
        ) for repeat in parent.list_repeats]
        alternative = CrisprCandidate(
            list(parent.list_repeats), gapped, list(parent.list_spacers),
            list(parent.list_repeat_starts),
        )
        self.assertFalse(self.structural._filter_by_column(parent))
        self.assertTrue(self.structural._filter_by_column(alternative))
        expected_children, expected_records, _ = self.repair(dna, [parent])
        self.assertEqual(1, len(expected_children))
        for members in ([parent], [alternative], [parent, alternative],
                        [alternative, parent]):
            children, records, stats = self.repair(dna, members)
            self.assertEqual(expected_children, children)
            self.assertEqual(expected_records, records)
            self.assertEqual(1, stats["parents_attempted"])
            self.assertEqual(len(members) - 1, stats["duplicate_geometries"])
            locus = Locus(source_geometry(parent)[0][0], source_geometry(parent)[1][-1])
            refined = DetectionRefinement(
                {locus: [fuzzy(candidate) for candidate in members]},
                PARAMETERS, False, dna,
            )
            expected_geometry = source_geometry(expected_children[0])
            self.assertEqual(1, sum(source_geometry(candidate) == expected_geometry
                                    for candidate in refined.output()[tuple(locus)]))

    def test_passing_short_core_does_not_prevent_long_parent_repair(self):
        dna, parent = synthetic()
        starts, ends = source_geometry(parent)
        short = build(dna, None, starts[:3], ends[:3])
        locus = Locus(starts[0], ends[-1])
        refinement = DetectionRefinement(
            {locus: [fuzzy(short), fuzzy(parent)]}, PARAMETERS, False, dna,
        )
        candidates = refinement.output()[tuple(locus)]
        self.assertTrue(any(len(c.list_repeats) == 3 for c in candidates))
        repaired = [c for c in candidates if len(c.list_repeats) == 7]
        self.assertEqual(1, len(repaired))
        self.assertEqual(CORE, repaired[0].consensus)
        self.assertNotEqual(candidate_score_fingerprint(parent),
                            candidate_score_fingerprint(repaired[0]))
        old_row = [0.91, parent, []]
        cache = CandidateScoreCache([{tuple(locus): [old_row]}, {}, {}, {}, {}])
        old_row[1] = repaired[0]
        scored = []

        def scorer(by_locus):
            scored.extend(by_locus[tuple(locus)])
            return {tuple(locus): [[0.73, repaired[0], []]]}

        fresh = cache.merge_after_enhancement(
            [{tuple(locus): [old_row]}, {}, {}, {}, {}], scorer,
        )
        self.assertEqual([repaired[0]], scored)
        self.assertEqual(0.73, fresh[tuple(locus)][0][0])
        self.assertEqual(1, cache.recomputed_count)

    def test_no_change_for_valid_mixed_duplicate_short_or_conserved_edge(self):
        cases = [synthetic(right=False), synthetic(core="ACGT" * 5),
                 synthetic(identical_spacers=True), synthetic(count=4)]
        dna, parent = synthetic()
        # Conserved outside column with variation one position inside cannot
        # be discarded merely to hide the true interior divergence.
        starts, ends = source_geometry(parent)
        chars = list(dna)
        for end in ends:
            chars[end] = "C"
        dna = "".join(chars)
        cases.append((dna, build(dna, None, starts, [end + 1 for end in ends])))
        for dna, parent in cases:
            with self.subTest(repeats=parent.list_repeats):
                children, _, _ = self.repair(dna, [parent])
                self.assertEqual([], children)

    def test_restrictions_missing_and_malformed_source_fail_closed(self):
        dna, parent = synthetic()
        for marker in (False, None, 1, "True"):
            restricted = copy.deepcopy(parent)
            restricted.source_path_classification_eligible = marker
            self.assertEqual([], self.repair(dna, [restricted])[0])
            self.assertEqual([], self.repair(dna, [parent, restricted])[0])
            normalized = DetectionRefinement._normalize_fuzzy_candidate(fuzzy(restricted))
            self.assertIs(marker, normalized.source_path_classification_eligible)
        self.assertEqual([], self.repair(None, [parent])[0])
        with self.assertRaisesRegex(ValueError, "source DNA"):
            self.repair(123, [parent])
        with self.assertRaisesRegex(ValueError, "source sequence|reconstruct"):
            self.repair(dna[:80], [parent])
        malformed = copy.deepcopy(parent)
        malformed.list_spacers[0] = "N" * len(malformed.list_spacers[0])
        for members in ([parent, malformed], [malformed, parent]):
            with self.assertRaisesRegex(ValueError, "spacer does not reconstruct"):
                self.repair(dna, members)

    def test_duplicate_fuzzy_geometry_keeps_restriction_before_repair(self):
        dna, parent = synthetic()
        restricted = fuzzy(parent)
        restricted.source_path_classification_eligible = False
        locus = Locus(*[source_geometry(parent)[i][j] for i, j in ((0, 0), (1, -1))])
        for members in ([fuzzy(parent), restricted], [restricted, fuzzy(parent)]):
            refinement = DetectionRefinement({locus: members}, PARAMETERS, False, dna)
            self.assertEqual(0, refinement.boundary_repair_stats["children_generated"])
            candidate = refinement.output()[tuple(locus)][0]
            self.assertFalse(candidate.source_path_classification_eligible)

    def test_insufficient_canonical_boundary_observations_abstain(self):
        dna, parent = synthetic()
        starts, ends = source_geometry(parent)
        bases = list(dna)
        for end in ends[:2]:
            bases[end - 2] = "N"
        dna = "".join(bases)
        parent = build(dna, None, starts, ends)
        children, records, _ = self.repair(dna, [parent])
        self.assertEqual([], children)
        self.assertGreater(records[0]["variant_rejections"]["weak_source_boundary"], 0)

    def test_ambiguous_minimum_never_uses_iteration_order_or_scores(self):
        # Exercise fail-closed decision logic if two distinct geometric trims
        # survive the independent filters. Real terminal-only edge evidence
        # usually makes this impossible, but future predicate changes must not
        # accidentally turn ambiguity into an arbitrary winner.
        dna, parent = synthetic(left=True)
        real_evidence = self.structural.evidence
        passed = real_evidence(build(dna, None,
            [s + 1 for s in source_geometry(parent)[0]],
            [e - 1 for e in source_geometry(parent)[1]]))
        calls = [0]

        def evidence(candidate):
            calls[0] += 1
            return real_evidence(candidate) if calls[0] == 1 else passed

        with patch("components.repeat_boundary_repair._weak_removed_columns", return_value=True), \
                patch("components.repeat_boundary_repair._source_phase_boundary_evidence", return_value={"passes": True}), \
                patch.object(self.structural, "evidence", side_effect=evidence):
            children, records, _ = self.repair(dna, [parent])
        self.assertEqual([], children)
        self.assertEqual("ambiguous_minimum_trim", records[0]["reason"])

    def test_budget_deduplication_and_idempotence(self):
        dna, base = synthetic(count=18)
        starts, ends = source_geometry(base)
        parents = [build(dna, None, starts[:count], ends[:count])
                   for count in range(7, 19)]
        first, records, stats = self.repair(dna, parents + [parents[0]])
        second, reversed_records, _ = self.repair(dna, list(reversed(parents)))
        self.assertEqual(first, second)
        self.assertEqual(MAX_PARENTS_PER_LOCUS, stats["parents_attempted"])
        self.assertEqual(15 * MAX_PARENTS_PER_LOCUS, stats["hypotheses_considered"])
        self.assertGreater(stats["parent_budget_skipped"], 0)
        self.assertEqual(stats["parent_budget_skipped"],
                         stats["canonical_eligibility_unexamined"])
        self.assertEqual(MAX_PARENTS_PER_LOCUS, stats["qualifying_parents"])
        self.assertEqual([r["parent_geometry"] for r in records],
                         [r["parent_geometry"] for r in reversed_records])
        self.assertEqual([], self.repair(dna, parents + first)[0])

    def test_phase_collapse_preserves_repair_and_cannot_launder_restriction(self):
        dna, parent = synthetic(count=28)
        starts, ends = source_geometry(parent)
        bases = list(dna)
        for index, end in enumerate(ends):
            bases[end - 1] = "ACGT"[index % 4]
        dna = "".join(bases)
        parent = build(dna, None, starts, ends)
        child = self.repair(dna, [parent])[0][0]
        collapsed, _, stats = collapse_source_backed_phase_equivalent_candidates(
            dna, {(0, len(dna)): [parent, child]},
        )
        self.assertEqual(1, stats["groups_collapsed"])
        canonical = collapsed[(0, len(dna))][0]
        self.assertEqual(child.boundary_repair_records, canonical.boundary_repair_records)
        restricted = copy.deepcopy(parent)
        restricted.source_path_classification_eligible = False
        collapsed, _, stats = collapse_source_backed_phase_equivalent_candidates(
            dna, {(0, len(dna)): [restricted, child]},
        )
        self.assertEqual(0, stats["groups_collapsed"])
        self.assertIs(restricted, collapsed[(0, len(dna))][0])
        self.assertFalse(restricted.source_path_classification_eligible)

    def test_audit_survives_rebuilt_objects_and_report_serialization(self):
        dna, parent = synthetic()
        child, records, stats = self.repair(dna, [parent])
        rebuilt = build(dna, None, *source_geometry(child[0]))
        pipeline = Pipeline.__new__(Pipeline)
        pipeline.repeat_boundary_repair_records = records
        pipeline.repeat_boundary_repair_stats = stats
        pipeline.categories = [{(0, len(dna)): [[0.9, rebuilt, []]]}, {}, {}, {}, {}]
        audit = pipeline._repeat_boundary_repair_audit()
        self.assertEqual("Bona-fide", audit["records"][0]["final_geometry_outcomes"][0]["category"])
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "source.fasta"
            fasta.write_text(">source\n" + dna + "\n")
            model = ReportModelBuilder(
                file_path=fasta, parameters={"param_repeat_boundary_repair": audit},
                flags={}, categories=[{}, {}, {}, {}, {}], non_array_data={},
                list_feature_names=[],
            ).build()
            restored = json.loads(json.dumps(model))
        self.assertEqual("generated", restored["run"]["parameters"]["param_repeat_boundary_repair"]["records"][0]["status"])


if __name__ == "__main__":
    unittest.main()
