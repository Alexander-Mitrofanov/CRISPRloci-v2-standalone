import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest

from components.components_detection import (
    FuzzySearch,
    FuzzySourceOccurrence,
    FuzzySourcePathProvenance,
    FUZZY_SOURCE_PATH_SCHEMA,
)
from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
)
from components.internal_repeat_completion import (
    InternalSearchContext,
    InternalSearchPolicy,
    analyze_gap_graph,
    analyze_internal_geometry,
    propose_internal_repairs,
    replay_internal_attestation,
    BROKEN_REASON,
)
from components.internal_repeat_transaction import (
    complete_internal_repeats,
    replay_internal_completion_basis,
    METHOD,
)
from components.module_detection_refinement import DetectionRefinement
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.repeat_boundary_repair import propose_terminal_repairs, source_geometry
from components.report_model import ReportModelBuilder
from tests.test_repeat_boundary_repair import PARAMETERS, fuzzy, Locus


MOTIF = "TGTTTCAGACATGCCCGGTTTAGACGGGATTAAGAC"


def lineage(candidate, broken=None):
    geometry = tuple(
        FuzzySourceOccurrence(
            FuzzySearch._source_occurrence_id(start, end), start, end, repeat
        )
        for (start, end), repeat in zip(
            zip(*source_geometry(candidate)), candidate.list_repeats
        )
    )
    spacers = tuple(candidate.list_spacers)
    if broken is None:
        broken = any(len(spacer) > 140 for spacer in spacers)
    record = FuzzySourcePathProvenance(
        FUZZY_SOURCE_PATH_SCHEMA,
        FuzzySearch._source_family_id(MOTIF, geometry, spacers),
        MOTIF,
        geometry,
        spacers,
        tuple(item.occurrence_id for item in geometry),
        "local_chain_family_parent" if broken else "complete_fuzzy_search",
        (0, len(geometry)),
    )
    candidate.source_path_provenance = (record,)
    candidate.source_path_classification_eligible = not broken
    if broken:
        candidate.source_path_classification_reason = BROKEN_REASON
    return candidate


def synthetic(count=9, missing=(3,), variants=None, identical=False, spacer_length=32):
    rng = random.Random(617)
    variants = variants or {}
    source = "ACGT" * 8
    intervals = []
    for i in range(count):
        start = len(source)
        source += variants.get(i, MOTIF)
        intervals.append((start, len(source)))
        source += (
            "AGCT" * (spacer_length // 4)
            if identical
            else "".join(rng.choice("ACGT") for _ in range(spacer_length))
        )
    selected = [
        interval for index, interval in enumerate(intervals) if index not in missing
    ]
    return source, lineage(build(source, None, *zip(*selected))), intervals


def enumerate_graph(nodes, start, end, lower, upper):
    nodes = sorted(
        {tuple(node["interval"]): node for node in nodes}.values(),
        key=lambda node: node["interval"],
    )
    results = []

    def walk(previous, path, cost):
        if path and lower <= end - previous <= upper:
            results.append((path, cost))
        for index, node in enumerate(nodes):
            if lower <= node["interval"][0] - previous <= upper:
                walk(node["interval"][1], path + (index,), cost + node["edit_cost"])

    walk(start, (), 0)
    return nodes, results


class InternalGraphTests(unittest.TestCase):
    def graph(self, nodes, start=0, end=50, lower=2, upper=16, policy=None):
        return analyze_gap_graph(
            InternalSearchContext("A" * (end + 1), policy),
            nodes,
            start,
            end,
            lower,
            upper,
        )

    def test_dp_matches_independent_exhaustive_small_graphs(self):
        rng = random.Random(891)
        for _ in range(60):
            nodes = [
                {
                    "interval": (start, start + rng.randint(2, 5)),
                    "edit_cost": rng.randrange(4),
                }
                for start in sorted(rng.sample(range(3, 43), 9))
            ]
            ordered, paths = enumerate_graph(nodes, 0, 50, 2, 16)
            result = self.graph(list(reversed(nodes)))
            self.assertTrue(result["complete"])
            self.assertEqual(len(paths), int(result["path_count"]))
            self.assertEqual(
                sorted({len(path) for path, _ in paths}), result["copy_counts"]
            )
            for count_result in result["count_results"]:
                subset = [
                    (path, cost)
                    for path, cost in paths
                    if len(path) == count_result["copies"]
                ]
                minimum = min(cost for _, cost in subset)
                self.assertEqual(minimum, count_result["minimum_cost"])
                self.assertEqual(
                    sum(cost == minimum for _, cost in subset),
                    int(count_result["minimum_path_count"]),
                )
            if result["reason"] == "unique_minimum_tiling":
                best = min(paths, key=lambda item: item[1])[0]
                self.assertEqual(
                    tuple(ordered[i]["interval"] for i in best),
                    result["selected_intervals"],
                )

    def test_higher_cost_copy_count_still_blocks_and_duplicate_geometry_does_not_multiply(
        self,
    ):
        nodes = [
            {"interval": (5, 7), "edit_cost": 0},
            {"interval": (11, 13), "edit_cost": 5},
        ]
        result = self.graph(nodes, end=20, lower=4, upper=14)
        self.assertEqual([1, 2], result["copy_counts"])
        self.assertEqual("copy_count_ambiguity", result["reason"])
        duplicate = self.graph(nodes + [dict(nodes[0])], end=20, lower=4, upper=14)
        self.assertEqual(result, duplicate)

    def test_exponential_exact_counts_unique_then_tied_and_big_integers(self):
        def layered(count, tie=False):
            return [
                {
                    "interval": (10 * i + shift, 10 * i + shift + 2),
                    "edit_cost": 0 if shift == 0 or tie else 1,
                }
                for i in range(1, count + 1)
                for shift in (0, 1)
            ]

        result = self.graph(layered(10), end=112, lower=7, upper=11)
        self.assertEqual("1024", result["path_count"])
        self.assertEqual("unique_minimum_tiling", result["reason"])
        self.assertFalse(result["diagnostic_paths_complete"])
        tie = self.graph(layered(10, True), end=112, lower=7, upper=11)
        self.assertEqual("1024", tie["minimum_path_count"])
        self.assertEqual("minimum_cost_geometry_ambiguity", tie["reason"])
        self.assertEqual(2, len(tie["paths"]))
        big = self.graph(layered(70, True), end=712, lower=7, upper=11)
        self.assertEqual(str(2**70), json.loads(json.dumps(big))["path_count"])
        small_diagnostics = self.graph(
            layered(10, True),
            end=112,
            lower=7,
            upper=11,
            policy=replace(InternalSearchPolicy(), max_gap_paths=2),
        )
        self.assertEqual(tie, small_diagnostics)

    def test_actual_work_ceilings_exact_limit_and_one_beyond(self):
        nodes = [
            {"interval": (5, 7), "edit_cost": 0},
            {"interval": (11, 13), "edit_cost": 1},
        ]
        baseline = self.graph(nodes, end=20, lower=4, upper=14)
        for name, count in baseline["graph_work"].items():
            for scope in ("gap", "parent", "source"):
                with self.subTest(name=name, scope=scope):
                    policy = replace(
                        InternalSearchPolicy(), **{"max_" + scope + "_" + name: count}
                    )
                    self.assertTrue(
                        self.graph(nodes, end=20, lower=4, upper=14, policy=policy)[
                            "complete"
                        ]
                    )
                    policy = replace(policy, **{"max_" + scope + "_" + name: count - 1})
                    rejected = self.graph(
                        nodes, end=20, lower=4, upper=14, policy=policy
                    )
                    self.assertFalse(rejected["complete"])
                    self.assertEqual(scope + "_" + name + "_budget", rejected["reason"])


class InternalSourceTests(unittest.TestCase):
    def setUp(self):
        self.structural = build_advanced_fuzzy_search_filter(PARAMETERS)

    def proposals(self, source, parents, policy=None):
        return propose_internal_repairs(
            InternalSearchContext(source, policy), parents, self.structural
        )

    def test_frozen_fixed_parent_geometric_replays_do_not_invent_eligibility(self):
        raw = (
            Path(__file__).parent / "data/internal_repeat_completion.json"
        ).read_bytes()
        self.assertEqual(
            "c68d40f661aee10898199222a8ea293e306e277d2ddce61445bdc67fabee363d",
            hashlib.sha256(raw).hexdigest(),
        )
        for fixture in json.loads(raw):
            source = fixture["source_fragment"]
            self.assertEqual(
                fixture["source_fragment_sha256"],
                hashlib.sha256(source.encode()).hexdigest(),
            )
            offset = fixture["source_offset_zero_based"]
            intervals = [
                (start - offset, end - offset)
                for start, end in fixture["parent"]["repeat_intervals"]
            ]
            parent = build(source, None, *zip(*intervals))
            fixed = copy.deepcopy(parent)
            corrections, _, _ = propose_terminal_repairs(
                source, [parent], self.structural
            )
            if corrections:
                parent = corrections[0]
            child, record = analyze_internal_geometry(
                InternalSearchContext(source), parent, self.structural
            )
            self.assertEqual(fixed, build(source, None, *zip(*intervals)))
            self.assertEqual(
                "not_authorized_by_geometry", record["eligibility_decision"]
            )
            if fixture["accession"] == "NC_005140.1":
                self.assertIsNone(child)
                self.assertEqual("minimum_cost_geometry_ambiguity", record["reason"])
                gap = record["gaps"][0]
                self.assertEqual(
                    ([2], "529", 2, "9"),
                    (
                        gap["copy_counts"],
                        gap["path_count"],
                        gap["minimum_cost"],
                        gap["minimum_path_count"],
                    ),
                )
                continue
            self.assertEqual(
                fixture["expected_occurrences_inside_parent_envelope"],
                len(child.list_repeats),
            )
            self.assertTrue(
                set(zip(*source_geometry(parent))) <= set(zip(*source_geometry(child)))
            )
            self.assertTrue(self.structural.evidence(child).passes)
            if len(parent.list_repeats) == 19:
                self.assertEqual("390625", record["gaps"][0]["path_count"])
            else:
                self.assertEqual(45, len(record["gaps"]))
            # Source-only fixture is never silently relabelled with an original
            # eligibility marker; a terminal-attested diagnostic is distinct.
            if not corrections:
                self.assertEqual([], self.proposals(source, [parent])[0])

    def test_exact_substitution_internal_insertion_deletion_and_multiple_gaps(self):
        variants = [
            MOTIF,
            MOTIF[:12] + "A" + MOTIF[13:],
            MOTIF[:14] + "T" + MOTIF[14:],
            MOTIF[:14] + MOTIF[15:],
        ]
        for variant in variants:
            with self.subTest(variant=variant):
                source, parent, intervals = synthetic(
                    missing=(3, 6), variants={3: variant}
                )
                original = copy.deepcopy(parent)
                children, records = self.proposals(source, [parent])
                self.assertEqual(1, len(children), records)
                self.assertEqual(intervals, list(zip(*source_geometry(children[0]))))
                self.assertEqual(original, parent)
                self.assertEqual(
                    original.source_path_provenance, parent.source_path_provenance
                )
                self.assertFalse(children[0].source_path_classification_eligible)
                self.assertTrue(
                    all(
                        gap["gap"][1] - gap["gap"][0] < 140
                        for gap in records[0]["gaps"]
                    )
                )
                replay = replay_internal_attestation(
                    source, json.loads(json.dumps(records[0])), self.structural
                )
                self.assertEqual(source_geometry(children[0]), source_geometry(replay))

    def test_substitution_deletion_boundary_tie_abstains_regardless_order_or_score(
        self,
    ):
        source, parent, intervals = synthetic(
            missing=(3,), variants={3: "G" + MOTIF[1:]}
        )
        duplicate = copy.deepcopy(parent)
        parent.score = 0.9
        duplicate.score = 0.51
        for candidates in ([parent], [parent, duplicate], [duplicate, parent]):
            children, records = self.proposals(source, candidates)
            self.assertEqual([], children)
            self.assertEqual("minimum_cost_geometry_ambiguity", records[0]["reason"])
            self.assertGreaterEqual(int(records[0]["gaps"][0]["minimum_path_count"]), 2)

    def test_unconnected_hit_unresolved_second_gap_duplicate_spacers_and_divergence(
        self,
    ):
        cases = [
            synthetic(missing=(3, 6), variants={6: "AC" * 18}),
            synthetic(missing=(3,), identical=True),
            synthetic(
                missing=(3,), variants={0: "AC" * 18, 1: "GT" * 18, 2: "AG" * 18}
            ),
        ]
        source, parent, intervals = synthetic()
        insertion = intervals[3][0]
        shifted_source = source[:insertion] + "A" * 100 + source[insertion:]
        shifted = tuple(
            tuple(x + (100 if x >= insertion else 0) for x in coordinates)
            for coordinates in source_geometry(parent)
        )
        cases.append(
            (shifted_source, lineage(build(shifted_source, None, *shifted)), None)
        )
        for source, parent, _ in cases:
            original = copy.deepcopy(parent)
            children, records = self.proposals(source, [parent])
            self.assertEqual([], children)
            self.assertEqual(original, parent)
            self.assertEqual(
                original.source_path_classification_eligible,
                parent.source_path_classification_eligible,
            )

    def test_unknown_malformed_missing_and_stale_lineage_never_launder(self):
        source, parent, _ = synthetic(missing=(3, 4))
        self.assertFalse(parent.source_path_classification_eligible)
        self.assertEqual(1, len(self.proposals(source, [parent])[0]))
        mutations = [
            lambda p: setattr(p, "source_path_classification_eligible", True),
            lambda p: setattr(p, "source_path_classification_eligible", None),
            lambda p: setattr(p, "source_path_classification_reason", "unknown"),
            lambda p: setattr(
                p, "source_path_restriction_reasons", (BROKEN_REASON, "unknown")
            ),
            lambda p: setattr(p, "source_path_provenance", ()),
            lambda p: setattr(
                p,
                "source_path_provenance",
                (p.source_path_provenance[0]._replace(family_id="stale"),),
            ),
            lambda p: setattr(
                p,
                "source_path_provenance",
                (p.source_path_provenance[0]._replace(source_index_span=(1, 3)),),
            ),
        ]
        for mutation in mutations:
            bad = copy.deepcopy(parent)
            mutation(bad)
            for candidates in ([bad], [parent, bad], [bad, parent]):
                self.assertEqual([], self.proposals(source, candidates)[0])
        with self.assertRaises(ValueError):
            analyze_internal_geometry(
                InternalSearchContext(source[:50]), parent, self.structural
            )
        bad = copy.deepcopy(parent)
        bad.list_spacers[0] = "A" * len(bad.list_spacers[0])
        self.assertEqual([], self.proposals(source, [parent, bad])[0])

    def test_caps_canonical_source_and_cache_counters(self):
        source, parent, _ = synthetic()
        context = InternalSearchContext(source)
        child, record = analyze_internal_geometry(context, parent, self.structural)
        gap = record["gaps"][0]
        cases = [
            ("max_gap_bases", gap["gap"][1] - gap["gap"][0]),
            ("max_parent_gap_bases", record["searched_gap_bases"]),
            ("max_parent_interval_tests", record["interval_tests"]),
            ("max_source_interval_tests", record["interval_tests"]),
            ("max_gap_nodes", len(gap["nodes"])),
        ]
        for name, count in cases:
            exact = InternalSearchContext(
                source, replace(InternalSearchPolicy(), **{name: count})
            )
            self.assertIsNotNone(
                analyze_internal_geometry(exact, parent, self.structural)[0], name
            )
            small = InternalSearchContext(
                source, replace(InternalSearchPolicy(), **{name: count - 1})
            )
            self.assertIsNone(
                analyze_internal_geometry(small, parent, self.structural)[0], name
            )
        before = context.stats["distinct_interval_tests"]
        replay, replay_record = analyze_internal_geometry(
            context, parent, self.structural
        )
        self.assertEqual(record, replay_record)
        self.assertEqual(before, context.stats["distinct_interval_tests"])
        self.assertEqual(1, context.stats["gap_cache_hits"])
        source = source[:40] + "N" + source[41:]
        invalid = build(source, None, *source_geometry(parent))
        with self.assertRaisesRegex(ValueError, "canonical"):
            analyze_internal_geometry(
                InternalSearchContext(source), invalid, self.structural
            )

    def test_pending_only_refinement_preserves_restrictions_on_duplicate_normalization(
        self,
    ):
        source, parent, _ = synthetic(missing=(3, 4))
        core = build(
            source, None, parent.list_repeat_starts[:3], source_geometry(parent)[1][:3]
        )
        refinement = DetectionRefinement(
            {Locus(0, len(source)): [fuzzy(core), fuzzy(parent)]},
            PARAMETERS,
            False,
            source,
        )
        self.assertEqual(1, len(refinement.pending_internal_repairs))
        child_geometry = source_geometry(refinement.pending_internal_repairs[0])
        self.assertFalse(
            any(
                source_geometry(candidate) == child_geometry
                for values in refinement.output().values()
                for candidate in values
            )
        )
        bad = fuzzy(parent)
        bad.source_path_classification_reason = "additional_external_restriction"
        for candidates in ([fuzzy(parent), bad], [bad, fuzzy(parent)]):
            refinement = DetectionRefinement(
                {Locus(0, len(source)): candidates}, PARAMETERS, False, source
            )
            self.assertEqual([], refinement.pending_internal_repairs)

    def test_boundary_then_internal_composition_and_genome_context_across_refinement(
        self,
    ):
        from tests.test_repeat_boundary_repair import synthetic as boundary_source

        source, full = boundary_source(count=9)
        keep = [index for index in range(9) if index != 2]
        geometry = source_geometry(full)
        parent = lineage(
            build(
                source,
                None,
                [geometry[0][i] for i in keep],
                [geometry[1][i] for i in keep],
            )
        )
        corrected, records, _ = propose_terminal_repairs(
            source, [parent], self.structural
        )
        self.assertEqual(1, len(corrected))
        pending, audits = self.proposals(source, corrected)
        self.assertEqual(1, len(pending), audits)
        self.assertEqual(9, len(pending[0].list_repeats))
        self.assertEqual(
            "terminal_repair", audits[0]["parent_authorizations"][0]["kind"]
        )
        replay_internal_attestation(source, audits[0], self.structural)

        first_source, first, _ = synthetic()
        second_source, second, _ = synthetic()
        source = first_source + second_source
        shifted = tuple(
            tuple(x + len(first_source) for x in axis)
            for axis in source_geometry(second)
        )
        second = lineage(build(source, None, *shifted))
        probe = InternalSearchContext(source)
        analyze_internal_geometry(probe, first, self.structural)
        context = InternalSearchContext(
            source,
            replace(
                InternalSearchPolicy(),
                max_source_interval_tests=probe.stats["distinct_interval_tests"],
            ),
        )
        first_refinement = DetectionRefinement(
            {Locus(0, len(first_source)): [fuzzy(first)]},
            PARAMETERS,
            False,
            source,
            context,
        )
        self.assertEqual(1, len(first_refinement.pending_internal_repairs))
        repeated = DetectionRefinement(
            {Locus(0, len(first_source)): [fuzzy(first)]},
            PARAMETERS,
            False,
            source,
            context,
        )
        self.assertEqual(1, len(repeated.pending_internal_repairs))
        next_refinement = DetectionRefinement(
            {Locus(len(first_source), len(source)): [fuzzy(second)]},
            PARAMETERS,
            False,
            source,
            context,
        )
        self.assertEqual([], next_refinement.pending_internal_repairs)
        self.assertEqual(
            "source_interval_budget",
            next_refinement.internal_repair_records[0]["reason"],
        )
        self.assertEqual(
            probe.stats["distinct_interval_tests"],
            context.stats["distinct_interval_tests"],
        )

    def test_component_parent_budget_is_deterministic_and_geometry_deduplicated(self):
        source, whole, _ = synthetic(count=18, missing=(7,))
        starts, ends = source_geometry(whole)
        parents = [
            lineage(build(source, None, starts[:count], ends[:count]))
            for count in range(10, 18)
        ]
        policy = replace(InternalSearchPolicy(), max_parents_per_component=2)
        first, records = self.proposals(source, parents + [parents[0]], policy)
        second, reverse = self.proposals(source, list(reversed(parents)), policy)
        self.assertEqual(
            [source_geometry(x) for x in first], [source_geometry(x) for x in second]
        )
        self.assertEqual(2, len(first))
        self.assertEqual(
            6, sum(record["reason"] == "component_parent_budget" for record in records)
        )
        self.assertEqual(records, reverse)


class InternalTransactionTests(unittest.TestCase):
    setUp = InternalSourceTests.setUp
    proposals = InternalSourceTests.proposals

    def prepared(self):
        source, parent, _ = synthetic(missing=(3, 4))
        pending, records = self.proposals(source, [parent])
        core = build(
            source, None, parent.list_repeat_starts[:3], source_geometry(parent)[1][:3]
        )
        old = [0.91, core, [], {"classification_basis": {"method": "model_score_band"}}]
        categories = [{(1, 2): [old]}, {}, {}, {}, {(3, 4): [[0.88, parent, []]]}]
        return source, parent, pending, categories, old

    def scorer(self, score=0.8):
        return lambda mapping: {
            key: [[score, candidate, []] for candidate in candidates]
            for key, candidates in mapping.items()
        }

    def test_fresh_score_band_single_event_immutable_idempotent_report_roundtrip(self):
        source, parent, pending, categories, old = self.prepared()
        original = copy.deepcopy(categories)
        updated, audit = complete_internal_repeats(
            source, categories, pending, self.structural, self.scorer()
        )
        self.assertEqual(original, categories)
        self.assertFalse(parent.source_path_classification_eligible)
        accepted = [row for i in (0, 2) for rows in updated[i].values() for row in rows]
        self.assertEqual(1, len(accepted))
        row = accepted[0]
        self.assertEqual(9, len(row[1].list_repeats))
        self.assertEqual(0.8, row[0])
        self.assertTrue(
            any(
                item[0] == old[0] and item[1] is old[1]
                for rows in updated[1].values()
                for item in rows
            )
        )
        self.assertEqual(METHOD, row[3]["classification_basis"]["method"])
        again, _ = complete_internal_repeats(
            source, updated, pending, self.structural, self.scorer()
        )
        self.assertEqual(updated, again)
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "source.fasta"
            fasta.write_text(">source\n" + source + "\n")
            model = ReportModelBuilder(
                file_path=fasta,
                parameters=PARAMETERS,
                flags={},
                categories=updated,
                non_array_data={},
                list_feature_names=[],
            ).build()
            restored = json.loads(json.dumps(model))
        self.assertIn(METHOD, json.dumps(restored))
        basis = row[3]["classification_basis"]["internal_repeat_completion"]
        for field, value in (
            ("completed_score", 0.9),
            ("source_sha256", "stale"),
            ("core_geometry", source_geometry(parent)),
        ):
            bad = copy.deepcopy(basis)
            bad[field] = value
            with self.assertRaises(ValueError):
                replay_internal_completion_basis(source, row[1], bad, row[0], 0)

    def test_score_failure_cross_key_neighbor_and_stale_score_cannot_commit(self):
        source, parent, pending, categories, old = self.prepared()
        updated, _ = complete_internal_repeats(
            source, categories, pending, self.structural, self.scorer(0.74)
        )
        self.assertEqual(categories, updated)
        neighbor = build(
            source,
            None,
            parent.list_repeat_starts[-3:],
            source_geometry(parent)[1][-3:],
        )
        categories[2][(500, 600)] = [[0.7, neighbor, []]]
        updated, audit = complete_internal_repeats(
            source, categories, pending, self.structural, self.scorer()
        )
        self.assertEqual(categories, updated)
        self.assertTrue(
            any(item["reason"] == "other_accepted_overlap" for item in audit)
        )

        def stale(mapping):
            return {
                key: [[0.9, copy.deepcopy(candidate), []] for candidate in candidates]
                for key, candidates in mapping.items()
            }

        with self.assertRaisesRegex(ValueError, "identity"):
            complete_internal_repeats(
                source, categories, pending, self.structural, stale
            )

    def test_scorer_mutation_and_tampered_attestation_abort(self):
        source, parent, pending, categories, old = self.prepared()

        def mutate(mapping):
            old[0] = 0.99
            return self.scorer()(mapping)

        with self.assertRaisesRegex(ValueError, "snapshot"):
            complete_internal_repeats(
                source, categories, pending, self.structural, mutate
            )
        old[0] = 0.91
        bad = copy.deepcopy(pending)
        bad[0].internal_repair_record["gaps"][0]["minimum_cost"] = 100
        updated, audit = complete_internal_repeats(
            source, categories, bad, self.structural, self.scorer()
        )
        self.assertEqual(categories, updated)
        self.assertEqual("invalid_attestation", audit[0]["reason"])


if __name__ == "__main__":
    unittest.main()
