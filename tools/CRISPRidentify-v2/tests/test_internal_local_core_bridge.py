import copy
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from components.components_evaluated_arrays_enhancement import (
    build_source_backed_candidate_from_boundaries as build,
)
from components.internal_completion_authority import (
    capture_current_authority_inventory,
    current_by_geometry,
)
from components.internal_local_core_bridge import (
    _restore,
    _digest,
    collect_bridge_inventory,
    compose_local_core,
    replay_local_core,
    validate_declared_family_against_run,
)
from components.internal_repeat_completion import (
    InternalSearchContext,
    InternalSearchPolicy,
    analyze_gap_graph,
)
from components.internal_repeat_transaction import complete_internal_repeats
from components.module_evaluation import build_advanced_fuzzy_search_filter
from components.repeat_boundary_repair import source_geometry
from components.report_model import ReportModelBuilder
from components.report_provenance_validation import (
    ReportProvenanceError,
    ReportProvenanceValidationContext,
    validate_classification_basis,
)
from tests import test_internal_repeat_completion as source_tests


class InternalRevisionTests(unittest.TestCase):
    def test_all_restricted_child_duplicates_veto_both_routes_and_preserve_evidence(
        self,
    ):
        t = source_tests.InternalTransactionTests()
        t.setUp()
        for route in ("pending", "current"):
            for marker, reason in (
                (False, "external_unknown_restriction"),
                (None, "pending_internal_completion"),
                (1, "pending_internal_completion"),
                (True, "external_unknown_restriction"),
            ):
                for reverse in (False, True):
                    with self.subTest(route=route, marker=marker, reverse=reverse):
                        source, parent, pending, categories, old = t.prepared()
                        restricted = copy.deepcopy(pending[0])
                        restricted.source_path_classification_eligible = marker
                        restricted.source_path_classification_reason = reason
                        if route == "pending":
                            pending.append(restricted)
                            if reverse:
                                pending.reverse()
                        else:
                            categories[4][(5000, 6000)] = [[0.88, restricted, []]]
                            if reverse:
                                categories[4] = dict(
                                    reversed(list(categories[4].items()))
                                )
                        original = copy.deepcopy(categories)
                        updated, audit = complete_internal_repeats(
                            source, categories, pending, t.structural, t.scorer()
                        )
                        self.assertEqual(original, updated)
                        self.assertFalse(
                            any(item["status"] == "committed" for item in audit)
                        )
                        again, _ = complete_internal_repeats(
                            source, updated, pending, t.structural, t.scorer()
                        )
                        self.assertEqual(updated, again)
                        if route == "current":
                            self.assertTrue(
                                any(
                                    row[1] is restricted
                                    for mapping in updated
                                    for rows in mapping.values()
                                    for row in rows
                                )
                            )
                        with tempfile.TemporaryDirectory() as tmp:
                            fasta = Path(tmp) / "source.fna"
                            fasta.write_text(">source\n" + source + "\n")
                            report = ReportModelBuilder(
                                file_path=fasta,
                                parameters={},
                                flags={},
                                categories=updated,
                                non_array_data={},
                                list_feature_names=[],
                            ).build()
                            self.assertEqual(
                                1,
                                json.loads(json.dumps(report))["summary"][
                                    "accepted_array_count"
                                ],
                            )

    def test_pending_unknown_lineage_and_scoring_reason_mutation_fail_closed(self):
        t = source_tests.InternalTransactionTests()
        t.setUp()
        source, parent, pending, categories, old = t.prepared()
        malformed = copy.deepcopy(pending[0])
        malformed.source_path_provenance = None
        updated, _ = complete_internal_repeats(
            source, categories, pending + [malformed], t.structural, t.scorer()
        )
        self.assertEqual(categories, updated)

        def mutate(mapping):
            old[1].source_path_classification_reason = "new_unknown_restriction"
            return t.scorer()(mapping)

        with self.assertRaisesRegex(ValueError, "snapshot"):
            complete_internal_repeats(source, categories, pending, t.structural, mutate)

    def test_diagnostic_completeness_after_actual_slice(self):
        nodes = [
            {"interval": (5, 7), "edit_cost": 1},
            {"interval": (6, 8), "edit_cost": 1},
        ]
        for cap in (1, 2):
            context = InternalSearchContext(
                "A" * 30, replace(InternalSearchPolicy(), max_gap_paths=cap)
            )
            result = analyze_gap_graph(context, nodes, 0, 20, 4, 14)
            self.assertTrue(result["complete"])
            self.assertEqual("minimum_cost_geometry_ambiguity", result["reason"])
            self.assertEqual(cap, len(result["paths"]))
            self.assertEqual(cap == 2, result["diagnostic_paths_complete"])
        # A more costly copy-count alternative is still part of the complete
        # family even when only minimum witnesses are printed.
        result = analyze_gap_graph(
            InternalSearchContext("A" * 30),
            [
                {"interval": (5, 7), "edit_cost": 0},
                {"interval": (11, 13), "edit_cost": 5},
            ],
            0,
            20,
            4,
            14,
        )
        self.assertEqual([1, 2], result["copy_counts"])
        self.assertFalse(result["diagnostic_paths_complete"])


class NativeLocalCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture = json.load(
            gzip.open(Path(__file__).parent / "data/internal_local_core_native.json.gz")
        )
        cls.fixture = fixture
        cls.source = fixture["source"]
        cls.structural = build_advanced_fuzzy_search_filter(fixture["parameters"])
        assert (
            hashlib.sha256(cls.source.encode()).hexdigest() == fixture["source_sha256"]
        )
        categories = [{}, {}, {}, {}, {}]
        for row in fixture["rows"]:
            item = row["candidate"]
            candidate = _restore(
                cls.source,
                {
                    "geometry": item["geometry"],
                    "authority": {k: v for k, v in item.items() if k != "geometry"},
                },
            )
            categories[row["category"]].setdefault(tuple(row["key"]), []).append(
                [row["score"], candidate, []]
            )
        pending = []
        for item in fixture["pending"]:
            pending.append(
                _restore(
                    cls.source,
                    {
                        "geometry": item["geometry"],
                        "authority": {k: v for k, v in item.items() if k != "geometry"},
                    },
                )
            )
        cls.categories = categories
        cls.pending = pending
        cls.record = pending[0].internal_repair_record
        cls.target = source_geometry(pending[0])
        cls.core = source_geometry(next(iter(categories[0].values()))[0][1])
        cls.inventory = collect_bridge_inventory(
            cls.target,
            fixture["boundary_records"],
            current_by_geometry(categories),
            pending,
        )
        cls.proof = compose_local_core(
            cls.source,
            cls.core,
            cls.record,
            cls.inventory,
            cls.structural,
            InternalSearchContext(cls.source),
        )

    def compose(self, inventory=None, record=None, core=None):
        return compose_local_core(
            self.source,
            self.core if core is None else core,
            self.record if record is None else record,
            self.inventory if inventory is None else inventory,
            self.structural,
            InternalSearchContext(self.source),
        )

    def run_evidence(self, categories=None):
        return {
            "boundary_records": self.fixture["boundary_records"],
            "internal_records": [
                candidate.internal_repair_record for candidate in self.pending
            ],
            "current_authority_inventory": capture_current_authority_inventory(
                self.source,
                self.categories if categories is None else categories,
                self.pending,
            ),
        }

    def test_complete_native_family_map_transaction_report_and_exact_cache(self):
        self.assertEqual(131, len(self.inventory["current_rows"]))
        self.assertEqual(6, len(self.proof["geometry_family"]))
        self.assertAlmostEqual(3752.086690116019, self.proof["phase"]["winner_margin"])
        self.assertTrue(
            all(item["parent_offset"] == 107 for item in self.proof["maps"])
        )
        self.assertTrue(
            all(item["witness_offset"] == 108 for item in self.proof["maps"])
        )
        self.assertTrue(
            all(
                tuple(map(tuple, item["child_geometry"])) != self.target
                for item in self.proof["witnesses"]
            )
        )
        context = InternalSearchContext(self.source)
        replay_local_core(
            self.source, self.core, self.record, self.proof, self.structural, context
        )
        comparisons = context.stats["local_core_phase_comparisons"]
        replay_local_core(
            self.source,
            self.core,
            self.record,
            copy.deepcopy(self.proof),
            self.structural,
            context,
        )
        self.assertEqual(comparisons, context.stats["local_core_phase_comparisons"])
        self.assertEqual(1, context.stats["local_core_proof_cache_hits"])
        categories = copy.deepcopy(self.categories)
        old = next(iter(categories[0].values()))[0]
        old.append({"classification_basis": {"method": "model_score_band"}})
        original = copy.deepcopy(categories)

        def score(mapping):
            return {
                key: [[0.8584082272961254, candidate, []] for candidate in values]
                for key, values in mapping.items()
            }

        updated, audit = complete_internal_repeats(
            self.source,
            categories,
            self.pending,
            self.structural,
            score,
            boundary_records=self.fixture["boundary_records"],
        )
        self.assertEqual(original, categories)
        accepted = [row for i in (0, 2) for rows in updated[i].values() for row in rows]
        self.assertEqual(1, len(accepted))
        self.assertEqual(self.target, source_geometry(accepted[0][1]))
        evidence = accepted[0][3]["classification_basis"]["internal_repeat_completion"]
        self.assertEqual(
            [1],
            [item["category_index"] for item in evidence["previous_completed_rows"]],
        )
        self.assertEqual({}, evidence["previous_completed_rows"][0]["authority"])
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "source.fna"
            fasta.write_text(">source\n" + self.source + "\n")
            parameters = dict(
                self.fixture["parameters"],
                param_repeat_boundary_repair={
                    "records": self.fixture["boundary_records"]
                },
                param_internal_repeat_completion={
                    "records": [
                        candidate.internal_repair_record for candidate in self.pending
                    ],
                    "current_authority_inventory": self.run_evidence(categories)[
                        "current_authority_inventory"
                    ],
                },
            )
            model = ReportModelBuilder(
                file_path=fasta,
                parameters=parameters,
                flags={},
                categories=updated,
                non_array_data={},
                list_feature_names=[],
            ).build()
            restored = json.loads(json.dumps(model))
        self.assertEqual(1, restored["summary"]["accepted_array_count"])
        # A cached valid source replay cannot bypass a changed run-level census,
        # including when replaying just the displaced core's retained role.
        run = self.run_evidence(categories)
        validation = ReportProvenanceValidationContext.from_source(
            self.source,
            self.fixture["source_sha256"],
            categories=updated,
            internal_completion_evidence=run,
        )
        rows = [
            (i, row)
            for i, mapping in enumerate(updated)
            for group in mapping.values()
            for row in group
            if len(row) > 3
            and "internal_repeat_completion" in row[3].get("classification_basis", {})
        ]
        self.assertEqual(
            {False, True},
            {
                row[3]["classification_basis"].get("retained_core", False)
                for _, row in rows
            },
        )

        def validate(category, row):
            return validate_classification_basis(
                self.source,
                self.fixture["source_sha256"],
                row[1],
                row[3]["classification_basis"],
                row_score=row[0],
                category_index=category,
                allowed_methods={
                    "source_backed_internal_repeat_completion",
                    "model_score_band",
                },
                validation_context=validation,
            )

        for category, row in rows:
            validate(category, row)
        run["current_authority_inventory"]["rows"].pop()
        for category, row in rows:
            with self.assertRaisesRegex(ReportProvenanceError, "current-family member"):
                validate(category, row)
        again, _ = complete_internal_repeats(
            self.source,
            updated,
            self.pending,
            self.structural,
            score,
            boundary_records=self.fixture["boundary_records"],
        )
        self.assertEqual(updated, again)

    def test_every_unknown_current_role_restricts_and_missing_inventory_fails(self):
        geometries = [
            self.target,
            self.core,
            tuple(map(tuple, self.record["parent_geometry"])),
        ]
        witness = self.proof["witnesses"][0]
        geometries.extend(
            [
                tuple(map(tuple, witness[name]))
                for name in ("parent_geometry", "child_geometry")
            ]
        )
        for geometry in geometries:
            bad = copy.deepcopy(self.inventory)
            bad["current_rows"].append(
                {
                    "geometry": geometry,
                    "authority": {
                        "source_path_classification_eligible": False,
                        "source_path_classification_reason": "external_unknown_restriction",
                    },
                }
            )
            with self.assertRaises(ValueError):
                self.compose(bad)
        bad = copy.deepcopy(self.inventory)
        bad["complete"] = False
        with self.assertRaisesRegex(ValueError, "incomplete"):
            self.compose(bad)
        bad = copy.deepcopy(self.inventory)
        bad["pending_records"][0]["authority"][
            "source_path_classification_eligible"
        ] = None
        with self.assertRaises(ValueError):
            self.compose(bad)

    def test_full_phase_abstention_other_winner_and_outside_rule_cannot_use_local_corridor(
        self,
    ):
        for reason in (
            "source_score_tie",
            "weak_source_score_margin",
            "weak_source_boundary_evidence",
        ):
            with patch(
                "components.internal_local_core_bridge.compare_source_backed_phase_candidates",
                return_value={"status": "abstained", "reason": reason},
            ):
                with self.assertRaisesRegex(ValueError, "did not select"):
                    self.compose()
        other = build(self.source, None, *self.proof["geometry_family"][0])
        if source_geometry(other) == self.target:
            other = build(self.source, None, *self.proof["geometry_family"][1])
        with patch(
            "components.internal_local_core_bridge.compare_source_backed_phase_candidates",
            return_value={
                "status": "selected",
                "reason": "unique_source_phase_winner",
                "winner_candidate": other,
            },
        ):
            with self.assertRaisesRegex(ValueError, "did not select"):
                self.compose()
        bad = copy.deepcopy(self.inventory)
        shifted = (
            tuple(s + 13 for s in self.target[0]),
            tuple(e + 13 for e in self.target[1]),
        )
        bad["current_rows"].append({"geometry": shifted, "authority": {}})
        with self.assertRaisesRegex(ValueError, "outside the bounded"):
            self.compose(bad)

    def test_family_cap_order_stale_witness_and_declared_member_tampering(self):
        bad = copy.deepcopy(self.inventory)
        for index in range(33):
            ends = list(self.target[1])
            ends[index] += 1
            bad["current_rows"].append(
                {"geometry": (self.target[0], tuple(ends)), "authority": {}}
            )
        with self.assertRaisesRegex(ValueError, "32 geometries"):
            self.compose(bad)
        ordered = copy.deepcopy(self.inventory)
        for key in ("current_rows", "pending_records", "boundary_records"):
            ordered[key].reverse()
        # Collection canonicalizes input order; compare decision/map sets.
        proof = self.compose(ordered)
        self.assertEqual(self.proof["geometry_family"], proof["geometry_family"])
        self.assertEqual(
            {_digest(item) for item in self.proof["maps"]},
            {_digest(item) for item in proof["maps"]},
        )
        bad = copy.deepcopy(self.inventory)
        bad["boundary_records"][0]["source_sha256"] = "stale"
        with self.assertRaises(ValueError):
            self.compose(bad)
        for mutation in ("family", "inventory", "map", "phase", "source"):
            proof = copy.deepcopy(self.proof)
            if mutation == "family":
                proof["geometry_family"] = proof["geometry_family"][1:]
            elif mutation == "inventory":
                proof["inventory"]["boundary_records"].pop()
            elif mutation == "map":
                proof["maps"][0]["parent_offset"] += 1
            elif mutation == "phase":
                proof["phase"]["winner_margin"] += 1
            else:
                proof["source_sha256"] = "stale"
            with self.assertRaises(ValueError):
                replay_local_core(
                    self.source,
                    self.core,
                    self.record,
                    proof,
                    self.structural,
                    InternalSearchContext(self.source),
                )

    def test_self_consistent_omitted_family_fails_independent_run_manifest(self):
        bad = copy.deepcopy(self.inventory)
        removable = next(
            item["child_geometry"]
            for item in self.proof["witnesses"]
            if not any(
                mapping["witness_identity"] == item["identity"]
                for mapping in self.proof["maps"]
            )
        )
        bad["boundary_records"] = [
            record
            for record in bad["boundary_records"]
            if tuple(map(tuple, record["child_geometry"]))
            != tuple(map(tuple, removable))
        ]
        proof = self.compose(bad)
        with self.assertRaisesRegex(ValueError, "omits or changes"):
            validate_declared_family_against_run(
                proof, self.target, self.run_evidence()
            )

    def test_current_only_tied_phase_cannot_be_omitted_from_report_family(self):
        categories = copy.deepcopy(self.categories)
        ends = list(self.target[1])
        ends[2] += 1
        contrary = build(self.source, None, self.target[0], ends)
        categories[1].setdefault((0, 1), []).append([0.8584082272961254, contrary, []])
        complete = collect_bridge_inventory(
            self.target,
            self.fixture["boundary_records"],
            current_by_geometry(categories),
            self.pending,
        )
        self.assertEqual(132, len(complete["current_rows"]))
        with self.assertRaisesRegex(ValueError, "source_score_tie"):
            self.compose(complete)
        # The original, self-consistent 131-row proof still selects the target.
        # Its omission must fail against the independently captured 132 rows.
        proof = self.compose()
        with self.assertRaisesRegex(ValueError, "current-family member"):
            validate_declared_family_against_run(
                proof, self.target, self.run_evidence(categories)
            )

    def test_current_census_preserves_duplicates_authority_presence_and_order(self):
        for marker in (None, False, True, 1):
            for role in (self.target, self.core):
                categories = copy.deepcopy(self.categories)
                duplicate = build(self.source, None, *role)
                duplicate.source_path_classification_eligible = marker
                duplicate.source_path_classification_reason = (
                    "external_unknown_restriction"
                )
                duplicate.source_path_restriction_reasons = ["unknown"]
                duplicate.source_path_provenance = None
                categories[4].setdefault((0, 1), []).append([0.1, duplicate, []])
                run = self.run_evidence(categories)
                with self.assertRaisesRegex(ValueError, "current-family member"):
                    validate_declared_family_against_run(self.proof, self.target, run)
                snapshot = copy.deepcopy(run)
                duplicate.source_path_restriction_reasons.append(
                    "changed_after_capture"
                )
                self.assertEqual(snapshot, run)
        # Even an unrestricted duplicate cannot silently disappear.
        run = self.run_evidence()
        run["current_authority_inventory"]["rows"].append(
            copy.deepcopy(run["current_authority_inventory"]["rows"][0])
        )
        with self.assertRaisesRegex(ValueError, "current-family member"):
            validate_declared_family_against_run(self.proof, self.target, run)
        for field, value in (
            ("source_path_classification_eligible", None),
            ("source_path_classification_eligible", True),
            ("source_path_classification_reason", None),
            ("source_path_provenance", None),
        ):
            run = self.run_evidence()
            item = next(
                item
                for item in run["current_authority_inventory"]["rows"]
                if item["authority"] == {}
            )
            item["authority"][field] = value
            with self.assertRaisesRegex(ValueError, "current-family member"):
                validate_declared_family_against_run(self.proof, self.target, run)
        run = self.run_evidence()
        proof = copy.deepcopy(self.proof)
        for key in ("boundary_records", "internal_records"):
            run[key] = list(reversed(run[key]))
        run["current_authority_inventory"]["rows"].reverse()
        run["current_authority_inventory"]["target_geometries"].reverse()
        proof["inventory"]["current_rows"].reverse()
        validate_declared_family_against_run(
            proof, self.target, json.loads(json.dumps(run))
        )

    def test_missing_incomplete_stale_or_changed_current_census_fails_closed(self):
        for mutation in (
            "absent",
            "null",
            "incomplete",
            "source",
            "scope",
            "removed",
            "geometry",
            "policy",
            "rows",
            "targets",
        ):
            run = self.run_evidence()
            current = run["current_authority_inventory"]
            if mutation == "absent":
                run.pop("current_authority_inventory")
            elif mutation == "null":
                run["current_authority_inventory"] = None
            elif mutation == "incomplete":
                current["complete"] = False
            elif mutation == "source":
                current["source_sha256"] = "stale"
            elif mutation == "scope":
                current["target_geometries"] = [self.core]
            elif mutation == "removed":
                current["rows"].pop()
            elif mutation == "geometry":
                current["rows"][0]["geometry"] = self.core
            elif mutation == "rows":
                current.pop("rows")
            elif mutation == "targets":
                current.pop("target_geometries")
            else:
                current["policy"] = "stale"
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                validate_declared_family_against_run(self.proof, self.target, run)

    def test_pipeline_captures_immutable_complete_census_before_transaction(self):
        from components.pipeline import Pipeline

        pipeline = Pipeline.__new__(Pipeline)
        pipeline.parameters = self.fixture["parameters"]
        pipeline.genome_sequence = self.source
        pipeline.categories = copy.deepcopy(self.categories)
        pipeline.pending_internal_repairs = self.pending
        expected = self.run_evidence(pipeline.categories)["current_authority_inventory"]

        def mutate(source, categories, *args, **kwargs):
            for mapping in categories:
                mapping.clear()
            return categories, []

        with (
            patch("components.pipeline.complete_internal_repeats", side_effect=mutate),
            patch.object(pipeline, "_score_candidates_complete"),
            patch.object(pipeline, "_scored_rows_from_categories", return_value={}),
        ):
            pipeline._complete_pending_internal_repeats()
        captured = pipeline._internal_repeat_audit()["current_authority_inventory"]
        self.assertEqual(_digest(expected), _digest(captured))
        self.assertEqual(131, len(captured["rows"]))
        captured["rows"].clear()
        self.assertEqual(
            131,
            len(
                pipeline._internal_repeat_audit()["current_authority_inventory"]["rows"]
            ),
        )

    def test_conflicting_local_maps_abstain_even_if_future_witness_corridor_changes(
        self,
    ):
        # Fault-inject the preceding corridor decision to exercise the final
        # agreement guard independently of today's narrow trim domain.
        parent = tuple(map(tuple, self.record["parent_geometry"]))
        witnesses = [
            {
                "identity": str(offset),
                "parent_geometry": (
                    parent[0][offset : offset + 3],
                    parent[1][offset : offset + 3],
                ),
                "child_geometry": (
                    parent[0][offset : offset + 3],
                    parent[1][offset : offset + 3],
                ),
            }
            for offset in (107, 108)
        ]
        resolved = {
            "current": {},
            "witnesses": witnesses,
            "geometry_family": self.proof["geometry_family"],
            "phase": self.proof["phase"],
        }
        with (
            patch(
                "components.internal_local_core_bridge.resolve_bridge_family",
                return_value=resolved,
            ),
            patch("components.internal_local_core_bridge._corridor", return_value=0),
        ):
            with self.assertRaisesRegex(ValueError, "unique agreeing"):
                self.compose()

    def test_local_map_cannot_use_inserted_copy_and_accepted_neighbor_still_blocks(
        self,
    ):
        bad = copy.deepcopy(self.record)
        parent = tuple(map(tuple, bad["parent_geometry"]))
        bad["parent_geometry"] = (
            parent[0][:108] + parent[0][109:],
            parent[1][:108] + parent[1][109:],
        )
        with self.assertRaisesRegex(ValueError, "retained-anchor"):
            self.compose(record=bad)
        categories = copy.deepcopy(self.categories)
        neighbor = build(
            self.source, None, self.target[0][200:203], self.target[1][200:203]
        )
        categories[2][("other", "key")] = [[0.7, neighbor, []]]

        def score(mapping):
            return {
                key: [[0.8584082272961254, c, []] for c in values]
                for key, values in mapping.items()
            }

        updated, audit = complete_internal_repeats(
            self.source,
            categories,
            self.pending,
            self.structural,
            score,
            boundary_records=self.fixture["boundary_records"],
        )
        self.assertEqual(categories, updated)
        self.assertTrue(
            any(item.get("reason") == "other_accepted_overlap" for item in audit)
        )


if __name__ == "__main__":
    unittest.main()
