"""Whole generated paths gain contextual Possible evidence, never flank votes."""
import copy
import hashlib
import unittest
from unittest.mock import patch

from components.candidate_origins import mark_repair_origin
from components.contextual_spacer_rescue import (
    METHOD, OBSERVATION_KEY, _choose, _geometry, _scope_checks,
    contextual_spacer_member_intervals, contextual_spacer_summary,
    rescue_contextual_spacer_arrays, validate_contextual_spacer_basis,
    validate_contextual_spacer_audit, validate_contextual_recorded_metrics, _references, _qualified,
)
from components.initial_competition import initial_competition_groups
from components.localized_spacer_context import (
    capture_pending_context, digest, localized_block,
)
from components.module_evaluation import build_advanced_fuzzy_search_filter
from tests.contextual_spacer_fixture import fixture, load_fixture, local_candidate
from tests.test_candidate_evidence import PARAMETERS


def flatten(categories, index):
    return [row for values in categories[index].values() for row in values]


def mark_repair(candidate, source):
    mark_repair_origin(candidate, {"policy": "terminal-only-trim-v1", "status": "generated",
        "coordinate_system": "zero_based_half_open", "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "child_geometry": [candidate.list_repeat_starts, [b for _a, b in _geometry(candidate)]]})


def replay(source, output, audit):
    for row in flatten(output, 2):
        if len(row) == 4 and row[3].get("classification_basis", {}).get("method") == METHOD:
            validate_contextual_spacer_basis(source, row[1], row[3]["classification_basis"], row[0], 2, row, output)
    report = {"source": {"sha256": hashlib.sha256(source.encode()).hexdigest()},
              "arrays": [{"classification_basis": row[3]["classification_basis"]}
                         for row in flatten(output, 2) if len(row) == 4]}
    validate_contextual_spacer_audit(report, audit, provenance_categories=output, source_sequence=source)


class ContextualSpacerTests(unittest.TestCase):
    def run_fixture(self, **kwargs):
        source, categories, pending, scorer = fixture(**kwargs)
        output, audit = rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, scorer, lambda: ())
        return source, output, audit

    def test_actual_top_tie_all_paths_and_original_failure_remain(self):
        source, output, audit = self.run_fixture(ordinals=[258, 260, 261, 262, 263])
        self.assertEqual(audit["stats"], {"pending_observations": 5, "fresh_candidates_scored": 15, "selections_committed": 1})
        self.assertEqual(len(flatten(output, 4)), 5)
        row = flatten(output, 2)[0]
        self.assertGreater(row[0], .75)
        summary = contextual_spacer_summary(row)
        self.assertEqual(summary["confidence_category"], "Possible")
        self.assertIs(summary["accepted_core_present"], False)
        self.assertEqual(summary["member_count"], 2)
        self.assertEqual(summary["variable_repeat_indexes"], [16])
        self.assertEqual({len(p) for p in contextual_spacer_member_intervals(row)}, {17})
        self.assertEqual(len({tuple(map(tuple, p)) for p in contextual_spacer_member_intervals(row)}), 2)
        records = audit["observations"]
        self.assertNotEqual(records[0]["flanks"]["right"]["fresh"]["score"], records[1]["flanks"]["right"]["fresh"]["score"])
        self.assertEqual(records[0]["fresh"]["score"], records[1]["fresh"]["score"])
        for record in records:
            self.assertEqual([r["predicate"] for r in record["validation"]["structural_filter"]["predicate_results"]
                              if r["status"] == "failed"], ["identical_spacer_burden"])
            self.assertEqual(record["block"]["spacer_run"], [10, 14])
            self.assertEqual([len(record["flanks"][r]["fresh"]["repeat_intervals"]) for r in ("left", "right")], [11, 3])
        replay(source, output, audit)

    def test_actual_negative_phases_and_ineligible_parents_do_not_acquire_exception(self):
        source, output, audit = self.run_fixture(ordinals=[0, 250, 251, 252, 253, 259, 266, 267, 268, 279])
        self.assertEqual(len(flatten(output, 2)), 0)
        self.assertGreater(len(audit["observations"]), 0)
        self.assertTrue(all(not r["qualified"] for r in audit["observations"]))
        self.assertTrue(all(len(r["fresh"]["repeat_intervals"]) >= 3 for r in audit["observations"]))
        replay(source, output, audit)

    def test_aliases_permutations_preserve_geometry_set_and_every_observation(self):
        sets = []
        for ordinals in ([258, 260], [260, 258], [258, 260, 258]):
            source, output, audit = self.run_fixture(ordinals=ordinals)
            row = flatten(output, 2)[0]
            sets.append(set(contextual_spacer_summary(row)["admissible_geometry_ids"]))
            self.assertEqual(len(flatten(output, 4)), len(ordinals))
            self.assertEqual(len({r[3][OBSERVATION_KEY]["observation_id"] for r in flatten(output, 4)}), len(ordinals))
            replay(source, output, audit)
        self.assertEqual(sets, [sets[0]]*3)

    def test_mutated_pending_and_scoring_callbacks_fail_before_output(self):
        source, categories, pending, scorer = fixture()
        duplicate = pending + pending
        with self.assertRaisesRegex(ValueError, "duplicate pending"):
            rescue_contextual_spacer_arrays(source, categories, duplicate, PARAMETERS, scorer, lambda: ())
        pending[0].candidate.consensus += "A"
        with self.assertRaisesRegex(ValueError, "representation or provenance"):
            rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, scorer, lambda: ())
        for attribute in ("consensus", "list_gaped_repeats", "competition_origins"):
            source, categories, pending, scorer = fixture()
            def corrupt(requests):
                rows = scorer(requests)
                flank = next(iter(requests.values()))[1]
                if attribute == "consensus":
                    flank.consensus += "A"
                elif attribute == "list_gaped_repeats":
                    flank.list_gaped_repeats = ["T" + g[1:] for g in flank.list_gaped_repeats]
                else:
                    mark_repair(flank, source)
                return rows
            with self.subTest(attribute=attribute), self.assertRaises(ValueError):
                rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, corrupt, lambda: ())
            self.assertEqual(categories, [{}, {}, {}, {}, {}])

    def test_capture_is_owned_and_block_localization_is_not_a_threshold_override(self):
        data = load_fixture()
        source = data["source_window"]["sequence"]
        source_hash = hashlib.sha256(source.encode()).hexdigest()
        candidate = local_candidate(data["cases"]["260"]["observed_feasibility_result"]["assessments"]["full"], data["source_window"]["absolute_start"])
        structural = build_advanced_fuzzy_search_filter(PARAMETERS)
        pending = capture_pending_context(source, source_hash, (0, len(source)), candidate, structural)
        candidate.source_path_classification_eligible = False
        pending.validate_current(source_hash)
        self.assertTrue(pending.candidate.source_path_classification_eligible)
        base = pending.candidate
        for modification in ("outside_duplicate", "two_runs", "uniform", "different_strings"):
            changed = copy.deepcopy(base)
            if modification == "outside_duplicate":
                changed.list_spacers[0] = changed.list_spacers[10]
            elif modification == "two_runs":
                changed.list_spacers[:4] = [changed.list_spacers[0]] * 4
            elif modification == "uniform":
                changed.list_spacers[:] = [changed.list_spacers[10]] * 16
            else:
                changed.list_spacers[:] = ["A" * 30 + format(i, "04b") for i in range(16)]
            self.assertIsNone(localized_block(changed, structural))
        self.assertIsNone(capture_pending_context(source, source_hash, (0, len(source)), candidate, structural))

    def test_higher_rank_failure_is_audited_but_equal_rank_failure_vetoes(self):
        source, output, audit = self.run_fixture(boundary_tie=True)
        records = copy.deepcopy(audit["observations"])
        records[1]["qualified"] = False
        decision, members = _choose(audit["source_sha256"], records, [])
        self.assertEqual(decision["reason"], "unsafe_same_rank_peer")
        self.assertFalse(members)
        records[1]["fresh"]["score"] += .01
        decision, members = _choose(audit["source_sha256"], records, [])
        self.assertEqual(decision["action"], "selected")
        self.assertEqual(len(members), 1)
        self.assertEqual(decision["observation_ids"], [r["observation_id"] for r in records])
        records[1]["qualified"] = True
        records[1]["fresh"]["score"] = records[0]["fresh"]["score"]
        records[1]["fresh"]["repeat_intervals"][-1][0] += 1
        decision, _ = _choose(audit["source_sha256"], records, [])
        self.assertEqual(decision["reason"], "different_start_tied_geometries")

    def test_union_envelope_and_direct_all_member_witness_are_checked(self):
        source, output, audit = self.run_fixture(boundary_tie=True)
        records = audit["observations"]
        ends = [r["fresh"]["repeat_intervals"][-1][1] for r in records]
        barrier = {"category_index": 0, "locus_key": [0, 1], "interval": [min(ends), max(ends)], "representation_sha256": "synthetic"}
        decision, _ = _choose(audit["source_sha256"], records, [barrier])
        self.assertEqual(decision["reason"], "accepted_envelope_collision")
        altered = copy.deepcopy(records)
        from components.occurrence_protection import occurrence_snapshot
        witness = copy.deepcopy(flatten(output, 4)[0])
        # An exact source-valid terminal occurrence outside both admitted paths.
        witness[1].list_repeat_starts = [len(source)-1]
        witness[1].list_repeats = [source[-1]]
        witness[1].list_repeats_gaped = [source[-1]]
        witness[1].list_gaped_repeats = [source[-1]]
        witness[1].list_spacers = []
        altered[0]["ordinary_references"] = [occurrence_snapshot(witness)]
        decision, _ = _choose(audit["source_sha256"], altered, [])
        self.assertEqual(decision["reason"], "original_occurrence_conflict")
        self.assertEqual(decision["ordinary_witness_guard"]["direct_checks"], 2)

    def test_current_repair_owner_fence_and_historical_scope_never_clear(self):
        source, categories, pending, scorer = fixture()
        candidate = copy.deepcopy(pending[0].candidate)
        mark_repair(candidate, source)
        pending = [capture_pending_context(source, hashlib.sha256(source.encode()).hexdigest(), (0, len(source)), candidate,
                    build_advanced_fuzzy_search_filter(PARAMETERS))]
        full = scorer({"f": [candidate]})["f"][0]
        # Construct ordinary flanks independently, not inherited repair origins.
        data = load_fixture()
        result = data["cases"]["260"]["observed_feasibility_result"]
        flanks = [local_candidate(result["assessments"][role], data["source_window"]["absolute_start"]) for role in ("left", "right")]
        owners = scorer({"f": flanks})["f"]
        categories[1][(0, len(source))] = owners
        out, audit = rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selections_committed"], 0)
        self.assertEqual(audit["observations"][0]["current_scope"]["evidence"]["decision"], "repair_scope_conflict")
        replay(source, out, audit)
        _groups, conflicts = initial_competition_groups(owners + [full])
        categories = [{}, {}, {}, {(0, len(source)): list(conflicts)}, {}]
        out, audit = rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selections_committed"], 0)
        self.assertFalse(audit["observations"][0]["current_scope"]["blocked"])
        self.assertEqual(audit["observations"][0]["initial_competition"][0]["decision"], "repair_scope_conflict")
        replay(source, out, audit)

    def test_single_owner_reference_is_bound_to_retained_ordinary_row(self):
        source, categories, pending, scorer = fixture()
        candidate = copy.deepcopy(pending[0].candidate)
        mark_repair(candidate, source)
        pending = [capture_pending_context(source, hashlib.sha256(source.encode()).hexdigest(), (0, len(source)), candidate,
                    build_advanced_fuzzy_search_filter(PARAMETERS))]
        data = load_fixture()
        result = data["cases"]["260"]["observed_feasibility_result"]
        owner = local_candidate(result["assessments"]["left"], data["source_window"]["absolute_start"])
        categories[1][(0, len(source))] = scorer({"x": [owner]})["x"]
        output, audit = rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selections_committed"], 1)
        self.assertEqual(len(audit["observations"][0]["ordinary_references"]), 1)
        replay(source, output, audit)
        missing = list(output)
        missing[1] = {}
        with self.assertRaises(ValueError):
            replay(source, missing, audit)

    def test_complete_nonselected_replay_and_typed_mutations_fail(self):
        source, output, audit = self.run_fixture(ordinals=[258, 260, 261])
        for mutation in ("count", "ordinal", "nonselected_score", "nonselected_validation", "member_omission", "source", "retained"):
            out, evidence = copy.deepcopy(output), copy.deepcopy(audit)
            if mutation == "count":
                evidence["stats"]["selections_committed"] = True
            elif mutation == "ordinal":
                evidence["observations"][0]["ordinal"] = False
            elif mutation == "nonselected_score":
                evidence["observations"][-1]["fresh"]["score"] += .001
            elif mutation == "nonselected_validation":
                evidence["observations"][-1]["flanks"]["right"]["validation"]["qualified"] = False
            elif mutation == "member_omission":
                row = flatten(out, 2)[0]
                basis = row[3]["classification_basis"]
                basis["contextual_spacer_context"]["members"].pop()
                basis["evidence_sha256"] = digest(basis["contextual_spacer_context"])
            elif mutation == "source":
                evidence["source_sha256"] = "a" * 64
            else:
                next(iter(out[4].values())).pop()
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                replay(source, out, evidence)

    def test_existing_repair_bridge_connects_two_frozen_ordinary_owners(self):
        from tests.test_initial_competition import row
        left = row(.7, [(0, 20), (40, 60), (80, 100)])
        right = row(.7, [(250, 270), (290, 310), (330, 350)])
        proposal = row(.8, [(0, 20), (40, 60), (80, 100), (140, 160)], repair=True)
        bridge = row(.65, [(150, 170), (200, 220), (250, 270)], repair=True)
        self.assertFalse(_scope_checks([left, right], [proposal])[id(proposal[1])]["blocked"])
        decision = _scope_checks([left, right], [proposal], [bridge])[id(proposal[1])]
        self.assertTrue(decision["blocked"])
        self.assertEqual(decision["evidence"]["decision"], "repair_scope_conflict")
        self.assertEqual(len(decision["evidence"]["ordinary_components"]), 2)

    def test_exact_equal_rank_negative_peer_and_recorded_metrics_replay(self):
        source, categories, pending, scorer = fixture(boundary_tie=True)
        target = pending[1].candidate
        def tagged(requests):
            result = scorer(requests)
            for rows in result.values():
                for row in rows:
                    if row[1] is target:
                        row.append({"tags": ["tandem"]})
            return result
        output, audit = rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, tagged, lambda: ())
        self.assertEqual(audit["decisions"][0]["reason"], "unsafe_same_rank_peer")
        self.assertEqual(audit["stats"]["selections_committed"], 0)
        replay(source, output, audit)
        for record in audit["observations"]:
            evidence = validate_contextual_recorded_metrics(source, record, audit["structural_parameters"])
            self.assertEqual(evidence["source_and_stored_metrics"], "passed")
        corrupted = copy.deepcopy(audit["observations"][0])
        corrupted["flanks"]["right"]["fresh"]["features"] = copy.deepcopy(corrupted["fresh"]["features"])
        with self.assertRaises(ValueError):
            validate_contextual_recorded_metrics(source, corrupted, audit["structural_parameters"])

    def test_whole_audit_cannot_omit_a_current_ordinary_owner(self):
        source, categories, pending, scorer = fixture()
        candidate = copy.deepcopy(pending[0].candidate)
        mark_repair(candidate, source)
        pending = [capture_pending_context(source, hashlib.sha256(source.encode()).hexdigest(), (0, len(source)), candidate,
                    build_advanced_fuzzy_search_filter(PARAMETERS))]
        data = load_fixture()
        result = data["cases"]["260"]["observed_feasibility_result"]
        owners = scorer({"x": [local_candidate(result["assessments"][role], data["source_window"]["absolute_start"])
                                for role in ("left", "right")]})["x"]
        categories[1][(0, len(source))] = owners
        output, audit = rescue_contextual_spacer_arrays(source, categories, pending, PARAMETERS, scorer, lambda: ())
        tampered = copy.deepcopy(audit)
        tampered["ordinary_scope_rows"].pop()
        record = tampered["observations"][0]
        full = flatten(output, 4)[0]
        scope = _scope_checks(owners[:1], [full])[id(full[1])]
        record["current_scope"] = scope
        record["ordinary_references"] = _references(record["initial_competition"], scope)
        record["qualified"] = _qualified(record)
        report = {"source": {"sha256": audit["source_sha256"]}, "arrays": []}
        with self.assertRaisesRegex(ValueError, "inventory is incomplete"):
            validate_contextual_spacer_audit(report, tampered, provenance_categories=output, source_sequence=source)

    def test_projected_quarantine_requires_exact_prequarantine_disposition(self):
        source, output, audit = self.run_fixture()
        row = flatten(output, 2)[0]
        # This unit tests only the audit join contract. OutputMaker validates
        # the biological demotion record independently before supplying it.
        key = next(iter(output[2]))
        disposition = {"action": "quarantined", "source_category_index": 2,
            "locus_key": list(key), "row_ordinal": 0, "score": row[0],
            "source_architecture": {"repeats": [{"source_interval": {"start": a, "end": b}}
                                                  for a, b in _geometry(row[1])]}}
        report = {"source": {"sha256": audit["source_sha256"]}, "arrays": []}
        validate_contextual_spacer_audit(report, audit, provenance_categories=output,
                                        likely_tandem_records=[disposition], source_sequence=source)
        disposition["row_ordinal"] = True
        with self.assertRaises(ValueError):
            validate_contextual_spacer_audit(report, audit, provenance_categories=output,
                                            likely_tandem_records=[disposition], source_sequence=source)

    def test_refinement_capture_is_side_channel_and_keeps_default_fallback_exact(self):
        from components.module_detection_refinement import DetectionRefinement
        from tests.test_repeat_boundary_repair import fuzzy, Locus
        from components.localized_spacer_context import candidate_capture
        data = load_fixture()
        source = data["source_window"]["sequence"]
        source_hash = hashlib.sha256(source.encode()).hexdigest()
        paths = [local_candidate(data["cases"][str(i)]["observed_feasibility_result"]["assessments"]["full"],
                                 data["source_window"]["absolute_start"]) for i in (258, 260)]
        outcomes = []
        for enabled in (False, True):
            mapping = {Locus(0, len(source)): [fuzzy(p) for p in paths]}
            with patch("components.module_detection_refinement.propose_terminal_repairs", return_value=([], [], {})), \
                    patch("components.module_detection_refinement.propose_internal_repairs", return_value=([], [])):
                refined = DetectionRefinement(mapping, PARAMETERS, False, source_dna=source,
                                              capture_contextual_spacers=enabled)
            outcomes.append([[candidate_capture(candidate, source_hash) for candidate in values]
                             for values in refined.output().values()])
            self.assertEqual(len(refined.pending_contextual_spacers), 2 if enabled else 0)
            for pending in refined.pending_contextual_spacers:
                pending.validate_current(source_hash)
                self.assertEqual(len(pending.candidate.list_repeats), 17)
        self.assertEqual(outcomes[0], outcomes[1])
        with self.assertRaises(TypeError):
            DetectionRefinement({}, PARAMETERS, False, capture_contextual_spacers=1)

    def test_empty_route_has_no_scoring_cds_or_owner_inventory_work(self):
        source, categories, _pending, _scorer = fixture()
        def forbidden(*args, **kwargs):
            raise AssertionError("empty contextual lane performed work")
        with patch("components.contextual_spacer_rescue.validate_ordinary_witness_inventory", side_effect=forbidden):
            output, audit = rescue_contextual_spacer_arrays(source, categories, [], PARAMETERS, forbidden, forbidden)
        self.assertIs(output, categories)
        replay(source, output, audit)


if __name__ == "__main__":
    unittest.main()
