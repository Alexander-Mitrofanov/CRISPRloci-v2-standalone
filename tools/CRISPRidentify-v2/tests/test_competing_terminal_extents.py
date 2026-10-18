"""Source-real, observational terminal disclosure; no model or detector runs."""
import copy
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from components.candidate_evidence import stage_candidate_evidence, _identity
from components.candidate_origins import inherit_candidate_origins, mark_repair_origin
from components.components_detection_refinement import CrisprCandidate
from components.competing_terminal_extents import (
    KEY, build_competing_terminal_extents, competing_terminal_extent_summaries,
    exact_terminal_relation, validate_competing_terminal_extents,
)
from components.initial_competition import initial_competition_groups
from components.likely_tandem_transaction import apply_likely_tandem_quarantine
from components.output_transaction import OutputTransaction
from components.report_extensions import attach_report_extensions
from components.report_model import ReportModelBuilder
from tests.test_candidate_evidence import PARAMETERS as BASE_PARAMETERS


PARAMETERS = dict(BASE_PARAMETERS, param_reconstruction_policy="conservative")


def fixture():
    """Actual local strings/vectors, explicitly synthetic local ordinary origins.

    The full absolute historical row/origin attestations stay in the fixture.
    They are not transplanted onto the translated 550 bp source window.
    No prediction or biological-negative clearance is claimed.
    """
    data = json.loads((Path(__file__).parent / "fixtures/competing_terminal_cp113104.json").read_text())
    source = data["source_window"]["sequence"]
    categories = [{}, {}, {}, {}, {}]
    names = None
    for index, view in zip((2, 3), data["synthetic_local_construction_views"]):
        candidate = CrisprCandidate(list_repeat_starts=view["repeat_starts"],
            list_repeats=view["repeats"], list_repeats_gaped=view["gapped_repeats"],
            list_spacers=view["spacers"])
        sets = view["stored_feature_sets"]
        names = [[f["name"] for f in group["features"]] for group in sets]
        vectors = [np.asarray([[f["value"] for f in group["features"]]], dtype=float) for group in sets]
        categories[index][tuple(view["local_search_key"])] = [[view["stored_own_score"], candidate, vectors]]
    return source, categories, names


def report_fixture(source, categories, names, *, parameters=None, provenance=None):
    parameters = PARAMETERS if parameters is None else parameters
    with tempfile.TemporaryDirectory() as temporary:
        fasta = Path(temporary) / "fixture.fa"
        fasta.write_text(">fixture\n" + source + "\n")
        return ReportModelBuilder(str(fasta), categories, {}, names, parameters=parameters,
                                  provenance_categories=provenance).build()


def context(source, categories, names, *, provenance=None, parameters=None):
    return dict(source_sequence=source, accepted_categories=categories,
                provenance_categories=categories if provenance is None else provenance,
                parameters=PARAMETERS if parameters is None else parameters, feature_names=names)


def build_fixture():
    source, categories, names = fixture()
    report = report_fixture(source, categories, names)
    inputs = context(source, categories, names)
    extension = build_competing_terminal_extents(report, **inputs)
    report["extensions"] = {KEY: extension}
    return source, categories, names, report, inputs


def _row(categories, category):
    return next(iter(categories[category].values()))[0]


def source_candidate(source, intervals):
    repeats = [source[a:b] for a, b in intervals]
    width = max(map(len, repeats))
    return CrisprCandidate(list_repeat_starts=[a for a, _ in intervals],
        list_repeats=repeats, list_repeats_gaped=[r + "-" * (width - len(r)) for r in repeats],
        list_spacers=[source[a[1]:b[0]] for a, b in zip(intervals, intervals[1:])])


class CompetingTerminalTests(unittest.TestCase):
    def test_real_five_six_observation_and_diagnostics_do_not_change_science(self):
        source, categories, names = fixture()
        original = copy.deepcopy(categories)
        report = report_fixture(source, categories, names)
        before = copy.deepcopy(report)
        extension = build_competing_terminal_extents(report, **context(source, categories, names))
        validate_competing_terminal_extents(report, extension, **context(source, categories, names))
        self.assertEqual(report, before)
        self.assertEqual(extension["counts"], dict(accepted_anchors=1, competing_geometries=1,
            candidate_observations=1, links=1, additional_accepted_arrays=0, additional_accepted_repeats=0))
        observation = next(iter(extension["observations"].values()))
        self.assertEqual(observation["score"], .5540780107828767)
        self.assertEqual(observation["repeat_count"], 6)
        self.assertEqual(observation["source_span"], [115, 448])
        self.assertEqual(observation["legacy_locus_key"], [15, 420])
        self.assertEqual(observation["negative_evidence"]["accepted_family_tandem_review_status"], "not_recorded")
        link = extension["links"][0]
        self.assertEqual(link["relation"], dict(anchor_to_candidate_indexes=list(range(5)),
            added_left_indexes=[], added_right_indexes=[5]))
        diagnostic = link["observations"][0]
        self.assertEqual(diagnostic["status"], "unaccepted_observation")
        self.assertEqual(diagnostic["original_occurrences"], {"status": "no_recorded_witness", "checks": []})
        self.assertIn("stored_score_below_current_donor_floor", diagnostic["routing_exclusions"])
        self.assertIsNone(diagnostic["identity_diagnostic"]["admission_cutoff"])
        for i in (2, 3):
            self.assertEqual(_row(categories, i)[1].__dict__, _row(original, i)[1].__dict__)
            np.testing.assert_array_equal(_row(categories, i)[2][0], _row(original, i)[2][0])

    def test_exact_terminal_relation_is_not_internal_phase_or_transitive(self):
        a = [[10, 20], [30, 40]]
        for b, mapping in (([[0, 5]] + a, [1, 2]), (a + [[50, 60]], [0, 1]),
                           ([[0, 5]] + a + [[50, 60]], [1, 2])):
            self.assertEqual(exact_terminal_relation(a, b)["anchor_to_candidate_indexes"], mapping)
        for b in (a, [[9, 20], [30, 40], [50, 60]], [[10, 20], [22, 25], [30, 40]],
                  [[0, 5], [10, 20], [31, 40]], [[0, 5], [30, 40]]):
            self.assertIsNone(exact_terminal_relation(a, b))
        with self.assertRaises(ValueError):
            exact_terminal_relation(a, [[False, 5], [10, 20], [30, 40]])

    def test_geometry_aliases_and_repeated_object_locations_preserve_full_metadata(self):
        source, categories, names = fixture()
        key = next(iter(categories[3]))
        row = categories[3][key][0]
        # Additional metadata must be preserved even when it is not part of the
        # existing sidecar candidate ID. These records are not authority votes.
        row.append({"classification_basis": {"method": "model_score_band"}, "audit_note": "synthetic original"})
        alias = copy.deepcopy(row)
        candidate = alias[1]
        aligned = [r[:i + 1] + " " + r[i + 1:] for i, r in enumerate(candidate.list_repeats_gaped)]
        alias[1] = CrisprCandidate(candidate.list_repeats, aligned, candidate.list_spacers, candidate.list_repeat_starts)
        inherit_candidate_origins(alias[1], [candidate])
        alias[3]["audit_note"] = "synthetic alignment alias"
        lower = copy.deepcopy(row)
        lower[0] = .51
        categories[3][key] = [row, alias, row, lower]
        report = report_fixture(source, categories, names)
        extension = build_competing_terminal_extents(report, **context(source, categories, names))
        self.assertEqual(extension["counts"]["candidate_observations"], 4)
        self.assertEqual(extension["counts"]["competing_geometries"], 1)
        self.assertEqual({r["row_metadata"]["audit_note"] for r in extension["observations"].values()},
                         {"synthetic original", "synthetic alignment alias"})
        self.assertEqual(len({r["observation_id"] for r in extension["observations"].values()}), 4)
        self.assertEqual(len({r["candidate_id"] for r in extension["observations"].values()}), 3)
        categories[3][key].reverse()
        permuted = build_competing_terminal_extents(report, **context(source, categories, names))
        self.assertEqual(extension["counts"], permuted["counts"])
        self.assertEqual(extension["links"][0]["source_geometry_id"], permuted["links"][0]["source_geometry_id"])

    def test_scope_eligibility_and_structure_are_not_broadened(self):
        for change in ("cross_key", "ineligible", "low", "failed_structure"):
            source, categories, names = fixture()
            key = next(iter(categories[3]))
            row = categories[3][key][0]
            if change == "cross_key":
                categories[3] = {(0, len(source)): [row]}
            elif change == "ineligible":
                row[1].source_path_classification_eligible = False
            elif change == "low":
                categories[4] = categories[3]
                categories[3] = {}
            params = dict(PARAMETERS, **({"param_min_repeats": 7} if change == "failed_structure" else {}))
            report = report_fixture(source, categories, names, parameters=params)
            extension = build_competing_terminal_extents(report, **context(source, categories, names, parameters=params))
            self.assertEqual(extension["observations"], {}, change)

    def test_known_negative_is_blocked_but_missing_is_not_passed(self):
        source, categories, names = fixture()
        row = _row(categories, 3)
        row.append({"classification_basis": {"method": "model_score_band"}, "tags": ["tandem"]})
        report = report_fixture(source, categories, names)
        extension = build_competing_terminal_extents(report, **context(source, categories, names))
        diagnostic = extension["links"][0]["observations"][0]
        self.assertEqual(diagnostic["status"], "blocked_observation")
        self.assertEqual(diagnostic["block_reasons"], ["candidate_recorded_negative"])
        self.assertEqual(next(iter(extension["observations"].values()))["recorded_negative_reasons"], ["tandem_tag"])

    def test_collision_checks_both_accepted_categories(self):
        for index, score in ((0, .8), (2, .61)):
            source, categories, names = fixture()
            # Source-real slices sufficient for publication, not claimed as a
            # biological true array. Their span overlaps the added tail only.
            candidate = source_candidate(source, [(425, 448), (460, 483), (495, 518)])
            categories[index][(425, 518)] = [[score, candidate, copy.deepcopy(_row(categories, 2)[2])]]
            report = report_fixture(source, categories, names)
            extension = build_competing_terminal_extents(report, **context(source, categories, names))
            diagnostic = extension["links"][0]["observations"][0]
            self.assertIn("overlaps_other_accepted_array", diagnostic["block_reasons"])
            self.assertEqual(len(diagnostic["collisions"]), 1)

    def test_original_s1_witness_and_multi_owner_conflicts_are_blocked(self):
        from tests.test_candidate_evidence import _fixture
        source, sample = _fixture()
        all_starts = list(range(0, len(source), 57))
        intervals = [(a, a + 37) for a in all_starts]
        key = (0, len(source))
        for mode in ("occurrence", "scope"):
            anchor = [.61, source_candidate(source, intervals[:3]), []]
            if mode == "occurrence":
                ordinary = [.59, source_candidate(source, intervals[:6]), []]
                proposed = [.55, source_candidate(source, [intervals[i] for i in (0, 1, 2, 4, 5, 6)]), []]
                ordinary_rows = [ordinary]
            else:
                ordinary_rows = [anchor, [.6, source_candidate(source, intervals[5:8]), []]]
                proposed = [.55, source_candidate(source, intervals), []]
            candidate = proposed[1]
            mark_repair_origin(candidate, {"policy": "terminal-only-trim-v1", "status": "generated",
                "coordinate_system": "zero_based_half_open", "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
                "child_geometry": [candidate.list_repeat_starts, [a + len(r) for a, r in zip(
                    candidate.list_repeat_starts, candidate.list_repeats)]]})
            _groups, conflicts = initial_competition_groups(ordinary_rows + [proposed])
            self.assertEqual(len(conflicts), 1)
            # Staged source-valid historical witness observations, not a claim
            # that this synthetic final category placement reruns the ranking.
            categories = [{}, {key: ordinary_rows}, {key: [anchor]}, {key: list(conflicts)}, {}]
            report = report_fixture(source, categories, [])
            extension = build_competing_terminal_extents(report, **context(source, categories, []))
            diagnostic = extension["links"][0]["observations"][0]
            self.assertEqual(diagnostic["status"], "blocked_observation")
            expected = "original_occurrence_conflict" if mode == "occurrence" else "candidate_repair_scope_conflict"
            self.assertIn(expected, diagnostic["block_reasons"])
            validate_competing_terminal_extents(report, extension, **context(source, categories, []))

    def test_collision_visible_only_in_actual_coupled_family_union(self):
        from tests.test_reconstruction_transaction import fixture as family_fixture
        from components.reconstruction_transaction import reconstruct_accepted_arrays
        from components.reconstruction_boundary_family import boundary_family_summary
        from components.feature_schema import COMPACT_FEATURE_SCHEMA
        source, categories, scorer = family_fixture()
        key = next(iter(categories[1]))
        donor = categories[1][key][0][1]
        repeats = donor.list_repeats[:]
        repeats[-1] = repeats[-1][:-6]
        short = CrisprCandidate(repeats, [r + ("-" * 6 if i == 5 else "") for i, r in enumerate(repeats)],
                                donor.list_spacers[:], donor.list_repeat_starts[:])
        categories[1][key].append(scorer({key: [short]})[key][0])
        # The validated family's short actual display ends at431; its other
        # actual member ends437. This terminal hypothesis begins at431.
        motif = source[431:460]
        source += "ACGTACG"
        other_intervals = []
        for spacer in ("AC" * 16, "GT" * 16, "AG" * 16):
            start = len(source)
            source += motif
            other_intervals.append((start, len(source)))
            source += spacer
        categories[2][(431, len(source))] = [[.61, source_candidate(source, other_intervals),
                                             copy.deepcopy(categories[0][key][0][2])]]
        categories[3][(431, len(source))] = [[.55, source_candidate(source, [(431, 460)] + other_intervals),
                                             copy.deepcopy(categories[1][key][0][2])]]
        output, _audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        accepted_family = output[0][key][0]
        family = boundary_family_summary(accepted_family)
        self.assertIsNotNone(family)
        self.assertEqual(family["union_span"][1], 437)
        self.assertEqual(accepted_family[1].list_repeat_starts[-1] + len(accepted_family[1].list_repeats[-1]), 431)
        names = [COMPACT_FEATURE_SCHEMA]
        report = report_fixture(source, output, names)
        extension = build_competing_terminal_extents(report, **context(source, output, names))
        diagnostic = extension["links"][0]["observations"][0]
        self.assertEqual(diagnostic["collisions"][0]["source_span"], [32, 437])
        self.assertEqual(diagnostic["status"], "blocked_observation")

    def test_current_anchor_locations_match_post_quarantine_sidecar(self):
        source, categories, names = fixture()
        repeat = "GTTTCAGTATCCTCGCTGTTGGATCGTTGCCTAGACC"
        prefix = (repeat + "A" * 20) * 3
        shift = len(prefix)
        moved = [{}, {}, {}, {}, {}]
        for i in (2, 3):
            key = next(iter(categories[i]))
            row = categories[i][key][0]
            geometry = [(a + shift, a + shift + len(r)) for a, r in zip(
                row[1].list_repeat_starts, row[1].list_repeats)]
            moved[i][tuple(x + shift for x in key)] = [[row[0], source_candidate(prefix + source, geometry), row[2]]]
        early = source_candidate(prefix + source, [(i * 57, i * 57 + 37) for i in range(3)])
        moved[2][(0, shift)] = [[.62, early, copy.deepcopy(_row(categories, 2)[2])]]
        source = prefix + source
        final, _, stats = apply_likely_tandem_quarantine(moved)
        self.assertEqual(sum(len(v) for v in final[2].values()), 1)
        report = report_fixture(source, final, names, provenance=moved)
        extension = build_competing_terminal_extents(report, **context(source, final, names, provenance=moved))
        anchor = next(iter(extension["anchors"].values()))["observation"]
        self.assertEqual(anchor["category_position"], 0)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            with OutputTransaction() as publication:
                stage_candidate_evidence(publication, path, source_header="fixture", source_sequence=source,
                    categories=final, parameters=PARAMETERS, feature_names=names, canonical_report=report)
                publication.publish()
            retained = next(json.loads(line) for line in gzip.decompress(path.read_bytes()).splitlines()
                            if json.loads(line).get("category") == "Possible")
        self.assertEqual(anchor["observation_id"], retained["observation_id"])
        self.assertEqual(retained, {k: anchor[k] for k in retained})

    def test_rehashed_advertised_negative_and_unknown_rule_fail_recorded_replay(self):
        _source, _categories, _names, report, _inputs = build_fixture()
        extension = report["extensions"][KEY]
        oid = next(iter(extension["observations"]))
        for mutation in (lambda r: r["negative_evidence"].update(accepted_family_tandem_review_status="passed"),
                         lambda r: r["negative_evidence"].update(tags=["tandem"]),
                         lambda r: r.update(classification_basis_method="invented"),
                         lambda r: r["feature_values"].update(repeat_len=100),
                         lambda r: r.update(score_origin="freshly_scored"),
                         lambda r: r["row_metadata"].update(initial_competition={"rule": "invented"})):
            altered = copy.deepcopy(extension)
            record = altered["observations"][oid]
            mutation(record)
            record["snapshot_sha256"] = _identity({k: v for k, v in record.items() if k != "snapshot_sha256"})
            with self.assertRaises(ValueError):
                validate_competing_terminal_extents(report, altered)

    def test_source_score_origin_and_typed_relationship_tampering_fail(self):
        source, categories, names, report, inputs = build_fixture()
        extension = report["extensions"][KEY]
        oid = next(iter(extension["observations"]))
        mutations = [
            lambda e: e["observations"][oid].update(score=True),
            lambda e: e["observations"][oid]["repeat_intervals"][0].__setitem__(0, 114),
            lambda e: e["observations"][oid]["candidate_origins"].update(ordinary_birth_sha256="a" * 64),
            lambda e: e["observations"][oid]["feature_sets"][0]["features"][0].update(value=100.0),
            lambda e: e["links"][0]["relation"]["anchor_to_candidate_indexes"].__setitem__(0, False),
            lambda e: e.update(coverage="unknown"),
            lambda e: e["observations"].clear(),
        ]
        for mutate in mutations:
            altered = copy.deepcopy(extension)
            mutate(altered)
            with self.assertRaises(ValueError):
                validate_competing_terminal_extents(report, altered, **inputs)
        # A coherent snapshot-digest rewrite still cannot defeat actual-row
        # membership, even for an otherwise plausible finite own score.
        altered = copy.deepcopy(extension)
        r = altered["observations"][oid]
        r["score"] = .56
        r["snapshot_sha256"] = _identity({k: v for k, v in r.items() if k != "snapshot_sha256"})
        with self.assertRaises(ValueError):
            validate_competing_terminal_extents(report, altered, **inputs)
        bad_source = dict(inputs, source_sequence="A" + source[1:])
        with self.assertRaises(ValueError):
            build_competing_terminal_extents(report, **bad_source)
        _row(categories, 3)[1].list_repeats_gaped[0] = "A" + _row(categories, 3)[1].list_repeats_gaped[0][1:]
        with self.assertRaises(ValueError):
            build_competing_terminal_extents(report, **inputs)

    def test_source_snapshots_equal_optional_sidecar_and_projection_does_not_hide_them(self):
        source, categories, names, report, inputs = build_fixture()
        extension = report["extensions"][KEY]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "evidence.gz"
            with OutputTransaction() as transaction:
                stage_candidate_evidence(transaction, path, source_header="fixture", source_sequence=source,
                    categories=categories, parameters=PARAMETERS, feature_names=names, canonical_report=report)
                transaction.publish()
            records = [json.loads(x) for x in gzip.decompress(path.read_bytes()).splitlines()]
        retained = next(r for r in records if r.get("category") == "Possible Discarded")
        observed = extension["observations"][retained["observation_id"]]
        self.assertEqual(retained, {k: observed[k] for k in retained})
        projected = [dict(categories[0]), {}, dict(categories[2]), {}, {}]
        after = build_competing_terminal_extents(report, **context(source, projected, names, provenance=categories))
        self.assertEqual(extension, after)
        quarantine, _, _ = apply_likely_tandem_quarantine(categories)
        self.assertIs(quarantine[3][next(iter(categories[3]))][0], _row(categories, 3))

    def test_extension_only_api_is_pure_and_legacy_or_missing_inventory_omits(self):
        source, categories, names = fixture()
        report = report_fixture(source, categories, names)
        before = copy.deepcopy(report)
        options = context(source, categories, names)
        context_without_source = {k: v for k, v in options.items() if k != "source_sequence"}
        attached = attach_report_extensions(report, source_sequence=source,
            competing_terminal_extent_context=context_without_source)
        self.assertEqual(report, before)
        self.assertEqual({k: v for k, v in attached.items() if k != "extensions"}, report)
        summary = next(iter(competing_terminal_extent_summaries(attached).values()))[0]
        self.assertEqual((summary["accepted_repeat_count"], summary["repeat_count"]), (5, 6))
        self.assertEqual(summary["status"], "unaccepted_observation")
        self.assertEqual(competing_terminal_extent_summaries(report), {})
        for changed in (dict(options, parameters=dict(PARAMETERS, param_reconstruction_policy="legacy")),
                        dict(options, provenance_categories=None)):
            self.assertIsNone(build_competing_terminal_extents(report, **changed))


if __name__ == "__main__":
    unittest.main()
