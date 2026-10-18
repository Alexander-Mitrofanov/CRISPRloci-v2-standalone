"""Atomic category updates and source/model-bound reporting for reconstruction."""

import copy
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest

import numpy as np

from components.compact_model import predict_compact_scores
from components.components_detection_refinement import CrisprCandidate
from components.feature_schema import COMPACT_FEATURE_SCHEMA
from components.initial_competition import (
    initial_competition_groups,
    with_competition_evidence,
)
from components.likely_tandem_transaction import apply_likely_tandem_quarantine
from components.reconstruction_transaction import (
    METHOD, reconstruct_accepted_arrays, validate_reconstruction_basis,
)
from components.report_model import ReportModelBuilder, ReportModelError
from tests.test_candidate_evidence import PARAMETERS


def fixture(*, equal_period=False):
    motif = "CGGTTTCAGTATCCTCGCTGTTGGATCGTTGCCTAGACCC"
    rng = random.Random(5077)
    source = "ACGT" * 8
    intervals = []
    for index in range(6):
        start = len(source)
        source += motif
        intervals.append((start, len(source)))
        spacer_length = 35 if equal_period else 31 + index
        source += "".join(rng.choice("ACGT") for _ in range(spacer_length))

    def candidate(count):
        selected = intervals[:count]
        return CrisprCandidate(
            list_repeats=[source[a:b] for a, b in selected],
            list_repeats_gaped=[source[a:b] for a, b in selected],
            list_spacers=[source[a[1]:b[0]] for a, b in zip(selected, selected[1:])],
            list_repeat_starts=[a for a, _ in selected],
        )

    def scored(candidate):
        # Model-real synthetic vectors test transaction binding. They are not
        # presented as independent feature extraction or biological evidence.
        n = len(candidate.list_repeats)
        features = np.array([[40, n, 1 if n == 3 else .97, .45, 34.5,
                              .254305, 0, 0, -8.2, 15.9]])
        score = float(predict_compact_scores(features)[0])
        return [score, candidate, [features]]

    categories = [{(0, len(source)): [scored(candidate(3))]},
                  {(0, len(source)): [scored(candidate(6))]}, {}, {}, {}]
    def score_batch(requests):
        return {key: [scored(c) for c in values] for key, values in requests.items()}
    return source, categories, score_batch


class ReconstructionTransactionTests(unittest.TestCase):
    def _initial_tie_fixture(self):
        source, categories, scorer = fixture()
        key = next(iter(categories[0]))
        core = categories[0][key][0]
        donor = categories[1][key][0]
        starts = [donor[1].list_repeat_starts[index] for index in (0, 1, 3)]
        length = len(core[1].list_repeats[0])
        tied_candidate = CrisprCandidate(
            list_repeat_starts=starts,
            list_repeats=[source[start:start + length] for start in starts],
            list_repeats_gaped=[source[start:start + length] for start in starts],
            list_spacers=[source[left + length:right] for left, right in zip(starts, starts[1:])],
        )
        tied = [core[0], tied_candidate, copy.deepcopy(core[2])]
        observation = initial_competition_groups([core, tied])[0][0].evidence
        categories[0][key][0] = with_competition_evidence(core, observation)
        categories[1][key].append(tied)
        return source, categories, scorer, observation

    def test_initial_tie_is_historical_evidence_not_reconstruction_authority(self):
        # Two tied short parses are both compatible with one longer donor.
        # Their initial ambiguity must neither forbid this selection nor be
        # silently relabelled as a fresh score/geometry observation.
        source, categories, scorer, observation = self._initial_tie_fixture()
        original_core = next(iter(categories[0].values()))[0]
        fresh_vector = None
        for index, values in ((2, (.99, .95, .9)), (8, (-20., -1.)), (3, (.3, .6))):
            for value in values:
                vector = original_core[2][0].copy()
                vector[0, index] = value
                score = float(predict_compact_scores(vector)[0])
                if score >= .75 and abs(score - original_core[0]) > 1e-12:
                    fresh_vector = vector
                    break
            if fresh_vector is not None:
                break
        self.assertIsNotNone(fresh_vector)

        def fresh_scorer(requests):
            result = scorer(requests)
            for rows in result.values():
                for row in rows:
                    if row[1] is original_core[1]:
                        row[2] = [fresh_vector.copy()]
                        row[0] = float(predict_compact_scores(fresh_vector)[0])
            return result

        output, audit = reconstruct_accepted_arrays(
            source, categories, PARAMETERS, fresh_scorer, lambda: ())
        selected = next(iter(output[0].values()))[0]
        basis = selected[3]["classification_basis"]
        proof = basis["reconstruction"]
        self.assertEqual(audit["stats"]["selections_committed"], 1)
        self.assertEqual(proof["anchor_initial_competition"], observation)
        self.assertEqual(proof["anchor_initial_competition"]["candidate_score"], original_core[0])
        self.assertNotEqual(proof["anchor_fresh"]["score"], original_core[0])
        self.assertNotIn("initial_competition", selected[3])
        round_trip = json.loads(json.dumps(basis))
        self.assertEqual(round_trip, validate_reconstruction_basis(
            source, selected[1], round_trip, selected[0], 0, selected, output))
        projected = [output[0], {}, output[2], {}, {}]
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "initial-tie.fa"
            fasta.write_text(">initial-tie\n" + source + "\n")
            report = ReportModelBuilder(
                fasta, projected, {}, [COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS, param_accepted_only=True),
                provenance_categories=output).build()
        published = report["arrays"][0]["classification_basis"]["reconstruction"]
        self.assertEqual(published["anchor_initial_competition"], observation)

    def test_initial_tie_metadata_cannot_be_stale_or_removed_with_a_new_digest(self):
        source, categories, scorer, _observation = self._initial_tie_fixture()
        output, _audit = reconstruct_accepted_arrays(
            source, categories, PARAMETERS, scorer, lambda: ())
        selected = next(iter(output[0].values()))[0]
        basis = selected[3]["classification_basis"]
        for mutation in ("score", "geometry", "remove", "malformed"):
            corrupt = copy.deepcopy(basis)
            observation = corrupt["reconstruction"]["anchor_initial_competition"]
            if mutation == "score":
                observation["candidate_score"] = .51
            elif mutation == "geometry":
                observation["candidate_geometry"][0][0] += 1
            elif mutation == "remove":
                del corrupt["reconstruction"]["anchor_initial_competition"]
            else:
                observation["tied_candidates"] = []
            corrupt["evidence_sha256"] = hashlib.sha256(json.dumps(
                corrupt["reconstruction"], sort_keys=True, separators=(",", ":"),
                allow_nan=False).encode()).hexdigest()
            with self.subTest(mutation=mutation), self.assertRaisesRegex(ValueError, "competition"):
                validate_reconstruction_basis(
                    source, selected[1], corrupt, selected[0], 0, selected, output)
        malformed = copy.deepcopy(categories)
        core = next(iter(malformed[0].values()))[0]
        core[3]["initial_competition"]["candidate_score"] = .51
        with self.assertRaisesRegex(ValueError, "competition"):
            reconstruct_accepted_arrays(source, malformed, PARAMETERS, scorer, lambda: ())

    def test_occurrence_witness_survives_proof_only_projection_and_replay(self):
        from components.candidate_origins import mark_repair_origin
        from components.initial_competition import validate_ordinary_witness_inventory
        source, categories, scorer = fixture()
        key = next(iter(categories[0]))
        core, donor = categories[0][key][0], categories[1][key][0]
        mark_repair_origin(donor[1], {
            "policy": "terminal-only-trim-v1", "status": "generated",
            "coordinate_system": "zero_based_half_open",
            "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry": [donor[1].list_repeat_starts,
                [a + len(r) for a, r in zip(donor[1].list_repeat_starts, donor[1].list_repeats)]],
        })
        groups, conflicts = initial_competition_groups([core, donor])
        self.assertFalse(conflicts)
        # Keep anchor's prior v1/unannotated observation. Only donor contributes
        # v2; its metadata is removed from the live selected row after commit.
        annotated = next(r for r in groups[0].rows if r[1] is donor[1])
        categories[1][key] = [annotated]
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        selected = output[0][key][0]
        self.assertEqual(audit["stats"]["selections_committed"], 1)
        basis = selected[3]["classification_basis"]
        self.assertIsNone(basis["reconstruction"]["anchor_initial_competition"])
        self.assertIsNotNone(basis["reconstruction"]["donor_initial_competition"])
        self.assertNotIn("initial_competition", selected[3])
        self.assertEqual(len(validate_ordinary_witness_inventory(output, source)), 1)
        self.assertEqual(basis, validate_reconstruction_basis(
            source, selected[1], basis, selected[0], 0, selected, output))
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "protected.fa"
            fasta.write_text(">source\n" + source + "\n")
            projected = [output[0], {}, {}, {}, {}]
            report = ReportModelBuilder(fasta, projected, {}, [COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS, param_accepted_only=True), provenance_categories=output).build()
            self.assertEqual(report["arrays"][0]["classification_basis"], basis)
        for field in ("source_span_sha256", "representation_sha256", "score"):
            corrupt = copy.deepcopy(basis)
            ref = corrupt["reconstruction"]["donor_initial_competition"]["occurrence_protection"]["ordinary_reference"]
            ref[field] = .51 if field == "score" else "0" * 64
            corrupt["evidence_sha256"] = hashlib.sha256(json.dumps(corrupt["reconstruction"],
                sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
            with self.subTest(field=field), self.assertRaises(ValueError):
                validate_reconstruction_basis(source, selected[1], corrupt, selected[0], 0, selected, output)
        without_original = [output[0], {}, {}, {}, {}]
        with self.assertRaisesRegex(ValueError, "ordinary occurrence witness"):
            validate_ordinary_witness_inventory(without_original, source)

    def test_proof_only_deferred_summary_matches_retained_repair(self):
        from components.candidate_origins import mark_repair_origin
        from components.initial_competition import validate_ordinary_witness_inventory
        source, categories, scorer = fixture()
        key = next(iter(categories[0]))
        core, donor = categories[0][key][0], categories[1][key][0]
        starts = [donor[1].list_repeat_starts[i] for i in (0,1,3)]
        length = len(core[1].list_repeats[0])
        proposed = CrisprCandidate(list_repeat_starts=starts,
            list_repeats=[source[a:a+length] for a in starts],
            list_repeats_gaped=[source[a:a+length] for a in starts],
            list_spacers=[source[a+length:b] for a,b in zip(starts, starts[1:])])
        repaired = scorer({key:[proposed]})[key][0]
        mark_repair_origin(proposed, {
            "policy":"terminal-only-trim-v1", "status":"generated",
            "coordinate_system":"zero_based_half_open",
            "source_sha256":hashlib.sha256(source.encode()).hexdigest(),
            "child_geometry":[starts, [a+length for a in starts]],
        })
        groups, conflicts = initial_competition_groups([core, donor, repaired])
        self.assertEqual(len(conflicts), 1)
        categories[0][key] = [with_competition_evidence(core, groups[0].evidence)]
        categories[1][key].extend(conflicts)
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selections_committed"], 1)
        selected = output[0][key][0]
        self.assertNotIn("initial_competition", selected[3])
        proof = selected[3]["classification_basis"]["reconstruction"]
        self.assertEqual(len(proof["anchor_initial_competition"]["occurrence_protection"]["deferred_candidates"]), 1)
        self.assertEqual(len(validate_ordinary_witness_inventory(output, source)), 1)
        corrupt = copy.deepcopy(output)
        current = corrupt[0][key][0]
        current[3]["classification_basis"]["reconstruction"]["anchor_initial_competition"]["occurrence_protection"]["deferred_candidates"][0]["candidate"]["score"] = .91
        with self.assertRaisesRegex(ValueError, "retained deferred repair"):
            validate_ordinary_witness_inventory(corrupt, source)

    def test_atomic_replacement_keeps_own_score_and_archives_core(self):
        source, categories, scorer = fixture()
        old = next(iter(categories[0].values()))[0]
        alternative = next(iter(categories[1].values()))[0]
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        new = next(iter(output[0].values()))[0]
        self.assertFalse(audit["experimental"])
        self.assertFalse(audit["biological_accuracy_established"])
        self.assertIs(old, next(iter(categories[0].values()))[0])
        self.assertIs(alternative, next(iter(categories[1].values()))[0])
        self.assertIs(new[1], alternative[1])
        self.assertEqual(new[0], alternative[0])
        self.assertNotEqual(new[0], old[0])
        self.assertIn(old, output[1][(0, len(source))])
        self.assertEqual(METHOD, new[3]["classification_basis"]["method"])
        self.assertEqual(1, audit["stats"]["selections_committed"])
        self.assertEqual(2, audit["stats"]["fresh_candidates_scored"])

    def test_report_proof_replays_and_rejects_source_or_score_changes(self):
        source, categories, scorer = fixture()
        output, _ = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        row = next(iter(output[0].values()))[0]
        basis = row[3]["classification_basis"]
        self.assertEqual(basis, validate_reconstruction_basis(
            source, row[1], basis, row[0], 0, row, output))
        for changed_source, score in (("T" + source[1:], row[0]), (source, row[0] + .01)):
            with self.assertRaises(ValueError):
                validate_reconstruction_basis(changed_source, row[1], basis, score, 0, row, output)
        no_core = list(output)
        no_core[1] = {}
        with self.assertRaisesRegex(ValueError, "preserved original core"):
            validate_reconstruction_basis(source, row[1], basis, row[0], 0, row, no_core)

    def test_applicable_coding_screen_proof_replays_before_and_after_json(self):
        # Six 40-base repeats separated by distinct 35-base spacers trigger
        # the source-CDS screen (five adjacent periods of 75). An assessed
        # empty CDS inventory clears that scoped screen; it is not missing.
        # This must replay in memory as well as after JSON serialization,
        # because the tandem decision dataclass contains tuple coordinates.
        source, categories, scorer = fixture(equal_period=True)
        coding_requests = []
        def assessed_empty_cds():
            coding_requests.append(True)
            return ()
        output, audit = reconstruct_accepted_arrays(
            source, categories, PARAMETERS, scorer, assessed_empty_cds)
        self.assertEqual(coding_requests, [True])
        self.assertEqual(audit["stats"]["selections_committed"], 1)
        row = next(iter(output[0].values()))[0]
        basis = row[3]["classification_basis"]
        coding = basis["reconstruction"]["reconstruction_validation"]["coding_tandem"]
        self.assertEqual(coding["status"], "passed")
        self.assertTrue(coding["source_cds_assessed"])
        self.assertEqual(coding["considered_source_cds"], [])
        self.assertEqual(len(coding["decision"]["period_unit_intervals"]), 0)
        for recorded_basis in (basis, json.loads(json.dumps(basis))):
            with self.subTest(serialized=recorded_basis is not basis):
                self.assertEqual(recorded_basis, validate_reconstruction_basis(
                    source, row[1], recorded_basis, row[0], 0, row, output))

    def test_altered_negative_attestation_fails_even_with_recomputed_digest(self):
        source, categories, scorer = fixture()
        output, _ = reconstruct_accepted_arrays(source, categories, PARAMETERS, scorer, lambda: ())
        row = next(iter(output[0].values()))[0]
        basis = copy.deepcopy(row[3]["classification_basis"])
        basis["reconstruction"]["reconstruction_validation"]["likely_tandem"]["status"] = "inapplicable"
        basis["evidence_sha256"] = hashlib.sha256(json.dumps(
            basis["reconstruction"], sort_keys=True, separators=(",", ":"),
            allow_nan=False).encode()).hexdigest()
        with self.assertRaisesRegex(ValueError, "checks differ"):
            validate_reconstruction_basis(source, row[1], basis, row[0], 0, row, output)

    def test_eligibility_reasons_round_trip_for_both_core_and_reconstruction(self):
        source, categories, scorer = fixture()
        core = next(iter(categories[0].values()))[0]
        donor = next(iter(categories[1].values()))[0]
        core[1].source_path_classification_eligible = True
        core[1].source_path_classification_reason = "source_core_revalidated"
        donor[1].source_path_classification_eligible = True
        donor[1].source_path_classification_reason = "source_extent_revalidated"
        output, audit = reconstruct_accepted_arrays(
            source, categories, PARAMETERS, scorer, lambda: ())
        self.assertEqual(audit["stats"]["selections_committed"], 1)
        selected = next(iter(output[0].values()))[0]
        basis = selected[3]["classification_basis"]
        proof = basis["reconstruction"]
        for key in ("anchor_original", "anchor_fresh"):
            self.assertEqual(proof[key]["eligibility_reason"], "source_core_revalidated")
        for key in ("donor_original", "reconstruction"):
            self.assertEqual(proof[key]["eligibility_reason"], "source_extent_revalidated")
        self.assertEqual(basis, validate_reconstruction_basis(
            source, selected[1], basis, selected[0], 0, selected, output))

    def test_selection_claims_fail_replay_even_with_recomputed_digest(self):
        # This contract verifies the chosen path and anti-fusion evidence. It
        # does not claim to replay ranking over the full earlier donor universe.
        source, categories, scorer = fixture()
        output, _ = reconstruct_accepted_arrays(
            source, categories, PARAMETERS, scorer, lambda: ())
        row = next(iter(output[0].values()))[0]
        original_basis = row[3]["classification_basis"]
        corruptions = (
            (("selection", "relationship", "kind"), "internal_repeat_insertion", "relationship differs"),
            (("selection", "relationship", "correspondence", "anchor_to_candidate_indexes"),
             [99, 98, 97], "relationship differs"),
            (("selection", "selected", "geometry_id"), "0" * 64, "location differs"),
            (("selection", "anchor", "repeat_intervals"), [[0, 1]], "location differs"),
            (("selection", "selected", "score"), 0.123, "location differs"),
            (("selection", "reconstruction_score"), 0.123, "score differs"),
            (("anchor_original", "repeat_intervals"), [[0, 1]], "fresh scoring changed"),
            (("donor_original", "gapped_repeats"), ["ACGT"], "fresh scoring changed"),
            (("anchor_fresh", "eligibility_reason"), "changed_authorization", "fresh scoring changed"),
        )
        for path, value, expected_error in corruptions:
            with self.subTest(path=path):
                basis = copy.deepcopy(original_basis)
                destination = basis["reconstruction"]
                for key in path[:-1]:
                    destination = destination[key]
                destination[path[-1]] = value
                basis["evidence_sha256"] = hashlib.sha256(json.dumps(
                    basis["reconstruction"], sort_keys=True, separators=(",", ":"),
                    allow_nan=False).encode()).hexdigest()
                with self.assertRaisesRegex(ValueError, expected_error):
                    validate_reconstruction_basis(
                        source, row[1], basis, row[0], 0, row, output)

    def test_accepted_only_report_replays_with_full_provenance_categories(self):
        source, categories, scorer = fixture()
        provenance, _ = reconstruct_accepted_arrays(
            source, categories, PARAMETERS, scorer, lambda: ())
        final, _records, _stats = apply_likely_tandem_quarantine(provenance)
        projected = [final[0], {}, final[2], {}, {}]
        selected = next(iter(projected[0].values()))[0]
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "reconstruction.fa"
            fasta.write_text(">reconstruction\n" + source + "\n", encoding="ascii")
            kwargs = {
                "file_path": str(fasta), "categories": projected,
                "non_array_data": {}, "list_feature_names": [COMPACT_FEATURE_SCHEMA],
                "parameters": dict(PARAMETERS, param_accepted_only=True), "flags": {},
            }
            report = ReportModelBuilder(
                **kwargs, provenance_categories=provenance).build()
            self.assertEqual(report["summary"]["accepted_array_count"], 1)
            self.assertEqual(report["arrays"][0]["repeat_count"], 6)
            self.assertEqual(report["arrays"][0]["certainty_score"], selected[0])
            self.assertEqual(report["arrays"][0]["classification_basis"]["method"], METHOD)
            # Projected output alone intentionally cannot attest its archived
            # core: supplying provenance is a functional reporting requirement.
            with self.assertRaisesRegex(ReportModelError, "preserved original core"):
                ReportModelBuilder(**kwargs).build()

    def test_scoring_mutation_and_missing_results_abort_before_commit(self):
        for mode in ("mutation", "missing"):
            source, categories, scorer = fixture()
            original = next(iter(categories[0].values()))[0]
            def broken(request):
                result = scorer(request)
                if mode == "mutation":
                    next(iter(request.values()))[0].list_repeat_starts[0] += 1
                else:
                    result = {}
                return result
            with self.assertRaises(ValueError):
                reconstruct_accepted_arrays(source, categories, PARAMETERS, broken, lambda: ())
            self.assertIs(original, next(iter(categories[0].values()))[0])

    def test_no_donor_avoids_model_and_coding_work(self):
        source, categories, _ = fixture()
        categories[1] = {}
        def unexpected(*args):
            self.fail("unnecessary external work")
        output, audit = reconstruct_accepted_arrays(source, categories, PARAMETERS, unexpected, unexpected)
        self.assertEqual(0, audit["stats"]["fresh_candidates_scored"])
        self.assertEqual(categories, output)


if __name__ == "__main__":
    unittest.main()
