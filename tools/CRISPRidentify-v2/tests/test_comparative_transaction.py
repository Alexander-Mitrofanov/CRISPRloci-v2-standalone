"""Round-trip source and model evidence for a competing accepted geometry."""

import copy
import hashlib
import unittest
from unittest.mock import patch

from components.refinement import comparative_transaction as transaction
from components.refinement.potential_degenerate import geometry
from components.rescue.architecture_rescue_coding_context import SourceCDS
from tests.test_candidate_evidence import PARAMETERS
from tests.test_reconstruction_transaction import fixture


class CompetingGeometryTests(unittest.TestCase):
    def test_ranked_primary_and_rival_are_replayed_in_their_published_roles(self):
        dna, categories, score = fixture()
        key, rows = next(iter(categories[0].items()))
        parent = rows[0]
        parent.append({"classification_basis": {
            "method": "model_score_band", "final_category_before_source_transactions": "Bona-fide",
        }})
        uniform = geometry(categories[1][key][0][1])[:4]
        a, b = uniform[-1]
        compact = uniform[:-1] + ((a, b - 1),)
        details = {
            "alternative_repeat_intervals": [list(p) for p in uniform],
            "spacer_lengths": [31, 32, 33], "repeat_edit_distances": [0, 0, 0, 1],
            "alternative_metrics": {
                "spacer_lengths": [31, 32, 33], "repeat_edit_distances": [0, 0, 0, 0],
                "interpretation": "uniform_repeat_width",
            },
        }
        choice = (compact, "retained_indel_allele_extension", details)
        params = dict(PARAMETERS, param_max_edit_distance=6)
        digest = hashlib.sha256(dna.encode()).hexdigest()
        with patch.object(transaction, "_selection", return_value=choice):
            result = transaction.apply(dna, categories, params, score, lambda: ())
            cat = next(i for i in (0, 2) if key in result[i])
            entry = result[cat][key][0]
            basis = entry[3]["classification_basis"]
            self.assertEqual(geometry(entry[1]), uniform)
            self.assertEqual(basis["alternative"]["candidate"]["repeat_intervals"], [list(p) for p in compact])
            self.assertEqual(transaction.validate(dna, basis, entry, cat, result, digest), basis)
            changed = copy.deepcopy(basis)
            changed["alternative"]["candidate"] = changed["selected"]
            with self.assertRaisesRegex(ValueError, "alternative intervals"):
                transaction.validate(dna, changed, entry, cat, result, digest)

    def test_original_witness_and_fresh_score_are_required_for_replay(self):
        dna, categories, score = fixture()
        key, rows = next(iter(categories[0].items()))
        parent = rows[0]
        parent.append(
            {
                "classification_basis": {
                    "method": "model_score_band",
                    "final_category_before_source_transactions": "Bona-fide",
                }
            }
        )
        path = geometry(categories[1][key][0][1])
        selected = (
            path,
            "regular_core_divergent_terminal",
            {"equal_edit_alternative_boundary": True},
        )
        params = dict(PARAMETERS, param_max_edit_distance=6)
        with patch.object(transaction, "_selection", return_value=selected):
            result = transaction.apply(dna, categories, params, score, lambda: ())
            self.assertIs(categories[0][key][0], parent)
            self.assertIn(parent, result[1][key])
            cat = next(i for i in (0, 2) if key in result[i])
            entry = result[cat][key][0]
            basis = entry[3]["classification_basis"]
            digest = hashlib.sha256(dna.encode()).hexdigest()
            self.assertEqual(
                transaction.validate(dna, basis, entry, cat, result, digest), basis
            )
            for field, bad in (
                ("rule", "other"),
                ("geometry_status", "certain"),
                ("previous_provenance_sha256", "0" * 64),
                ("alternative", {}),
            ):
                changed = copy.deepcopy(basis)
                changed[field] = bad
                with self.subTest(field=field), self.assertRaises(ValueError):
                    transaction.validate(dna, changed, entry, cat, result, digest)
            without_witness = list(result)
            without_witness[1] = {}
            with self.assertRaisesRegex(ValueError, "parent witness"):
                transaction.validate(dna, basis, entry, cat, without_witness, digest)

    def test_coding_veto_for_a_regular_split_in_non_in_frame_array(self):
        dna, cats, _ = fixture()
        parent = next(iter(cats[0].values()))[0][1]
        old = geometry(parent)
        lo, hi = old[0][1], old[1][0]
        details = {"split_spacer_interval": [lo, hi]}

        def scope(start, end):
            cds = SourceCDS(start, end, False, 90.0)
            return transaction._cds_scope(
                dna, parent, "periodic_internal_insertion", details, lambda: (cds,)
            )

        self.assertIsNone(scope(lo, hi))
        self.assertEqual(len(scope(lo, hi - 1)), 1)

    def test_below_band_support_only_applies_to_designated_possible_cores(self):
        self.assertEqual(
            transaction._category(2, 0.54, 0.42, "periodic_internal_insertion"), 2
        )
        self.assertIsNone(
            transaction._category(0, 0.95, 0.42, "periodic_internal_insertion")
        )
        self.assertIsNone(
            transaction._category(2, 0.54, 0.42, "regular_core_divergent_terminal")
        )
        self.assertIsNone(
            transaction._category(2, 0.40, 0.85, "periodic_internal_insertion")
        )
