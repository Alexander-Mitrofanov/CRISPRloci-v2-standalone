import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from components.reconstruction.family_terminal_extension import terminal_projection
from components.reconstruction.reconstruction_boundary_family import (
    digest,
    validate_boundary_family_audit,
)
from components.reconstruction.reconstruction_transaction import (
    reconstruct_accepted_arrays,
    validate_boundary_family_basis,
)
from components.refinement.potential_degenerate import geometry
from components.reporting.report_model import ReportModelBuilder
from components.scoring.feature_schema import COMPACT_FEATURE_SCHEMA
from tests.test_candidate_evidence import PARAMETERS
from tests.test_reconstruction_boundary_family import fixture


class TerminalFamilyTests(unittest.TestCase):
    @staticmethod
    def run_family(*, fresh_high=True, divergent=False):
        source, categories, scorer = fixture()
        key = next(iter(categories[0]))
        core = categories[0][key][0]
        terminal = core[1].list_repeats[0]
        if divergent:
            substitutions = str.maketrans("ACGT", "CGTA")
            terminal = terminal[:14] + terminal[14:18].translate(substitutions) + terminal[18:]
        source += terminal

        def score(request):
            batch = scorer(request)
            for rows in batch.values():
                for r in rows:
                    # A fixed synthetic feature vector with a real model score
                    # exercises the family transaction's per-member confidence gate.
                    if fresh_high and len(r[1].list_repeats) == 7:
                        r[0], r[2] = core[0], core[2]
            return batch

        output, audit = reconstruct_accepted_arrays(
            source,
            categories,
            dict(PARAMETERS, param_degenerated=True, param_max_edit_distance=6),
            score,
            lambda: (),
        )
        return source, output, audit, output[0][key][0]

    def test_complete_family_survives_native_proof_and_report(self):
        source, output, audit, row = self.run_family(divergent=True)
        basis = row[3]["classification_basis"]
        proof = basis["boundary_family"]
        ext = proof["terminal_extension"]
        originals = {
            m["location"]["geometry_id"]: m["fresh"]["repeat_intervals"]
            for m in ext["parent_members"]
        }
        self.assertEqual(
            {m["parent_geometry_id"] for m in proof["members"]}, set(originals)
        )
        for m in proof["members"]:
            self.assertEqual(
                m["fresh"]["repeat_intervals"],
                originals[m["parent_geometry_id"]] + [ext["added_interval"]],
            )
        self.assertEqual(
            validate_boundary_family_basis(
                source, row[1], basis, row[0], 0, row, output
            ),
            basis,
        )
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / "source.fa"
            path.write_text(">source\n" + source + "\n")
            report = ReportModelBuilder(
                path,
                [output[0], {}, {}, {}, {}],
                {},
                [COMPACT_FEATURE_SCHEMA],
                parameters=dict(PARAMETERS, param_accepted_only=True),
                provenance_categories=output,
            ).build()
        self.assertEqual(
            report["arrays"][0]["boundary_family"]["variable_repeat_indexes"], [3, 4]
        )
        self.assertEqual(report["arrays"][0]["repeat_count"], 7)
        repeats = report["arrays"][0]["repeats"]
        self.assertTrue(all("status" not in repeat for repeat in repeats[:-1]))
        self.assertEqual(repeats[-1]["status"], "potential_degenerate_repeat")
        self.assertEqual(repeats[-1]["location"], "right")
        self.assertEqual(repeats[-1]["ordinary_edit_limit"], 6)
        self.assertEqual(repeats[-1]["placement_method"], "local_spacer_period")
        self.assertEqual(repeats[-1]["edit_distance"], 4)
        validate_boundary_family_audit(report, audit)
        audit.pop("terminal_family_extensions")
        with self.assertRaisesRegex(ValueError, "terminal extensions differ"):
            validate_boundary_family_audit(report, audit)

    def test_weak_extension_retains_the_entire_original_family(self):
        _, _, audit, row = self.run_family(fresh_high=False)
        self.assertNotIn("terminal_family_extensions", audit)
        self.assertEqual(len(row[1].list_repeats), 6)
        self.assertNotIn(
            "terminal_extension", row[3]["classification_basis"]["boundary_family"]
        )

    def test_one_disagreeing_parent_vetoes_the_extension(self):
        source, cats, _scorer = fixture()
        members = [r[1] for r in next(iter(cats[1].values()))]
        source += members[0].list_repeats[0]

        def propose(_source, member, _ordinary):
            return geometry(member) + (
                ((len(source) - len(members[0].list_repeats[0]), len(source)),)
                if member is members[0]
                else ()
            )

        with patch(
            "components.reconstruction.family_terminal_extension.plan_path",
            side_effect=propose,
        ):
            self.assertIsNone(terminal_projection(source, members, 6))

    def test_changed_extension_cannot_be_relabelled_with_a_new_digest(self):
        source, output, _, row = self.run_family()
        basis = row[3]["classification_basis"]
        proof = basis["boundary_family"]
        proof["terminal_extension"]["added_interval"][0] -= 1
        basis["evidence_sha256"] = digest(proof)
        with self.assertRaisesRegex(ValueError, "terminal placement differs"):
            validate_boundary_family_basis(
                source, row[1], basis, row[0], 0, row, output
            )
