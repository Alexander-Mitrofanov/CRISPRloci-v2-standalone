import copy
import hashlib
import json
from unittest import TestCase

from components.reconstruction.reconstruction_boundary_family import geometry_id
from components.reporting.tied_occurrences import (
    AUDIT_SCHEMA,
    COORDINATES,
    EXTENSION,
    audit_text,
    tied_reconstruction_occurrences,
    validate_tied_reconstruction_occurrences,
)


def location(geometry, category):
    return dict(
        category_index=category,
        geometry_id=geometry_id(geometry),
        source_span=[geometry[0][0], geometry[-1][1]],
        repeat_count=len(geometry),
        repeat_intervals=geometry,
        score=0.8,
        locus_key=[0, 75],
    )


class TiedOccurrenceTests(TestCase):
    def setUp(self):
        self.source = "ACGT" * 20
        core = [[18, 22], [36, 40]]
        a = [[3, 7], *core, [51, 55], [68, 72]]
        b = [[4, 8], *core, [52, 56], [68, 72]]
        self.array = {
            "id": "anchor",
            "category": "Bona-fide",
            "certainty_score": 0.8,
            "legacy_interval_key": [0, 75],
            "repeats": [
                {
                    "source_interval": {"start": x + 1, "end": y},
                    "sequence": self.source[x:y],
                }
                for x, y in core
            ],
        }
        self.report = {
            "source": {
                "id": "chr",
                "length": len(self.source),
                "sha256": hashlib.sha256(self.source.encode()).hexdigest(),
            },
            "arrays": [self.array],
            "extensions": {},
        }
        policy = {
            "minimum_anchor_score": 0.75,
            "minimum_candidate_score": 0.60,
            "score_tie_tolerance": 1e-12,
        }
        self.peers = [
            {"candidate": location(p, 1), "vetoes": [], "fresh_score": 0.8}
            for p in (a, b)
        ]
        self.evidence = {
            "policy": policy,
            "freshly_scored_peer_rows": 2,
            "direct_peer_rows": 2,
            "direct_peer_validation_complete": True,
            "anchor_validation": {
                k: True
                for k in (
                    "source_valid",
                    "structural_valid",
                    "negatives_passed",
                    "score_fresh",
                    "qualified",
                )
            },
            "anchor_score": 0.8,
            "best_reconstruction_score": 0.8,
            "best_repeat_count": len(a),
            "same_rank_peer_geometries": 2,
            "same_rank_peers": self.peers,
        }
        self.decision = {
            "action": "abstained",
            "reason": "ambiguous_best_geometry",
            "coordinate_system": COORDINATES,
            "selected": None,
            "rule": "conservative-core-reconstruction",
            "tied_geometry_ids": sorted(geometry_id(p) for p in (a, b)),
            "evidence": self.evidence,
            "anchor": location(core, 0),
        }
        self.audit = {
            "schema": AUDIT_SCHEMA,
            "coordinate_system": COORDINATES,
            "source_sha256": self.report["source"]["sha256"],
            "policy": policy,
            "decisions": [self.decision],
        }

    def build(self):
        return tied_reconstruction_occurrences(self.report, self.audit, self.source)

    def test_intersection_skips_contested_intervals_and_preserves_accepted_arrays(self):
        original = copy.deepcopy(self.report)
        output = self.build()
        (row,) = output["anchors"]
        self.assertEqual(row["common_repeat_intervals"], [[18, 22], [36, 40], [68, 72]])
        self.assertEqual(row["additional_repeat_intervals"], [[68, 72]])
        self.assertEqual(row["tied_geometry_ids"], self.decision["tied_geometry_ids"])
        ledger = audit_text(self.audit)
        self.assertEqual(
            output["reconstruction_evidence"]["sha256"],
            hashlib.sha256(ledger.encode()).hexdigest(),
        )
        self.assertEqual(self.report, original)
        self.report["extensions"][EXTENSION] = output
        validate_tied_reconstruction_occurrences(
            self.report, json.loads(ledger), self.source
        )

    def test_stale_source_public_anchor_and_modified_intersection_fail_replay(self):
        self.report["extensions"][EXTENSION] = self.build()
        with self.assertRaisesRegex(ValueError, "source binding"):
            validate_tied_reconstruction_occurrences(
                self.report, self.audit, self.source[::-1]
            )
        self.array["id"] = "different"
        with self.assertRaisesRegex(ValueError, "failed replay"):
            validate_tied_reconstruction_occurrences(
                self.report, self.audit, self.source
            )
        self.array["id"] = "anchor"
        self.report["extensions"][EXTENSION]["anchors"][0][
            "common_repeat_intervals"
        ].append([51, 55])
        with self.assertRaisesRegex(ValueError, "failed replay"):
            validate_tied_reconstruction_occurrences(
                self.report, self.audit, self.source
            )

    def test_missing_peer_invalid_geometry_and_veto_fail_closed(self):
        peer = self.peers.pop()
        with self.assertRaisesRegex(ValueError, "tie set"):
            self.build()
        self.peers.append(peer)
        peer["vetoes"] = ["overlapping_accepted_array"]
        with self.assertRaisesRegex(ValueError, "vetoed"):
            self.build()
        peer["vetoes"] = []
        peer["candidate"]["repeat_intervals"][-1][1] = 81
        with self.assertRaisesRegex(ValueError, "outside source"):
            self.build()

    def test_absent_audit_replaced_anchor_other_abstentions_and_empty_support(self):
        self.assertIsNone(
            tied_reconstruction_occurrences(self.report, None, self.source)
        )
        self.array["category"] = "Low score"
        self.assertIsNone(self.build())
        self.array["category"] = "Bona-fide"
        self.decision["reason"] = "same_rank_peer_veto"
        self.assertIsNone(self.build())
        self.decision["reason"] = "ambiguous_best_geometry"
        changed = [[4, 8], [18, 22], [36, 40], [52, 56], [69, 73]]
        self.peers[1]["candidate"] = location(changed, 1)
        self.decision["tied_geometry_ids"] = sorted(
            p["candidate"]["geometry_id"] for p in self.peers
        )
        self.assertIsNone(self.build())

    def test_wrong_digest_and_unvalidated_rank_do_not_publish(self):
        original = self.peers[0]["candidate"]["geometry_id"]
        self.peers[0]["candidate"]["geometry_id"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "geometry ID"):
            self.build()
        self.peers[0]["candidate"]["geometry_id"] = original
        self.evidence["direct_peer_validation_complete"] = False
        with self.assertRaisesRegex(ValueError, "direct-peer"):
            self.build()


class TiedOccurrencePublicationTests(TestCase):
    def test_machine_profile_publishes_native_tie_and_matching_ledger(self):
        import tempfile
        from pathlib import Path

        from components.reconstruction.reconstruction_transaction import (
            reconstruct_accepted_arrays,
        )
        from components.refinement.candidates import CrisprCandidate
        from components.reporting.stage import OutputMaker
        from components.scoring.feature_schema import COMPACT_FEATURE_SCHEMA
        from tests.test_candidate_evidence import PARAMETERS
        from tests.test_extent_family_exports import annotations
        from tests.test_reconstruction_transaction import fixture

        source, original, scorer = fixture()
        (key,) = original[1]
        first = original[1][key][0][1]
        starts = list(first.list_repeat_starts)
        starts[-1] += 1
        repeats = list(first.list_repeats)
        repeats[-1] = repeats[-1][1:]
        gapped = list(repeats)
        gapped[-1] = "-" + gapped[-1]
        second = CrisprCandidate(
            list_repeat_starts=starts,
            list_repeats=repeats,
            list_repeats_gaped=gapped,
            list_spacers=[
                source[a + len(r) : b] for a, r, b in zip(starts, repeats, starts[1:])
            ],
        )
        original[1][key].extend(scorer({key: [second]})[key])
        categories, audit = reconstruct_accepted_arrays(
            source, original, PARAMETERS, scorer, lambda: ()
        )
        self.assertEqual(audit["decisions"][0]["reason"], "ambiguous_best_geometry")
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            fasta = root / "source.fa"
            fasta.write_text(">source\n" + source + "\n")
            output = root / "out"
            output.mkdir()
            OutputMaker(
                file_path=fasta,
                parameters=PARAMETERS,
                flags={"output_profile": "machine", "flag_strand": True},
                result_path=str(output),
                pickle_result_path=str(root / "pickle"),
                json_result_path="",
                categories=categories,
                non_array_data=annotations(),
                list_features=[COMPACT_FEATURE_SCHEMA],
                header=">source\n",
                reconstruction_evidence=audit,
            )
            report = json.loads((output / "report.json").read_text())
            observed = report["extensions"][EXTENSION]
            ledger = (
                output / observed["reconstruction_evidence"]["filename"]
            ).read_bytes()
            self.assertEqual(
                hashlib.sha256(ledger).hexdigest(),
                observed["reconstruction_evidence"]["sha256"],
            )
            self.assertEqual(len(observed["anchors"][0]["common_repeat_intervals"]), 5)
            self.assertEqual(
                len(observed["anchors"][0]["additional_repeat_intervals"]), 2
            )
            self.assertEqual(report["summary"]["accepted_array_count"], 1)
            self.assertEqual(report["arrays"][0]["repeat_count"], 3)
            validate_tied_reconstruction_occurrences(report, json.loads(ledger), source)
