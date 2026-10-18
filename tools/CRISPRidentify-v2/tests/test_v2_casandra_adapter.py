import hashlib
import importlib.util
import json
import shutil
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "components" / "casandra_adapter.py"
)


def load_module():
    external = types.ModuleType("components.external_tools")
    bounded = types.ModuleType("components.bounded_subprocess")

    class ExternalToolError(RuntimeError):
        pass

    external.ExternalToolError = ExternalToolError

    class BoundedSubprocessTimeoutError(RuntimeError):
        pass

    class BoundedSubprocessTruncationError(RuntimeError):
        pass

    def run_bounded(*args, **kwargs):
        raise AssertionError("default bounded runner was not injected")

    bounded.BoundedSubprocessTimeoutError = BoundedSubprocessTimeoutError
    bounded.BoundedSubprocessTruncationError = (
        BoundedSubprocessTruncationError
    )
    bounded.run_bounded = run_bounded
    package = types.ModuleType("components")
    package.__path__ = []
    spec = importlib.util.spec_from_file_location(
        "v2_casandra_adapter_under_test", MODULE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(
        sys.modules,
        {
            "components": package,
            "components.bounded_subprocess": bounded,
            "components.external_tools": external,
            spec.name: module,
        },
    ):
        spec.loader.exec_module(module)
    return module


casandra = load_module()


def digest_bytes(value):
    return hashlib.sha256(value).hexdigest()


def json_bytes(value):
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode("utf-8")


def prediction(protein_id, profile, cas_type, subtype, positive, negative):
    is_cas = profile is not None
    profile_families = {
        "profile-cas3": ["cas3"],
        "profile-cas5": ["cas5"],
    }.get(profile, [])
    return {
        "best_positive_profile": profile,
        "classification": {
            "class": "1" if is_cas else None,
            "subtype": subtype,
            "type": cas_type,
        },
        "evidence": {
            "decision_threshold": 0.0,
            "model_id": casandra.CASANDRA_MODEL_ID,
            "positive_profile_required": True,
            "profile_hits": 1 if is_cas else 0,
            "report_evalue": 1000.0,
            "type_v_subtype_masked": False,
            "raw_classification": {
                "class": "1" if is_cas else None,
                "subtype": subtype,
                "type": cas_type,
            },
            "hierarchy_conflict_masked": False,
        },
        "hard_negative_profile_score": negative,
        "is_cas": is_cas,
        "profile_families": profile_families,
        "positive_profile_score": positive,
        "score_margin": positive - negative,
        "sequence_id": protein_id,
        "sequence_length": 2,
        "variant": None,
    }


class FixtureRunner:
    def __init__(self, input_path, installation, mutation=None, returncode=0):
        self.input_path = Path(input_path)
        self.installation = installation
        self.mutation = mutation
        self.returncode = returncode
        self.commands = []

    def __call__(self, command, **kwargs):
        self.commands.append((list(command), dict(kwargs)))
        if self.returncode:
            return types.SimpleNamespace(
                returncode=self.returncode,
                stdout="",
                stderr="deliberate fixture failure",
            )
        output = Path(command[command.index("--output") + 1])
        output.mkdir(parents=True)
        (output / "extraction").mkdir()
        source_sha256 = hashlib.sha256(self.input_path.read_bytes()).hexdigest()
        dna = "ATGAAATAA" * 3 + "A" * 5001 + "ATGAAATAA"
        identifiers = ["p1", "p2", "p3", "p4"]
        predictions = [
            prediction("p1", "profile-cas3", "I", "I-A", 15.0, -2.0),
            prediction("p2", None, None, None, -200.0, -200.0),
            prediction("p3", "profile-cas5", "I", "I-A", 9.0, -3.0),
            prediction("p4", "profile-cas3", "I", None, 7.0, -4.0),
        ]
        intervals = ((1, 9), (10, 18), (19, 27), (5029, 5037))
        proteins = []
        for protein_id, predicted, interval in zip(
                identifiers, predictions, intervals):
            start, end = interval
            nucleotide = dna[start - 1:end]
            protein_sequence = "MK"
            proteins.append({
                "caller": {
                    "closed": False,
                    "mask": True,
                    "mode": "single",
                    "program": "Pyrodigal",
                    "requested_mode": "auto",
                    "translation_table": 11,
                    "version": "3.7.1",
                },
                "contig_id": "record",
                "end_1based_inclusive": end,
                "genome_source_sha256": source_sha256,
                "nucleotide_sequence": nucleotide,
                "nucleotide_sha256": digest_bytes(nucleotide.encode("ascii")),
                "prediction": predicted,
                "protein_id": protein_id,
                "protein_length": 2,
                "protein_sequence": protein_sequence,
                "protein_sha256": digest_bytes(
                    protein_sequence.encode("ascii")
                ),
                "start_1based": start,
                "strand": "+",
                "translation_table": 11,
            })

        accepted_direct = {
            "abstained_subtype": False,
            "class": "1",
            "confidence": {
                "class_vote_fraction": 1.0,
                "subtype_vote_fraction": 1.0,
                "type_vote_fraction": 1.0,
            },
            "raw_subtype_winner": "I-A",
            "subtype": "I-A",
            "subtype_supporting_proteins": ["p1", "p3"],
            "type": "I",
            "votes": {
                "class": {"1": 24.0},
                "subtype": {"I-A": 24.0},
                "type": {"I": 24.0},
            },
        }
        rejected_direct = {
            "abstained_subtype": False,
            "class": "1",
            "confidence": {
                "class_vote_fraction": 1.0,
                "subtype_vote_fraction": 0.0,
                "type_vote_fraction": 1.0,
            },
            "raw_subtype_winner": None,
            "subtype": None,
            "subtype_supporting_proteins": [],
            "type": "I",
            "votes": {
                "class": {"1": 7.0},
                "subtype": {},
                "type": {"I": 7.0},
            },
        }
        construction = {
            "crispr_array_evidence_used": False,
            "max_cas_gap_nt": 5000,
            "max_non_cas_genes": 5,
        }
        cassettes = [{
            "cas_gene_count": 2,
            "cas_protein_ids": ["p1", "p3"],
            "cassette_id": "casandra|contig=record|cassette=0001|loc=1-27",
            "classification": accepted_direct,
            "construction": dict(construction),
            "context_gene_count": 3,
            "context_protein_ids": ["p1", "p2", "p3"],
            "contig_id": "record",
            "end_1based_inclusive": 27,
            "final_classification": {
                "class": "1",
                "confidence": 1.0,
                "direct_profile_result": json.loads(
                    json.dumps(accepted_direct)
                ),
                "evidence": {
                    "raw_classification": {
                        "class": "1", "type": "I", "subtype": "I-A",
                    },
                    "hierarchy_conflict_masked": False,
                },
                "method": "direct_profile_aggregation",
                "subtype": "I-A",
                "type": "I",
            },
            "genome_evidence_gate": {
                "accepted": True,
                "rule": "class1_requires_two_cas_genes",
            },
            "start_1based": 1,
        }]
        rejected = [{
            "cas_gene_count": 1,
            "cas_protein_ids": ["p4"],
            "cassette_id": (
                "casandra|contig=record|cassette=0002|loc=5029-5037"
            ),
            "classification": rejected_direct,
            "construction": dict(construction),
            "context_gene_count": 1,
            "context_protein_ids": ["p4"],
            "contig_id": "record",
            "end_1based_inclusive": 5037,
            "final_classification": {
                "class": "1",
                "confidence": 1.0,
                "direct_profile_result": json.loads(
                    json.dumps(rejected_direct)
                ),
                "evidence": {
                    "raw_classification": {
                        "class": "1", "type": "I", "subtype": None,
                    },
                    "hierarchy_conflict_masked": False,
                },
                "method": "direct_profile_type_aggregation",
                "subtype": None,
                "type": "I",
            },
            "genome_evidence_gate": {
                "accepted": False,
                "rule": "no_supported_subtype",
            },
            "start_1based": 5029,
        }]
        run = {
            "bundle_id": self.installation.bundle_id,
            "bundle_role": self.installation.bundle_role,
            "cas_cassettes": 1,
            "cas_proteins": 3,
            "classified_cassettes": 1,
            "cpu_threads": 3,
            "crispr_array_prediction": False,
            "genes": 4,
            "hierarchy_conflict_masked_cassettes": 0,
            "input": {
                "mode": "pyrodigal_de_novo",
                "path": str(self.input_path),
                "sha256": source_sha256,
            },
            "offline_inference": True,
            "program": "CasAndra",
            "raw_cassette_candidates": 2,
            "rejected_cassette_candidates": 1,
            "schema_version": casandra.CASANDRA_SCHEMA_VERSION,
        }
        fixture = {
            "proteins": proteins,
            "predictions": predictions,
            "cassettes": cassettes,
            "rejected": rejected,
            "run": run,
        }
        if self.mutation is not None:
            self.mutation(fixture)

        payloads = {
            "run.json": json_bytes(run),
            "proteins.jsonl": b"".join(json_bytes(row) for row in proteins),
            "protein_predictions.jsonl": b"".join(
                json_bytes(row) for row in predictions
            ),
            "cassettes.jsonl": b"".join(
                json_bytes(row) for row in cassettes
            ),
            "rejected_cassette_candidates.jsonl": b"".join(
                json_bytes(row) for row in rejected
            ),
            "extraction/audit.json": json_bytes({
                "schema_version": 3,
                "selected_mode": "single",
            }),
            "extraction/manifest.json": json_bytes({
                "caller": {
                    "mode": "single",
                    "program": "Pyrodigal",
                    "requested_mode": "auto",
                    "translation_table": 11,
                    "version": "3.7.1",
                },
                "schema_version": 3,
            }),
        }
        for name, content in payloads.items():
            path = output / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
        manifest = {
            "bundle_manifest_sha256": (
                self.installation.model_manifest_sha256
            ),
            "files": {
                name: {
                    "sha256": digest_bytes(content),
                    "size": len(content),
                }
                for name, content in payloads.items()
            },
            "schema_version": casandra.CASANDRA_SCHEMA_VERSION,
        }
        (output / "manifest.json").write_bytes(json_bytes(manifest))
        if fixture.get("tamper_run_after_manifest"):
            (output / "run.json").write_text("{}\n", encoding="utf-8")
        return types.SimpleNamespace(
            returncode=0,
            stdout=json.dumps(run, indent=2, sort_keys=True),
            stderr="",
        )


def apply_real_mixed_type_vi_cassette(fixture):
    """Recreate NC_008709.1 cassette 0058's independent vote winners."""

    specifications = (
        (0, "profile-cas3", "1", "III", None, 12.178770065307617),
        (1, "profile-cas3", "1", "I", None, 10.81540584564209),
        (2, "profile-cas5", "2", "VI", "VI-B", 12.440103530883789),
    )
    families = {"profile-cas3": ["cas3"], "profile-cas5": ["cas5"]}
    for index, profile, cas_class, cas_type, subtype, score in specifications:
        predicted = fixture["predictions"][index]
        classification = {
            "class": cas_class,
            "type": cas_type,
            "subtype": subtype,
        }
        predicted.update({
            "best_positive_profile": profile,
            "classification": classification,
            "hard_negative_profile_score": -2.0,
            "is_cas": True,
            "positive_profile_score": score,
            "profile_families": families[profile],
            "score_margin": score + 2.0,
        })
        predicted["evidence"].update({
            "hierarchy_conflict_masked": False,
            "profile_hits": 1,
            "raw_classification": dict(classification),
        })

    type_iii = specifications[0][-1]
    type_i = specifications[1][-1]
    type_vi = specifications[2][-1]
    class_1 = type_iii + type_i
    total = class_1 + type_vi
    direct = {
        "abstained_subtype": False,
        "class": "1",
        "confidence": {
            "class_vote_fraction": class_1 / total,
            "subtype_vote_fraction": 1.0,
            "type_vote_fraction": type_vi / total,
        },
        "raw_subtype_winner": "VI-B",
        "subtype": "VI-B",
        "subtype_supporting_proteins": ["p3"],
        "type": "VI",
        "votes": {
            "class": {"1": class_1, "2": type_vi},
            "subtype": {"VI-B": type_vi},
            "type": {"I": type_i, "III": type_iii, "VI": type_vi},
        },
    }
    cassette = fixture["cassettes"][0]
    cassette.update({
        "cas_gene_count": 3,
        "cas_protein_ids": ["p1", "p2", "p3"],
        "classification": direct,
        "final_classification": {
            "class": "2",
            "confidence": 1.0,
            "direct_profile_result": json.loads(json.dumps(direct)),
            "evidence": {
                "raw_classification": {
                    "class": "1", "type": "VI", "subtype": "VI-B",
                },
                "hierarchy_conflict_masked": True,
            },
            "method": "direct_profile_aggregation",
            "subtype": "VI-B",
            "type": "VI",
        },
        "genome_evidence_gate": {
            "accepted": True,
            "rule": "class2_multi_gene",
        },
    })
    fixture["run"]["cas_proteins"] = 4
    fixture["run"]["hierarchy_conflict_masked_cassettes"] = 1


class CasAndraAdapterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fasta = self.root / "input.fasta"
        self.dna = "ATGAAATAA" * 3 + "A" * 5001 + "ATGAAATAA"
        self.fasta.write_text(">record fixture\n{}\n".format(self.dna))
        self.installation = casandra.CasAndraInstallation(
            package_root=self.root / "casandra",
            version=casandra.CASANDRA_VERSION,
            bundle_id=casandra.CASANDRA_BUNDLE_ID,
            bundle_role=casandra.CASANDRA_BUNDLE_ROLE,
            model_id=casandra.CASANDRA_MODEL_ID,
            model_manifest_sha256=(
                casandra.CASANDRA_MODEL_MANIFEST_SHA256
            ),
            profile_metadata={
                "profile-cas3": {
                    "families": ["cas3"],
                    "intended_use": "primary_search",
                },
                "profile-cas5": {
                    "families": ["cas5"],
                    "intended_use": "primary_search",
                },
            },
        )

    def predictor(self, runner):
        return casandra.CasAndraPredictor(
            input_fasta=self.fasta,
            source_id="record",
            source_dna=self.dna,
            cpu_count=3,
            gene_mode="auto",
            installation=self.installation,
            command_runner=runner,
        )

    def test_cli_boundary_and_normalized_contract_include_orphan_hits(self):
        runner = FixtureRunner(self.fasta, self.installation)
        result = self.predictor(runner).run()

        command, options = runner.commands[0]
        self.assertEqual(
            [sys.executable, "-P", "-m", "casandra", "predict-genome"],
            command[:5],
        )
        self.assertEqual("auto", command[command.index("--gene-mode") + 1])
        self.assertEqual("3", command[command.index("--threads") + 1])
        self.assertEqual(3600.0, options["timeout"])
        self.assertEqual(1024 * 1024, options["stdout_limit"])
        self.assertEqual(64 * 1024, options["stderr_limit"])
        self.assertEqual("strict", options["errors"])
        self.assertIn("PYTHONPYCACHEPREFIX", options["env"])
        self.assertTrue(options["env"]["PYTHONPYCACHEPREFIX"].endswith("pycache"))

        self.assertEqual("computed", result["status"])
        self.assertEqual("CasAndra", result["method"]["provider"])
        self.assertEqual(casandra.CASANDRA_VERSION,
                         result["method"]["provider_version"])
        self.assertEqual("auto", result["method"]["gene_mode"])
        self.assertEqual("single", result["method"]["selected_gene_mode"])
        self.assertEqual(3600.0, result["method"]["timeout_seconds"])
        self.assertEqual(3, result["gene_count"])
        self.assertEqual(
            ["p1", "p3", "p4"],
            [gene["protein_id"] for gene in result["genes"]],
        )
        self.assertEqual(
            0,
            result["run_statistics"]["hierarchy_conflict_masked_proteins"],
        )
        self.assertEqual(
            0,
            result["run_statistics"]["hierarchy_conflict_masked_cassettes"],
        )
        self.assertEqual(
            {
                "raw_classification": {
                    "class": "1", "type": "I", "subtype": "I-A",
                },
                "hierarchy_conflict_masked": False,
            },
            result["genes"][0]["classification_evidence"],
        )
        system = result["systems"][0]
        self.assertEqual(1, system["cassette_id"])
        self.assertEqual(
            "casandra|contig=record|cassette=0001|loc=1-27",
            system["provider_cassette_id"],
        )
        self.assertEqual("CAS-I-A", system["type"])
        self.assertEqual(
            "positive_bitscore_weighted_subtype_vote_fraction",
            system["classification"]["confidence_kind"],
        )
        self.assertEqual(["cas3"], system["genes"][0]["profile_families"])
        self.assertEqual("cas3", system["genes"][0]["annotation"])
        self.assertEqual(
            [True, False, True],
            [gene["is_cas"] for gene in system["genes"]],
        )
        self.assertEqual("unknown", system["genes"][1]["annotation"])
        self.assertFalse(system["experimental"])
        self.assertIsNone(system["experimental_reason"])
        self.assertFalse(system["evidence"]["type_vii_experimental"])
        self.assertEqual(
            {
                "raw_classification": {
                    "class": "1", "type": "I", "subtype": "I-A",
                },
                "hierarchy_conflict_masked": False,
            },
            system["evidence"]["provider_final_classification_evidence"],
        )
        self.assertNotIn("probabilities", system)
        self.assertNotIn("bitscore", system["genes"][0])

    def test_real_mixed_type_vi_cassette_is_strictly_normalized(self):
        result = self.predictor(FixtureRunner(
            self.fasta,
            self.installation,
            mutation=apply_real_mixed_type_vi_cassette,
        )).run()

        system = result["systems"][0]
        self.assertEqual("CAS-VI-B", system["type"])
        self.assertEqual(
            {"class": "2", "type": "VI", "subtype": "VI-B"},
            {
                key: system["classification"][key]
                for key in ("class", "type", "subtype")
            },
        )
        self.assertEqual(
            {
                "raw_classification": {
                    "class": "1", "type": "VI", "subtype": "VI-B",
                },
                "hierarchy_conflict_masked": True,
            },
            system["evidence"]["provider_final_classification_evidence"],
        )
        self.assertEqual(
            "class2_multi_gene",
            system["evidence"]["genome_evidence_gate"],
        )
        self.assertEqual(
            1,
            result["run_statistics"]["hierarchy_conflict_masked_cassettes"],
        )

    def test_final_cassette_conflict_flag_is_recomputed(self):
        def lie_about_flag(fixture):
            apply_real_mixed_type_vi_cassette(fixture)
            fixture["cassettes"][0]["final_classification"]["evidence"][
                "hierarchy_conflict_masked"
            ] = False

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "final cassette hierarchy_conflict_masked disagrees"):
            self.predictor(FixtureRunner(
                self.fasta, self.installation, mutation=lie_about_flag
            )).run()

    def test_raw_final_cassette_must_match_direct_votes(self):
        def forge_raw_route(fixture):
            apply_real_mixed_type_vi_cassette(fixture)
            evidence = fixture["cassettes"][0]["final_classification"][
                "evidence"
            ]
            evidence["raw_classification"]["class"] = "2"
            evidence["hierarchy_conflict_masked"] = False

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "raw final class/type disagrees with direct profile votes"):
            self.predictor(FixtureRunner(
                self.fasta, self.installation, mutation=forge_raw_route
            )).run()

    def test_empty_positive_profile_family_evidence_is_rejected(self):
        def remove_families(fixture):
            fixture["predictions"][3]["profile_families"] = []

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=remove_families
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "families disagree"):
            self.predictor(runner).run()

    def test_profile_family_disagreement_is_rejected(self):
        def replace_family(fixture):
            fixture["predictions"][0]["profile_families"] = ["not-cas3"]

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=replace_family
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "families disagree"):
            self.predictor(runner).run()

    def test_protein_decision_is_recomputed(self):
        def corrupt_decision(fixture):
            prediction = fixture["predictions"][0]
            prediction["is_cas"] = False
            prediction["profile_families"] = []
            prediction["classification"] = {
                "class": None,
                "type": None,
                "subtype": None,
            }
            prediction["evidence"]["raw_classification"] = dict(
                prediction["classification"]
            )
            prediction["evidence"]["hierarchy_conflict_masked"] = False

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_decision
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "frozen protein decision"):
            self.predictor(runner).run()

    def test_protein_decision_constants_are_pinned(self):
        def corrupt_threshold(fixture):
            fixture["predictions"][0]["evidence"][
                "decision_threshold"
            ] = -1.0

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_threshold
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "protein decision threshold"):
            self.predictor(runner).run()

    def test_frozen_subtype_enum_is_enforced(self):
        def invent_subtype(fixture):
            fixture["predictions"][0]["classification"]["subtype"] = (
                "I-BOGUS"
            )

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=invent_subtype
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "subtype is outside"):
            self.predictor(runner).run()

    def test_class_type_mapping_is_enforced(self):
        def mismatch_class(fixture):
            fixture["predictions"][0]["classification"]["class"] = "2"

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=mismatch_class
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "class and type are inconsistent"):
            self.predictor(runner).run()

    def test_inconsistent_raw_subtype_is_masked_and_verified(self):
        def mask_conflicting_subtype(fixture):
            prediction = fixture["predictions"][3]
            prediction["evidence"]["raw_classification"] = {
                "class": "1",
                "type": "I",
                "subtype": "V-M",
            }
            prediction["evidence"]["hierarchy_conflict_masked"] = True

        result = self.predictor(FixtureRunner(
            self.fasta, self.installation, mutation=mask_conflicting_subtype
        )).run()
        self.assertEqual("computed", result["status"])
        self.assertIsNone(
            result["genes"][2]["classification"]["subtype"]
        )
        masked_gene = result["genes"][2]
        self.assertEqual(
            {"class": "1", "type": "I", "subtype": "V-M"},
            masked_gene["classification_evidence"]["raw_classification"],
        )
        self.assertTrue(
            masked_gene["classification_evidence"][
                "hierarchy_conflict_masked"
            ]
        )
        self.assertEqual(
            1,
            result["run_statistics"]["hierarchy_conflict_masked_proteins"],
        )

    def test_emitted_classification_must_match_raw_normalization(self):
        def lie_about_normalization(fixture):
            prediction = fixture["predictions"][0]
            prediction["evidence"]["raw_classification"] = {
                "class": "1",
                "type": "I",
                "subtype": "V-M",
            }
            prediction["evidence"]["hierarchy_conflict_masked"] = True

        with self.assertRaisesRegex(
                casandra.ExternalToolError, "hierarchy-normalized raw"):
            self.predictor(FixtureRunner(
                self.fasta,
                self.installation,
                mutation=lie_about_normalization,
            )).run()

    def test_hierarchy_conflict_flag_is_recomputed(self):
        def lie_about_conflict_flag(fixture):
            fixture["predictions"][0]["evidence"][
                "hierarchy_conflict_masked"
            ] = True

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "hierarchy_conflict_masked disagrees"):
            self.predictor(FixtureRunner(
                self.fasta,
                self.installation,
                mutation=lie_about_conflict_flag,
            )).run()

    def test_raw_classification_is_required(self):
        def omit_raw_classification(fixture):
            del fixture["predictions"][0]["evidence"]["raw_classification"]

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "raw prediction.*classification is malformed"):
            self.predictor(FixtureRunner(
                self.fasta,
                self.installation,
                mutation=omit_raw_classification,
            )).run()

    def test_raw_classification_requires_exact_fields(self):
        def omit_raw_field(fixture):
            del fixture["predictions"][0]["evidence"]["raw_classification"][
                "subtype"
            ]

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "raw prediction.*classification is malformed"):
            self.predictor(FixtureRunner(
                self.fasta,
                self.installation,
                mutation=omit_raw_field,
            )).run()

    def test_raw_classification_enum_is_enforced(self):
        def invent_raw_subtype(fixture):
            fixture["predictions"][0]["evidence"]["raw_classification"][
                "subtype"
            ] = "I-BOGUS"

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "raw prediction.*subtype is outside"):
            self.predictor(FixtureRunner(
                self.fasta,
                self.installation,
                mutation=invent_raw_subtype,
            )).run()

    def test_negative_prediction_requires_null_raw_classification(self):
        def add_negative_raw_claim(fixture):
            prediction = fixture["predictions"][1]
            prediction["evidence"]["raw_classification"] = {
                "class": None,
                "type": None,
                "subtype": "I-A",
            }
            prediction["evidence"]["hierarchy_conflict_masked"] = True

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "negative prediction.*raw Cas classification"):
            self.predictor(FixtureRunner(
                self.fasta,
                self.installation,
                mutation=add_negative_raw_claim,
            )).run()

    def test_type_vii_without_subtype_is_consistent_and_experimental(self):
        def make_type_vii(fixture):
            for index in (0, 2):
                classification = fixture["predictions"][index]["classification"]
                classification["type"] = "VII"
                classification["subtype"] = None
                fixture["predictions"][index]["evidence"][
                    "raw_classification"
                ] = dict(classification)
                fixture["predictions"][index]["evidence"][
                    "hierarchy_conflict_masked"
                ] = False
            direct = {
                "abstained_subtype": False,
                "class": "1",
                "confidence": {
                    "class_vote_fraction": 1.0,
                    "subtype_vote_fraction": 0.0,
                    "type_vote_fraction": 1.0,
                },
                "raw_subtype_winner": None,
                "subtype": None,
                "subtype_supporting_proteins": [],
                "type": "VII",
                "votes": {
                    "class": {"1": 24.0},
                    "subtype": {},
                    "type": {"VII": 24.0},
                },
            }
            cassette = fixture["cassettes"][0]
            cassette["classification"] = direct
            cassette["final_classification"] = {
                "class": "1",
                "confidence": 1.0,
                "direct_profile_result": json.loads(json.dumps(direct)),
                "evidence": {
                    "raw_classification": {
                        "class": "1", "type": "VII", "subtype": None,
                    },
                    "hierarchy_conflict_masked": False,
                },
                "method": "direct_profile_type_aggregation",
                "subtype": None,
                "type": "VII",
            }

        result = self.predictor(FixtureRunner(
            self.fasta, self.installation, mutation=make_type_vii
        )).run()
        system = result["systems"][0]
        self.assertEqual("CAS-VII", system["type"])
        self.assertIsNone(system["classification"]["subtype"])
        self.assertEqual(
            "positive_bitscore_weighted_type_vote_fraction",
            system["classification"]["confidence_kind"],
        )
        self.assertTrue(system["experimental"])
        self.assertEqual(
            casandra.TYPE_VII_EXPERIMENTAL_REASON,
            system["experimental_reason"],
        )
        self.assertEqual(
            system["experimental_reason"],
            system["evidence"]["experimental_reason"],
        )
        self.assertTrue(system["evidence"]["type_vii_experimental"])

    def test_type_vii_final_claim_without_member_evidence_is_rejected(self):
        def fake_type_vii(fixture):
            final = fixture["cassettes"][0]["final_classification"]
            final.update({
                "type": "VII",
                "subtype": None,
                "method": "direct_profile_type_aggregation",
            })

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=fake_type_vii
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "hierarchy-normalized raw route"):
            self.predictor(runner).run()

    def test_weighted_votes_are_recomputed(self):
        def corrupt_votes(fixture):
            fixture["cassettes"][0]["classification"]["votes"]["type"][
                "I"
            ] = 23.0

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_votes
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "recomputed frozen-v4 semantics"):
            self.predictor(runner).run()

    def test_direct_route_confidence_is_recomputed(self):
        def corrupt_confidence(fixture):
            fixture["cassettes"][0]["final_classification"][
                "confidence"
            ] = 0.75

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_confidence
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "direct-route confidence"):
            self.predictor(runner).run()

    def test_direct_route_method_is_enforced(self):
        def corrupt_method(fixture):
            fixture["cassettes"][0]["final_classification"]["method"] = (
                "type_iii_ordered_profile_architecture_extratrees"
            )

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_method
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "direct cassette route"):
            self.predictor(runner).run()

    def test_genome_gate_rule_is_recomputed(self):
        def corrupt_gate(fixture):
            fixture["cassettes"][0]["genome_evidence_gate"]["rule"] = (
                "class2_multi_gene"
            )

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_gate
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "genome evidence gate"):
            self.predictor(runner).run()

    def test_every_positive_gene_must_belong_to_a_candidate(self):
        def omit_rejected_candidate(fixture):
            fixture["rejected"].clear()

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=omit_rejected_candidate
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "do not partition all"):
            self.predictor(runner).run()

    def test_reported_protein_must_match_table_11_translation(self):
        def corrupt_translation(fixture):
            protein = fixture["proteins"][0]
            protein["protein_sequence"] = "MM"
            protein["protein_sha256"] = digest_bytes(b"MM")

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=corrupt_translation
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "table-11 translation"):
            self.predictor(runner).run()

    def test_vendored_source_tree_matches_adapter_pin(self):
        package_root = (
            ROOT
            / "vendor"
            / "casandra"
            / "src"
            / "casandra"
        )
        self.assertEqual(
            casandra.CASANDRA_SOURCE_TREE_SHA256,
            casandra._casandra_source_aggregate(package_root),
        )

    def test_tampered_source_tree_fails_preflight_integrity(self):
        source_root = (
            ROOT
            / "vendor"
            / "casandra"
            / "src"
            / "casandra"
        )
        tampered_root = self.root / "tampered-casandra"
        tampered_root.mkdir()
        for source_path in source_root.glob("*.py"):
            shutil.copy2(source_path, tampered_root / source_path.name)
        init_path = tampered_root / "__init__.py"
        init_path.write_text(
            init_path.read_text(encoding="utf-8") + "\n# tampered\n",
            encoding="utf-8",
        )
        distribution = types.SimpleNamespace(
            version=casandra.CASANDRA_VERSION,
            locate_file=lambda relative: tampered_root,
        )
        package_spec = types.SimpleNamespace(
            origin=str(tampered_root / "__init__.py")
        )
        with (
            mock.patch.object(
                casandra.importlib.metadata,
                "distribution",
                return_value=distribution,
            ),
            mock.patch.object(
                casandra.importlib.util,
                "find_spec",
                return_value=package_spec,
            ),
            self.assertRaisesRegex(
                casandra.ExternalToolError, "source does not match"
            ),
        ):
            casandra.inspect_casandra_installation()

    def test_timeout_is_fail_closed_without_legacy_fallback(self):
        def time_out(command, **kwargs):
            raise casandra.BoundedSubprocessTimeoutError("deliberate timeout")

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "exceeded its .* timeout; no legacy"):
            self.predictor(time_out).run()

    def test_stdout_truncation_is_fail_closed_without_legacy_fallback(self):
        def truncate(command, **kwargs):
            raise casandra.BoundedSubprocessTruncationError(
                "deliberate truncation"
            )

        with self.assertRaisesRegex(
                casandra.ExternalToolError,
                "stdout exceeded.*no legacy"):
            self.predictor(truncate).run()

    def test_timeout_must_be_positive_finite_non_bool(self):
        for timeout in (0, -1, float("inf"), True, "3600"):
            with self.subTest(timeout=timeout):
                with self.assertRaisesRegex(ValueError, "timeout_seconds"):
                    casandra.CasAndraPredictor(
                        input_fasta=self.fasta,
                        source_id="record",
                        source_dna=self.dna,
                        installation=self.installation,
                        timeout_seconds=timeout,
                    )

    def test_nonzero_exit_is_fail_closed_without_legacy_fallback(self):
        runner = FixtureRunner(
            self.fasta, self.installation, returncode=9
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "no legacy Cas backend was used"):
            self.predictor(runner).run()
        self.assertEqual(1, len(runner.commands))

    def test_tampered_manifested_output_is_rejected(self):
        def tamper(fixture):
            fixture["tamper_run_after_manifest"] = True

        runner = FixtureRunner(self.fasta, self.installation, mutation=tamper)
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "manifest validation"):
            self.predictor(runner).run()

    def test_source_coordinate_mismatch_is_rejected(self):
        def shift_first_gene(fixture):
            fixture["proteins"][0]["start_1based"] = 2

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=shift_first_gene
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "do not reconstruct"):
            self.predictor(runner).run()

    def test_wrong_provider_bundle_is_rejected(self):
        def replace_bundle(fixture):
            fixture["run"]["bundle_id"] = "unreviewed-model"

        runner = FixtureRunner(
            self.fasta, self.installation, mutation=replace_bundle
        )
        with self.assertRaisesRegex(
                casandra.ExternalToolError, "bundle_id"):
            self.predictor(runner).run()

    def test_preflight_error_is_actionable(self):
        with mock.patch.object(
                casandra,
                "inspect_casandra_installation",
                side_effect=casandra.ExternalToolError("model hash mismatch"),
        ):
            error = casandra.casandra_installation_error()
        self.assertIn("--cas True requires", error)
        self.assertIn("scripts/setup_local.sh", error)


if __name__ == "__main__":
    unittest.main()

