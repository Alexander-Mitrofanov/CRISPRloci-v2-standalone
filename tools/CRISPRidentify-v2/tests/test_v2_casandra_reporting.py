import copy
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
V2 = ROOT

TYPE_VII_REASON = (
    "CasAndra v4 Type VII calls lack independent whole-genome cassette validation"
)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cas_output = _load(
    "v2_cas_output_reporting_contract",
    V2 / "components" / "cas_output.py",
)
cas_html = _load(
    "v2_cas_html_reporting_contract",
    V2 / "components" / "cas_html.py",
)


def _gene(
        protein_id, start, end, *, is_cas, annotation="unknown",
        profile=None, score=None):
    classification = {
        "class": "1" if is_cas else None,
        "type": "I" if is_cas else None,
        "subtype": "I-F" if is_cas else None,
    }
    return {
        "protein_id": protein_id,
        "source_id": "record",
        "start": start,
        "end": end,
        "strand": "+",
        "is_cas": is_cas,
        "annotation": annotation,
        "classification": dict(classification),
        "classification_evidence": {
            "raw_classification": dict(classification),
            "hierarchy_conflict_masked": False,
        },
        "profile": profile,
        "score": score,
        "score_kind": (
            "CasAndra positive-profile score" if score is not None else None
        ),
    }


def _casandra_result():
    first = _gene(
        "casandra|protein=1", 100, 150, is_cas=True,
        annotation="cas7", profile="profile_7", score=42.5,
    )
    second = _gene(
        "casandra|protein=2", 210, 300, is_cas=True,
        annotation="cas5", profile="profile_5", score=31.25,
    )
    orphan = _gene(
        "casandra|protein=orphan", 400, 450, is_cas=True,
        annotation="cas1", profile="profile_1", score=19.0,
    )
    orphan["classification"]["subtype"] = None
    orphan["classification_evidence"] = {
        "raw_classification": {
            "class": "1", "type": "I", "subtype": "V-M",
        },
        "hierarchy_conflict_masked": True,
    }
    context = _gene(
        "casandra|protein=context", 160, 200, is_cas=False,
    )
    return {
        "status": "computed",
        "coordinate_system": {
            "end_inclusive": True,
            "indexing": 1,
            "reference": "input_forward_strand",
        },
        "method": {
            "provider": "CasAndra",
            "provider_version": "0.1.0.dev0+crispridentify.3",
            "bundle_id": "casandra-test-bundle",
            "bundle_role": "deployment refit",
            "pipeline_schema_version": 4,
            "gene_caller": "Pyrodigal 3.7.1",
            "gene_mode": "single",
            "translation_table": 11,
            "model_manifest_sha256": "a" * 64,
        },
        "gene_count": 3,
        "genes": [first, second, orphan],
        "protein_count": 12,
        "run_statistics": {
            "hierarchy_conflict_masked_cassettes": 0,
            "raw_cassette_candidates": 1,
            "rejected_cassette_candidates": 0,
            "hierarchy_conflict_masked_proteins": 1,
        },
        "system_count": 1,
        "systems": [{
            "cassette_id": 1,
            "provider_cassette_id": (
                "casandra|contig=record|cassette=<1&>|loc=100-300"
            ),
            "source_id": "record",
            "start": 100,
            "end": 300,
            "type": "CAS-I-F",
            "classification": {
                "class": "1",
                "type": "I",
                "subtype": "I-F",
                "method": "direct_profile_aggregation",
                "confidence": 0.75,
                "confidence_kind": "provider_reported",
            },
            "evidence": {
                "genome_evidence_gate": "class1_requires_two_cas_genes",
                "provider_direct_profile_result": {
                    "class": "1", "type": "I", "subtype": "I-F",
                },
                "provider_final_classification_evidence": {
                    "raw_classification": {
                        "class": "1", "type": "I", "subtype": "I-F",
                    },
                    "hierarchy_conflict_masked": False,
                },
                "type_vii_experimental": False,
                "experimental_reason": None,
            },
            "experimental": False,
            "experimental_reason": None,
            "genes": [first, context, second],
        }],
    }


def _type_vii_result():
    payload = _casandra_result()
    system = payload["systems"][0]
    system["type"] = "CAS-VII"
    system["classification"].update({
        "type": "VII",
        "subtype": None,
        "method": "direct_profile_type_aggregation",
        "confidence_kind": (
            "positive_bitscore_weighted_type_vote_fraction"
        ),
    })
    system.update({
        "experimental": True,
        "experimental_reason": TYPE_VII_REASON,
    })
    system["evidence"].update({
        "genome_evidence_gate": "class1_requires_two_cas_genes",
        "provider_direct_profile_result": {
            "class": "1", "type": "VII", "subtype": None,
        },
        "provider_final_classification_evidence": {
            "raw_classification": {
                "class": "1", "type": "VII", "subtype": None,
            },
            "hierarchy_conflict_masked": False,
        },
        "type_vii_experimental": True,
        "experimental_reason": TYPE_VII_REASON,
    })
    return payload


def _mixed_type_vi_result():
    payload = _casandra_result()
    system = payload["systems"][0]
    system["type"] = "CAS-VI-B"
    system["classification"].update({
        "class": "2",
        "type": "VI",
        "subtype": "VI-B",
    })
    system["evidence"].update({
        "genome_evidence_gate": "class2_multi_gene",
        "provider_direct_profile_result": {
            "class": "1", "type": "VI", "subtype": "VI-B",
        },
        "provider_final_classification_evidence": {
            "raw_classification": {
                "class": "1", "type": "VI", "subtype": "VI-B",
            },
            "hierarchy_conflict_masked": True,
        },
    })
    payload["run_statistics"][
        "hierarchy_conflict_masked_cassettes"
    ] = 1
    return payload


def _model():
    return {
        "source": {"id": "record", "length": 1_000},
        "arrays": [{
            "id": "CRISPR-20-60-bf1",
            "source_interval": {"start": 20, "end": 60},
        }],
        "cas": _casandra_result(),
    }


class CasAndraProjectionTests(unittest.TestCase):
    def test_normalization_preserves_provider_claims_and_display_id(self):
        normalized = cas_output.normalize_cas_data(
            _casandra_result(), source_id="record", source_length=1_000
        )

        system = normalized["systems"][0]
        self.assertEqual(1, system["cassette_id"])
        self.assertEqual(
            "casandra|contig=record|cassette=<1&>|loc=100-300",
            system["provider_cassette_id"],
        )
        self.assertEqual("CAS-I-F", system["type"])
        self.assertEqual("I-F", system["classification"]["subtype"])
        self.assertEqual("provider_reported", system["confidence_kind"])
        self.assertEqual(1, system["missing_gene_count"])
        self.assertNotIn("probabilities", system)
        self.assertEqual(3, normalized["gene_count"])
        orphan = next(
            gene for gene in normalized["genes"]
            if gene["protein_id"] == "casandra|protein=orphan"
        )
        self.assertEqual(
            {"class": "1", "type": "I", "subtype": "V-M"},
            orphan["classification_evidence"]["raw_classification"],
        )
        self.assertTrue(
            orphan["classification_evidence"]["hierarchy_conflict_masked"]
        )
        self.assertEqual(
            1,
            normalized["run_statistics"][
                "hierarchy_conflict_masked_proteins"
            ],
        )

    def test_mixed_cassette_evidence_and_qc_count_are_preserved(self):
        normalized = cas_output.normalize_cas_data(
            _mixed_type_vi_result(),
            source_id="record",
            source_length=1_000,
        )

        system = normalized["systems"][0]
        self.assertEqual("CAS-VI-B", system["type"])
        self.assertEqual("2", system["classification"]["class"])
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
            1,
            normalized["run_statistics"][
                "hierarchy_conflict_masked_cassettes"
            ],
        )

    def test_report_replays_raw_final_hierarchy_normalization(self):
        payload = _mixed_type_vi_result()
        payload["systems"][0]["classification"]["class"] = "1"

        with self.assertRaisesRegex(
                cas_output.CasReportError,
                "hierarchy-normalized raw provider claim"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_report_recomputes_cassette_conflict_flag(self):
        payload = _mixed_type_vi_result()
        payload["systems"][0]["evidence"][
            "provider_final_classification_evidence"
        ]["hierarchy_conflict_masked"] = False

        with self.assertRaisesRegex(
                cas_output.CasReportError,
                "hierarchy_conflict_masked disagrees"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_report_links_raw_final_claim_to_direct_evidence(self):
        payload = _mixed_type_vi_result()
        evidence = payload["systems"][0]["evidence"][
            "provider_final_classification_evidence"
        ]
        evidence["raw_classification"]["class"] = "2"
        evidence["hierarchy_conflict_masked"] = False

        with self.assertRaisesRegex(
                cas_output.CasReportError,
                "raw final class/type disagrees with direct"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_masked_cassette_count_must_match_system_evidence(self):
        payload = _mixed_type_vi_result()
        payload["run_statistics"][
            "hierarchy_conflict_masked_cassettes"
        ] = 0

        with self.assertRaisesRegex(
                cas_output.CasReportError,
                "masked cassette count disagrees"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_inconsistent_positive_context_gene_fails_closed(self):
        payload = _casandra_result()
        payload["systems"][0]["genes"][0] = dict(
            payload["systems"][0]["genes"][0]
        )
        payload["systems"][0]["genes"][0]["score"] = 99.0

        with self.assertRaisesRegex(
                cas_output.CasReportError, "disagrees with top-level gene"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_classification_evidence_shape_is_strict(self):
        payload = _casandra_result()
        payload["genes"][2]["classification_evidence"][
            "hierarchy_conflict_masked"
        ] = "yes"

        with self.assertRaisesRegex(
                cas_output.CasReportError,
                "hierarchy_conflict_masked must be a boolean"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_hierarchy_conflict_run_count_is_recomputed(self):
        payload = _casandra_result()
        payload["run_statistics"][
            "hierarchy_conflict_masked_proteins"
        ] = 0

        with self.assertRaisesRegex(
                cas_output.CasReportError, "hierarchy conflict count"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_context_classification_evidence_must_match_top_gene(self):
        payload = _casandra_result()
        context_gene = copy.deepcopy(payload["systems"][0]["genes"][0])
        context_gene["classification_evidence"]["raw_classification"][
            "subtype"
        ] = "V-F"
        context_gene["classification_evidence"][
            "hierarchy_conflict_masked"
        ] = True
        payload["systems"][0]["genes"][0] = context_gene

        with self.assertRaisesRegex(
                cas_output.CasReportError, "disagrees with top-level gene"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_type_vii_without_subtype_is_preserved(self):
        payload = _type_vii_result()
        normalized = cas_output.normalize_cas_data(
            payload, source_id="record", source_length=1_000
        )
        self.assertEqual("CAS-VII", normalized["systems"][0]["type"])
        self.assertIsNone(
            normalized["systems"][0]["classification"]["subtype"]
        )
        system = normalized["systems"][0]
        self.assertTrue(system["experimental"])
        self.assertEqual(TYPE_VII_REASON, system["experimental_reason"])
        self.assertTrue(system["evidence"]["type_vii_experimental"])
        self.assertEqual(
            TYPE_VII_REASON, system["evidence"]["experimental_reason"]
        )

    def test_experimental_flag_must_match_type_vii_classification(self):
        payload = _casandra_result()
        system = payload["systems"][0]
        system["experimental"] = True
        system["experimental_reason"] = TYPE_VII_REASON
        system["evidence"]["type_vii_experimental"] = True
        system["evidence"]["experimental_reason"] = TYPE_VII_REASON

        with self.assertRaisesRegex(
                cas_output.CasReportError, "disagrees with classification"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_experimental_evidence_must_match_system(self):
        payload = _type_vii_result()
        payload["systems"][0]["evidence"]["experimental_reason"] = "different"

        with self.assertRaisesRegex(
                cas_output.CasReportError, "evidence experimental reason"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )

    def test_native_ids_must_be_unique_and_display_ids_sequential(self):
        payload = _casandra_result()
        duplicate = copy.deepcopy(payload["systems"][0])
        duplicate["cassette_id"] = 3
        duplicate["start"] = 500
        duplicate["end"] = 600
        duplicate["provider_cassette_id"] = "second"
        for gene in duplicate["genes"]:
            gene["start"] += 400
            gene["end"] += 300
        payload["systems"].append(duplicate)
        payload["system_count"] = 2
        payload["run_statistics"]["raw_cassette_candidates"] = 2

        with self.assertRaisesRegex(
                cas_output.CasReportError, "display IDs must be sequential"):
            cas_output.normalize_cas_data(
                payload, source_id="record", source_length=1_000
            )


class CasAndraArtifactTests(unittest.TestCase):
    def test_csvs_are_provider_neutral_and_keep_orphan_positives(self):
        with tempfile.TemporaryDirectory() as temporary:
            cas_output.write_cas_outputs(
                temporary,
                _casandra_result(),
                source_id="record",
                source_length=1_000,
            )
            with (Path(temporary) / "Cas_systems.csv").open(
                    newline="", encoding="utf-8") as stream:
                systems = list(csv.DictReader(stream))
            with (Path(temporary) / "Cas_genes.csv").open(
                    newline="", encoding="utf-8") as stream:
                genes = list(csv.DictReader(stream))

        self.assertEqual("CAS-I-F", systems[0]["type"])
        self.assertEqual("direct_profile_aggregation",
                         systems[0]["classification_method"])
        self.assertEqual("provider_reported", systems[0]["confidence_kind"])
        self.assertEqual("1", systems[0]["non_cas_context_gene_count"])
        self.assertEqual("False", systems[0]["experimental"])
        self.assertEqual("", systems[0]["experimental_reason"])
        self.assertNotIn("top_probability", systems[0])
        orphan = next(
            row for row in genes
            if row["protein_id"] == "casandra|protein=orphan"
        )
        self.assertEqual("", orphan["cassette_id"])
        self.assertEqual("True", orphan["is_cas"])
        self.assertNotIn("bitscore", orphan)

    def test_type_vii_csv_exposes_experimental_reason(self):
        with tempfile.TemporaryDirectory() as temporary:
            cas_output.write_cas_outputs(
                temporary,
                _type_vii_result(),
                source_id="record",
                source_length=1_000,
            )
            with (Path(temporary) / "Cas_systems.csv").open(
                    newline="", encoding="utf-8") as stream:
                system = next(csv.DictReader(stream))

        self.assertEqual("True", system["experimental"])
        self.assertEqual(TYPE_VII_REASON, system["experimental_reason"])

    def test_malformed_projection_does_not_publish_csvs(self):
        payload = _casandra_result()
        payload["systems"][0]["end"] = 2_000
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(
                    cas_output.CasReportError, "outside the input sequence"):
                cas_output.write_cas_outputs(
                    temporary,
                    payload,
                    source_id="record",
                    source_length=1_000,
                )
            self.assertFalse((Path(temporary) / "Cas_systems.csv").exists())
            self.assertFalse((Path(temporary) / "Cas_genes.csv").exists())

    def test_html_names_provider_and_score_semantics_without_legacy_claims(self):
        rendered = cas_html.render_cas_section(_model())

        self.assertIn("CasAndra 0.1.0.dev0+crispridentify.3", rendered)
        self.assertIn("Provider-reported confidence: 0.750", rendered)
        self.assertIn("provider_reported", rendered)
        self.assertIn("Non-Cas context genes", rendered)
        self.assertIn("CasAndra positive-profile score", rendered)
        self.assertIn(
            "casandra|contig=record|cassette=&lt;1&amp;&gt;|loc=100-300",
            rendered,
        )
        self.assertNotIn("HMM2019", rendered)
        self.assertNotIn("ExtraTrees", rendered)
        self.assertNotIn("Top vote share", rendered)


    def test_html_visibly_labels_experimental_type_vii(self):
        model = _model()
        model["cas"] = _type_vii_result()
        rendered = cas_html.render_cas_section(model)

        self.assertIn("Experimental Type VII call.", rendered)
        self.assertIn(TYPE_VII_REASON, rendered)

if __name__ == "__main__":
    unittest.main()

