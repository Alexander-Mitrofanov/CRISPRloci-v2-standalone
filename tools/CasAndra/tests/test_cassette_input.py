import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from casandra.api import DEFAULT_MODEL_BUNDLE
from casandra.architecture_runtime import load_architecture_model
from casandra.cassette_input import (
    classify_ordered_predictions,
    classify_protein_cassette,
)
from casandra.model import ProteinModel
from casandra.prediction import predict_query_records
from casandra.runtime import parse_fasta


class UnexpectedArchitectureModel:
    def classify(self, *_args, **_kwargs):
        raise AssertionError("Type I cassette must not use the Type III architecture model")


class RecordingArchitectureModel:
    def __init__(self):
        self.ordered_genes = []

    def classify(self, cassette, genes_by_id, _predictions_by_id):
        self.ordered_genes = [
            (protein_id, genes_by_id[protein_id]["start_1based"])
            for protein_id in cassette["cas_protein_ids"]
        ]
        return "III-A", 0.91


def _prediction(sequence_id, *, is_cas, cas_type=None, subtype=None, score=30.0):
    return {
        "sequence_id": sequence_id,
        "is_cas": is_cas,
        "positive_profile_score": score,
        "best_positive_profile": "Cas3" if is_cas else None,
        "classification": {
            "class": "1" if is_cas else None,
            "type": cas_type if is_cas else None,
            "subtype": subtype if is_cas else None,
        },
    }


def test_supplied_cassette_preserves_order_and_reports_non_cas_records():
    proteins = [("cas3", "MAAA"), ("other", "MDDD"), ("cas8", "MCCC")]
    predictions = [
        _prediction("cas3", is_cas=True, cas_type="I", subtype="I-E", score=50.0),
        _prediction("other", is_cas=False),
        _prediction("cas8", is_cas=True, cas_type="I", subtype="I-E", score=40.0),
    ]

    result = classify_ordered_predictions(
        proteins, predictions, UnexpectedArchitectureModel()
    )

    assert result["input_protein_ids"] == ["cas3", "other", "cas8"]
    assert result["cas_protein_ids"] == ["cas3", "cas8"]
    assert result["non_cas_protein_ids"] == ["other"]
    assert result["classification"]["type"] == "I"
    assert result["classification"]["subtype"] == "I-E"
    assert result["classification"]["method"] == "direct_profile_aggregation"
    assert result["coordinates_available"] is False


def test_supplied_cassette_rejects_prediction_reordering():
    proteins = [("first", "MAAA"), ("second", "MCCC")]
    predictions = [
        _prediction("second", is_cas=False),
        _prediction("first", is_cas=False),
    ]

    try:
        classify_ordered_predictions(
            proteins, predictions, UnexpectedArchitectureModel()
        )
    except ValueError as error:
        assert "record order" in str(error)
    else:
        raise AssertionError("reordered predictions were accepted")


def test_type_iii_cassette_routes_ordered_proteins_through_architecture_model():
    proteins = [("cas10", "MAAA"), ("cas7", "MCCC"), ("cas5", "MDDD")]
    predictions = [
        _prediction("cas10", is_cas=True, cas_type="III", subtype=None, score=80.0),
        _prediction("cas7", is_cas=True, cas_type="III", subtype=None, score=50.0),
        _prediction("cas5", is_cas=True, cas_type="III", subtype=None, score=40.0),
    ]
    model = RecordingArchitectureModel()

    result = classify_ordered_predictions(proteins, predictions, model)

    assert result["input_protein_ids"] == ["cas10", "cas7", "cas5"]
    assert model.ordered_genes == [("cas10", 1), ("cas7", 16), ("cas5", 31)]
    assert result["classification"]["type"] == "III"
    assert result["classification"]["subtype"] == "III-A"
    assert (
        result["classification"]["method"]
        == "type_iii_ordered_profile_architecture_extratrees"
    )
    assert result["classification"]["confidence"] == 0.91


def test_canonical_spyogenes_type_iia_cassette_uses_architecture_not_fixture_override():
    fixture = (
        Path(__file__).resolve().parent
        / "fixtures"
        / "canonical_type_iia_cassette.faa"
    )
    proteins = list(parse_fasta(fixture))
    protein_model = ProteinModel.from_directory(DEFAULT_MODEL_BUNDLE / "protein")
    predictions = predict_query_records(
        protein_model,
        [
            {"sequence_id": sequence_id, "sequence": sequence}
            for sequence_id, sequence in proteins
        ],
        threads=1,
    )
    architecture_model = load_architecture_model(
        DEFAULT_MODEL_BUNDLE / "cassette_architecture"
    )

    result = classify_ordered_predictions(proteins, predictions, architecture_model)

    assert result["input_protein_ids"] == [
        "SPY_RS04360_cas9",
        "SPY_RS04365_cas1",
        "SPY_RS04370_cas2",
        "SPY_RS04375_csn2",
    ]
    assert result["cas_protein_ids"] == [
        "SPY_RS04360_cas9",
        "SPY_RS04365_cas1",
        "SPY_RS04370_cas2",
    ]
    assert result["non_cas_protein_ids"] == ["SPY_RS04375_csn2"]
    assert result["direct_profile_classification"]["subtype"] == "II-C"
    assert result["classification"]["type"] == "II"
    assert result["classification"]["subtype"] == "II-A"
    assert (
        result["classification"]["method"]
        == "type_ii_ordered_profile_architecture_extratrees"
    )
    assert result["classification"]["direct_profile_result"] == result[
        "direct_profile_classification"
    ]
    # The cap15 bundle includes the September full361 architecture refit.
    assert result["classification"]["confidence"] == 0.5682433452220029


def test_cassette_output_binds_run_and_manifest_to_same_bundle(tmp_path, monkeypatch):
    def fake_predict(_model, _input_path, output_path, *, threads):
        prediction = _prediction(
            "cas3", is_cas=True, cas_type="I", subtype="I-E", score=50.0
        )
        prediction.update(
            {
                "sequence_length": 4,
                "cas_family": "Cas3",
                "evidence": {"model_id": "test-model"},
            }
        )
        output_path.write_text(json.dumps(prediction) + "\n", encoding="utf-8")
        return {
            "model_id": "test-model",
            "input_records": 1,
            "cas_predictions": 1,
            "cpu_threads": threads,
        }

    monkeypatch.setattr("casandra.cassette_input.predict_proteins", fake_predict)
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    bundle_manifest_path = bundle / "manifest.json"
    bundle_manifest_path.write_text('{"trusted":true}\n', encoding="utf-8")
    input_path = tmp_path / "cassette.faa"
    input_path.write_text(">cas3\nMKTW\n", encoding="ascii")
    output = tmp_path / "result"

    run = classify_protein_cassette(
        bundle,
        input_path,
        output,
        protein_model=SimpleNamespace(),
        architecture_model=UnexpectedArchitectureModel(),
        bundle_manifest={"role": "test-role"},
        bundle_config={"bundle_id": "test-bundle"},
        threads=2,
    )

    expected_bundle_digest = hashlib.sha256(
        bundle_manifest_path.read_bytes()
    ).hexdigest()
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert run["bundle_manifest_sha256"] == expected_bundle_digest
    assert manifest["bundle_manifest_sha256"] == expected_bundle_digest
    assert run["program"] == "CasAndra"
    assert run["program_version"] == "0.3.0.dev2"
    assert run["offline_inference"] is True
    assert run["crispr_array_prediction"] is False
    assert set(manifest["files"]) == {"proteins.jsonl", "cassette.json", "run.json"}
