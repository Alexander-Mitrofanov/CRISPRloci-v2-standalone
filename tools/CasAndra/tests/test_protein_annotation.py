import hashlib
import json
from types import SimpleNamespace

import pytest

from casandra.protein_annotation import annotate_protein_fasta


def _fake_predict(_model, _input_path, output_path, *, threads):
    rows = [
        {
            "sequence_id": "cas2_without_type",
            "sequence_length": 4,
            "is_cas": True,
            "best_positive_profile": "C25_Cas2_1",
            "cas_family": "Cas2",
            "result": "Cas2",
            "classification": {"class": None, "type": None, "subtype": None},
            "positive_profile_score": 42.0,
            "hard_negative_profile_score": 1.0,
            "score_margin": 41.0,
            "evidence": {"model_id": "test-model"},
        },
        {
            "sequence_id": "negative",
            "sequence_length": 4,
            "is_cas": False,
            "best_positive_profile": None,
            "cas_family": None,
            "result": "no cas",
            "classification": {"class": None, "type": None, "subtype": None},
            "positive_profile_score": -10.0,
            "hard_negative_profile_score": 2.0,
            "score_margin": -12.0,
            "evidence": {"model_id": "test-model"},
        },
    ]
    output_path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )
    return {
        "schema_version": 1,
        "model_id": "test-model",
        "input_records": 2,
        "cas_predictions": 1,
        "cpu_threads": threads,
    }


def test_annotation_output_is_atomic_provenance_bound_and_checksummed(
    tmp_path, monkeypatch
):
    monkeypatch.setattr("casandra.protein_annotation.predict_proteins", _fake_predict)
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "manifest.json").write_text('{"trusted":true}\n', encoding="utf-8")
    input_path = tmp_path / "proteins.faa"
    input_path.write_text(
        ">cas2_without_type\nMKTW\n>negative\nACDX\n", encoding="ascii"
    )
    output = tmp_path / "annotation"
    model = SimpleNamespace(
        metadata={"C25_Cas2_1": {"families": ["cas2", "cas"]}}
    )

    run = annotate_protein_fasta(
        bundle,
        input_path,
        output,
        protein_model=model,
        bundle_manifest={"role": "test-role"},
        bundle_config={"bundle_id": "test-bundle"},
        threads=2,
    )

    assert sorted(path.name for path in output.iterdir()) == [
        "manifest.json",
        "protein_predictions.jsonl",
        "run.json",
    ]
    assert run["analysis"] == "annotate_cas_genes"
    assert run["program_version"] == "0.3.0.dev2"
    assert run["bundle_id"] == "test-bundle"
    assert run["model_id"] == "test-model"
    assert run["protein_records"] == 2
    assert run["cas_proteins"] == 1
    assert run["result_contract"] == {"positive": "cas_family", "negative": "no cas"}
    assert run["wall_seconds"] >= 0
    assert run["inputs"][0]["sha256"] == hashlib.sha256(input_path.read_bytes()).hexdigest()
    prediction_rows = [
        json.loads(line)
        for line in (output / "protein_predictions.jsonl").read_text().splitlines()
    ]
    assert [row["result"] for row in prediction_rows] == ["Cas2", "no cas"]

    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["bundle_manifest_sha256"] == hashlib.sha256(
        (bundle / "manifest.json").read_bytes()
    ).hexdigest()
    assert set(manifest["files"]) == {"protein_predictions.jsonl", "run.json"}
    for name, record in manifest["files"].items():
        path = output / name
        assert record == {
            "size": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }


def test_annotation_rejects_positive_without_curated_family(tmp_path, monkeypatch):
    def invalid_predict(model, input_path, output_path, *, threads):
        summary = _fake_predict(model, input_path, output_path, threads=threads)
        rows = [json.loads(line) for line in output_path.read_text().splitlines()]
        rows[0]["cas_family"] = None
        output_path.write_text(
            "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
        )
        return summary

    monkeypatch.setattr("casandra.protein_annotation.predict_proteins", invalid_predict)
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "manifest.json").write_text("{}\n", encoding="utf-8")
    input_path = tmp_path / "proteins.faa"
    input_path.write_text(
        ">cas2_without_type\nMKTW\n>negative\nACDX\n", encoding="ascii"
    )
    output = tmp_path / "annotation"
    model = SimpleNamespace(metadata={"C25_Cas2_1": {"families": ["cas2"]}})

    with pytest.raises(RuntimeError, match="curated Cas-family"):
        annotate_protein_fasta(
            bundle,
            input_path,
            output,
            protein_model=model,
            bundle_manifest={},
            bundle_config={"bundle_id": "test-bundle"},
            threads=1,
        )
    assert not output.exists()


@pytest.mark.parametrize("malformation", ["missing", "contradictory"])
def test_annotation_rejects_missing_or_contradictory_result(
    tmp_path, monkeypatch, malformation
):
    def invalid_predict(model, input_path, output_path, *, threads):
        summary = _fake_predict(model, input_path, output_path, threads=threads)
        rows = [json.loads(line) for line in output_path.read_text().splitlines()]
        if malformation == "missing":
            rows[0].pop("result")
        else:
            rows[0]["result"] = "no cas"
        output_path.write_text(
            "".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8"
        )
        return summary

    monkeypatch.setattr("casandra.protein_annotation.predict_proteins", invalid_predict)
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "manifest.json").write_text("{}\n", encoding="utf-8")
    input_path = tmp_path / "proteins.faa"
    input_path.write_text(
        ">cas2_without_type\nMKTW\n>negative\nACDX\n", encoding="ascii"
    )
    output = tmp_path / "annotation"
    model = SimpleNamespace(metadata={"C25_Cas2_1": {"families": ["cas2"]}})

    with pytest.raises(RuntimeError, match="result disagrees"):
        annotate_protein_fasta(
            bundle,
            input_path,
            output,
            protein_model=model,
            bundle_manifest={},
            bundle_config={"bundle_id": "test-bundle"},
            threads=1,
        )
    assert not output.exists()
