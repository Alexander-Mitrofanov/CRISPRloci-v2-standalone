import copy
import json
from importlib.resources import files
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
import xgboost as xgb

from repeattyper import Record, classify
from repeattyper.cli import main
from repeattyper.errors import InputError, ModelError, OutputError
from repeattyper.features import FEATURE_NAMES, matrix, reverse_complement
from repeattyper.model import load_model
from repeattyper.input import read_fasta
from repeattyper.output import publish

DIGEST = "0" * 64
HERE = Path(__file__).parent


def test_exact_upstream_scores_and_labels():
    spec = spec_from_file_location("upstream_xgb", HERE / "upstream/xgb.py")
    upstream = module_from_spec(spec)
    spec.loader.exec_module(upstream)
    records, _ = read_fasta(HERE / "upstream/repeats.fasta")
    import random
    rng = random.Random(42)
    records += [Record(f"random{i}", "".join(rng.choices("ACGT", k=length)))
                for i, length in enumerate([1, 2, 3, 4, 10, 20, 40, 60, 100, 200] * 10)]
    records += [Record(f"rc{i}", reverse_complement(row.sequence)) for i, row in enumerate(records)]
    data = files("repeattyper").joinpath("data")
    original = upstream.XGB(SimpleNamespace(threads=1, kmer=4,
        xgb=str(data.joinpath("xgb_repeats.model")), typedict=str(data.joinpath("type_dict.tab")),
        repeats=[row.sequence for row in records]))
    original.predict_repeats()
    booster, _ = load_model()
    scores = booster.predict(xgb.DMatrix(matrix([row.sequence for row in records]),
                                       feature_names=FEATURE_NAMES), iteration_range=(0, 41))
    np.testing.assert_array_equal(scores, original.z_pred)
    report = classify(records, input_sha256=DIGEST)
    assert [row["winner"] for row in report["results"]] == original.z_type
    np.testing.assert_array_equal([row["native_score"] for row in report["results"]], original.z_max)
    count = len(records) // 2
    for a, b in zip(report["results"][:count], report["results"][count:], strict=True):
        assert (a["winner"], a["native_score"]) == (b["winner"], b["native_score"])


def test_empty_ambiguous_duplicates_uncertainty():
    assert classify([], input_sha256=DIGEST)["results"] == []
    report = classify([Record("a", "ACGT"), Record("b", "ACGT"), Record("n", "ACGTN")],
                      input_sha256=DIGEST, threshold=1)
    assert report["counts"] == {"records": 3, "distinct_classified_sequences": 1,
                                "predicted": 0, "uncertain": 2, "unsupported": 1}
    assert report["results"][0]["prediction"] is None
    assert report["results"][0]["native_score"] == report["results"][1]["native_score"]
    assert report["results"][2]["winner"] is None
    assert report["results"][2]["native_score"] is None


@pytest.mark.parametrize("sequence", ["", "A" * 201, "ACTU", "<script>", "acgt"])
def test_invalid_api_sequences(sequence):
    with pytest.raises(InputError):
        classify([Record("a", sequence)], input_sha256=DIGEST)


@pytest.mark.parametrize("threshold", [float("nan"), float("inf"), -1, 2, True])
def test_invalid_threshold(threshold):
    with pytest.raises(InputError):
        classify([], input_sha256=DIGEST, threshold=threshold)


@pytest.mark.parametrize("text", [">a\nACGT\n>a\nACTG\n", ">a\n", "ACTG", ">\nACGT", ">a\nACGTU\n"])
def test_invalid_fasta(tmp_path, text):
    path = tmp_path / "input.fa"
    path.write_text(text)
    assert main(["classify", "--input", str(path), "--output", str(tmp_path / "out")]) == 2
    assert not (tmp_path / "out").exists()


def test_fasta_normalization_cli_atomic_publication(tmp_path):
    path = tmp_path / "input.fa"
    path.write_text(">a description\naCgT\n>other\nACGT\n")
    output = tmp_path / "out"
    assert main(["classify", "--input", str(path), "--output", str(output)]) == 0
    report = json.loads((output / "report.json").read_text())
    assert report["counts"]["records"] == 2
    assert report["results"][0]["sequence"] == "ACGT"
    assert (output / "predictions.tsv").read_text().count("\n") == 3
    original = (output / "report.json").read_bytes()
    assert main(["classify", "--input", str(path), "--output", str(output)]) == 4
    assert (output / "report.json").read_bytes() == original
    link = tmp_path / "link"
    link.symlink_to(path)
    assert main(["classify", "--input", str(link), "--output", str(tmp_path / "other")]) == 2


def test_model_tampering_rejected(monkeypatch):
    from repeattyper import model
    data = files("repeattyper").joinpath("data")
    provenance = json.loads(data.joinpath("provenance.json").read_text())
    corrupt = copy.deepcopy(provenance)
    corrupt["files"]["xgb_repeats.json"] = DIGEST
    monkeypatch.setattr(model.json, "loads", lambda _: corrupt)
    with pytest.raises(ModelError, match="checksum"):
        model.load_model()


def test_bounded_input_and_records(tmp_path):
    path = tmp_path / "big.fa"
    path.write_bytes(b">a\n" + b"A" * (2 * 1024 * 1024))
    with pytest.raises(InputError, match="MiB"):
        read_fasta(path)
    with pytest.raises(InputError, match="10000"):
        classify([Record(str(i), "A") for i in range(10001)], input_sha256=DIGEST)


def test_output_symlink(tmp_path):
    (tmp_path / "dir").mkdir()
    (tmp_path / "link").symlink_to(tmp_path / "dir", target_is_directory=True)
    with pytest.raises(OutputError):
        publish({}, tmp_path / "link/out")


def test_short_repeat_features():
    values = matrix(["A", "AC", "ACG", "ACGT"])
    assert values.shape == (5, 138)
    assert np.isfinite(values).all()
