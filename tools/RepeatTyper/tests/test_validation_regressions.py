"""Regressions at the input, immutable-model and publication boundaries."""
import hashlib
import json
import shutil
from importlib.resources import files

import numpy as np
import pytest

from repeattyper import Record, __version__, classify, inference, model, output
from repeattyper.cli import main
from repeattyper.errors import InputError, ModelError, OutputError

DIGEST = "0" * 64


@pytest.mark.parametrize("records", [None, "ACGT", {}, [None], [Record(None, "ACGT")],
                                   [Record("a", None)], [Record("a", ["A"])],
                                   [Record("a", "ACGT"), Record("a", "ACGT")]])
def test_malformed_api_records_are_input_errors(records):
    with pytest.raises(InputError):
        classify(records, input_sha256=DIGEST)


@pytest.mark.parametrize("digest", [None, 42, b"0" * 64, "G" * 64, "0" * 63])
def test_malformed_digests_are_input_errors(digest):
    with pytest.raises(InputError, match="SHA-256"):
        classify([], input_sha256=digest)


@pytest.mark.parametrize("threshold", [10**10000, -(10**10000), "0.75", None],
                         ids=["huge-positive", "huge-negative", "string", "null"])
def test_malformed_numeric_thresholds_are_input_errors(threshold):
    with pytest.raises(InputError, match="Threshold"):
        classify([], input_sha256=DIGEST, threshold=threshold)


@pytest.fixture
def bundled_data_copy(tmp_path, monkeypatch):
    data = tmp_path / "data"
    shutil.copytree(files("repeattyper").joinpath("data"), data)
    monkeypatch.setattr(model, "files", lambda _: tmp_path)
    return data


def test_provenance_cannot_change_its_own_expected_checksums(bundled_data_copy):
    data = bundled_data_copy
    provenance = json.loads((data / "provenance.json").read_bytes())
    changed_model = (data / "xgb_repeats.json").read_bytes() + b"\n"
    (data / "xgb_repeats.json").write_bytes(changed_model)
    provenance["files"]["xgb_repeats.json"] = hashlib.sha256(changed_model).hexdigest()
    provenance["upstream_commit"] = "0" * 40
    (data / "provenance.json").write_text(json.dumps(provenance))
    with pytest.raises(ModelError, match="provenance checksum"):
        classify([], input_sha256=DIGEST)


@pytest.mark.parametrize("filename", ["provenance.json", "xgb_repeats.json", "type_dict.tab"])
def test_missing_bundled_resources_have_model_exit_code(bundled_data_copy, tmp_path, filename):
    (bundled_data_copy / filename).unlink()
    fasta = tmp_path / "repeats.fa"
    fasta.write_text(">a\nACGT\n")
    destination = tmp_path / "results"
    assert main(["classify", "--input", str(fasta), "--output", str(destination)]) == 3
    assert not destination.exists()


def test_model_is_loaded_from_the_bytes_that_were_verified(bundled_data_copy, monkeypatch):
    resource = bundled_data_copy / "xgb_repeats.json"
    real_read = type(resource).read_bytes
    reads = []

    def read_once(path):
        if path == resource:
            reads.append(path)
            if len(reads) > 1:
                return b"changed after checksum validation"
        return real_read(path)

    monkeypatch.setattr(type(resource), "read_bytes", read_once)
    assert model.load_model()[1]["id"] == "crisprcastyper-1.9.0-repeat-model"
    assert reads == [resource]


@pytest.mark.parametrize("scores", [
    [[1.0] + [0.0] * 36] * 2,
    np.zeros((1, 37), dtype=np.float32),
    np.zeros((2, 37), dtype=np.float32),
    np.ones((2, 37), dtype=np.float32),
    np.full((2, 37), np.nan, dtype=np.float32),
    np.full((2, 37), np.inf, dtype=np.float32),
    np.full((2, 37), -1, dtype=np.float32),
    np.full((2, 37), "invalid", dtype=object),
    np.eye(2, 37, dtype=bool),
])
def test_invalid_model_scores_are_model_errors(monkeypatch, scores):
    booster = type("Booster", (), {"predict": lambda *args, **kwargs: scores})()
    monkeypatch.setattr(inference, "load_model", lambda: (booster, {"labels": ["I-A"] * 37}))
    with pytest.raises(ModelError, match="invalid score matrix"):
        classify([Record("a", "ACGT")], input_sha256=DIGEST)


@pytest.mark.parametrize("objective,classes", [("multi:softmax", "37"), ("multi:softprob", "36")])
def test_incompatible_objective_is_rejected(monkeypatch, objective, classes):
    configuration = {"learner": {"objective": {"name": objective},
                                  "learner_model_param": {"num_class": classes}}}
    booster = type("Booster", (), {"load_model": lambda *args: None,
        "num_features": lambda *args: 138, "attr": lambda *args: "41",
        "save_config": lambda *args: json.dumps(configuration)})()
    monkeypatch.setattr(model.xgb, "Booster", lambda *args: booster)
    with pytest.raises(ModelError, match="contract mismatch"):
        model.load_model()


@pytest.mark.parametrize("create_existing", ["directory", "file", "symlink"])
def test_concurrent_output_creation_is_never_overwritten(tmp_path, monkeypatch, create_existing):
    report = classify([], input_sha256=DIGEST)
    destination = tmp_path / "output"
    target = tmp_path / "unrelated"
    target.write_text("preserve target")
    real_rename = output._rename_new_directory
    inodes = []

    def race(source, destination):
        if create_existing == "directory":
            destination.mkdir()
        elif create_existing == "file":
            destination.write_text("preserve output")
        else:
            destination.symlink_to(target)
        inodes.append(destination.lstat().st_ino)
        return real_rename(source, destination)

    monkeypatch.setattr(output, "_rename_new_directory", race)
    with pytest.raises(OutputError):
        output.publish(report, destination)
    assert destination.lstat().st_ino == inodes[0]
    assert target.read_text() == "preserve target"
    if create_existing == "directory":
        assert list(destination.iterdir()) == []
    elif create_existing == "file":
        assert destination.read_text() == "preserve output"
    else:
        assert destination.is_symlink()
    assert not list(tmp_path.glob(".repeattyper-*"))


def test_patch_version_keeps_model_policy_and_schema():
    report = classify([Record("a", "ACGTN")], input_sha256=DIGEST)
    assert __version__ == "2.0.1"
    assert report["tool"] == {"name": "RepeatTyper", "version": __version__}
    assert report["schema"] == {"name": "repeattyper.classification", "version": "1.0.0"}
    assert report["policy"]["threshold"] == 0.75
    assert report["results"][0]["winner"] is None
    assert report["results"][0]["native_score"] is None
