"""Publication consistency, compatibility, and failure-boundary contracts."""

import csv
import hashlib
import json
import os
from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path

import pytest

from crisprtracrrna import pipeline, reporting
from crisprtracrrna.config import Settings
from crisprtracrrna.errors import ToolError, TracrError


def predictions():
    return [
        {
            "id": "model_1",
            "interval": {"start": 2, "end": 4},
            "strand": "-",
            "sequence": "ACG",
            "best_score": 20,
            "best_evalue": 1e-8,
        },
        {
            "id": "short_candidate_1",
            "interval": {"start": 5, "end": 7},
            "strand": "+",
            "sequence": "TTA",
            "ranking": {"score": 3.25},
            "candidate_origin": "anti_repeat",
            "admission_branch": "short_rescue",
            "evidence_class": "interaction_cm_terminator",
            "boundary_kind": "transcript_hypothesis",
            "array": {"id": "array,1", "category": "Bona-fide"},
            "anti_repeat": {"identity": 0.9, "coverage": 0.8},
            "interaction": {"energy_kcal_mol": -10.5},
            "model_hit": {"score": 22.5},
            "terminator": {"evalue": 0.001},
            "cas_context": {"id": "cassette_1", "distance_bases": 25},
            "short_rescue_admission": {"accepted": True},
        },
    ]


@pytest.fixture
def mocked_run(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs,quoted"
    inputs.mkdir()
    (inputs / "a.fasta").write_text(">shared;id%=x first header\nACGTTTACGT\n>empty\nAAAAAAAAAA\n")
    (inputs / "b.fasta").write_text(">shared;id%=x second header\nACGTTTACGT\n")

    def analyze(source, work, settings):
        work.mkdir()
        rows = predictions() if source.id == "f000001.s000001" else []
        if source.id == "f000002.s000001":
            rows = predictions()[:1]
        return {
            "predictions": rows,
            # Accepted short-rescue proposals alias their prediction dictionaries.
            "short_rescue_proposals": rows[1:],
            "arrays": [],
            "cas": {"status": "not_requested"},
        }

    monkeypatch.setattr(pipeline, "doctor", lambda settings: {"ready": True})
    monkeypatch.setattr(pipeline, "analyze_record", analyze)
    monkeypatch.setattr(pipeline, "command_provenance", lambda *args: {"command": []})
    model = tmp_path / "model.cm"
    model.write_text("mock covariance model\n")
    return inputs, tmp_path / "output", Settings(mode="model", models=(str(model),))


def test_retained_records_match_final_report_and_export_identities(mocked_run, monkeypatch):
    inputs, output, settings = mocked_run
    publications = []

    def publish(path, value):
        publications.append(path.relative_to(output).as_posix())
        reporting.atomic_json(path, value)

    monkeypatch.setattr(pipeline, "atomic_json", publish)
    returned = pipeline.run(inputs, output, settings)
    report = json.loads((output / "report.json").read_text())
    assert report["records"] == returned["records"]
    assert report["summary"] == {"sequence_count": 3, "prediction_count": 3}
    assert publications[-2:] == ["manifest.json", "report.json"]
    assert not (output / "failure.json").exists()
    identities = []
    for record in report["records"]:
        retained = json.loads((output / "work" / record["key"] / "record.json").read_text())
        assert retained == record
        identities.extend(prediction["export_id"] for prediction in record["predictions"])
    assert identities == [
        "f000001.s000001.p000001",
        "f000001.s000001.p000002",
        "f000002.s000001.p000001",
    ]
    first = report["records"][0]
    assert first["short_rescue_proposals"][0] == first["predictions"][1]
    assert report["records"][1]["predictions"] == []
    with (output / "summary.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["prediction_id"] for row in rows] == identities
    assert rows[0]["source_file"] == str(inputs / "a.fasta")
    assert rows[0]["sequence_id"] == rows[2]["sequence_id"] == "shared;id%=x"
    fasta = (output / "predictions.fasta").read_text()
    gff = (output / "predictions.gff3").read_text()
    for identity in identities:
        assert f">{identity} " in fasta
        assert f"ID={identity};" in gff
    assert "source=shared%3Bid%25%3Dx" in fasta
    assert "f000001.s000001:shared%3Bid%25%3Dx" in gff
    assert "f000002.s000001:shared%3Bid%25%3Dx" in gff
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["algorithm"] == "sha256"
    assert set(manifest["files"]) == {
        "configuration.json",
        "summary.csv",
        "predictions.fasta",
        "predictions.gff3",
        *(f"work/{record['key']}/record.json" for record in report["records"]),
    }
    for name, digest in manifest["files"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest


def test_export_bytes_preserve_model_defaults_and_complete_evidence(tmp_path):
    """These literal bytes were checked against the frozen pre-refactor exporter."""
    record = {
        "key": "r1",
        "source": {"file": 'some,"file.fasta', "id": "seq;1%=x", "length": 10},
        "predictions": predictions(),
    }
    original = deepcopy(record)
    pipeline.write_exports(tmp_path, [record])
    expected_csv = (
        "source_file,sequence_id,prediction_id,candidate_origin,admission_branch,evidence_class,"
        "boundary_kind,start,end,strand,sequence,score,evalue,array_id,array_category,"
        "anti_repeat_identity,anti_repeat_coverage,interaction_energy_kcal_mol,model_score,"
        "terminator_evalue,cas_cassette_id,cas_distance_bases\r\n"
        '"some,""file.fasta",seq;1%=x,r1.p000001,covariance_model,model_mode,cm_only,'
        "merged_model_region,2,4,-,ACG,20,1e-08,,,,,,,,,\r\n"
        '"some,""file.fasta",seq;1%=x,r1.p000002,anti_repeat,short_rescue,'
        "interaction_cm_terminator,transcript_hypothesis,5,7,+,TTA,3.25,,"
        '"array,1",Bona-fide,0.9,0.8,-10.5,22.5,0.001,cassette_1,25\r\n'
    )
    expected_fasta = (
        ">r1.p000001 source=seq%3B1%25%3Dx start=2 end=4 strand=-\nACG\n"
        ">r1.p000002 source=seq%3B1%25%3Dx start=5 end=7 strand=+\nTTA\n"
    )
    expected_gff = (
        "##gff-version 3\n##sequence-region r1:seq%3B1%25%3Dx 1 10\n"
        "r1:seq%3B1%25%3Dx\tCRISPRtracrRNA\tncRNA\t2\t4\t.\t-\t.\t"
        "ID=r1.p000001;inference=ab%20initio%20prediction;product=tracrRNA%20candidate\n"
        "r1:seq%3B1%25%3Dx\tCRISPRtracrRNA\tncRNA\t5\t7\t.\t+\t.\t"
        "ID=r1.p000002;inference=ab%20initio%20prediction;product=tracrRNA%20candidate\n"
    )
    assert (tmp_path / "summary.csv").read_bytes() == expected_csv.encode()
    assert (tmp_path / "predictions.fasta").read_bytes() == expected_fasta.encode()
    assert (tmp_path / "predictions.gff3").read_bytes() == expected_gff.encode()
    for ordinal, prediction in enumerate(record["predictions"], 1):
        assert prediction.pop("export_id") == f"r1.p{ordinal:06d}"
    assert record == original


@pytest.mark.parametrize(
    ("stage", "completed"),
    [
        ("work", 0),
        ("settings", 0),
        ("configuration.json", 0),
        ("record.json", 1),
        ("exports", 3),
        ("manifest.json", 3),
        ("report.json", 3),
    ],
)
def test_publication_failures_leave_failure_marker_with_persisted_count(
    mocked_run, monkeypatch, stage, completed
):
    inputs, output, settings = mocked_run
    primary = OSError(f"injected {stage} failure")

    def fail(*args, **kwargs):
        raise primary

    original_mkdir = Path.mkdir

    def mkdir(path, *args, **kwargs):
        if path == output / "work":
            raise primary
        return original_mkdir(path, *args, **kwargs)

    def publish(path, value):
        if path.name == stage and (stage != "record.json" or path.parent.name == "f000001.s000002"):
            raise primary
        reporting.atomic_json(path, value)

    if stage == "work":
        monkeypatch.setattr(Path, "mkdir", mkdir)
    elif stage == "settings":
        monkeypatch.setattr(Settings, "public", fail)
    elif stage == "exports":
        monkeypatch.setattr(pipeline, "write_exports", fail)
    else:
        monkeypatch.setattr(pipeline, "atomic_json", publish)
    with pytest.raises(OSError) as caught:
        pipeline.run(inputs, output, settings)
    assert caught.value is primary
    assert not (output / "report.json").exists()
    assert len(list(output.glob("work/*/record.json"))) == completed
    assert json.loads((output / "failure.json").read_text()) == {
        "status": "failed",
        "error": str(primary),
        "completed_sequences": completed,
    }


@pytest.mark.parametrize("error_type", [ToolError, KeyboardInterrupt, SystemExit])
def test_failure_marker_error_never_masks_primary_exception(
    mocked_run, monkeypatch, caplog, error_type
):
    inputs, output, settings = mocked_run
    primary = error_type("primary scientific failure or interruption")

    def analyze(*args):
        raise primary

    def publish(path, value):
        if path.name == "failure.json":
            raise OSError("secondary storage failure")
        reporting.atomic_json(path, value)

    monkeypatch.setattr(pipeline, "analyze_record", analyze)
    monkeypatch.setattr(pipeline, "atomic_json", publish)
    with pytest.raises(error_type) as caught:
        pipeline.run(inputs, output, settings)
    assert caught.value is primary
    assert not (output / "report.json").exists()
    assert "Could not publish failure marker" in caplog.text
    assert "secondary storage failure" in caplog.text


def test_preflight_and_existing_output_remain_untouched(mocked_run, monkeypatch):
    inputs, output, settings = mocked_run
    monkeypatch.setattr(pipeline, "doctor", lambda _: {"ready": False, "errors": ["not ready"]})
    with pytest.raises(TracrError, match="Preflight failed"):
        pipeline.run(inputs, output, settings)
    assert not output.exists()
    output.mkdir()
    retained = output / "report.json"
    retained.write_bytes(b"existing successful run\n")
    monkeypatch.setattr(pipeline, "doctor", lambda _: {"ready": True})
    with pytest.raises(TracrError, match="Refusing to overwrite"):
        pipeline.run(inputs, output, settings)
    assert list(output.iterdir()) == [retained]
    assert retained.read_bytes() == b"existing successful run\n"


@pytest.mark.parametrize("value", [object(), float("nan"), float("inf"), -float("inf")])
def test_atomic_serialization_failure_preserves_previous_file(tmp_path, value):
    destination = tmp_path / "report.json"
    destination.write_bytes(b"previous publication\n")
    with pytest.raises((TypeError, ValueError)):
        reporting.atomic_json(destination, {"value": value})
    assert destination.read_bytes() == b"previous publication\n"
    assert list(tmp_path.iterdir()) == [destination]


@pytest.mark.parametrize("stage", ["write", "replace"])
def test_atomic_io_failure_cleans_partial_file_and_preserves_destination(
    tmp_path, monkeypatch, stage
):
    destination = tmp_path / "report.json"
    destination.write_bytes(b"previous publication\n")
    primary = OSError(f"injected {stage} failure")
    original_open = Path.open

    @contextmanager
    def failing_open(path, *args, **kwargs):
        with original_open(path, *args, **kwargs) as handle:

            class FailingWriter:
                def write(self, data):
                    handle.write(data[:3])
                    raise primary

            yield FailingWriter()

    def replace(*args):
        raise primary

    with monkeypatch.context() as patch:
        if stage == "write":
            patch.setattr(Path, "open", failing_open)
        else:
            patch.setattr(Path, "replace", replace)
        with pytest.raises(OSError) as caught:
            reporting.atomic_json(destination, {"value": 1})
    assert caught.value is primary
    assert destination.read_bytes() == b"previous publication\n"
    assert list(tmp_path.iterdir()) == [destination]


def test_atomic_cleanup_error_is_logged_without_hiding_primary(tmp_path, monkeypatch, caplog):
    destination = tmp_path / "report.json"
    destination.write_bytes(b"previous publication\n")
    primary = OSError("primary replace failure")

    def replace(*args):
        raise primary

    def unlink(*args, **kwargs):
        raise OSError("secondary cleanup failure")

    with monkeypatch.context() as patch:
        patch.setattr(Path, "replace", replace)
        patch.setattr(Path, "unlink", unlink)
        with pytest.raises(OSError) as caught:
            reporting.atomic_json(destination, {"value": 1})
    assert caught.value is primary
    assert destination.read_bytes() == b"previous publication\n"
    assert "Could not remove temporary JSON file" in caplog.text
    assert "secondary cleanup failure" in caplog.text


def test_atomic_attempts_use_distinct_siblings_and_leave_existing_tmp_untouched(
    tmp_path, monkeypatch
):
    destination = tmp_path / "report.json"
    old_temporary = tmp_path / "report.json.tmp"
    old_temporary.write_bytes(b"belongs to another attempt\n")
    attempts = []
    original_replace = Path.replace

    def replace(temporary, target):
        attempts.append(temporary)
        assert temporary.parent == destination.parent
        assert target == destination
        return original_replace(temporary, target)

    monkeypatch.setattr(Path, "replace", replace)
    reporting.atomic_json(destination, {"value": "first"})
    reporting.atomic_json(destination, {"value": "café"})
    assert len(set(attempts)) == 2
    assert old_temporary.read_bytes() == b"belongs to another attempt\n"
    assert destination.read_bytes() == b'{\n  "value": "caf\\u00e9"\n}\n'
    assert set(tmp_path.iterdir()) == {destination, old_temporary}


@pytest.mark.skipif(os.name != "posix", reason="POSIX umask permissions")
def test_atomic_publication_respects_normal_file_creation_permissions(tmp_path):
    reference = tmp_path / "normal.json"
    reference.write_text("{}\n")
    destination = tmp_path / "atomic.json"
    reporting.atomic_json(destination, {})
    assert destination.stat().st_mode & 0o777 == reference.stat().st_mode & 0o777


def test_pipeline_publication_helpers_remain_reexports():
    assert pipeline.atomic_json is reporting.atomic_json
    assert pipeline.write_exports is reporting.write_exports
