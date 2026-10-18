"""Observed detection handoff, real installed CCTK and independent publication checks."""

import copy
from dataclasses import replace
import hashlib
import io
import json
import os
import shutil
import sys
import zipfile
import pytest
from crispr_integration.array_compare import run_array_compare_pipeline
from crispr_integration.runner import RunnerConfig
from crispr_web.array_exchange import array_input, native_input
from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.fasta import FastaValidationError
from crispr_web.worker import Worker
from .test_viral_search import configured_settings
from .test_sequence_exchange import detection_ports


def sha(value):
    return hashlib.sha256(
        value.encode() if isinstance(value, str) else value
    ).hexdigest()


@pytest.fixture
def native():
    executable = os.environ.get("CCTK_EXECUTABLE") or shutil.which("cctk")
    if not executable:
        pytest.skip("Install pinned CCTK or set CCTK_EXECUTABLE")
    return (executable,)


def rich_packet(tmp_path):
    doc = detection_ports(tmp_path)["arrays.json"]
    first = doc["arrays"][0]
    first["strand"] = None
    first["repeats"] = None
    # A separate source with a colliding array/occurrence ID must remain separate.
    source = copy.deepcopy(doc["sources"][0])
    source["id"] = "other"
    doc["sources"].append(source)
    second = copy.deepcopy(first)
    second["source_id"] = "other"
    second["strand"] = "-"
    doc["arrays"].append(second)
    # Preserve a duplicate occurrence at its own coordinate and a deletion slot.
    end = first["spacers"][-1]["end"]
    duplicate = copy.deepcopy(first["spacers"][0])
    duplicate.update(
        id="gap000004", ordinal=3, start=end + 1, end=end + len(duplicate["sequence"])
    )
    first["interval"]["end"] = duplicate["end"]
    first["spacers"].append(duplicate)
    first["spacers"].append(
        {
            "ordinal": 4,
            "sequence": None,
            "status": "deletion",
            "boundary_position": duplicate["end"],
        }
    )
    return doc


@pytest.fixture
def pipeline(tmp_path, native):
    doc = rich_packet(tmp_path)
    text, mapping = array_input(json.dumps(doc), max_bytes=1_000_000)
    folder = tmp_path / "job/input/records"
    folder.mkdir(parents=True)
    native_doc, records = native_input(doc)
    for i, (key, seq) in enumerate(records.items()):
        (folder / f"{i}.fasta").write_text(f">{key}\n{seq}\n")
    mapfile = folder.parent.parent / "array-map.json"
    mapfile.write_text(json.dumps(mapping))
    output = tmp_path / "out"
    run_array_compare_pipeline(
        folder,
        output,
        runner_config=RunnerConfig(("unused",), ("unused",), ("unused",), cctk=native),
        map_path=mapfile,
        timeout_seconds=60,
    )
    claim = ClaimedJob(
        "a" * 32,
        "array_compare",
        1,
        0,
        "bona_fide_possible",
        1,
        True,
        True,
        (),
        array_map_sha256=sha(mapfile.read_bytes()),
    )
    return (
        claim,
        output,
        {key: (len(seq), sha(seq)) for key, seq in records.items()},
        folder,
    )


def valid(data):
    claim, out, bindings, folder = data
    return _validated_pipeline_status(claim, out, bindings, input_root=folder)


def test_real_duplicate_scoped_unknown_and_deletion_comparison(pipeline):
    assert valid(pipeline) == "completed"
    report = json.loads((pipeline[1] / "array-compare/comparison.json").read_text())
    assert report["edges"][0]["nshared"] == 2 and report["edges"][0]["jaccard"] == 1
    assert (
        report["arrays"][0]["strand"] is None and report["arrays"][1]["strand"] == "-"
    )
    assert len(report["arrays"][0]["spacers"]) == 4
    assert report["arrays"][0]["spacers"][3]["id"] != "gap000004"


@pytest.mark.parametrize(
    "mutation",
    [
        "sequence",
        "nshared",
        "jaccard",
        "components",
        "strand",
        "order",
        "symbol",
        "input_hash",
        "map",
        "fasta",
    ],
)
def test_independent_tamper_rejected(pipeline, mutation):
    claim, out, bindings, folder = pipeline
    path = out / "array-compare/comparison.json"
    doc = json.loads(path.read_text())
    if mutation == "sequence":
        doc["arrays"][0]["spacers"][0]["sequence"] = "ACGT"
    elif mutation == "nshared":
        doc["edges"][0]["nshared"] = 3
    elif mutation == "jaccard":
        doc["edges"][0]["jaccard"] = 0.1
    elif mutation == "components":
        doc["components"] = [["array000000"], ["array000001"]]
    elif mutation == "strand":
        doc["arrays"][0]["strand"] = "+"
    elif mutation == "order":
        doc["arrays"].reverse()
    elif mutation == "symbol":
        doc["arrays"][0]["spacers"][0]["symbol"] = "s_" + "0" * 64
    elif mutation == "input_hash":
        doc["input_sha256"] = "0" * 64
    elif mutation == "map":
        (folder.parent.parent / "array-map.json").write_text("{}")
    elif mutation == "fasta":
        next(folder.iterdir()).write_text(">changed\nACGT\n")
    path.write_text(json.dumps(doc))
    assert valid(pipeline) is None


@pytest.mark.parametrize(
    "change",
    [
        lambda d: d["arrays"][0]["spacers"][0].update(sequence=None),
        lambda d: d["arrays"][0]["spacers"][0].update(
            sequence="N" * len(d["arrays"][0]["spacers"][0]["sequence"])
        ),
        lambda d: d["arrays"][0]["spacers"][0].update(start=1),
        lambda d: d["arrays"][0]["spacers"][1].update(ordinal=1),
        lambda d: d["arrays"][0].update(strand="forward"),
        lambda d: d["sources"][0].update(length=1),
    ],
)
def test_malformed_observations_rejected(tmp_path, change):
    doc = detection_ports(tmp_path)["arrays.json"]
    change(doc)
    with pytest.raises(FastaValidationError):
        array_input(json.dumps(doc), max_bytes=1_000_000)


def comparison_settings(factory, tmp_path):
    settings = configured_settings(
        factory,
        tmp_path,
        integration_command=(sys.executable, "-m", "crispr_integration"),
    )
    path = settings.scientific_release_manifest
    doc = json.loads(path.read_text())
    doc["schema_version"] = "9.0.0"
    doc["scientific_artifacts"]["cctk"] = {
        **doc["scientific_artifacts"]["crisprspacer"],
        "name": "cctk",
        "display_name": "CCTK",
        "version": "2.0.1",
    }
    path.write_text(json.dumps(doc))
    return replace(settings, scientific_release_manifest_sha256=sha(path.read_bytes()))


def test_http_detection_handoff_worker_and_protected_exports(
    settings_factory, tmp_path, native
):
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app

    packet = detection_ports(tmp_path)["arrays.json"]
    settings = comparison_settings(settings_factory, tmp_path)
    settings.runner_config.write_text(
        json.dumps(
            {
                "crispridentify": ["unused"],
                "spacerplacer": ["unused"],
                "crispr_evor": ["unused"],
                "cctk": list(native),
            }
        )
    )
    with TestClient(create_app(settings)) as client:
        assert "array_compare" in client.get("/api/v1/config").json()["modes"]
        response = client.post(
            "/api/v1/jobs",
            json={"mode": "array_compare", "sequence": json.dumps(packet)},
        )
        assert response.status_code == 202, response.text
        data = response.json()
        jid = data["job"]["job_id"]
        headers = {"Authorization": "Bearer " + data["access_token"]}
        assert Worker(settings, worker_id="cctk-http").run_once()
        result = client.get("/api/v1/jobs/" + jid, headers=headers).json()
        assert result["status"] == "completed", result
        assert result["summary"]["array_compare"]["counts"]["isolates"] == 1
        exports = {a["name"]: a for a in result["artifacts"]}
        for name in [
            "array-comparison.json",
            "array-sharing.tsv",
            "arrays.json",
            "result.json",
            "results.zip",
        ]:
            assert client.get(exports[name]["download_url"]).status_code in {401, 404}
            response = client.get(exports[name]["download_url"], headers=headers)
            assert response.status_code == 200
            if name == "arrays.json":
                assert response.json() == packet
            if name == "result.json":
                assert response.json()["complete"] is True
            if name == "results.zip":
                with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
                    assert "results/array-comparison.json" in archive.namelist()
                    assert not any(
                        "array-compare/" in n or "array-compare-input/" in n
                        for n in archive.namelist()
                    )


def test_unverified_repeat_annotations_still_require_consistent_bounds(tmp_path):
    doc = detection_ports(tmp_path)["arrays.json"]
    unit = doc["arrays"][0]["repeats"][0]
    unit.update(sequence=None, status="not_source_verified")
    assert array_input(json.dumps(doc), max_bytes=1_000_000)
    unit["end"] = doc["sources"][0]["length"] + 1
    with pytest.raises(FastaValidationError):
        array_input(json.dumps(doc), max_bytes=1_000_000)


def test_integration_rejects_extra_fields_and_rewritten_interpretation(pipeline):
    from crispr_integration.array_compare import validate_report

    out = pipeline[1]
    native = json.loads((out / "array-compare-input/arrays.json").read_text())
    report = json.loads((out / "array-compare/comparison.json").read_text())
    for changes in [{"extra": "field"}, {"interpretation": "inferred phylogeny"}]:
        with pytest.raises(ValueError):
            validate_report(
                {**report, **changes}, native, report["input_sha256"], report["options"]
            )


def test_public_comparison_preserves_semantics_null_and_zero_annotations(tmp_path):
    from crispr_web.summaries.array_compare import public_report
    from crispr_web.publication.compact.array_compare import array_compare_data
    mapping = rich_packet(tmp_path)
    mapping["arrays"][0]["model_score"] = 0.0
    mapping["arrays"][1]["model_score"] = None
    rows = [{"source_id": row["source_id"], "array_id": row["id"], "spacers": row["spacers"], "observed_spacer_count": sum(s.get("status") != "deletion" for s in row["spacers"])} for row in mapping["arrays"]]
    for folder in ("array-compare", "array-compare-input"):
        (tmp_path / folder).mkdir()
    (tmp_path / "array-compare-input/map.json").write_text(json.dumps({"document": mapping, "submitted_sha256": "a" * 64}))
    native = {"tool": {"name": "CCTK", "version": "2.0.1"}, "arrays": rows, "options": {}, "edges": [], "components": [], "isolates": [], "interpretation": "Observed sharing"}
    (tmp_path / "array-compare/comparison.json").write_text(json.dumps(native))
    doc = public_report(tmp_path)
    assert doc["coordinate_system"] == mapping["coordinate_system"]
    assert doc["score_semantics"] == mapping["score_semantics"]
    assert [row["model_score"] for row in doc["arrays"]][:2] == [0.0, None]
    assert "claimed producers were not independently verified" in doc["origin_semantics"]
    class Reader:
        def json(self, name):
            assert name == "array-comparison.json"
            return doc
    compact = array_compare_data(Reader())
    assert compact["score_semantics"] == mapping["score_semantics"]
    assert compact["arrays"][0]["model_score"] == 0.0
    assert compact["arrays"][1]["model_score"] is None
