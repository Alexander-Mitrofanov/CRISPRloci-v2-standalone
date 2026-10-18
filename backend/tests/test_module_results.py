"""Every module publishes compact science, independently of browser preview caps."""

import json
from pathlib import Path
import zipfile

import pytest

from crispr_web.exchange import INPUT_KINDS
from crispr_web.models import JobSubmission
from crispr_web.publication.compact.common import Reader
from crispr_web.publication.module_result import (
    MODULES,
    module_result,
    write_module_result,
)
from crispr_web.service import JobService
from crispr_web.worker import Worker

from .test_annotations import annotation_output as annotation_output
from .test_leader import leader_output as leader_output
from .test_sequence_exchange import detection_ports

FIXTURES = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).resolve().parents[2]


def write(root, name, doc):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))


def load(path):
    return json.loads(path.read_text())


def exported(root, mode, status="completed"):
    assert write_module_result(root, mode, status, max_bytes=20_000_000) == []
    doc = load(root / "result.json")
    assert doc["module"] == mode and doc["status"] == status
    assert doc["complete"] is True
    assert set(doc) == {"schema", "module", "status", "complete", "data"}
    assert doc["schema"] == {"name": "crisprloci.module-result", "version": "1.0.0"}
    return doc["data"]


def test_every_submission_mode_has_a_compact_exporter():
    assert MODULES == set(INPUT_KINDS)
    schema = load(ROOT / "docs/schemas/module-result-v1.json")
    assert set(schema["properties"]["module"]["enum"]) == MODULES
    examples = ROOT / "docs/examples/module-results"
    assert {path.stem for path in examples.glob("*.json")} == MODULES
    for mode in MODULES:
        doc = load(examples / f"{mode}.json")
        assert doc["module"] == mode and doc["complete"] is True


def test_detection_keeps_source_evidence_and_uncertainty(tmp_path):
    doc = detection_ports(tmp_path)["arrays.json"]
    write(tmp_path, "arrays.json", doc)
    compact = exported(tmp_path, "detection")
    assert compact == {key: value for key, value in doc.items() if key != "schema"}
    assert compact["arrays"][0]["boundary"]["exact"] is False
    assert (
        len((tmp_path / "result.json").read_bytes())
        < (FIXTURES / "repeats/identify.json").stat().st_size / 5
    )


@pytest.mark.parametrize("kind", ["isolated", "rna", "context"])
def test_repeat_folds_are_deduplicated_without_losing_occurrences(tmp_path, kind):
    doc = load(FIXTURES / f"repeats/{kind}.json")
    write(tmp_path, "repeat/repeat-analysis.json", doc)
    mode = "repeat_context" if kind == "context" else "repeats"
    data = exported(tmp_path, mode)
    assert len(data["instances"]) == len(doc["instances"])
    folds = {fold["id"]: fold for fold in data["folds"]}
    for row, original in zip(data["instances"], doc["instances"]):
        fold = folds[row["fold_id"]]
        assert fold["rna_sequence"] == original["rna_sequence"]
        assert fold["status"] == original["isolated"]["status"]
        assert (
            row["comparison"]["support_change"]
            == original["comparison"]["support_change"]
        )
        for key in ("structure", "mfe_kcal_mol", "reason"):
            if key in original["isolated"]:
                assert fold[key] == original["isolated"][key]
    assert data["classification"] == doc["classification"]
    assert '"pairs"' not in json.dumps(data) and '"unpaired"' not in json.dumps(data)
    if kind == "context":
        assert len(folds) < len(data["instances"])
        assert len(data["contexts"]) == len(doc["contexts"])
        assert (tmp_path / "result.json").stat().st_size < (
            FIXTURES / "repeats/context.json"
        ).stat().st_size / 10
    if kind == "isolated":
        assert any(fold.get("mfe_kcal_mol") == 0 for fold in folds.values())
        assert any(fold.get("reason") == "ambiguous_bases" for fold in folds.values())
        assert data["instances"][0]["start"] is None


def test_tracr_uses_full_export_and_retains_missing_values(annotation_output):
    from crispr_web.publication.annotations import write_annotation_exports
    from crispr_web.summaries.annotations import annotation_summary

    root = annotation_output[0]
    cas, tracr, _ = annotation_summary(
        root, max_rows=0, accepted_categories=[], warnings=[]
    )
    assert tracr["truncated"] is True
    write_annotation_exports(root, {"casandra": cas, "tracrrna": tracr})
    data = exported(root, "tracrrna")
    assert len(data["candidates"]) == 1
    assert data["candidates"][0]["score"] == 0
    assert data["candidates"][0]["evalue"] is None
    assert data["candidates"][0]["strand"] == "-"
    assert "casandra" not in data and "detection" not in data


def test_leader_and_loci_preserve_both_sides_and_deduplicate_sources(leader_output):
    from crispr_web.publication.annotations import write_annotation_exports
    from crispr_web.summaries.leader import leader_summary

    root = leader_output[0]
    preview = leader_summary(root, max_rows=1)
    assert preview["truncated"] is True
    write_annotation_exports(root, {"crisprleader": preview})
    data = exported(root, "leader")
    assert data["prediction_performed"] is False
    assert len(data["arrays"]) == 1 and len(data["contexts"]) == 2
    assert {row["side"] for row in data["contexts"]} == {"left", "right"}
    assert all("array_segments" not in row for row in data["contexts"])
    source = data["sources"][0]
    # Compose already-normalized evidence; scientific stage validation has its
    # own independent tests. Empty Cas/tracr evidence must stay empty.
    write(root, "arrays.json", load(ROOT / "docs/examples/module-exchange/arrays.json"))
    write(
        root,
        "annotations.json",
        {
            "casandra": {
                "status": "completed",
                "truncated": False,
                "gene_count": 5,
                "genes": [],
                "cassettes": [],
            },
            "tracrrna": {
                "status": "completed",
                "truncated": False,
                "mode": "complete",
                "model_type": "II",
                "sources": [source],
                "candidates": [],
            },
        },
    )
    loci = exported(root, "loci")
    assert loci["sources"] == [source]
    assert loci["cas"]["genes"] == [] and loci["cas"]["evaluated_gene_count"] == 5
    assert loci["tracrrna"]["candidates"] == []
    assert all(
        "sources" not in loci[key] for key in ("detection", "cas", "tracrrna", "flanks")
    )
    doc = load(root / "annotations.json")
    doc["tracrrna"]["sources"][0]["sha256"] = "0" * 64
    write(root, "annotations.json", doc)
    with pytest.raises(ValueError, match="inconsistent locus"):
        module_result(root, "loci", "completed")


@pytest.mark.parametrize("empty", [False, True])
def test_cas_has_no_unrequested_placeholders_and_preserves_zero(tmp_path, empty):
    gene = {
        "id": "gene_1",
        "source_id": "contig",
        "start": 2,
        "end": 4,
        "strand": "-",
        "type": "II",
        "subtype": None,
        "profile": None,
        "score_margin": 0.0,
    }
    write(
        tmp_path,
        "annotations.json",
        {
            "casandra": {
                "status": "completed",
                "truncated": False,
                "gene_count": 10,
                "genes": [] if empty else [gene],
                "cassettes": [],
            },
            "tracrrna": {"status": "not_requested"},
            "detection": None,
        },
    )
    data = exported(tmp_path, "cas")
    assert data["genes"] == ([] if empty else [gene])
    assert data["evaluated_gene_count"] == 10
    assert not {"tracrrna", "detection", "model", "truncated"} & data.keys()


@pytest.mark.parametrize("empty", [False, True])
def test_protospacer_preserves_exactness_and_no_hit(tmp_path, empty):
    doc = load(FIXTURES / "protospacer/analysis.json")
    if empty:
        doc.update(matches=[], outcome="no_matches")
    else:
        doc["matches"].append({**doc["matches"][0], "relative_strand": "-"})
    write(tmp_path, "protospacer/analysis.json", doc)
    data = exported(tmp_path, "protospacer")
    assert data["outcome"] == doc["outcome"] and data["search_complete"] is True
    assert (
        data["policy"]["pam"] == "not_evaluated"
        and data["policy"]["substitutions"] == 0
    )
    assert len(data["spacers"]) == (0 if empty else 1)
    assert len(data["matches"]) == len(doc["matches"])
    assert all("spacer_sequence" not in hit for hit in data["matches"])


def test_viral_preserves_heuristic_limits_and_occurrence_identity(tmp_path):
    doc = load(ROOT / "frontend/tests/fixtures/viral.json")["summary"]["viral_search"]
    doc["queries"][0]["subject_limit_saturated"] = True
    write(tmp_path, "viral/result.json", doc)
    data = exported(tmp_path, "viral_search")
    assert data["retrieval_completeness"] == "heuristic"
    assert data["queries"][0]["subject_limit_saturated"] is True
    assert data["policy"]["host_inference"] == "not_performed"
    queries = {row["id"]: row for row in data["queries"]}
    subjects = {row["accession"] for row in data["subjects"]}
    for hit, raw in zip(data["matches"], doc["matches"]):
        assert queries[hit["query_id"]]["sequence"] == raw["query_sequence"]
        assert hit["subject_accession"] in subjects
        assert hit["relative_strand"] == raw["relative_strand"]
        assert "query_alignment" not in hit and "raw_hsp_row" not in hit
    assert all("taxid" in subject for subject in data["subjects"])


def test_repeat_map_keeps_ties_conflicts_and_incomplete_annotations(tmp_path):
    doc = load(ROOT / "frontend/tests/fixtures/repeat-map.json")["summary"][
        "repeat_map"
    ]
    write(tmp_path, "repeat-mapping.json", doc)
    data = exported(tmp_path, "repeat_map")
    exact = data["results"][0]
    assert exact["best_distance"] == 0 and exact["best_hit_count"] == 25
    assert len(exact["hits"]) == 20 and exact["hits_truncated"] is True
    assert exact["annotations"]["family"]["status"] == "conflicting"
    assert exact["annotations"]["motif"]["missing_reference_count"] == 1
    references = {row["reference_id"] for row in data["references"]}
    for row in data["results"]:
        for hit in row["hits"]:
            assert hit["reference_id"] in references
            assert "family" not in hit and "alignment_tie_policy" not in hit
    assert len(data["results"]) == len(doc["results"])


@pytest.mark.parametrize(
    "mode,case,prediction",
    [
        ("orientation", "forwardprediction", "Forward"),
        ("orientation", "reverseprediction", "Reverse"),
        ("orientation", "ndprediction", "ND"),
        ("orientation", "noeligible", None),
        ("reconstruction", "source", None),
        ("reconstruction", "noeligible", None),
        ("detection", "source", None),
    ],
)
def test_worker_registers_compact_result_and_bundles_it(
    settings_factory, mode, case, prediction
):
    settings = settings_factory(max_summary_rows=0)
    service = JobService(settings)
    _, queued = service.submit(
        JobSubmission(sequence=f">{case}\nACGT\n", mode=mode), "198.51.100.1"
    )
    assert Worker(settings, store=service.store, worker_id="compact").run_once()
    job = service.get(queued["job_id"])
    assert job["status"] in {"completed", "completed_no_eligible_groups"}
    artifact = next(row for row in job["artifacts"] if row["name"] == "result.json")
    path, _ = service.artifact_file(job["job_id"], artifact["artifact_id"])
    doc = load(path)
    assert doc["complete"] is True, doc
    assert doc["module"] == mode
    data = doc["data"]
    if case == "noeligible":
        assert data["groups"] == [] and data["excluded_arrays"]
    elif mode in {"orientation", "reconstruction"}:
        assert len(data["groups"]) == 1
        group = data["groups"][0]
        assert group["members"] and "ccdb_file" not in json.dumps(group)
        if prediction:
            assert group["orientation"]["prediction"] == prediction
            assert group["selected_hypothesis"] == (
                "reverse" if prediction == "Reverse" else "input"
            )
            assert group["selection_basis"] == (
                "unresolved_default" if prediction == "ND" else "supported"
            )
    bundle, _ = service.result_bundle(job["job_id"])
    with zipfile.ZipFile(bundle) as archive:
        assert json.loads(archive.read("results/result.json")) == doc


@pytest.mark.parametrize(
    "failure",
    [
        "missing",
        "symlink",
        "duplicate_key",
        "nonfinite",
        "private",
        "size",
        "truncated",
    ],
)
def test_unavailable_data_is_never_published_as_empty_science(tmp_path, failure):
    doc = load(FIXTURES / "protospacer/analysis.json")
    path = tmp_path / "protospacer/analysis.json"
    write(tmp_path, "protospacer/analysis.json", doc)
    if failure == "missing":
        path.unlink()
    elif failure == "symlink":
        path.rename(tmp_path / "other.json")
        path.symlink_to(tmp_path / "other.json")
    elif failure == "duplicate_key":
        path.write_text('{"matches":[],"matches":[]}')
    elif failure == "nonfinite":
        doc["matches"][0]["target_start"] = float("nan")
        write(tmp_path, "protospacer/analysis.json", doc)
    elif failure == "private":
        doc["reference"]["name"] = str(tmp_path / "private.db")
        write(tmp_path, "protospacer/analysis.json", doc)
    elif failure == "truncated":
        write(
            tmp_path,
            "annotations.json",
            {"casandra": {"status": "completed", "truncated": True}},
        )
    assert write_module_result(
        tmp_path,
        "cas" if failure == "truncated" else "protospacer",
        "completed",
        max_bytes=400 if failure == "size" else 20_000_000,
    )
    fallback = load(tmp_path / "result.json")
    assert (
        fallback["complete"] is False
        and fallback["reason"] == "compact_export_unavailable"
    )
    assert "data" not in fallback
    assert not list(tmp_path.glob(".*.tmp"))


def test_reader_accounts_bytes_once_and_bounds_csv_as_well(tmp_path):
    write(tmp_path, "one.json", {"value": "x" * 20})
    read = Reader(tmp_path)
    read.remaining = (tmp_path / "one.json").stat().st_size
    assert read.json("one.json") == read.json("one.json")
    assert read.remaining == 0
    with pytest.raises(ValueError):
        read.bytes("one.json")


def test_cancel_propagates_and_failed_jobs_do_not_publish(tmp_path):
    class Cancelled(RuntimeError):
        pass

    def cancel():
        raise Cancelled()

    write(
        tmp_path,
        "protospacer/analysis.json",
        load(FIXTURES / "protospacer/analysis.json"),
    )
    with pytest.raises(Cancelled):
        write_module_result(
            tmp_path,
            "protospacer",
            "completed",
            max_bytes=20_000_000,
            checkpoint=cancel,
        )
    for status in ("failed", "cancelled"):
        assert (
            write_module_result(tmp_path, "protospacer", status, max_bytes=20_000_000)
            == []
        )
    assert not (tmp_path / "result.json").exists()
