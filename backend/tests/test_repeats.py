"""Frozen ViennaRNA output crosses the independent worker/publication boundary."""

from dataclasses import replace
from pathlib import Path
import json
import shutil

import pytest

from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.fasta import FastaValidationError
from crispr_web.models import JobSubmission
from crispr_web.publication.repeats import write_repeat_exports
from crispr_web.service import JobService
from crispr_web.summaries.repeats import repeat_summary
from crispr_web.worker import Worker

FIXTURE = Path(__file__).parent / "fixtures/repeats"


def setup_output(tmp_path, kind="context"):
    root, inputs = tmp_path / "output", tmp_path / "input"
    inputs.mkdir()
    (root / "repeat").mkdir(parents=True)
    fasta = (
        FIXTURE
        / {"context": "source.fasta", "isolated": "repeats.fasta", "rna": "rna.fasta"}[
            kind
        ]
    )
    # Backend inputs contain one normalized record per file.
    for i, record in enumerate(fasta.read_text().split(">")[1:]):
        (inputs / f"{i:06d}.fasta").write_text(">" + record)
    document = json.loads((FIXTURE / f"{kind}.json").read_text())
    if kind == "context":
        native = root / "crispridentify/source/report.json"
        native.parent.mkdir(parents=True)
        shutil.copyfile(FIXTURE / "identify.json", native)
    mode = "repeat_context" if kind == "context" else "repeats"
    pipeline = {
        "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
        "status": "completed",
        "workflow": mode,
        "stages": {
            "crisprrepeat": {"status": "completed", **document["counts"]},
            "crispridentify": {"status": "completed"},
        },
    }

    def flush():
        (root / "repeat/repeat-analysis.json").write_text(json.dumps(document))
        (root / "integration_result.json").write_text(json.dumps(pipeline))

    flush()
    claimed = ClaimedJob(
        "a" * 32,
        mode,
        1,
        0,
        "bona_fide_possible",
        1,
        True,
        True,
        (),
        molecule="RNA" if kind == "rna" else "DNA",
    )
    bindings = {r["id"]: (r["length"], r["sha256"]) for r in document["sources"]}
    return root, inputs, document, pipeline, claimed, bindings, flush


def validate(data):
    root, inputs, _, _, claimed, bindings, _ = data
    return _validated_pipeline_status(claimed, root, bindings, input_root=inputs)


@pytest.mark.parametrize("kind", ["context", "isolated", "rna"])
def test_real_numeric_report_validates_and_exports(tmp_path, kind):
    data = setup_output(tmp_path, kind)
    assert validate(data) == "completed"
    root, _, doc, _, _, _, _ = data
    summary = repeat_summary(root, max_rows=1)
    assert len(summary["instances"]) == 1
    assert summary["truncated"] == (len(doc["instances"]) > 1)
    write_repeat_exports(root, {"repeats": summary})
    assert json.loads((root / "repeat-analysis.json").read_text()) == doc
    assert (root / "repeat-analysis.fasta").read_text().count(">") == len(
        doc["instances"]
    )
    if kind == "context":
        row = summary["instances"][0]
        assert row["context"]["effective_window"] == 150
        assert len(row["context"]["unpaired_probabilities"]) == len(row["rna_sequence"])


@pytest.mark.parametrize(
    "mutation",
    [
        "sequence",
        "offset",
        "probability",
        "support",
        "orientation",
        "omission",
        "counts",
        "mfe",
        "accessibility",
        "window",
        "source",
        "classify",
        "stage",
        "context_duplicate",
        "context_span",
        "symlink",
        "count_type",
        "offset_type",
        "index_type",
        "terminal_type",
    ],
)
def test_rejects_corrupt_context(tmp_path, mutation):
    data = setup_output(tmp_path)
    root, _, doc, pipeline, _, _, flush = data
    row = doc["instances"][0]
    if mutation == "sequence":
        row["rna_sequence"] = "A" + row["rna_sequence"][1:]
    elif mutation == "offset":
        row["context_start"] += 1
    elif mutation == "probability":
        row["isolated"]["pairs"][0][2] = 1.01
    elif mutation == "support":
        row["comparison"]["context_pair_support"] = 0.987654321
    elif mutation == "orientation":
        row["orientation"] = "reverse"
    elif mutation == "omission":
        doc["instances"].pop()
    elif mutation == "counts":
        pipeline["stages"]["crisprrepeat"]["instances"] += 1
    elif mutation == "mfe":
        row["isolated"]["structure"] = "." * len(row["rna_sequence"])
    elif mutation == "accessibility":
        row["isolated"]["unpaired"][0] = 0.5
    elif mutation == "window":
        doc["contexts"][0]["window"] = 100
    elif mutation == "source":
        doc["sources"][0]["sha256"] = "a" * 64
    elif mutation == "classify":
        doc["classification"]["strand"] = {"status": "completed", "strand": "+"}
    elif mutation == "stage":
        pipeline["stages"]["crispridentify"]["status"] = "failed"
    elif mutation == "context_duplicate":
        doc["contexts"][1] = doc["contexts"][0]
    elif mutation == "context_span":
        doc["contexts"][0]["pairs"].append([0, 90, 0.01])
    elif mutation == "count_type":
        doc["counts"]["arrays"] = True
    elif mutation == "offset_type":
        row["context_start"] = False
    elif mutation == "index_type":
        row["repeat_index"] = True
    elif mutation == "terminal_type":
        row["terminal"] = 1
    flush()
    if mutation == "symlink":
        path = root / "repeat/repeat-analysis.json"
        target = tmp_path / "outside.json"
        path.rename(target)
        path.symlink_to(target)
    assert validate(data) is None


def test_rna_policy_and_claim_roundtrip(settings_factory):
    service = JobService(settings_factory())
    service.store.initialize()
    _, job = service.submit(
        JobSubmission(mode="repeats", molecule="RNA", sequence=">r\nGGGUUUCCC\n"),
        "127.0.0.1",
    )
    claimed = service.store.claim_next("repeat-test")
    assert claimed.molecule == "RNA"
    command = Worker(service.settings).command_for(claimed, service.settings.data_root / "output")
    assert command[command.index("--molecule") + 1] == "RNA"
    assert job["options"]["molecule"] == "RNA"
    with pytest.raises(FastaValidationError, match="RNA input"):
        service.submit(
            JobSubmission(mode="detection", molecule="RNA", sequence=">r\nGGGUUUCCC\n"),
            "127.0.0.2",
        )
    with pytest.raises(FastaValidationError):
        service.submit(
            JobSubmission(mode="repeats", sequence=">r\n" + "A" * 201), "127.0.0.3"
        )
    assert replace(claimed, molecule="DNA").molecule == "DNA"
