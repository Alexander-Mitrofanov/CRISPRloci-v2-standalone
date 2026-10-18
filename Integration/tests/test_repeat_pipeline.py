from pathlib import Path
import json
import shutil

import pytest

from crispr_integration.repeat_validation import validate_repeat_report
from crispr_integration.repeats import run_repeat_pipeline
from crispr_integration.runner import RunnerConfig

FIXTURE = Path(__file__).resolve().parents[2] / "backend/tests/fixtures/repeats"


@pytest.mark.parametrize("kind", ["context", "isolated", "rna"])
def test_real_numeric_report(kind, tmp_path):
    report = json.loads((FIXTURE / f"{kind}.json").read_text())
    fasta = FIXTURE / {"context": "source.fasta", "isolated": "repeats.fasta", "rna": "rna.fasta"}[kind]
    reports = None
    if kind == "context":
        reports = tmp_path / "identify"
        reports.mkdir()
        shutil.copyfile(FIXTURE / "identify.json", reports / "report.json")
    validate_repeat_report(report, fasta, reports, "RNA" if kind == "rna" else "DNA", ["Bona-fide", "Possible"])
    report["instances"][0]["rna_sequence"] += "A"
    with pytest.raises(ValueError, match="RNA mapping"):
        validate_repeat_report(report, fasta, reports, "RNA" if kind == "rna" else "DNA", ["Bona-fide", "Possible"])


@pytest.mark.parametrize("mutation", [
    "missing_context", "duplicate_context", "orientation", "context_binding", "offset",
    "nonfinite_support", "support_type", "energy", "unpaired", "missing_score", "consensus", "instance_id",
    "count_type", "offset_type", "index_type", "terminal_type",
])
def test_repeat_context_rejects_tampered_native_evidence(tmp_path, mutation):
    report = json.loads((FIXTURE / "context.json").read_text())
    reports = tmp_path / "identify"
    reports.mkdir()
    shutil.copyfile(FIXTURE / "identify.json", reports / "report.json")
    row = report["instances"][0]
    if mutation == "missing_context":
        report["contexts"].pop()
    elif mutation == "duplicate_context":
        report["contexts"][1] = report["contexts"][0]
    elif mutation == "orientation":
        report["contexts"][0]["orientation"] = "unknown"
    elif mutation == "context_binding":
        row["context_id"] = report["contexts"][1]["id"]
    elif mutation == "offset":
        row["context_start"] += 1
    elif mutation == "nonfinite_support":
        row["comparison"]["support_change"] = float("nan")
    elif mutation == "support_type":
        row["comparison"]["support_change"] = True
    elif mutation == "energy":
        row["isolated"]["mfe_kcal_mol"] = float("inf")
    elif mutation == "unpaired":
        row["isolated"]["unpaired"][0] = 1.0
    elif mutation == "missing_score":
        row["comparison"]["context_pair_support"] = None
    elif mutation == "consensus":
        report["consensus_folds"].pop()
    elif mutation == "count_type":
        report["counts"]["arrays"] = True
    elif mutation == "offset_type":
        row["context_start"] = False
    elif mutation == "index_type":
        row["repeat_index"] = True
    elif mutation == "terminal_type":
        row["terminal"] = 1
    else:
        row["id"] = "forged"
    with pytest.raises((ValueError, KeyError)):
        validate_repeat_report(report, FIXTURE / "source.fasta", reports, "DNA", ["Bona-fide", "Possible"])


@pytest.mark.parametrize("context", [False, True])
def test_isolated_runner_commands(tmp_path, monkeypatch, context):
    from crispr_integration import runner

    calls = []
    job = tmp_path / "output"

    def detect(source, output, **kwargs):
        assert kwargs["tools"] == ()
        path = Path(output) / "crispridentify/source/report.json"
        path.parent.mkdir(parents=True)
        (Path(output) / "logs").mkdir()
        shutil.copyfile(FIXTURE / "identify.json", path)
        return {
            "schema": {"name": "crispr-tool-integration-run", "version": "1.0.0"},
            "status": "completed",
            "stages": {"crispridentify": {"status": "completed"}},
        }

    def execute(stage, command, logs, *, timeout_seconds):
        calls.append(command)
        assert stage == "crisprrepeat" and command[0] == "/isolated/crisprrepeat"
        assert timeout_seconds == 60
        native = Path(command[command.index("--output") + 1])
        native.mkdir()
        shutil.copyfile(FIXTURE / ("context.json" if context else "isolated.json"), native / "repeat-analysis.json")
        return {"status": "completed"}

    monkeypatch.setattr(runner, "run_pipeline", detect)
    monkeypatch.setattr(runner, "_run_process", execute)
    config = RunnerConfig(("identify",), ("spacerplacer",), ("evor",), crisprrepeat=("/isolated/crisprrepeat",))
    result = run_repeat_pipeline(
        FIXTURE / ("source.fasta" if context else "repeats.fasta"),
        job,
        runner_config=config,
        context=context,
        timeout_seconds=60,
    )
    assert result["status"] == "completed"
    assert ("--reports" in calls[0]) is context
    assert calls[0][calls[0].index("--orientation") + 1] == "both"
