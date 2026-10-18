"""Private option transport and argv use only existing scientific flags."""
from dataclasses import replace
import json

import pytest

from crispr_integration import cli
from crispr_integration.runner import RunnerConfig
from crispr_integration.scientific_options import (
    TOOL_OPTIONS_ENV, load_tool_options_environment, tool_arguments, validate_tool_options,
)


def test_exact_boolean_syntax_and_numeric_arguments():
    config = RunnerConfig(("identify",), ("spacerplacer",), ("evor",), tool_options={
        "crispridentify": {"--fast_run": False, "--min_repeats": 4},
        "crispr_tracrrna": {"--model-rescue": False},
        "crisprrepeat": {"--no-lonely-pairs": True},
    })
    assert tool_arguments(config, "crispridentify") == ["--fast_run", "False", "--min_repeats", "4"]
    assert tool_arguments(config, "crispr_tracrrna") == ["--no-model-rescue"]
    assert tool_arguments(config, "crisprrepeat") == ["--no-lonely-pairs"]
    disabled = replace(config, tool_options={"crisprrepeat": {"--no-lonely-pairs": False}})
    assert tool_arguments(disabled, "crisprrepeat") == []


def test_cli_transports_options_without_new_cli_flags(monkeypatch, tmp_path):
    options = {"crispridentify": {"--min_repeats": 4}}
    monkeypatch.setenv(TOOL_OPTIONS_ENV, json.dumps(options))
    received = {}

    def pipeline(*args, **kwargs):
        received.update(kwargs)
        return {"status": "completed"}

    monkeypatch.setattr(cli, "run_pipeline", pipeline)
    assert cli.main(["run", str(tmp_path / "input"), str(tmp_path / "output"), "--tools", "identify"]) == 0
    assert received["runner_config"].tool_options == options
    assert load_tool_options_environment("detection") == {}


@pytest.mark.parametrize("raw", [
    "null", "[]", "NaN", "{" + " " * 8192 + "}",
    '{"crispridentify":{"--result_folder":"/tmp/x"}}',
    '{"crispridentify":{"--min_len_rep":50,"--max_len_rep":20}}',
    '{"crispr_tracrrna":{"--gene-mode":"meta"}}',
])
def test_transport_is_bounded_and_rejects_wrong_workflow(monkeypatch, raw):
    monkeypatch.setenv(TOOL_OPTIONS_ENV, raw)
    with pytest.raises(ValueError):
        load_tool_options_environment("detection")


def test_context_only_flags_cannot_silently_apply_to_isolated_repeats():
    with pytest.raises(ValueError, match="does not apply"):
        validate_tool_options({"crisprrepeat": {"--span": 60}}, "repeats")


@pytest.mark.parametrize(("selection", "options", "expected"), [
    ("casandra", {"casandra": {"--gene-mode": "meta"}}, ["--gene-mode", "meta"]),
    ("tracrrna", {"crispr_tracrrna": {"--merge-gap": 12}}, ["--merge-gap", "12"]),
    ("loci", {"crispr_tracrrna": {"--model-rescue": False, "--similarity": 0.9}}, ["--no-model-rescue", "--similarity", "0.9"]),
])
def test_annotation_builders_apply_flags_to_actual_stage(monkeypatch, tmp_path, selection, options, expected):
    from crispr_integration.annotations import run_annotation_pipeline
    import crispr_integration.runner as runner

    class CommandCaptured(BaseException):
        pass

    captured = []

    def capture(stage, command, *args, **kwargs):
        captured.extend(command)
        raise CommandCaptured

    source = tmp_path / "input.fasta"
    source.write_text(">sample\nACGT\n")
    config = RunnerConfig(("identify",), ("spacerplacer",), ("evor",), tool_options=options)
    monkeypatch.setattr(runner, "_run_process", capture)
    with pytest.raises(CommandCaptured):
        run_annotation_pipeline(source, tmp_path / "output", runner_config=config, selection=selection)
    assert captured[-len(expected):] == expected
