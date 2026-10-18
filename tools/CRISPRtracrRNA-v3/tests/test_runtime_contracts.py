"""Input boundaries and process provenance stay explicit and reproducible."""

import hashlib
import json
import os
import sys
from pathlib import Path

import pytest

from crisprtracrrna import cli
from crisprtracrrna.config import Settings, load_settings
from crisprtracrrna.errors import ToolError, TracrError
from crisprtracrrna.models import Interval
from crisprtracrrna.process import Runner, _stderr_tail
from crisprtracrrna.sequence import subsequence


@pytest.mark.parametrize("name", ["timeout", "similarity", "coverage", "model_rescue_evalue"])
@pytest.mark.parametrize("value", [True, False, None, "0.5", float("nan"), float("inf")])
def test_numeric_settings_reject_invalid_values(name, value):
    with pytest.raises(TracrError, match=name):
        Settings(**{name: value})


@pytest.mark.parametrize(
    "name", ["threads", "tail_length", "intarna_batch_size", "upstream_length", "merge_gap"]
)
@pytest.mark.parametrize("value", [True, False, 1.0, None, "1"])
def test_integer_settings_reject_coercible_values(name, value):
    with pytest.raises(TracrError, match=name):
        Settings(**{name: value})


@pytest.mark.parametrize("weights", [None, [], "array", 4, {"array": 0.5}])
def test_weights_require_the_complete_mapping(weights):
    with pytest.raises(TracrError, match="weights"):
        Settings(weights=weights)


@pytest.mark.parametrize("value", [True, False, None, "0.5", float("nan"), float("inf"), -1])
def test_weights_reject_invalid_numbers(value):
    weights = {**Settings().weights, "array": value}
    with pytest.raises(TracrError, match="weights"):
        Settings(weights=weights)


@pytest.mark.parametrize("name", ["cmscan", "fasta36", "blastn", "intarna", "erpin", "erpin_model"])
@pytest.mark.parametrize("value", ['""', "4", "true", "[]"])
def test_config_path_shapes_fail_before_resolution(tmp_path, name, value):
    config = tmp_path / "settings.toml"
    config.write_text(f"{name} = {value}\n")
    with pytest.raises(TracrError, match=name):
        load_settings(config, {})


@pytest.mark.parametrize("name", ["models", "identify_command", "casandra_command"])
@pytest.mark.parametrize("value", ['[""]', '["valid", 4]', '"tool"', "4"])
def test_config_arrays_reject_malformed_items(tmp_path, name, value):
    config = tmp_path / "settings.toml"
    config.write_text(f"{name} = {value}\n")
    with pytest.raises(TracrError, match=name):
        load_settings(config, {})


def test_config_paths_overrides_literal_arguments_and_redaction(tmp_path):
    interpreter = tmp_path / "venv/bin/python"
    interpreter.parent.mkdir(parents=True)
    interpreter.symlink_to(sys.executable)
    config = tmp_path / "settings.toml"
    config.write_text(
        'threads = 2\nmodels = ["data/custom.cm"]\nerpin_model = "data/custom.epn"\n'
        'cmscan = "venv/bin/python"\n'
        'identify_command = ["venv/bin/python", "scripts/identify.py", "--flag", '
        '"--output=data/results", "--data/path", "$(literal)"]\n'
        '[environment]\nTOKEN = "private-value"\n'
    )
    settings = load_settings(config, {"threads": 4, "coverage": None, "similarity": 0.9})
    assert settings.threads == 4
    assert settings.similarity == 0.9
    assert settings.coverage == Settings().coverage
    assert settings.models == (str(tmp_path / "data/custom.cm"),)
    assert settings.erpin_model == str(tmp_path / "data/custom.epn")
    assert settings.cmscan == str(interpreter)
    assert settings.identify_command == (
        str(interpreter),
        str(tmp_path / "scripts/identify.py"),
        "--flag",
        "--output=data/results",
        "--data/path",
        "$(literal)",
    )
    assert settings.public()["environment"] == ["TOKEN"]
    assert "private-value" not in json.dumps(settings.public())


@pytest.mark.parametrize("factory", [Interval, Interval.from_public])
@pytest.mark.parametrize("start,end", [(0.0, 2), (1, 2.0), (False, 2), (1, True), ("1", 2)])
def test_coordinates_reject_nonintegers_before_arithmetic(factory, start, end):
    with pytest.raises(ValueError, match="integers"):
        factory(start, end)


@pytest.mark.parametrize("strand", ["?", "", None, 0])
def test_subsequence_rejects_unknown_strands(strand):
    with pytest.raises(TracrError, match="strand"):
        subsequence("ACGTA", Interval(1, 4), strand)


def test_valid_coordinate_conversions_and_strands_are_unchanged():
    interval = Interval.from_public(2, 4)
    assert interval == Interval(1, 4)
    assert interval.public() == {"start": 2, "end": 4, "length": 3}
    assert subsequence("ACGTA", interval, "+") == "CGT"
    assert subsequence("ACGTA", interval, "-") == "ACG"


@pytest.mark.parametrize("command", [[], (), "echo", None, [""], [4], ["echo", None]])
def test_runner_rejects_malformed_argv(tmp_path, command):
    with pytest.raises(ToolError, match="Command"):
        Runner(tmp_path).run(command)


@pytest.mark.parametrize(
    "arguments,expected",
    [(["-c", ""], ""), (["-c", "import sys; print(repr(sys.argv[1]))", ""], "''\n")],
)
def test_runner_preserves_valid_empty_argument_strings(tmp_path, arguments, expected):
    runner = Runner(tmp_path)
    assert runner.run([sys.executable, *arguments]).read_text() == expected
    assert runner.records[0]["argv"][1:] == arguments


def _tool(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"#!/bin/sh\nprintf '%s\\n' '{text}'\n")
    path.chmod(0o755)
    return path


@pytest.mark.parametrize("lookup", ["explicit", "relative_path", "empty_path"])
def test_runner_resolves_and_hashes_the_child_executable(tmp_path, monkeypatch, lookup):
    caller = tmp_path / "caller"
    child = caller / "child"
    relative = "bin/tool" if lookup == "relative_path" else "tool"
    _tool(caller / relative, "wrong caller executable")
    expected = _tool(child / relative, "child executable")
    monkeypatch.chdir(caller)
    runner = Runner(
        tmp_path / "logs", environment={"PATH": "bin" if lookup == "relative_path" else ""}
    )
    command = "./tool" if lookup == "explicit" else "tool"
    output = runner.run([command], cwd=Path("child"))
    assert output.read_text().strip() == "child executable"
    record = runner.records[0]
    assert record["argv"] == [str(expected)]
    assert record["cwd"] == str(child)
    assert record["executable_sha256"] == hashlib.sha256(expected.read_bytes()).hexdigest()
    assert json.loads((runner.directory / "commands.json").read_text())[0] == record


def test_runner_preserves_path_search_order_and_default_child_cwd(tmp_path):
    first = _tool(tmp_path / "bin/tool", "first")
    _tool(tmp_path / "tool", "second")
    runner = Runner(tmp_path, environment={"PATH": f"bin{os.pathsep}"})
    assert runner.run(("tool",)).read_text().strip() == "first"
    assert runner.records[0]["argv"] == [str(first)]


def test_runner_keeps_the_executable_symlink(tmp_path):
    interpreter = tmp_path / "venv/bin/python"
    interpreter.parent.mkdir(parents=True)
    interpreter.symlink_to(sys.executable)
    runner = Runner(tmp_path / "logs")
    output = runner.run([str(interpreter), "-c", "import sys; print(sys.executable)"])
    assert output.read_text().strip() == str(interpreter)
    assert runner.records[0]["argv"][0] == str(interpreter)
    assert str(interpreter) != str(interpreter.resolve())


def test_failed_process_keeps_full_multibyte_evidence_and_bounded_diagnostic(tmp_path):
    runner = Runner(tmp_path)
    diagnostic = "\U0001f9ec" * 20_000 + "FINAL_MARKER"
    script = (
        "import sys; sys.stderr.buffer.write("
        "('\\U0001f9ec' * 20000 + 'FINAL_MARKER').encode('utf-8')); sys.exit(9)"
    )
    with pytest.raises(ToolError, match="exited 9") as failure:
        runner.run([sys.executable, "-c", script])
    assert str(failure.value).split("\n", 1)[1] == diagnostic[-3000:]
    assert (tmp_path / "0001.stderr.log").read_text() == diagnostic
    assert runner.records[0]["returncode"] == 9


def test_stderr_tail_reads_a_bounded_window_at_a_utf8_boundary(tmp_path, monkeypatch):
    path = tmp_path / "stderr"
    diagnostic = "\U0001f9ec" * 8000 + "FINAL"
    path.write_text(diagnostic, encoding="utf-8")
    original_open = Path.open
    reads = []

    class BoundedReader:
        def __init__(self):
            self.handle = original_open(path, "rb")

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.handle.close()

        def seek(self, *args):
            return self.handle.seek(*args)

        def read(self, size=-1):
            assert 0 <= size <= 12000
            reads.append(size)
            return self.handle.read(size)

    monkeypatch.setattr(Path, "open", lambda *args, **kwargs: BoundedReader())
    assert _stderr_tail(path) == diagnostic[-3000:]
    assert reads == [12000]


def test_cli_legacy_aliases_and_explicit_booleans(tmp_path, monkeypatch):
    captured = []

    def run(source, output, settings):
        captured.append((source, output, settings))
        return {"summary": {"sequence_count": 1, "prediction_count": 0}}

    monkeypatch.setattr(cli, "run", run)
    source, output = tmp_path / "input.fa", tmp_path / "result"
    assert (
        cli.main(
            [
                "--file",
                str(source),
                "--output_folder",
                str(output),
                "--run_type",
                "model_run",
                "--model_type",
                "v",
                "--perform_type_v_anti_repeat_analysis",
                "False",
                "--anti_repeat_similarity_threshold",
                "0.8",
                "--anti_repeat_coverage_threshold",
                "0.75",
                "--no-model-rescue",
                "--no-mask-arrays",
            ]
        )
        == 0
    )
    assert len(captured) == 1
    actual_source, actual_output, settings = captured[0]
    assert (actual_source, actual_output) == (source, output)
    assert (settings.mode, settings.model_type) == ("model", "V")
    assert settings.type_v_motifs is False
    assert settings.model_rescue is False
    assert settings.mask_arrays is False
    assert (settings.similarity, settings.coverage) == (0.8, 0.75)


def test_cli_doctor_applies_overrides_and_bare_boolean_flags(tmp_path, monkeypatch):
    captured = []

    def doctor(settings):
        captured.append(settings)
        return {"ready": True}

    monkeypatch.setattr(cli, "doctor", doctor)
    config = tmp_path / "settings.toml"
    config.write_text('model_type = "V"\nthreads = 2\n')
    assert cli.main(["doctor", "--config", str(config), "--threads", "4", "--type-v-motifs"]) == 0
    assert captured[0].threads == 4
    assert captured[0].type_v_motifs is True
