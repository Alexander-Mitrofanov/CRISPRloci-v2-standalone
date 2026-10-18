import csv
import json
import os
import sys
import time
from pathlib import Path

import pytest

from crisprtracrrna.cli import boolean, main
from crisprtracrrna.config import Settings
from crisprtracrrna.errors import ToolError, TracrError
from crisprtracrrna.pipeline import run, write_exports
from crisprtracrrna.process import Runner


def test_shell_metacharacters_are_literals(tmp_path):
    work = tmp_path / "path with spaces"
    runner = Runner(work)
    literal = "$(touch should_not_exist); echo unsafe"
    output = runner.run([sys.executable, "-c", "import sys; print(sys.argv[1])", literal])
    assert output.read_text().strip() == literal
    assert not (work / "should_not_exist").exists()
    assert json.loads((work / "commands.json").read_text())[0]["returncode"] == 0


def test_external_failure_is_not_a_negative_prediction(tmp_path):
    runner = Runner(tmp_path)
    with pytest.raises(ToolError, match="exited 7"):
        runner.run(
            [sys.executable, "-c", 'import sys; print("diagnostic", file=sys.stderr); sys.exit(7)']
        )
    with pytest.raises(ToolError, match="not found"):
        runner.run(["crisprtracr-deliberately-missing-binary"])


@pytest.mark.skipif(os.name != "posix", reason="POSIX process-group cancellation")
def test_timeout_kills_descendant_process(tmp_path):
    pid_file = tmp_path / "pid"
    script = (
        "import subprocess,time,pathlib; "
        'p=subprocess.Popen(["' + sys.executable + '","-c","import time;time.sleep(30)"]); '
        f"pathlib.Path({str(pid_file)!r}).write_text(str(p.pid)); time.sleep(30)"
    )
    with pytest.raises(ToolError, match="timed out"):
        Runner(tmp_path / "logs", timeout=0.5).run([sys.executable, "-c", script])
    pid = int(pid_file.read_text())
    # A killed descendant can briefly remain a zombie until the host's init reaps it.
    stat = tmp_path / "unused"
    for _ in range(20):
        stat = Path(f"/proc/{pid}/stat")
        if not stat.exists() or stat.read_text().split()[2] == "Z":
            break
        time.sleep(0.05)
    assert not stat.exists() or stat.read_text().split()[2] == "Z"


def test_exports_quote_commas_and_distinguish_sources(tmp_path):
    record = {
        "key": "s1",
        "source": {"file": "some,file.fasta", "id": "seq.1", "length": 10},
        "predictions": [
            {
                "interval": {"start": 2, "end": 4},
                "strand": "-",
                "sequence": "ACG",
                "best_evalue": 1e-8,
                "best_score": 20,
            }
        ],
    }
    write_exports(tmp_path, [record])
    with (tmp_path / "summary.csv").open() as handle:
        row = next(csv.DictReader(handle))
    assert row["source_file"] == "some,file.fasta"
    assert row["strand"] == "-"
    assert "\tncRNA\t2\t4\t.\t-\t" in (tmp_path / "predictions.gff3").read_text()
    write_exports(tmp_path, [])
    assert len((tmp_path / "summary.csv").read_text().splitlines()) == 1


def test_explicit_false_and_invalid_cli(tmp_path):
    assert boolean("False") is False
    assert main(["run", "--input", str(tmp_path / "missing"), "--threads", "0"]) == 1


def test_pipeline_failure_never_publishes_complete_report(tmp_path, monkeypatch):
    from crisprtracrrna import pipeline

    source = tmp_path / "x.fasta"
    source.write_text(">x\nACTG\n")
    monkeypatch.setattr(pipeline, "doctor", lambda settings: {"ready": True})

    def fail(*args):
        raise ToolError("failed scientific computation")

    monkeypatch.setattr(pipeline, "analyze_record", fail)
    output = tmp_path / "failed"
    with pytest.raises(ToolError):
        run(source, output, Settings())
    assert not (output / "report.json").exists()
    assert json.loads((output / "failure.json").read_text())["status"] == "failed"
    with pytest.raises(TracrError, match="overwrite"):
        run(source, output, Settings())


@pytest.mark.parametrize("timeout", [None, 0.05])
def test_command_log_failure_does_not_hide_process_failure(tmp_path, monkeypatch, caplog, timeout):
    original_write = Path.write_text

    def fail_command_log(path, *args, **kwargs):
        if path.name == "commands.json":
            raise OSError("command journal unavailable")
        return original_write(path, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", fail_command_log)
    runner = Runner(tmp_path, timeout=timeout or 10)
    script = "import sys; sys.exit(7)" if timeout is None else "import time; time.sleep(10)"
    expected = "exited 7" if timeout is None else "timed out"
    with pytest.raises(ToolError, match=expected):
        runner.run([sys.executable, "-c", script])
    assert "Could not publish command records" in caplog.text
    if timeout is None:
        assert runner.records[0]["returncode"] == 7
    else:
        assert runner.records[0]["timed_out"] is True


def test_successful_process_requires_its_command_log(tmp_path, monkeypatch):
    original_write = Path.write_text

    def fail_command_log(path, *args, **kwargs):
        if path.name == "commands.json":
            raise OSError("command journal unavailable")
        return original_write(path, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", fail_command_log)
    with pytest.raises(ToolError, match="Cannot publish command records"):
        Runner(tmp_path).run([sys.executable, "-c", ""])
