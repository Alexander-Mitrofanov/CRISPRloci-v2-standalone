import hashlib
import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from unittest.mock import patch

from crisprtracrrna.config import Settings, load_settings
from crisprtracrrna.errors import ToolError, TracrError
from crisprtracrrna.process import Runner
from crisprtracrrna import pipeline, reporting

checks = []
with tempfile.TemporaryDirectory(prefix="tracrrna-runtime-review-") as td:
    base = Path(td)
    real = base / "actual" / "child"
    real.mkdir(parents=True)
    link = base / "linked"
    link.symlink_to(real, target_is_directory=True)
    tools = base / "actual" / "bin"
    tools.mkdir()
    script = tools / "probe"
    script.write_text("#!/bin/sh\nprintf '%s\\n' 'chosen child tool'\n")
    script.chmod(0o755)
    unusable = real / "blocked"
    unusable.mkdir()
    (unusable / "probe").write_text("not executable\n")
    child = link / ".." / "child"
    env = {**os.environ, "PATH": f"blocked:{base / 'missing'}:../bin:"}
    expected = subprocess.run(["probe"], cwd=child, env=env, capture_output=True, check=True).stdout
    runner = Runner(base / "logs", environment={"PATH": env["PATH"]})
    result = runner.run(["probe"], cwd=child)
    assert result.read_bytes() == expected == b"chosen child tool\n"
    launched = Path(runner.records[0]["argv"][0])
    assert launched.is_absolute()
    assert launched.read_bytes() == script.read_bytes()
    assert runner.records[0]["executable_sha256"] == hashlib.sha256(script.read_bytes()).hexdigest()
    assert str(link) in str(launched), "symlink identity should remain"
    checks.append("Mixed absolute/relative/empty PATH through symlink cwd and '..': native subprocess and Runner agree; SHA-256 identifies launched tool")

    timeout_runner = Runner(base / "timeout", timeout=0.1)
    child_script = "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)"
    started = time.monotonic()
    try:
        timeout_runner.run([sys.executable, "-c", child_script])
    except ToolError as exc:
        assert "timed out after 0.1s" in str(exc)
    else:
        raise AssertionError("Timeout did not fail")
    assert time.monotonic() - started < 5
    record = timeout_runner.records[0]
    assert record["timed_out"] is True
    assert json.loads((base / "timeout" / "commands.json").read_text())[0] == record
    checks.append("SIGTERM-ignoring child is killed after bounded timeout; persisted timeout record matches memory")

    out = base / "atomic"
    out.mkdir()
    destination = out / "report.json"
    old_temporary = out / "report.json.tmp"
    old_temporary.write_text("prior unrelated evidence")
    reporting.atomic_json(destination, {"writer": -1, "payload": "initial"})
    barrier = threading.Barrier(8)
    def publish(writer):
        barrier.wait()
        for iteration in range(10):
            reporting.atomic_json(destination, {"writer": writer, "iteration": iteration, "payload": str(writer) * 20000})
            actual = json.loads(destination.read_text())
            assert actual["payload"] == str(actual["writer"]) * 20000
    with ThreadPoolExecutor(max_workers=8) as executor:
        list(executor.map(publish, range(8)))
    assert sorted(p.name for p in out.iterdir()) == ["report.json", "report.json.tmp"]
    assert old_temporary.read_text() == "prior unrelated evidence"
    checks.append("80 concurrent atomic replacements across eight writers remain parseable and self-consistent, clean their own temporary files, and preserve old fixed sibling")

    source = base / "input.fa"
    source.write_text(">first\nAAAA\n>second\nCCCC\n")
    model = base / "model.cm"
    model.write_text("mock")
    settings = Settings(mode="model", models=(str(model),))
    def analyze(source, work, settings):
        work.mkdir()
        return {"predictions": [], "arrays": [], "cas": {"status": "not_requested"}}
    for marker_fails in (False, True):
        output = base / f"late-failure-{marker_fails}"
        original = OSError("second retained record cannot be published")
        def publish_record(path, value):
            if path.name == "record.json" and path.parent.name.endswith("s000002"):
                raise original
            if marker_fails and path.name == "failure.json":
                raise OSError("failure marker storage also unavailable")
            return reporting.atomic_json(path, value)
        with patch.object(pipeline, "doctor", return_value={"ready": True}), patch.object(pipeline, "analyze_record", side_effect=analyze), patch.object(pipeline, "atomic_json", side_effect=publish_record):
            try:
                pipeline.run(source, output, settings)
            except OSError as exc:
                assert exc is original
            else:
                raise AssertionError("Expected retained-record failure")
        assert len(list(output.glob("work/*/record.json"))) == 1
        assert not (output / "report.json").exists()
        if not marker_fails:
            assert json.loads((output / "failure.json").read_text())["completed_sequences"] == 1
    checks.append("Second-record write failures count exactly one retained record, never publish report.json and preserve original exception even when failure.json also fails")

    configuration = base / "config.toml"
    configuration.write_text('timeout = true\n')
    try:
        load_settings(configuration, {})
    except TracrError as exc:
        assert "timeout" in str(exc)
    else:
        raise AssertionError("Boolean timeout accepted")
    assert load_settings(configuration, {"timeout": 2}).timeout == 2
    checks.append("Boolean numeric TOML setting rejected with named TracrError; valid CLI-style numeric override still wins")

print(json.dumps({"status": "passed", "checks": checks}, indent=2))
