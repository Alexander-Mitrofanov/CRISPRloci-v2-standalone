from pathlib import Path
import hashlib
import json
import shutil

import pytest

from crispr_integration.leader import add_leader_context
from crispr_integration.runner import RunnerConfig

FIXTURE = Path(__file__).resolve().parents[2] / "backend/tests/fixtures/leader"


def test_isolated_context_commands_and_reused_identify_contract(tmp_path, monkeypatch):
    from crispr_integration import runner

    job = tmp_path / "job"
    reports = job / "crispridentify/source"
    reports.mkdir(parents=True)
    (job / "logs").mkdir()
    shutil.copyfile(FIXTURE / "identify.json", reports / "report.json")
    calls = []

    def execute(stage, command, logs, *, timeout_seconds):
        calls.append(command)
        assert command[0] == "/isolated/crisprleader"
        assert timeout_seconds == 60 and logs == job / "logs"
        if command[1] == "extract":
            out = Path(command[command.index("--output") + 1])
            out.mkdir()
            shutil.copyfile(FIXTURE / "result.json", out / "result.json")
            raw = (out / "result.json").read_bytes()
            (out / "manifest.json").write_text(
                json.dumps(
                    {
                        "schema_version": "1.0.0",
                        "files": {"result.json": {"sha256": hashlib.sha256(raw).hexdigest(), "size_bytes": len(raw)}},
                    }
                )
            )
        return {"status": "completed"}

    monkeypatch.setattr(runner, "_run_process", execute)
    config = RunnerConfig(("identify",), ("spacerplacer",), ("evor",), crisprleader=("/isolated/crisprleader",))
    result = {"workflow": "leader", "stages": {}}
    add_leader_context(job, FIXTURE / "source.fasta", result, config, 60, 37)
    assert [command[1] for command in calls] == ["adapt", "extract", "validate"]
    assert result["stages"]["crisprleader"]["context_count"] == 2
    assert result["stages"]["crisprleader"]["prediction_performed"] is False
    with pytest.raises(ValueError, match="between 1 and 5000"):
        add_leader_context(job, FIXTURE / "source.fasta", result, config, 60, 0)
