import os
from pathlib import Path
import sys

import pytest

from crisprspacer.blast_process import blast_environment, run_blast
from crisprspacer.io import InputError, LimitError


def test_environment_drops_inherited_loaders_and_disables_network_features(tmp_path, monkeypatch):
    monkeypatch.setenv("LD_PRELOAD", "/private/injected.so")
    monkeypatch.setenv("BLASTDB", "/private/db")
    env = blast_environment(tmp_path)
    assert "LD_PRELOAD" not in env and "BLASTDB" not in env
    assert env["BLAST_USAGE_REPORT"] == "false"
    assert env["NCBI_CONFIG__BLAST__DATA_LOADERS"] == "none"
    assert env["NCBI_DONT_USE_NCBIRC"] == "1"


def test_native_failure_and_atomic_stream_paths(tmp_path):
    with pytest.raises(InputError, match="status 7"):
        run_blast([sys.executable, "-c", "raise SystemExit(7)"], tmp_path, tmp_path/"out", tmp_path/"err", 2)
    with pytest.raises(InputError, match="new files"):
        run_blast([sys.executable, "-c", "pass"], tmp_path, tmp_path/"out", tmp_path/"err", 2)


@pytest.mark.skipif(os.name != "posix", reason="Process group containment is POSIX-specific")
def test_deadline_stops_child_process_group(tmp_path):
    script = "import subprocess,sys,time; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)']); print(p.pid,flush=True); time.sleep(30)"
    with pytest.raises(LimitError, match="deadline"):
        run_blast([sys.executable, "-c", script], tmp_path, tmp_path/"out", tmp_path/"err", 1)
    pid = int((tmp_path/"out").read_text())
    proc = Path(f"/proc/{pid}/stat")
    if proc.exists():
        assert proc.read_text().split()[2] == "Z"  # terminated, awaiting init's reaper


def test_output_flood_is_a_resource_error(tmp_path):
    script = "import sys; sys.stdout.write('x'*1000000);sys.stdout.flush()"
    with pytest.raises(LimitError, match="byte budget"):
        run_blast([sys.executable, "-c", script], tmp_path, tmp_path/"out", tmp_path/"err", 5, 1000)
