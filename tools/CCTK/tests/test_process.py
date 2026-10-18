import argparse
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest
from cctk.errors import CCTKError
from cctk.process import run
from cctk import sequence_operations as seq


def dead(pid):
    # SIGKILL delivery is asynchronous for a descendant not owned by this test.
    deadline = time.monotonic() + 1
    while time.monotonic() < deadline:
        if _dead_now(pid):
            return True
        time.sleep(0.01)
    return _dead_now(pid)


def _dead_now(pid):
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
        return text.split()[2] == "Z"
    except FileNotFoundError:
        return True


def test_shell_metacharacters_are_literal():
    result = run(
        [
            sys.executable,
            "-c",
            "import sys; print(sys.argv[1])",
            "$(echo unsafe); filename with spaces",
        ]
    )
    assert result.stdout.strip() == "$(echo unsafe); filename with spaces"


def test_external_nonzero_and_noisy_output():
    with pytest.raises(CCTKError, match="exit 3"):
        run([sys.executable, "-c", "raise SystemExit(3)"])
    with pytest.raises(CCTKError, match="output exceeds"):
        run([sys.executable, "-c", 'import os; os.write(1,b"x"*(40*1024*1024))'])


def test_timeout_cleans_descendant(tmp_path):
    path = tmp_path / "pid"
    source = 'import subprocess,sys,time; p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(90)"]); open(sys.argv[1],"w").write(str(p.pid)); time.sleep(90)'
    with pytest.raises(CCTKError, match="deadline"):
        run([sys.executable, "-c", source, path], timeout=0.25)
    assert dead(int(path.read_text()))


def test_successful_leader_does_not_leave_descendant(tmp_path):
    path = tmp_path / "pid"
    source = 'import subprocess,sys; p=subprocess.Popen([sys.executable,"-c","import time; time.sleep(90)"]); open(sys.argv[1],"w").write(str(p.pid))'
    run([sys.executable, "-c", source, path])
    assert dead(int(path.read_text()))


def test_real_cli_sigterm_cleans_external_group(tmp_path):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "a.fna").write_text(">a\nACGTACGT\n")
    pidfile = tmp_path / "child.pid"
    fake = tmp_path / "minced"
    fake.write_text(
        f'#!{sys.executable}\nimport subprocess,sys,time\np=subprocess.Popen([sys.executable,"-c","import time; time.sleep(90)"])\nopen({str(pidfile)!r},"w").write(str(p.pid))\ntime.sleep(90)\n'
    )
    fake.chmod(0o700)
    from cctk import minced

    minced.build_parser(argparse.ArgumentParser())
    # Current command uses the original operator executable option.
    result = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "cctk",
            "minced",
            "-i",
            str(genomes),
            "-o",
            str(tmp_path / "out"),
            "-m",
            "--minced-path",
            str(fake),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    try:
        until = time.monotonic() + 5
        while not pidfile.exists() and result.poll() is None and time.monotonic() < until:
            time.sleep(0.02)
        assert pidfile.exists(), result.communicate(timeout=1)[1]
        result.send_signal(signal.SIGTERM)
        result.communicate(timeout=3)
        assert result.returncode != 0
        assert dead(int(pidfile.read_text()))
    finally:
        if result.poll() is None:
            result.kill()
            result.wait()


def test_small_extraction_batches(monkeypatch):
    monkeypatch.setattr(seq, "run_blastcmd", lambda ns, db, fs, loc: (ns, ["ACGT"]))
    assert seq.pool_MP_blastdbcmd([], "db", 10, 8) == []
    assert seq.pool_MP_blastdbcmd([(0, "%s %s %s\\n", "a 1-4 plus ")], "db", 10, 8) == ["ACGT"]


def test_blast_argument_vector_and_extraction(monkeypatch, tmp_path):
    from cctk import spacerblast

    calls = []

    def fake(argv, **kwargs):
        calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0, "", "")

    monkeypatch.setattr(spacerblast, "run", fake)
    args = argparse.Namespace(
        spacer_file="query with spaces.fa",
        blast_db_path="db;literal",
        num_threads=1,
        max_target_seqs=10,
        evalue=0.001,
        other_blast_options="-word_size 7",
    )
    assert spacerblast.run_blastn(args) == []
    assert calls[0][0][2] == "query with spaces.fa"
    assert calls[0][0][4] == "db;literal"
    monkeypatch.setattr(
        seq, "run", lambda argv, **kwargs: subprocess.CompletedProcess(argv, 0, ">a\nACGT\n", "")
    )
    assert seq.run_blastcmd([0], "db path", "%s %s %s\\n", "a 1-4 plus ") == ([0], ["ACGT"])


def test_external_bad_utf8_is_domain_error():
    with pytest.raises(CCTKError, match="UTF-8"):
        run([sys.executable, "-c", "import os; os.write(1,bytes([255]))"])


def test_optional_execution_from_threads():
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(run, [sys.executable, "-c", "print('safe')"]) for _ in range(8)]
        assert [f.result(timeout=5).stdout for f in futures] == ["safe\n"] * 8
