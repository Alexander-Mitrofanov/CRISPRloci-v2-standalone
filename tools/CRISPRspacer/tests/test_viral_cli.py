"""Exercise the public operator/search boundaries with generated DNA only."""
import json
from pathlib import Path
import shutil

import pytest

from crisprspacer.cli import main


@pytest.mark.skipif(any(shutil.which(name) is None for name in ("makeblastdb", "blastdbcmd", "blastn")),
                    reason="NCBI BLAST+ required")
def test_public_viral_workflow(tmp_path, capsys):
    examples = Path(__file__).resolve().parents[1] / "examples"
    reference, result = tmp_path / "reference", tmp_path / "result"
    assert main(["viral-index", "--sources", str(examples / "synthetic-viral-targets.fasta"),
                 "--name", "synthetic", "--source-release", "test", "--source-url",
                 "https://example.org/synthetic", "--output", str(reference)]) == 0
    assert main(["validate-viral-reference", str(reference)]) == 0
    assert main(["viral-search", "--spacers", str(examples / "synthetic-viral-spacers.fasta"),
                 "--reference", str(reference), "--output", str(result)]) == 0
    assert main(["validate-viral-result", str(result), "--reference", str(reference)]) == 0
    report = json.loads((result / "result.json").read_text())
    assert report["outcome"] == "matches"
    assert len(report["candidate_viruses"]) == 3
    assert all(hit["substitutions"] == 0 for hit in report["matches"])
    capsys.readouterr()
    assert main(["viral-search", "--spacers", str(examples / "synthetic-viral-spacers.fasta"),
                 "--reference", str(reference), "--expected-reference-sha256", "0" * 64,
                 "--output", str(tmp_path / "bad")]) == 2
    assert not (tmp_path / "bad").exists()
    assert json.loads(capsys.readouterr().err)["status"] == "error"
