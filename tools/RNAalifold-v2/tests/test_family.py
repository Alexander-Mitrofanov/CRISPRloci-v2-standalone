import json
from pathlib import Path
import subprocess
import sys

import pytest

from rna_family_tool.core import columns, fasta, pair_evidence, representatives, structure_pairs

FIXTURE = Path(__file__).parent / "data/family.fasta"


def test_duplicate_sequence_weights_and_crosswalk():
    records = fasta(FIXTURE.read_bytes(), molecule="DNA")
    sources, groups = representatives(records)
    assert len(sources) == 4 and len(groups) == 3
    assert groups[0]["member_ids"] == ["homolog_a", "homolog_c"]
    assert sources[0]["representative_id"] == sources[2]["representative_id"]


@pytest.mark.parametrize("payload", [b">a\nACGU\n>b\nACGU\n", b">a\n"+b"N"*40+b"\n>b\n"+b"G"*40,
                                     b">a\n"+b"A"*40+b"\n>a\n"+b"G"*40,
                                     b">\n"+b"A"*40+b"\n>b\n"+b"G"*40])
def test_invalid_inputs_rejected(payload):
    with pytest.raises((ValueError, IndexError)):
        fasta(payload)


def test_identical_only_family_rejected():
    with pytest.raises(ValueError, match="distinct"):
        representatives({"a": "ACGU"*10, "b": "ACGU"*10})


def test_gap_aware_evidence_and_compensatory_counts():
    rows = ["GAAAC", "CAAAG", "GAAA-", "AAAAC"]
    evidence = pair_evidence(rows, 1, 5, 0.5)
    assert evidence["canonical_sequences"] == 2
    assert evidence["noncanonical_sequences"] == 1
    assert evidence["gapped_sequences"] == 1
    assert evidence["compensatory_sequence_pairs"] == 1
    assert columns(["A-", "GU"])[0] == "NU"


@pytest.mark.parametrize("value", ["((.", ")..", "[.]"])
def test_unbalanced_structures_rejected(value):
    with pytest.raises(ValueError):
        structure_pairs(value)


def test_real_native_alignment_fold_and_duplicate_weight(tmp_path):
    pytest.importorskip("RNA")
    from rna_family_tool.cli import run
    import RNA
    result = run(FIXTURE, tmp_path / "result", molecule="DNA", mafft="mafft")
    assert result["counts"] == {"sources": 4, "unique_sequences": 3, "alignment_columns": 36,
                               "mfe_pairs": 12, "retained_probabilities": 14}
    assert result["pairs"][0]["compensatory_sequence_pairs"] == 2
    model = RNA.md()
    model.temperature = 37
    compound = RNA.fold_compound([row["sequence"] for row in result["alignment"]], model)
    structure, energy = compound.mfe()
    assert result["fold"]["structure"] == structure
    assert result["fold"]["mfe_score_kcal_mol"] == pytest.approx(energy)
    assert json.loads((tmp_path / "result/report.json").read_text())["counts"] == result["counts"]
    with pytest.raises(ValueError, match="already exist"):
        run(FIXTURE, tmp_path / "result", molecule="DNA")


def test_native_unpaired_family_is_valid(tmp_path):
    pytest.importorskip("RNA")
    from rna_family_tool.cli import run
    source = tmp_path / "unpaired.fa"
    source.write_text(">a\n" + "A"*40 + "\n>b\n" + "A"*39 + "C\n")
    result = run(source, tmp_path / "unpaired", molecule="RNA", mafft="mafft")
    assert set(result["fold"]["structure"]) == {"."}
    assert result["counts"]["mfe_pairs"] == 0
    assert result["pairs"] == []


def test_command_failure_never_publishes_report(tmp_path):
    from rna_family_tool.cli import run
    with pytest.raises(subprocess.CalledProcessError):
        run(FIXTURE, tmp_path / "failed", molecule="DNA", mafft="/bin/false")
    assert not (tmp_path / "failed/report.json").exists()


def test_cli_rejects_existing_output(tmp_path):
    result = subprocess.run([sys.executable, "-m", "rna_family_tool.cli", "--input", str(FIXTURE),
                             "--output", str(tmp_path), "--molecule", "DNA"], capture_output=True)
    assert result.returncode == 2


def test_both_mafft_children_remove_ambient_helpers_and_keep_temp_paths(tmp_path, monkeypatch):
    import os
    from types import SimpleNamespace
    from rna_family_tool import cli
    monkeypatch.setenv("MAFFT_BINARIES", str(tmp_path / "untrusted-helpers"))
    monkeypatch.setenv("TMPDIR", str(tmp_path / "job-temp"))
    observed = []
    def execute(command, **kwargs):
        observed.append(command)
        assert "MAFFT_BINARIES" not in kwargs["env"]
        assert kwargs["env"]["TMPDIR"] == str(tmp_path / "job-temp")
        if "--version" in command:
            return SimpleNamespace(stdout="", stderr="v7.526")
        kwargs["stdout"].write(b">rna_000001\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n>rna_000002\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC\n")
    monkeypatch.setattr(cli.subprocess, "run", execute)
    groups = [{"id": "rna_000001", "sequence": "A" * 30}, {"id": "rna_000002", "sequence": "A" * 29 + "C"}]
    cli.align(groups, tmp_path, mafft="/pinned/mafft", threads=1, timeout_seconds=60)
    assert len(observed) == 2
    assert os.environ["MAFFT_BINARIES"] == str(tmp_path / "untrusted-helpers")
