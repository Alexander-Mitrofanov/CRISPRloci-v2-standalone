from types import SimpleNamespace
import numpy as np
import pytest
from phigaro.bounded import (
    read_fasta,
    parse_genes,
    read_hits,
    regions_for,
    predict,
    process,
)
from phigaro.finder.v2 import scan_phages


def test_terminal_region_numpy_and_no_mutation():
    scores = np.array([0.0, 47.0, 48.0])
    before = scores.copy()
    assert scan_phages(scores, 45.39, 46) == [(1, 2)]
    assert np.array_equal(scores, before)


def test_no_peak_or_short_run():
    assert scan_phages([45.5, 45.7, 0], 45.39, 46) == []
    assert scan_phages([], 45.39, 46) == []


def test_terminal_gene_is_retained_in_interval():
    genes = [
        dict(
            id=f"g{i}",
            source_id="s",
            start=i * 9 + 1,
            end=(i + 1) * 9,
            gc=0.5,
            best_hit=dict(profile="VOG0001", evalue=0, score=100),
        )
        for i in range(50)
    ]
    regions = regions_for({"s": "A" * 450}, genes)
    assert len(regions) == 1
    assert regions[0]["end"] == max(
        g["end"] for g in genes if g["id"] in regions[0]["gene_ids"]
    )
    assert regions[0]["sequence"] == "A" * (regions[0]["end"] - regions[0]["start"] + 1)


@pytest.mark.parametrize(
    "text", [">x\nACGT\n>x\nACGT\n", "ACGT\n", ">x\nAUCG\n", ">\nACGT\n", ">x\n"]
)
def test_bad_fasta_rejected(tmp_path, text):
    p = tmp_path / "in.fasta"
    p.write_text(text)
    with pytest.raises(ValueError):
        read_fasta(p)


def test_translation_and_partial_mapping(tmp_path):
    p = tmp_path / "p.faa"
    p.write_text(">record000001_1 # 1 # 9 # -1 # partial=01;gc_cont=0.222\nMK*\n")
    rows = parse_genes(p, {"source": "TTATTTCAT"}, {"record000001": "source"})
    assert rows[0]["partial_5prime"] is True and rows[0]["partial_3prime"] is False
    assert rows[0]["translation"] == "MK"
    p.write_text(p.read_text().replace("MK*", "MR*"))
    with pytest.raises(ValueError):
        parse_genes(p, {"source": "TTATTTCAT"}, {"record000001": "source"})


def test_native_gc_is_separate_from_source_gc(tmp_path):
    p = tmp_path / "p.faa"
    p.write_text(">record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.333\nMK*\n")
    row = parse_genes(p, {"s": "ATGAAATAA"}, {"record000001": "s"})[0]
    assert row["gc"] == 0.333 and row["source_gc"] == 1 / 9


@pytest.mark.parametrize("sequence", ["ATGGCNTAA", "ATGAARTAA"])
def test_native_ambiguous_codons_are_x_and_preserve_source_hash(tmp_path, sequence):
    from hashlib import sha256

    p = tmp_path / "p.faa"
    p.write_text(">record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.333\nMX*\n")
    row = parse_genes(p, {"s": sequence}, {"record000001": "s"})[0]
    assert row["translation"] == "MX"
    assert row["sequence_sha256"] == sha256(sequence.encode()).hexdigest()


def test_repeated_terminal_stop_is_rejected(tmp_path):
    p = tmp_path / "p.faa"
    p.write_text(">record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.333\nMK**\n")
    with pytest.raises(ValueError):
        parse_genes(p, {"s": "ATGAAATAA"}, {"record000001": "s"})


@pytest.mark.parametrize(
    "sequence,protein", [("AAAAAATAA", "MK"), ("ATGAAAAAA", "MKK")]
)
def test_complete_cds_requires_start_and_stop(tmp_path, sequence, protein):
    p = tmp_path / "p.faa"
    p.write_text(
        f">record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.333\n{protein}*\n"
    )
    with pytest.raises(ValueError):
        parse_genes(p, {"s": sequence}, {"record000001": "s"})


def test_hmm_unknown_and_nonfinite_rejected(tmp_path):
    p = tmp_path / "h.tsv"
    gene = dict(id="g", best_hit=None)
    for name, score in [("bad", "1"), ("VOG0001", "nan")]:
        p.write_text(f"g - {name} - 1e-20 {score} 0 1 1 0 1 1 1 1 1 1 1 1 desc\n")
        with pytest.raises(ValueError):
            read_hits(p, [gene], {"VOG0001"})


def test_model_hash_mismatch_prevents_launch(tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">s\nACGT\n")
    models = tmp_path / "x.hmm"
    models.write_text("NAME VOG0001\n")
    with pytest.raises(ValueError):
        predict(
            SimpleNamespace(
                input=source,
                models_hmm=models,
                models_sha256="0" * 64,
                output=tmp_path / "out",
            )
        )
    assert not (tmp_path / "out").exists()


def test_native_command_is_shell_free_and_propagates_failure(tmp_path, monkeypatch):
    import subprocess

    def fail(command, **kwargs):
        assert (
            command == ["tool with spaces", "input with spaces"]
            and "shell" not in kwargs
            and kwargs["check"] is True
        )
        raise subprocess.CalledProcessError(7, command)

    from phigaro import bounded

    monkeypatch.setattr(bounded, "_run_native", fail)
    with pytest.raises(subprocess.CalledProcessError):
        process(["tool with spaces", "input with spaces"], tmp_path, "native", 1)
