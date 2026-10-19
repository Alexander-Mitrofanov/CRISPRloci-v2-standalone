"""Independent controls for source grouping, numeric parity and native evidence."""

from itertools import permutations

import numpy as np
import pytest

from phigaro.bounded import PARAMETERS, parse_genes, read_fasta, read_hits, regions_for
from phigaro.finder.v2 import calc_scores, scan_phages


def reference_scores(values, window):
    half = window // 2
    padded = np.concatenate([np.zeros(half), values, np.zeros(half)])
    positions = np.arange(2 * half + 1) / (2 * half)
    kernel = 1 - abs(positions - 0.5) / 0.5
    return [np.sum(padded[i : i + len(kernel)] * kernel) for i in range(len(values))]


@pytest.mark.parametrize("window", [2, 3, 32, 33, 64, 65])
@pytest.mark.parametrize("size", [0, 1, 2, 31, 32, 33, 10000])
def test_triangular_scores_are_bitwise_identical(window, size):
    values = np.random.default_rng(size).normal(size=size)
    before = values.copy()
    assert np.array_equal(calc_scores(values, window), reference_scores(values, window))
    assert np.array_equal(values, before)


def test_mixed_scale_and_nonfinite_score_reductions_match():
    values = [1e100, -1e100, 1e-100, 0.7, -2.2, np.nan, np.inf] * 8
    with np.errstate(invalid="ignore"):
        assert np.array_equal(
            calc_scores(values, 32), reference_scores(values, 32), equal_nan=True
        )


def test_peak_boundaries_and_separate_runs():
    assert scan_phages([45.39, 46.0, 45.39, 45.38, 46.0], 45.39, 46.0) == [
        (0, 2),
        (4, 4),
    ]
    assert scan_phages([45.39] * 4, 45.39, 46.0) == []


def test_region_grouping_cannot_merge_record_boundaries():
    # Individually short supported runs cannot gain evidence from another record.
    genes = [
        dict(
            id=f"{source}_{i}",
            source_id=source,
            start=i * 9 + 1,
            end=(i + 1) * 9,
            gc=0.5,
            best_hit=dict(profile="VOG0001", evalue=0, score=100),
        )
        for source in ("a", "b")
        for i in range(8)
    ]
    assert regions_for({"a": "A" * 72, "b": "A" * 72}, genes[::-1]) == []


def test_source_and_gene_order_are_deterministic_with_overlapping_genes():
    genes = [
        dict(
            id=f"{source}_{i}",
            source_id=source,
            start=i * 9 + 1,
            end=(i + 1) * 9,
            gc=0.5,
            best_hit=dict(
                profile="VOG0001", evalue=PARAMETERS["max_evalue"], score=100
            ),
        )
        for source in ("a", "b")
        for i in range(50)
    ]
    genes[30]["end"] = 450
    sequences = {"b": "C" * 450, "a": "A" * 450, "unassessed": "T"}
    regions = regions_for(sequences, genes[::-1])
    assert regions == regions_for(sequences, genes)
    assert [r["source_id"] for r in regions] == ["b", "a"]
    assert [r["id"] for r in regions] == ["region000001", "region000002"]
    assert all(
        r["end"] == max(g["end"] for g in genes if g["id"] in r["gene_ids"])
        for r in regions
    )
    assert all(r["supporting_gene_ids"] == r["gene_ids"] for r in regions)
    assert all(
        r["sequence"] == sequences[r["source_id"]][r["start"] - 1 : r["end"]]
        for r in regions
    )


def test_fasta_wrapping_preserves_normalized_sequence_and_record_order(tmp_path):
    sequence = "acgtryswkmbdhvn" * 3000
    source = tmp_path / "wrapped.fa"
    source.write_text(
        ">s_one comment\n" + "\n".join(sequence) + "\n>s:two\nacgt\n", encoding="ascii"
    )
    assert list(read_fasta(source).items()) == [
        ("s_one", sequence.upper()),
        ("s:two", "ACGT"),
    ]


@pytest.mark.parametrize(
    "text", ["".join(f">s{i}\nA\n" for i in range(21)), ">s\n" + "A" * 2_000_001]
)
def test_fasta_record_and_base_bounds(tmp_path, text):
    source = tmp_path / "large.fa"
    source.write_text(text, encoding="ascii")
    with pytest.raises(ValueError, match="limited to 20 records"):
        read_fasta(source)


@pytest.mark.parametrize(
    "text",
    ["MK*\n", "MK*\n>record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.111\nMK*\n"],
)
def test_orphan_native_protein_is_rejected(tmp_path, text):
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(text, encoding="ascii")
    with pytest.raises(ValueError, match="before a Prodigal header"):
        parse_genes(proteins, {"s": "ATGAAATAA"}, {"record000001": "s"})


def test_duplicate_native_gene_is_rejected(tmp_path):
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(
        ">record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.111\nMK*\n" * 2,
        encoding="ascii",
    )
    with pytest.raises(ValueError, match="duplicate predicted genes"):
        parse_genes(proteins, {"s": "ATGAAATAA"}, {"record000001": "s"})


def test_hmmer_best_hit_ties_are_independent_of_native_row_order(tmp_path):
    hits = tmp_path / "hits.tsv"
    rows = [
        ("VOG0002", "1e-20", "100"),
        ("VOG0001", "1e-20", "100"),
        ("VOG0003", "1e-20", "90"),
        ("VOG0004", "1e-10", "200"),
    ]
    for order in permutations(rows):
        gene = dict(id="g", best_hit=None)
        hits.write_text(
            "".join(
                f"g - {profile} - {evalue} {score} 0 1 1 0 1 1 1 1 1 1 1 1 desc\n"
                for profile, evalue, score in order
            )
        )
        read_hits(hits, [gene], {row[0] for row in rows})
        assert gene["best_hit"] == dict(profile="VOG0001", evalue=1e-20, score=100.0)


def test_input_hash_binds_the_exact_parsed_fasta_snapshot(tmp_path, monkeypatch):
    import hashlib
    from types import SimpleNamespace
    from phigaro import bounded

    source, models = tmp_path / 'input.fa', tmp_path / 'models.hmm'
    original = b'>original\nacgt\n'
    source.write_bytes(original)
    models.write_text('NAME VOG0001\n')
    model_digest = hashlib.sha256(models.read_bytes()).hexdigest()
    original_digest = bounded.file_digest

    def change_input_after_parsing(path):
        if path == models:
            source.write_bytes(b'>changed\nTTTTTT\n')
        return original_digest(path)

    monkeypatch.setattr(bounded, 'file_digest', change_input_after_parsing)
    result = bounded.predict(SimpleNamespace(input=source, models_hmm=models,
        models_sha256=model_digest, output=tmp_path / 'result'))
    assert result['input_sha256'] == hashlib.sha256(original).hexdigest()
    assert result['sources'][0]['id'] == 'original'
    assert result['sources'][0]['sha256'] == hashlib.sha256(b'ACGT').hexdigest()


def test_fasta_snapshot_preserves_universal_newlines_and_rejects_embedded_controls(tmp_path):
    import hashlib
    from phigaro.bounded import read_fasta_snapshot

    source = tmp_path / 'input.fa'
    payload = b'>source\r\nac\rgt\n'
    source.write_bytes(payload)
    assert read_fasta_snapshot(source) == ({'source': 'ACGT'}, hashlib.sha256(payload).hexdigest())
    source.write_bytes(b'>source\nAC\vGT\n')
    with pytest.raises(ValueError):
        read_fasta_snapshot(source)
