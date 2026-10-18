"""Regressions for the modified offline/source-bound workflow."""
from pathlib import Path

import pandas as pd
import pytest

from defense_context.pipeline import call_genes, collect, read_fasta
from defense_finder_posttreat.df_systems import build_defense_finder_systems


def test_shared_protein_does_not_duplicate_system_summary():
    rows = pd.DataFrame([
        dict(sys_id=system, type="RM", subtype="RM_II", activity="Defense",
             hit_id="shared", hit_pos=2, gene_name="profile")
        for system in ("system_a", "system_b")
    ])
    result = build_defense_finder_systems(rows)
    assert result.sys_id.to_list() == ["system_a", "system_b"]
    assert result.genes_count.to_list() == [1, 1]


def test_duplicate_ids_and_empty_sources_rejected(tmp_path):
    path = tmp_path / "bad.fna"
    for text in (">same\nACGT\n>same\nACGT\n", ">empty\n", "ACGT\n"):
        path.write_text(text)
        with pytest.raises(ValueError):
            read_fasta(path)


def test_contigs_keep_independent_order_and_partial_evidence():
    records = read_fasta(Path(__file__).parent / "data/context_fragment.fna")
    seq = next(iter(records.values()))
    genes = call_genes({"first": seq[:10000], "second": seq[:10000]})
    assert {g["replicon"] for g in genes} == {"r000001", "r000002"}
    assert sum(g["position"] == 1 for g in genes) == 2
    assert all(type(g["partial_begin"]) is bool and type(g["partial_end"]) is bool for g in genes)


def test_native_extract_suffix_is_not_lost(tmp_path):
    from defense_context.pipeline import GROUPS
    for group in GROUPS:
        folder = tmp_path / group
        (folder / "hmmer_results").mkdir(parents=True)
        (folder / "best_solution.tsv").write_text("# No Systems found\n")
    (tmp_path / GROUPS[0] / "hmmer_results/profile.res_hmm_extract").write_text(
        "# native header\ngene\treplicon\t1\t100\tprofile\t1e-12\t45\t0.9\t0.8\t1\t80\n")
    _, _, hits = collect(tmp_path, [{"id": "gene", "replicon": "replicon", "position": 1, "protein": "M" * 100}])
    assert len(hits) == 1
    assert hits[0]["profile"] == "profile"


@pytest.mark.parametrize("reverse", [False, True])
def test_ambiguous_cds_hashes_bind_exact_submitted_bases(reverse):
    """Actual Pyrodigal calls preserve IUPAC hashes and translate ambiguity as X."""
    import hashlib
    import pyrodigal
    from defense_context.pipeline import COMPLEMENT

    records = read_fasta(Path(__file__).parent / "data/context_fragment.fna")
    sequence = next(iter(records.values()))
    sequence = sequence[:300] + "RYSWKMBDHVN" + sequence[311:]
    if reverse:
        sequence = sequence.translate(COMPLEMENT)[::-1]
    native = sorted(pyrodigal.GeneFinder(meta=True).find_genes(sequence.encode("ascii")),
                    key=lambda g: (g.begin, g.end, g.strand))
    genes = call_genes({"source": sequence})
    changed = []
    for gene, raw in zip(genes, native, strict=True):
        cds = sequence[gene["start"] - 1:gene["end"]]
        if gene["strand"] == -1:
            cds = cds.translate(COMPLEMENT)[::-1]
        assert gene["cds_sha256"] == hashlib.sha256(cds.encode()).hexdigest()
        assert gene["protein"] == raw.translate(include_stop=False)
        if set(cds) - set("ACGTN"):
            changed.append(gene)
            assert gene["cds_sha256"] != hashlib.sha256(raw.sequence().encode()).hexdigest()
            assert "X" in gene["protein"]
    assert changed, "Fixture must include a called CDS traversing ambiguity"
