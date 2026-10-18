import json

import pytest

from crisprrepeat.analysis import analyze
from crisprrepeat.folding import Folder, Settings, comparison, structure_pairs
from crisprrepeat.inputs import load_arrays, read_fasta, reverse_complement, sha
from crisprrepeat.references import References, align


def genomic_fixture(root, negative=False):
    repeat = "GCGCGAAAACGCGC"
    sequence = "AAA" + repeat + "CCTAGATTTGGT" + repeat + "TTT"
    positions = [(3, 3 + len(repeat)), (3 + len(repeat) + 12, len(sequence) - 3)]
    if negative:
        sequence = reverse_complement(sequence)
        positions = sorted((len(sequence) - b, len(sequence) - a) for a, b in positions)
    fasta = root / "genome.fa"
    fasta.write_text(">contig\n" + sequence + "\n")
    report = {
        "schema": {"name": "CRISPRidentify-report", "version": "1.2.0"},
        "source": {
            "id": "contig",
            "length": len(sequence),
            "sha256": sha(sequence),
            "sha256_scope": "uppercase_sequence",
        },
        "arrays": [
            {
                "id": "native",
                "category": "Bona-fide",
                "repeat_count": 2,
                "source_interval": {"start": 4, "end": len(sequence) - 3},
                "source_sequence": sequence[3:-3],
                "consensus_repeat": sequence[positions[0][0] : positions[0][1]],
                "repeats": [
                    {"index": i, "sequence": sequence[a:b], "source_interval": {"start": a + 1, "end": b}}
                    for i, (a, b) in enumerate(positions, 1)
                ],
            }
        ],
    }
    path = root / "report.json"
    path.write_text(json.dumps(report))
    return fasta, path, report


def test_whole_window_matches_global_ensemble_and_native_accessibility():
    folder = Folder(Settings(window=150, span=80))
    sequence = "GCGCGAAAACGCGC"
    isolated = folder.isolated(sequence)
    local = folder.context(sequence)
    p = {(i, j): x for i, j, x in isolated["pairs"]}
    q = {(i, j): x for i, j, x in local["pairs"]}
    assert p.keys() == q.keys()
    assert max(abs(p[key] - q[key]) for key in p) < 1e-6
    assert local["unpaired"] == pytest.approx(isolated["unpaired"], abs=1e-6)
    result = comparison(isolated, local)
    assert result["support_change"] == pytest.approx(0, abs=1e-6)


def test_native_context_competition_changes_same_repeat_support():
    folder = Folder()
    repeat = "GCGCGAAAACGCGC"
    alone = folder.isolated(repeat)
    context = folder.context(repeat + "AAAA" + "GCGCGUUUUCGCGC")
    result = comparison(alone, context)
    assert result["status"] == "completed"
    assert abs(result["support_change"]) > 0.01
    assert any(i < len(repeat) <= j for i, j, p in context["pairs"] if p > 0.01)


def test_unpaired_empty_motif_is_missing_not_perfect():
    folder = Folder()
    isolated = folder.isolated("AAAAAAAAAAAA")
    assert isolated["structure"] == "." * 12
    result = comparison(isolated, folder.context("AAAAAAAAAAAA"))
    assert result["reason"] == "no_reference_pairs"
    assert result["support_change"] is None
    assert result["context_pair_support"] is None


def test_span_counts_endpoints_and_preserves_isolated_score():
    folder = Folder(Settings(window=30, span=5))
    isolated = folder.isolated("GCGCGAAAACGCGC")
    result = comparison(isolated, folder.context("GCGCGAAAACGCGC"))
    assert result["reason"] == "motif_exceeds_context_span"
    assert result["isolated_pair_support"] is not None
    assert result["context_pair_support"] is None


def test_source_reconstruction_orientation_and_tampering(tmp_path):
    fasta, report_path, raw = genomic_fixture(tmp_path)
    report = analyze(fasta, reports=report_path)
    assert report["counts"] == {"sources": 1, "arrays": 1, "instances": 4, "contexts": 2}
    contexts = {c["id"]: c for c in report["contexts"]}
    for row in report["instances"]:
        context = contexts[row["context_id"]]
        offset = row["context_start"]
        assert context["rna_sequence"][offset : offset + len(row["rna_sequence"])] == row["rna_sequence"]
    raw["arrays"][0]["repeats"][0]["sequence"] = "T" * 14
    report_path.write_text(json.dumps(raw))
    with pytest.raises(ValueError, match="reconstruct"):
        load_arrays(report_path, read_fasta(fasta), ["Bona-fide"])


def test_reverse_source_produces_same_oriented_analysis(tmp_path):
    forward, reverse = tmp_path / "f", tmp_path / "r"
    forward.mkdir()
    reverse.mkdir()
    f, fp, _ = genomic_fixture(forward)
    r, rp, _ = genomic_fixture(reverse, True)
    a = analyze(f, reports=fp, orientation="forward")
    b = analyze(r, reports=rp, orientation="reverse")
    assert a["contexts"][0]["pairs"] == b["contexts"][0]["pairs"]
    assert a["contexts"][0]["rna_sequence"] == b["contexts"][0]["rna_sequence"]


def test_rna_dna_normalization_and_ambiguity(tmp_path):
    fasta = tmp_path / "repeat.fa"
    fasta.write_text(">rna\nGCGCGAAAACGCGC\n>unknown\nNNNNNN\n")
    result = analyze(fasta, molecule="RNA", orientation="forward")
    assert result["instances"][0]["start"] is None
    assert result["instances"][1]["isolated"]["status"] == "unsupported"
    fasta.write_text(">mixed\nACTU\n")
    with pytest.raises(ValueError):
        read_fasta(fasta, "RNA", repeats=True)


def test_reference_oracle_and_ties_preserve_origins(tmp_path):
    assert align("ACGT", "AGT")["edit_distance"] == 1
    sequence = "GCGCGAAAACGCGC"
    path = tmp_path / "references.json"
    path.write_text(
        json.dumps(
            {
                "schema": {"name": "CRISPRrepeat-reference", "version": "1.0.0"},
                "release": "synthetic-test",
                "license": "CC0",
                "source_url": "https://example.org/fixture",
                "build_date": "2026-09-27",
                "entries": [
                    {
                        "id": "fixture",
                        "sequence": sequence,
                        "sha256": sha(sequence),
                        "occurrences": [
                            {
                                "accession": "synthetic",
                                "source_url": "https://example.org/fixture",
                                "evidence": "synthetic control",
                            }
                        ],
                    }
                ],
            }
        )
    )
    hits = References(path).search(sequence)
    assert hits["hits"][0]["match"] == "exact"
    assert hits["hits"][0]["occurrences"][0]["accession"] == "synthetic"
    assert References(path).search("AAAAAAAAAAAAAAAAAAAA")["status"] == "no_hit"


def test_deterministic_and_dot_bracket():
    a = Folder().isolated("GCGCGAAAACGCGC")
    b = Folder().isolated("GCGCGAAAACGCGC")
    assert a == b
    assert a["motif_pairs"] == structure_pairs(a["structure"])
    assert len(a["unpaired"]) == len(a["structure"])
