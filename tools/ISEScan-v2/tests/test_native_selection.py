"""Independent controls for native selection, source binding and execution reuse."""

from pathlib import Path
import random
import sys
from unittest.mock import patch

import pytest

NATIVE = Path(__file__).resolve().parents[1] / "isescan_v2/native"
sys.path.insert(0, str(NATIVE))
import constants  # noqa: E402
import isPredict  # noqa: E402
import pred  # noqa: E402
import tools as native_tools  # noqa: E402


def orf_hit(begin, end, strand="+", family="IS1_1", evalue=1e-20):
    return (("seq000001", begin, end, strand), family, evalue, evalue, 0,
            {"orfhits4tpase": []})


@pytest.mark.parametrize("prepare", [isPredict.prepare4hmmsearch, isPredict.prepare4phmmer])
@pytest.mark.parametrize("table,updated,retry", [
    (None, False, True), (b"", False, True), (b"#", False, True),
    (b"# [ok]", False, True), (b"header\n# [ok]\n", False, False),
    (b"header\n# [ok]\n", True, True), (b"header\n# [ok]\ntruncated", False, True),
])
def test_hmmer_reuses_only_completed_unchanged_tables(tmp_path, prepare, table, updated, retry):
    protein = tmp_path / "protein file.faa"
    protein.write_text(">seq000001_1_3_+\nM\n")
    result = tmp_path / "results" / "org" / "model file.protein file.faa"
    if table is not None:
        result.parent.mkdir(parents=True)
        result.write_bytes(table)
    args, outputs = prepare("model file", [(str(protein), "org", updated)],
                            str(tmp_path / "results"), 2)
    assert outputs == [str(result)]
    assert args == ([("model file", str(protein), str(result), 2)] if retry else [])


def test_empty_and_missing_proteomes_never_become_search_targets(tmp_path):
    empty = tmp_path / "empty.faa"
    empty.touch()
    assert isPredict.prepare4hmmsearch("m", [(str(empty), "o", False),
        (str(tmp_path / "missing"), "o", True)], str(tmp_path), 1) == ([], [])


@pytest.mark.parametrize("use_intersection", [False, True])
@pytest.mark.parametrize("seed", range(10))
def test_interval_sweep_matches_independent_all_pairs_oracle(monkeypatch, use_intersection, seed):
    monkeypatch.setattr(constants, "intersected2remove", use_intersection)
    monkeypatch.setattr(constants, "min4intersect", 3)
    monkeypatch.setattr(constants, "overlap2removeRedundancy", 0.5)
    rng = random.Random(seed)
    boundaries = [(start, start + rng.randrange(50)) for start in
                  [rng.randrange(1, 500) for _ in range(70)]]
    boundaries.extend([(1, 10), (10, 19), (20, 20), (20, 20), (1, 500)])
    expected = set()
    for i, (a, b) in enumerate(boundaries):
        for j in range(i + 1, len(boundaries)):
            c, d = boundaries[j]
            intersection = min(b, d) - max(a, c) + 1
            qualifies = intersection >= 3 if use_intersection else (
                intersection / min(b - a + 1, d - c + 1) >= 0.5)
            if qualifies:
                expected.update((i, j))
    assert pred._intersected_ids(boundaries) == expected


def test_sparse_intervals_do_not_invoke_pairwise_scientific_measure():
    with patch.object(native_tools, "chooseMeasure", wraps=native_tools.chooseMeasure) as measure:
        assert pred._intersected_ids([(i, i + 10) for i in range(1, 100001, 20)]) == set()
        assert measure.call_count == 0


def test_hit_selection_keeps_copy_priority_then_evalue_then_original_tie():
    def hit(label, copies, evalue):
        return {"label": label, "bd": (100, 199), "hmmhit": ("IS1_1", evalue),
                "occurence": {"ncopy4is": copies}}
    hits = [hit("single", 1, 1e-50), hit("first", 2, 1e-20), hit("tie", 5, 1e-20)]
    assert pred.parallel4overlappedHits(("seq000001", hits)) == [hits[1]]


def test_orf_selection_keeps_majority_family_and_original_equal_evalue_hit():
    hits = [orf_hit(100, 199, family="IS1_1", evalue=1e-50),
            orf_hit(100, 199, family="IS5_1"), orf_hit(100, 199, "-", family="IS5_1")]
    result = pred.parall4orfhits(("seq000001", hits))
    assert len(result) == 1
    assert result[0][:4] == hits[1][:4]
    assert result[0][4] == 2
    assert result[0][5] is hits[1][5]


def test_orf_coordinate_dedup_ignores_strand_and_preserves_rank():
    hits = [(0, "", "source_with_underscores_10_20_+"),
            (1, "", "source_with_underscores_10_20_-"),
            (2, "", "source_with_underscores_21_30_+")]
    assert pred.refine_hmm_hits(hits) == [hits[0], hits[2]]


def copy_pair(query="seq000001_IS1_1_1_1000_100_199_+", start=100, end=199):
    return {"qseqid": query, "sseqid": "seq000001", "qstart": 100, "qend": 199,
            "sstart": start, "send": end, "length": 100}


def copy_search(pairs, orfhits=None):
    with (patch.object(pred, "writeOrfExt2fileOnStream", return_value=">query\nAAA\n"),
          patch.object(native_tools, "doBlastn2seqOnStream", return_value=("blast table", "")),
          patch.object(native_tools, "getBlastResult4dnaOnStream", return_value=pairs)):
        return pred.getFullIS4seqOnStream(("seq000001", orfhits or [orf_hit(100, 199)],
                                         ("org", "source.fasta", "A" * 1000), 1))


def test_copy_search_moves_self_first_preserving_other_ties_and_orf_identity():
    query = orf_hit(100, 199)
    first, self_pair, last = copy_pair(start=400, end=499), copy_pair(), copy_pair(start=799, end=700)
    pairs = copy_search([first, self_pair, last], [query, orf_hit(100, 199, family="IS5_1")])
    assert next(iter(pairs.values())) == [self_pair, first, last]
    assert all(pair["orfhit"] is query for pair in pairs[self_pair["qseqid"]])


@pytest.mark.parametrize("pairs,message", [
    ([copy_pair(start=400, end=499)], "lacks self-alignment"),
    ([copy_pair(query="seq000001_IS1_1_1_1000_300_399_+")], "ORF is absent"),
    ([dict(copy_pair(), sseqid="other_source")], "subject differs"),
    ([copy_pair(query="other_source_IS1_1_1_1000_100_199_+")], "query differs"),
])
def test_copy_search_rejects_unbound_or_incomplete_evidence(pairs, message):
    with pytest.raises(RuntimeError, match=message):
        copy_search(pairs)


@pytest.mark.parametrize("fails", [False, True])
def test_blast_subject_is_closed_and_removed_on_success_or_failure(fails):
    subjects = []
    def blast(nthreads, query, subject, **kwargs):
        path = Path(subject)
        subjects.append(path)
        assert path.read_text() == ">seq000001\nAAA\n"
        if fails:
            raise RuntimeError("native failure")
        return "", ""
    with (patch.object(pred, "writeOrfExt2fileOnStream", return_value=">query\nAAA\n"),
          patch.object(native_tools, "doBlastn2seqOnStream", side_effect=blast)):
        if fails:
            with pytest.raises(RuntimeError, match="native failure"):
                pred.getFullIS4seqOnStream(("seq000001", [], ("o", "f", "AAA"), 1))
        else:
            assert pred.getFullIS4seqOnStream(("seq000001", [], ("o", "f", "AAA"), 1)) == {}
    assert subjects and all(not path.exists() for path in subjects)


def test_added_copy_inherits_first_exact_source_orf_and_has_no_invented_tpase():
    original = orf_hit(100, 199)
    hits = {"seq000001": [original, orf_hit(100, 199, family="IS5_1")]}
    pairs = {"seq000001": {"query": [copy_pair(), copy_pair(start=799, end=700)]}}
    with patch.object(pred, "removeOverlappedOrfhits", side_effect=lambda value: value):
        result = pred.addNonORFcopy(pairs, hits)
    copy = result["seq000001"][-1]
    assert copy[0] == ("seq000001", 700, 799, "-")
    assert copy[1:5] == original[1:5]
    assert copy[5] == {"orfhits4tpase": []}


def blast_row(strand="+", qlen=800):
    return (f"seq000001_IS1_1_100_899_200_399_{strand} seq000001 100 800 0 0 "
            f"1 800 100 899 0 100 800 {qlen} 1000\n")


def test_blast_query_length_mismatch_is_an_error():
    with pytest.raises(RuntimeError, match="different lengths for query sequence"):
        native_tools.getBlastResult4dnaOnStream(blast_row(qlen=801))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_blast_copy_coordinates_bind_both_query_strands_to_source(strand):
    hit = native_tools.getBlastResult4dnaOnStream(blast_row(strand))[0]
    assert (hit["qstart"], hit["qend"], hit["orfBegin"], hit["orfEnd"]) == (100, 899, 200, 399)


@pytest.mark.parametrize("strand", ["+", "-", ""])
def test_single_file_exports_original_source_tpase_and_missing_values(tmp_path, strand):
    seqid, dna = "seq000001", "AAAACCCCTTTT"
    raw = [( (seqid, 5, 8, strand), "IS1_1", 1e-20, 1e-15, 0)] if strand else []
    hit = {"orf": (seqid, 1, 12, "+"), "bd": (1, 12), "type": "p", "tirs": [],
           "hmmhit": ("IS1_1", 1e-30, 1e-30, 1, {"orfhits4tpase": raw}),
           "occurence": {"ncopy4is": 1, "sim4is": 1, "ncopy4orf": 1, "sim4orf": 1}}
    protein_id = f"{seqid}_5_8_{strand}"
    pred.outputIS4multipleSeqOneFile({seqid: [hit]}, {seqid: ("o", "f", dna)},
        {seqid: ("f", {protein_id: "MP"} if strand else {})}, {}, "o/f", str(tmp_path))
    import csv
    with (tmp_path / "o/f.tsv").open() as table:
        row = next(csv.DictReader(table, delimiter="\t"))
    assert (row["orfBegin"], row["orfEnd"], row["strand"], row["E-value4copy"]) == (
        ("5", "8", strand, "1e-15") if strand else ("0", "0", "", "-1"))
    assert row["E-value"] == "1e-30" and row["type"] == "p"
    if strand:
        assert (tmp_path / "o/f.orf.faa").read_text() == f">{protein_id}\nMP\n"
        assert (tmp_path / "o/f.orf.fna").read_text() == f">{protein_id}\n{'GGGG' if strand == '-' else 'CCCC'}\n"
    else:
        assert (tmp_path / "o/f.orf.faa").read_text() == ""
        assert (tmp_path / "o/f.orf.fna").read_text() == ""


@pytest.mark.parametrize("source,family", [
    ("source_with_underscores", "IS1_1"), ("gi|228288719|ref|NC_012624.1|", "IS200/IS605_1"),
])
def test_query_source_binding_preserves_compound_ids_and_is200_family(source, family):
    query = f"{source}_{family}_1_1000_100_199_-"
    assert pred._query_orf(query, source) == (100, 199, "-")
