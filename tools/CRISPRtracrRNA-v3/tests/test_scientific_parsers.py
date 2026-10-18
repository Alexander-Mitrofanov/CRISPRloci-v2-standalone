import random
import re

import pytest

from crisprtracrrna.config import Settings
from crisprtracrrna.errors import ContractError
from crisprtracrrna.evidence import best_hit, candidate_evidence, merge_model_hits, type_v_motifs
from crisprtracrrna.models import AntiRepeat, Array, Hit, Interval, Sequence
from crisprtracrrna.search import (
    make_candidates,
    parse_alignment_table,
    parse_cm_table,
    parse_erpin,
    parse_intarna,
)


def test_cm_multiple_sequences_reverse_strand_and_multiple_models(tmp_path):
    path = tmp_path / "hits.tbl"
    path.write_text(
        "# Infernal output\n"
        "m - seq1 - cm 1 40 80 41 - no 1 .5 0 22.5 1e-8 ! description\n"
        "m - seq2 - cm 1 40 3 42 + no 1 .5 0 18.0 1e-6 ! description\n"
    )
    hits = parse_cm_table(path, {"seq1": 100, "seq2": 100}, "a.cm")
    assert hits[0].strand == "-"
    assert hits[0].interval == Interval(40, 80)
    assert hits[1].query == "seq2"
    assert best_hit([hits[1], hits[0]]) == hits[0]
    assert best_hit(hits, forward_only=True) == hits[1]
    assert best_hit([]) is None
    path.write_text("# No hits\n")
    assert parse_cm_table(path, {"seq1": 100}, "a") == []


@pytest.mark.parametrize(
    "row",
    [
        "bad",
        "m - bad - cm 1 40 1 40 + no 1 .5 0 22 1e-8 ! x",
        "m - seq1 - cm 1 40 1 400 + no 1 .5 0 22 1e-8 ! x",
        "m - seq1 - cm 1 40 1 40 + no 1 .5 0 nan 1e-8 ! x",
        "m - seq1 - cm 1 40 40 1 + no 1 .5 0 22 1e-8 ! x",
        "m - seq1 - cm 1 40 1 40 - no 1 .5 0 22 1e-8 ! x",
    ],
)
def test_cm_malformed_not_silently_no_hits(tmp_path, row):
    path = tmp_path / "bad.tbl"
    path.write_text(row)
    with pytest.raises(ContractError):
        parse_cm_table(path, {"seq1": 100}, "a")


def test_alignment_coverage_uses_query_span_and_rejects_any_array_overlap(tmp_path):
    source = Sequence("seq", "seq", "A" * 100)
    array = Array("a", Interval(40, 60), "A" * 10, "+", "Possible", 0.8)
    queries = {"repeat": (array, "A" * 10, "+")}
    path = tmp_path / "hits.tsv"
    path.write_text("repeat seq 70 15 0 0 1 6 10 15 1 5\nrepeat seq 100 10 0 0 1 10 38 47 1 5\n")
    hits = parse_alignment_table(path, queries, source, [array], Settings())
    assert len(hits) == 2  # same locus, two strand hypotheses
    assert hits[0].coverage == 0.6
    assert hits[0].identity == 0.7  # inclusive threshold
    assert hits[0].strand == "-"
    assert hits[0].interval == Interval(9, 15)
    path.write_text("seq repeat 70 15 0 0 10 15 1 6 1 5\n")
    assert parse_alignment_table(path, queries, source, [array], Settings(), True) == hits


@pytest.mark.parametrize("identity", ["100", "10"])
@pytest.mark.parametrize(
    ("column", "value"),
    [
        (3, "bad"),
        (3, "0"),
        (3, "-1"),
        (3, "9"),
        (4, "bad"),
        (4, "1.5"),
        (4, "-1"),
        (4, "11"),
        (5, "bad"),
        (5, "1.5"),
        (5, "-1"),
        (5, "11"),
        (10, "bad"),
        (10, "nan"),
        (10, "inf"),
        (10, "-1"),
        (11, "bad"),
        (11, "nan"),
        (11, "inf"),
    ],
)
def test_alignment_validates_all_numeric_fields_before_filtering(tmp_path, identity, column, value):
    source = Sequence("seq", "seq", "A" * 100)
    array = Array("a", Interval(40, 60), "A" * 10, "+", "Possible", 0.8)
    queries = {"repeat": (array, array.repeat, "+")}
    cols = f"repeat seq {identity} 10 0 0 1 10 1 10 1e-5 5".split()
    cols[column] = value
    path = tmp_path / "hits.tsv"
    path.write_text(" ".join(cols) + "\n")
    with pytest.raises(ContractError, match="Invalid anti-repeat alignment table"):
        parse_alignment_table(path, queries, source, [array], Settings())


def test_alignment_rejects_combined_counts_exceeding_length(tmp_path):
    source = Sequence("seq", "seq", "A" * 100)
    array = Array("a", Interval(40, 60), "A" * 10, "+", "Possible", 0.8)
    queries = {"repeat": (array, array.repeat, "+")}
    path = tmp_path / "hits.tsv"
    path.write_text("repeat seq 10 10 6 5 1 10 1 10 1e-5 5\n")
    with pytest.raises(ContractError, match="Invalid alignment counts"):
        parse_alignment_table(path, queries, source, [array], Settings())


def test_erpin_every_hit_not_only_first_and_score_not_evalue(tmp_path):
    path = tmp_path / "erpin.txt"
    path.write_text(
        "Training set\n>s\nFW 1 1..20 6.0 1e-2\nATG\n"
        "FW 2 21..40 12.0 1e-6\nATG\nRC 3 80..61 14.0 1e-8\n"
    )
    hits = parse_erpin(path, {"s": 100})["s"]
    assert len(hits) == 3
    assert best_hit(hits, forward_only=True).score == 12
    assert hits[2].strand == "-"
    assert hits[2].interval == Interval(60, 80)


def test_intarna_named_columns_and_empty_output(tmp_path):
    path = tmp_path / "intarna.csv"
    path.write_text(
        "E;end2;start2;hybridDP;start1;end1;extra\n-10;8;1;(((.....&.....)));2;9;future\n"
    )
    result = parse_intarna(path, 10, 10)
    assert result["energy_kcal_mol"] == -10
    assert result["target_interval"] == {"start": 2, "end": 9, "length": 8}
    path.write_text("start1;end1;start2;end2;E;hybridDP\n")
    assert parse_intarna(path, 10, 10) is None
    path.write_text("unrecognized\n")
    with pytest.raises(ContractError):
        parse_intarna(path, 10, 10)


def test_missing_evidence_gets_no_consistency_bonus_and_significance_ranks_correctly():
    source = Sequence("s", "s", "ACGT" * 100)
    array = Array("a", Interval(200, 220), "ACGT", None, "Possible", 0.8)
    anti = AntiRepeat(array, "ACGT", "+", Interval(10, 14), "+", 0.9, 0.8, True)
    candidate = make_candidates(source, [anti], Settings())[0]
    cas = {"cassettes": []}
    empty = candidate_evidence(candidate, source, [], [], None, cas, Settings())
    for feature in [
        "interaction",
        "model",
        "terminator",
        "orientation",
        "anti_tail",
        "tail_terminator",
    ]:
        assert empty["ranking"]["features"][feature] == 0
    weak = Hit(candidate.id, Interval(80, 100), "+", 5, 0.1)
    strong = Hit(candidate.id, Interval(80, 100), "+", 15, 1e-6)
    weak_result = candidate_evidence(candidate, source, [], [weak], None, cas, Settings())
    strong_result = candidate_evidence(candidate, source, [], [weak, strong], None, cas, Settings())
    assert strong_result["ranking"]["score"] > weak_result["ranking"]["score"]
    assert strong_result["terminator"]["evalue"] == 1e-6


def test_model_merge_preserves_strand_and_model_identity():
    source = Sequence("s", "s", "ACGT" * 100)
    hits = [
        Hit("s", Interval(10, 20), "+", 10, 0.01, "a"),
        Hit("s", Interval(15, 25), "-", 11, 0.001, "b"),
        Hit("s", Interval(19, 30), "+", 12, 0.001, "c"),
    ]
    rows = merge_model_hits(source, hits, 0)
    assert len(rows) == 2
    assert len(rows[0]["hits"]) == 2
    assert rows[1]["strand"] == "-"
    assert rows[1]["sequence"] == "TACGTACGTA"


def test_type_v_masks_distance_and_inclusive_coordinates():
    seq = "TTT" + "A" * 12 + "CTTTC" + "CCTCCAAAAAAGGAGG"
    result = type_v_motifs(seq)
    assert result["paired_motif"]["intervening_bases"] == 12
    assert result["paired_motif"]["CTTTC"] == {"start": 16, "end": 20, "length": 5}
    assert result["mask_motif"]["sequence"] == "CCTCCAAAAAAGGAGG"
    assert type_v_motifs("CC|CCAAAAAAGGAGG")["mask_motif"] is None


@pytest.mark.parametrize("gap", [11, 12, 28, 29])
def test_type_v_paired_motif_gap_boundaries(gap):
    result = type_v_motifs("TTT" + "A" * gap + "CTTTC")["paired_motif"]
    if 12 <= gap <= 28:
        assert result == {
            "TTT": Interval(0, 3).public(),
            "CTTTC": Interval(3 + gap, 8 + gap).public(),
            "intervening_bases": gap,
        }
    else:
        assert result is None


def test_type_v_pair_selection_matches_reference_with_overlaps_and_competing_matches():
    randomizer = random.Random(20261013)
    sequences = [
        "TTTT" + "A" * 11 + "CTTTC",
        "TTT" + "A" * 12 + "CTTTC" + "A" * 4 + "CTTTC",
        "TTT" + "A" * 12 + "TTT" + "A" * 12 + "CTTTC",
        "TTT" * 500 + "A" * 29 + "CTTTC",
        "CTTTC" * 100 + "TTT" * 100,
    ]
    for _ in range(200):
        sequence = "".join(randomizer.choices("ACGTN", k=200))
        for motif in ("TTT", "TTT", "CTTTC", "CTTTC"):
            position = randomizer.randrange(len(sequence) - len(motif) + 1)
            sequence = sequence[:position] + motif + sequence[position + len(motif) :]
        sequences.append(sequence)
    for sequence in sequences:
        expected = None
        for first in re.finditer("(?=TTT)", sequence):
            for second in re.finditer("(?=CTTTC)", sequence):
                gap = second.start() - first.start() - 3
                if 12 <= gap <= 28:
                    expected = {
                        "TTT": Interval(first.start(), first.start() + 3).public(),
                        "CTTTC": Interval(second.start(), second.start() + 5).public(),
                        "intervening_bases": gap,
                    }
                    break
            if expected:
                break
        assert type_v_motifs(sequence)["paired_motif"] == expected
