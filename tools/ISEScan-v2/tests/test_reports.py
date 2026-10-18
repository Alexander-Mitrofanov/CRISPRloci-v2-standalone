import csv
import io
import pytest
from isescan_v2.report import FIELDS, document, fasta

SOURCE = b">example\nAAAACCCCTTTT\n"


def row(**changes):
    data = dict(
        zip(
            FIELDS,
            [
                "seq000001",
                "IS1",
                "IS1_1",
                "1",
                "12",
                "12",
                "2",
                "1",
                "4",
                "9",
                "12",
                "8",
                "4",
                "4",
                "0",
                "0",
                "0",
                "",
                "1",
                "1e-20",
                "-1",
                "p",
                "0",
                "AAAA:AAAA",
            ],
        )
    )
    data.update(changes)
    stream = io.StringIO()
    writer = csv.DictWriter(
        stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerow(data)
    return stream.getvalue().encode()


@pytest.mark.parametrize("evalue", ["-1", "-1.0"])
def test_missing_orf_stays_missing_and_tir_is_not_tsd(evalue):
    report = document(
        SOURCE, row(**{"E-value4copy": evalue}),
        {"threads": 1, "include_partial": True}, {},
    )
    item = report["elements"][0]
    assert item["transposase_orf"] is None
    assert item["evalue_this_copy"] is None
    assert item["evalue_best_copy"] == 1e-20
    assert item["target_site_duplication"] is None
    assert item["terminal_inverted_repeats"]["identical_positions"] == 4
    assert item["native_type"] == "p" and item["contig_edge"] is True


@pytest.mark.parametrize("evalue", ["-1", "-1.0"])
def test_present_orf_rejects_missing_per_copy_evalue(evalue):
    with pytest.raises(ValueError, match="Invalid numeric evidence"):
        document(SOURCE, row(
            orfBegin="5", orfEnd="8", strand="+", orfLen="4",
            **{"E-value4copy": evalue},
        ), {}, {})


@pytest.mark.parametrize("evalue", ["0", "1e-10"])
def test_present_orf_preserves_nonnegative_per_copy_evalue(evalue):
    report = document(SOURCE, row(
        orfBegin="5", orfEnd="8", strand="+", orfLen="4",
        **{"E-value4copy": evalue},
    ), {}, {})
    item = report["elements"][0]
    assert item["transposase_orf"] == {"start": 5, "end": 8, "strand": "+"}
    assert item["transposase_missing_reason"] is None
    assert item["evalue_this_copy"] == float(evalue)
    assert item["evalue_best_copy"] == 1e-20


@pytest.mark.parametrize(
    "changes",
    [
        {"isEnd": "13"},
        {"isLen": "11"},
        {"type": "inactive"},
        {"tir": "AAAT:AAAA"},
        {"irId": "3"},
        {"E-value": "nan"},
        {"E-value4copy": "-2"},
        {"orfBegin": "0", "orfEnd": "2"},
        {"start1": "0"},
        {"strand": "+"},
    ],
)
def test_rejects_forged_or_inconsistent_native_evidence(changes):
    with pytest.raises((ValueError, KeyError)):
        document(SOURCE, row(**changes), {"threads": 1, "include_partial": True}, {})


@pytest.mark.parametrize(
    "data", [b">a\nAAA\n>a\nCCC\n", b"AAA", b">a\nAAZ\n", b">a\n", b">a;touch\nAAA\n"]
)
def test_invalid_source_is_rejected(data):
    with pytest.raises(ValueError):
        fasta(data)


def test_empty_table_is_completed_zero_calls():
    report = document(
        SOURCE,
        ("\t".join(FIELDS) + "\n").encode(),
        {"threads": 1, "include_partial": True},
        {},
    )
    assert report["elements"] == []


def test_rejects_alignment_columns_with_two_gaps():
    # Ungapped sequences and reported metrics otherwise agree with source DNA.
    with pytest.raises(ValueError, match="differ from source DNA"):
        document(SOURCE, row(tir="AAAA-:AAAA-", irLen="5", nGaps="1"), {}, {})
