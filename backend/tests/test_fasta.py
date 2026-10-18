from __future__ import annotations

import tracemalloc

import pytest

from crispr_web.fasta import FastaValidationError, _iter_text_lines, normalize_fasta


LINE_BOUNDARY_VALUES = (
    "\r\n",
    "\n",
    "\r",
    "\v",
    "\f",
    "\x1c",
    "\x1d",
    "\x1e",
    "\x85",
    "\u2028",
    "\u2029",
)
LINE_BOUNDARIES = [
    pytest.param(boundary, id=identifier)
    for boundary, identifier in zip(
        LINE_BOUNDARY_VALUES,
        (
            "crlf",
            "lf",
            "cr",
            "vertical-tab",
            "form-feed",
            "file-separator",
            "group-separator",
            "record-separator",
            "next-line",
            "line-separator",
            "paragraph-separator",
        ),
        strict=True,
    )
]


def normalized(text, **overrides):
    limits = {
        "max_total_bases": 100,
        "max_record_bases": 80,
        "max_records": 4,
        "max_header_characters": 100,
    }
    limits.update(overrides)
    return normalize_fasta(text, **limits)


def test_normalizes_iupac_multirecord_and_hashes_stable_input():
    result = normalized(">sample one\nacgt rysw\n>sample-2\nNNBDHV\n")
    assert result.source_ids == ("sample", "sample-2")
    assert result.total_bases == 14
    assert result.combined == ">sample\nACGTRYSW\n>sample-2\nNNBDHV\n"
    assert len(result.sha256) == 64
    assert [record.filename for record in result.records] == [
        "001_sample.fasta",
        "002_sample-2.fasta",
    ]


@pytest.mark.parametrize("identifier", ["repeat_", "repeat-", "x" * 63 + "_suffix"])
def test_normalization_is_idempotent_for_reusable_identifiers(identifier):
    first = normalized(f">{identifier}\nACGT\n")
    second = normalized(first.combined)
    assert second.combined == first.combined
    assert second.source_ids == first.source_ids


def test_rejects_unicode_before_case_normalization():
    # U+017F LATIN SMALL LETTER LONG S uppercases to ASCII S. Validation must
    # inspect the literal submitted symbol before normalizing ASCII case.
    with pytest.raises(FastaValidationError, match="non-IUPAC"):
        normalized(">confusable\nA\u017fGT\n")


def test_invalid_symbol_precedes_record_limit_after_streaming_scan():
    with pytest.raises(FastaValidationError, match="non-IUPAC"):
        normalized(">confusable\n\u017fAAAA\n", max_record_bases=1)


def test_incremental_line_scanner_matches_splitlines_boundaries():
    text = "alpha\r\nbeta\ngamma\rdelta\vepsilon\fzeta\x1ceta\x1dtheta\x1eiota\x85kappa\u2028lambda\u2029"
    assert list(_iter_text_lines(text)) == text.splitlines()


@pytest.mark.parametrize("boundary", LINE_BOUNDARIES)
def test_each_python_line_boundary_normalizes_fasta(boundary):
    result = normalized(f">isolate_A{boundary}ac{boundary}>isolate_B{boundary}gt")
    assert [record.sequence for record in result.records] == ["AC", "GT"]


def test_mixed_python_line_boundaries_normalize_one_record():
    lines = [">mixed", "A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B"]
    text = "".join(
        line
        + (LINE_BOUNDARY_VALUES[index] if index < len(LINE_BOUNDARY_VALUES) else "")
        for index, line in enumerate(lines)
    )
    result = normalized(text)
    assert result.records[0].sequence == "ACGTRYSWKMB"


def test_pathological_short_line_input_is_rejected_with_bounded_memory():
    # Eager splitlines plus a list of every sequence fragment used more than
    # 100 MiB at the production request boundary. The incremental parser must
    # stop as soon as the configured record limit is crossed.
    text = ">many-lines\n" + "AA\n" * 200_000
    tracemalloc.start()
    try:
        with pytest.raises(FastaValidationError, match="base limit"):
            normalized(
                text,
                max_total_bases=1_000,
                max_record_bases=100,
            )
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 4_000_000


@pytest.mark.parametrize(
    "text, message",
    [
        (">same\nACGT\n>same description\nACGT\n", "duplicate"),
        (">a/b\nACGT\n>a?b\nACGT\n", "collide"),
        (">bad\nACGTZ\n", "non-IUPAC"),
        (">empty\n>next\nACGT\n", "empty"),
    ],
)
def test_rejects_ambiguous_or_invalid_provenance(text, message):
    with pytest.raises(FastaValidationError, match=message):
        normalized(text)


def test_enforces_per_record_total_and_record_count_limits():
    with pytest.raises(FastaValidationError, match="base limit"):
        normalized(">a\nAAAA\n>b\nAAAA\n", max_record_bases=3)
    with pytest.raises(FastaValidationError, match="total limit"):
        normalized(">a\nAAAA\n>b\nAAAA\n", max_total_bases=7)
    with pytest.raises(FastaValidationError, match="2-record"):
        normalized(">a\nA\n>b\nA\n>c\nA\n", max_records=2)
