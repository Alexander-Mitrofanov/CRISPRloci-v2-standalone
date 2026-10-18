import json
from pathlib import Path

import pytest

from crispr_integration.parser import ReportFormatError, load_identified_arrays


FIXTURES = Path(__file__).parent / "fixtures"
SOURCE_SHA256 = "5c6a650ff081787232ffc2057ac5b923e9337eb04b2653eb4a7d702ea3d20ba7"


def _canonical_document(arrays):
    return {
        "schema": {"name": "CRISPRidentify-report", "version": "1.1.0"},
        "source": {
            "id": "x",
            "length": 12,
            "sha256": SOURCE_SHA256,
            "sha256_scope": "uppercase_sequence",
        },
        "arrays": arrays,
    }


def test_current_canonical_report_contract():
    arrays = load_identified_arrays(FIXTURES / "current_report.json")

    assert len(arrays) == 3
    assert arrays[0].source_id == "genome_A"
    assert arrays[0].sequence_orientation == "source"
    assert arrays[0].strand == "+"
    assert arrays[1].strand == "-"
    assert arrays[2].strand is None


def test_supported_same_major_aliases_and_oriented_sequences_are_not_guessed():
    arrays = load_identified_arrays(FIXTURES / "future_report.json")

    assert [array.category for array in arrays] == ["Bona-fide", "Possible"]
    assert all(array.sequence_orientation == "biological" for array in arrays)
    assert arrays[0].strand is None


def test_unknown_canonical_schema_major_is_rejected(tmp_path):
    document = json.loads((FIXTURES / "current_report.json").read_text(encoding="utf-8"))
    document["schema"]["version"] = "2.0.0"
    report = tmp_path / "report.json"
    report.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ReportFormatError, match="unsupported CRISPRidentify-report schema major 2"):
        load_identified_arrays(report)


def test_canonical_shape_without_schema_is_rejected(tmp_path):
    document = _canonical_document([])
    del document["schema"]
    report = tmp_path / "report.json"
    report.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ReportFormatError, match="missing a schema object"):
        load_identified_arrays(report)


def test_canonical_source_provenance_is_required(tmp_path):
    document = _canonical_document([])
    del document["source"]["sha256_scope"]
    report = tmp_path / "report.json"
    report.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(ReportFormatError, match="source.sha256_scope"):
        load_identified_arrays(report)


def test_report_is_bound_to_uppercase_normalized_fasta(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">genome_A example\nacgtacgtacgt\n", encoding="utf-8")

    arrays = load_identified_arrays(
        FIXTURES / "current_report.json",
        expected_fasta_path=fasta,
    )

    assert len(arrays) == 3


def test_stale_same_id_report_hash_is_rejected(tmp_path):
    document = json.loads((FIXTURES / "current_report.json").read_text(encoding="utf-8"))
    document["source"]["sha256"] = "0" * 64
    report = tmp_path / "report.json"
    report.write_text(json.dumps(document), encoding="utf-8")
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">genome_A\nACGTACGTACGT\n", encoding="utf-8")

    with pytest.raises(ReportFormatError, match="SHA-256 differs from submitted FASTA"):
        load_identified_arrays(report, expected_fasta_path=fasta)


def test_misassociated_report_source_id_is_rejected(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">different_genome\nACGTACGTACGT\n", encoding="utf-8")

    with pytest.raises(ReportFormatError, match="does not match any submitted FASTA record"):
        load_identified_arrays(
            FIXTURES / "current_report.json",
            expected_fasta_path=fasta,
        )


def test_report_source_length_mismatch_is_rejected(tmp_path):
    document = json.loads((FIXTURES / "current_report.json").read_text(encoding="utf-8"))
    document["source"]["length"] = 11
    report = tmp_path / "report.json"
    report.write_text(json.dumps(document), encoding="utf-8")
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">genome_A\nACGTACGTACGT\n", encoding="utf-8")

    with pytest.raises(ReportFormatError, match="length 11 != submitted"):
        load_identified_arrays(report, expected_fasta_path=fasta)


def test_missing_report_for_submitted_multifasta_record_is_rejected(tmp_path):
    fasta = tmp_path / "input.fasta"
    fasta.write_text(
        ">genome_A\nACGTACGTACGT\n>genome_B\nACGTACGTACGA\n",
        encoding="utf-8",
    )

    with pytest.raises(ReportFormatError, match="missing submitted FASTA source.*genome_B"):
        load_identified_arrays(
            FIXTURES / "current_report.json",
            expected_fasta_path=fasta,
        )


def test_legacy_complete_summary_and_empty_deletion_spacer():
    arrays = load_identified_arrays(FIXTURES / "legacy_complete.json")

    assert [array.source_id for array in arrays] == ["legacy_A", "legacy_B"]
    assert arrays[1].strand == "-"
    assert arrays[1].omitted_empty_spacers == 1
    assert len(arrays[1].spacers) == 2


def test_malformed_spacers_fail_instead_of_falling_back_to_display_output(tmp_path):
    report = tmp_path / "report.json"
    report.write_text(json.dumps(_canonical_document([{
            "id": "a", "category": "Bona-fide", "consensus_repeat": "ACGT", "spacers": "not-a-list"
        }])), encoding="utf-8")

    with pytest.raises(ReportFormatError, match="spacers must be a list"):
        load_identified_arrays(report)


def test_array_with_only_deleted_spacers_is_preserved_for_adapter_filtering(tmp_path):
    report = tmp_path / "report.json"
    report.write_text(json.dumps(_canonical_document([{
            "id": "a", "category": "Low score", "consensus_repeat": "ACGT", "spacers": [{"sequence": ""}]
        }])), encoding="utf-8")

    arrays = load_identified_arrays(report)
    assert arrays[0].spacers == ()
    assert arrays[0].omitted_empty_spacers == 1
