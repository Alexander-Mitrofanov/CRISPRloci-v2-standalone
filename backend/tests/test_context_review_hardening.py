"""Independent review regressions for native CDSs and safe context exports."""

import importlib
import json
from copy import deepcopy
from pathlib import Path

import pytest

from crispr_integration.prophage_evidence import native_genes as adapter_genes
from crispr_web.execution.validation_prophage_evidence import native_genes as api_genes
from crispr_web.publication.context_text import atomic_text
from crispr_integration.source_spacer_evidence import (
    context_document as adapter_context,
)
from crispr_web.execution.validation_source_spacer_evidence import (
    context_document as api_context,
)


@pytest.mark.parametrize("parse", [adapter_genes, api_genes])
@pytest.mark.parametrize(
    "coding,protein",
    [("AAAAAATAA", "MK*"), ("ATGAAAAAA", "MKK"), ("ATGAAATAA", "MK***")],
)
def test_complete_cds_rejects_invalid_initiation_stop_or_terminal_residues(
    parse, coding, protein
):
    fasta = (
        f">record000001_1 # 1 # 9 # 1 # partial=00;gc_cont=0.0\n{protein}\n".encode()
    )
    with pytest.raises(ValueError):
        parse(fasta, b"", {"source": coding + "A" * 20000})


@pytest.mark.parametrize("parse", [adapter_genes, api_genes])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize(
    "coding,protein,p5,p3",
    [
        ("GTGAAATAA", "MK*", False, False),
        ("GGGAAATAA", "GK*", True, False),
        ("ATGAAA", "MK", False, True),
        ("ATGGCNTAA", "MX*", False, False),
        ("ATGAARTAA", "MX*", False, False),
    ],
)
def test_partial_orientation_and_alternative_initiation_remain_valid(
    parse, strand, coding, protein, p5, p3
):
    flags = [str(int(p5)), str(int(p3))]
    if strand == -1:
        coding = coding.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[
            ::-1
        ]
        flags.reverse()
    fasta = f">record000001_1 # 1 # {len(coding)} # {strand} # partial={''.join(flags)};gc_cont=0.0\n{protein}\n".encode()
    gene = parse(fasta, b"", {"source": coding + "A" * 20000})[0]
    assert gene["translation"] == protein.removesuffix("*")
    assert gene["partial_5prime"] is p5 and gene["partial_3prime"] is p3


EXPORTS = [
    ("locus_compare", "locus-comparison.json"),
    ("locus_compare", "locus-alignments.tsv"),
    ("locus_compare", "locus-proteins.faa"),
    ("prophage_context", "prophage-context.json"),
    ("prophage_context", "prophage-regions.tsv"),
    ("prophage_context", "prophage-spacer-matches.tsv"),
    ("insertion_context", "insertion-context.json"),
    ("insertion_context", "insertion-elements.tsv"),
    ("defense_context", "defense-context.json"),
    ("defense_context", "defense-systems.tsv"),
    ("rna_family", "rna-family.json"),
    ("rna_family", "rna-family-alignment.fasta"),
    ("rna_family", "rna-family-consensus.dbn"),
    ("rna_family", "rna-family-pairs.tsv"),
]


@pytest.mark.parametrize("mode,name", EXPORTS)
def test_context_exports_never_follow_precreated_destination_links(
    tmp_path, monkeypatch, mode, name
):
    module = importlib.import_module(f"crispr_web.publication.{mode}")
    doc = {
        "loci": [],
        "links": [],
        "regions": [],
        "spacer_matches": [],
        "elements": [],
        "systems": [],
        "alignment": [],
        "pairs": [],
        "fold": {"consensus": "ACGU", "structure": "...."},
    }
    monkeypatch.setattr(module, "public_report", lambda *args, **kwargs: doc)
    if hasattr(module, "_json"):
        monkeypatch.setattr(module, "_json", lambda *args, **kwargs: doc)
    root = tmp_path / "output"
    root.mkdir()
    victim = tmp_path / "protected.txt"
    victim.write_text("must remain unchanged")
    (root / name).symlink_to(victim)
    getattr(module, f"write_{mode}_exports")(root, {mode: {"available": True}})
    assert victim.read_text() == "must remain unchanged"
    assert not (root / name).is_symlink()
    assert (root / name).is_file()


def test_atomic_export_does_not_follow_temporary_link(tmp_path):
    victim = tmp_path / "protected.txt"
    victim.write_text("unchanged")
    (tmp_path / ".report.tsv.tmp").symlink_to(victim)
    with pytest.raises(FileExistsError):
        with atomic_text(tmp_path, "report.tsv") as stream:
            stream.write("new data")
    assert victim.read_text() == "unchanged"
    assert not (tmp_path / "report.tsv").exists()


def test_atomic_export_byte_budget_preserves_previous_complete_file(tmp_path):
    target = tmp_path / "report.tsv"
    target.write_text("previous complete export")
    with pytest.raises(ValueError, match="byte limit"):
        with atomic_text(tmp_path, "report.tsv", max_bytes=3) as stream:
            stream.write("éé")
    assert target.read_text() == "previous complete export"
    assert not (tmp_path / ".report.tsv.tmp").exists()


@pytest.mark.parametrize("build", [adapter_context, api_context])
def test_complete_spacer_crosswalk_preserves_native_identity_source_and_exclusions(
    build,
):
    case = json.loads(
        (Path(__file__).parent / "fixtures/prophage-spacer/fixture.json").read_text()
    )
    doc = build(case["search"], case["reports"], case["records"], case["regions"])
    assert doc["source_inventory_complete"] is True
    assert doc["source_arrays"][0]["native_id"] == "native_array_1"
    query = doc["source_spacers"][0]
    assert (
        query["id"] == "spacer-000001"
        and query["array_id"] == doc["source_arrays"][0]["id"]
    )
    assert (
        query["sequence"]
        == case["records"][query["source_id"]][query["start"] - 1 : query["end"]]
    )
    assert all(match["query_id"] == query["id"] for match in doc["spacer_matches"])
    assert doc["search_counts"]["excluded_array_matches"] == 1


def test_bounded_source_preview_prioritizes_shown_matches(tmp_path, monkeypatch):
    from crispr_web.summaries import prophage_context as module

    document = {
        "source_inventory_complete": True,
        "genes": [],
        "regions": [],
        "spacer_matches": [{"query_id": "q2"}],
        "source_spacers": [
            {"id": "q1", "array_id": "a1"},
            {"id": "q2", "array_id": "a2"},
        ],
        "source_arrays": [{"id": "a1"}, {"id": "a2"}],
        "counts": {"source_arrays": 2, "source_spacers": 2},
    }
    monkeypatch.setattr(module, "_safe_file", lambda *args: tmp_path / "report.json")
    monkeypatch.setattr(module, "public_report", lambda *args: deepcopy(document))
    preview = module.prophage_context_summary(tmp_path, max_rows=1)
    assert preview["source_spacers"] == [{"id": "q2", "array_id": "a2"}]
    assert preview["source_arrays"] == [{"id": "a2"}]
    assert preview["source_inventory_complete"] is False
    assert (
        preview["source_arrays_truncated"] is True
        and preview["source_spacers_truncated"] is True
    )
    assert preview["counts"] == document["counts"]
    assert document["source_inventory_complete"] is True


@pytest.mark.parametrize("build", [adapter_context, api_context])
def test_native_palindrome_and_other_array_exclusion_are_complete(build):
    case = json.loads(
        (
            Path(__file__).parent / "fixtures/context-review/palindrome-search.json"
        ).read_text()
    )
    doc = build(case["search"], case["reports"], case["records"], case["regions"])
    assert doc["search_counts"]["excluded_array_matches"] == 2
    assert doc["search_counts"]["candidate_matches"] == 1
    assert doc["spacer_matches"][0]["strand"] == "both"
    assert [row["status"] for row in doc["source_spacers"]] == [
        "eligible",
        "category_excluded",
    ]


@pytest.mark.parametrize("build", [adapter_context, api_context])
def test_source_spacer_policy_covers_short_sources_palindromes_and_excluded_queries(build):
    case = json.loads((Path(__file__).parent / "fixtures/context-review/palindrome-search.json").read_text())
    doc = build(case["search"], case["reports"], case["records"], case["regions"])
    policy = doc["search_policy"]
    assert max(map(len, case["records"].values())) < 20_000
    assert policy == {
        "scope": "same-record", "query_categories": ["Bona-fide", "Possible"], "query_alphabet": "ACGT", "min_query_length": 16,
        "max_query_length": 1000, "method": "exhaustive_full_length_exact", "substitutions": 0, "gaps": 0,
        "strands": "both", "palindromic_strand": "both", "array_exclusion": "overlap_with_all_reported_arrays", "array_padding": 0,
        "search_complete": True, "results_truncated": False, "pam": "not_evaluated", "targeting_activity": "not_assessed",
        "short_source_matching": "included_independently_of_prophage_eligibility",
    }


@pytest.mark.parametrize("build", [adapter_context, api_context])
def test_completed_no_query_search_still_publishes_scope(build):
    case = json.loads((Path(__file__).parent / "fixtures/prophage-spacer/fixture.json").read_text())
    case["reports"][0]["arrays"] = []
    case["search"].update(queries=[], matches=[], outcome="no_eligible_queries", coverage=[])
    case["search"]["counts"] = {"query_occurrences": 0, "eligible_queries": 0, "targets": 1, "alignment_occurrences": 0, "candidate_matches": 0, "excluded_array_matches": 0}
    doc = build(case["search"], case["reports"], case["records"], [])
    assert doc["spacer_matches"] == []
    assert doc["search_policy"]["search_complete"] is True
    assert doc["search_policy"]["targeting_activity"] == "not_assessed"
