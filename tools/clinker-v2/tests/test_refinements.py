"""Independent controls for scientific equivalence and refinement defects."""
import copy
import json
from collections import Counter
from itertools import product

import pytest
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from clinker.align import Globaligner, Alignment, Link, compare_pairs, compute_identity
from clinker.bounded import SCHEMA, compare
from clinker.classes import Cluster, Gene, Locus


def named_gene(label):
    return Gene(label=label, names={"locus_tag": label}, translation="MA")


def test_function_groups_merge_only_groups_with_linked_genes():
    a, b, c, d, isolated = map(named_gene, ("a", "b", "c", "d", "isolated"))
    one = Cluster("one", [Locus("one", [a, b, isolated])])
    two = Cluster("two", [Locus("two", [c, d])])
    alignment = Alignment(query=one, target=two)
    for x, y in ((a, c), (b, d), (a, d)):
        alignment.add_link(x, y, 1, 1)
    aligner = Globaligner()
    aligner.add_alignment(alignment)
    aligner.build_gene_groups(
        functions={"first": ["a"], "second": ["b"], "isolated": ["isolated"]},
        colours={"first": "#FF0000", "second": "#00FF00"},
    )
    assert [(group.label, set(group.genes)) for group in aligner.groups] == [
        ("first", {a.uid, b.uid, c.uid, d.uid}), ("isolated", {isolated.uid})
    ]
    assert aligner.groups[0].colour == "#FF0000"


def test_link_deserialization_does_not_mutate_input():
    original = Link(query=named_gene("a"), target=named_gene("b"), identity=.5, similarity=.75)
    raw = original.to_dict()
    frozen = copy.deepcopy(raw)
    assert Link.from_dict(raw).to_dict() == frozen
    assert raw == frozen
    assert Link.from_dict(raw).to_dict() == frozen


@pytest.mark.parametrize("strand", (1, -1))
def test_ranged_locus_fallback_extracts_original_source_coordinates(strand):
    record = SeqRecord(Seq("CCCATGGCTTAAAAA"), id="source")
    record.features = [SeqFeature(FeatureLocation(3, 12, strand), type="CDS",
                                  qualifiers={"locus_tag": ["g"]})]
    frozen = copy.deepcopy(record)
    locus = Locus.from_seqrecord(record, set_origin=False, start=3, end=12)
    assert locus.genes[0].sequence == str(frozen.features[0].extract(frozen.seq))
    assert locus.genes[0].translation == str(frozen.features[0].extract(frozen.seq).translate()).removesuffix("*")
    assert (locus.genes[0].start, locus.genes[0].end) == (3, 12)
    assert record.seq == frozen.seq
    assert record.features == frozen.features


def test_adjacency_multiset_count_matches_independent_oracle():
    one = [(1, 2)] * 100 + [(2, 3)] * 3 + [(3, 4)]
    two = [(1, 2)] * 50 + [(2, 3)] * 7 + [(4, 3)]
    expected = sum((Counter(one) & Counter(two)).values())
    assert compare_pairs(one, two) == expected == 53


@pytest.mark.parametrize("one,two", [
    ("M-A.AXJ", "MGA.AXZ"), ("--..", "--.."),
    *product("ACDEFGHIKLMNPQRSTVWYBXZJUO", repeat=2),
])
def test_identity_and_similarity_match_independent_residue_oracle(one, two):
    groups = ("GAVLI", "FYW", "CM", "ST", "KRH", "DENQ", "P")
    denominator = sum(not (a == b and a in "-.") for a, b in zip(one, two))
    if not denominator:
        with pytest.raises(ValueError):
            compute_identity((one, two))
        return
    matches = sum(a == b and a not in "-." for a, b in zip(one, two))
    similar = sum(a != b and any(a in group and b in group for group in groups) for a, b in zip(one, two))
    assert compute_identity((one, two)) == (matches / denominator, (matches + similar) / denominator)


def test_synteny_matrix_scores_each_unordered_pair_once(monkeypatch):
    aligner = Globaligner()
    clusters = [Cluster(str(i), []) for i in range(3)]
    aligner.add_clusters(*clusters)
    scores = {frozenset((a.uid, b.uid)): i for i, (a, b) in enumerate(
        ((clusters[0], clusters[1]), (clusters[0], clusters[2]), (clusters[1], clusters[2])), 1)}
    calls = []
    def synteny(one, two, i):
        calls.append((one, two, i))
        return scores[frozenset((one, two))] * i
    monkeypatch.setattr(aligner, "synteny", synteny)
    assert aligner.matrix(i=2).tolist() == [[0, 2, 4], [2, 0, 6], [4, 6, 0]]
    assert len(calls) == 3


def repeated_inputs(coding="ATGGCTTAA", translation="MA", count=3):
    gene = {"start": 1, "end": len(coding), "strand": "+", "phase": 0,
            "translation_table": 11, "partial_5prime": False, "partial_3prime": False,
            "translation": translation}
    document = {"schema": SCHEMA, "coordinate_system": "1-based-inclusive-source-forward", "loci": [
        {"id": source, "source_id": source, "start": 1, "end": len(coding),
         "genes": [{**gene, "id": f"g{i:03d}"} for i in range(count)]} for source in ("a", "b")]}
    return f">a\n{coding}\n>b\n{coding}\n".encode(), json.dumps(document).encode()


def test_duplicate_proteins_reuse_alignment_without_losing_gene_evidence(monkeypatch):
    original = Align.PairwiseAligner
    calls = []
    class CountingAligner:
        def __init__(self, **kwargs):
            object.__setattr__(self, "aligner", original(**kwargs))
        def __setattr__(self, name, value):
            setattr(self.aligner, name, value)
        def align(self, one, two):
            calls.append((one, two))
            return self.aligner.align(one, two)
    monkeypatch.setattr(Align, "PairwiseAligner", CountingAligner)
    report, _ = compare(*repeated_inputs())
    assert report["evaluated_protein_pairs"] == 9
    assert report["alignment_cells"] == 36
    assert len(report["links"]) == 9
    assert len(calls) == 1
    assert {(link["query_gene"], link["target_gene"]) for link in report["links"]} == set(product(
        ("g000", "g001", "g002"), repeat=2))
    assert all(link["identity"] == 1 and link["score"] == 9 for link in report["links"])


def test_ambiguous_codons_keep_biopython_resolution_and_explicit_uncertainty():
    report, _ = compare(*repeated_inputs("ATGGCNRAYSARWTANNNTANTAA", "MABZJXX", 1))
    assert report["loci"][0]["genes"][0]["translation"] == "MABZJXX"
    assert report["links"][0]["aligned_query"] == "MABZJXX"
    assert report["links"][0]["identity"] == 1


def test_aggregate_work_accounting_matches_explicit_pairs_at_cell_boundary(monkeypatch):
    from clinker import bounded
    fasta, raw = repeated_inputs()
    document = json.loads(raw)
    records = {name: "ATGGCTTAA" for name in ("a", "b")}
    lengths = [len(gene["translation"]) for gene in document["loci"][0]["genes"]]
    expected_cells = sum(a * b for a, b in product(lengths, lengths))
    monkeypatch.setitem(bounded.LIMITS, "alignment_cells", expected_cells)
    _, pairs, cells = bounded.validate_features(document, records)
    assert (pairs, cells) == (9, expected_cells)
    monkeypatch.setitem(bounded.LIMITS, "alignment_cells", expected_cells - 1)
    with pytest.raises(ValueError, match="work limit"):
        bounded.validate_features(document, records)


def test_total_cross_locus_pair_limit_is_enforced_before_alignment(monkeypatch):
    from clinker import bounded
    fasta, raw = repeated_inputs(count=15)
    document = json.loads(raw)
    gene_locus = document["loci"][0]
    document["loci"] = [{**gene_locus, "id": f"locus{i}"} for i in range(8)]
    def unexpected_alignment(**kwargs):
        pytest.fail("Work limit must be checked before constructing the aligner")
    monkeypatch.setattr(bounded.Align, "PairwiseAligner", unexpected_alignment)
    with pytest.raises(ValueError, match="work limit"):
        bounded.compare(fasta, json.dumps(document).encode())


def test_complete_session_deserialization_can_be_repeated():
    one = Cluster("one", [Locus("one", [named_gene("a")])])
    two = Cluster("two", [Locus("two", [named_gene("b")])])
    alignment = Alignment(query=one, target=two)
    alignment.add_link(one.loci[0].genes[0], two.loci[0].genes[0], .5, .75)
    aligner = Globaligner()
    aligner.add_alignment(alignment)
    raw = aligner.to_dict()
    frozen = copy.deepcopy(raw)
    for _ in range(2):
        restored = Globaligner.from_dict(raw)
        assert restored.to_dict() == frozen
        assert raw == frozen


@pytest.mark.parametrize("delimiter", (None, ",", "\t"))
def test_native_empty_links_format_as_valid_negative_evidence(delimiter):
    from clinker.formatters import format_links
    assert format_links([], headers=False, delimiter=delimiter) == ""
    assert "Identity" in format_links([], headers=True, delimiter=delimiter)
