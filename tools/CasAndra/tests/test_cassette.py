from casandra.cassette import classify_cassette, construct_cassettes


def prediction(sequence_id, is_cas, subtype=None, score=10):
    cas_type = subtype.split("-")[0] if subtype else None
    return {
        "sequence_id": sequence_id,
        "is_cas": is_cas,
        "positive_profile_score": score,
        "classification": {"class": 2 if cas_type == "V" else 1, "type": cas_type, "subtype": subtype},
    }


def gene(protein_id, order, start):
    return {
        "protein_id": protein_id,
        "contig_id": "ctg",
        "gene_order": order,
        "start_1based": start,
        "end_1based_inclusive": start + 99,
    }


def test_weighted_cassette_vote_and_abstention():
    rows = [prediction("a", True, "V-A", 20), prediction("b", True, "V-B", 10)]
    assert classify_cassette(rows, 0.6)["subtype"] == "V-A"
    assert classify_cassette(rows, 0.8)["subtype"] is None


def test_nearby_cas_hits_bridge_limited_noncas_genes():
    genes = [gene("a", 1, 1), gene("x", 2, 200), gene("b", 3, 400)]
    predictions = [
        prediction("a", True, "V-A"),
        prediction("x", False),
        prediction("b", True, "V-A"),
    ]
    cassettes = construct_cassettes(genes, predictions, max_non_cas_genes=1)
    assert len(cassettes) == 1
    assert cassettes[0]["cas_gene_count"] == 2
    assert cassettes[0]["context_gene_count"] == 3


def test_distant_hits_form_separate_cassettes():
    genes = [gene("a", 1, 1), gene("b", 2, 30_000)]
    predictions = [prediction("a", True, "V-A"), prediction("b", True, "V-A")]
    assert len(construct_cassettes(genes, predictions)) == 2
