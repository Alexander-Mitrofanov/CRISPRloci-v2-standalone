from casandra.cassette_filter import passes_genome_evidence_gate


def test_class1_singleton_is_rejected_but_two_gene_subtyped_cassette_passes():
    base = {"final_classification": {"class": "1", "subtype": "I-C"}, "cas_protein_ids": ["a"]}
    assert not passes_genome_evidence_gate({**base, "cas_gene_count": 1}, {}, {})[0]
    assert passes_genome_evidence_gate({**base, "cas_gene_count": 2}, {}, {})[0]


def test_strong_long_class2_single_effector_passes():
    cassette = {
        "final_classification": {"class": "2", "subtype": "V-A"},
        "cas_gene_count": 1,
        "cas_protein_ids": ["p"],
    }
    genes = {"p": {"protein_length": 1000}}
    predictions = {"p": {"positive_profile_score": 100}}
    assert passes_genome_evidence_gate(cassette, genes, predictions)[0]


def test_bundle_policy_overrides_single_effector_thresholds():
    cassette = {
        "final_classification": {"class": "2", "subtype": "V-A"},
        "cas_gene_count": 1,
        "cas_protein_ids": ["p"],
    }
    genes = {"p": {"protein_length": 500}}
    predictions = {"p": {"positive_profile_score": 60}}
    strict = {
        "class_2_single_effector_minimum_aa": 600,
        "class_2_single_effector_minimum_positive_profile_score": 70,
    }
    assert not passes_genome_evidence_gate(cassette, genes, predictions, strict)[0]
