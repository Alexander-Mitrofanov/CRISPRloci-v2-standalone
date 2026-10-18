from casandra.cassette_architecture import cassette_features


def test_architecture_features_preserve_order_and_scores():
    cassette = {"cas_protein_ids": ["b", "a"]}
    genes = {"a": {"start_1based": 1}, "b": {"start_1based": 100}}
    predictions = {
        "a": {"best_positive_profile": "P1", "positive_profile_score": 100, "classification": {"type": "V", "subtype": "V-A"}},
        "b": {"best_positive_profile": "P2", "positive_profile_score": 50, "classification": {"type": "V", "subtype": "V-A"}},
    }
    features = cassette_features(cassette, genes, predictions)
    assert features["bigram:P1>P2"] == 1
    assert features["subtype:V-A:score"] == 1.5
