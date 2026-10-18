from casandra.cassette_hybrid import classify_hybrid_cassette


class FixedArchitectureModel:
    def __init__(self, subtype, confidence):
        self.subtype = subtype
        self.confidence = confidence
        self.calls = 0

    def classify(self, *_args):
        self.calls += 1
        return self.subtype, self.confidence


def test_non_type_iii_keeps_direct_profile_result():
    cassette = {
        "classification": {
            "class": "2",
            "type": "V",
            "subtype": "V-A",
            "confidence": {"subtype_vote_fraction": 0.9},
        }
    }
    result = classify_hybrid_cassette(cassette, None, {}, {})
    assert result["subtype"] == "V-A"
    assert result["method"] == "direct_profile_aggregation"


def test_type_ii_uses_ordered_architecture_and_preserves_direct_profile_result():
    direct = {
        "class": "2",
        "type": "II",
        "subtype": "II-C",
        "confidence": {"subtype_vote_fraction": 1.0},
    }
    cassette = {"classification": direct}
    model = FixedArchitectureModel("II-A", 0.54)

    result = classify_hybrid_cassette(cassette, model, {}, {})

    assert model.calls == 1
    assert result == {
        "class": "2",
        "type": "II",
        "subtype": "II-A",
        "method": "type_ii_ordered_profile_architecture_extratrees",
        "confidence": 0.54,
        "direct_profile_result": direct,
    }


def test_type_iii_keeps_its_explicit_architecture_method():
    direct = {
        "class": "1",
        "type": "III",
        "subtype": None,
        "confidence": {"subtype_vote_fraction": 0.0},
    }
    model = FixedArchitectureModel("III-A", 0.91)

    result = classify_hybrid_cassette({"classification": direct}, model, {}, {})

    assert result["subtype"] == "III-A"
    assert result["method"] == "type_iii_ordered_profile_architecture_extratrees"
    assert result["direct_profile_result"] is direct
