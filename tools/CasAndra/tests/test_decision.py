from casandra.decision import is_cas_hit


def test_no_hit_zero_margin_is_not_cas():
    assert not is_cas_hit({"best_positive_profile": None, "margin": 0.0}, 0.0)


def test_actual_positive_hit_at_threshold_is_cas():
    assert is_cas_hit({"best_positive_profile": "Cas12a", "margin": 0.0}, 0.0)
