from casandra.api import Analyzer, DEFAULT_MODEL_BUNDLE


def test_analyzer_preloads_verified_default_bundle():
    analyzer = Analyzer.from_bundle(DEFAULT_MODEL_BUNDLE)
    summary = analyzer.inspect()
    assert summary["integrity"] == "verified"
    assert summary["profile_hmms"] > 0
    thread_settings = {
        key: value
        for key, value in analyzer.architecture_model.estimator.get_params().items()
        if key == "n_jobs" or key.endswith("__n_jobs")
    }
    assert thread_settings and set(thread_settings.values()) == {1}
