def test_runtime_fast_imports():
    from casandra.runtime_fast import scan_profiles_fast
    assert callable(scan_profiles_fast)

