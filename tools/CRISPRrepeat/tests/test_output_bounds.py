"""Serialization budgets apply to expanded JSON, not shared Python objects."""

import json

import pytest

from crisprrepeat import exports


def test_repeated_reference_metadata_aborts_before_full_serialization(tmp_path, monkeypatch):
    shared = {"accession": "example", "evidence": "x" * 80}
    report = {"instances": [], "references": [shared] * 10000}
    monkeypatch.setattr(exports, "MAX_OUTPUT_BYTES", 512)
    # A full JSON allocation would defeat the limit before it could be checked.
    monkeypatch.setattr(exports.json, "dumps", lambda *a, **k: pytest.fail("Unbounded JSON allocation"))
    output = tmp_path / "limited"
    with pytest.raises(ValueError, match="output byte limit"):
        exports.write_outputs(output, report)
    assert list(output.iterdir()) == []


def test_streamed_json_byte_limit_includes_newline(tmp_path, monkeypatch):
    report = {"instances": [], "provenance": "Δοκιμή"}
    expected = (json.dumps(report, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    monkeypatch.setattr(exports, "MAX_OUTPUT_BYTES", len(expected))
    output = tmp_path / "exact"
    exports.write_outputs(output, report)
    assert (output / "repeat-analysis.json").read_bytes() == expected
    assert not (output / ".report.tmp").exists()
    monkeypatch.setattr(exports, "MAX_OUTPUT_BYTES", len(expected) - 1)
    with pytest.raises(ValueError, match="output byte limit"):
        exports.write_outputs(tmp_path / "over", report)


def test_nonfinite_output_leaves_no_completed_or_private_report(tmp_path):
    output = tmp_path / "nonfinite"
    with pytest.raises(ValueError):
        exports.write_outputs(output, {"instances": [], "bad": float("nan")})
    assert list(output.iterdir()) == []
