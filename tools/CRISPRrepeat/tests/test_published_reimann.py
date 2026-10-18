"""Pinned-engine regression on attributed transcripts, not a cleavage classifier."""

import json
from pathlib import Path
import runpy
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "reimann_2017.py"
BENCHMARK = runpy.run_path(str(SCRIPT))


@pytest.fixture(scope="module")
def result():
    return BENCHMARK["run_benchmark"]()


def runs_for(result, case_id):
    case = next(case for case in result["cases"] if case["id"] == case_id)
    return {run["preset"]: run for run in case["runs"]}


def supports(run):
    return {row["id"]: row["computed_comparison"]["context_pair_support"] for row in run["repeats"]}


def test_published_sequence_binding_rejects_changed_transcript(tmp_path):
    fixture, _ = BENCHMARK["read_fixture"](BENCHMARK["DEFAULT_FIXTURE"])
    assert [case["length"] for case in fixture["cases"]] == [190, 121]
    fixture["cases"][0]["rna_sequence"] = "A" + fixture["cases"][0]["rna_sequence"][1:]
    path = tmp_path / "altered.json"
    path.write_text(json.dumps(fixture))
    with pytest.raises(ValueError, match="checksum mismatch"):
        BENCHMARK["read_fixture"](path)


def test_same_repeat_exhibits_published_context_contrast(result):
    runs = runs_for(result, "CRISPR2_VIII_190nt")
    for run in runs.values():
        assert run["model"]["version"] == "2.7.0"
        r6, r7 = run["repeats"]
        a, b = r6["computed_comparison"], r7["computed_comparison"]
        assert a["status"] == b["status"] == "completed"
        assert a["motif_origin"] == b["motif_origin"] == "isolated_mfe_hypothesis"
        assert a["motif_pairs"] == b["motif_pairs"]
        assert a["isolated_pair_support"] == b["isolated_pair_support"]
        assert a["context_pair_support"] < b["context_pair_support"]
        assert r6["experimental_cleavage"] == {"Cas6-1": False, "Cas6-2a": False}
        assert r7["experimental_cleavage"] == {"Cas6-1": True, "Cas6-2a": True}


def test_window_sensitivity_with_lonely_pair_policy_held_constant(result):
    runs = runs_for(result, "CRISPR2_VIII_190nt")
    local, whole = runs["local_no_lonely_pairs"], runs["whole_fragment_no_lonely_pairs"]
    assert local["model"]["no_lonely_pairs"] is whole["model"]["no_lonely_pairs"] is True
    assert local["effective_context"]["window"] == 150
    assert whole["effective_context"]["window"] == 190
    assert local["effective_context"]["span"] == whole["effective_context"]["span"] == 80
    # Numeric regression anchors for this engine/fixture, never biological cutoffs.
    assert supports(local) == pytest.approx({"R6": 0.5600637441352296, "R7": 0.9519081398447753}, abs=1e-6)
    assert supports(whole) == pytest.approx({"R6": 0.0049783170709506675, "R7": 0.9447940141408562}, abs=1e-6)
    assert supports(runs["default_local"]) == pytest.approx(
        {"R6": 0.5600921628370721, "R7": 0.9478926029214516}, abs=1e-6)


def test_experimental_counterexample_is_not_hidden_or_classified(result):
    challenge = runs_for(result, "CRISPR2_I_121nt")
    contrast = runs_for(result, "CRISPR2_VIII_190nt")
    for preset, run in challenge.items():
        row = run["repeats"][0]
        comparison = row["computed_comparison"]
        assert row["id"] == "R3"
        assert (row["start"], row["end"]) == (42, 79)
        assert row["experimental_cleavage"] == {"Cas6-1": True, "Cas6-2a": True}
        assert row["published_model_score"] == 0.0001
        assert comparison["context_pair_support"] < supports(contrast[preset])["R6"]
        assert "cleavage_prediction" not in comparison
        assert run["effective_context"]["window"] == 121
    assert challenge["default_local"]["effective_context"]["requested_window"] == 150


def test_benchmark_cli_records_provenance_and_finite_json(tmp_path):
    output = tmp_path / "benchmark.json"
    subprocess.run([sys.executable, str(SCRIPT), "--output", str(output)], check=True, cwd=tmp_path)
    report = json.loads(output.read_text())
    assert len(report["fixture_sha256"]) == len(report["script_sha256"]) == 64
    assert report["source"]["doi"] == "10.1093/nar/gkw786"
    assert report["elapsed_seconds"] >= 0
    assert len(report["limitations"]) >= 5
    assert len(report["cases"]) == 2
