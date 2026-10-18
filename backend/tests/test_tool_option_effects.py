"""Exercise real mapping CLI overrides and independent requested-setting gates."""

from dataclasses import replace
import json

import pytest

from crispr_integration.repeat_map import run_repeat_map_pipeline
from crispr_integration.runner import RunnerConfig
from crispr_web.summaries.repeat_map import public_report

from .test_repeat_map import NATIVE, map_pipeline, map_reference, validated
from .test_repeats import setup_output, validate

# Imported fixtures retain the real frozen-reference preparation used by map tests.
__all__ = ["map_pipeline", "map_reference"]


@pytest.mark.parametrize("max_distance,top_k", [(0, 1), (1, 3), (2, 100)])
def test_mapping_flags_change_real_hits_and_still_check_all_ties(map_pipeline, tmp_path, max_distance, top_k):
    claim, _, bindings, source, reference = map_pipeline
    options = {"crisprmap": {"--max-distance": max_distance, "--top-k": top_k}}
    output = tmp_path / "custom-map"
    run_repeat_map_pipeline(
        source, output,
        runner_config=RunnerConfig(("unused",), ("unused",), ("unused",), crisprmap=NATIVE, tool_options=options),
        reference_path=reference, reference_sha256=claim.repeat_reference_sha256, timeout_seconds=60,
    )
    custom = replace(claim, tool_options=options)
    assert validated((custom, output, bindings, source, reference)) == "completed"
    assert validated((claim, output, bindings, source, reference)) is None
    report = json.loads((output / "repeat-map/report.json").read_text())
    assert public_report(report)["policy"] == {
        "max_distance": max_distance, "top_k": top_k,
        "orientation_scope": "relative_to_reference_not_transcription",
    }
    rows = {row["id"]: row for row in report["results"]}
    assert rows["near"]["status"] == ("no_hit" if max_distance == 0 else "near")
    assert len(rows["exact"]["hits"]) == min(top_k, 25)
    # Display limits cannot erase conflicting evidence from undisplayed ties.
    assert rows["exact"]["best_hit_count"] == 25
    assert rows["exact"]["annotations"]["family"]["status"] == "conflicting"


def test_repeat_gate_binds_requested_model_instead_of_default(tmp_path):
    data = list(setup_output(tmp_path, "isolated"))
    data[4] = replace(data[4], tool_options={"crisprrepeat": {"--temperature": 20}})
    data[3]["tool_options"] = data[4].tool_options
    data[-1]()
    # An output that ignored the requested option must fail validation.
    assert validate(data) is None
    data[2]["model"]["temperature"] = 20
    data[-1]()
    assert validate(data) == "completed"


def test_repeat_gate_binds_context_window_and_span(tmp_path):
    data = list(setup_output(tmp_path, "context"))
    data[4] = replace(data[4], tool_options={"crisprrepeat": {"--window": 200, "--span": 100}})
    data[3]["tool_options"] = data[4].tool_options
    data[2]["model"].update(window=200, span=100)
    data[-1]()
    assert validate(data) is None
    for row in data[2]["contexts"]:
        row.update(requested_window=200, requested_span=100,
                   window=min(len(row["rna_sequence"]), 200), span=min(len(row["rna_sequence"]), 100))
    data[-1]()
    assert validate(data) == "completed"
