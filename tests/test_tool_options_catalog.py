"""Keep independently validated API/worker options and the public menu aligned."""

import json
from pathlib import Path

from crispr_integration.scientific_options import OPTIONS as INTEGRATION_OPTIONS
from crispr_integration.scientific_options import TOOL_MODES as INTEGRATION_MODES
from crispr_web.scientific_options import OPTIONS, TOOL_MODES

ROOT = Path(__file__).resolve().parents[1]


def test_public_catalog_matches_both_independent_allowlists():
    catalog = json.loads((ROOT / "frontend/src/contracts/tool-options-v1.json").read_text())
    tools = {tool["id"]: tool for tool in catalog["tools"] if tool["flags"]}
    assert set(tools) == set(OPTIONS) == set(INTEGRATION_OPTIONS)
    for key, tool in tools.items():
        assert set(tool["modes"]) == set(TOOL_MODES[key]) == set(INTEGRATION_MODES[key])
        assert {field["flag"] for field in tool["flags"]} == set(OPTIONS[key]) == set(INTEGRATION_OPTIONS[key])
        for field in tool["flags"]:
            api = OPTIONS[key][field["flag"]]
            worker = INTEGRATION_OPTIONS[key][field["flag"]]
            assert vars(api) == vars(worker)
            assert field["default"] == api.default
            assert field["type"] == ("enum" if api.kind == "select" else api.kind)
            assert field.get("min") == api.minimum
            assert field.get("max") == api.maximum
            assert tuple(field.get("choices", ())) == api.choices
            assert tuple(field.get("modes", ())) == api.modes
