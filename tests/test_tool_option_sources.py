"""Check menu claims against locked upstream parsers without importing tool runtimes."""

import ast
import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
CATALOG = json.loads((ROOT / "frontend/src/contracts/tool-options-v1.json").read_text())
LOCK = json.loads((ROOT / "tools/tools.lock.json").read_text())
SOURCES = {
    "clinker": ("clinker-v2", "clinker/bounded.py"),
    "phigaro": ("Phigaro-v2", "phigaro/bounded.py"),
    "isescan": ("ISEScan-v2", "isescan_v2/__main__.py"),
    "defensefinder": ("DefenseFinder-v2", "defense_context/pipeline.py"),
    "rnaalifold": ("RNAalifold-v2", "rna_family_tool/cli.py"),
    "cctk": ("CCTK", "src/cctk/cli.py"),
    "spacepharer": ("SpacePHARER", "python/spacepharer/cli.py"),
    "crispridentify": ("CRISPRidentify-v2", "CRISPRidentify.py"),
    "casandra": ("CasAndra", "src/casandra/cli.py"),
    "crispr_tracrrna": ("CRISPRtracrRNA-v3", "src/crisprtracrrna/cli.py"),
    "spacerplacer": ("SpacerPlacer", "spacerplacer.py"),
    "crispr_evor": ("evOr", "spacerplacer.py"),
    "crisprleader": ("CRISPRleader-v2", "src/crisprleader/cli.py"),
    "crisprrepeat": ("CRISPRrepeat", "src/crisprrepeat/cli.py"),
    "repeattyper": ("RepeatTyper", "src/repeattyper/cli.py"),
    "crisprmap": ("CRISPRmap-v2", "src/crisprmap/cli.py"),
    "crisprspacer": ("CRISPRspacer", "src/crisprspacer/cli.py"),
}


def locked_tree(tool_name, relative_path):
    tool = next(row for row in LOCK["tools"] if row["name"] == tool_name)
    checkout = ROOT / tool["path"]
    if (ROOT / "tools/bundled-sources.json").is_file():
        sys.path.insert(0, str(ROOT / "scripts"))
        from bundled_tools import verify_bundled

        assert verify_bundled(ROOT, tool)
        return ast.parse((checkout / relative_path).read_text()), tool
    if not (checkout / ".git").exists():
        pytest.skip(f"Independent {tool_name} checkout is not provisioned")
    # A checkout can be on a newer development revision. Its locked Git object,
    # not the working tree or HEAD, defines the deployed flag contract.
    source = subprocess.run(
        ["git", "-C", str(checkout), "show", f"{tool['commit']}:{relative_path}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return ast.parse(source), tool


def parser_options(tree):
    result = {}
    for call in ast.walk(tree):
        if not (
            isinstance(call, ast.Call)
            and isinstance(call.func, ast.Attribute)
            and call.func.attr == "add_argument"
        ):
            continue
        keywords = {item.arg: item.value for item in call.keywords}
        for argument in call.args:
            if isinstance(argument, ast.Constant) and isinstance(argument.value, str):
                if argument.value.startswith("--"):
                    result.setdefault(argument.value, []).append(keywords)
    # Resolve the CLI's simple literal loop-generated long flags without imports.
    for loop in (node for node in ast.walk(tree) if isinstance(node, ast.For)):
        if not isinstance(loop.target, ast.Name) or not isinstance(loop.iter, (ast.Tuple, ast.List)):
            continue
        for call in (node for node in ast.walk(loop) if isinstance(node, ast.Call)):
            if not isinstance(call.func, ast.Attribute) or call.func.attr != "add_argument":
                continue
            for argument in call.args:
                if (isinstance(argument, ast.BinOp) and isinstance(argument.op, ast.Add)
                    and isinstance(argument.left, ast.Constant) and argument.left.value == "--"
                    and isinstance(argument.right, ast.Name) and argument.right.id == loop.target.id):
                    for item in loop.iter.elts:
                        if isinstance(item, ast.Constant) and isinstance(item.value, str):
                            result.setdefault("--" + item.value, []).append({field.arg: field.value for field in call.keywords})
    return result


def settings_defaults(tree):
    settings = next(node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == "Settings")
    return {
        field.target.id: ast.literal_eval(field.value)
        for field in settings.body
        if isinstance(field, ast.AnnAssign)
        and isinstance(field.target, ast.Name)
        and isinstance(field.value, (ast.Constant, ast.UnaryOp))
    }


def argument_default(arguments, flag, defaults):
    if "default" in arguments:
        value = ast.literal_eval(arguments["default"])
        if value is not None:
            return value
    action = arguments.get("action")
    if isinstance(action, ast.Constant):
        if action.value == "store_true":
            return False
        if action.value == "store_false":
            return True
    # tracrRNA intentionally leaves parser defaults unset, then merges its
    # versioned Settings dataclass; booleans can explicitly default to None.
    return defaults[flag.removeprefix("--").replace("-", "_")]


def test_catalog_covers_each_locked_tool_once():
    ids = [tool["id"] for tool in CATALOG["tools"]]
    assert len(ids) == len(set(ids))
    assert set(ids) == set(SOURCES)
    assert {source[0] for source in SOURCES.values()} == {tool["name"] for tool in LOCK["tools"]}


@pytest.mark.parametrize("tool", CATALOG["tools"], ids=lambda tool: tool["id"])
def test_editable_and_managed_flags_are_shipped_and_defaults_are_real(tool):
    name, source_path = SOURCES[tool["id"]]
    tree, locked = locked_tree(name, source_path)
    arguments = parser_options(tree)
    assert tool["revision"] == locked["commit"]
    assert locked["commit"] in tool["documentation"]
    defaults = {}
    if tool["id"] == "crispr_tracrrna":
        config_tree, _ = locked_tree(name, "src/crisprtracrrna/config.py")
        defaults = settings_defaults(config_tree)
    for entry in tool["flags"]:
        flag = entry["flag"]
        assert flag in arguments, f"{name} does not ship {flag}"
        definitions = arguments[flag]
        assert any(
            argument_default(definition, flag, defaults) == entry["default"]
            for definition in definitions
        ), f"{name} {flag} menu default differs from the locked CLI"
        if "choices" in entry:
            assert any(
                "choices" in definition
                and list(ast.literal_eval(definition["choices"])) == entry["choices"]
                for definition in definitions
            ), f"{name} {flag} menu choices differ from the locked CLI"
    for group in tool["managed"]:
        for flag in group["flags"]:
            assert flag in arguments, f"{name} does not ship managed flag {flag}"
