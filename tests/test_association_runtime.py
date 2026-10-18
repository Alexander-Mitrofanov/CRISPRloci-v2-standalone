"""Production SpacePHARER build recipe and source/runtime receipt boundaries."""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import provision_tools
import manage_release
from runtime.common import source_inventory, environment_inventory, load_lock


def test_vendored_source_recipe_receipt_and_old_identity_stability(tmp_path):
    for name in ("SpacePHARER", "RepeatTyper"):
        root = tmp_path / name
        recipe = root / "lib/mmseqs/lib/zstd/build/cmake/lib/CMakeLists.txt"
        recipe.parent.mkdir(parents=True)
        recipe.write_text("source recipe")
        artifact = root / "build/generated"
        artifact.parent.mkdir()
        artifact.write_text("generated")
        actual = source_inventory(root)
        assert (str(recipe.relative_to(root)) in actual) == (name == "SpacePHARER")
        assert "build/generated" not in actual
        if name == "SpacePHARER":
            previous = actual.copy()
            recipe.write_text("changed")
            assert source_inventory(root) != previous


def test_native_recipe_compiles_tests_installs_and_receipts_executable(
    tmp_path, monkeypatch
):
    prefix = tmp_path / "env"
    (prefix / "bin").mkdir(parents=True)
    build = tmp_path / "build"
    (build / "src").mkdir(parents=True)
    (build / "src/spacepharer-native").write_text("native executable")
    (build / "src/spacepharer-native").chmod(0o755)
    commands = []
    monkeypatch.setattr(
        provision_tools,
        "run",
        lambda args, **kwargs: commands.append([str(v) for v in args]),
    )
    monkeypatch.setattr(
        provision_tools.subprocess, "check_output", lambda *a, **kw: "2.1.0\n"
    )
    result = provision_tools.build_spacepharer_runtime(
        tmp_path / "source", prefix, build, "2.1.0", {}
    )
    assert "-DHAVE_SSE4_1=1" in commands[0] and "-DHAVE_AVX2=0" in commands[0]
    assert any("--output-on-failure" in args for args in commands)
    native = result["spacepharer-native"]
    assert native["version"] == "2.1.0"
    inventory = environment_inventory(prefix)
    assert inventory["bin/spacepharer-native"]["sha256"] == native["sha256"]
    assert inventory["bin/spacepharer-native"]["executable"]
    assert manage_release.scientific_manifest_version({"spacepharer": {}}) == "8.0.0"


def test_current_toolset_preserves_spacepharer_pin():
    path = Path(__file__).resolve().parents[1] / "tools/tools.lock.json"
    doc = load_lock(path)
    assert len(doc["tools"]) >= 12
    tool = next(t for t in doc["tools"] if t["name"] == "SpacePHARER")
    assert tool["version"] == "2.1.1" and len(tool["commit"]) == 40
    assert json.loads(path.read_text()) == doc
