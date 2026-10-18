"""Control generators must preserve exact invariants and per-item stability."""

import importlib.util
import random
import sys
from collections import Counter
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
spec = importlib.util.spec_from_file_location("council_benchmark", SCRIPTS / "council_benchmark.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_exact_dinucleotide_and_ambiguity_invariants():
    for sequence in ("", "A", "AAAA", "ACGNTACRGT", "AACCGGTTACGTGACTAGCTAGCTGATCGATCGGAT"):
        for value in range(30):
            shuffled = module.dinucleotide_shuffle(sequence, random.Random(value))
            assert len(shuffled) == len(sequence)
            assert Counter(shuffled) == Counter(sequence)
            assert Counter(zip(shuffled, shuffled[1:], strict=False)) == Counter(
                zip(sequence, sequence[1:], strict=False)
            )
            assert all(shuffled[i] == base for i, base in enumerate(sequence) if base not in "ACGT")


def test_order_independent_seed_and_mono_invariants():
    sequences = {"a": "ACGT" * 15, "b": "AGCTTAGCRNNACGT" * 4}

    def generate(items):
        return {
            key: module.dinucleotide_shuffle(value, random.Random(module.seed(key, "test", 1)))
            for key, value in items
        }

    assert generate(sequences.items()) == generate(reversed(list(sequences.items())))
    assert generate(sequences.items()) == {
        key: value
        for key, value in generate([("new", "AGCT"), *sequences.items()]).items()
        if key != "new"
    }
    assert module.seed("a", "test", 1) != module.seed("a", "test", 2)
    source = "ACGTNRNACGT" * 3
    result = module.mono_shuffle(source, random.Random(10))
    assert Counter(result) == Counter(source)
    assert all(result[i] == base for i, base in enumerate(source) if base not in "ACGT")


def test_derangement_is_bijective_cross_locus_and_excludes_identical_and_rc():
    repeats = ["AACC", "ACGC", "ATAC", "CGGA", "GGTT", "AGCG"]
    pool = [
        {"identity": str(i), "label": str(i), "repeat": repeat} for i, repeat in enumerate(repeats)
    ]
    assignments = module.donor_assignment(pool)
    assert assignments == module.donor_assignment(pool[::-1])
    assert set(assignments) == set(assignments.values()) == {str(i) for i in range(len(pool))}
    for recipient, donor in assignments.items():
        assert recipient != donor
        assert repeats[int(donor)] not in {
            repeats[int(recipient)],
            module.reverse_complement(repeats[int(recipient)]),
        }


def test_default_projection_never_hides_experimental_candidates():
    from council_audit import projection

    original = {"id": "standard", "admission_branch": "standard", "export_id": "p1"}
    unexpected = {"id": "short", "admission_branch": "short_rescue", "export_id": "p2"}
    assert projection([original]) != projection([original, unexpected])
    assert projection([original], experimental=True) == projection(
        [original, unexpected], experimental=True
    )


def test_auditor_accepts_valid_zero_prediction_exports(tmp_path):
    from council_audit import audit_run

    from crisprtracrrna.adapters import sha256_file
    from crisprtracrrna.models import Sequence
    from crisprtracrrna.pipeline import atomic_json, write_exports

    source = tmp_path / "source.fasta"
    source.write_text(">zero\nAAAA\n")
    output = tmp_path / "output"
    output.mkdir()
    record = {
        "key": "f000001.s000001",
        "source": {
            "file": str(source),
            "file_sha256": sha256_file(source),
            "id": "zero",
            "length": 4,
            "sha256": Sequence("zero", "zero", "AAAA").sha256,
        },
        "predictions": [],
    }
    write_exports(output, [record])
    atomic_json(
        output / "manifest.json", {"files": {p.name: sha256_file(p) for p in output.iterdir()}}
    )
    atomic_json(
        output / "report.json",
        {"status": "complete", "records": [record], "summary": {"prediction_count": 0}},
    )
    assert audit_run(output) == {"predictions": 0, "proposals": 0, "evidence_files": 3}


def test_amended_settings_honor_explicit_config_and_reject_stock(tmp_path):
    import pytest

    selected_config = tmp_path / "selected.toml"
    selected_binary = module.ROOT / "native/fasta36-seeded-20260927/fasta36-council-seed20260927"
    selected_config.write_text(f'fasta36 = "{selected_binary}"\n')
    selected = module.execution_settings(selected_config, True, pinned=True)
    assert selected.threads == 1 and selected.short_anti_repeat_rescue
    assert Path(selected.fasta36) == selected_binary
    stock_config = tmp_path / "stock.toml"
    stock_config.write_text('fasta36 = "fasta36"\n')
    with pytest.raises(ValueError):
        module.execution_settings(stock_config, True, pinned=True)


def test_conditional_passes_requested_configuration(tmp_path, monkeypatch):
    import pytest

    (tmp_path / "controls").mkdir()
    (tmp_path / "controls/manifest.json").write_text('{"repeat_controls": []}')
    monkeypatch.setattr(module, "ROUND", tmp_path)
    selected = tmp_path / "selected.toml"

    class StopBeforeNative(Exception):
        pass

    def capture(config, experimental, pinned=False):
        assert config == selected and experimental and pinned
        raise StopBeforeNative

    monkeypatch.setattr(module, "execution_settings", capture)
    with pytest.raises(StopBeforeNative):
        module.conditional("test", selected)


def test_pinned_auditor_rejects_stock_metadata_and_wrong_stage_profile():
    import pytest
    from council_audit import require_stage_profile

    for name in ("current_default_pinned_native", "current_experimental_pinned_native"):
        with pytest.raises(AssertionError):
            require_stage_profile(name, {"profile": "default"})
    with pytest.raises(AssertionError):
        require_stage_profile(
            "current_experimental_pinned_native",
            {"profile": "default", "native_dependency_profile": "council_seed20260927_threads1"},
        )
    assert (
        require_stage_profile(
            "current_experimental_pinned_native",
            {
                "profile": "experimental",
                "native_dependency_profile": "council_seed20260927_threads1",
            },
        )
        == "experimental"
    )


def test_historical_retention_is_independent_of_new_baseline():
    from council_results import retention

    historical = {"old": {"detected": True}, "new": {"detected": False}}
    pinned = {"old": {"detected": False}, "new": {"detected": True}}
    result = retention(historical, {"archived": pinned, "current": pinned, "experimental": pinned})
    assert all(value["lost"] == ["old"] and not value["retained"] for value in result.values())
