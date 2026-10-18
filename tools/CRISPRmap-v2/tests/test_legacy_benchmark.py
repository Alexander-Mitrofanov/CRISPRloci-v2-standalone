"""Test scientific split units and abstention accounting with independent fixtures."""
import importlib.util
from pathlib import Path

from crisprmap.reference import build_reference

spec = importlib.util.spec_from_file_location("run_benchmark", Path(__file__).parents[1] / "benchmarks/run_benchmark.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def reference_fixture():
    records = [
        {"id": "a", "sequence": "AACGT"}, {"id": "a_rc", "sequence": "ACGTT"},
        {"id": "b", "sequence": "AATGT"}, {"id": "b_missing", "sequence": "AATGT"},
        {"id": "c", "sequence": "CCCCG"}, {"id": "d", "sequence": "ATATA"},
    ]
    labels = {"a": "F1", "a_rc": "F1", "b": "F2", "b_missing": None, "c": "F3", "d": None}
    return build_reference(records, {ident: {"family": label} for ident, label in labels.items()},
                           {"namespace": "synthetic-test", "source_url": "local:synthetic", "source_sha256": "0" * 64})


def test_holdout_removes_every_occurrence_and_reverse_complement():
    reference = reference_fixture()
    for requested in range(1, 4):
        queries, truths, training, group_count, selected = benchmark.make_holdout(reference, requested)
        assert group_count == 4
        assert len(queries) == len(truths) == len(selected) == requested
        assert not selected.intersection(benchmark.group_key(row["sequence"]) for row in training["entries"])
        assert sum(len(group) for group in truths) + len(training["entries"]) == len(reference["entries"])


def test_family_absent_excludes_conflicting_occurrences_and_never_assigns():
    result = benchmark.family_absent(reference_fixture(), max_distance=3)
    assert result["family_count"] == 2  # F2 has only a partially labelled query group, excluded from scored truth.
    assert result["query_count"] == 2
    assert result["query_assignment_count"] == 0
    assert result["per_family"]["F1"]["training_records"] == 4  # a and its reverse complement both removed.


def test_missing_original_labels_are_not_accuracy_negatives():
    rows = [{"annotations": {"family": {"value": None}, "motif": {"value": None}}}]
    result = benchmark.label_comparison(rows, [{"family": None, "motif": None}])
    assert result["family"]["known_original_labels"] == 0
    assert result["family"]["agreement_among_assigned_known_labels"] is None
    assert result["family"]["raw_value_matches_including_null"] == 1
    assert benchmark.formatted_rate(None) == "unavailable (zero eligible cases)"


def test_benchmark_rejects_valid_but_relabelled_reference(tmp_path, monkeypatch, capsys):
    import hashlib
    import json
    import pytest
    raw = b">Taxon_A_G_1_M1_F2\tNC_123456_1_5_9\nACGT\n"
    (tmp_path / "reference.fa").write_bytes(raw)
    monkeypatch.setitem(benchmark.legacy.ASSETS, "reference.fa", {
        "url": "https://example.invalid/pinned.fa", "bytes": len(raw),
        "sha256": hashlib.sha256(raw).hexdigest(),
    })
    expected = benchmark.legacy.build_legacy_reference(tmp_path)
    relabelled = build_reference(expected["entries"], {"Taxon_A_G_1_M1_F2": {"family": "F99", "motif": "M1", "accession": "NC_123456"}}, expected["provenance"])
    (tmp_path / "reference.json").write_text(json.dumps(relabelled))
    with pytest.raises(SystemExit) as error:
        benchmark.main(["--data-dir", str(tmp_path), "--output-dir", str(tmp_path / "results")])
    assert error.value.code == 2
    assert "rebuilt from pinned original bytes and labels" in capsys.readouterr().err
    assert not (tmp_path / "results").exists()
