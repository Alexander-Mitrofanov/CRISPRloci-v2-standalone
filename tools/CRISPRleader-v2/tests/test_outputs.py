"""Synthetic engineering regressions, independent of scientific prediction models."""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from urllib.parse import unquote

import pytest

import crisprleader.export as export_module
from crisprleader.errors import InputError
from crisprleader.export import write_bundle
from crisprleader.models import AnalysisConfig, SourceRecord
from crisprleader.pipeline import extract
from crisprleader.validation import read_json, validate_result

_DNA = "ACGTRYSWKMBDHVNACGTACGATTCGAGCTAACGTGCAAT"
_SOURCE_ID = "seq|one:synthetic%"
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _oriented(sequence: str, strand: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1] if strand == "-" else sequence


def _slice(segments: list[dict[str, int]]) -> str:
    return "".join(_DNA[part["start"]:part["end"]] for part in segments)


@pytest.fixture
def context(tmp_path: Path) -> tuple[dict, dict[str, SourceRecord], Path]:
    source = SourceRecord(_SOURCE_ID, _DNA, "circular", "explicit synthetic fixture topology")
    fasta = tmp_path / "source.fa"
    fasta.write_text(f">{_SOURCE_ID} synthetic fixture\n{_DNA}\n")
    arrays = []
    for identifier, segments, anchors in (
        ("near-start", [{"start": 4, "end": 18}],
         [[{"start": 4, "end": 8}], [{"start": 14, "end": 18}]]),
        ("wrapped", [{"start": len(_DNA) - 5, "end": len(_DNA)}, {"start": 0, "end": 9}],
         [[{"start": len(_DNA) - 5, "end": len(_DNA) - 1}], [{"start": 5, "end": 9}]]),
    ):
        arrays.append({
            "id": identifier, "source_id": source.id, "segments": segments,
            "terminal_repeats": [{"segments": anchor, "sequence": _slice(anchor)} for anchor in anchors],
            "consensus_repeat": _slice(anchors[0]), "strand": "?", "strand_constraint": False,
        })
    annotation = tmp_path / "arrays.json"
    annotation.write_text(json.dumps({
        "schema_version": "1.0.0", "coordinate_system": "zero_based_half_open",
        "sources": [{"id": source.id, "length": source.length, "sha256": source.sha256,
                     "topology": source.topology, "topology_provenance": source.topology_provenance}],
        "arrays": arrays,
    }))
    report = extract(fasta, annotation, config=AnalysisConfig(flank_length=9))
    return report, {source.id: source}, tmp_path


def test_context_extraction_is_explicit_and_source_bound(context: tuple) -> None:
    report, sources, _ = context
    assert report["analysis_mode"] == "context_extraction"
    assert report["execution_status"] == "completed"
    assert report["reference"] is None
    assert report["predictions"] == []
    assert len(report["candidates"]) == 4
    assert {candidate["strand"] for candidate in report["candidates"]} == {"+", "-"}
    assert any(len(candidate["segments"]) == 2 for candidate in report["candidates"])
    assert any(len(array["segments"]) == 2 for array in report["arrays"])
    for candidate in report["candidates"]:
        assert candidate["sequence"] == _oriented(_slice(candidate["segments"]), candidate["strand"])
    validate_result(report, sources)


def test_export_formats_reconstruct_circular_arcs_and_preserve_identity(context: tuple) -> None:
    report, _, temporary = context
    output = write_bundle(report, temporary / "export")
    fasta: dict[str, str] = {}
    identifier = ""
    for line in (output / "sequences.fasta").read_text().splitlines():
        if line.startswith(">"):
            identifier = line[1:].split()[0]
            assert identifier not in fasta
            fasta[identifier] = ""
            source_field = line.split("source=", 1)[1]
            assert unquote(source_field) == _SOURCE_ID
        else:
            fasta[identifier] += line
    gff_features: dict[str, list[tuple[list[str], dict[str, str]]]] = defaultdict(list)
    by_part = {}
    for line in (output / "features.gff3").read_text().splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        assert len(columns) == 9
        assert unquote(columns[0]) == _SOURCE_ID
        attributes = dict(part.split("=", 1) for part in columns[8].split(";"))
        gff_features[attributes["feature_id"]].append((columns, attributes))
        by_part[attributes["ID"]] = columns
    candidates = {(row["array_id"], row["strand"]): row for row in report["candidates"]}
    assert set(gff_features) == set(fasta)
    assert len(fasta) == len(candidates)
    saw_wrapped = False
    for feature_id, parts in gff_features.items():
        parts.sort(key=lambda item: int(item[1]["part"].split("/")[0]))
        strand = parts[0][0][6]
        array_id = unquote(parts[0][1]["array_id"])
        candidate = candidates[(array_id, strand)]
        gff_segments = [{"start": int(cols[3]) - 1, "end": int(cols[4])} for cols, _ in parts]
        assert gff_segments == candidate["segments"]
        assert fasta[feature_id] == _oriented(_slice(gff_segments), strand)
        assert fasta[feature_id] == candidate["sequence"]
        assert all(attrs["evidence_kind"] == "extracted_flank" for _, attrs in parts)
        assert all(cols[2] == "sequence_feature" for cols, _ in parts)
        if len(parts) == 2:
            saw_wrapped = True
            assert gff_segments[0]["end"] == len(_DNA)
            assert gff_segments[1]["start"] == 0
            assert {attrs["part"] for _, attrs in parts} == {"1/2", "2/2"}
    assert saw_wrapped
    bed_rows = [line.split("\t") for line in (output / "features.bed").read_text().splitlines()]
    assert len(bed_rows) == len(by_part)
    for chromosome, start, end, name, score, strand in bed_rows:
        assert chromosome == _SOURCE_ID  # BED has no GFF percent-decoding convention.
        gff = by_part[name]
        assert (int(start), int(end), strand) == (int(gff[3]) - 1, int(gff[4]), gff[6])
        assert score == "0"
    assert read_json(output / "result.json") == report
    manifest = json.loads((output / "manifest.json").read_text())
    for name, expected in manifest["files"].items():
        content = (output / name).read_bytes()
        assert expected == {"size_bytes": len(content), "sha256": hashlib.sha256(content).hexdigest()}


def test_preexisting_outputs_and_dangling_symlinks_are_never_replaced(context: tuple) -> None:
    report, _, temporary = context
    directory = temporary / "existing"
    directory.mkdir()
    sentinel = directory / "sentinel.txt"
    sentinel.write_text("keep this file")
    symlink = temporary / "dangling"
    symlink.symlink_to(temporary / "missing")
    for target in (directory, symlink):
        with pytest.raises(InputError, match="exists"):
            write_bundle(report, target)
    assert sentinel.read_text() == "keep this file"
    assert symlink.is_symlink()


def test_concurrent_empty_target_is_not_replaced(context: tuple, monkeypatch: pytest.MonkeyPatch) -> None:
    report, _, temporary = context
    target = temporary / "concurrent"
    publish = export_module._publish_no_replace
    concurrent_inode: list[int] = []

    def create_competing_target(source: Path, destination: Path) -> None:
        destination.mkdir()
        concurrent_inode.append(destination.stat().st_ino)
        publish(source, destination)

    monkeypatch.setattr(export_module, "_publish_no_replace", create_competing_target)
    with pytest.raises(FileExistsError):
        write_bundle(report, target)
    assert target.stat().st_ino == concurrent_inode[0]
    assert list(target.iterdir()) == []
    assert list(temporary.glob(".concurrent.*")) == []


def test_export_failure_removes_partial_staging_directory(context: tuple) -> None:
    report, _, temporary = context
    broken = deepcopy(report)
    broken["provenance"]["invalid_number"] = float("nan")
    with pytest.raises(ValueError):
        write_bundle(broken, temporary / "partial")
    assert not (temporary / "partial").exists()
    assert list(temporary.glob(".partial.*")) == []


def test_fifo_json_is_rejected_without_waiting_for_a_writer(tmp_path: Path) -> None:
    fifo = tmp_path / "input.json"
    os.mkfifo(fifo)
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(Path(export_module.__file__).resolve().parents[1])
    probe = subprocess.run(
        [sys.executable, "-c", "from crisprleader.errors import InputError; "
         "from crisprleader.validation import read_json; import sys\n"
         "try: read_json(sys.argv[1])\n"
         "except InputError as error: print(error); sys.exit(0)\n"
         "sys.exit(1)", str(fifo)],
        env=environment, capture_output=True, text=True, timeout=3, check=False,
    )
    assert probe.returncode == 0, probe.stderr
    assert "regular file" in probe.stdout


@pytest.mark.parametrize("mutate", [
    lambda doc: doc["sources"][0].update(id=[]),
    lambda doc: doc["sources"][0].update(sha256="0" * 64),
    lambda doc: doc["sources"][0].update(length=True),
    lambda doc: doc["sources"][0].update(topology_provenance=None),
    lambda doc: doc["arrays"][0].update(source_id="unknown"),
    lambda doc: doc["arrays"][0].update(id={}),
    lambda doc: doc["arrays"].append(deepcopy(doc["arrays"][0])),
    lambda doc: doc["arrays"][0].update(terminal_repeats=[None, None]),
    lambda doc: doc["arrays"][0].update(consensus_repeat="ACGU"),
    lambda doc: doc["candidates"][0].update(id=[]),
    lambda doc: doc["candidates"][0].update(array_id="unknown"),
    lambda doc: doc["candidates"][0].update(side=[]),
    lambda doc: doc["candidates"][0].update(strand="-"),
    lambda doc: doc["candidates"][0].update(requested_length=True),
    lambda doc: doc["candidates"][0].update(truncated_reason="contig_boundary"),
    lambda doc: doc["candidates"][0].update(observed_repeat="AAAA"),
    lambda doc: doc["candidates"][0].update(consensus_repeat="AAAA"),
    lambda doc: doc["candidates"].pop(),
    lambda doc: doc["candidates"].append(deepcopy(doc["candidates"][0])),
    lambda doc: doc.update(arrays=[False]),
    lambda doc: doc.update(reference={"asserted_bundle": "none"}),
    lambda doc: doc.update(execution_status="running"),
])
def test_malformed_identities_and_geometry_fail_as_input_errors(context: tuple, mutate) -> None:
    report, sources, _ = context
    broken = deepcopy(report)
    mutate(broken)
    with pytest.raises(InputError):
        validate_result(broken, sources)


def test_reconstructed_but_displaced_flank_is_rejected(context: tuple) -> None:
    report, sources, _ = context
    broken = deepcopy(report)
    candidate = broken["candidates"][0]
    candidate["segments"] = [{"start": 10, "end": 19}]
    candidate["sequence"] = _DNA[10:19]  # Slice agrees; adjacency intentionally does not.
    with pytest.raises(InputError, match="adjacent window"):
        validate_result(broken, sources)


def test_source_matching_but_interior_terminal_anchor_is_rejected(context: tuple) -> None:
    report, sources, _ = context
    broken = deepcopy(report)
    anchor = broken["arrays"][0]["terminal_repeats"][0]
    anchor["segments"] = [{"start": 5, "end": 9}]
    anchor["sequence"] = _DNA[5:9]
    with pytest.raises(InputError, match="flush"):
        validate_result(broken, sources)


@pytest.mark.parametrize("status", ["unsupported", "ambiguous", "no_supported_hit"])
@pytest.mark.parametrize("selection", [
    {"selected_candidate_id": "not-a-candidate"},
    {"core_start": 999, "core_end": 1000},
    {"core": {}},
])
def test_no_call_cannot_select_candidate_or_core(context: tuple, status: str, selection: dict) -> None:
    report, sources, _ = context
    document = deepcopy(report)
    document["analysis_mode"] = "prediction"
    document["predictions"] = [{
        "source_id": array["source_id"], "array_id": array["id"], "status": status,
        "selected_candidate_id": None, "core_start": None, "core_end": None, "core": None,
        "reasons": ["synthetic validator regression"], "evidence": [], "motifs": [], "details": {},
    } for array in document["arrays"]]
    document["predictions"][0].update(selection)
    with pytest.raises(InputError, match="cannot select"):
        validate_result(document, sources)


@pytest.mark.parametrize("contents", [
    '{"a":1,"a":2}', '{"a":NaN}', '{"a":1e999}', '[]', '{broken',
])
def test_result_reader_rejects_malformed_json(tmp_path: Path, contents: str) -> None:
    path = tmp_path / "invalid.json"
    path.write_text(contents)
    with pytest.raises(InputError):
        read_json(path)


@pytest.mark.parametrize("limit", [-1, 0, True, 1.5, None])
def test_result_reader_rejects_invalid_budget_before_opening(tmp_path: Path, limit) -> None:
    # A negative file read size would otherwise remove the promised read bound.
    with pytest.raises(InputError, match="byte budget"):
        read_json(tmp_path / "not-opened.json", limit=limit)


@pytest.mark.parametrize("config", [None, [], {}, {"flank_length": True},
                                    {"flank_length": 0}, {"flank_length": 5001}])
def test_result_requires_bounded_effective_window(context: tuple, config) -> None:
    report, sources, _ = context
    broken = deepcopy(report)
    broken["effective_config"] = config
    with pytest.raises(InputError, match="effective_config"):
        validate_result(broken, sources)


def test_effective_window_must_match_both_candidate_hypotheses(context: tuple) -> None:
    report, sources, _ = context
    broken = deepcopy(report)
    # The contexts still reconstruct perfectly, but the declared run configuration is false.
    broken["effective_config"]["flank_length"] = 8
    with pytest.raises(InputError, match="differs from effective_config"):
        validate_result(broken, sources)


@pytest.mark.parametrize("field,value", [
    ("source_sha256", "0" * 64), ("source_length", 1), ("source_length", True),
    ("topology", "unknown"), ("extraction_policy", "inferred-leader-boundary"),
    ("window_is_biological_boundary", True), ("window_is_biological_boundary", 0),
])
def test_context_provenance_cannot_disagree_with_geometry(context: tuple, field, value) -> None:
    report, sources, _ = context
    broken = deepcopy(report)
    broken["candidates"][0]["metadata"][field] = value
    with pytest.raises(InputError):
        validate_result(broken, sources)


def test_extra_candidate_metadata_survives_without_changing_bound_claims(context: tuple) -> None:
    report, sources, _ = context
    annotated = deepcopy(report)
    annotated["candidates"][0]["metadata"]["operator_note"] = {"confidence": None}
    validate_result(annotated, sources)
