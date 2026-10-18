"""Portable sequence contracts, source binding, and bounded job admission."""

import hashlib
import json
from pathlib import Path

import pytest

from crispr_web.exchange import INPUT_KINDS, SequenceSet, sequence_input, sequence_set
from crispr_web.fasta import FastaValidationError, normalize_fasta
from crispr_web.models import JobSubmission
from crispr_web.publication.array_data import array_documents
from crispr_web.publication.sequence_ports import write_document
from crispr_web.service import JobService
from crispr_web.worker import Worker

FIXTURE = Path(__file__).parent / "fixtures/repeats"


def detection_ports(tmp_path, mutate=None):
    source = normalize_fasta(
        (FIXTURE / "source.fasta").read_text(),
        max_records=10,
        max_total_bases=10000,
        max_record_bases=10000,
        max_header_characters=200,
    )
    bindings = {
        row.identifier: (
            row.sequence,
            hashlib.sha256(row.sequence.encode()).hexdigest(),
        )
        for row in source.records
    }
    doc = json.loads((FIXTURE / "identify.json").read_text())
    if mutate:
        mutate(doc)
    root = tmp_path / "results"
    path = root / "crispridentify/source/report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc))
    return array_documents(
        root,
        bindings,
        categories=["Bona-fide", "Possible"],
        mode="detection",
        max_files=100,
        max_bytes=1_000_000,
    )


def test_real_detection_ports_are_compact_complete_and_reusable(tmp_path):
    exports = detection_ports(tmp_path)
    native = json.loads((FIXTURE / "identify.json").read_text())
    compact = exports["arrays.json"]
    assert len(json.dumps(compact)) < len(json.dumps(native)) / 5
    assert compact["arrays"][0]["boundary"]["exact"] is False
    assert "feature_sets" not in json.dumps(compact)
    for kind, modes, count in [
        ("repeat", ["repeats", "repeat_map", "repeat_type"], 3),
        ("spacer", ["viral_search", "spacer_association"], 2),
    ]:
        packet = exports[kind + "s.json"]
        assert len(packet["records"]) == count
        for mode in modes:
            text, metadata = sequence_input(
                json.dumps(packet), mode=mode, molecule="DNA", max_bytes=100000
            )
            assert text.count(">") == count
            assert len(metadata["origins"]) == count
        for row, raw in zip(
            sorted(packet["records"], key=lambda row: row["origin"]["ordinal"]),
            native["arrays"][0][kind + "s"],
        ):
            assert row["sequence"] == raw["sequence"]
            assert row["origin"]["start"] == raw["source_interval"]["start"]


@pytest.mark.parametrize("tamper", [None, "native_id", "digest", "duplicate_source", "directory"])
def test_locus_ports_restore_original_ids_from_bound_record_keys(tmp_path, tamper):
    sequence = (FIXTURE / "source.fasta").read_text().splitlines()
    sequence = "".join(line for line in sequence if not line.startswith(">"))
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    bindings = {source_id: (sequence, digest) for source_id in ("first", "second")}
    records = []
    for index, source_id in enumerate(bindings, 1):
        key = f"f000001.s{index:06d}"
        native = json.loads((FIXTURE / "identify.json").read_text())
        native["source"]["id"] = "wrong" if tamper == "native_id" else key
        directory = "f000002.s000001" if tamper == "directory" else key
        path = tmp_path / "tracrrna/work" / directory / "identify/source/report.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(native))
        records.append({
            "key": key,
            "source": {"id": "first" if tamper == "duplicate_source" else source_id,
                       "length": len(sequence), "sha256": "0" * 64 if tamper == "digest" else digest},
        })
    (tmp_path / "tracrrna/report.json").write_text(json.dumps({"status": "complete", "records": records}))
    options = dict(categories=["Bona-fide", "Possible"], mode="loci", max_files=100, max_bytes=1_000_000)
    if tamper:
        with pytest.raises(ValueError, match="(locus|array source)"):
            array_documents(tmp_path, bindings, **options)
        return
    documents = array_documents(tmp_path, bindings, **options)
    assert {source["id"] for source in documents["arrays.json"]["sources"]} == set(bindings)
    for kind in ("repeats", "spacers"):
        packet = documents[kind + ".json"]
        assert packet["complete"] is True
        assert {row["origin"]["source_id"] for row in packet["records"]} == set(bindings)
        assert len({row["id"] for row in packet["records"]}) == len(packet["records"])


def test_unknown_orientation_and_missing_sequences_are_not_invented(tmp_path):
    def mutate(doc):
        doc["arrays"][0]["orientation"] = {"strand": None}
        doc["arrays"][0]["repeats"][0].pop("sequence")

    exports = detection_ports(tmp_path, mutate)
    assert exports["arrays.json"]["arrays"][0]["strand"] is None
    assert exports["arrays.json"]["arrays"][0]["repeats"][0]["sequence"] is None
    assert exports["repeats.json"]["complete"] is False
    with pytest.raises(FastaValidationError, match="Incomplete"):
        sequence_input(
            json.dumps(exports["repeats.json"]),
            mode="repeats",
            molecule="DNA",
            max_bytes=100000,
        )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda p: p.update(complete=False),
        lambda p: p.update(kind="spacer"),
        lambda p: p.update(records=[]),
        lambda p: p.update(command=["run"]),
        lambda p: p["schema"].update(version="2.0.0"),
        lambda p: p["records"].append(p["records"][0]),
        lambda p: p["records"][0].update(id="x\n>injected"),
        lambda p: p["records"][0].update(sequence="ACGU"),
        lambda p: p.update(complete=1),
    ],
)
def test_import_rejects_incompatible_or_ambiguous_packets(mutation):
    packet = sequence_set("repeat", "DNA", [{"id": "r", "sequence": "ACGT"}])
    mutation(packet)
    with pytest.raises(FastaValidationError):
        sequence_input(
            json.dumps(packet), mode="repeats", molecule="DNA", max_bytes=100000
        )


@pytest.mark.parametrize(
    "text",
    [
        '{"schema":1,"schema":2}',
        '{"x":NaN}',
        '{"x":' + "[" * 1100 + "0" + "]" * 1100 + "}",
    ],
)
def test_strict_json_rejects_duplicate_nonfinite_and_deep_input(text):
    with pytest.raises(FastaValidationError):
        sequence_input(text, mode="repeats", molecule="DNA", max_bytes=100000)


def test_admission_uses_existing_limits_and_preserves_origins(
    tmp_path, settings_factory
):
    settings = settings_factory(max_records=2)
    service = JobService(settings)
    service.store.initialize()
    packet = detection_ports(tmp_path)["repeats.json"]
    submission = JobSubmission(mode="repeats", sequence=json.dumps(packet))
    with pytest.raises(FastaValidationError, match="record limit"):
        service.submit(submission, "127.0.0.1")
    packet["records"] = packet["records"][:2]
    _, job = service.submit(
        submission.model_copy(update={"sequence": json.dumps(packet)}), "127.0.0.1"
    )
    root = service.job_root(job["job_id"])
    assert job["input"]["source_ids"] == [row["id"] for row in packet["records"]]
    metadata = json.loads((root / "handoff.json").read_text())
    assert (
        metadata["origins"][packet["records"][0]["id"]]
        == packet["records"][0]["origin"]
    )
    assert "sequence" not in metadata


def test_completed_worker_registers_ports_and_excludes_them_on_failure(
    settings_factory,
):
    service = JobService(settings_factory())
    service.store.initialize()
    _, job = service.submit(
        JobSubmission(mode="detection", sequence=">source\nACGT\n"), "127.0.0.1"
    )
    assert Worker(service.settings, worker_id="exchange-test").run_once()
    finished = service.get(job["job_id"])
    assert finished["status"] == "completed"
    artifact = next(
        a for a in finished["artifacts"] if a["name"] == "input-sequences.json"
    )
    path, _ = service.artifact_file(job["job_id"], artifact["artifact_id"])
    assert (
        SequenceSet.model_validate_json(path.read_text()).records[0].sequence == "ACGT"
    )


def test_file_cap_removes_partial_json_and_cancellation_propagates(tmp_path):
    with pytest.raises(ValueError):
        write_document(tmp_path, "repeats.json", {"data": "A" * 100}, max_bytes=10)
    assert not list(tmp_path.iterdir())

    def cancelled():
        raise RuntimeError("cancelled")

    with pytest.raises(RuntimeError, match="cancelled"):
        write_document(
            tmp_path, "repeats.json", {}, max_bytes=100, checkpoint=cancelled
        )
    assert not list(tmp_path.iterdir())


def test_module_ports_and_published_schema_do_not_drift():
    root = Path(__file__).resolve().parents[2]
    assert (
        json.loads((root / "frontend/src/contracts/module-ports-v1.json").read_text())
        == INPUT_KINDS
    )
    schema = json.loads((root / "docs/schemas/sequence-set-v1.json").read_text())
    schema.pop("$schema")
    assert schema == SequenceSet.model_json_schema()


def test_prophage_validated_detection_exports_all_typed_ports(tmp_path, settings_factory):
    from types import SimpleNamespace
    from crispr_web.publication.sequence_ports import write_sequence_exports
    from crispr_web.array_exchange import array_input
    from crispr_integration.source_spacer_evidence import context_document
    case = json.loads((Path(__file__).parent / "fixtures/context-review/palindrome-search.json").read_text())
    # Complete native-shaped retained evidence, including a low-score array excluded from handoffs.
    context_document(case["search"], case["reports"], case["records"], case["regions"])
    job = tmp_path / "job"
    root = job / "results"
    (job / "input").mkdir(parents=True)
    (job / "input/normalized.fasta").write_text("".join(f">{key}\n{value}\n" for key, value in case["records"].items()))
    for i, report in enumerate(case["reports"]):
        report_path = root / f"crispridentify/record{i}/report.json"
        report_path.parent.mkdir(parents=True)
        report_path.write_text(json.dumps(report))
    claim = SimpleNamespace(mode="prophage_context", molecule="DNA", category_policy="bona_fide_possible")
    warnings = write_sequence_exports(root, job, claim, settings_factory())
    assert warnings == []
    arrays = json.loads((root / "arrays.json").read_text())
    assert len(arrays["arrays"]) == 1
    assert arrays["arrays"][0]["id"] == "native_array_0"
    assert arrays["arrays"][0]["strand"] is None
    array_input(json.dumps(arrays), max_bytes=1_000_000)
    for name, modes, wanted in [("repeats.json", ["repeats", "repeat_map", "repeat_type"], 2), ("spacers.json", ["viral_search", "spacer_association"], 1)]:
        doc = json.loads((root / name).read_text())
        assert doc["complete"] is True and len(doc["records"]) == wanted
        assert all(row["origin"]["array_id"] == "native_array_0" for row in doc["records"])
        for mode in modes:
            fasta, origins = sequence_input(json.dumps(doc), mode=mode, molecule="DNA", max_bytes=1_000_000)
            assert fasta.count(">") == wanted
            assert len(origins["origins"]) == wanted
