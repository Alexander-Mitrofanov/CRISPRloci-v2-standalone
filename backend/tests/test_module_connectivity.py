"""Producer exports -> consumer admission, across every public module and role.

Scientific execution is covered by each independent module's native/fixture
tests. This matrix tests the previously under-covered connections themselves.
"""

from dataclasses import replace
import hashlib
import json
from types import SimpleNamespace

import pytest

from crispr_web.array_exchange import array_input
from crispr_web.exchange import INPUT_KINDS, association_groups, sequence_input, sequence_set
from crispr_web.execution.input_integrity import _verify_job_inputs
from crispr_web.fasta import FastaValidationError, normalize_fasta
from crispr_web.models import JobSubmission, ScientificMode
from crispr_web.publication.sequence_ports import write_sequence_exports
from crispr_web.service import JobService
from .test_array_compare import comparison_settings, rich_packet
from .test_sequence_exchange import FIXTURE, detection_ports


DETECTORS = {"detection", "orientation", "reconstruction", "leader", "repeat_context", "prophage_context", "loci"}
EDGES = [
    (producer, consumer)
    for producer, kind in INPUT_KINDS.items()
    for consumer, expected in INPUT_KINDS.items()
    if ("spacer" if producer == "array_compare" else kind) == expected
]


def normalize(text, molecule="DNA"):
    return normalize_fasta(text, max_records=1000, max_record_bases=10000,
                           max_total_bases=100000, max_header_characters=200, molecule=molecule)


@pytest.fixture
def admission_settings(settings_factory, tmp_path):
    """Admission-only reference identities; these are never scientific fixtures."""
    settings = comparison_settings(settings_factory, tmp_path)
    path = settings.scientific_release_manifest
    manifest = json.loads(path.read_text())
    for name, version in [("spacepharer", "2.1.0"), ("repeattyper", "2.0.0"), ("crisprmap", "2.0.0a1")]:
        manifest["scientific_artifacts"][name] = {
            **manifest["scientific_artifacts"]["crisprspacer"], "name": name,
            "display_name": name, "version": version,
        }
    manifest["schema_version"] = "10.0.0"
    for name in ("clinker", "phigaro", "isescan", "defensefinder", "rnaalifold"):
        manifest["scientific_artifacts"][name] = {
            **manifest["scientific_artifacts"]["crisprspacer"], "name": name, "display_name": name, "version": "1.0.0"}
    manifest["scientific_artifacts"]["crisprspacer"]["version"] = "0.3.0"
    path.write_text(json.dumps(manifest))
    reference = tmp_path / "reference.json"
    reference.write_text("{}")
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    return replace(settings, scientific_release_manifest_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                   spacer_database_path=tmp_path, spacer_database_sha256=digest,
                   viral_reference_path=tmp_path, viral_reference_sha256=digest,
                   association_reference_path=tmp_path, association_reference_sha256=digest,
                   repeat_reference_path=reference, repeat_reference_sha256=digest)


def export_input(tmp_path, settings, mode, *, rna=False):
    root = tmp_path / "results"
    job_root = tmp_path / "producer"
    (job_root / "input").mkdir(parents=True)
    root.mkdir(exist_ok=True)
    if mode == "array_compare":
        packet = rich_packet(tmp_path)
        text, mapping = array_input(json.dumps(packet), max_bytes=1000000)
        (job_root / "array-map.json").write_text(json.dumps(mapping))
    elif mode in DETECTORS:
        detection_ports(tmp_path)
        text = (FIXTURE / "source.fasta").read_text()
        if mode == "loci":
            original = normalize(text).records[0]
            report = json.loads((FIXTURE / "identify.json").read_text())
            report["source"]["id"] = "f000001.s000001"
            native = root / "tracrrna/work/f000001.s000001/identify/source/report.json"
            native.parent.mkdir(parents=True)
            native.write_text(json.dumps(report))
            (root / "tracrrna/report.json").write_text(json.dumps({
                "status": "complete", "records": [{"key": "f000001.s000001", "source": {
                    "id": original.identifier, "length": len(original.sequence),
                    "sha256": hashlib.sha256(original.sequence.encode()).hexdigest(),
                }}],
            }))
    elif mode == "rna_family":
        text = ">family_a\n" + "ACGT" * 10 + "\n>family_b\n" + "ACGT" * 9 + "ACGA\n"
    else:
        text = ">occurrence_\n" + ("ACGU" if rna else "ACGT") * 6 + "\n"
    molecule = "RNA" if rna else "DNA"
    normalized = normalize(text, molecule)
    (job_root / "input/normalized.fasta").write_text(normalized.combined)
    claim = SimpleNamespace(mode=mode, molecule=molecule, category_policy="bona_fide_possible")
    assert write_sequence_exports(root, job_root, claim, settings) == []
    return root


def submit(service, consumer, text):
    _, job = service.submit(JobSubmission(mode=consumer, sequence=text,
        association_grouping="per_source" if consumer == "spacer_association" else None), "192.0.2.1")
    root = service.job_root(job["job_id"])
    assert _verify_job_inputs(service.settings, root, service.store.get_job(job["job_id"]))
    return root


@pytest.mark.parametrize("producer,consumer", EDGES)
def test_every_input_export_round_trips_to_compatible_consumers(tmp_path, admission_settings, producer, consumer):
    root = export_input(tmp_path, admission_settings, producer)
    text = (root / "input-sequences.json").read_text()
    packet = json.loads(text)
    if consumer == "locus_compare" and len(packet["records"]) < 2:
        with pytest.raises(FastaValidationError, match="at least two"):
            submit(JobService(admission_settings), consumer, text)
        return
    # Non-detection spacer inputs need an explicit user group instead of a
    # fabricated source annotation.
    if consumer == "spacer_association" and producer != "array_compare":
        _, job = JobService(admission_settings).submit(JobSubmission(mode=consumer, sequence=text,
            association_grouping="single_set"), "192.0.2.1")
        assert job["input"]["source_ids"] == [r["id"] for r in packet["records"]]
    else:
        accepted = submit(JobService(admission_settings), consumer, text)
        assert normalize((accepted / "input/normalized.fasta").read_text()).source_ids == tuple(r["id"] for r in packet["records"])


@pytest.mark.parametrize("producer", sorted(DETECTORS))
@pytest.mark.parametrize("filename,consumer", [("repeats.json", "repeats"), ("repeats.json", "repeat_map"),
    ("repeats.json", "repeat_type"), ("spacers.json", "viral_search"),
    ("spacers.json", "spacer_association"), ("arrays.json", "array_compare")])
def test_every_detector_output_is_directly_admitted(tmp_path, admission_settings, producer, filename, consumer):
    root = export_input(tmp_path, admission_settings, producer)
    submit(JobService(admission_settings), consumer, (root / filename).read_text())


def test_comparison_spacers_keep_sources_through_another_handoff(tmp_path, admission_settings):
    root = export_input(tmp_path, admission_settings, "array_compare")
    text = (root / "spacers.json").read_text()
    packet = json.loads(text)
    fasta, metadata = sequence_input(text, mode="spacer_association", molecule="DNA", max_bytes=1000000)
    records = normalize(fasta)
    groups = association_groups(records.source_ids, "per_source", metadata)
    assert len(packet["records"]) == 5  # All occurrences, excluding the deletion slot.
    assert [len(g["members"]) for g in groups["groups"]] == [3, 2]
    accepted = submit(JobService(admission_settings), "spacer_association", text)
    output = tmp_path / "second"
    output.mkdir()
    assert write_sequence_exports(output, accepted, SimpleNamespace(mode="spacer_association", molecule="DNA"), admission_settings) == []
    assert json.loads((output / "input-sequences.json").read_text()) == packet


def test_rna_output_only_connects_to_rna_repeat_analysis(tmp_path, admission_settings):
    root = export_input(tmp_path, admission_settings, "repeats", rna=True)
    text = (root / "input-sequences.json").read_text()
    assert sequence_input(text, mode="repeats", molecule="RNA", max_bytes=1000000)
    for mode in ["repeat_map", "repeat_type", "repeats"]:
        with pytest.raises(FastaValidationError, match="Select RNA"):
            sequence_input(text, mode=mode, molecule="DNA", max_bytes=1000000)


def test_matrix_covers_all_modes():
    assert set(INPUT_KINDS) == {mode.value for mode in ScientificMode}
    assert {p for p, _ in EDGES} == set(INPUT_KINDS)


@pytest.mark.parametrize("mode", sorted(INPUT_KINDS))
@pytest.mark.parametrize("kind", ["target", "repeat", "spacer", "rna"])
def test_incompatible_roles_are_rejected_for_every_mode(mode, kind):
    text = json.dumps(sequence_set(kind, "DNA", [{"id": "record_", "sequence": "ACGT"}]))
    if INPUT_KINDS[mode] == kind:
        assert sequence_input(text, mode=mode, molecule="DNA", max_bytes=10000)
    else:
        with pytest.raises(FastaValidationError):
            sequence_input(text, mode=mode, molecule="DNA", max_bytes=10000)


def test_canonical_ids_do_not_lose_origins_during_grouping(tmp_path, admission_settings):
    packet = sequence_set("spacer", "DNA", [{
        "id": key, "sequence": "ACGT", "origin": {
            "source_id": "s" * 240, "source_sha256": "a" * 64,
            "array_id": "array1", "ordinal": index, "start": index * 4, "end": index * 4 + 3,
        },
    } for index, key in enumerate(["repeat_", "repeat-"], 1)])
    root = submit(JobService(admission_settings), "spacer_association", json.dumps(packet))
    groups = json.loads((root / "association-groups.json").read_text())
    assert groups["groups"] == [{"id": "source1", "members": ["repeat_", "repeat-"]}]


def test_http_download_to_submission_round_trip_has_independent_capabilities(settings_factory):
    from fastapi.testclient import TestClient
    from crispr_web.api import create_app
    from crispr_web.worker import Worker

    settings = settings_factory()
    with TestClient(create_app(settings)) as client:
        first = client.post("/api/v1/jobs", json={"mode": "detection", "sequence": ">source_\nACGT\n"})
        assert first.status_code == 202
        original = first.json()
        jid = original["job"]["job_id"]
        headers = {"Authorization": "Bearer " + original["access_token"]}
        assert Worker(settings, worker_id="connection").run_once()
        result = client.get("/api/v1/jobs/" + jid, headers=headers).json()
        assert result["status"] == "completed"
        port = next(a for a in result["artifacts"] if a["name"] == "input-sequences.json")
        assert client.get(port["download_url"]).status_code in {401, 404}
        downloaded = client.get(port["download_url"], headers=headers)
        assert downloaded.status_code == 200
        second = client.post("/api/v1/jobs", json={"mode": "detection", "sequence": downloaded.text})
        assert second.status_code == 202
        next_job = second.json()
        assert next_job["access_token"] != original["access_token"]
        assert client.get("/api/v1/jobs/" + next_job["job"]["job_id"], headers=headers).status_code == 404
        assert Worker(settings, worker_id="connection").run_once()
        next_headers = {"Authorization": "Bearer " + next_job["access_token"]}
        final = client.get("/api/v1/jobs/" + next_job["job"]["job_id"], headers=next_headers).json()
        assert final["status"] == "completed"
        next_port = next(a for a in final["artifacts"] if a["name"] == "input-sequences.json")
        assert client.get(next_port["download_url"], headers=next_headers).json() == downloaded.json()
