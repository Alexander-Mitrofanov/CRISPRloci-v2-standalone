"""Provisioning gates and bounded admissions for the five context analyses."""
from dataclasses import replace
import json

import pytest
from fastapi.testclient import TestClient

from crispr_web.api import create_app
from crispr_web.context_modes import MODE_ARTIFACTS, context_mode_available
from crispr_web.exchange import sequence_set
from crispr_web.fasta import FastaValidationError
from crispr_web.models import JobSubmission
from crispr_web.service import JobService, ReferenceUnavailableError
from .test_module_connectivity import admission_settings as admission_settings


@pytest.mark.parametrize("mode,requirements", MODE_ARTIFACTS.items())
def test_every_required_identity_is_needed(mode, requirements):
    assert context_mode_available(mode, dict.fromkeys(requirements))
    for missing in requirements:
        assert not context_mode_available(mode, {key: {} for key in requirements if key != missing})


@pytest.mark.parametrize("mode", MODE_ARTIFACTS)
def test_unprovisioned_context_modes_are_not_admitted(settings_factory, mode):
    with pytest.raises(ReferenceUnavailableError, match="not provisioned"):
        JobService(settings_factory()).submit(JobSubmission(mode=mode, sequence=">a\nACGT\n"), "192.0.2.1")


def test_only_provisioned_modes_are_advertised(settings_factory, admission_settings):
    manifest_bytes = admission_settings.scientific_release_manifest.read_bytes()
    with TestClient(create_app(settings_factory())) as client:
        assert not set(MODE_ARTIFACTS).intersection(client.get("/api/v1/config").json()["modes"])
    admission_settings.scientific_release_manifest.write_bytes(manifest_bytes)
    with TestClient(create_app(admission_settings)) as client:
        assert set(MODE_ARTIFACTS).issubset(client.get("/api/v1/config").json()["modes"])


@pytest.mark.parametrize("mode", ["locus_compare", "rna_family"])
def test_comparisons_need_multiple_sources(admission_settings, mode):
    with pytest.raises(FastaValidationError, match="at least two"):
        JobService(admission_settings).submit(JobSubmission(mode=mode, sequence=">a\n" + "A"*40), "192.0.2.1")


@pytest.mark.parametrize("mode,sequence", [
    ("locus_compare", ">a\n" + "A"*50001 + "\n>b\n" + "C"*40),
    ("locus_compare", "".join(f">a{i}\nACGT\n" for i in range(9))),
    ("rna_family", ">a\n" + "A"*501 + "\n>b\n" + "C"*40),
    ("rna_family", ">a\n" + "A"*29 + "\n>b\n" + "C"*40),
    ("rna_family", ">a\n" + "A"*40 + "\n>b\n" + "A"*40),
    ("rna_family", ">a\n" + "N"*40 + "\n>b\n" + "A"*40),
    ("rna_family", "".join(f">a{i}\n" + "ACGT"*10 + "\n" for i in range(51))),
    ("defense_context", "".join(f">a{i}\n" + "A"*700000 + "\n" for i in range(3))),
])
def test_context_specific_bounds_survive_larger_operator_defaults(admission_settings, mode, sequence):
    settings = replace(admission_settings, max_total_bases=3_000_000, max_record_bases=2_000_000, max_records=100)
    with pytest.raises(FastaValidationError):
        JobService(settings).submit(JobSubmission(mode=mode, sequence=sequence), "192.0.2.1")


@pytest.mark.parametrize("molecule,base", [("DNA", "T"), ("RNA", "U")])
def test_rna_family_fasta_and_typed_port_admit_correct_alphabet(admission_settings, molecule, base):
    packet = sequence_set("rna", molecule, [{"id": "a", "sequence": "ACG"+base*37}, {"id": "b", "sequence": "GCA"+base*37}])
    _, job = JobService(admission_settings).submit(JobSubmission(mode="rna_family", molecule=molecule, sequence=json.dumps(packet)), "192.0.2.1")
    assert job["input"]["record_count"] == 2


def test_repeat_port_cannot_claim_rna_family_homology(admission_settings):
    packet = sequence_set("repeat", "RNA", [{"id": "a", "sequence": "ACGU"*10}])
    with pytest.raises(FastaValidationError):
        JobService(admission_settings).submit(JobSubmission(mode="rna_family", molecule="RNA", sequence=json.dumps(packet)), "192.0.2.1")


@pytest.mark.parametrize("molecule,base", [("DNA", "T"), ("RNA", "U")])
def test_rna_family_worker_preserves_admitted_molecule(admission_settings, molecule, base):
    from crispr_web.worker import Worker

    service = JobService(admission_settings)
    _, job = service.submit(
        JobSubmission(mode="rna_family", molecule=molecule,
                      sequence=f">a\n{'ACG' + base * 37}\n>b\n{'GCA' + base * 37}\n"),
        "192.0.2.1",
    )
    claim = service.store.claim_next("rna-family-command-regression")
    assert claim is not None and claim.job_id == job["job_id"]
    worker = Worker(admission_settings, store=service.store, worker_id="rna-family-command-regression")
    command = worker.command_for(claim, service.job_root(claim.job_id) / "output")
    assert command[command.index("--tools") + 1] == "rna-family"
    assert command.count("--molecule") == 1
    assert command[command.index("--molecule") + 1] == molecule
