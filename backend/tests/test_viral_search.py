"""Real synthetic BLAST pipeline, independent replay, mutation and public export gates."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from crispr_integration.runner import RunnerConfig
from crispr_integration.viral import run_viral_pipeline
from crispr_web.db import ClaimedJob
from crispr_web.execution.validation_pipeline import _validated_pipeline_status
from crispr_web.publication.viral import write_viral_exports
from crispr_web.summaries.viral import viral_summary
from crispr_web.models import JobSubmission
from crispr_web.service import JobService, ReferenceUnavailableError


def configured_settings(settings_factory, tmp_path, **overrides):
    names = (
        "crispridentify",
        "spacerplacer",
        "crispr_evor",
        "casandra",
        "crispr_tracrrna",
        "crisprspacer",
    )
    artifacts = {
        name: {
            "name": name,
            "display_name": name,
            "version": "0.2.0",
            "kind": "git_source",
            "path": "/private/" + name,
            "sha256": "a" * 64,
            "sha256_scope": "source_inventory_v1",
            "repository": "https://github.com/example/" + name + ".git",
            "commit": "b" * 40,
        }
        for name in names
    }
    artifacts["integration"] = {
        "name": "integration",
        "display_name": "Integration",
        "version": "0.6.0",
        "kind": "wheel",
        "path": "/private/integration.whl",
        "sha256": "c" * 64,
        "sha256_scope": "file_bytes",
    }
    document = {
        "schema": "crispr-web-scientific-release",
        "schema_version": "5.0.0",
        "release_id": "d" * 64,
        "scientific_artifacts": artifacts,
    }
    payload = json.dumps(document).encode()
    manifest = tmp_path / "spacer-release.json"
    manifest.write_bytes(payload)
    return settings_factory(
        scientific_release_manifest=manifest,
        scientific_release_manifest_sha256=hashlib.sha256(payload).hexdigest(),
        **overrides,
    )


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "tools/CRISPRspacer"
NATIVE = (
    sys.executable,
    "-c",
    f"import sys;sys.path.insert(0,{str(SOURCE / 'src')!r});from crisprspacer.cli import main;sys.exit(main())",
)


@pytest.fixture(scope="module")
def viral_reference(tmp_path_factory):
    if any(shutil.which(n) is None for n in ("blastn", "makeblastdb", "blastdbcmd")):
        pytest.skip("Requires BLAST+ for synthetic reference")
    root = tmp_path_factory.mktemp("viral-reference")
    reference = root / "reference"
    subprocess.run(
        [
            *NATIVE,
            "viral-index",
            "--sources",
            str(SOURCE / "examples/synthetic-viral-targets.fasta"),
            "--name",
            "synthetic",
            "--source-release",
            "test",
            "--source-url",
            "https://example.org/synthetic",
            "--output",
            str(reference),
        ],
        check=True,
        capture_output=True,
    )
    return reference, hashlib.sha256(
        (reference / "reference.json").read_bytes()
    ).hexdigest()


@pytest.fixture
def viral_pipeline(tmp_path, viral_reference):
    reference, digest = viral_reference
    source = tmp_path / "input"
    source.mkdir()
    text = (SOURCE / "examples/synthetic-viral-spacers.fasta").read_text().splitlines()
    records = {text[i][1:]: text[i + 1] for i in range(0, len(text), 2)}
    for i, (key, value) in enumerate(records.items()):
        (source / f"{i:06d}.fasta").write_text(f">{key}\n{value}\n")
    out = tmp_path / "output"
    result = run_viral_pipeline(
        source,
        out,
        runner_config=RunnerConfig(
            ("unused",),
            ("unused",),
            ("unused",),
            crisprspacer=NATIVE,
            blastn=(shutil.which("blastn"),),
        ),
        reference_path=reference,
        reference_sha256=digest,
        timeout_seconds=60,
    )
    assert result["status"] == "completed"
    claim = ClaimedJob(
        "a" * 32,
        "viral_search",
        1,
        0,
        "bona_fide_possible",
        1,
        True,
        True,
        (),
        viral_reference_sha256=digest,
    )
    bindings = {
        k: (len(v), hashlib.sha256(v.encode()).hexdigest()) for k, v in records.items()
    }
    return claim, out, bindings, source, reference


def validated(data):
    claim, out, bindings, source, reference = data
    return _validated_pipeline_status(
        claim, out, bindings, input_root=source, viral_reference=reference
    )


def rehash(out):
    root = out / "viral"
    p = root / "manifest.json"
    m = json.loads(p.read_text())
    for name in m["files"]:
        data = (root / name).read_bytes()
        m["files"][name] = {
            "sha256": hashlib.sha256(data).hexdigest(),
            "bytes": len(data),
        }
    p.write_text(json.dumps(m))


def test_real_reverse_pipeline_and_safe_complete_exports(viral_pipeline):
    assert validated(viral_pipeline) == "completed"
    out = viral_pipeline[1]
    preview = viral_summary(out, max_rows=1)
    assert (
        preview["counts"]["matched_subject_accessions"] == 3
        and preview["matches_truncated"]
    )
    write_viral_exports(out, {"viral_search": preview})
    public = json.loads((out / "viral-matches.json").read_text())
    assert (
        len(public["matches"]) == 3 and "blast" not in public and "limits" not in public
    )
    assert "/tmp/" not in json.dumps(public)
    assert len((out / "viral-matches.tsv").read_text().splitlines()) == 4


@pytest.mark.parametrize(
    "mutation",
    [
        "strand",
        "position",
        "title",
        "omitted_hit",
        "omitted_raw",
        "counts",
        "queries",
        "reference",
        "policy",
        "source",
        "manifest",
        "stage",
        "symlink",
        "boolean_position",
    ],
)
def test_mutations_fail_independent_validation(viral_pipeline, mutation):
    claim, out, bindings, source, reference = viral_pipeline
    p = out / "viral/result.json"
    doc = json.loads(p.read_text())
    if mutation == "strand":
        doc["matches"][0]["relative_strand"] = "-"
    elif mutation == "position":
        doc["matches"][0]["subject_start"] += 1
    elif mutation == "boolean_position":
        doc["matches"][0]["query_start"] = True
    elif mutation == "title":
        doc["matches"][0]["subject_title"] = "Forged virus"
    elif mutation == "omitted_hit":
        doc["matches"].pop()
    elif mutation == "omitted_raw":
        (out / "viral/blast-hsps.tsv").write_text("")
    elif mutation == "counts":
        doc["counts"]["raw_hsps"] += 1
    elif mutation == "queries":
        doc["queries"][0]["returned_subjects"] += 1
    elif mutation == "reference":
        doc["reference"]["reference_sha256"] = "a" * 64
    elif mutation == "policy":
        doc["policy"]["max_substitutions"] = 1
    elif mutation == "source":
        next(source.iterdir()).write_text(">forged\nACGT\n")
    elif mutation == "stage":
        stage = out / "integration_result.json"
        d = json.loads(stage.read_text())
        d["stages"]["crisprspacer"]["reference_sha256"] = "b" * 64
        stage.write_text(json.dumps(d))
    p.write_text(json.dumps(doc))
    rehash(out)
    if mutation == "manifest":
        (out / "viral/matches.tsv").write_text("corrupt")
    if mutation == "symlink":
        raw = out / "viral/blast-hsps.tsv"
        outside = out / "outside"
        raw.rename(outside)
        raw.symlink_to(outside)
    assert validated(viral_pipeline) is None
    assert (
        validated(
            (
                replace(claim, viral_reference_sha256="b" * 64),
                out,
                bindings,
                source,
                reference,
            )
        )
        is None
    )


def test_unconfigured_search_is_unavailable_before_admission(settings_factory):
    service = JobService(settings_factory())
    with pytest.raises(ReferenceUnavailableError):
        service.submit(
            JobSubmission(mode="viral_search", sequence=">q\nACGT\n"), "127.0.0.1"
        )


def test_viral_reference_and_policy_survive_queue_claim(
    settings_factory, tmp_path, viral_reference
):
    reference, digest = viral_reference
    # The scientific manifest is a test fixture; identity version follows the pinned tool.
    settings = configured_settings(settings_factory, tmp_path)
    manifest = settings.scientific_release_manifest
    data = json.loads(manifest.read_text())
    data["scientific_artifacts"]["crisprspacer"]["version"] = "0.3.0"
    manifest.write_text(json.dumps(data))
    settings = replace(
        settings,
        scientific_release_manifest_sha256=hashlib.sha256(
            manifest.read_bytes()
        ).hexdigest(),
        viral_reference_path=reference,
        viral_reference_sha256=digest,
    )
    service = JobService(settings)
    _, job = service.submit(
        JobSubmission(
            mode="viral_search",
            sequence=">q\nACGTACGTACGTACGTACGT\n",
            viral_max_mismatches=0,
        ),
        "127.0.0.1",
    )
    assert job["options"]["viral_max_mismatches"] == 0
    claim = service.store.claim_next("test-worker")
    assert claim.viral_reference_sha256 == digest and claim.viral_max_mismatches == 0
