"""Release guards for independently pinned context engines and model archives."""

import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from runtime import context_tools
from runtime.common import load_lock
import provision_tools


def archive(path, entries):
    with tarfile.open(path, "w:gz") as handle:
        for name, data in entries:
            member = tarfile.TarInfo(name)
            member.size = len(data)
            handle.addfile(member, io.BytesIO(data))


def pin(monkeypatch, path, entries):
    expected = b"".join(data for _, data in sorted(entries))
    monkeypatch.setattr(context_tools, "PVOG_ARCHIVE", context_tools.digest(path))
    monkeypatch.setattr(context_tools, "PVOG_HMM", hashlib.sha256(expected).hexdigest())
    monkeypatch.setattr(context_tools, "PVOG_COUNT", len(entries))
    monkeypatch.setattr(context_tools, "PVOG_BYTES", len(expected))
    return expected


def test_pvog_single_pass_reorders_and_preserves_suffix_profiles(tmp_path, monkeypatch):
    entries = [
        ("AllvogHMMprofiles/VOG2b.hmm", b"profile2b\n"),
        ("AllvogHMMprofiles/VOG1.hmm", b"profile1\n"),
    ]
    source = tmp_path / "models.tar.gz"
    archive(source, entries)
    expected = pin(monkeypatch, source, entries)
    target = tmp_path / "models.hmm"
    context_tools.concatenate_pvogs(source, target)
    assert target.read_bytes() == expected
    assert not list(tmp_path.glob(".pvogs-*"))


@pytest.mark.parametrize(
    "name",
    [
        "../escape.hmm",
        "/absolute.hmm",
        "AllvogHMMprofiles/../escape.hmm",
        "AllvogHMMprofiles/unexpected.txt",
    ],
)
def test_reference_tar_paths_are_not_extracted(tmp_path, monkeypatch, name):
    entries = [(name, b"model")]
    source = tmp_path / "models.tar.gz"
    archive(source, entries)
    pin(monkeypatch, source, entries)
    with pytest.raises(ValueError):
        context_tools.concatenate_pvogs(source, tmp_path / "models.hmm")
    assert not (tmp_path / "models.hmm").exists()


def test_duplicate_model_members_rejected(tmp_path, monkeypatch):
    entries = [("AllvogHMMprofiles/VOG1.hmm", b"a")] * 2
    source = tmp_path / "models.tar.gz"
    archive(source, entries)
    pin(monkeypatch, source, entries)
    with pytest.raises(ValueError, match="Duplicate"):
        context_tools.concatenate_pvogs(source, tmp_path / "models.hmm")


def test_model_digest_failure_never_publishes_destination(tmp_path, monkeypatch):
    entries = [("AllvogHMMprofiles/VOG1.hmm", b"a")]
    source = tmp_path / "models.tar.gz"
    archive(source, entries)
    pin(monkeypatch, source, entries)
    monkeypatch.setattr(context_tools, "PVOG_HMM", "0" * 64)
    with pytest.raises(ValueError, match="identity"):
        context_tools.concatenate_pvogs(source, tmp_path / "models.hmm")
    assert not (tmp_path / "models.hmm").exists()


def test_native_model_receipt_is_required_and_bound(tmp_path):
    executable = tmp_path / "hmmsearch"
    executable.write_bytes(b"native")
    model_dir = tmp_path / "models"
    model_dir.mkdir()
    tool = {
        "name": "DefenseFinder-v2",
        "native_tools": {
            "hmmsearch": {
                "version": "3.4",
                "path": str(executable),
                "sha256": context_tools.digest(executable),
            },
            "defense-finder-models": {
                "path": str(model_dir),
                "commit": context_tools.DEFENSE_COMMIT,
                "manifest_sha256": context_tools.DEFENSE_MANIFEST,
            },
        },
    }
    context_tools.verify_context_native_receipt(tool, tmp_path)
    executable.write_bytes(b"mutated")
    with pytest.raises(ValueError, match="bytes changed"):
        context_tools.verify_context_native_receipt(tool, tmp_path)
    tool["native_tools"].pop("hmmsearch")
    with pytest.raises(ValueError, match="receipt missing"):
        context_tools.verify_context_native_receipt(tool, tmp_path)


def test_context_runtime_dependencies_do_not_merge_sp_model_or_enable_updates():
    assert provision_tools.CONTEXT_CONDA["ISEScan-v2"] == [
        "python=3.11",
        "pip",
        "numpy=1.26.4",
        "scipy=1.13.1",
        "biopython=1.85",
        "fastcluster=1.2.6",
        "fraggenescan=1.32",
        "hmmer=3.4",
        "blast=2.16.0",
    ]
    assert provision_tools.CONTEXT_CONDA["RNAalifold-v2"][-1] == "mafft=7.526"
    assert (
        provision_tools.CONTEXT_REQUIREMENTS["DefenseFinder-v2"]
        == "requirements-context.lock"
    )
    assert "SpacerPlacer" not in provision_tools.CONTEXT_TOOLS
    assert "evOr" not in provision_tools.CONTEXT_TOOLS
    for name, filename in provision_tools.CONTEXT_REQUIREMENTS.items():
        assert (Path(__file__).resolve().parents[1] / "tools" / name / filename).is_file()


def test_lock_loader_accepts_all_seventeen_pinned_tools(tmp_path):
    root = Path(__file__).resolve().parents[1]
    lock = json.loads((root / "tools/tools.lock.json").read_text())
    for name in provision_tools.CONTEXT_TOOLS:
        if name not in {item["name"] for item in lock["tools"]}:
            lock["tools"].append(
                {
                    "name": name,
                    "path": f"tools/{name}",
                    "repository": f"https://github.com/Alexander-Mitrofanov/{name}.git",
                    "commit": "a" * 40,
                    "version": "1.0.0",
                    "entrypoint": name,
                }
            )
    path = tmp_path / "tools.lock.json"
    path.write_text(json.dumps(lock))
    assert len(load_lock(path)["tools"]) == 17


def test_native_context_gate_cannot_pass_fabricated_zero_as_positive():
    from runtime.context_gate import check_summary

    summary = {
        "insertion_context": {
            "available": True,
            "counts": {"elements": 0},
            "elements": [],
        }
    }
    check_summary("insertion_context", "insertion-negative", summary)
    with pytest.raises(RuntimeError, match="control count"):
        check_summary("insertion_context", "insertion-positive", summary)


def test_rna_gate_requires_distinct_sources_and_preserves_no_activity_claim():
    from runtime.context_gate import check_summary

    summary = {
        "rna_family": {
            "available": True,
            "counts": {"sources": 4, "unique_sequences": 3, "mfe_pairs": 12},
            "engine": {"version": "2.7.0"},
            "policy": {"activity": "not_assessed"},
        }
    }
    check_summary("rna_family", "rna-positive", summary)
    summary["rna_family"]["counts"]["unique_sequences"] = 4
    with pytest.raises(RuntimeError, match="deduplication"):
        check_summary("rna_family", "rna-positive", summary)


def test_gate_rejects_unfit_comparisons_and_rna_for_genomic_workflows():
    from runtime.context_gate import rejected_submissions

    cases = rejected_submissions({key: [] for key in ("clinker", "phigaro", "isescan", "defensefinder", "rnaalifold")})
    assert len(cases) == 5
    assert {case["mode"] for case in cases} == {"locus_compare", "prophage_context", "insertion_context", "defense_context", "rna_family"}
    assert all(case.get("molecule") == "RNA" for case in cases if case["mode"] != "locus_compare")


def test_context_activation_requires_positive_negative_and_rejection_controls():
    from runtime.context_gate import require_coverage

    gate = {"cases": [{"mode": "insertion_context", "fixture": "insertion-positive", "status": "completed"}], "context_predispatch_rejections": 1}
    with pytest.raises(ValueError, match="missing native controls"):
        require_coverage({"isescan": {}}, gate)
    gate["cases"].append({"mode": "insertion_context", "fixture": "insertion-negative", "status": "completed"})
    require_coverage({"isescan": {}}, gate)
    gate["context_predispatch_rejections"] = 0
    with pytest.raises(ValueError, match="input rejections"):
        require_coverage({"isescan": {}}, gate)


def test_phigaro_gate_keeps_unassessed_short_input_separate_from_negative():
    from runtime.context_gate import check_summary

    doc = {"available": True, "counts": {"assessed_sources": 0, "genes": 0, "regions": 0}, "parameters": {"gene_mode": "single"}}
    check_summary("prophage_context", "prophage-ineligible", {"prophage_context": doc})
    with pytest.raises(RuntimeError, match="assessment"):
        check_summary("prophage_context", "prophage-negative", {"prophage_context": doc})
    doc["counts"].update(assessed_sources=1, genes=1)
    check_summary("prophage_context", "prophage-negative", {"prophage_context": doc})
    doc["parameters"]["gene_mode"] = "meta"
    with pytest.raises(RuntimeError, match="gene policy"):
        check_summary("prophage_context", "prophage-negative", {"prophage_context": doc})


def test_empty_gene_controls_are_short_and_keep_comparison_record_count(tmp_path):
    from runtime.context_gate import fixtures

    toolset = tmp_path / "toolset"
    genomic = toolset / "sources/CRISPRtracrRNA-v3/examples/complete.fasta"
    genomic.parent.mkdir(parents=True)
    genomic.write_text(">fixture\n" + "ACGT" * 2500 + "\n")
    defense = toolset / "sources/DefenseFinder-v2/tests/data/context_fragment.fna"
    defense.parent.mkdir(parents=True)
    defense.write_text(">fixture\nACGT\n")
    run_root = tmp_path / "run"
    run_root.mkdir()
    paths, cases = fixtures(toolset, run_root, {"clinker": (), "defensefinder": ()})
    # The pinned Pyrodigal 3.7.1 meta caller predicts one partial CDS in A*1000
    # or C*1000. These real-tested twelve-base controls produce no CDS with it.
    for name, sources in (("locus-negative", 2), ("defense-negative", 1)):
        records = [record.splitlines() for record in paths[name].read_text().split(">")[1:]]
        assert len(records) == sources
        assert all(len("".join(record[1:])) == 12 for record in records)
        assert all(set("".join(record[1:])) == set("ACGT") for record in records)
    assert ("locus_compare", "locus-negative", "II") in cases
    assert ("defense_context", "defense-negative", "II") in cases


@pytest.mark.parametrize(
    "mode,fixture,counts",
    [
        ("locus_compare", "locus-negative", {"loci": 2, "genes": 0, "links": 0}),
        (
            "defense_context",
            "defense-negative",
            {"sources": 1, "genes": 0, "systems": 0, "components": 0, "hits": 0},
        ),
    ],
)
def test_empty_gene_gate_still_rejects_partial_gene_calls(mode, fixture, counts):
    from runtime.context_gate import check_summary

    doc = {
        "available": True,
        "counts": counts,
        "policy": {"functional_activity": "not_assessed"},
    }
    check_summary(mode, fixture, {mode: doc})
    counts["genes"] = 1
    with pytest.raises(RuntimeError):
        check_summary(mode, fixture, {mode: doc})


def test_future_managed_launcher_removes_mafft_override_without_changing_caller(tmp_path, monkeypatch):
    import os
    import subprocess
    script = tmp_path / "probe.py"
    script.write_text("import json,os;print(json.dumps({k:os.environ.get(k) for k in ['MAFFT_BINARIES','TMPDIR','MAFFT_TMPDIR']}))")
    target = tmp_path / "runner"
    provision_tools.wrapper(target, [sys.executable, script], Path(sys.executable).parent)
    monkeypatch.setenv("MAFFT_BINARIES", str(tmp_path / "override"))
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    monkeypatch.setenv("MAFFT_TMPDIR", str(tmp_path / "mafft-temp"))
    result = subprocess.run([str(target)], check=True, text=True, capture_output=True)
    assert json.loads(result.stdout) == {"MAFFT_BINARIES": None, "TMPDIR": str(tmp_path), "MAFFT_TMPDIR": str(tmp_path / "mafft-temp")}
    assert os.environ["MAFFT_BINARIES"] == str(tmp_path / "override")
