"""Reconstruction, evidence, exports and paired prediction audits for round two."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from copy import deepcopy
from pathlib import Path
from urllib.parse import unquote

from crisprtracrrna.adapters import sha256_file
from crisprtracrrna.models import Interval
from crisprtracrrna.pipeline import atomic_json
from crisprtracrrna.sequence import read_fasta, reverse_complement

ROOT = Path(__file__).resolve().parents[1]
ROUND = ROOT / "benchmarks/council_round2"


def load(path):
    return json.loads(path.read_text())


def projection(predictions, experimental=False):
    output = []
    for original in predictions:
        if experimental and original.get("admission_branch") == "short_rescue":
            continue
        row = deepcopy(original)
        for key in ("admission_branch", "evidence_class"):
            row.pop(key, None)
        if experimental:
            row.pop("export_id", None)
        if row.get("interaction"):
            row["interaction"].pop("target_paired_positions", None)
            row["interaction"].pop("repeat_paired_positions", None)
        output.append(row)
    return output


def checked_interval(value, length):
    interval = Interval.from_public(value["start"], value["end"])
    assert interval.end <= length
    assert value["length"] == interval.length
    return interval


def orient(sequence, strand):
    return sequence if strand == "+" else reverse_complement(sequence)


def prediction_audit(p, source):
    interval = checked_interval(p["interval"], len(source))
    expected = orient(source[interval.start : interval.end], p["strand"])
    assert p["sequence"] == expected

    def genomic(local):
        return (
            Interval(interval.start + local.start, interval.start + local.end)
            if p["strand"] == "+"
            else Interval(interval.end - local.end, interval.end - local.start)
        )

    anti = p.get("anti_repeat")
    if anti:
        anti_local = checked_interval(anti["local_interval"], len(expected))
        anti_genomic = checked_interval(anti["interval"], len(source))
        assert genomic(anti_local) == anti_genomic
        assert anti["sequence"] == expected[anti_local.start : anti_local.end]
        assert anti["sequence"] == orient(
            source[anti_genomic.start : anti_genomic.end], p["strand"]
        )
    for hit in p.get("model_hits", []) + p.get("terminator_hits", []):
        local = checked_interval(hit["interval"], len(expected))
        global_interval = checked_interval(hit["genomic_interval"], len(source))
        assert genomic(local) == global_interval
        assert orient(expected[local.start : local.end], hit["strand"]) == orient(
            source[global_interval.start : global_interval.end], hit["genomic_strand"]
        )
    for signal in p.get("poly_t_signals", []):
        local = checked_interval(signal, len(expected))
        assert set(expected[local.start : local.end]) == {"T"}
    interaction = p.get("interaction")
    if interaction:
        target = checked_interval(interaction["target_interval"], len(expected))
        repeat = checked_interval(interaction["repeat_interval"], len(p["array"]["repeat"]))
        left, right = interaction["hybrid_dot_bracket"].split("&")
        assert len(left) == target.length and len(right) == repeat.length
        target_paired = [target.start + i + 1 for i, value in enumerate(left) if value == "("]
        repeat_paired = [repeat.start + i + 1 for i, value in enumerate(right) if value == ")"]
        assert target_paired == interaction["target_paired_positions"]
        assert repeat_paired == interaction["repeat_paired_positions"]
        assert len(target_paired) == len(repeat_paired) > 0
    admission = p.get("short_rescue_admission")
    if admission and admission["accepted"]:
        assert p["anti_repeat"]["complementary_hypothesis"]
        assert interaction is not None and interaction["energy_kcal_mol"] < 0
        actual = [
            i
            for i in interaction["target_paired_positions"]
            if anti_local.start <= i - 1 < anti_local.end
        ]
        assert actual == admission["anchored_paired_positions"]
        assert len(actual) >= 10 and len(actual) / anti_local.length >= 0.5
        assert all(
            pair["bases"] in {"AU", "UA", "GC", "CG", "GU", "UG"}
            for pair in admission["antiparallel_pairs"]
        )
        assert admission["supporting_model_hits"] or admission["supporting_terminator_hits"]
        for h in admission["supporting_model_hits"]:
            site = checked_interval(h["interval"], len(expected))
            assert (
                h["strand"] == "+"
                and h["evalue"] <= 0.01
                and site.start >= target.start
                and site.overlap(target) < 10
            )
        for h in admission["supporting_terminator_hits"]:
            site = checked_interval(h["interval"], len(expected))
            assert h["strand"] == "+" and h["evalue"] <= 0.001 and site.start >= target.end
    # Public/internal and reverse-complement interval round trips.
    reversed_interval = Interval(len(source) - interval.end, len(source) - interval.start)
    assert (
        orient(
            reverse_complement(source)[reversed_interval.start : reversed_interval.end],
            "-" if p["strand"] == "+" else "+",
        )
        == expected
    )


def audit_run(directory):
    report = load(directory / "report.json")
    assert report["status"] == "complete" and not (directory / "failure.json").exists()
    evidence = load(directory / "manifest.json")["files"]
    for name, digest in evidence.items():
        assert sha256_file(directory / name) == digest, (directory, name)
    with (directory / "summary.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    fasta_path = directory / "predictions.fasta"
    fasta = list(read_fasta(fasta_path)) if fasta_path.stat().st_size else []
    gff = [
        line.split("\t")
        for line in (directory / "predictions.gff3").read_text().splitlines()
        if not line.startswith("#")
    ]
    predictions = []
    proposals = 0
    for record in report["records"]:
        source = list(read_fasta(Path(record["source"]["file"])))
        original = next(s for s in source if s.id == record["source"]["id"])
        assert original.sha256 == record["source"]["sha256"]
        assert sha256_file(Path(record["source"]["file"])) == record["source"]["file_sha256"]
        for p in record["predictions"] + record.get("short_rescue_proposals", []):
            prediction_audit(p, original.bases)
        predictions.extend(record["predictions"])
        proposals += len(record.get("short_rescue_proposals", []))
    assert (
        len(rows)
        == len(fasta)
        == len(gff)
        == len(predictions)
        == report["summary"]["prediction_count"]
    )
    for row, fa, gf, p in zip(rows, fasta, gff, predictions, strict=True):
        assert row["prediction_id"] == fa.id == p["export_id"]
        attrs = dict(item.split("=", 1) for item in gf[8].split(";"))
        assert unquote(attrs["ID"]) == p["export_id"]
        assert row["sequence"] == fa.bases == p["sequence"]
        assert row["strand"] == gf[6] == p["strand"]
        assert int(row["start"]) == int(gf[3]) == p["interval"]["start"]
        assert int(row["end"]) == int(gf[4]) == p["interval"]["end"]
        assert row["admission_branch"] == p["admission_branch"]
        assert row["evidence_class"] == p["evidence_class"]
    return {
        "predictions": len(predictions),
        "proposals": proposals,
        "evidence_files": len(evidence),
    }


CONTROL_MANIFEST_SHA256 = "67efd2ede9053e2638de91c2646c08c157c347c13c833a8aa4a66b51725ebb55"


def require_stage_profile(name, manifest):
    conditional = name.startswith("conditional_")
    expected = (
        name
        if conditional
        else (
            "experimental"
            if name
            in {
                "natural_experimental",
                "current_experimental_pinned_native",
                "genomic_experimental",
            }
            else "default"
        )
    )
    assert manifest["profile"] == expected, "Wrong profile for requested stage"
    pinned = name.endswith("_pinned_native") or name.startswith(("genomic_", "conditional_"))
    if pinned:
        assert manifest.get("native_dependency_profile") == "council_seed20260927_threads1", (
            "Missing mandatory pinned native metadata"
        )
    return expected


def provenance_audit(name, manifest):
    require_stage_profile(name, manifest)
    protected = load(ROUND / "protection_manifest.json")["files"]
    protocol_path = ROUND / "protocol.json"
    plan_path = ROOT / "docs/COUNCIL_IMPROVEMENT_PLAN.md"
    assert (
        sha256_file(protocol_path) == protected[str(protocol_path)] == manifest["protocol_sha256"]
    )
    assert sha256_file(plan_path) == protected[str(plan_path)] == manifest["plan_sha256"]
    protocol = load(protocol_path)
    assert protocol["plan_sha256"] == manifest["plan_sha256"]
    assert sha256_file(ROUND / "controls/manifest.json") == CONTROL_MANIFEST_SHA256
    control_manifest = load(ROUND / "controls/manifest.json")
    assert control_manifest["protocol_sha256"] == manifest["protocol_sha256"]
    assert sha256_file(ROUND / "controls/generator.py") == control_manifest["generator_sha256"]
    conditional = name.startswith("conditional_")
    genomic = name.startswith("genomic_")
    expected_profile = (
        name
        if conditional
        else (
            "experimental"
            if name
            in {
                "natural_experimental",
                "current_experimental_pinned_native",
                "genomic_experimental",
            }
            else "default"
        )
    )
    experimental = expected_profile == "experimental" or conditional
    assert manifest["profile"] == expected_profile
    config = manifest["settings"]
    for key, value in {
        "mode": "complete",
        "model_type": "II",
        "threads": 1,
        "intarna_batch_size": 128,
        "short_anti_repeat_rescue": experimental,
        "similarity": 0.7,
        "coverage": 0.6,
        "model_rescue_evalue": 1e-5,
        "mask_arrays": True,
        "model_rescue": True,
    }.items():
        assert config[key] == value, (name, key)
    assert manifest["workers"] == 2 and manifest["threads_per_locus"] == 1
    if conditional:
        assert manifest["control_manifest_sha256"] == CONTROL_MANIFEST_SHA256
    entries = load(ROOT / "benchmarks/literature/loci.json")
    expected = {e["label"]: {"fasta": e["fasta"], "sha256": e["fasta_sha256"]} for e in entries}
    if genomic:
        expected = {
            e["label"]: {"fasta": e["fasta"], "sha256": e["sha256"]}
            for e in control_manifest["genomic_controls"]
        }
        assert set(expected) == set(protocol["controls"]["genomic_labels"])
    elif name == "small_default":
        expected = {k: v for k, v in expected.items() if k in {"Ain", "Kki", "Lsp1", "Mse", "Cca1"}}
    elif name == "small_repeat":
        expected = {k: v for k, v in expected.items() if k in {"Ain", "Kki", "Cca1"}}
    identities = [e["label"] for e in manifest["inputs"]]
    assert len(identities) == len(set(identities)) == len(expected)
    assert set(identities) == set(expected)
    for entry in manifest["inputs"]:
        wanted = expected[entry["label"]]
        assert entry["fasta"] == wanted["fasta"] and entry["fasta_sha256"] == wanted["sha256"]
        assert sha256_file(ROOT / entry["fasta"]) == wanted["sha256"]
    assert sha256_file(ROOT / "upstream-lock.json") == manifest["upstream_lock_sha256"]
    for path, digest in manifest["model_sha256"].items():
        assert sha256_file(Path(path)) == digest == protected[str(Path(path).resolve())]
    for data in manifest["native_sha256"].values():
        assert sha256_file(Path(data["path"])) == data["sha256"]
    # Source revisions are recorded per run. Every path must be represented;
    # archival source snapshots permit intentional later metadata-only fixes.
    assert set(
        str(path.relative_to(ROOT)) for path in (ROOT / "src/crisprtracrrna").glob("*.py")
    ) <= set(manifest["source_sha256"])
    require_pinned = name.endswith("_pinned_native") or name.startswith(
        ("genomic_", "conditional_")
    )
    if require_pinned:
        assert "native_dependency_profile" in manifest, "Missing mandatory pinned native metadata"
    if "native_dependency_profile" in manifest:
        assert manifest["native_dependency_profile"] == "council_seed20260927_threads1"
        amendment_path = ROUND / "protocol_amendment_v2.json"
        assert (
            sha256_file(amendment_path)
            == "9cb073a515af142d77a9a87e43f4c194ec3b41a04114c3a804db2898f374f0d5"
        )
        assert manifest["protocol_amendment_sha256"] == sha256_file(amendment_path)
        native_path = ROUND / "native_execution_manifest.json"
        assert manifest["native_execution_manifest_sha256"] == sha256_file(native_path)
        assert (
            sha256_file(native_path)
            == "80f2e839d61dac8babccdc6ded60361d41a00db83cb65ce13e08e1c261826e14"
        )
        native = load(native_path)
        frozen = load(ROUND / "dependency_freeze_sha256.json")["files"]
        assert native["binary_sha256"] == frozen[native["binary"]]
        assert native["configuration_sha256"] == frozen[str(ROOT / "local.reproducible.toml")]
        assert (
            native["build_manifest_sha256"]
            == frozen[str(ROOT / "native/fasta36-seeded-20260927/build_manifest.json")]
        )
        assert config["fasta36"] == native["binary"]
        assert manifest["native_sha256"]["fasta36"]["sha256"] == native["binary_sha256"]
        assert manifest["configuration_sha256"] == native["configuration_sha256"]
        assert sha256_file(ROOT / "local.reproducible.toml") == native["configuration_sha256"]
    return expected


def archived_provenance_audit():
    directory = ROUND / "archived_detector_pinned_native"
    manifest = load(directory / "run_manifest.json")
    native_path = ROUND / "native_execution_manifest.json"
    assert (
        sha256_file(native_path)
        == "80f2e839d61dac8babccdc6ded60361d41a00db83cb65ce13e08e1c261826e14"
    )
    native = load(native_path)
    frozen = load(ROUND / "dependency_freeze_sha256.json")["files"]
    assert native["binary_sha256"] == frozen[native["binary"]]
    assert native["configuration_sha256"] == frozen[str(ROOT / "local.reproducible.toml")]
    assert manifest["profile"] == "archived_default_pinned_native"
    assert manifest["native_execution_manifest_sha256"] == sha256_file(native_path)
    assert (
        manifest["protocol_amendment_sha256"]
        == "9cb073a515af142d77a9a87e43f4c194ec3b41a04114c3a804db2898f374f0d5"
    )
    assert manifest["protocol_sha256"] == sha256_file(ROUND / "protocol.json")
    assert manifest["configuration_sha256"] == native["configuration_sha256"]
    assert manifest["settings"]["fasta36"] == native["binary"]
    assert manifest["settings"]["threads"] == manifest["threads_per_locus"] == 1
    assert manifest["workers"] == 2
    assert manifest["executables"]["fasta36"]["sha256"] == native["binary_sha256"]
    assert sha256_file(Path(native["binary"])) == native["binary_sha256"]
    expected_source = load(ROUND / "baseline/source_sha256.json")
    assert manifest["baseline_source_archive_sha256"] == sha256_file(ROUND / "baseline/source.zip")
    source_root = Path(manifest["source_location"])
    assert set(manifest["source_sha256"]) == {
        "src/crisprtracrrna/" + p.name for p in source_root.glob("*.py")
    }
    for path, digest in manifest["source_sha256"].items():
        assert digest == expected_source[path] == sha256_file(source_root / Path(path).name)
    data_root = Path(manifest["resource_relocation"])
    for relative, digest in load(ROOT / "src/crisprtracrrna/data/manifest.json").items():
        assert sha256_file(data_root / relative) == digest
    expected_inputs = load(ROOT / "benchmarks/literature/loci.json")
    assert len(manifest["inputs"]) == len(manifest["results"]) == 79
    assert [e["label"] for e in manifest["inputs"]] == [e["label"] for e in expected_inputs]
    for entry, expected in zip(manifest["inputs"], expected_inputs, strict=True):
        assert entry["fasta"] == expected["fasta"]
        assert (
            entry["fasta_sha256"] == expected["fasta_sha256"] == sha256_file(ROOT / entry["fasta"])
        )
        report = load(directory / entry["label"] / "report.json")
        assert report["status"] == "complete"
        report_membership(report, entry, manifest["settings"])
    return manifest


def report_membership(report, entry, settings):
    assert report["configuration"] == settings
    source_path = ROOT / entry["fasta"]
    expected = list(read_fasta(source_path))
    assert len(report["records"]) == len(expected)
    for record, source in zip(report["records"], expected, strict=True):
        info = record["source"]
        assert Path(info["file"]).resolve() == source_path.resolve()
        assert info["file_sha256"] == entry["fasta_sha256"]
        assert info["id"] == source.id and info["sha256"] == source.sha256
        assert info["length"] == len(source.bases)
    if not settings.get("short_anti_repeat_rescue", False):
        assert all(not r.get("short_rescue_proposals") for r in report["records"])
        assert all(
            p.get("admission_branch") != "short_rescue"
            for r in report["records"]
            for p in r["predictions"]
        )


def stage(name, baseline=False, experimental=False, comparison_stage=None):
    directory = ROUND / name
    manifest = load(directory / "run_manifest.json")
    provenance_audit(name, manifest)
    if comparison_stage == "archived_detector_pinned_native":
        archived_provenance_audit()
    results = []
    for item in manifest["inputs"]:
        label = item["label"]
        assert sha256_file(ROOT / item["fasta"]) == item["fasta_sha256"]
        report = load(directory / label / "report.json")
        report_membership(report, item, manifest["settings"])
        row = {"label": label, **audit_run(directory / label)}
        if baseline:
            old = (
                load(ROUND / comparison_stage / label / "report.json")
                if comparison_stage
                else load(ROUND / "baseline/reports" / f"{label}.json")
            )
            before = [p for r in old["records"] for p in r["predictions"]]
            after = [p for r in report["records"] for p in r["predictions"]]
            same = projection(before, experimental) == projection(after, experimental)
            row["baseline_projection_equal"] = same
            assert same, label
            if experimental:
                old_ids = {p["id"]: p["export_id"] for p in before}
                row["export_id_shifts"] = sum(
                    old_ids[p["id"]] != p["export_id"] for p in after if p["id"] in old_ids
                )
        results.append(row)
    result = {
        "stage": name,
        "auditor_sha256": sha256_file(Path(__file__)),
        "comparison_stage": comparison_stage or ("historical_frozen" if baseline else None),
        "status": "pass",
        "runs": len(results),
        "totals": dict(
            sum(
                (
                    Counter({k: r[k] for k in ("predictions", "proposals", "evidence_files")})
                    for r in results
                ),
                Counter(),
            )
        ),
        "results": results,
    }
    atomic_json(ROUND / f"audit_{name}.json", result)
    print(json.dumps({key: value for key, value in result.items() if key != "results"}))


def conditional(name):
    directory = ROUND / f"conditional_{name}"
    manifest = load(directory / "run_manifest.json")
    provenance_audit("conditional_" + name, manifest)
    rows = []
    controls = load(ROUND / "controls/manifest.json")
    expected_controls = {}
    for item in controls["repeat_controls"]:
        if item["panel"] == name:
            expected_controls.setdefault(item["label"], {})[item["array_id"]] = item[
                "control_repeat"
            ]
    for entry in manifest["inputs"]:
        work = directory / entry["label"]
        result = load(work / "result.json")
        assert result["status"] == "complete"
        expected_repeats = expected_controls.get(entry["label"], {})
        assert result["control_repeats"] == expected_repeats
        assert result["arrays"] == len(expected_repeats)
        assert result["short_proposals"] == len(result["proposals"])
        admitted = sum(p["short_rescue_admission"]["accepted"] for p in result["proposals"])
        assert result["short_admitted"] == admitted
        source = next(read_fasta(work / "input.fasta"))
        assert source.sha256 == result["source_sha256"]
        original = next(read_fasta(ROOT / entry["fasta"]))
        assert source.sha256 == original.sha256
        files = load(work / "manifest.json")["files"]
        for path, digest in files.items():
            assert sha256_file(work / path) == digest
        for proposal in result["proposals"]:
            prediction_audit(proposal, source.bases)
        rows.append(
            {
                "label": entry["label"],
                "proposals": len(result["proposals"]),
                "short_admitted": result["short_admitted"],
                "evidence_files": len(files),
            }
        )
    assert len(rows) == manifest["complete"] == 79
    assert sum(r["short_admitted"] for r in rows) == manifest["short_admitted"]
    atomic_json(
        ROUND / f"audit_conditional_{name}.json",
        {"status": "pass", "runs": len(rows), "rows": rows},
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage")
    parser.add_argument("--baseline", action="store_true")
    parser.add_argument("--experimental", action="store_true")
    parser.add_argument("--conditional", action="store_true")
    parser.add_argument("--comparison-stage")
    args = parser.parse_args()
    if args.conditional:
        conditional(args.stage)
    else:
        stage(args.stage, args.baseline, args.experimental, args.comparison_stage)


if __name__ == "__main__":
    main()
