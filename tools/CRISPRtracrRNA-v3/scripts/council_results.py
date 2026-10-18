"""Summarize the frozen round without changing detector thresholds or historical data."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from statistics import median

from summarize_literature_benchmark import localization_match

from crisprtracrrna.adapters import sha256_file
from crisprtracrrna.pipeline import atomic_json

ROOT = Path(__file__).resolve().parents[1]
ROUND = ROOT / "benchmarks/council_round2"


def load(path):
    return json.loads(path.read_text())


def evidence_rows(stage):
    rows = []
    for entry in load(ROUND / stage / "run_manifest.json")["inputs"]:
        report = load(ROUND / stage / entry["label"] / "report.json")
        for record in report["records"]:
            for proposal in record.get("short_rescue_proposals", []):
                admission = proposal["short_rescue_admission"]
                rows.append(
                    {
                        "label": entry["label"],
                        "contig": record["source"]["id"],
                        "candidate_id": proposal["id"],
                        "interval": proposal["interval"],
                        "strand": proposal["strand"],
                        "array_id": proposal["array"]["id"],
                        "evidence_class": proposal["evidence_class"],
                        "admission": admission,
                        "prediction": proposal,
                    }
                )
    return rows


def localization_rows(directory, frozen=False):
    rows = {}
    for entry in load(ROOT / "benchmarks/literature/loci.json"):
        path = (
            directory / (entry["label"] + ".json")
            if frozen
            else directory / entry["label"] / "report.json"
        )
        report = load(path)
        ranked = [(r["source"]["id"], p) for r in report["records"] for p in r["predictions"]]
        ranked.sort(key=lambda cp: (-cp[1]["ranking"]["score"], cp[1]["id"]))
        matches = [
            (i, p)
            for i, (contig, p) in enumerate(ranked, 1)
            if entry["status"] == "exact_unique"
            and localization_match(p, entry["matches"][0], contig)
        ]
        rows[entry["label"]] = {
            "detected": bool(matches) if entry["status"] == "exact_unique" else None,
            "best_rank": matches[0][0] if matches else None,
            "candidate_count": len(ranked),
        }
    return rows


def retention(historical, stages):
    required = {label for label, row in historical.items() if row["detected"] is True}
    return {
        name: {
            "required_labels": sorted(required),
            "retained": sorted(label for label in required if stage[label]["detected"] is True),
            "lost": sorted(label for label in required if stage[label]["detected"] is not True),
        }
        for name, stage in stages.items()
    }


def main():
    default = load(ROUND / "current_default_pinned_native/metrics.json")
    experimental = load(ROUND / "current_experimental_pinned_native/metrics.json")
    historical = localization_rows(ROUND / "baseline/reports", frozen=True)
    archived = localization_rows(ROUND / "archived_detector_pinned_native")
    before = {row["label"]: row for row in default["rows"]}
    after = {row["label"]: row for row in experimental["rows"]}
    historical_retention = retention(
        historical,
        {"archived_pinned": archived, "default_pinned": before, "experimental_pinned": after},
    )
    assert len(next(iter(historical_retention.values()))["required_labels"]) == 68
    dependency_deltas = [
        {
            "label": label,
            "historical": historical[label],
            "archived_pinned": archived[label],
            "gained": archived[label]["detected"] is True
            and historical[label]["detected"] is not True,
            "lost": historical[label]["detected"] is True
            and archived[label]["detected"] is not True,
        }
        for label in historical
    ]
    atomic_json(ROUND / "historical_dependency_deltas.json", dependency_deltas)
    paired = [
        {
            "label": label,
            "default": before[label],
            "experimental": after[label],
            "added_candidates": after[label]["candidate_count"] - before[label]["candidate_count"],
            "gained": after[label]["detected"] is True and before[label]["detected"] is not True,
            "lost": before[label]["detected"] is True and after[label]["detected"] is not True,
        }
        for label in before
    ]
    atomic_json(ROUND / "paired_natural_results.json", paired)
    proposals = evidence_rows("current_experimental_pinned_native")
    truth = {e["label"]: e for e in load(ROOT / "benchmarks/literature/loci.json")}
    for row in proposals:
        ref = truth[row["label"]]
        row["reference_localized"] = (
            bool(localization_match(row["prediction"], ref["matches"][0], row["contig"]))
            if ref["status"] == "exact_unique"
            else None
        )
    atomic_json(ROUND / "natural_short_proposals.json", proposals)
    genomic_before = load(ROUND / "genomic_default/metrics.json")
    genomic_after = load(ROUND / "genomic_experimental/metrics.json")
    gb = {r["label"]: r for r in genomic_before["rows"]}
    ga = {r["label"]: r for r in genomic_after["rows"]}
    gp = evidence_rows("genomic_experimental")
    controls = {
        "baseline_total": genomic_before["candidate_count"],
        "experimental_total": genomic_after["candidate_count"],
        "ratio": (
            genomic_after["candidate_count"] / genomic_before["candidate_count"]
            if genomic_before["candidate_count"]
            else None
        ),
        "burden_increased_inputs": sum(
            ga[label]["candidate_count"] > gb[label]["candidate_count"] for label in gb
        ),
        "new_corroborated_short_admissions": [r for r in gp if r["admission"]["accepted"]],
        "paired": [
            {"label": label, "default": gb[label], "experimental": ga[label]} for label in gb
        ],
    }
    atomic_json(ROUND / "paired_genomic_controls.json", controls)
    conditional = {}
    for name in ("dinucleotide_1", "dinucleotide_2", "derangement"):
        manifest = load(ROUND / f"conditional_{name}/run_manifest.json")
        assert manifest["complete"] == 79
        rows = []
        for entry in manifest["inputs"]:
            result = load(ROUND / f"conditional_{name}" / entry["label"] / "result.json")
            admitted = [p for p in result["proposals"] if p["short_rescue_admission"]["accepted"]]
            by_array = dict.fromkeys(result["control_repeats"], 0)
            for proposal in admitted:
                by_array[proposal["array"]["id"]] += 1
            rows.append(
                {
                    "label": entry["label"],
                    "arrays": result["arrays"],
                    "bases": result["bases"],
                    "proposals": result["short_proposals"],
                    "admitted": len(admitted),
                    "admitted_per_mb": len(admitted) * 1e6 / result["bases"],
                    "admitted_per_array": by_array,
                    "rejections": dict(
                        Counter(
                            reason
                            for p in result["proposals"]
                            for reason in p["short_rescue_admission"]["rejection_reasons"]
                        )
                    ),
                    "admissions": admitted,
                }
            )
        counts = [r["admitted"] for r in rows]
        conditional[name] = {
            "complete": len(rows),
            "proposals": sum(r["proposals"] for r in rows),
            "admitted": sum(counts),
            "median": median(counts),
            "max": max(counts),
            "rows": rows,
            "batch_wall_seconds": manifest["batch_wall_seconds"],
        }
    atomic_json(ROUND / "conditional_controls_summary.json", conditional)
    frozen = load(ROUND / "controls/manifest.json")
    shuffle_stats = {}
    for name in ("dinucleotide_1", "dinucleotide_2", "derangement"):
        members = [r for r in frozen["repeat_controls"] if r["panel"] == name]
        shuffle_stats[name] = {
            "arrays": len(members),
            "unchanged": sum(r["unchanged"] for r in members),
            "distinct_control_repeats": len({r["control_repeat"] for r in members}),
            "same_as_other_panel": sum(
                any(
                    q["identity"] == r["identity"]
                    and q["panel"] != name
                    and q["control_repeat"] == r["control_repeat"]
                    for q in frozen["repeat_controls"]
                )
                for r in members
            ),
        }
    aggregate_evidence = {}
    for name, metrics in (("default", default), ("experimental", experimental)):
        aggregate_evidence[name] = {
            key: dict(sum((Counter(row.get(key, {})) for row in metrics["rows"]), Counter()))
            for key in ("origins", "admission_branches", "evidence_classes", "rejections")
        }
    geometry = []
    for row in experimental["rows"]:
        best = row.get("highest_ranked_match")
        if best is None:
            continue
        strand = truth[row["label"]]["matches"][0]["strand"]
        start_error, end_error = best["genomic_start_error"], best["genomic_end_error"]
        geometry.append(
            {
                "label": row["label"],
                "strand": strand,
                **best,
                "five_prime_error_along_transcription": start_error
                if strand == "+"
                else -end_error,
                "three_prime_error_along_transcription": end_error
                if strand == "+"
                else -start_error,
                "best_attainable_iou": row["best_iou"],
                "note": (
                    "Positive endpoint errors mean downstream in transcript orientation; "
                    "windows are not inferred endpoints."
                ),
            }
        )
    atomic_json(ROUND / "experimental_geometry.json", geometry)
    performance = load(ROUND / "performance/results.json")["results"]
    added_max = max(row["added_candidates"] for row in paired)
    ratio = experimental["batch_wall_seconds"] / default["batch_wall_seconds"]
    gates = {
        "original_historical_default_exact_projection": False,
        "natural_completion": default["complete"] == experimental["complete"] == 79,
        "dependency_controlled_default_exact_projection": load(
            ROUND / "audit_current_default_pinned_native.json"
        )["status"]
        == "pass",
        "experimental_preserves_standard_projection": load(
            ROUND / "audit_current_experimental_pinned_native.json"
        )["status"]
        == "pass",
        "interaction_map_equality": all(
            r["whole_map_equal"] and r["baseline_equal"] for r in performance
        ),
        "mse_grouped_faster": performance[0]["wall_seconds"]["128"]
        < performance[0]["wall_seconds"]["1"],
        "experimental_localized_min69": experimental["localized"] >= 69,
        "all_historical68_retained_in_all_pinned_stages": all(
            not r["lost"] for r in historical_retention.values()
        ),
        "all_pinned_default_detections_retained": not any(r["lost"] for r in paired),
        "top1_min53": experimental["top_k"]["1"] >= 53,
        "top5_min63": experimental["top_k"]["5"] >= 63,
        "coverage90_min64": experimental["coverage90"] >= 64,
        "candidate_burden_max5092": experimental["candidate_count"] <= 5092,
        "candidate_p95_max288": experimental["candidate_p95"] <= 288,
        "candidate_max928": experimental["candidate_max"] <= 928,
        "added_per_locus_max100": added_max <= 100,
        "experimental_wall_ratio_max2": ratio <= 2,
        "conditional_controls_complete": all(v["complete"] == 79 for v in conditional.values()),
        "conditional_median_max1": all(v["median"] <= 1 for v in conditional.values()),
        "conditional_max10": all(v["max"] <= 10 for v in conditional.values()),
    }
    result = {
        "historical_retention": historical_retention,
        "dependency_driven_gains": [r["label"] for r in dependency_deltas if r["gained"]],
        "dependency_driven_losses": [r["label"] for r in dependency_deltas if r["lost"]],
        "default": {k: v for k, v in default.items() if k != "rows"},
        "experimental": {k: v for k, v in experimental.items() if k != "rows"},
        "gained_labels": [r["label"] for r in paired if r["gained"]],
        "lost_labels": [r["label"] for r in paired if r["lost"]],
        "added_candidates_max_per_locus": added_max,
        "experimental_wall_ratio": ratio,
        "natural_rejections": dict(
            Counter(reason for r in proposals for reason in r["admission"]["rejection_reasons"])
        ),
        "aggregate_evidence": aggregate_evidence,
        "performance": performance,
        "shuffle_stats": shuffle_stats,
        "conditional": {
            k: {a: b for a, b in v.items() if a != "rows"} for k, v in conditional.items()
        },
        "genomic_control_review_required": controls["burden_increased_inputs"] > len(gb) / 2
        or (controls["ratio"] is not None and controls["ratio"] > 1.5)
        or bool(controls["new_corroborated_short_admissions"]),
        "genomic_control_new_short_admissions": len(controls["new_corroborated_short_admissions"]),
        "genomic_control_burden_ratio": controls["ratio"],
        "gates": gates,
        "protocol_sha256": sha256_file(ROUND / "protocol.json"),
        "control_manifest_sha256": sha256_file(ROUND / "controls/manifest.json"),
        "historical_equivalence_note": (
            "Original stock-native comparison FAILED at Lsp2, 16 versus 14 "
            "candidates; fixed-dependency comparisons are separately labeled and do not erase it."
        ),
        "limitations": [
            "Inspected development/regression loci; no independent generalization estimate",
            "Unmatched predictions are unverified, not established false positives",
            "Search windows are not native transcript endpoints",
            "Genomic controls disrupt Cas context; conditional controls retain inferred arrays",
            "No default promotion; no calibrated confidence, precision, specificity or FDR claim",
        ],
    }
    atomic_json(ROUND / "final_metrics.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
