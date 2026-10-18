"""Apply one explicit locus metric to every retained benchmark stage.

Counts candidate localization separately from transcript coverage and execution
errors. The positive-only reference set cannot estimate genome-wide precision.
"""

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks"


def localization_match(prediction, truth, contig):
    if contig != truth["contig"] or prediction["strand"] != truth["strand"]:
        return 0
    p = prediction["interval"]
    overlap = max(0, min(p["end"], truth["end"]) - max(p["start"], truth["start"]) + 1)
    shorter = min(p["end"] - p["start"] + 1, truth["end"] - truth["start"] + 1)
    return overlap if overlap >= 20 and overlap >= 0.5 * shorter else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", nargs="+", help="Stage paths relative to benchmarks/")
    args = parser.parse_args()
    entries = json.loads((BENCH / "literature/loci.json").read_text())
    stages = args.stage or [
        "baseline/complete",
        "baseline/genomic_models",
        "after_path_fix/complete",
        "after_rescue/complete",
        "final/complete",
        "final_fixed/complete",
    ]
    all_rows = []
    summary = {}
    for stage in stages:
        rows = []
        if not (BENCH / stage / "summary.json").exists():
            if args.stage:
                raise FileNotFoundError(BENCH / stage / "summary.json")
            continue
        for e in entries:
            row = {
                "stage": stage,
                "label": e["label"],
                "reference_id": e["id"],
                "mapping_status": e["status"],
                "status": "error",
                "detected": False if e["status"] == "exact_unique" else None,
                "best_rank": None,
                "reference_coverage": 0.0,
                "candidate_count": 0,
                "model_rescue_count": 0,
                "error": "",
                "best_origin": "",
            }
            report = BENCH / stage / e["label"] / "report.json"
            if report.exists():
                d = json.loads(report.read_text())
                row["status"] = "success"
                ranked = sorted(
                    [(rec["source"]["id"], p) for rec in d["records"] for p in rec["predictions"]],
                    key=lambda cp: (
                        -cp[1].get("ranking", {}).get("score", cp[1].get("best_score", 0))
                    ),
                )
                row["candidate_count"] = len(ranked)
                row["model_rescue_count"] = sum(
                    p.get("candidate_origin") == "covariance_model" for _, p in ranked
                )
                if e["status"] == "exact_unique":
                    truth = e["matches"][0]
                    for rank, (contig, p) in enumerate(ranked, 1):
                        overlap = localization_match(p, truth, contig)
                        if overlap:
                            row["detected"] = True
                            if row["best_rank"] is None:
                                row["best_rank"] = rank
                                row["best_origin"] = p.get("candidate_origin", "legacy_candidate")
                            row["reference_coverage"] = max(
                                row["reference_coverage"],
                                overlap / (truth["end"] - truth["start"] + 1),
                            )
            else:
                failure = report.with_name("failure.json")
                row["error"] = (
                    json.loads(failure.read_text())["error"]
                    if failure.exists()
                    else "No completed report"
                )
            rows.append(row)
        mapped = [r for r in rows if r["mapping_status"] == "exact_unique"]
        summary[stage] = {
            "attempted": len(rows),
            "successful": sum(r["status"] == "success" for r in rows),
            "mapped_reference_loci": len(mapped),
            "detected": sum(r["detected"] for r in mapped),
            "recall_including_execution_errors": sum(r["detected"] for r in mapped) / len(mapped),
            "top1": sum(r["best_rank"] == 1 for r in mapped),
            "top5": sum(r["best_rank"] is not None and r["best_rank"] <= 5 for r in mapped),
            "reference_coverage_at_least_90_percent": sum(
                r["reference_coverage"] >= 0.9 for r in mapped
            ),
            "total_candidates": sum(r["candidate_count"] for r in rows),
            "missed_or_error": [r["label"] for r in mapped if not r["detected"]],
        }
        all_rows.extend(rows)
    out = BENCH / "results"
    out.mkdir(exist_ok=True)
    (out / "genome_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "genome_results.json").write_text(json.dumps(all_rows, indent=2) + "\n")
    with (out / "genome_results.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(all_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)
    (out / "metric.json").write_text(
        json.dumps(
            {
                "truth": (
                    "Unique exact match of the published oriented sequence in the supplied locus, "
                    "both strands searched"
                ),
                "localization": (
                    "Same contig and strand, >=20 nt overlap and >=50% of the shorter interval"
                ),
                "coverage": (
                    "Maximum fraction of reference bases overlapped by a localized prediction"
                ),
                "errors": (
                    "Counted as non-detections in the fixed mapped denominator; "
                    "also reported separately"
                ),
                "precision": (
                    "Not estimated: reference loci are not exhaustively annotated "
                    "for alternative tracrRNAs"
                ),
                "history": (
                    "Initial exploratory run counters used a looser overlap rule. "
                    "This final metric uniformly rescores retained predictions from every stage."
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
