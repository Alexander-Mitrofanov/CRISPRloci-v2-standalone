"""Reproducible sequence and natural-locus benchmarks using the real v3 tools.

Sequence recognition is not genome sensitivity. Synthetic shuffled controls are
not validated biological negatives. Genome truth requires unique exact mapping.
"""

from __future__ import annotations

import argparse
import json
import random
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from summarize_literature_benchmark import localization_match
from summarize_sequence_benchmark import family

from crisprtracrrna.config import load_settings
from crisprtracrrna.models import Sequence
from crisprtracrrna.pipeline import run
from crisprtracrrna.process import Runner
from crisprtracrrna.search import scan_models

ROOT = Path(__file__).resolve().parents[1]


def sequence_benchmark(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    rows = [json.loads(line) for line in (ROOT / "benchmarks/literature/catalog.jsonl").open()]
    if args.sources:
        rows = [r for r in rows if r["source"] in args.sources]
    summary = {}
    results = []
    for model_type in ("II", "V"):
        members = [r for r in rows if family(r["subtype"]) == model_type]
        unique = {r["sequence_sha256"]: r["sequence"] for r in members}
        records = []
        seq_to_id = {}
        rng = random.Random(20260926)
        for i, (digest, sequence) in enumerate(sorted(unique.items())):
            identity = f"seq{i:06d}"
            seq_to_id[digest] = identity
            records.append(Sequence(identity, identity, sequence))
            bases = list(sequence)
            rng.shuffle(bases)
            records.append(Sequence(identity + "_shuffle", identity + "_shuffle", "".join(bases)))
        work = output / model_type
        work.mkdir()
        settings = load_settings(
            args.config, {"mode": "model", "model_type": model_type, "threads": 2}
        )
        runner = Runner(work / "logs", settings.timeout, settings.environment)
        started = time.monotonic()
        hits = scan_models(records, work, settings, runner)
        by_query = defaultdict(list)
        for h in hits:
            by_query[h.query].append(h)
        (work / "hits.json").write_text(
            json.dumps([{"query": h.query, **h.public()} for h in hits], indent=2)
        )
        (work / "sequence_ids.json").write_text(json.dumps(seq_to_id, indent=2))
        for row in members:
            identity = seq_to_id[row["sequence_sha256"]]
            eligible = [
                h for h in by_query[identity] if h.strand == "+" and h.interval.length >= 20
            ]
            results.append(
                {
                    "id": row["id"],
                    "source": row["source"],
                    "subtype": row["subtype"],
                    "label": row["label"],
                    "sequence_sha256": row["sequence_sha256"],
                    "detected": bool(eligible),
                    "best_score": max((h.score for h in eligible), default=None),
                    "max_reference_coverage": max(
                        (h.interval.length / row["length"] for h in eligible), default=0
                    ),
                    "wrong_strand_only": bool(by_query[identity])
                    and not any(h.strand == "+" for h in by_query[identity]),
                }
            )
        summary[model_type] = {
            "unique_sequences": len(unique),
            "elapsed_seconds": time.monotonic() - started,
            "recognized_unique_sequences": sum(
                any(h.strand == "+" and h.interval.length >= 20 for h in by_query[q])
                for q in seq_to_id.values()
            ),
            "shuffled_sequences_with_any_hit": sum(
                bool(by_query[q + "_shuffle"]) for q in seq_to_id.values()
            ),
            "settings": settings.public(),
        }
    grouped = defaultdict(Counter)
    for r in results:
        grouped[r["source"]]["total"] += 1
        grouped[r["source"]]["detected"] += r["detected"]
    summary["by_source"] = dict(grouped)
    summary["rule"] = (
        "At least one reported model hit on the correct strand spanning >=20 nt. "
        "Not full transcript boundary recovery."
    )
    (output / "sequence_results.json").write_text(json.dumps(results, indent=2) + "\n")
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


def genome_benchmark(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    entries = json.loads((ROOT / "benchmarks/literature/loci.json").read_text())
    settings = load_settings(args.config, {"mode": args.mode, "threads": 1})

    def analyze(entry):
        target = output / entry["label"]
        started = time.monotonic()
        try:
            report = run(ROOT / entry["fasta"], target, settings)
            predictions = [
                (rec["source"]["id"], p) for rec in report["records"] for p in rec["predictions"]
            ]
            match = entry["matches"][0] if entry["status"] == "exact_unique" else None
            eligible = []
            if match:
                for contig, p in predictions:
                    overlap = localization_match(p, match, contig)
                    if overlap:
                        eligible.append((p, overlap))
            ranked = sorted(
                predictions,
                key=lambda cp: -cp[1].get("ranking", {}).get("score", cp[1].get("best_score", 0)),
            )
            ranks = [i for i, (_, p) in enumerate(ranked, 1) if any(p is q for q, _ in eligible)]
            return {
                "id": entry["id"],
                "label": entry["label"],
                "mapping_status": entry["status"],
                "detected": bool(eligible) if match else None,
                "prediction_count": len(predictions),
                "accepted_array_count": sum(len(rec["arrays"]) for rec in report["records"]),
                "best_rank": min(ranks, default=None),
                "max_reference_coverage": max(
                    (o / (match["end"] - match["start"] + 1) for _, o in eligible), default=0
                )
                if match
                else None,
                "elapsed_seconds": time.monotonic() - started,
                "status": "success",
            }
        except Exception as exc:
            return {
                "id": entry["id"],
                "label": entry["label"],
                "mapping_status": entry["status"],
                "status": "error",
                "error": str(exc),
                "elapsed_seconds": time.monotonic() - started,
            }

    results = []
    with ThreadPoolExecutor(args.workers) as pool, (output / "results.jsonl").open("w") as out:
        for r in pool.map(analyze, entries):
            results.append(r)
            out.write(json.dumps(r) + "\n")
            out.flush()
            print(json.dumps(r), flush=True)
    (output / "summary.json").write_text(
        json.dumps(
            {
                "total": len(results),
                "success": sum(r["status"] == "success" for r in results),
                "mapped": sum(r["mapping_status"] == "exact_unique" for r in results),
                "detected": sum(r.get("detected") is True for r in results),
                "localization_rule": (
                    "Same contig and strand, >=20 nt overlap and >=50% of the shorter interval"
                ),
                "results": results,
            },
            indent=2,
        )
        + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=["sequence", "genome"])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=ROOT / "local.toml")
    parser.add_argument("--mode", choices=["complete", "model"], default="complete")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--sources", nargs="+", help="Optional source subset for sequence scans")
    args = parser.parse_args()
    (sequence_benchmark if args.kind == "sequence" else genome_benchmark)(args)


if __name__ == "__main__":
    main()
