"""Join every catalog row to retained model hits using sequence hashes.

Deduplicate within model families, distinguish unsupported source annotations,
and exclude Type III/VI accidentally screened by the first prefix-based runner.
Those raw scans remain intact; they do not enter the corrected denominators.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks"


def family(subtype):
    if subtype == "II" or subtype.startswith("II-") or subtype == "unknown":
        return "II"
    if subtype == "V" or subtype.startswith("V-"):
        return "V"
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", nargs="+", help="Scan paths relative to benchmarks/")
    args = parser.parse_args()
    lookup = {}
    batches = args.batch or [
        "final/sequences",
        "final/additional_sequences",
        "final/chyou_sequences",
        "final/briner_sequences",
    ]
    inputs = {}
    for batch in batches:
        for model in ("II", "V"):
            work = BENCH / batch / model
            ids = json.loads((work / "sequence_ids.json").read_text())
            if not ids:
                continue
            fasta = work / "cm-input.fasta"
            inputs[str(fasta.relative_to(ROOT))] = hashlib.sha256(fasta.read_bytes()).hexdigest()
            hits = defaultdict(list)
            for hit in json.loads((work / "hits.json").read_text()):
                hits[hit["query"]].append(hit)
            for digest, query in ids.items():
                eligible = [
                    h for h in hits[query] if h["strand"] == "+" and h["interval"]["length"] >= 20
                ]
                result = {
                    "detected": bool(eligible),
                    "best_score": max((h["score"] for h in eligible), default=None),
                    "best_evalue": min((h["evalue"] for h in eligible), default=None),
                    "max_model_span": max((h["interval"]["length"] for h in eligible), default=0),
                    "shuffle_any_hit": bool(hits[query + "_shuffle"]),
                    "scan_batch": batch,
                }
                if (model, digest) in lookup:
                    assert lookup[model, digest]["detected"] == result["detected"]
                else:
                    lookup[model, digest] = result
    rows = []
    unique = {}
    grouped = {
        name: defaultdict(Counter) for name in ("source", "subtype", "evidence", "source_subset")
    }
    for line in (BENCH / "literature/catalog.jsonl").open():
        ref = json.loads(line)
        model = family(ref["subtype"])
        key = (model, ref["sequence_sha256"])
        result = lookup.get(key) if model else None
        if model and result is None:
            raise ValueError(f"Expected scan missing for {ref['id']}")
        scope = "Type_II_reference" if model == "II" else "not_assessed_unsupported_subtype"
        if ref["subtype"] == "unknown":
            scope = "unknown_subtype_exploratory_Type_II_screen"
        elif model == "V":
            scope = (
                "Type_V_K_reference"
                if ref["subtype"].startswith("V-K")
                else "outside_V_K_model_scope_exploratory_screen"
            )
        row = {
            k: ref[k]
            for k in ("id", "source", "label", "subtype", "evidence", "length", "sequence_sha256")
        }
        row["source_subset"] = ref.get("source_subset", "")
        row.update(
            {
                "model_family_screen": model,
                "scope": scope,
                "detected": result["detected"] if result else None,
                "best_score": result["best_score"] if result else None,
                "best_evalue": result["best_evalue"] if result else None,
                "model_span_fraction": result["max_model_span"] / ref["length"] if result else None,
                "scan_batch": result["scan_batch"] if result else "",
            }
        )
        rows.append(row)
        if result:
            unique.setdefault(key, result)
        for group in grouped:
            counts = grouped[group][ref.get(group) or "not_applicable"]
            counts["records"] += 1
            counts["screened"] += result is not None
            counts["detected"] += bool(result and result["detected"])
            counts["not_assessed"] += result is None
    summary = {
        "catalog_records": len(rows),
        "screened_records": sum(r["detected"] is not None for r in rows),
        "not_assessed_records": sum(r["detected"] is None for r in rows),
        "by_model_pool": {
            model: {
                "unique_sequences": sum(k[0] == model for k in unique),
                "recognized_unique_sequences": sum(
                    v["detected"] for k, v in unique.items() if k[0] == model
                ),
                "shuffles_with_any_hit": sum(
                    v["shuffle_any_hit"] for k, v in unique.items() if k[0] == model
                ),
            }
            for model in ("II", "V")
        },
        "scan_input_sha256": inputs,
        "scan_batches": batches,
        **{f"by_{group}": counts for group, counts in grouped.items()},
        "recognition_rule": (
            "Reported hit on + strand, >=20 nt model span; not full boundary recovery"
        ),
        "controls": (
            "One mononucleotide shuffle per unique sequence, seed20260926; any strand hit. "
            "Not biological negatives. First batch retained for duplicate hashes. "
            "RNG stream depends on pool membership/order: exact historical control counts "
            "require the retained cm-input.fasta files identified by the input hashes."
        ),
        "history": (
            "Canonical exact major-type grouping. Initial raw summaries used startswith II/V, "
            "which also included III/VI. Retained hits are rejoined here; "
            "III/VI do not enter corrected denominators."
        ),
        "interpretation": (
            "Recognition of predicted sequences measures model agreement; "
            "no universal biological sensitivity or precision estimate."
        ),
    }
    out = BENCH / "results"
    with (out / "sequence_detection.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (out / "sequence_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "catalog_records",
                    "screened_records",
                    "not_assessed_records",
                    "by_model_pool",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
