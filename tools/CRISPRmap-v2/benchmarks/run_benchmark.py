#!/usr/bin/env python3
"""Compare with original outputs and audit historical annotation recovery.

Run from an installed checkout after scripts/fetch_legacy.py. This program
never downloads files, executes the original server, or estimates its runtime.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from crisprmap.io import read_fasta  # noqa: E402
from crisprmap.reference import build_reference, load_reference  # noqa: E402

spec = importlib.util.spec_from_file_location("fetch_legacy", ROOT / "scripts/fetch_legacy.py")
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def group_key(sequence):
    return min(sequence, reverse_complement(sequence))


def run_search(records, reference, max_distance, top_k=20):
    from crisprmap.core import search
    start = time.perf_counter()
    report = search(records, reference, max_distance=max_distance, top_k=top_k,
                    max_comparisons=20_000_000, deadline_seconds=600)
    return report["results"], time.perf_counter() - start


def label_comparison(rows, truth, fields=("family", "motif")):
    result = {}
    for field in fields:
        labelled = [(row, expected[field]) for row, expected in zip(rows, truth) if expected[field] is not None]
        assigned = [(row, expected) for row, expected in labelled if row["annotations"][field]["value"] is not None]
        result[field] = {
            "known_original_labels": len(labelled),
            "assigned_known_labels": len(assigned),
            "correct_assigned_known_labels": sum(row["annotations"][field]["value"] == expected for row, expected in assigned),
            "known_label_assignment_coverage": len(assigned) / len(labelled) if labelled else None,
            "agreement_among_assigned_known_labels": sum(row["annotations"][field]["value"] == expected for row, expected in assigned) / len(assigned) if assigned else None,
            "missing_original_labels": sum(expected[field] is None for expected in truth),
            "raw_value_matches_including_null": sum(row["annotations"][field]["value"] == expected[field] for row, expected in zip(rows, truth)),
        }
    return result


def original_comparison(data_dir, reference, max_distance):
    # Reverify raw files each run; a modified cached baseline must not be accepted.
    for name, asset in legacy.ASSETS.items():
        legacy.checked_bytes((data_dir / name).read_bytes(), asset)
    records, input_sha = read_fasta(data_dir / "original-input.fa")
    truth = legacy.parse_original_results((data_dir / "original-results.tsv").read_text())
    if len(records) != len(truth) or any(group_key(row["sequence"]) != group_key(expected["sequence"]) for row, expected in zip(records, truth)):
        raise ValueError("Original input/output pairing failed sequence or reverse-complement verification")
    rows, seconds = run_search(records, reference, max_distance)
    detail = []
    for record, row, expected in zip(records, rows, truth):
        detail.append({
            "input_id": record["id"], "original_row": expected["ordinal"],
            "original_database_member": expected["original_reference_id"] is not None,
            "v2_status": row["status"], "v2_best_distance": row["best_distance"],
            "original_family": expected["family"], "v2_family": row["annotations"]["family"]["value"],
            "v2_family_candidates": row["annotations"]["family"]["values"],
            "original_motif": expected["motif"], "v2_motif": row["annotations"]["motif"]["value"],
            "v2_family_evidence_status": row["annotations"]["family"]["status"],
            "v2_motif_evidence_status": row["annotations"]["motif"]["status"],
            "v2_reference_motif_values": row["annotations"]["motif"]["values"],
        })
    return {
        "kind": "published_original_output_replay", "original_version": "v1.0-2013",
        "original_example_job_id": "8888684", "query_count": len(records),
        "original_executable_rerun": False, "original_runtime_seconds": None,
        "runtime_speedup": None, "v2_wall_seconds_including_reference_validation": seconds,
        "input_sha256": input_sha, "output_sha256": legacy.ASSETS["original-results.tsv"]["sha256"],
        "status_counts": dict(Counter(row["status"] for row in rows)),
        "metrics": label_comparison(rows, truth), "rows": detail,
    }


def self_recovery(reference, max_distance):
    entries = reference["entries"]
    records = [{"id": f"self_{i}", "sequence": row["sequence"]} for i, row in enumerate(entries)]
    records += [{"id": f"rc_{i}", "sequence": reverse_complement(row["sequence"])} for i, row in enumerate(entries)]
    rows, seconds = run_search(records, reference, max_distance)
    counts = Counter(row["status"] for row in rows)
    annotation_states = {field: dict(Counter(row["annotations"][field]["status"] for row in rows)) for field in ("family", "motif")}
    # Family/motif value and conflict states must be invariant to query orientation.
    orientation_invariant = sum(a["annotations"] == b["annotations"] for a, b in zip(rows[:len(entries)], rows[len(entries):]))
    return {"kind": "reference_self_recovery_not_independent_validation", "query_count": len(records),
            "reference_record_count": len(entries), "status_counts": dict(counts),
            "exact_recovery_count": counts["exact"], "orientation_annotation_invariant_pairs": orientation_invariant,
            "annotation_status_counts": annotation_states,
            "v2_wall_seconds_including_reference_validation": seconds}


def make_holdout(reference, requested, fold=None):
    groups = defaultdict(list)
    for row in reference["entries"]:
        groups[group_key(row["sequence"])].append(row)
    keys = sorted(groups, key=lambda key: (hashlib.sha256(("crisprmap-v2-heldout-2026-v1:" + key).encode()).hexdigest(), key))
    selected = ({key for key in keys if int(hashlib.sha256(("crisprmap-v2-heldout-2026-v1:" + key).encode()).hexdigest(), 16) % 5 == fold}
                if fold is not None else set(keys[:min(requested, len(keys) - 1)]))
    train = [row for row in reference["entries"] if group_key(row["sequence"]) not in selected]
    records = [{"id": f"heldout_{i:04d}", "sequence": key} for i, key in enumerate(keys) if key in selected]
    truth = [groups[record["sequence"]] for record in records]
    metadata = {row["id"]: {field: row[field] for field in ("family", "motif", "accession")} for row in train}
    provenance = dict(reference["provenance"])
    provenance["benchmark_split"] = "SHA256-ranked canonical sequence groups, duplicate and reverse-complement occurrences removed together"
    train_ref = build_reference(train, metadata, provenance)
    train_groups = {group_key(row["sequence"]) for row in train}
    if selected & train_groups:
        raise AssertionError("Held-out leakage: an identical or reverse-complement query remains in reference")
    return records, truth, train_ref, len(groups), selected


def heldout_comparison(reference, requested, max_distance):
    records, truth_groups, rows, selected = [], [], [], set()
    seconds, training_records, training_hashes = 0.0, [], []
    for fold in (range(5) if requested == 0 else [None]):
        fold_records, fold_truth, train, group_count, fold_selected = make_holdout(reference, requested, fold)
        fold_rows, elapsed = run_search(fold_records, train, max_distance)
        records.extend(fold_records)
        truth_groups.extend(fold_truth)
        rows.extend(fold_rows)
        selected.update(fold_selected)
        training_records.append(len(train["entries"]))
        training_hashes.append(train["content_sha256"])
        seconds += elapsed
    summary = Counter()
    per_family = defaultdict(Counter)
    detail = []
    for record, row, group in zip(records, rows, truth_groups):
        labels = {entry["family"] for entry in group}
        known = sorted(label for label in labels if label is not None)
        complete_single_label = len(labels) == 1 and None not in labels
        candidates = row["annotations"]["family"]["values"]
        if complete_single_label:
            label = known[0]
            summary["unambiguous_label_queries"] += 1
            per_family[label]["queries"] += 1
            if candidates:
                summary["queries_with_candidates"] += 1
                per_family[label]["with_candidates"] += 1
            if label in candidates:
                summary["true_label_in_candidates"] += 1
                per_family[label]["true_label_in_candidates"] += 1
            if candidates == [label] and row["annotations"]["family"]["status"] == "candidate":
                summary["sole_candidate_agrees"] += 1
                per_family[label]["sole_candidate_agrees"] += 1
            if len(candidates) == 1 and row["annotations"]["family"]["status"] == "candidate":
                summary["sole_candidate_queries"] += 1
                per_family[label]["sole_candidate_queries"] += 1
        elif not known:
            summary["missing_historical_label_queries"] += 1
        else:
            summary["conflicting_or_partial_historical_label_queries"] += 1
        detail.append({"query_id": record["id"], "sequence_sha256": hashlib.sha256(record["sequence"].encode()).hexdigest(),
                       "historical_family_values": known, "historical_label_complete_unanimous": complete_single_label,
                       "v2_status": row["status"], "best_distance": row["best_distance"], "family_candidates": candidates,
                       "family_evidence_status": row["annotations"]["family"]["status"]})
    denominator = summary["unambiguous_label_queries"]
    family_metrics = {}
    for family, counts in sorted(per_family.items()):
        family_metrics[family] = {
            **dict(counts), "candidate_coverage": counts["with_candidates"] / counts["queries"],
            "true_label_candidate_recall": counts["true_label_in_candidates"] / counts["queries"],
            "agreement_among_sole_candidates": counts["sole_candidate_agrees"] / counts["sole_candidate_queries"] if counts["sole_candidate_queries"] else None,
        }
    macro = {
        "families_with_unambiguous_query_labels": len(family_metrics),
        "candidate_coverage": sum(row["candidate_coverage"] for row in family_metrics.values()) / len(family_metrics) if family_metrics else None,
        "true_label_candidate_recall": sum(row["true_label_candidate_recall"] for row in family_metrics.values()) / len(family_metrics) if family_metrics else None,
    }
    contributing = [row["agreement_among_sole_candidates"] for row in family_metrics.values() if row["agreement_among_sole_candidates"] is not None]
    macro["families_contributing_sole_candidate_agreement"] = len(contributing)
    macro["agreement_among_sole_candidates"] = sum(contributing) / len(contributing) if contributing else None
    return {
        "kind": "heldout_historical_family_candidate_agreement_not_biological_accuracy",
        "unit": "canonical_sequence_or_reverse_complement_group",
        "split_rule": ("Five folds: SHA256('crisprmap-v2-heldout-2026-v1:' + canonical_sequence) modulo 5; every group evaluated once" if requested == 0 else "First requested groups ordered by SHA256('crisprmap-v2-heldout-2026-v1:' + canonical_sequence)"),
        "homology_species_or_family_separated": False, "query_count": len(records),
        "total_sequence_groups": group_count, "training_records_per_fold": training_records,
        "training_reference_sha256_per_fold": training_hashes,
        "heldout_groups_sha256": hashlib.sha256("\n".join(sorted(selected)).encode()).hexdigest(),
        "exact_or_reverse_complement_overlap_count": 0,
        "status_counts": dict(Counter(row["status"] for row in rows)), "counts": dict(summary),
        "candidate_coverage_on_unambiguous_labels": summary["queries_with_candidates"] / denominator if denominator else None,
        "true_label_candidate_recall_on_unambiguous_labels": summary["true_label_in_candidates"] / denominator if denominator else None,
        "agreement_among_sole_candidate_queries": summary["sole_candidate_agrees"] / summary["sole_candidate_queries"] if summary["sole_candidate_queries"] else None,
        "per_family": family_metrics, "macro_per_family": macro,
        "v2_wall_seconds_including_reference_validation": seconds, "rows": detail,
    }


def family_absent(reference, max_distance):
    groups = defaultdict(list)
    for entry in reference["entries"]:
        groups[group_key(entry["sequence"])].append(entry)
    families = sorted({entry["family"] for entry in reference["entries"] if entry["family"] is not None})
    rows_by_family, elapsed_total = {}, 0.0
    for family in families:
        excluded = {key for key, entries in groups.items() if any(entry["family"] == family for entry in entries)}
        query_keys = sorted(key for key in excluded if {entry["family"] for entry in groups[key]} == {family})
        if not query_keys:
            continue
        training = [entry for entry in reference["entries"] if group_key(entry["sequence"]) not in excluded]
        metadata = {entry["id"]: {field: entry[field] for field in ("family", "motif", "accession")} for entry in training}
        provenance = {**reference["provenance"], "benchmark_absent_family": family}
        train = build_reference(training, metadata, provenance)
        queries = [{"id": f"absent_{family}_{i}", "sequence": key} for i, key in enumerate(query_keys)]
        rows, elapsed = run_search(queries, train, max_distance)
        elapsed_total += elapsed
        if any(family in row["annotations"]["family"]["values"] for row in rows):
            raise AssertionError("An excluded family leaked into absent-family reference")
        cross_family = sum(bool(row["annotations"]["family"]["values"]) for row in rows)
        rows_by_family[family] = {"queries": len(queries), "removed_canonical_groups": len(excluded),
                                  "training_records": len(training), "status_counts": dict(Counter(row["status"] for row in rows)),
                                  "cross_family_candidate_queries": cross_family,
                                  "complete_unanimous_cross_family_candidate_queries": sum(row["annotations"]["family"]["status"] == "candidate" for row in rows),
                                  "query_assignment_count": sum(row["annotations"]["family"]["value"] is not None for row in rows)}
    total = sum(row["queries"] for row in rows_by_family.values())
    candidates = sum(row["cross_family_candidate_queries"] for row in rows_by_family.values())
    return {"kind": "family_absent_historical_retrieval_stress_test_not_calibrated_open_set_validation",
            "design": "Remove every canonical group carrying one family, then query all complete unanimous groups of that family; repeat across all families",
            "query_count": total, "family_count": len(rows_by_family), "cross_family_candidate_queries": candidates,
            "cross_family_candidate_rate": candidates / total if total else None,
            "query_assignment_count": sum(row["query_assignment_count"] for row in rows_by_family.values()),
            "per_family": rows_by_family, "v2_wall_seconds_including_reference_validation": elapsed_total}


def controls(reference, max_distance):
    rng = random.Random(20260930)
    source = sorted({group_key(row["sequence"]) for row in reference["entries"]})
    records = []
    for i in range(100):
        sequence = rng.choice(source)
        records.append({"id": f"random_{i:03d}", "sequence": "".join(rng.choice("ACGT") for _ in sequence)})
        shuffled = list(sequence)
        rng.shuffle(shuffled)
        records.append({"id": f"shuffled_{i:03d}", "sequence": "".join(shuffled)})
    for base in "ACGT":
        records.append({"id": "homopolymer_" + base, "sequence": base * 30})
    rows, seconds = run_search(records, reference, max_distance)
    by_kind = defaultdict(Counter)
    for record, row in zip(records, rows):
        kind = record["id"].split("_")[0]
        by_kind[kind][row["status"]] += 1
        by_kind[kind]["query_count"] += 1
        by_kind[kind]["with_family_candidate_evidence"] += bool(row["annotations"]["family"]["values"])
        by_kind[kind]["query_assignment_count"] += row["annotations"]["family"]["value"] is not None
    return {"kind": "synthetic_engineering_controls_not_biological_specificity", "random_seed": 20260930,
            "query_count": len(records), "control_design": "100 length-matched uniform random; 100 mononucleotide-composition-preserving shuffles; 4 length-30 homopolymers",
            "status_counts_by_control": {kind: dict(values) for kind, values in sorted(by_kind.items())},
            "v2_wall_seconds_including_reference_validation": seconds}


def formatted_rate(value):
    return "unavailable (zero eligible cases)" if value is None else f"{value:.3f}"


def markdown_report(result):
    original, self_result, heldout, negative = (result[k] for k in ("original_comparison", "reference_self_recovery", "heldout", "synthetic_controls"))
    fam = original["metrics"]["family"]
    motif = original["metrics"]["motif"]
    rows = ["# Benchmark against recovered original CRISPRmap evidence", "", f"Generated: {result['generated_utc']}", "",
            "This is an original-output comparison, not a rerun of the original executable. The original service was offline and no local source release was recovered. Original runtime and speedup are unavailable.", "",
            f"Reference: `{result['reference_namespace']}`; {self_result['reference_record_count']:,} records. SHA256: `{result['reference_sha256']}`.", "",
            "## Original published 16-query example", "",
            "| Label | Original known labels | v2 assigned | Assigned agreement |", "| --- | ---: | ---: | ---: |",
            f"| Sequence family | {fam['known_original_labels']} | {fam['assigned_known_labels']} | {fam['correct_assigned_known_labels']}/{fam['assigned_known_labels']} |",
            f"| Structural motif | {motif['known_original_labels']} | {motif['assigned_known_labels']} | {motif['correct_assigned_known_labels']}/{motif['assigned_known_labels']} |", "",
            f"v2 statuses: `{json.dumps(original['status_counts'], sort_keys=True)}`. v2 elapsed: {original['v2_wall_seconds_including_reference_validation']:.3f} seconds.", "",
            "Known-label coverage is reported because abstention is not a correct positive assignment. Null-to-null matches are retained in JSON but are not counted as positive-label agreement. Original superclass inference is unsupported and excluded. Nonexact structural motifs are intentionally not assigned.", "",
            "The two positive-label differences are explicit: Ecoli_K12 matches reference occurrences with competing families F2/F4 and motifs M1/M2, so v2 abstains; the original selected F4/M2. Stepidermidis has original motif M28 but is nonexact, and v2 does not infer its structural motif. This tool does not reproduce the original structure-aware assignment capability.", "",
            "## Reference recovery", "",
            f"All {self_result['query_count']:,} reference and reverse-complement queries: {self_result['exact_recovery_count']:,} exact matches. Annotation evidence was orientation-invariant for {self_result['orientation_annotation_invariant_pairs']:,}/{self_result['reference_record_count']:,} pairs. This is an engineering self-recovery check, not independent validation.", "",
            "## Sequence-group holdout", "",
            f"Evaluated {heldout['query_count']} canonical sequence/reverse-complement groups. Split: {heldout['split_rule']}. Training record counts per fold: {heldout['training_records_per_fold']}. Each query group and every occurrence of its reverse complement were removed from its training fold. Closely related sequences may occur on both sides; species, homology and family separation were not imposed.", "",
            f"Counts: `{json.dumps(heldout['counts'], sort_keys=True)}`.", "",
            f"Candidate coverage for unambiguous historical families: {formatted_rate(heldout['candidate_coverage_on_unambiguous_labels'])}. True historical family present among candidates: {formatted_rate(heldout['true_label_candidate_recall_on_unambiguous_labels'])}. Agreement among single-candidate queries: {formatted_rate(heldout['agreement_among_sole_candidate_queries'])}.", "",
            f"Macro family metrics (equal family weight): `{json.dumps(heldout['macro_per_family'], sort_keys=True)}`.", "",
            "These are historical candidate-agreement metrics, not calibrated biological accuracy. The single-candidate metric requires complete unanimous candidate evidence; incomplete or conflicting reference evidence is excluded. Missing or contradictory query labels are excluded from those denominators and counted separately. The maximum edit distance was fixed before evaluation; no test-set tuning was performed. Per-family counts and per-query hashes are in the JSON.", "",
            "## Family-absent stress test", "",
            f"Removed each of {result['family_absent']['family_count']} historical families in turn, excluding complete duplicate/reverse-complement groups. Across {result['family_absent']['query_count']} unambiguous query groups, {result['family_absent']['cross_family_candidate_queries']} retrieved evidence for another family ({result['family_absent']['cross_family_candidate_rate']:.3f}). Query family assignments: {result['family_absent']['query_assignment_count']}; zero assigned labels follows from the nonexact-abstention design and is not evidence of biological specificity.", "",
            "Cross-family candidates are not validated family assignments and may reflect genuine similarity across historical partitions. This stress test quantifies candidate retrieval when the source family is absent; it does not calibrate biological novelty detection.", "",
            "## Synthetic controls and timing", "", f"Control outcomes: `{json.dumps(negative['status_counts_by_control'], sort_keys=True)}`.", "",
            "Synthetic sequences are engineering controls, not experimentally established non-CRISPR negatives. Random controls match length, and shuffles preserve mononucleotide composition; no dinucleotide preservation is claimed.", "",
            "| Workload | v2 wall seconds |", "| --- | ---: |"]
    for key in ("original_comparison", "reference_self_recovery", "heldout", "family_absent", "synthetic_controls"):
        rows.append(f"| {key} | {result[key]['v2_wall_seconds_including_reference_validation']:.3f} |")
    rows += ["", "Timing uses `time.perf_counter` around the in-process search including reference validation; it excludes download, file loading and interpreter startup. One measured run per workload; no confidence interval or original-tool speedup is claimed. Hardware can differ between runs.", "", "## Reproduce", "", "```bash", "python scripts/fetch_legacy.py", "python benchmarks/run_benchmark.py --output-dir benchmarks/results", "```", "", "The raw third-party assets stay in ignored data/. Public result files contain derived scientific comparison counts and query identifiers/hashes, not downloaded reference sequences. See docs/LEGACY_PROVENANCE.md and docs/SCIENCE.md for source and interpretation details.", ""]
    return "\n".join(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data" / "legacy" / "2013")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "benchmarks" / "results")
    parser.add_argument("--heldout-groups", type=int, default=0, help="0 evaluates every canonical group across five folds; positive selects that many groups")
    parser.add_argument("--max-distance", type=int, default=3)
    args = parser.parse_args(argv)
    if not 0 <= args.heldout_groups <= 2000:
        parser.error("--heldout-groups must be 0..2000")
    reference, _ = load_reference(args.data_dir / "reference.json")
    expected = legacy.build_legacy_reference(args.data_dir)
    if reference != expected:
        parser.error("reference.json differs from the reference rebuilt from pinned original bytes and labels")
    if reference["provenance"]["namespace"] != legacy.NAMESPACE:
        parser.error("This original-output benchmark requires the 2013 namespace")
    result = {"schema": "crisprmap.benchmark.v1", "generated_utc": datetime.now(timezone.utc).isoformat(),
              "software_version": __import__("crisprmap").__version__, "python": platform.python_version(),
              "platform": {"system": platform.system(), "machine": platform.machine()},
              "reference_namespace": reference["provenance"]["namespace"], "reference_sha256": reference["content_sha256"],
              "method": {"max_distance": args.max_distance, "distance_policy": "fixed_uncalibrated", "top_k": 20, "timing_repetitions": 1}}
    for name, function in [
        ("original_comparison", lambda: original_comparison(args.data_dir, reference, args.max_distance)),
        ("reference_self_recovery", lambda: self_recovery(reference, args.max_distance)),
        ("heldout", lambda: heldout_comparison(reference, args.heldout_groups, args.max_distance)),
        ("family_absent", lambda: family_absent(reference, args.max_distance)),
        ("synthetic_controls", lambda: controls(reference, args.max_distance)),
    ]:
        print(f"Running {name}...", flush=True)
        result[name] = function()
        print(f"Completed {name}: {result[name]['v2_wall_seconds_including_reference_validation']:.3f}s", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    legacy.atomic_write(args.output_dir / "benchmark.json", (json.dumps(result, indent=2, sort_keys=True) + "\n").encode())
    legacy.atomic_write(args.output_dir / "REPORT.md", markdown_report(result).encode())
    print(args.output_dir / "REPORT.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
