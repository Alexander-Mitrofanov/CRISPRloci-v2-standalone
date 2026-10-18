"""Diagnose natural-locus misses using retained evidence and cached native runs.

Run with PYTHONPATH=src. Isolating a known RNA is a truth-assisted diagnostic,
never a genome-recall experiment. Both repeat orientations are always queried.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from collections import defaultdict
from pathlib import Path

from summarize_literature_benchmark import localization_match
from summarize_sequence_benchmark import family

from crisprtracrrna.config import load_settings
from crisprtracrrna.process import Runner
from crisprtracrrna.sequence import reverse_complement

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks"


def fasta(path):
    result = {}
    key = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            key = line[1:].split()[0]
            result[key] = ""
        elif key:
            result[key] += line.strip()
    return result


def overlap(a, b):
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)


def distance(a, b):
    return max(0, b["start"] - a["end"] - 1, a["start"] - b["end"] - 1)


def parse_alignments(path, queries, truth, config, arrays, offset=0):
    result = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        c = line.split()
        if len(c) != 12:
            raise ValueError(f"Invalid alignment row in {path}: {line}")
        qs, qe, ss, se = map(int, c[6:10])
        if not 0 < qs <= len(queries[c[0]]) or not 0 < qe <= len(queries[c[0]]):
            raise ValueError(f"Invalid query coordinates in {path}")
        interval = {"start": min(ss, se) + offset, "end": max(ss, se) + offset}
        identity = float(c[2]) / 100
        coverage = (abs(qe - qs) + 1) / len(queries[c[0]])
        if not 0 <= identity <= 1:
            raise ValueError(f"Invalid identity in {path}")
        blocked = []
        if coverage < config["coverage"]:
            blocked.append("coverage")
        if identity < config["similarity"]:
            blocked.append("identity")
        if any(overlap(interval, a["interval"]) for a in arrays):
            blocked.append("array_overlap")
        result.append(
            {
                "query": c[0],
                "repeat_length": len(queries[c[0]]),
                "query_span": abs(qe - qs) + 1,
                "interval": interval,
                "identity": identity,
                "coverage": coverage,
                "bitscore": float(c[11]),
                "evalue": float(c[10]),
                "truth_overlap_nt": overlap(interval, truth),
                "truth_distance_nt": distance(interval, truth),
                "blocked_by": blocked,
            }
        )
    return result


def diagnostic_alignments(entry, work, output, settings, rerun):
    output.mkdir(parents=True, exist_ok=True)
    for name in ("repeats.fasta", "anti-search.fasta"):
        shutil.copyfile(work / name, output / name)
    forward = entry["sequence"]
    if entry["matches"][0]["strand"] == "-":
        forward = reverse_complement(forward)
    (output / "isolated-reference.fasta").write_text(">reference\n" + forward + "\n")
    (output / "published-repeats.fasta").write_text(
        ">published_f\n"
        + entry["repeat"]
        + "\n>published_r\n"
        + reverse_complement(entry["repeat"])
        + "\n"
    )
    runner = Runner(output / "commands", settings.timeout, settings.environment)
    blast = [
        settings.blastn,
        "-task",
        "blastn-short",
        "-query",
        "repeats.fasta",
        "-strand",
        "both",
        "-dust",
        "no",
        "-soft_masking",
        "false",
        "-evalue",
        "100",
        "-outfmt",
        "6",
    ]
    commands = {
        "isolated-fasta.tsv": [
            settings.fasta36,
            "-T",
            "1",
            "-m",
            "8",
            "repeats.fasta",
            "isolated-reference.fasta",
        ],
        "genomic-blast.tsv": [*blast, "-subject", "anti-search.fasta"],
        "isolated-blast.tsv": [*blast, "-subject", "isolated-reference.fasta"],
        "published-anti-search.tsv": [
            settings.fasta36,
            "-T",
            "1",
            "-m",
            "8",
            "published-repeats.fasta",
            "anti-search.fasta",
        ],
        "published-isolated-reference.tsv": [
            settings.fasta36,
            "-T",
            "1",
            "-m",
            "8",
            "published-repeats.fasta",
            "isolated-reference.fasta",
        ],
    }
    for name, command in commands.items():
        if rerun or not (output / name).exists():
            runner.run(command, cwd=output, stdout=output / name)


def align_summary(rows):
    relevant = [r for r in rows if r["truth_overlap_nt"]]
    passing_identity = [r for r in relevant if "identity" not in r["blocked_by"]]
    best = max(
        passing_identity or relevant, key=lambda r: (r["bitscore"], r["coverage"]), default=None
    )
    return {
        "total_alignments": len(rows),
        "truth_overlapping_alignments": len(relevant),
        "truth_overlapping_passing_default_gates": sum(not r["blocked_by"] for r in relevant),
        "max_truth_coverage_among_identity_passing": max(
            (r["coverage"] for r in passing_identity), default=0
        ),
        "best_truth_alignment": best,
        "truth_alignments": relevant,
    }


def make_diagnosis(entry, stage, settings, rerun):
    base = BENCH / stage / entry["label"]
    report = json.loads((base / "report.json").read_text())
    record = report["records"][0]
    work = base / "work" / record["key"]
    truth = entry["matches"][0]
    identify_path = next((work / "identify").glob("*/report.json"))
    identify = json.loads(identify_path.read_text())
    arrays = [
        {
            "id": a["id"],
            "interval": a["source_interval"],
            "repeat": a["consensus_repeat"],
            "category": a["category"],
            "score": a["certainty_score"],
            "published_repeat_exact_either_orientation": entry["repeat"]
            in (a["consensus_repeat"], reverse_complement(a["consensus_repeat"])),
            "truth_overlap_nt": overlap(a["source_interval"], truth),
        }
        for a in identify["arrays"]
        if a["id"] in record["arrays"]
    ]
    config = report["configuration"]
    queries = fasta(work / "repeats.fasta")
    out = BENCH / "diagnostics/final_failures" / entry["label"]
    diagnostic_alignments(entry, work, out, settings, rerun)
    summaries = {}
    for name, path, offset in [
        ("pipeline_fasta_genome", work / "anti-repeat.tsv", 0),
        ("diagnostic_blast_genome", out / "genomic-blast.tsv", 0),
        ("diagnostic_fasta_isolated", out / "isolated-fasta.tsv", truth["start"] - 1),
        ("diagnostic_blast_isolated", out / "isolated-blast.tsv", truth["start"] - 1),
    ]:
        summaries[name] = align_summary(
            parse_alignments(path, queries, truth, config, arrays, offset)
        )
    published_queries = fasta(out / "published-repeats.fasta")
    for name, path, offset in [
        ("published_repeat_fasta_genome", out / "published-anti-search.tsv", 0),
        (
            "published_repeat_fasta_isolated",
            out / "published-isolated-reference.tsv",
            truth["start"] - 1,
        ),
    ]:
        summaries[name] = align_summary(
            parse_alignments(path, published_queries, truth, config, arrays, offset)
        )
    matching_cms = [
        h for h in record["raw_genomic_model_hits"] if localization_match(h, truth, truth["contig"])
    ]
    wrong_strand_cms = [
        h
        for h in record["raw_genomic_model_hits"]
        if h["strand"] != truth["strand"] and overlap(h["interval"], truth) >= 20
    ]
    best_cm = min(matching_cms, key=lambda h: h["evalue"], default=None)
    nearest = min(
        (p for p in record["predictions"] if p["strand"] == truth["strand"]),
        key=lambda p: distance(p["interval"], truth),
        default=None,
    )
    nearest_summary = (
        None
        if nearest is None
        else {
            "id": nearest["id"],
            "interval": nearest["interval"],
            "strand": nearest["strand"],
            "distance_nt": distance(nearest["interval"], truth),
            "truth_overlap_nt": overlap(nearest["interval"], truth),
            "origin": nearest["candidate_origin"],
        }
    )
    primary = summaries["pipeline_fasta_genome"]
    blast = summaries["diagnostic_blast_genome"]
    if primary["truth_overlapping_passing_default_gates"]:
        reason = "candidate_placement_requires_review"
    elif blast["truth_overlapping_passing_default_gates"]:
        reason = "fasta_genomic_reporting_or_scoring_gap"
    elif any(x["max_truth_coverage_among_identity_passing"] > 0 for x in summaries.values()):
        reason = "short_repeat_similarity_below_coverage_gate"
    else:
        reason = "no_identity_passing_repeat_alignment"
    paths = [
        base / "report.json",
        identify_path,
        work / "anti-repeat.tsv",
        work / "repeats.fasta",
        out / "genomic-blast.tsv",
        out / "isolated-fasta.tsv",
        out / "isolated-blast.tsv",
        out / "published-repeats.fasta",
        out / "published-anti-search.tsv",
        out / "published-isolated-reference.tsv",
    ]
    return {
        "label": entry["label"],
        "reference_id": entry["id"],
        "truth": truth,
        "reference_length": len(entry["sequence"]),
        "published_repeat": entry["repeat"],
        "successful_execution": report["status"] == "complete",
        "arrays": arrays,
        "any_exact_repeat_either_orientation": any(
            a["published_repeat_exact_either_orientation"] for a in arrays
        ),
        "array_masks_truth": any(a["truth_overlap_nt"] for a in arrays),
        "published_repeat_substitution_rescues_genomic_alignment": bool(
            summaries["published_repeat_fasta_genome"]["truth_overlapping_passing_default_gates"]
        ),
        "settings": {
            k: config[k]
            for k in (
                "similarity",
                "coverage",
                "mask_arrays",
                "model_rescue",
                "model_rescue_evalue",
                "tail_length",
                "upstream_length",
            )
        },
        "anti_repeat_diagnosis": reason,
        "alignment_diagnostics": summaries,
        "genomic_model_localizing_hits": matching_cms,
        "genomic_model_wrong_strand_overlapping_hits": wrong_strand_cms,
        "best_localizing_model_hit": best_cm,
        "model_diagnosis": "localized_hit_above_rescue_evalue_cutoff"
        if best_cm
        else "no_localizing_hit_reported_by_existing_models",
        "candidate_count": len(record["predictions"]),
        "nearest_same_strand_candidate": nearest_summary,
        "ranking_failure": False,
        "evidence_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
        },
    }


def partial_rescue(entry, report):
    truth = entry["matches"][0]
    predictions = sorted(report["records"][0]["predictions"], key=lambda p: -p["ranking"]["score"])
    matches = [
        (rank, p, localization_match(p, truth, truth["contig"]))
        for rank, p in enumerate(predictions, 1)
    ]
    matches = [(rank, p, ov) for rank, p, ov in matches if ov]
    if not matches or max(ov / len(entry["sequence"]) for _, _, ov in matches) >= 0.9:
        return None
    rank, p, ov = max(matches, key=lambda r: r[2])
    return {
        "label": entry["label"],
        "truth": truth,
        "reference_length": len(entry["sequence"]),
        "prediction_interval": p["interval"],
        "origin": p["candidate_origin"],
        "boundary_kind": p["boundary_kind"],
        "truth_overlap_nt": ov,
        "reference_coverage": ov / len(entry["sequence"]),
        "rank": rank,
        "best_model_evalue": min(h["evalue"] for h in p["model_hits"]),
        "explanation": (
            "The CM identifies a conserved subregion. No passing anti-repeat supplies"
            " an oriented extension window; the rescue reports only model-supported "
            "coordinates and does not assert complete transcript boundaries."
        ),
    }


def sequence_thresholds():
    """Join canonical families across batches, retaining first-batch duplicate controls."""
    members = defaultdict(set)
    with (BENCH / "literature/catalog.jsonl").open() as handle:
        for line in handle:
            row = json.loads(line)
            model = family(row["subtype"])
            if model:
                members[model].add(row["sequence_sha256"])
    result = {}
    canonical = json.loads((BENCH / "results/sequence_summary.json").read_text())
    batches = canonical["scan_batches"]
    for model_type in ("II", "V"):
        lookup = {}
        inputs = {}
        for batch in batches:
            base = BENCH / batch / model_type
            ids = json.loads((base / "sequence_ids.json").read_text())
            if not ids:
                continue
            hits = defaultdict(list)
            for hit in json.loads((base / "hits.json").read_text()):
                hits[hit["query"]].append(hit)
            for digest, query in ids.items():
                if digest in members[model_type]:
                    lookup.setdefault(digest, (hits[query], hits[query + "_shuffle"]))
            for name in ("hits.json", "sequence_ids.json", "cm-input.fasta"):
                path = base / name
                inputs[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        if set(lookup) != members[model_type]:
            raise ValueError(f"Missing canonical sequence scans for {model_type}")
        tables = []
        for cutoff in (1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1, 10):
            tables.append(
                {
                    "threshold": cutoff,
                    "unique_sequences": len(lookup),
                    "recognized_correct_strand_span20": sum(
                        any(
                            h["evalue"] <= cutoff
                            and h["strand"] == "+"
                            and h["interval"]["length"] >= 20
                            for h in true_hits
                        )
                        for true_hits, _ in lookup.values()
                    ),
                    "shuffled_sequences_any_hit": sum(
                        any(h["evalue"] <= cutoff for h in shuffled_hits)
                        for _, shuffled_hits in lookup.values()
                    ),
                    "shuffled_sequences_any_strand_span20": sum(
                        any(
                            h["evalue"] <= cutoff and h["interval"]["length"] >= 20
                            for h in shuffled_hits
                        )
                        for _, shuffled_hits in lookup.values()
                    ),
                }
            )
        result[model_type] = {
            "thresholds": tables,
            "membership": (
                "Canonical catalog membership: II, II-* or unknown for the Type II screen; "
                "V or V-* for Type V. Type III and VI are excluded. Both retained scan "
                "batches are joined by sequence hash; the first batch wins duplicate hashes."
            ),
            "caveat": (
                "E-values depend on search space: this pooled sequence scan is not a "
                "calibration of the per-genome rescue E-value. Mononucleotide shuffled "
                "sequences are not verified biological negatives."
            ),
            "input_sha256": inputs,
        }
    return result


def write_outputs(rows, rescues, thresholds, stage, out):
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": stage,
        "methods": (
            "Retained CRISPRidentify, FASTA36 and cmscan outputs; matching requires "
            "same contig/strand, >=20 nt and >=50% shorter interval. Alternate BLAST "
            "uses identical identity/coverage gates and masks. Both repeat "
            "orientations are queried. Isolated-reference alignments use known truth "
            "and are diagnostic only."
        ),
        "limits": (
            "Absence of a cmscan table hit does not prove biological absence of a "
            "family. Scores and candidate counts are not calibrated precision "
            "estimates. No parameters or models were fitted to these references."
        ),
        "misses": rows,
        "partial_boundary_rescues": rescues,
        "model_threshold_diagnostic": thresholds,
        "pooled_sequence_threshold_diagnostic": sequence_thresholds(),
    }
    (out / "failure_diagnosis.json").write_text(json.dumps(payload, indent=2) + "\n")
    table = []
    for r in rows:
        cm = r["best_localizing_model_hit"]
        d = r["alignment_diagnostics"]
        iso = d["diagnostic_fasta_isolated"]["best_truth_alignment"]
        native = d["pipeline_fasta_genome"]
        table.append(
            {
                "label": r["label"],
                "reference_length": r["reference_length"],
                "accepted_arrays": len(r["arrays"]),
                "exact_repeat_either_orientation": r["any_exact_repeat_either_orientation"],
                "pipeline_alignments_at_truth": native["truth_overlapping_alignments"],
                "pipeline_passing_at_truth": native["truth_overlapping_passing_default_gates"],
                "blast_passing_at_truth": d["diagnostic_blast_genome"][
                    "truth_overlapping_passing_default_gates"
                ],
                "published_repeat_substitution_passing_at_truth": d[
                    "published_repeat_fasta_genome"
                ]["truth_overlapping_passing_default_gates"],
                "isolated_fasta_query_span": iso["query_span"] if iso else "",
                "isolated_fasta_repeat_length": iso["repeat_length"] if iso else "",
                "isolated_fasta_identity": iso["identity"] if iso else "",
                "isolated_fasta_coverage": iso["coverage"] if iso else "",
                "anti_repeat_diagnosis": r["anti_repeat_diagnosis"],
                "best_localizing_cm_evalue": cm["evalue"] if cm else "",
                "model_diagnosis": r["model_diagnosis"],
                "candidate_count": r["candidate_count"],
                "nearest_same_strand_candidate_distance": r["nearest_same_strand_candidate"][
                    "distance_nt"
                ]
                if r["nearest_same_strand_candidate"]
                else "",
            }
        )
    with (out / "failure_diagnosis.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(table[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(table)
    text = [
        "# Natural-locus failure analysis",
        "",
        (
            f"Stage: `{stage}`. All {len(rows)} misses ran successfully. "
            "None is solely a ranking failure: no returned candidate meets the fixed "
            "localization criterion. Every locus has an accepted array; all truth "
            "intervals lie outside array masks."
        ),
        "",
        (
            "Diagnostic searches keep both repeat orientations. Isolated-RNA searches"
            " are truth-assisted explanations, not additional benchmark detections. "
            "Native genome BLAST uses the same masked input, 70% identity and 60% "
            "query-span coverage gates."
        ),
        "",
        (
            "| Locus | Repeat evidence | Native FASTA alignments at truth | Isolated "
            "FASTA span; identity | CM at truth | Interpretation |"
        ),
        "|---|---|---:|---|---|---|",
    ]
    for r, t in zip(rows, table, strict=True):
        repeat = (
            "Exact published repeat (either strand)"
            if r["any_exact_repeat_either_orientation"]
            else "Inferred repeat differs from published repeat"
        )
        align = (
            f"{t['isolated_fasta_query_span']}/{t['isolated_fasta_repeat_length']}; "
            f"{100 * t['isolated_fasta_identity']:.1f}%"
            if t["isolated_fasta_identity"] != ""
            else "None"
        )
        cm = (
            f"E={t['best_localizing_cm_evalue']:g}, fails 1e-5 gate"
            if t["best_localizing_cm_evalue"] != ""
            else "No localizing hit reported"
        )
        text.append(
            f"| {r['label']} | {repeat} | {t['pipeline_alignments_at_truth']} "
            f"(passing: {t['pipeline_passing_at_truth']}) | {align} | {cm} | "
            f"{r['anti_repeat_diagnosis'].replace('_', ' ')} |"
        )
    text += [
        "",
        (
            "The isolated FASTA column reports the highest-scoring identity-passing "
            "alignment (or highest-scoring available alignment), not the maximum "
            "coverage across every alignment. Full per-alignment evidence, native "
            "BLAST comparisons, exact coordinates, configuration and source hashes "
            "are in `failure_diagnosis.json`."
        ),
        "",
        (
            "Eight loci have native FASTA alignments at the known RNA, rejected for "
            "insufficient repeat coverage. Esp2 and Bbo have no native FASTA genomic "
            "alignment at the RNA; BLAST and isolated-RNA diagnostics expose only "
            "subthreshold matches. Bbo additionally has a 26/38 nt isolated match at "
            "69.23% identity, failing the 70% identity gate. The accepted consensus "
            "differs from the published repeat in Bbo, Cca1, Ffr and Psp, but "
            "substituting the published repeat (both orientations) rescues none at "
            "default gates. Array inference differences therefore do not explain "
            "these failures by themselves. No correct-strand CM matches are confused "
            "with wrong-strand hits: Cca1, Cme1, Orh and Psp have correct-strand hits"
            " that fail the rescue threshold; the other six have no localizing hit in"
            " the raw default cmscan output."
        ),
        "",
        "## Partial boundaries from model rescue",
        "",
        "| Locus | Model-supported interval | Covered reference | Rank | E-value |",
        "|---|---|---|---:|---:|",
    ]
    for r in rescues:
        p = r["prediction_interval"]
        text.append(
            f"| {r['label']} | {p['start']}–{p['end']} ({r['truth']['strand']}) | "
            f"{r['truth_overlap_nt']}/{r['reference_length']} "
            f"({100 * r['reference_coverage']:.1f}%) | {r['rank']} | {r['best_model_evalue']:g} |"
        )
    text += [
        "",
        (
            "These four loci pass localization because a conserved subregion is "
            "detected, but do not pass the separate 90% reference-coverage metric. "
            "Their `boundary_kind=model_supported_region` is deliberate: extending to"
            " the complete transcript without anti-repeat or boundary evidence would "
            "invent precision."
        ),
        "",
        "## Model-threshold diagnostic",
        "",
        (
            "The following is a rescore of raw genomic CM hits, not a rerun or a "
            "recommended cutoff. The union with anti-repeat predictions is potential "
            "localization before rescue suppression/deduplication, not measured final"
            " pipeline recall. Counts before array exclusion/overlap deduplication "
            "expose the increasing evidence burden; unmatched hits are not "
            "established false positives."
        ),
        "",
        (
            "| Maximum E-value | Localized references by CM / 78 | Anti-repeat or raw"
            " CM union / 78 | Raw CM hits / 78 loci |"
        ),
        "|---:|---:|---:|---:|",
    ]
    for t in thresholds:
        text.append(
            f"| {t['threshold']} | {t['model_detected']} | "
            f"{t['union_with_anti_repeat_detected']} | {t['raw_hit_count']} |"
        )
    pooled = payload["pooled_sequence_threshold_diagnostic"].get("II")
    if pooled:
        text += [
            "",
            (
                "Canonical Type II/unknown sequence controls combine the retained batches "
                "by sequence hash, retaining the first batch for duplicates. Membership "
                "is exactly II, II-* or unknown; Type III is excluded. They provide a "
                "separate descriptive check. E-values depend on search-space size, "
                "so these values cannot "
                "calibrate the per-genome rescue threshold. Shuffles retain "
                "mononucleotide composition and are not validated biological negatives."
            ),
            "",
            (
                "| Maximum E-value | Correct-strand recognized / unique inputs | Shuffles"
                " with any reported hit / same denominator |"
            ),
            "|---:|---:|---:|",
        ]
        for t in pooled["thresholds"]:
            text.append(
                f"| {t['threshold']} | {t['recognized_correct_strand_span20']}/"
                f"{t['unique_sequences']} | {t['shuffled_sequences_any_hit']}/"
                f"{t['unique_sequences']} |"
            )
    text += [
        "",
        "## Follow-up decisions",
        "",
        (
            "Do not relabel these misses as ranking failures or silently lower gates."
            " If an alternate backend reports an otherwise eligible genomic "
            "alignment, evaluate it on the complete reference set and suitable "
            "negatives before changing the default. Shorter anti-repeat coverage and "
            "weaker CM support call for an explicitly calibrated sensitivity profile "
            "and independently held-out families. Missing native transcript "
            "boundaries require separate evidence. This positive-only, "
            "training-overlapping collection cannot establish specificity or select a"
            " new cutoff."
        ),
        "",
        (
            "Reproduce with `PYTHONPATH=src python "
            "scripts/diagnose_literature_failures.py`; add `--rerun-alignments` to "
            "rerun the cached native diagnostics."
        ),
    ]
    (out / "failure_analysis.md").write_text("\n".join(text) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", default="final_fixed/complete")
    parser.add_argument("--rerun-alignments", action="store_true")
    args = parser.parse_args()
    entries = json.loads((BENCH / "literature/loci.json").read_text())
    settings = load_settings(ROOT / "local.toml", {"threads": 1})
    rows, rescues = [], []
    thresholds = [
        {
            "threshold": e,
            "model_detected": 0,
            "union_with_anti_repeat_detected": 0,
            "raw_hit_count": 0,
        }
        for e in (1e-5, 1e-4, 1e-3, 1e-2, 0.1, 1, 10)
    ]
    for entry in entries:
        if entry["status"] != "exact_unique":
            continue
        report = json.loads((BENCH / args.stage / entry["label"] / "report.json").read_text())
        truth = entry["matches"][0]
        detected = any(
            localization_match(p, truth, r["source"]["id"])
            for r in report["records"]
            for p in r["predictions"]
        )
        anti_detected = any(
            localization_match(p, truth, r["source"]["id"])
            for r in report["records"]
            for p in r["predictions"]
            if p["candidate_origin"] == "anti_repeat"
        )
        for t in thresholds:
            hits = [
                h
                for r in report["records"]
                for h in r["raw_genomic_model_hits"]
                if h["evalue"] <= t["threshold"]
            ]
            localized = any(localization_match(h, truth, truth["contig"]) for h in hits)
            t["model_detected"] += localized
            t["union_with_anti_repeat_detected"] += anti_detected or localized
            t["raw_hit_count"] += len(hits)
        if not detected:
            rows.append(make_diagnosis(entry, args.stage, settings, args.rerun_alignments))
        else:
            rescue = partial_rescue(entry, report)
            if rescue:
                rescues.append(rescue)
    write_outputs(rows, rescues, thresholds, args.stage, BENCH / "results")
    print(
        json.dumps(
            {
                "misses": len(rows),
                "partial_boundary_rescues": len(rescues),
                "diagnoses": {r["label"]: r["anti_repeat_diagnosis"] for r in rows},
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
