"""Frozen round-two execution and controls; outputs never touch historical results.

The panels are development/regression evidence, not biological negatives or an
independent evaluation. Detector calls receive only FASTA and ordinary settings.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import re
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, replace
from datetime import UTC, datetime
from hashlib import sha256
from pathlib import Path

from summarize_literature_benchmark import localization_match

from crisprtracrrna.adapters import parse_identify, sha256_file
from crisprtracrrna.config import load_settings
from crisprtracrrna.evidence import candidate_evidence, evidence_class, short_rescue_admission
from crisprtracrrna.models import AntiRepeat, Array, Candidate, Interval, Sequence
from crisprtracrrna.pipeline import atomic_json, command_provenance, doctor, run
from crisprtracrrna.process import Runner
from crisprtracrrna.search import (
    interactions,
    make_candidates,
    scan_models,
    search_anti_repeats,
    terminators,
)
from crisprtracrrna.sequence import read_fasta, reverse_complement, write_fasta

ROOT = Path(__file__).resolve().parents[1]
ROUND = ROOT / "benchmarks/council_round2"
HISTORICAL = ROOT / "benchmarks/final_fixed/complete"
NAMESPACE = "crisprtracrrna-council-controls-v1"
MASTER_SEED = 20260927


def json_read(path):
    return json.loads(path.read_text())


def seed(source_hash, kind, replicate, pool="II"):
    value = f"{NAMESPACE}|{MASTER_SEED}|{pool}|{source_hash}|{kind}|{replicate}"
    return int(sha256(value.encode()).hexdigest(), 16)


def mono_shuffle(sequence, rng):
    positions = [i for i, base in enumerate(sequence) if base in "ACGT"]
    values = [sequence[i] for i in positions]
    rng.shuffle(values)
    result = list(sequence)
    for i, base in zip(positions, values, strict=True):
        result[i] = base
    return "".join(result)


def dinucleotide_shuffle(sequence, rng):
    """Randomized Euler trails preserve exact directed edge counts and endpoints.

    Ambiguities partition canonical runs and remain fixed. Each canonical run's
    endpoints remain fixed, so adjacent pairs across ambiguity boundaries do too.
    This is not claimed to sample uniformly from all possible shuffles.
    """
    result = list(sequence)
    for match in re.finditer("[ACGT]+", sequence):
        original = match.group()
        adjacency = defaultdict(list)
        for left, right in zip(original, original[1:], strict=False):
            adjacency[left].append(right)
        for values in adjacency.values():
            rng.shuffle(values)
        stack, path = [original[0]], []
        while stack:
            if adjacency[stack[-1]]:
                stack.append(adjacency[stack[-1]].pop())
            else:
                path.append(stack.pop())
        shuffled = "".join(reversed(path))
        assert len(shuffled) == len(original)
        assert Counter(zip(shuffled, shuffled[1:], strict=False)) == Counter(
            zip(original, original[1:], strict=False)
        )
        result[match.start() : match.end()] = shuffled
    return "".join(result)


def inputs():
    # Do not pass reference sequence, published repeat or truth to the detector.
    return [
        {"label": e["label"], "fasta": e["fasta"]}
        for e in json_read(ROOT / "benchmarks/literature/loci.json")
    ]


def frozen_context(label):
    work = next((HISTORICAL / label / "work").iterdir())
    source = next(read_fasta(work / "input.fasta"))
    arrays, _ = parse_identify(next((work / "identify").rglob("report.json")), source)
    return source, arrays


def donor_assignment(pool):
    """Deterministic bipartite matching on a frozen entire donor pool."""
    by_id = {x["identity"]: x for x in pool}
    choices = {}
    for identity, item in by_id.items():
        choices[identity] = sorted(
            (
                d
                for d in by_id
                if by_id[d]["label"] != item["label"]
                and by_id[d]["repeat"] not in {item["repeat"], reverse_complement(item["repeat"])}
            ),
            key=lambda d: sha256(f"{NAMESPACE}|{MASTER_SEED}|{identity}|{d}".encode()).hexdigest(),
        )
    assigned = {}

    def augment(recipient, visited):
        for donor in choices[recipient]:
            if donor in visited:
                continue
            visited.add(donor)
            if donor not in assigned or augment(assigned[donor], visited):
                assigned[donor] = recipient
                return True
        return False

    for identity in sorted(by_id):
        if not augment(identity, set()):
            raise ValueError("No full cross-locus repeat derangement exists")
    return {recipient: donor for donor, recipient in assigned.items()}


def prepare():
    output = ROUND / "controls"
    output.mkdir(exist_ok=False)
    protocol = json_read(ROUND / "protocol.json")
    entries = inputs()
    pool, contexts = [], {}
    for entry in entries:
        source, arrays = frozen_context(entry["label"])
        contexts[entry["label"]] = (source, arrays)
        for array in arrays:
            pool.append(
                {
                    "identity": entry["label"] + ":" + array.id,
                    "label": entry["label"],
                    "array_id": array.id,
                    "repeat": array.repeat,
                    "source_sha256": source.sha256,
                }
            )
    pool.sort(key=lambda x: x["identity"])
    assignment = donor_assignment(pool)
    by_id = {x["identity"]: x for x in pool}
    repeats = []
    for item in pool:
        for replicate in (1, 2):
            kind = "repeat_dinucleotide:" + item["array_id"]
            control_seed = seed(item["source_sha256"], kind, replicate)
            shuffled = dinucleotide_shuffle(item["repeat"], random.Random(control_seed))
            repeats.append(
                {
                    **item,
                    "panel": f"dinucleotide_{replicate}",
                    "replicate": replicate,
                    "seed": str(control_seed),
                    "control_kind": kind,
                    "control_repeat": shuffled,
                    "unchanged": shuffled == item["repeat"],
                    "donor": None,
                }
            )
        donor = by_id[assignment[item["identity"]]]
        repeats.append(
            {
                **item,
                "panel": "derangement",
                "control_kind": "repeat_derangement",
                "control_repeat": donor["repeat"],
                "donor": donor["identity"],
                "unchanged": False,
            }
        )
    genomic = []
    for entry in entries:
        if entry["label"] not in protocol["controls"]["genomic_labels"]:
            continue
        source, arrays = contexts[entry["label"]]
        protected = [False] * len(source.bases)
        for array in arrays:
            protected[array.interval.start : array.interval.end] = [True] * array.interval.length
        changed = list(source.bases)
        segments, start = [], 0
        while start < len(changed):
            if protected[start]:
                start += 1
                continue
            end = start + 1
            while end < len(changed) and not protected[end]:
                end += 1
            control_seed = seed(source.sha256, f"genomic_nonarray:{start}:{end}", 1)
            original = source.bases[start:end]
            shuffled = mono_shuffle(original, random.Random(control_seed))
            changed[start:end] = shuffled
            segments.append(
                {
                    "start_0": start,
                    "end_exclusive": end,
                    "seed": str(control_seed),
                    "unchanged": shuffled == original,
                }
            )
            start = end
        original_record = next(read_fasta(ROOT / entry["fasta"]))
        target = output / f"{entry['label']}.fasta"
        write_fasta(
            target, [Sequence(original_record.id, original_record.description, "".join(changed))]
        )
        assert Counter(source.bases) == Counter(changed)
        assert all(
            changed[i] == base
            for i, base in enumerate(source.bases)
            if base not in "ACGT" or protected[i]
        )
        genomic.append(
            {
                "label": entry["label"],
                "fasta": str(target.relative_to(ROOT)),
                "sha256": sha256_file(target),
                "original_sha256": source.sha256,
                "arrays": [asdict(a) for a in arrays],
                "segments": segments,
            }
        )
    manifest = {
        "namespace": NAMESPACE,
        "master_seed": MASTER_SEED,
        "created_at": datetime.now(UTC).isoformat(),
        "generator_sha256": sha256_file(Path(__file__)),
        "protocol_sha256": sha256_file(ROUND / "protocol.json"),
        "donor_pool": pool,
        "donor_pool_sha256": sha256(json.dumps(pool, sort_keys=True).encode()).hexdigest(),
        "repeat_controls": repeats,
        "genomic_controls": genomic,
    }
    atomic_json(output / "manifest.json", manifest)
    print(
        json.dumps(
            {
                "genomic_controls": len(genomic),
                "repeat_controls": len(repeats),
                "arrays": len(pool),
                "manifest_sha256": sha256_file(output / "manifest.json"),
            }
        ),
        flush=True,
    )


def provenance(settings, entries, profile):
    paths = sorted((ROOT / "src/crisprtracrrna").glob("*.py")) + [Path(__file__)]
    executables = doctor(settings)["executables"]
    result = {
        "created_at": datetime.now(UTC).isoformat(),
        "profile": profile,
        "settings": settings.public(),
        "source_sha256": {str(p.relative_to(ROOT)): sha256_file(p) for p in paths},
        "protocol_sha256": sha256_file(ROUND / "protocol.json"),
        "plan_sha256": sha256_file(ROOT / "docs/COUNCIL_IMPROVEMENT_PLAN.md"),
        "inputs": [{**e, "fasta_sha256": sha256_file(ROOT / e["fasta"])} for e in entries],
        "model_sha256": {str(p): sha256_file(p) for p in settings.model_paths()},
        "native_sha256": {
            name: {"path": path, "sha256": sha256_file(Path(path))}
            for name, path in executables.items()
        },
        "upstream_lock_sha256": sha256_file(ROOT / "upstream-lock.json"),
        "upstream_revisions": {
            "identify": command_provenance(settings.identify_command, settings.environment),
            "casandra": command_provenance(settings.casandra_command, settings.environment),
        },
        "python": platform.python_version(),
        "platform": platform.platform(),
        "workers": 2,
        "threads_per_locus": 1,
    }

    pinned = ROOT / "native/fasta36-seeded-20260927/fasta36-council-seed20260927"
    if Path(settings.fasta36).resolve() == pinned.resolve():
        result["native_dependency_profile"] = "council_seed20260927_threads1"
        result["protocol_amendment_sha256"] = sha256_file(ROUND / "protocol_amendment_v2.json")
        result["native_execution_manifest_sha256"] = sha256_file(
            ROUND / "native_execution_manifest.json"
        )
        result["configuration_sha256"] = sha256_file(ROOT / "local.reproducible.toml")
    return result


def execution_settings(config, experimental, pinned=False):
    settings = load_settings(
        config or ROOT / "local.toml", {"threads": 1, "short_anti_repeat_rescue": experimental}
    )
    if pinned:
        expected = ROOT / "native/fasta36-seeded-20260927/fasta36-council-seed20260927"
        if Path(settings.fasta36).resolve() != expected.resolve() or settings.threads != 1:
            raise ValueError(
                "Amended panels require --config local.reproducible.toml and native -T1"
            )
    return settings


def panel(name, experimental=False, labels=None, genomic=False, config=None):
    entries = inputs()
    if genomic:
        entries = [
            {"label": e["label"], "fasta": e["fasta"]}
            for e in json_read(ROUND / "controls/manifest.json")["genomic_controls"]
        ]
    if labels:
        entries = [e for e in entries if e["label"] in labels]
    output = ROUND / name
    output.mkdir(exist_ok=False)
    settings = execution_settings(
        config, experimental, pinned=genomic or name.endswith("_pinned_native")
    )
    manifest = provenance(settings, entries, "experimental" if experimental else "default")
    atomic_json(output / "run_manifest.json", manifest)
    started = time.monotonic()

    def execute(entry):
        before = time.monotonic()
        try:
            report = run(ROOT / entry["fasta"], output / entry["label"], settings)
            return {
                "label": entry["label"],
                "status": "complete",
                "elapsed_seconds": time.monotonic() - before,
                "prediction_count": report["summary"]["prediction_count"],
            }
        except Exception as exc:
            return {
                "label": entry["label"],
                "status": "failed",
                "error": str(exc),
                "elapsed_seconds": time.monotonic() - before,
            }

    results = []
    with ThreadPoolExecutor(2) as pool, (output / "results.jsonl").open("w") as handle:
        for row in pool.map(execute, entries):
            results.append(row)
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            print(json.dumps(row), flush=True)
    manifest.update(
        {
            "completed_at": datetime.now(UTC).isoformat(),
            "batch_wall_seconds": time.monotonic() - started,
            "results": results,
        }
    )
    atomic_json(output / "run_manifest.json", manifest)
    summarize(output, biological=not genomic)


def p95(values):
    return sorted(values)[math.ceil(0.95 * len(values)) - 1] if values else None


def summarize(output, biological=True):
    truth = {e["label"]: e for e in json_read(ROOT / "benchmarks/literature/loci.json")}
    manifest = json_read(output / "run_manifest.json")
    rows, all_lengths, total_bases = [], [], 0
    for entry in manifest["inputs"]:
        label = entry["label"]
        bases = sum(len(s.bases) for s in read_fasta(ROOT / entry["fasta"]))
        total_bases += bases
        row = {
            "label": label,
            "status": "failed",
            "candidate_count": 0,
            "detected": False,
            "best_rank": None,
            "coverage": 0.0,
            "bases": bases,
        }
        if not biological or truth[label]["status"] != "exact_unique":
            row["detected"] = None
            row["coverage"] = None
        report_path = output / label / "report.json"
        if report_path.exists():
            report = json_read(report_path)
            assert report["status"] == "complete"
            assert all(
                r["source"]["file_sha256"] == entry["fasta_sha256"] for r in report["records"]
            )
            ranked = [(r["source"]["id"], p) for r in report["records"] for p in r["predictions"]]
            ranked.sort(key=lambda cp: (-cp[1]["ranking"]["score"], cp[1]["id"]))
            all_lengths.extend(len(p["sequence"]) for _, p in ranked)
            row.update(
                {
                    "status": "complete",
                    "candidate_count": len(ranked),
                    "elapsed_seconds": report["wall_seconds"],
                    "distinct_intervals": len(
                        {
                            (c, p["strand"], p["interval"]["start"], p["interval"]["end"])
                            for c, p in ranked
                        }
                    ),
                    "origins": dict(Counter(p["candidate_origin"] for _, p in ranked)),
                    "admission_branches": dict(Counter(p["admission_branch"] for _, p in ranked)),
                    "evidence_classes": dict(Counter(p["evidence_class"] for _, p in ranked)),
                    "accepted_arrays": sum(len(r["arrays"]) for r in report["records"]),
                    "short_proposals": sum(
                        len(r["short_rescue_proposals"]) for r in report["records"]
                    ),
                    "rejections": dict(
                        Counter(
                            reason
                            for r in report["records"]
                            for p in r["short_rescue_proposals"]
                            for reason in p["short_rescue_admission"]["rejection_reasons"]
                        )
                    ),
                }
            )
            if biological and truth[label]["status"] == "exact_unique":
                ref = truth[label]["matches"][0]
                matching = [
                    (i, p, localization_match(p, ref, c))
                    for i, (c, p) in enumerate(ranked, 1)
                    if localization_match(p, ref, c)
                ]
                if matching:
                    row["detected"] = True
                    rank, pred, overlap = matching[0]
                    length = ref["end"] - ref["start"] + 1

                    def iou(p, overlap, length=length):
                        return overlap / (length + p["interval"]["length"] - overlap)

                    row.update(
                        {
                            "best_rank": rank,
                            "coverage": max(o / length for _, _, o in matching),
                            "highest_ranked_match": {
                                "id": pred["id"],
                                "rank": rank,
                                "coverage": overlap / length,
                                "iou": iou(pred, overlap),
                                "genomic_start_error": pred["interval"]["start"] - ref["start"],
                                "genomic_end_error": pred["interval"]["end"] - ref["end"],
                                "branch": pred["admission_branch"],
                            },
                            "best_iou": max(iou(p, o) for _, p, o in matching),
                        }
                    )
        row["candidates_per_mb"] = row["candidate_count"] * 1e6 / bases
        rows.append(row)
    counts = [r["candidate_count"] for r in rows]
    ordered = sorted(counts)
    summary = {
        "inputs": len(rows),
        "complete": sum(r["status"] == "complete" for r in rows),
        "candidate_count": sum(counts),
        "bases": total_bases,
        "candidates_per_mb": sum(counts) * 1e6 / total_bases,
        "candidate_median": ordered[len(ordered) // 2],
        "candidate_p95": p95(counts),
        "candidate_max": max(counts),
        "length_p95": p95(all_lengths),
        "batch_wall_seconds": manifest["batch_wall_seconds"],
        "short_admitted": sum(r.get("admission_branches", {}).get("short_rescue", 0) for r in rows),
        "short_proposals": sum(r.get("short_proposals", 0) for r in rows),
        "rows": rows,
    }
    if biological:
        mapped = [r for r in rows if truth[r["label"]]["status"] == "exact_unique"]
        summary.update(
            {
                "mapped": len(mapped),
                "localized": sum(r["detected"] is True for r in mapped),
                "top_k": {
                    k: sum(r["best_rank"] is not None and r["best_rank"] <= k for r in mapped)
                    for k in (1, 3, 5, 10)
                },
                "mrr": sum(1 / r["best_rank"] if r["best_rank"] else 0 for r in mapped)
                / len(mapped),
                "coverage90": sum(r["coverage"] >= 0.9 for r in mapped),
            }
        )
    atomic_json(output / "metrics.json", summary)
    return summary


def conditional(panel_name, config=None):
    output = ROUND / f"conditional_{panel_name}"
    output.mkdir(exist_ok=False)
    control_path = ROUND / "controls/manifest.json"
    controls = json_read(control_path)
    by_locus = defaultdict(dict)
    for item in controls["repeat_controls"]:
        if item["panel"] == panel_name:
            by_locus[item["label"]][item["array_id"]] = item["control_repeat"]
    settings = execution_settings(config, True, pinned=True)
    entries = inputs()
    manifest = provenance(settings, entries, "conditional_" + panel_name)
    manifest.update(
        {
            "control_manifest_sha256": sha256_file(control_path),
            "scope": (
                "conditional component control with frozen inferred arrays; "
                "no biological-negative claim"
            ),
        }
    )
    atomic_json(output / "run_manifest.json", manifest)
    started = time.monotonic()

    def execute(entry):
        work = output / entry["label"]
        work.mkdir()
        before = time.monotonic()
        try:
            source, arrays = frozen_context(entry["label"])
            arrays = [replace(a, repeat=by_locus[entry["label"]][a.id]) for a in arrays]
            fasta = work / "input.fasta"
            write_fasta(fasta, [source])
            runner = Runner(work / "logs", settings.timeout, settings.environment)
            proposals = []
            standard = search_anti_repeats(source, fasta, arrays, work, settings, runner, proposals)
            candidates = make_candidates(source, proposals, settings, "short_candidate")
            by_query = defaultdict(list)
            for hit in scan_models(
                [Sequence(c.id, c.id, c.sequence) for c in candidates], work, settings, runner
            ):
                by_query[hit.query].append(hit)
            terms = terminators(candidates, work, settings, runner)
            duplexes = interactions(candidates, work, settings, runner, anchored=True)
            proposal_rows = []
            for c in candidates:
                admission = short_rescue_admission(
                    c, by_query[c.id], terms.get(c.id, []), duplexes[c.id], settings
                )
                record = candidate_evidence(
                    c,
                    source,
                    by_query[c.id],
                    terms.get(c.id, []),
                    duplexes[c.id],
                    {"cassettes": []},
                    settings,
                )
                record["short_rescue_admission"] = admission
                record["admission_branch"] = "short_rescue"
                record["evidence_class"] = evidence_class(
                    duplexes[c.id] is not None,
                    bool(admission["supporting_model_hits"]),
                    bool(admission["supporting_terminator_hits"]),
                )
                proposal_rows.append(record)
            row = {
                "label": entry["label"],
                "status": "complete",
                "arrays": len(arrays),
                "bases": len(source.bases),
                "standard_proposals_not_scored": len(standard),
                "short_proposals": len(candidates),
                "short_admitted": sum(
                    r["short_rescue_admission"]["accepted"] for r in proposal_rows
                ),
                "proposals": proposal_rows,
                "elapsed_seconds": time.monotonic() - before,
                "source_sha256": source.sha256,
                "control_repeats": {a.id: a.repeat for a in arrays},
            }
            atomic_json(work / "result.json", row)
            atomic_json(
                work / "manifest.json",
                {
                    "files": {
                        str(p.relative_to(work)): sha256_file(p)
                        for p in sorted(work.rglob("*"))
                        if p.is_file()
                    }
                },
            )
            return {k: v for k, v in row.items() if k != "proposals"}
        except Exception as exc:
            row = {"label": entry["label"], "status": "failed", "error": str(exc)}
            atomic_json(work / "failure.json", row)
            return row

    results = []
    with ThreadPoolExecutor(2) as pool, (output / "results.jsonl").open("w") as handle:
        for row in pool.map(execute, entries):
            results.append(row)
            handle.write(json.dumps(row) + "\n")
            handle.flush()
            print(json.dumps(row), flush=True)
    counts = sorted(r["short_admitted"] for r in results if r["status"] == "complete")
    manifest.update(
        {
            "completed_at": datetime.now(UTC).isoformat(),
            "batch_wall_seconds": time.monotonic() - started,
            "results": results,
            "complete": len(counts),
            "short_admitted": sum(counts),
            "short_median": counts[len(counts) // 2] if counts else None,
            "short_max": max(counts, default=None),
        }
    )
    atomic_json(output / "run_manifest.json", manifest)


def reconstruct_candidate(p):
    array = p["array"]
    ai = array["interval"]
    a = Array(
        array["id"],
        Interval.from_public(ai["start"], ai["end"]),
        array["repeat"],
        array["predicted_strand"],
        array["category"],
        array["score"],
    )
    ar = p["anti_repeat"]
    anti = AntiRepeat(
        a,
        array["repeat"],
        array["tested_repeat_strand"],
        Interval.from_public(ar["interval"]["start"], ar["interval"]["end"]),
        p["strand"],
        ar["identity"],
        ar["coverage"],
        ar["complementary_hypothesis"],
    )
    return Candidate(
        p["id"],
        anti,
        Interval.from_public(p["interval"]["start"], p["interval"]["end"]),
        p["sequence"],
        Interval.from_public(ar["local_interval"]["start"], ar["local_interval"]["end"]),
    )


def performance():
    output = ROUND / "performance"
    output.mkdir(exist_ok=False)
    settings = load_settings(ROOT / "local.toml", {"threads": 1})
    results = []
    for label in ("Mse", "Ssu", "Sth1A"):
        report = json_read(ROUND / "baseline/reports" / f"{label}.json")
        predictions = [
            p
            for r in report["records"]
            for p in r["predictions"]
            if p["candidate_origin"] == "anti_repeat"
        ]
        candidates = [reconstruct_candidate(p) for p in predictions]
        maps, times, calls = {}, {}, {}
        for size in (1, 128):
            work = output / f"{label}_{size}"
            work.mkdir()
            runner = Runner(work / "logs", settings.timeout, settings.environment)
            before = time.monotonic()
            maps[size] = interactions(
                candidates, work, replace(settings, intarna_batch_size=size), runner
            )
            times[size] = time.monotonic() - before
            commands = json_read(work / "logs/commands.json")
            calls[size] = len(commands)
            atomic_json(work / "interactions.json", maps[size])
        old_maps = {p["id"]: p["interaction"] for p in predictions}
        projected = {
            k: None
            if v is None
            else {
                key: value
                for key, value in v.items()
                if key not in {"target_paired_positions", "repeat_paired_positions"}
            }
            for k, v in maps[128].items()
        }
        row = {
            "label": label,
            "candidates": len(candidates),
            "whole_map_equal": maps[1] == maps[128],
            "baseline_equal": old_maps == projected,
            "wall_seconds": times,
            "calls": calls,
            "grouped_fraction_of_serial": times[128] / times[1],
        }
        results.append(row)
        print(json.dumps(row), flush=True)
    atomic_json(
        output / "results.json",
        {
            "settings": settings.public(),
            "results": results,
            "native_executable_sha256": sha256_file(
                Path(doctor(settings)["executables"]["IntaRNA"])
            ),
        },
    )


def reverse_stress():
    output = ROUND / "reverse_stress"
    output.mkdir(exist_ok=False)
    entry = next(e for e in inputs() if e["label"] == "Ain")
    source = next(read_fasta(ROOT / entry["fasta"]))
    fasta = output / "Ain_reverse.fasta"
    write_fasta(fasta, [Sequence(source.id, source.description, reverse_complement(source.bases))])
    settings = load_settings(ROOT / "local.toml", {"threads": 1})
    report = run(fasta, output / "run", settings)
    original = json_read(ROUND / "small_default/Ain/report.json")
    forward = {
        (p["interval"]["start"], p["interval"]["end"], p["strand"])
        for r in original["records"]
        for p in r["predictions"]
    }
    reverse = {
        (
            len(source.bases) - p["interval"]["end"] + 1,
            len(source.bases) - p["interval"]["start"] + 1,
            "-" if p["strand"] == "+" else "+",
        )
        for r in report["records"]
        for p in r["predictions"]
    }
    atomic_json(
        output / "comparison.json",
        {
            "forward_count": len(forward),
            "reverse_mapped_count": len(reverse),
            "shared_intervals": len(forward & reverse),
            "forward_only": sorted(forward - reverse),
            "reverse_only": sorted(reverse - forward),
            "note": (
                "Upstream array inference and native alignment ordering may be asymmetric; "
                "source reconstruction remains mandatory."
            ),
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=["prepare", "panel", "conditional", "performance", "reverse"]
    )
    parser.add_argument("--name")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--experimental", action="store_true")
    parser.add_argument("--genomic", action="store_true")
    parser.add_argument("--labels", nargs="+")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "panel":
        panel(args.name, args.experimental, args.labels, args.genomic, args.config)
    elif args.action == "conditional":
        conditional(args.name, args.config)
    elif args.action == "performance":
        performance()
    else:
        reverse_stress()


if __name__ == "__main__":
    main()
