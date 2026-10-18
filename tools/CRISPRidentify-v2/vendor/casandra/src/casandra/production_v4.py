"""Fast CPU-only protein inference with a frozen reporting cutoff."""

from __future__ import annotations

import json
from pathlib import Path

from casandra.classification_v4 import normalize_classification
from casandra.decision import is_cas_hit
from casandra.production_v2 import load_model
from casandra.runtime import parse_fasta
from casandra.runtime_fast import scan_profiles_fast


MODEL_ID = "casandra-production-profile-hmm-type-v-v4-fast"


def _profile_families(metadata, profile_id):
    if profile_id is None:
        return []
    raw_families = metadata.get(profile_id, {}).get("families", [])
    families = []
    for value in raw_families:
        family = str(value).strip()
        if family and family not in families:
            families.append(family)
    return families


def predict_query_records_v4(model_dir: Path, queries: list[dict[str, str]], threads: int = 6):
    manifest, card, metadata, hmms = load_model(model_dir)
    report_evalue = float(card["report_evalue"])
    scans = scan_profiles_fast(queries, hmms, metadata, threads, report_evalue)
    threshold = float(card["detection_threshold"])
    records = []
    for query, scan in zip(queries, scans, strict=True):
        is_cas = is_cas_hit(scan, threshold)
        raw_classification = (
            dict(scan["classification"])
            if is_cas else
            {"class": None, "type": None, "subtype": None}
        )
        classification, hierarchy_conflict_masked = normalize_classification(
            raw_classification
        )
        profile_families = (
            _profile_families(metadata, scan["best_positive_profile"])
            if is_cas else []
        )
        records.append({
            "sequence_id": query["sequence_id"],
            "sequence_length": len(query["sequence"]),
            "is_cas": is_cas,
            "positive_profile_score": scan["positive_score"],
            "hard_negative_profile_score": scan["negative_score"],
            "score_margin": scan["margin"],
            "best_positive_profile": scan["best_positive_profile"],
            "profile_families": profile_families,
            "classification": classification,
            "variant": None,
            "evidence": {
                "model_id": manifest["model_id"],
                "profile_hits": scan["hit_count"],
                "decision_threshold": threshold,
                "positive_profile_required": True,
                "report_evalue": report_evalue,
                "type_v_subtype_masked": False,
                "raw_classification": raw_classification,
                "hierarchy_conflict_masked": hierarchy_conflict_masked,
            },
        })
    return records


def predict_proteins_v4(model_dir: Path, fasta_path: Path, output_path: Path, *, threads=6):
    queries = [{"sequence_id": name, "sequence": seq} for name, seq in parse_fasta(fasta_path)]
    records = predict_query_records_v4(model_dir, queries, threads)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for row in records:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    return {
        "model_id": MODEL_ID, "input_records": len(records),
        "cas_predictions": sum(x["is_cas"] for x in records),
        "output": str(output_path), "cpu_threads": threads,
    }

