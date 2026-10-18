"""Canonical Cas protein prediction API."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from casandra.decision import is_cas_hit
from casandra.model import ProteinModel, validate_threads
from casandra.runtime import parse_fasta
from casandra.runtime_fast import scan_profiles_fast


def _model(value: ProteinModel | Path) -> ProteinModel:
    return value if isinstance(value, ProteinModel) else ProteinModel.from_directory(value)


def profile_cas_family(
    metadata: dict[str, dict[str, object]], profile_id: object
) -> str | None:
    """Return the curated primary family of a raw HMM profile as a display label."""

    if not isinstance(profile_id, str) or not profile_id:
        return None
    profile = metadata.get(profile_id)
    families = profile.get("families") if isinstance(profile, dict) else None
    if not isinstance(families, list) or not families:
        return None
    primary = families[0]
    if not isinstance(primary, str) or not primary.strip():
        return None
    normalized = primary.strip()
    if any(character.isspace() or ord(character) < 32 for character in normalized):
        return None
    return normalized[0].upper() + normalized[1:]


def predict_query_records(
    model: ProteinModel | Path,
    queries: Iterable[dict[str, str]],
    *,
    threads: int = 6,
) -> list[dict[str, object]]:
    """Predict already-loaded amino-acid records without writing intermediate files."""
    loaded = _model(model)
    threads = validate_threads(threads)
    values = list(queries)
    report_evalue = float(loaded.card["report_evalue"])
    overlay = getattr(loaded, 'specificity', None)
    scans = scan_profiles_fast(
        values,
        list(loaded.hmms),
        loaded.metadata,
        threads,
        report_evalue,
        **({'capture_evidence': True} if overlay is not None else {}),
    )
    threshold = float(loaded.card["detection_threshold"])
    records: list[dict[str, object]] = []
    for query, scan in zip(values, scans, strict=True):
        is_cas = is_cas_hit(scan, threshold)
        cas_family = (
            profile_cas_family(loaded.metadata, scan["best_positive_profile"])
            if is_cas
            else None
        )
        records.append(
            {
                "sequence_id": query["sequence_id"],
                "sequence_length": len(query["sequence"]),
                "is_cas": is_cas,
                "positive_profile_score": scan["positive_score"],
                "hard_negative_profile_score": scan["negative_score"],
                "score_margin": scan["margin"],
                "best_positive_profile": scan["best_positive_profile"],
                "cas_family": cas_family,
                "result": cas_family if is_cas else "no cas",
                "classification": (
                    scan["classification"]
                    if is_cas
                    else {"class": None, "type": None, "subtype": None}
                ),
                "variant": None,
                "evidence": {
                    "model_id": loaded.model_id,
                    "profile_hits": scan["hit_count"],
                    "decision_threshold": threshold,
                    "positive_profile_required": True,
                    "report_evalue": report_evalue,
                    "type_v_subtype_masked": False,
                },
            }
        )
    if overlay is not None:
        return overlay.apply_records(values, records, scans, threads=threads)
    return records


def predict_proteins(
    model: ProteinModel | Path,
    fasta_path: Path,
    output_path: Path,
    *,
    threads: int = 6,
) -> dict[str, object]:
    """Classify an amino-acid FASTA and write versioned JSON Lines output."""
    loaded = _model(model)
    queries = [
        {"sequence_id": name, "sequence": sequence}
        for name, sequence in parse_fasta(Path(fasta_path))
    ]
    records = predict_query_records(loaded, queries, threads=threads)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for row in records:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    return {
        "schema_version": 1,
        "model_id": loaded.model_id,
        "input_records": len(records),
        "cas_predictions": sum(bool(row["is_cas"]) for row in records),
        "output": str(output_path),
        "cpu_threads": validate_threads(threads),
    }
