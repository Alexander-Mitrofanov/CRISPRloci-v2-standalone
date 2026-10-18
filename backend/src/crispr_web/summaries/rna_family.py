"""Bounded RNA-family preview preserving alignment and uncertainty."""
from .common import _json, _safe_file

WARNINGS = ["Homology, biological orientation and RNA boundaries are supplied by the user, not inferred.",
            "Exact duplicates contribute one sequence weight; related nonidentical sequences are not phylogenetically reweighted.",
            "Comparative folding scores include covariance contributions. Descriptive covariation is not a significance test or activity evidence."]


def public_report(doc):
    return {"schema": {"name": "crisprloci.rna-family", "version": "1.0.0"},
            "policy": {key: doc["policy"][key] for key in ("homology", "orientation", "boundaries", "deduplication", "weights",
                "temperature_celsius", "energy_parameters", "covariance_factor", "noncanonical_factor", "probability_threshold",
                "pair_coordinates", "covariation", "activity")},
            "engine": {key: doc["engine"][key] for key in ("name", "version")},
            "aligner": {key: doc["aligner"][key] for key in ("name", "version", "strategy", "threads")},
            "molecule": doc["input"]["molecule"],
            "counts": {key: doc["counts"][key] for key in ("sources", "unique_sequences", "alignment_columns", "mfe_pairs", "retained_probabilities")},
            "sources": [{key: row[key] for key in ("id", "length", "sha256", "rna_sha256", "representative_id")} for row in doc["sources"]],
            "alignment": [{key: row[key] for key in ("id", "member_ids", "sequence")} for row in doc["alignment"]],
            "columns": [{"column": row["column"], "counts": {key: row["counts"][key] for key in "ACGU-"}} for row in doc["columns"]],
            "fold": {key: doc["fold"][key] for key in ("consensus", "structure", "mfe_score_kcal_mol", "ensemble_free_energy_kcal_mol", "temperature_celsius")},
            "pairs": [{key: row[key] for key in ("left", "right", "probability", "canonical_sequences", "noncanonical_sequences",
                        "gapped_sequences", "pair_types", "compensatory_sequence_pairs")} for row in doc["pairs"]],
            "probabilities": [{key: row[key] for key in ("left", "right", "probability")} for row in doc["probabilities"]],
            "warnings": WARNINGS}


def rna_family_summary(root, *, max_rows, checkpoint=None):
    path = root / "rna-family/report.json"
    if _safe_file(root, path) is None:
        return {"available": False, "alignment": [], "pairs": [], "warnings": []}
    doc = public_report(_json(root, path, max_bytes=20_000_000, checkpoint=checkpoint))
    doc.pop("probabilities")
    doc.pop("columns")
    limit = max(0, min(max_rows, 500))
    doc["pairs_truncated"] = len(doc["pairs"]) > limit
    doc["pairs"] = doc["pairs"][:limit]
    doc["available"] = True
    return doc
