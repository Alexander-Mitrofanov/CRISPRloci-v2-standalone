"""Sparse ordered-profile features for CPU-only cassette classification."""

from __future__ import annotations


def cassette_features(cassette, genes_by_id, predictions_by_id):
    features: dict[str, float] = {}
    ordered_ids = sorted(
        cassette["cas_protein_ids"], key=lambda value: int(genes_by_id[value]["start_1based"])
    )
    tokens = []
    for protein_id in ordered_ids:
        prediction = predictions_by_id[protein_id]
        profile = prediction["best_positive_profile"]
        if profile is None:
            continue
        tokens.append(str(profile))
        key = f"profile:{profile}"
        features[key] = features.get(key, 0.0) + 1.0
        for head in ("type", "subtype"):
            label = prediction["classification"].get(head)
            if label is not None:
                key = f"{head}:{label}:score"
                features[key] = features.get(key, 0.0) + max(
                    1.0, float(prediction["positive_profile_score"])
                ) / 100.0
    features["architecture:cas_gene_count"] = float(len(tokens))
    for left, right in zip(tokens, tokens[1:]):
        key = f"bigram:{left}>{right}"
        features[key] = features.get(key, 0.0) + 1.0
    return features


def dominant_cassette(cassettes):
    return max(
        cassettes,
        key=lambda row: (
            int(row["cas_gene_count"]),
            float(row["classification"]["confidence"]["subtype_vote_fraction"]),
        ),
        default=None,
    )
