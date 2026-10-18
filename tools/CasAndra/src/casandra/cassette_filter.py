"""Configurable high-precision genome-context evidence gate."""

from __future__ import annotations

from collections.abc import Mapping


DEFAULT_POLICY = {
    "supported_subtype_required": True,
    "class_1_minimum_cas_genes": 2,
    "class_2_minimum_cas_genes_or_single_effector": 2,
    "class_2_single_effector_minimum_aa": 400,
    "class_2_single_effector_minimum_positive_profile_score": 50.0,
}


def passes_genome_evidence_gate(
    cassette,
    genes_by_id,
    predictions_by_id,
    policy: Mapping[str, object] | None = None,
):
    configured = {**DEFAULT_POLICY, **dict(policy or {})}
    final = cassette["final_classification"]
    if configured["supported_subtype_required"] and final["subtype"] is None:
        return False, "no_supported_subtype"
    count = int(cassette["cas_gene_count"])
    cas_class = str(final["class"])
    if cas_class == "1":
        minimum = int(configured["class_1_minimum_cas_genes"])
        return count >= minimum, f"class1_requires_{minimum}_cas_genes"
    if cas_class == "2":
        minimum = int(configured["class_2_minimum_cas_genes_or_single_effector"])
        if count >= minimum:
            return True, "class2_multi_gene"
        protein_id = cassette["cas_protein_ids"][0]
        length = int(genes_by_id[protein_id]["protein_length"])
        score = float(predictions_by_id[protein_id]["positive_profile_score"])
        minimum_length = int(configured["class_2_single_effector_minimum_aa"])
        minimum_score = float(
            configured["class_2_single_effector_minimum_positive_profile_score"]
        )
        accepted = length >= minimum_length and score >= minimum_score
        return accepted, "class2_single_effector_thresholds"
    return count >= 2, "unresolved_class_requires_two_cas_genes"

