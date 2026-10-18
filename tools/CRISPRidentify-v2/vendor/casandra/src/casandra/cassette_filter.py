"""High-precision genome-context evidence gate for candidate Cas cassettes."""

from __future__ import annotations


def passes_genome_evidence_gate(cassette, genes_by_id, predictions_by_id):
    final = cassette["final_classification"]
    if final["subtype"] is None and final["type"] != "VII":
        return False, "no_supported_subtype"
    count = int(cassette["cas_gene_count"])
    cas_class = str(final["class"])
    if cas_class == "1":
        return (count >= 2, "class1_requires_two_cas_genes")
    if cas_class == "2":
        if count >= 2:
            return True, "class2_multi_gene"
        protein_id = cassette["cas_protein_ids"][0]
        length = int(genes_by_id[protein_id]["protein_length"])
        score = float(predictions_by_id[protein_id]["positive_profile_score"])
        accepted = length >= 400 and score >= 50.0
        return accepted, "class2_single_effector_requires_length400_score50"
    return (count >= 2, "unresolved_class_requires_two_cas_genes")
