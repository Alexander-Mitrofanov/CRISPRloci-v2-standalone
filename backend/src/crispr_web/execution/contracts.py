"""Shared worker execution identifiers."""

MODE_TO_TOOL = {
    "locus_compare": "locus-compare",
    "prophage_context": "prophage-context",
    "insertion_context": "insertion-context",
    "defense_context": "defense-context",
    "rna_family": "rna-family",
    "array_compare": "array-compare",
    "repeat_map": "repeat-map",
    "repeat_type": "repeat-type",
    "spacer_association": "spacer-association",
    "viral_search": "viral-search",
    "protospacer": "protospacer",
    "repeats": "repeats",
    "repeat_context": "repeat-context",
    "detection": "identify",
    "cas": "casandra",
    "tracrrna": "tracrrna",
    "loci": "loci",
    "leader": "leader",
    "orientation": "crispr-evor",
    "reconstruction": "spacerplacer",
}

SUCCESS_STATUSES = frozenset({"completed", "completed_no_eligible_groups"})
