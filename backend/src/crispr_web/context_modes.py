"""Release requirements and input bounds for connected context analyses."""

MODE_ARTIFACTS = {
    "locus_compare": ("clinker", "casandra"),
    "prophage_context": ("phigaro", "crispridentify", "crisprspacer"),
    "insertion_context": ("isescan", "casandra"),
    "defense_context": ("defensefinder", "casandra"),
    "rna_family": ("rnaalifold",),
}


def context_mode_available(mode, identities):
    return all(name in identities for name in MODE_ARTIFACTS.get(mode, ()))
