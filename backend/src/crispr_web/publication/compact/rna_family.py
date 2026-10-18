"""Complete family evidence in the module-result exchange envelope."""


def rna_family_data(read):
    doc = read.json("rna-family.json")
    return {key: doc[key] for key in ("molecule", "policy", "engine", "aligner", "counts", "sources", "alignment", "columns", "fold", "pairs", "probabilities", "warnings")}
