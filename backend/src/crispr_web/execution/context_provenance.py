"""Method citations for the independently provisioned context workflows."""

CONTEXT_CITATIONS = {
    "locus_compare": [{"tool": "clinker", "role": "Source-bound protein alignment and gene-neighborhood comparison",
                       "doi": "10.1093/bioinformatics/btab007", "url": "https://doi.org/10.1093/bioinformatics/btab007"}],
    "prophage_context": [{"tool": "Phigaro", "role": "Predicted prophage intervals joined to source-spacer sequence matches; activity not assessed",
                          "doi": "10.1093/bioinformatics/btaa250", "url": "https://doi.org/10.1093/bioinformatics/btaa250"},
                         {"tool": "CRISPRidentify", "role": "Source-bound CRISPR arrays supplying the spacer queries",
                          "doi": "10.1093/nar/gkaa1158", "url": "https://doi.org/10.1093/nar/gkaa1158"}],
    "insertion_context": [{"tool": "ISEScan", "role": "Insertion-sequence interval and terminal-repeat evidence; gene disruption not established",
                           "doi": "10.1093/bioinformatics/btx433", "url": "https://doi.org/10.1093/bioinformatics/btx433"}],
    "defense_context": [{"tool": "DefenseFinder", "role": "Defense-system and component predictions within the explicitly pinned model scope",
                         "doi": "10.1038/s41467-022-30269-9", "url": "https://doi.org/10.1038/s41467-022-30269-9"}],
    "rna_family": [{"tool": "RNAalifold", "role": "Comparative consensus folding of a user-justified homologous RNA alignment",
                    "doi": "10.1186/1471-2105-9-474", "url": "https://doi.org/10.1186/1471-2105-9-474"},
                   {"tool": "ViennaRNA", "role": "Unmodified comparative folding engine and pairing probabilities",
                    "doi": "10.1186/1748-7188-6-26", "url": "https://doi.org/10.1186/1748-7188-6-26"},
                   {"tool": "MAFFT", "role": "Multiple alignment of distinct RNA-family sequences",
                    "doi": "10.1093/molbev/mst010", "url": "https://doi.org/10.1093/molbev/mst010"}],
}
