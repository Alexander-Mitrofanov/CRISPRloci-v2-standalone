# Scientific design basis

The research council reviewed these primary sources before the owner approved
implementation. This record motivates the target design; it is not evidence of
accuracy for the incomplete v2 predictor. The full council plan remains in the
parent project's internal `improvenments/CRISPRleader-v2` folder.

| Primary source | Design constraint retained |
| --- | --- |
| Alkhnbashi et al. 2016, [Characterizing leader sequences of CRISPR loci](https://doi.org/10.1093/bioinformatics/btw454) | Distinguish reference-supported conserved cores from full functional leader boundaries. |
| Buyukyoruk et al. 2025, [Discovery of Diverse CRISPR Leader Motifs](https://doi.org/10.1089/crispr.2024.0093) | Motif interpretation has taxon/subtype scope; fixed upstream windows are not boundary ground truth. |
| Nuñez et al. 2016, [CRISPR Immunological Memory Requires a Host Factor for Specificity](https://doi.org/10.1016/j.molcel.2016.04.027) | Preserve system-specific experimental scope. |
| Santiago-Frangos et al. 2021, [Distribution and phasing of sequence motifs that facilitate CRISPR adaptation](https://doi.org/10.1016/j.cub.2021.05.068) | Preserve positional context without asserting a universal motif requirement. |
| Kim et al. 2019, [CRISPR DNA elements controlling site-specific spacer integration and proper repeat length](https://doi.org/10.1093/nar/gkz677) | A junction model needs the observed terminal repeat, not a consensus substitution. |
| Van Orden et al. 2020, [CRISPR type II-A subgroups exhibit phylogenetically distinct mechanisms](https://doi.org/10.1074/jbc.RA120.013554) | Subtype alone does not establish a motif's applicability. |
| Zhang et al. 2013, [Processing-independent CRISPR RNAs limit natural transformation](https://doi.org/10.1016/j.molcel.2013.05.001) | Upstream promoter evidence cannot be a universal leader requirement. |
| Wheeler and Eddy 2013, [nhmmer](https://doi.org/10.1093/bioinformatics/btt403) | A DNA profile route requires evaluation; it is not implemented here. |
| Grant, Bailey and Noble 2011, [FIMO](https://doi.org/10.1093/bioinformatics/btr064) | Declare the PWM null and complete testing scope; numerical significance is not biological function. |

The original CRISPRleader paper reports 195 reference clusters. The audited
legacy repository snapshot contained 193 HMM/FASTA pairs. That discrepancy and
undocumented repository-wide reuse terms preclude silently treating the legacy
inventory as a complete new reference release. No legacy models, matrices,
biological tables or source code were copied into this package.

All model curation, biological implementation, calibration and independent
benchmark gates listed in STATUS.md remain unfinished. The generic numerical
scanner is independently authored and tested on synthetic matrices only.
