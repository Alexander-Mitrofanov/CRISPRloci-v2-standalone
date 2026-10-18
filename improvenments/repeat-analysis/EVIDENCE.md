# Repeat-analysis evidence and council record

Reviewed 27 September 2026. Targeted primary-source review and local source
inspection; no scientific execution or accuracy benchmark was performed.
Read the adopted [implementation decision](DECISION.md).

## Biological and historical evidence

- **B1 — Reimann et al., 2017**, [Structural constraints and enzymatic promiscuity
  in the Cas6-dependent generation of crRNAs](https://doi.org/10.1093/nar/gkw786).
  Specific Synechocystis substrates support investigating sequence-context
  effects. Fragment length can matter beyond immediately adjacent spacers.
  Exact tested substrates and study-specific settings are needed for
  reproduction; this paper does not validate a general processing classifier.
- **B2 — Alkhnbashi et al., 2021**, [CRISPRloci](https://doi.org/10.1093/nar/gkab456),
  with [official help](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=CRISPRloci).
  Establishes the isolated/array-context comparison and repeat-only historical
  scope. Its local window/span is a reproduction baseline, not an optimal
  parameter claim. The historical label "structure accuracy" is inappropriate
  as experimental accuracy for the proposed descriptive report.
- **B3 — Shao et al., 2016**, [A Non-Stem-Loop CRISPR RNA Is Processed by Dual
  Binding Cas6](https://doi.org/10.1016/j.str.2016.02.009).
  Protein-stabilized processing is a counterexample to making isolated hairpin
  formation a requirement for a functional repeat.
- **B4 — Deltcheva et al., 2011**, [CRISPR RNA maturation by trans-encoded small
  RNA and host factor RNase III](https://doi.org/10.1038/nature09886).
  Supports retaining tracrRNA interaction evidence separately from isolated
  repeat folding for applicable systems.
- **B5 — Swarts et al., 2017**, [Structural Basis for Guide RNA Processing and
  Seed-Dependent DNA Targeting by CRISPR-Cas12a](https://doi.org/10.1016/j.molcel.2017.03.016).
  Cas12a's protein-bound repeat structure and processing limit interpretation
  of generic pseudoknot-free RNA models. Do not extrapolate across all Cas12.

## Numerical methods

- **N1 — [RNAfold official manual](https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/man/RNAfold.html)**.
  MFE and ensemble outputs are different objects. PostScript `ubox` entries
  encode square roots of pair probabilities. Pin all effective settings.
- **N2 — [RNAplfold official manual](https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/man/RNAplfold.html)**.
  Local pairs are averaged over windows. Native unpaired-region probabilities,
  span/window and reporting cutoff require explicit handling. The live manual
  inspected describes 2.7.2; the repository pins 2.7.0, whose actual behavior
  must be verified during implementation.
- **N3 — Stephan Bernhart, 2007**, [Variations on RNA folding — Locally stable
  structures and RNA hybridization](https://www.tbi.univie.ac.at/downloads/team/degree/bernhart_diss.pdf),
  section 5.3.2, printed pp. 96–97. Explains why pair-specific window averages
  are not a normalized distribution of partners. The council's restriction on
  row sums, entropy and summed external probability follows this distinction.
- **N4 — Lorenz et al., 2011**, [ViennaRNA Package 2.0](https://doi.org/10.1186/1748-7188-6-26).
  Scientific basis for using established thermodynamic tools; does not validate
  the new integration or prove cellular folding from sequence alone.

## Reference and classifier candidates

- **R1 — Lange et al., 2013**, [CRISPRmap](https://doi.org/10.1093/nar/gkt606),
  [official documentation](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=CRISPRmap),
  [original repeat FASTA](https://rna.informatik.uni-freiburg.de/CRISPRmap/oridata/CRISPRmap_data.fa).
  Family, structure motif and Cas subtype are different classifications.
  Historical versions have different vocabularies. The original downloadable
  collection is v1 material, not automatically the later v2 reference. Some
  accessions lack version suffixes; origin reconstruction, model recovery and
  redistribution terms remain implementation gates.
- **R2 — Alkhnbashi et al., 2014**, [CRISPRstrand](https://doi.org/10.1093/bioinformatics/btu459).
  Published strand inference is scientific precedent, not validation of the
  currently bundled CNN. Obtain training provenance or explicitly limit claims
  about independence of an external test.
- **R3 — Russel et al., 2020**, [CRISPRCasTyper](https://doi.org/10.1089/crispr.2020.0059),
  [official repository / RepeatTyper](https://github.com/Russel88/CRISPRCasTyper).
  Separate repeat subtyping is available. Code is MIT licensed; the exact
  model/data provenance still needs review. Models changed over time, so pin
  weights and label vocabulary independently. The repository's suggested score
  cutoff is a candidate baseline, not this project's calibrated threshold.
- **R4 — Nethery et al., 2021**, [CRISPRclassify](https://doi.org/10.1089/crispr.2021.0021),
  [official repository](https://github.com/CRISPRlab/CRISPRclassify).
  `classifyRepeats()` is an appropriate separate implementation comparator;
  shared training/reference data still require an overlap audit. Package/data
  redistribution terms and exact model artifacts were not resolved here.

## Current source facts

These are observations of the pinned checkout, not tests of its runtime.

| ID | Source | Finding |
| --- | --- | --- |
| C1 | [report_model.py](../../tools/CRISPRidentify-v2/components/report_model.py), `_public_interval` and repeat construction | Native one-based inclusive repeat intervals and observed source sequences are available; canonical arrays retain source span, consensus and architecture. |
| C2 | [Integration parser](../../Integration/crispr_integration/parser.py), [backend identify validator](../../backend/src/crispr_web/execution/validation_identify.py) | Evolutionary `ArrayRecord` omits repeat geometry. Backend report binding does not independently reconstruct every repeat. |
| C3 | [portable_strand.py](../../tools/CRISPRidentify-v2/components/portable_strand.py), [model manifest](../../tools/CRISPRidentify-v2/trained_models/strand/crisprstrand_model_r_v1.json) | Checksummed CNN artifact, independent forward/reverse scores, 38-base truncation after orientation, zero encoding of noncanonical input and legacy ties choosing reverse. |
| C4 | [environment.yml](../../tools/CRISPRidentify-v2/environment.yml), [core lock](../../tools/CRISPRidentify-v2/locks/conda-linux-64-core.lock) | ViennaRNA 2.7.0 is pinned. Availability of both commands and exact options must be smoke-tested in the chosen release. |
| C5 | [compact_feature_extraction.py](../../tools/CRISPRidentify-v2/components/compact_feature_extraction.py), [verified repeat FASTA](../../tools/CRISPRidentify-v2/tools/blasting/Verified_repeats_dataset1.fa) | Detector MFE is an energy feature, not the P4 report. Anonymous `db1_N` repeat identifiers cannot establish genomic origins. |
| C6 | [runner](../../Integration/crispr_integration/runner.py), [publication](../../backend/src/crispr_web/publication/annotations.py), [public contract](../../frontend/src/contracts/public-api-v1.json) | Worker adapters and approved artifact publication provide insertion points. Existing public modes lack typed repeat-only input. |

## Council and scope

| Agent task | Responsibility | Adopted finding |
| --- | --- | --- |
| `repeat_biology` | Primary mechanisms, metric semantics, controls | Per-instance context evidence; explicit motif-pair support; no generalized processing/functional claims. |
| `repeat_architecture` | Read-only source/runtime/contract audit | Dedicated source reconstruction at both validation boundaries; focused worker adapter and bounded publication. |
| `repeat_reference_methods` | Families, strand/subtype methods, references and evaluation | Retrieval first; separate assignments, pinned training/reference provenance, selective predictions and family-aware splits. |
| Coordinator | Historical/numerical cross-check and final decision | Stage P4a–P4c; reject invalid RNAplfold row-sum interpretations; distinguish deliverable hypotheses from validated predictions. |

The architecture reviewer explicitly corrected an initial suggestion to sum
local pair values into outside-repeat probabilities and enforce global row-sum
bounds. All three reviews inform the decision; agent agreement is not evidence
of biological performance.

The biology and reference reviewers also checked the written decision. Their
final corrections clarify motif/model eligibility, incomparable pair sets,
shared comparator training data and abstention for absent family labels.

Only internal Markdown planning files were changed. No scientific code, tool
checkout, runtime, reference database, job or production service was changed.
Document links and the decision's cross-references are checked before delivery.
Application/browser/scientific tests are deferred to the implementation gates
because this task changes no executable code or UI.
