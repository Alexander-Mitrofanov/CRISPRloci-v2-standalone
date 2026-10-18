# Leader council evidence and review record

Targeted primary-literature and source audit, 27 September 2026. Computational
methods motivate the proposed design; no improvement in v2 accuracy has yet
been demonstrated. Experimental findings apply to their tested systems. This
is not an exhaustive systematic review.

## Scientific evidence

| Source | What it supports | Limit imposed on the plan |
| --- | --- | --- |
| Alkhnbashi et al. (2016), Bioinformatics 32:i576–i585, [Characterizing leader sequences of CRISPR loci](https://doi.org/10.1093/bioinformatics/btw454) | Repeat-related families and conservation/profile evidence for conserved cores. | A conserved core is not automatically the full functional leader or promoter; reference scope and boundary validation matter. |
| Buyukyoruk et al. (2025), The CRISPR Journal 8:137–148, [Discovery of Diverse CRISPR Leader Motifs…](https://doi.org/10.1089/crispr.2024.0093) | Taxon/subtype-associated motifs and putative biological interpretations. Verified journal publication, online 8 January 2025. | Fixed upstream extraction is not boundary ground truth; many functional interpretations are predictions and there is no universal leader motif. |
| Nuñez et al. (2016), Molecular Cell 62:824–833, [CRISPR Immunological Memory Requires a Host Factor for Specificity](https://doi.org/10.1016/j.molcel.2016.04.027) | IHF-mediated leader recognition and integration specificity in the studied E. coli system. | A mechanism-specific experiment does not establish a universal promoter or full-leader boundary rule. |
| Santiago-Frangos et al. (2021), Current Biology, [Distribution and phasing of sequence motifs that facilitate CRISPR adaptation](https://doi.org/10.1016/j.cub.2021.05.068) | Sequence position and phasing can matter; includes type I-F experimental evidence. | Preserve positional context; absence of a canonical IHF motif does not exclude a leader. |
| Kim et al. (2019), NAR 47, [CRISPR DNA elements controlling site-specific spacer integration and proper repeat length by a Type II CRISPR–Cas system](https://doi.org/10.1093/nar/gkz677) | Leader-proximal and leader–repeat junction evidence in a studied II-A system. | Include the observed terminal repeat when a model spans the junction; do not infer an upstream endpoint from that motif alone. |
| Van Orden et al. (2020), JBC, [CRISPR type II-A subgroups exhibit phylogenetically distinct mechanisms for prespacer insertion](https://doi.org/10.1074/jbc.RA120.013554) | Mechanistic variation even within a subtype. | Record organism/lineage scope rather than promoting every motif to all subtype members. |
| Zhang et al. (2013), Molecular Cell, [Processing-independent CRISPR RNAs limit natural transformation in Neisseria meningitidis](https://doi.org/10.1016/j.molcel.2013.05.001) | Repeat-embedded promoters in the examined systems. | Upstream promoter evidence cannot be a universal leader requirement. |
| Wheeler and Eddy (2013), Bioinformatics 29:2487–2489, [nhmmer: DNA homology search with profile HMMs](https://doi.org/10.1093/bioinformatics/btt403) | A DNA-specific profile-search implementation to evaluate. | Search significance/alignments alone do not validate leader boundaries; inspect alphabet and profile format. |
| Grant, Bailey and Noble (2011), Bioinformatics 27:1017–1018, [FIMO: scanning for occurrences of a given motif](https://doi.org/10.1093/bioinformatics/btr064) | Statistical PWM scanning baseline. | Define background and actual testing scope; a motif p/q-value is not a leader/function probability. |

Author documentation confirms [nhmmer's DNA scope](https://www.ebi.ac.uk/jdispatcher/blog/2024-07-23-nhmmer)
and [FIMO's reporting behavior](https://meme-suite.org/meme/doc/fimo.html?man_type=web).
The latter warns that capped reporting can affect completeness and q-value
calculation. Exact external versions, packaging terms and settings will be
selected and pinned during authorized implementation, not assumed from a
historical paper.

Additional experimental curation candidates include
[McGinn and Marraffini (2016)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5096952/),
[Rollie et al. (2018)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5815122/) and
[Kieper et al. (2019)](https://academic.oup.com/femsle/article/366/11/fnz129/5525085).
Their accession/coordinate mapping and assay-specific labels still require
curation. A mutant that reduces adaptation is a control for that assay, not
necessarily a proven negative for conserved-core or transcription prediction.

## Source inventories and reuse gates

| Source snapshot | Observed inventory | Consequence |
| --- | --- | --- |
| [BackofenLab/CRISPRleader](https://github.com/BackofenLab/CRISPRleader/tree/bde50f7939bbe400a4812e4f2959f9340b743bf5), commit `bde50f7939bbe400a4812e4f2959f9340b743bf5` | 193 HMMs and 193 cluster FASTAs: 51 archaeal, 142 bacterial. Python 2.7-era code/bundled programs. No repository-wide code/data license found; license paths identified for bundled MView do not cover the entire collection. | Compare against the paper's 195-cluster total (51+144). Reconcile, do not invent missing models. Do not assume downloadable code/models/binaries can be redistributed under the new project's license. |
| [WiedenheftLab/2024_Buyukyoruk_et_al_CRISPRjournal](https://github.com/WiedenheftLab/2024_Buyukyoruk_et_al_CRISPRjournal/tree/c1fba682bc971576afcda481c86c9d15eb526924), commit `c1fba682bc971576afcda481c86c9d15eb526924` | Motif/sequence/control resources and analysis scripts. No license file found in the audited tree. | Useful reproducibility reference; asset-specific reuse/redistribution basis remains to be established. Repository name does not determine publication status. |

These are audit observations, not a determination that reuse is prohibited.
Do not assume an article's license governs a separate source repository or all
third-party assets. The historical Freiburg download endpoint could not be
retrieved in the reference audit; this does not prove the data are unavailable
elsewhere. No authors were contacted and no permissions were obtained here.

Preferred fallback: independently implement published methods and build fresh
models from a fixed, manageable genome panel with reviewed source terms and
curated factual locus/evidence records. Record accession **versions**, sequence
hashes and known restrictions. [NCBI's policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)
does not automatically transfer third-party rights; do not label all downloaded
data CC0. Generate alignments/models independently and preserve the build recipe.
Do not copy undocumented HMMs/PWMs or publication figures/tables into the new
repository. A user-supplied reference option is supplementary, not the default
completion path.

## Reference bundle and benchmark ledger

The proposed manifest must record:

- Bundle schema/ID/version, file hashes, build-code commit, dependencies, source
  freeze date, taxonomy snapshot, parameters/seeds and attribution/rights ledger.
- Source accession.version, sequence digest, retrieval/source URL, genome
  release date, taxon, applicable terms and any exclusion.
- Locus ID, coordinates/strand, observed terminal repeat, repeat family,
  subtype/context evidence, truncation and boundary meaning: core, functional
  interval, promoter/TSS, junction or merely an extracted window.
- Evidence tier (experimental, reviewed comparative, computational candidate),
  DOI plus exact evidence locator, curator/reviewer, exact/interval/unknown label.
- Model family/type, training locus IDs, applicability, build parameters,
  calibration records, cutoff/policy version and supported output claims.
- Split assignment grouped by source genome, close homology and repeat/leader
  family, with exclusions and leakage audit covering preprocessing as well as
  final model fitting.

Separate engineering fixtures from biological benchmarks. Original training
clusters evaluated against their own supplied HMMs are regression examples,
not independent validation. Experimental labels can still overlap model
training sequences; audit both label independence and sequence leakage.

Measure false calls per tested flank/genome, detection precision/recall when
labels justify them, prediction coverage, abstention, known-strand accuracy,
boundary error and interval overlap. Stratify supported-family versus
unseen-family results, lineage and fragmentation. Report uncertainty with
genome/family-aware resampling where feasible. Boundary tolerance windows
must be predeclared and compare equivalent label meanings.

Permissible scope-specific acceptance targets will be fixed before final test
evaluation. The council has not asserted a universal 95% accuracy threshold,
full-family coverage or a universal exact leader boundary. Insufficient labels
require a narrower validated scope or explicit prerelease status.

## Code boundary checked in the current server

- [Identify's existing flank routine](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/components/components_non_array_computations.py:501)
  extracts fixed flanks; it is not an existing conserved-core predictor.
- [Canonical Identify output](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/components/report_model.py:1225)
  declares one-based inclusive source-forward coordinates with separate
  orientation and source digest.
- [The evolutionary adapter](/home/alex/CodexProjects/CRISPRloci-v2/Integration/crispr_integration/parser.py:366)
  does not retain the genomic interval in its array record. The new adapter
  therefore reads canonical input directly without importing server modules.

## Final council review

All three reviewers support the standalone plan subject to the owner confirming
it. Scientific review required explicit motif domains and observed-repeat
junction scanning. Reference review required an actual default biological
bundle and passing predeclared scientific gates, not merely a recipe/report.
These corrections are incorporated in PLAN.md. Architecture review establishes
the frozen-contract, separate-worker ownership and later server-release boundary,
and required explicit circular GFF3/BED exports and separation of broken-reference
setup failures from scientific lack of family coverage. Both were incorporated.

At the council stage only planning documents were created. The owner subsequently
confirmed the plan, and implementation began. This literature review does not
establish accuracy of the incomplete predictor or substitute for the independent
biological evaluation. See README.md for current implementation and blocker status.
