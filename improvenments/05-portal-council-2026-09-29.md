# CRISPRloci portal methods council — 29 September 2026

The subsequent [connected scientific pipelines plan](06-connected-scientific-pipelines.md)
refines these recommendations around complete biological investigations and
reusable intermediate analyses, following the project owner's clarification.

The proposed direction is **one evidence-linked investigation per locus**:
identify its components, interpret its architecture, examine its targets and
mobile-element context, and compare it with related loci. Existing methods
provide a substantial foundation. The largest gains now come from connecting
their evidence and adding a few methods that answer distinct biological questions.

This is a planning assessment by three specialist agents and a coordinating
agent, with a second critique round. It is not an external expert review or a
biological benchmark. It revisits the September 27 roadmap against newer code;
earlier “future work” lists must not be treated as current implementation status.

## Inspection and current functionality

Inspection covered API request models and service capabilities, Integration
grouping, scientific tool entry points, summaries, result components, reference
search, pinned tool identities and scientific contracts. The public
[portal](https://alexander-mitrofanov.github.io/crispr-evolution-web/) loaded,
reported **Service ready**, advertised the methods below, and displayed its
locus, both spacer-search and database forms. The visible request cap was
2,000,000 bases and retention was 72 hours after completion.

No fresh scientific jobs were submitted. This confirms source capabilities and
read-only interface availability, not end-to-end execution of all deployed
methods or identity equality between the deployed release and this workspace.

| Area | Implemented in the inspected source | Important boundary |
| --- | --- | --- |
| Locus annotation | CRISPRidentify arrays, CasAndra Cas/cassette classification, CRISPRtracrRNA candidates, CRISPRleader adjacent windows | “Complete locus” currently combines these four tools, not every downstream analysis. |
| Array detection and strand | Categorized arrays and bundled CRISPRstrand inference | Strand prediction already exists; it is distinct from evOr acquisition-order evidence. |
| Cas and tracrRNA | Standalone Cas mode; Type II or V-K tracrRNA model screening | Candidates and model labels do not establish activity or full transcript boundaries. |
| Leaders | Both source-verified adjacent context windows | No biological leader-boundary prediction. |
| Repeats | Supplied DNA/RNA folds and per-instance isolated/local array-context folding | No calibrated repeat-family/subtype classifier or activity inference. |
| Evolution | SpacerPlacer reconstruction; evOr input/reverse comparison and histories | Exact orientation-neutral consensus-repeat grouping plus transitive spacer overlap; this is not a general locus-orthology model. |
| Spacer searches | Target DNA against the fixed 2017 spacer collection; supplied spacers against viral RefSeq with 0–2 substitutions | Exact first direction; full-length ungapped BLAST filtering in the second. No web PAM or array-overlap assessment; heuristic retrieval can miss matches. |
| Reference exploration | CRISPRcasDb-derived accession, repeat, spacer and Cas searches with occurrence details | Current code also supports same-length, both-orientation sequence neighborhoods with at most one literal substitution. This is not broad homology search. |
| Presentation and reproducibility | Basic array/Cas/tracrRNA genomic context strip, evolution evidence/history views, bounded previews, downloads, provenance, recovery/cancellation | Views remain largely mode-specific; context strips are not yet a comparative locus browser. |

Source anchors: [mode definitions](../frontend/src/science.js),
[submission contract](../backend/src/crispr_web/models.py),
[mode-specific result tabs](../frontend/src/components/results/ResultsView.vue),
[existing genomic strip](../frontend/src/components/results/annotations/AnnotationResults.vue),
[grouping](../Integration/crispr_integration/adapter.py),
[catalogue similarity](../backend/src/crispr_web/catalog_similarity.py),
[CRISPRstrand implementation](../tools/CRISPRidentify-v2/components/portable_strand.py),
and [eight-tool lock](../tools/tools.lock.json).

## Recommended additions

### 1. Source-bound self-targeting and direct array-to-target analysis

**Question:** which spacers match the supplied genome or associated replicons,
and which have external viral matches?

The pinned CRISPRspacer already has `prepare`, `self` and supplied-target
`search` commands. Its self workflow excludes matches overlapping any reported
array, retains exclusions, and supports explicitly declared assembly membership
and circular traversal. Its general search also retains target flanks. These
capabilities are absent from the public job modes inspected here.

Expose this existing science through independent web validation, then let users
send detected spacer occurrences into the existing viral workflow while
preserving array ID, ordinal, sequence and source identity. Default to
same-record candidates until explicit assembly metadata exists. Multi-FASTA
membership alone must not imply one organism; distinguish chromosome/plasmid
labels and their provenance when available. Undetected arrays remain an
overlap uncertainty, and matches are not evidence of active autoimmunity.

**Decision:** first scientific integration. This reuses a pinned tool but still
needs adapters, output validation, publication, UI and release gates.
[Existing tool contract](../tools/CRISPRspacer/README.md),
[CLI implementation](../tools/CRISPRspacer/src/crisprspacer/cli.py).

### 2. Repeat classification and observed Cas architecture

**Question:** what kind of locus is this, which machinery is observed, and do
the repeat and Cas evidence agree?

Evaluate standalone **RepeatTyper** for repeat-derived subtype evidence and
**CRISPRmap** for sequence/structure-family relationships. They answer different
questions. Orphan arrays and disagreement with nearby Cas classifications are
particularly valuable cases. Retain reference/model identity, raw scores and
unsupported/uncertain outcomes. Do not treat the classifier score as a calibrated
portal-wide confidence. [RepeatTyper documentation](https://github.com/Russel88/CRISPRCasTyper),
[CRISPRmap paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3783184/).

Add an interpretation layer over CasAndra: observed adaptation, processing and
interference components relative to a **versioned architecture**; partial genes,
contig edges, competing cassette associations and derived configurations.
CCTyper's component reporting is a useful comparator, not a reason to replace
CasAndra or run a second complete detector by default.

The 2025 classification recognizes seven types and 46 subtypes and includes
noncanonical architectures. A missing canonical component need not mean a
defective system. Preserve native installed-model labels and separately version
any newer vocabulary mapping; the publication does not establish installed-model
coverage. [Updated classification](https://www.nature.com/articles/s41564-025-02180-8).

**Decision:** early priority, after an independent family/clade-held-out
evaluation. Report observed architecture rather than a “functional immunity” score.

### 3. Protospacer context and qualified PAM evidence

**Question:** what lies next to a match, and is any known targeting rule
applicable to this particular effector?

Start with inspectable target flanks, coordinates, strand transformations,
alignment differences and contig-edge clipping. Where justified, attach an
experimentally characterized effector's motif rule with its source, orientation
convention and applicability limits. Subtype or protein similarity alone must
not silently transfer that rule. With uncertain array direction or array–Cas
assignment, retain alternatives and unknown compatibility.

Evaluate **Spacer2PAM** later for candidate motif inference from sufficiently
diverse, independently supported matches. Its documented workflow uses remote
BLAST/Entrez and optional PHASTER; it needs adaptation to the portal's offline,
pinned-reference worker design. Motif enrichment is not an experimentally
established interference requirement. [Official implementation](https://github.com/grybnicky/Spacer2PAM),
[method paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8990532/).

**Decision:** target-flank evidence belongs in the first target integration;
known-rule checks are conditional; inferred motifs need a later scientific gate.
Do not apply generic DNA-PAM logic to Type III/VI. Relevant Type III evidence
includes transcription and PFS/anti-tag relationships.
[Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6664175/).

### 4. Observed array comparison and locus synteny

**Question:** which loci are comparable, how do their observed spacer repertoires
differ, and is the surrounding Cas/gene arrangement conserved?

Derive and expose a spacer-sharing graph from the existing validated arrays,
including isolated arrays and the shared spacers supporting each edge.
The current adapter publishes overlap groups, not an explicit graph.
Transitive groups can connect A–B–C
without A sharing a spacer with C; show this rather than hiding it in a group ID.
Add aligned observed arrays with consistent spacer colors and explicit exact
versus clustered identity policies, linked to the existing history explorer.

**CCTK network/CRISPRdiff** provides useful methods and interoperable formats.
It accepts arrays from other detectors, so no additional detection stage is
necessary. Prefer native accessible rendering of validated portal data, with
CCTK export/benchmark interoperability; do not claim to execute CCTK when using
an independently implemented descriptive view.
[CCTK paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10457644/),
[official repository](https://github.com/Alan-Collins/CRISPR_comparison_toolkit).

Evaluate **clinker** for protein/gene-cluster comparisons around candidate loci.
It accepts annotated genomic inputs and would require source-bound CDS
translations and locus windows. Render validated alignment relationships in the
portal instead of serving tool-generated HTML. Shared spacers or gene similarity
alone do not prove orthology, functional equivalence or horizontal transfer.
[clinker repository](https://github.com/gamcil/clinker).

Later, accept an independent tree with an explicit sample/array mapping and
compare histories using the existing tools where supported or **CCTK Constrain**
as a separate sensitivity analysis. Constrain requires fully resolved trees;
reject unsupported polytomies rather than silently producing ancestry.
[Constrain contract](https://crispr-comparison-toolkit.readthedocs.io/en/latest/content/constrain.html).

**Decision:** high-value second delivery group. Broader comparative relationships
must remain separate from current validated evolutionary eligibility rules.

### 5. Prophage, plasmid and other mobile-element context

**Question:** does a same-source target fall in a predicted integrated virus,
or is the CRISPR locus itself associated with a mobile element?

**geNomad** is a strong candidate because it supplies virus/plasmid predictions,
provirus intervals and gene markers. This could join self-match and locus tracks
with relevant genomic context. Whole-contig/assembly analysis may be necessary;
a tiny locus window is not equivalent input. [Official workflow and outputs](https://portal.nersc.gov/genomad/quickstart.html).

Its current license states accredited-academic, internal, noncommercial research
restrictions. Public-service applicability must be resolved before adoption;
this assessment does not establish permission. [License](https://github.com/apcamargo/genomad/blob/main/LICENSE).
**VIBRANT** is an alternative to benchmark for viruses/proviruses, with a GPL-3+
license and an older documented runtime/model; it does not replace the plasmid
branch. Importing externally generated, source-bound MGE annotations is also
useful. [VIBRANT repository](https://github.com/AnantharamanLab/VIBRANT).

**Decision:** second delivery group, contingent on licensing, lineage-appropriate
region/boundary benchmarks and measured resources. A predicted prophage is not
proof of active infection.

### 6. Biological leader prediction and specialist mechanistic depth

**Question:** where is the conserved leader, and what evidence supports locus
direction and specialized functions?

The current CRISPRleader window extractor leaves a central biological gap.
Start curated-reference and boundary-benchmark work now. The historical
CRISPRleader is a useful baseline but its Python 2.7 implementation is not a
ready modern replacement. Require family support, conserved-core intervals,
clipping, divergent/unsupported cases and independent validation before enabling
leader claims. [Historical tool](https://github.com/BackofenLab/CRISPRleader).

Display detector strand, leader evidence and evOr acquisition-order evidence
alongside one another. They answer related but distinct questions; disagreement
should be retained, not resolved by an invented aggregate confidence.

A later specialist extension could interpret **Type III accessory/signalling
proteins** using curated profiles, after checking overlap with CasAndra.
Published research code is a starting point rather than a production adapter.
[Primary study](https://academic.oup.com/nar/article/52/12/7129/7684601),
[authors' code](https://github.com/vihoikka/hoikkala_etal_typeIII_effectors).

**Decision:** high research priority for leaders, evidence-gated delivery;
targeted accessory depth later.

## Council discussion and resolved order

| Perspective | Challenge | Resolution |
| --- | --- | --- |
| Locus biology | Leaders are central; postponing delivery must not make flank extraction the permanent endpoint. | Begin leader reference/benchmark research alongside early integrations; enable prediction only after validation. |
| Target ecology | Comparative work needs cohorts; flank evidence immediately improves existing self/viral matches. | Move transparent target context into the first target integration; separate it from motif inference. |
| Comparative genomics | Another visualization package does not itself produce a joined portal. | Extend existing native views and contracts; use external methods only where they add distinct analysis. |
| All reviewers | Proximity, score agreement and missing canonical genes can produce misleading biological claims. | Preserve conflicting evidence, explicit array–cassette ambiguity and versioned architecture expectations. |
| Coordinator | Existing reports and some READMEs describe older five-tool/exact-search states. | Base inventory on current code and lock; distinguish standalone tool support, web exposure and deployed execution verification. |

Proposed delivery order:

1. **Shared locus evidence model:** stable source/array/spacer/Cas IDs; explicit
   sample/assembly/topology; selected-locus zoom and detail; consistent exports.
   Preserve the existing genomic strip and history explorer as starting points.
   Annotated-genome/GenBank and protein input are useful enablers, not separate
   biological methods; define translation and coordinate provenance first.
2. **Early analysis depth:** existing self workflow, detected-spacer target handoff,
   target flanks, repeat typing and observed component architecture.
3. **Comparative and ecological context:** observed arrays, Cas-neighborhood
   synteny, MGE annotation, then independently supplied tree comparisons.
4. **Evidence-gated inference:** validated leaders, candidate PAM discovery and
   focused accessory-system analysis. Leader research starts during step 1.

The narrative is: **What is present? → What architecture is supported? → What
does its spacer repertoire match? → What genomic context explains those matches?
→ How does this locus differ and evolve across samples?**

## Scope boundaries and release requirements

Defer generic guide design, broad defense mining, de novo anti-CRISPR discovery,
and raw-read metagenomic reconstruction from the core portal roadmap. Targeted
Acr/Aca context may later help specific self-match/provirus investigations.
**iPHoP** fits a future virus-centered host-association workflow, but its CRISPR
component must not be counted as independent corroboration of the same spacer
hits. [iPHoP repository](https://github.com/simroux/iphop).
Optional expression/small-RNA tracks could add experimental context later;
sequence annotation and folding cannot establish activity.

Validate new methods with biological positives, hard negatives, family/clade
holdouts, reverse complements, contig edges, circular origins, missed arrays,
orphan/derived systems, conflicting assignments and duplicate/correlated target
hits. Measure abstention and reference coverage as well as errors. Freeze
evaluation datasets and acceptance criteria before final benchmarking; measure
RAM, storage and runtime rather than promising unsupported capacity.

Retain isolated immutable scientific releases, operator-managed reference
updates, independent Integration/backend validation, bounded publication,
cancellation, retention and per-job capabilities. A richer workspace must not
silently create indefinite storage or weaken recovery-link authorization.
Every implemented contract change needs appropriate `scripts/check.sh` coverage,
explicit browser-flow checks and new real-science release cases. Do not confuse
software contract tests with biological validation.

This task changed only this planning document. It did not alter application or
scientific code, dependencies, references, runtimes or services. No deployment,
new analysis jobs or application regression tests were performed. Local document
links and the private-data scanner were checked. The workspace root does not
expose a usable Git repository, so no root commit/diff or hook installation was
performed; existing independent tool checkouts were not modified.
