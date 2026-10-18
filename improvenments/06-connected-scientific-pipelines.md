# Connected scientific pipelines for CRISPRloci

The council recommends completing biological investigations by adding reusable
analysis steps between the existing scientific tools. The first deliverable
should connect a detected locus to its interpretable spacer targets. The second
should connect differences between comparable loci to their potential targets
and evolutionary histories. Comparative RNA and regulatory-region analysis can
then reuse the same sequence selection, alignment and evidence machinery.

This refines the [previous council assessment](05-portal-council-2026-09-29.md)
following the project owner's request to prioritize connected scientific
pipelines. Three agents reconsidered locus biology, target ecology and comparative
analysis, then critiqued the combined proposal. The pipelines below are proposed
portal designs. Cited studies establish precedents for analytical steps; they do
not validate our proposed combinations or measure how frequently all biologists
perform them. Reviewed 29 September 2026.

## Reusable scientific analyses to add

| Analysis step | What the researcher currently has to assemble | Reusable result and downstream consumers |
| --- | --- | --- |
| Sequence collection and exact deduplication | Extract repeats, spacers, proteins or flanks; retain identifiers after transformations; collapse duplicate search work | A sequence set with every source occurrence and transformation. Feeds target search, alignments, family comparisons and histories. |
| Candidate biological associations | Decide which array, cassette and RNA candidates could belong together | An array–cassette–RNA relationship table with evidence and alternatives. Feeds architecture interpretation, effector-specific target rules and comparative grouping. |
| Reference recruitment and locus comparison groups | Find related examples, distinguish separate loci within a genome and choose a defensible comparison set | A versioned selection with inclusion/exclusion reasons and separate repeat, Cas and synteny evidence. Feeds cohort comparisons and conserved-region analysis. |
| Multiple alignment and conservation analysis | Align selected sequences; inspect mutations, gaps, conserved blocks and logos; map differences back to genomes | Alignments, per-column coverage/conservation and source-coordinate differences. Feeds repeat variation, leader research, RNA analysis and protein comparison. |
| Genomic interval intersections | Join match coordinates to genes, arrays, prophages, replicons and target flanks | Feature-overlap and distance tables retaining strand and coordinate provenance. Feeds self-match interpretation and target reports. |
| Matrices and explicit aggregation | Convert annotation and match tables into comparable counts, presence matrices, distances and networks | Sample/locus × spacer, Cas-component and target matrices with units, denominators and assessed/unknown states. Feeds biological comparisons and selected downstream statistics. |
| Observed change and history reconciliation | Match spacer identities across arrays, tools, samples and inferred events | An occurrence/sequence/cluster/event crosswalk. Feeds target annotations on observed differences and inferred gains or losses. |

These are scientific analyses with explicit assumptions, not just new visual
components. Some are deterministic joins; others, especially locus grouping,
alignment, clustering and inference, require method-specific validation.
Researcher-selected subsets and exported intermediate tables should be usable
by another analysis without retyping identifiers or recreating provenance.

## Pipeline 1 Interpret the targets of a detected locus

**Question:** Which sequences does this array match, what genomic context
explains those matches, and which associations merit follow-up?

1. Accept genomic assemblies and explicit contig-to-assembly declarations. Use
   current array/Cas annotation. Record candidate array–cassette relationships;
   ambiguous relationships remain valid input for sequence matching.
   Keep proximity, repeat compatibility and architecture evidence separate;
   permit several partners. Cross-contig distance is unknown, and assembly
   membership alone does not establish a functional relationship.
2. Extract a spacer occurrence inventory: sample, contig, array, ordinal,
   source coordinates, sequence and orientation evidence. Search each exact
   sequence once while retaining all original occurrences.
3. Search the same source assembly, supplied target sequences and frozen public
   references through the existing CRISPRspacer capabilities. Preserve skipped
   queries, search bounds, incomplete retrieval and no accepted match outcomes.
4. Rejoin matches to all originating arrays. Normalize target accession/version
   and exact sequence identity. Retain duplicate accession relationships and
   supplied genome/segment/taxonomy mappings without turning each viral segment
   into an independent species or host association.
5. Intersect target intervals with source-bound annotations: arrays, genes,
   coding strand, declared replicons and available MGE regions. Preserve both
   flanks and any missing context at contig boundaries. Initial gene/MGE tracks
   can come from validated imports; new predictors need not block the pipeline.
6. Produce an array-by-target support matrix, per-target coordinate map and
   evidence table. Count distinct spacer sequences separately from their copies
   and matching sites. Show multi-target ambiguity, eligible/searched denominators
   and unavailable annotations.

**Completed product:** the researcher can select a candidate target, inspect
which arrays support it, see where those matches occur, identify contradictory
evidence and export the selected sequences and evidence for follow-up. This
closes the path from genome annotation to a target investigation.

The public viral workflow already aggregates distinct spacer support per
accession in [CRISPRspacer](../tools/CRISPRspacer/src/crisprspacer/viral_search.py).
Extend that existing analysis with array/sample identities, target relationships,
feature intersections and matrices. The standalone tool already supports self
and supplied-target search; their portal exposure is missing.
[Current commands](../tools/CRISPRspacer/src/crisprspacer/cli.py).

This connected investigation has precedent in CRISPRTarget's combination of
alignments, recognition-context inspection and genomic comparisons.
[CRISPRTarget paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3737339/).
Self-match interpretation similarly benefits from joining candidate targets to
coding regions, strand and prophage intervals, as in this
[self-targeting study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7708048/).

**Useful branches:** array overlaps prompt review of array evidence; unexplained
same-assembly matches can be inspected against MGE and architecture evidence;
well-supported external matches can enter conditional PAM analysis. A surviving
self-match does not establish active autoimmunity or explain it by an Acr.

## Pipeline 2 Compare loci and explain their differences

**Question:** How do corresponding loci differ across isolates, and which
candidate targets are associated with those differences?

1. Accept a sample manifest and local locus-quality evidence. Distinguish samples,
   assemblies and multiple loci per genome. Contig truncation and unresolved
   genomic regions must remain unknown observations.
2. Form candidate homologous-locus groups using repeat relationships, Cas
   architecture and conserved flanking genes/order. Keep all contributing
   evidence visible. A syntenic locus may be comparable after complete spacer
   turnover; shared spacers are not a prerequisite for every comparison.
3. Build observed spacer repertoires, arrays and Cas-component matrices. Align
   ordered arrays where correspondence is supported. Report unalignable
   repertoires instead of inventing spacer-column homology.
4. Calculate defined repertoire distances and sharing networks. Report which
   units and denominator each metric counts. Show direct shared-spacer edges,
   including isolated arrays, separately from transitive group membership.
5. For eligible groups, reuse SpacerPlacer/evOr. Keep the existing conservative
   evolutionary grouping separate from broader locus-comparison groups. Add
   independent organism-tree or genomic-relatedness context when available.
6. Join Pipeline 1's exact-sequence target evidence to observed differences and
   reconstructed events. Preserve the mapping from each model spacer cluster
   to its sequence variants; a representative sequence's matches must not be
   assigned to every member or to an unknown ancestral sequence.

**Completed product:** aligned loci, a table of observed differences, applicable
reconstructed histories and an event-to-target evidence table. A researcher can
ask whether a set of newly observed spacers has candidate matches to the same
viral group, or whether related isolates have different target repertoires.
These are testable sequence-based hypotheses, not evidence of changed resistance.

CCTK connects repertoire networks, array comparison and historical analysis.
CRISPRscope's associated study combines genomic relatedness with spacer sharing,
diversity and target analyses. These provide precedents for the operations;
their study-specific rate assumptions are not portal defaults.
[CCTK](https://pmc.ncbi.nlm.nih.gov/articles/PMC10457644/),
[CRISPRscope study](https://link.springer.com/article/10.1186/s40168-022-01328-6).

**Optional longitudinal branch:** attach real sampling dates and lineage/subject
metadata to observed differences. First observed presence is not an acquisition
date; strain replacement, missing sequence and uncertain array direction remain
alternative explanations. Date labels alone do not justify a turnover rate.

**Optional group comparisons:** metadata-aware contrasts can follow descriptive
matrices, but require a defined sampling unit, contrast and treatment of lineage,
subject, site and multiple testing. Pairwise distances are not independent
biological replicates. Leave results descriptive when the design cannot support
the requested inference.

## Pipeline 3 Compare conserved RNA and regulatory regions

**Question:** Which repeat, tracrRNA or adjacent-region hypotheses have
comparative support, and which variants warrant closer examination?

1. Reuse candidate locus relationships to select related examples from uploaded
   genomes or frozen references. Selection rules differ for proteins, repeats,
   RNAs and whole loci; a repeat-similarity hit alone is not locus orthology.
2. Extract source-bound Cas proteins, repeat instances, both array flanks and
   supported tracrRNA regions. Retain unresolved orientation hypotheses.
3. Align credible sequence families. Report coverage, gaps, duplicate sampling,
   conserved blocks, substitutions and indels with links back to original
   coordinates. Show contributing sequences rather than only a sequence logo.
4. Overlay existing repeat-context folding and native tracrRNA evidence. The
   installed tracrRNA tool already retains anti-repeat alignment, IntaRNA duplex,
   covariance-model and terminator evidence; another duplex engine is not the
   first missing piece. Comparative consensus RNA structure is an additional
   method to evaluate for suitable alignments.
   Keep transcript orientation, acquisition polarity and source-coordinate strand
   distinct, and compare compatible transcript-orientation hypotheses. Do not
   pool folding metrics into a processing-confidence score or infer cleavage
   positions from them.
5. Export a candidate dossier: conserved and variable regions, supporting and
   conflicting evidence, uncertain boundaries and representative sequences.

**Completed product:** a biologist can select related candidates for transcript
boundary mapping, RNA-processing studies or regulatory follow-up from one
inspectable evidence collection. This is deeper than listing a fold or an RNA
candidate independently.

Comparative tracrRNA research combines repeat pairing, RNA clustering, multiple
alignment and consensus structure prediction; this directly motivates the
proposed connections. Its study systems do not establish the installed portal's
model coverage. [Comparative tracrRNA study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6546382/).
The current tool's boundary semantics remain applicable.
[Installed tracrRNA contract](../tools/CRISPRtracrRNA-v3/README.md).

**Leader branch:** adjacent-sequence conservation can prioritize regions for
study, but must not be presented as a validated leader-boundary predictor.
The original leader work motivates a dedicated reference/model and boundary
evaluation rather than assuming an ordinary alignment solves the problem.
[CRISPRleader paper](https://academic.oup.com/bioinformatics/article/32/17/i576/2450784).

## How the pipelines improve one another

The strongest reusable relationship is
**sample → locus → array → spacer occurrence → exact sequence → target match**.
Candidate Cas/RNA associations attach to the locus/array; comparative groups
attach to locus occurrences; reconstructed events attach through an explicit
spacer crosswalk. These connections let one result become another analysis's
input without confusing observed and inferred objects.

For example, related loci can supply more comparative RNA evidence and a
better-defined effector-associated spacer set. The target workflow can then
collect flanks from supported matches. A conditional motif method can analyse
those contexts, and its result can be compared across related effector groups.
A published Cas9 pipeline already connects Cas clustering, associated repeats
and spacers, orientation and dereplicated phage/plasmid targets to PAM inference.
This is precedent for synergy, not a claim that its rules generalize to every
CRISPR type. [Connected PAM study](https://www.nature.com/articles/s41467-026-69098-5).

Preserve dependencies instead of creating circular support. A reference chosen
for its motif cannot independently validate that motif; PAM-derived orientation
cannot independently confirm the same PAM. Matches for motif discovery must be
selected without requiring the motif under investigation. Show independent
target contexts and related-reference redundancy. Use held-out families or
target groups for evaluation; label unvalidated combinations exploratory.

## Council decisions and implementation sequence

The locus-biology reviewer prioritized candidate array–cassette associations and
small comparative sequence sets because they improve downstream interpretation.
The target and comparative reviewers favored completing the target workflow
first, followed by linking cohort differences to targets. The combined decision
is to make explicit, possibly ambiguous biological associations part of the
shared input contract, while allowing sequence searches to proceed without
resolved effector assignments.

1. Deliver **Pipeline 1** with source-preserving sequence collection, relationships,
   target reconciliation, interval joins and support matrices. Permit imported
   annotation tracks so a new MGE predictor is not a prerequisite.
2. Deliver **Pipeline 2**, reusing those same objects and adding locus comparison,
   observed-array analysis and sequence-to-event mapping.
3. Develop **Pipeline 3** from the shared sequence/alignment methods. Expose
   existing RNA evidence early; release new comparative inference after its
   validation. Leader research can begin alongside the first two pipelines.

The initial sequence collection, interval joins and matrices should be treated
as first-class methods with downloadable results and explicit policies. External
packages earn a place where they implement a needed algorithm: alignment,
comparative RNA structure, gene synteny or motif analysis. A package name does
not itself define a complete scientific pipeline.

## Contract and validation requirements

Keep exact sequences, biological occurrences and analytical clusters separate.
Keep observed presence/absence/unknown separate from model gains/losses. Local
locus callability, such as a represented neighborhood, is necessary before a
missing detection can support absence; overall genome completeness is insufficient.
Every matrix states whether it counts sequences, copies, loci, samples, sites or
inferred events, and identifies skipped, incomplete and unassessed cells.

Downstream calculations must consume full independently validated artifacts,
never capped UI previews. Current summaries deliberately truncate sequences and
matches: [sequence previews](../backend/src/crispr_web/summaries/detection_sequences.py),
[viral previews](../backend/src/crispr_web/summaries/viral.py).
Joining separate jobs requires authorized access to each input artifact and
must respect retention; a workflow does not grant indefinite persistence.

Use explicit sequence and interval transformations, stable identifiers, tool and
reference versions, selection policies and dependency records. Reusing a valid
intermediate requires matching input identity and effective scientific settings.
Run scientific operations only in isolated workers, preserve independent backend
validation and provision references through operator-controlled releases.

Validate scientific joins with duplicated sequences, divergent variants,
negative-strand intervals, circular origins, incomplete loci, mixed systems,
multimapping targets and censored searches. Separately evaluate homology,
alignment, motif and historical inferences using suitable biological controls.
No deployment or application implementation is included in this planning task.
Only repository planning notes changed; local links and the private-data scanner
were checked. Application tests and scientific benchmarks were not run.
