# Repeat analysis — council implementation decision

27 September 2026. Design decision for roadmap **P4**, parity items **F04–F06**.

Implementation follow-up: the user requested an independent tool. P4a and repeat-only folding now live in [CRISPRrepeat](https://github.com/Alexander-Mitrofanov/CRISPRrepeat), integrated through the queued worker and web results. See [implementation scope and validation](../../docs/repeat-analysis.md). Predictive classifiers and a curated public reference database remain deferred under the gates below.
The original decision below records an internal council of three agent reviewers and a coordinator; this is not independent experimental validation. The first implementation has engineering and selected published-substrate tests, without generalized predictive accuracy claims.

## Decision and scientific purpose

Build a **repeat evidence workflow** with two inputs: a detected array with its
source sequence, and user-supplied repeats. Start with observed repeat instances
and their array context; add local reference retrieval; gate inferred strand,
subtype and family assignments on their own evaluations. Use a focused
Integration adapter with pinned ViennaRNA tools and existing detector assets.

The first scientific question is: **which repeat instances change support for
a proposed RNA structural motif when their surrounding array sequence is
included?** The second is: **which documented repeat families and subtype
associations resemble this sequence, and where does the evidence disagree?**

Context-dependent processing has experimental precedent in particular Cas6
systems [B1]. This motivates prioritizing instances for follow-up, not a
universal processing predictor. Weak folding must never remove a detected
array. Repeat evidence must not claim immune activity, establish host identity,
or override CasAndra, tracrRNA evidence or evOr.

Source IDs below resolve in [EVIDENCE.md](EVIDENCE.md). All algorithm choices,
metrics and delivery gates below are council recommendations, not results of
a benchmark performed during this task.

## 1. First release: each repeat in its observed context

1. Read canonical detector records and reconstruct the array and every repeat
   directly from immutable input sequence. Analyze actual instances, including
   substitutions, indels and terminal variants. Fold the consensus separately;
   it is not a substitute for an observed repeat and may not occur in the array.
2. Convert each selected transcription-orientation hypothesis to an explicit
   RNA analysis sequence. Compute isolated MFE/dot-bracket and partition-function
   pair probabilities using RNAfold. Deduplicate identical folds within the job
   while preserving every occurrence and its coordinates.
3. Run RNAplfold on the observed, oriented array span. Retain per-instance
   intrarepeat pairs, individual pairs to spacers/other repeats, and native
   single-nucleotide unpaired probabilities. A genome span is an analysis
   substrate, not an experimentally established full precursor transcript.
4. Start evaluation with the historical **150-nt window / 80-nt maximum span**
   configuration [B2]. Freeze the effective temperature, energy parameter file,
   dangling-end model, lonely-pair policy and executable versions. The initial
   reference temperature is 37°C, explicitly a computational setting, not an
   estimate of the organism's growth temperature. Compare lonely-pair policies
   and other settings during development before freezing a public preset.
5. Show repeat instances in array order with an isolated/context comparison,
   a motif-pair overlay and a matrix of competing interactions. Keep an
   inspectable table and full downloads behind any graphical summary.

RNAplfold scans the whole submitted array but uses local windows. Label this
**local folding in array context**; it does not include every possible
long-range interaction. Do not manufacture a single full-array MFE dot-bracket
structure from local pair averages. An isolated fold can have a structure
drawing; the context view should show pair support and interactions.

Reimann's experimental reproduction uses the exact tested fragments, including
synthetic additions, and its own reported settings; it is a separate benchmark
preset. Longer-fragment effects make an arbitrary small focal window an
inadequate replacement for this analysis [B1]. Whole-fragment RNAfold or changed
local windows can later provide explicitly labeled sensitivity analyses.

## 2. Define the comparison before implementation

For a repeat instance, let **M** be an explicitly stored, nonempty set of
repeat-local base pairs. Its origin is either a mapped curated motif, a family
consensus with stated evidence, or the instance's isolated MFE hypothesis.
An MFE hypothesis must not be called a known functional motif.

Let **P(i,j)** be the isolated RNAfold pair probability and **Q(g(i),g(j))**
the RNAplfold local pair support after mapping into the oriented array.

```text
isolated_pair_support = mean(P(i,j) for (i,j) in M)
context_pair_support  = mean(Q(g(i),g(j)) for (i,j) in M)
support_change        = context_pair_support - isolated_pair_support
```

The two supports lie in [0,1]; the change lies in [-1,1]. They describe average
support for the listed pairs. They are not the probability the entire motif
forms, experimental structure accuracy, or a cleavage probability. The methods
use different ensembles, so the difference is a descriptive model comparison,
not a causal biological effect size. Show both supports and the pair count,
not only their difference; two weakly supported motifs are not positive evidence.

For an unpaired MFE, **M is empty and the comparison is null**, with reason
`no_reference_pairs`. Folding still completed successfully. Keep its ensemble
and accessibility output. For a mapped motif with missing positions, report
pair coverage; withhold a directly comparable full-motif score. Any partial
score must carry a separate label and its actual denominator. Never silently
drop gapped endpoints or score an unmappable pair as zero.

Freeze M before examining the context result. Pair sets derived independently
from different instances need not be homologous; keep them distinct in the
instance table. Cross-instance motif comparisons require a common mapped motif
and explicit coverage. Check that each motif pair fits the effective context
window/span. A pair forbidden by those settings makes the full comparison
unavailable; it is not evidence of biological loss of that pairing. Preserve
valid isolated outputs when only context support and its delta are unavailable.
Do not pool distinct pair sets as support for one shared motif. Also check
model eligibility: noncanonical/protein-mediated interactions and pseudoknot
motifs may be incompletely represented. Label that limitation and withhold
unsupported motif comparisons rather than infer dysfunction from a low score.

**Numerical rule:** RNAplfold averages pairs over different sets of windows;
its row sums can exceed one [N3]. Do not calculate `1 - sum(Q)`, partner entropy,
total-variation distances or a summed "external pairing probability" from this
matrix. Use native `-u 1` accessibility and individual competitor pairs. Native
segment accessibility is required if a multi-base footprint is later added;
products of singleton probabilities are not equivalent.

Compute metrics before display sparsification. Use numeric library output or
a strictly parsed private tool output with verified probability semantics;
RNAfold dot-plot entries contain square roots of probabilities [N1]. Record
numeric/reporting cutoffs. If exact motif-pair values are unavailable below a
cutoff, supply an error bound or mark the metric unavailable; absence from a
truncated plot is not a measured zero. Test against complete small matrices.

## 3. Input, orientation and interpretation contract

| Contract area | Required behavior |
| --- | --- |
| Molecule and normalization | Preserve original sequence/hash and DNA/RNA declaration. Produce uppercase RNA and, where needed, DNA model views with a transformation log. Reject mixed T/U in the initial repeat-only contract. Ambiguous bases remain visible; the first folding implementation reports unsupported affected substrates rather than inventing bases. |
| Coordinates | Normalize internally to zero-based half-open source-forward intervals. Explicitly convert the detector's one-based inclusive intervals. Store repeat-local to array-local to source mappings, with reverse-complement transformations. Repeat-only input has no invented genomic coordinates. |
| Source boundaries | Use all source bases in the observed array span, not concatenated consensus copies. Record clipping and uncertain transcript ends; flag terminal instances for boundary sensitivity. Circular genomic coordinates do not imply circular RNA. Support wrapping only with an explicit segmented mapping and tests; otherwise return unsupported. |
| Transcription direction | Keep supplied orientation, detector evidence, future leader evidence and evOr acquisition-order evidence separate. For unresolved DNA orientation, retain both hypotheses within bounds. RNA declared as the transcribed 5′→3′ sequence is folded as supplied. Do not select strand by the more stable fold. |
| Missing context | A repeat-only record supports isolated folding and reference analysis; array-context evidence is unsupported with reason `context_not_supplied`. |
| Candidate categories | Retain source quality/category and relationships between alternative candidates. Repeat analysis does not require enough related arrays for evolutionary reconstruction. Analyze uncertain categories only under an explicit request/policy, retaining the distinction from the default cohort. |
| Stage state | Give folding, context, strand, retrieval, family and subtype separate `not_requested`, `completed`, `no_hit`, `unsupported`, `failed`, `truncated` states and reasons. Null metrics have their own reasons. These do not replace durable job states. |

Protein-stabilized non-hairpin Cas6 repeats, tracrRNA-dependent processing and
Cas12a pseudoknots establish limits of a generic RNA-only model [B3–B5]. Show
applicable Cas/tracrRNA evidence alongside folding, preserving uncertain locus
associations. Do not generalize a processing rule across all Cas12 subtypes.
Equilibrium folding also does not model transcription kinetics or protein binding.

## 4. Repeat-only references and classification

**Select local reference retrieval as the next capability.** Build an immutable,
operator-managed reference release with original sequences, accessions and
versions, array locations, sequence orientation, organism/taxonomic snapshot,
annotation provenance, licenses, checksums and build recipe. Deduplicate search
sequences while retaining all occurrence records. Unresolved historical origins
remain unresolved; a consensus reference need not be an exact genomic substring.

Search both orientations where meaningful, reporting exact matches separately
from near matches. Retain alignment, identity, query/reference coverage,
mismatches, gaps, scoring parameters and all equally supported matches, subject
to explicit reporting caps. Use exhaustive short-sequence alignment as the
small-reference oracle and an initial implementation candidate; choose any
index after measuring recall and cost. Source organisms are reference origins,
not a host identification for the submitted repeat.

**Reuse the portable CRISPRstrand model as a candidate strand method.** The
current asset is a CNN with two independent orientation scores and 38-base
truncation, not the published 2014 SVM merely under a new packaging [C3, R2].
Preserve its raw outputs and preprocessing; never normalize them into a
posterior by assumption. Add an interpretation layer that can report unresolved
ties/palindromes and calibrated uncertainty, without changing legacy detector
behavior. Until evaluated, raw method evidence cannot establish a confident
orientation. Training provenance is an unresolved scientific gate.

**Select standalone RepeatTyper as the first subtype candidate to benchmark**
against a nearest-reference baseline and CRISPRclassify [R3–R4]. Do not add a
second complete locus pipeline. Pin the model independently from its executable
and evaluate only labels actually represented in that model. Do not import an
upstream score threshold as calibrated local confidence. Preserve tied or
conflicting subtype associations, especially for repeats shared by different
systems; nearby Cas context is separate evidence, not an automatic tiebreaker.

**Keep CRISPRmap family/motif assignment separate** [R1]. Its historical
download and later family vocabulary differ. A nearest sequence match is not a
validated family assignment. Restore a user-controlled CRISPRmap handoff with
clear destination and selected sequence; no automatic external submission.
Local assignment remains gated on recoverable references/model provenance,
redistribution rights and comparison against the correct historical release.

## 5. Repository implementation boundaries

| Layer | Implementation work |
| --- | --- |
| Scientific adapter | Add a focused repeat module/contract under `Integration/crispr_integration/`; keep stage ordering in the existing runner. Read canonical detector geometry, not the evolutionary `ArrayRecord` projection. Expose a repeat-only CLI path and optional detected-array stage. |
| Runtime | Verify RNAfold and RNAplfold in the existing pinned ViennaRNA environment, then add operator-owned shell-free commands and a real scientific smoke test. Record binaries, parameters and strand/model identities in the scientific release. |
| API/worker | Add typed repeat input and optional analysis settings through the API, queue serialization and worker invocation. The API only validates/queues; contained workers execute. Preserve current modes, checkpoints, deadlines and cancellation. |
| Independent validation | Add Integration validation and an independent backend `execution/validation_repeats.py`. Both reconstruct sequences and mappings from immutable inputs. Validate IDs/counts, hashes, finite energies, pair bounds, motif arithmetic, dot-bracket syntax, reference bindings and truncation. Enforce row-sum bounds only for global marginals, never RNAplfold averages. |
| Publication | Publish approved repeat JSON, instance TSV, analysis FASTA, isolated dot-bracket and pair tables. Keep private tool files out of the public bundle. Sparse previews declare cutoff/caps; exact metric inputs remain reproducible from approved artifacts. |
| Frontend | Extend `src/contracts/public-api-v1.json` and focused results components. Provide instance selection, pair-support changes, matrix inspection, orientation alternatives and separate reference/subtype evidence. Render supplied scientific results rather than inferring biology in the browser. |

The existing detector validates native construction, but today's backend
identify gate primarily binds report-level metadata; it does not independently
reconstruct each repeat [C1–C2]. Adding detailed source checks is prerequisite,
not optional hardening after the feature is shown.

Bound record count/length, array length, repeat instances, orientation branches,
pair-table size, output bytes, CPU, memory and elapsed time before allocation.
Measure numeric limits in a development pilot. Cache only within the job using
sequence, orientation, full settings and tool/parameter identity. On a cap,
report affected instances and stage state; never silently shorten context or
return a score from incomplete probabilities. A failed optional analysis must
not turn missing evidence into a negative result or erase valid detection.

Integration changes require updated contracts, fixtures, version/wheel and
release identities. Reference builds are operator actions; running jobs remain
bound to their release. Existing containment, capabilities, retention and
frontend-only public export rules continue to apply.

The current Integration process helper uses `stdin=DEVNULL` and has no child
working-directory parameter. Add narrow support for an owned input handle and
private per-analysis working directory, or use a focused wrapper: RNAplfold
requires stdin and writes sidefiles. Keep these inside attempt storage under
the worker's existing containment. Do not use shell redirection or raw FASTA
headers as filenames. Verify actual pinned-version behavior for short arrays,
effective window sizes and output formats before accepting production results.

## 6. Demonstrate scientific value before stronger claims

| Gate | Evidence and pass condition |
| --- | --- |
| Source and numeric correctness | All hand-checked reconstruction, reverse-strand, indel, empty-motif and cutoff fixtures agree; tampered sequences/coordinates fail independently. Repeated pinned runs agree within declared numeric tolerances. No unexplained provenance or arithmetic discrepancies. |
| Mechanistic reproduction | Recover exact published Reimann fragments, motifs, assay conditions and outcome labels. Compare isolated-only and context predictions, reproduce reported context-sensitive examples and publish disagreements. This is reproduction, not an independent test. |
| Incremental biological value | Assemble independent experimentally labeled processing/structure data, with study/locus/family separation. Compare context support to isolated support and simple length/composition baselines. Report ranking/discrimination, effect sizes, coverage and uncertainty. If data are insufficient, release descriptive evidence only and make no processing-prediction claim. |
| Mechanism and sensitivity | Include non-hairpin Cas6, tracrRNA-dependent cases, Cas12a and orphan arrays. Inspect temperature/window/span/lonely-pair choices, fragment ends and reverse orientations. Spacer shuffle/swap controls test context response, not biological accuracy. Preserve parameter-sensitive and unsupported cases. |
| Retrieval | Compare to an exhaustive oracle, including reverse complements, duplicates, ties, indels, no-hit and caps. Verify every reported source/alignment and retain ambiguous origins. |
| Strand/subtype/family | Partition by canonical repeat homology and genome/strain, with separate family/clade holdouts. Keep duplicates, reverse complements and overlapping fragments together. Labels must be independent of the method being tested. Include absent families/subtypes and length/composition-matched controls. |
| Selective predictions | Report error versus accepted coverage, per-label results/confusion, calibration and unsupported-family abstention, with uncertainty grouped by locus/family. Set useful-coverage/error targets on development data with scientific owners, freeze them before the held-out run, and release only supported strata. |
| Engineering/release | Run relevant `scripts/check.sh` gates and browser acceptance for new UI/input/download/no-hit flows; add real-science release tests and bounded runtime measurements. Verify miniature reference releases offline. |

A fixed pretrained model is only "family-held-out" when its training inventory
allows overlap exclusion. Otherwise label it an external test with unknown
historical overlap, or retrain on a transparent frozen partition. Do not turn
random repeat-row splits or synthetic tests into claims of biological accuracy.

Separate known-family generalization from unknown-family rejection: a family
absent from the reference vocabulary is tested for abstention, not correct
assignment to its unavailable label. Subtype generalization across unseen
repeat families is another distinct evaluation. Freeze homology-clustering
criteria using development data and report sensitivity to those choices.

The first implementation artifacts should be accession-level manifests for
experimental fragments, reference occurrences, classification cases, orientation
cases and grouped splits. Record source hashes, label evidence, licenses and
training/reference overlap before tuning parameters. These inventories are
planned work; no independent biological benchmark set was curated in this task.

## 7. Delivery order and council resolution

| Milestone | Deliverable | Scientific claim allowed |
| --- | --- | --- |
| P4a | Source contract, per-instance isolated/context analysis, exact exports and evidence UI; reproduction and sensitivity report | Descriptive structural hypotheses and context comparisons |
| P4b | Typed repeat-only workflow, immutable reference retrieval, raw audited strand evidence and optional validated strand interpretation | Reference similarity; orientation only within evaluated scope |
| P4c | Selected subtype classifier and separately validated local family method; historical handoff available earlier | Calibrated/abstaining associations within the released vocabulary |

P0/P1 contracts and evidence surfaces are prerequisites. The reference-release
part of P6 precedes P4b; leader work P3 is optional additional evidence, not a
blocker for folding. F04–F06 remain incomplete until their actual acceptance
criteria pass; a partial milestone must not be labeled full historical parity.

The biology reviewer favored full-array local analysis over adding an arbitrary
third focal-window ensemble to the initial workflow. The reference reviewer
favored inspectable retrieval before predictive labels. The architecture
reviewer favored a focused adapter that reuses current assets. The coordinator
adopted these choices and required method-specific numerical validation:
RNAplfold row sums and summed outside-repeat probabilities were explicitly
removed from the initial architecture recommendation after review.

**Go decision:** proceed with P4a as the first implementation slice, and prepare
reference provenance for P4b alongside it. Predictive subtype/family labels and
processing claims remain gated. This decision completes the requested planning;
it does not assert that any analysis has been implemented or validated.
