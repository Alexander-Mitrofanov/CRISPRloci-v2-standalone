# Scientific and delivery plan

Status: confirmed by the project owner on 27 September 2026; implementation authorized.
The linked [evidence register](EVIDENCE.md) separates published observations,
source availability and council design choices.

## 1. Product and scope

Create CRisprleader v2 as a standalone scientific Python package and command-line
tool in its own repository. Its main question is: **given a genome and a CRISPR
array, which adjacent sequence has evidence for a conserved leader core, where
is that supported core, and what remains uncertain?**

The first release will provide:

- Validated genome/array inputs and a direct canonical CRISPRidentify adapter.
- Source-bound flanks on both array ends, including negative strand, fragmented
  records and explicitly circular topology.
- Reference-family leader-core inference and separately reported motif evidence.
- Explicit conflicts, unsupported families, partial sequence and unresolved
  boundaries, rather than a forced prediction for every array.
- Reproducible JSON, TSV, FASTA, GFF3/BED outputs; local reference verification;
  documentation, examples, CI and a scientific benchmark report.

The initial predictor covers validated reference families. Universal de novo
leader discovery, promoter/TSS prediction, experimental activity prediction,
new array detection and automatic changes to other tools are outside this
release. Motifs may provide useful auxiliary subtype evidence but must not
replace or silently revise CasAndra classifications.

## 2. Inputs and biological quantities

Primary input is genomic FASTA plus a versioned array-annotation JSON. Each
source record includes ID, length, uppercase-sequence SHA-256 and topology with
its provenance. Each array includes source-forward interval/segments, observed
terminal-repeat anchors, consensus repeat, source category and optional strand
evidence. Accept supplied annotations without requiring CRISPRidentify to run.
Optional Cas context is source-bound metadata, not a required dependency.

The Identify adapter reads its canonical report and original FASTA directly.
It must validate schema, coordinates, sequences and digests without importing
Identify or the web server. The current Integration evolutionary `ArrayRecord`
omits genomic intervals and must not become the new input model. Legacy text
labeled “Leader region” is excluded as inference input.

Keep these distinct in output:

| Quantity | Meaning |
| --- | --- |
| Extracted flank | Sequence available for analysis, with a requested window and clipping details. |
| Predicted conserved core | Interval supported by a stated reference/model and validated boundary policy. |
| Leader-side evidence | Evidence favoring one array end, separate from supplied transcription strand. |
| Motif occurrence | Sequence/position match with model, background and statistical scope. |
| Putative functional association | Literature-supported interpretation within its stated taxon/system scope. |
| Full functional leader / promoter / TSS | Unavailable unless independently supported; not inferred from the flank window or one motif. |

A later sequence-only flank mode can report local matches, but cannot invent
genomic coordinates. It is not required to replace genome-based inference.

## 3. Proposed inference pipeline

**A. Validate and normalize.** Bind annotations to the supplied genome. Preserve
repeat boundaries and category provenance; do not silently repair an array or
shift its boundary to make a leader model fit. Convert to zero-based, half-open,
source-forward coordinates internally and in native JSON. Adapt explicitly
from Identify's one-based inclusive convention.

**B. Extract candidate flanks.** Evaluate both array ends when feasible, even
when a computational strand call is supplied, so disagreement can be detected.
Preserve the supplied call; an explicit user strand constraint controls
selection but does not erase contradictory evidence. Orient candidates into
their proposed leader-to-array direction. Window policy is bounded and records
requested/available length; it can extend beyond 200 nt where justified by
reference requirements. It never defines a biological boundary by itself.

Represent circular features as ordered source-forward segments and reverse
complement the extracted arc once when needed. Unknown topology does not allow
wrapping. Never reuse the same nucleotide twice or cross the array to fill a
window. Record ambiguous bases, neighboring features/arrays and truncation;
do not assume every overlapping annotation disproves a leader.

**C. Retrieve reference families.** Shortlist plausible families using repeat
similarity and reference metadata. Cas subtype/taxonomy are soft context, not
universal exclusion rules. Record the candidate set and retrieval policy.
Benchmark shortlist recall against broader searches on small datasets so
retrieval cannot hide missed families. Preserve tied/competing families.

**D. Score core evidence.** Evaluate a DNA profile-HMM route, using a suitable
DNA search implementation such as `nhmmer`, against a simpler sequence-alignment
baseline. Inspect profile alphabet/format before importing legacy files. Raw
bitscores, E-values, coverage and alignment intervals retain their own meanings.
Family-specific calibration and explicit terminal-repeat anchoring determine
whether a supported core interval can be reported. A local match is not
automatically a full leader, and local alignment/envelope endpoints are not
automatically validated core boundaries. Partial core matches remain partial.

**E. Score motif evidence separately.** Evaluate a reproducible PWM scanner
against FIMO or use an appropriately licensed pinned implementation. Match
context includes orientation, distance to the observed leader-repeat junction,
taxon/subtype scope, background model and motif provenance. Multiple-testing
scope must include the tested motifs, candidate sides and positions; reporting
caps cannot silently change significance. A short AT-rich motif alone cannot
produce a leader call or prove integration/transcription activity.

Each motif declares its allowed domain: leader-only, leader–repeat junction,
or repeat. Junction models include the observed genomic terminal repeat,
not just the consensus repeat, with reversible genome mapping. Do not clip
junction-spanning motifs at the flank boundary.

**F. Make transparent decisions.** Start with auditable evidence rules, not an
unvalidated combined machine-learning score. Fit cutoffs on development data,
freeze the policy, and evaluate on untouched test partitions. Do not average
heterogeneous/correlated scores into a probability. Preserve independent raw
evidence and explanations for support, conflict and abstention.

**G. Validate and export.** Reconstruct every predicted genomic sequence from
its original segments and strand using a separate validation path. Export only
supported coordinate features as such; keep candidate windows/motifs explicitly
labeled. Stable IDs join all files. Missing fields remain null/absent according
to the documented schema; zero remains a legitimate measured value.

GFF3 converts to one-based inclusive coordinates; BED uses zero-based half-open
coordinates. Circular origin-spanning features export as explicitly linked
segments, never a flattened min–max interval or a BED interval with end below
start. Independent reconstruction and round-trip fixtures cover each format.

## 4. Result contract and commands

Separate execution state from scientific status and sequence completeness.
Proposed dimensions:

- Execution: completed, invalid input, failed.
- Decision: supported, ambiguous, no model-supported hit, unsupported.
- Reasons: orientation conflict, reference unavailable, motif-only evidence,
  insufficient sequence, unresolved boundary, or declared unsupported input.
- Flank/core completeness: complete, truncated, unknown/not applicable.

These states do not assert biological absence. The report retains each tested
orientation, source/array IDs, terminal-repeat anchor, sequence digests,
candidate segments, model hits, motif hits, raw scores, decision policy version,
reference hashes, source licenses/citations and effective settings. Keep
promoter and full-leader claims unavailable unless supported independently.

A required reference pack that is missing, corrupt or incompatible causes an
actionable setup/execution failure. A valid pack with no coverage for an input
family can produce a completed scientific unsupported result. Do not turn a
broken installation into a biological no-call.

Proposed commands are `crisprleader predict`, `crisprleader adapt identify`,
`crisprleader validate` and `crisprleader references verify`. Maintainer
reference build/acquisition is a separate documented workflow. Ordinary
prediction runs offline and never download or update models/tools.

External programs run through a narrow process adapter: fixed argument lists,
operator-selected paths, private work directories, deadlines, cancellation,
process-group termination, bounded logs and independently validated output.
Do not load pickle or execute code from a reference bundle.

## 5. Data and scientific acceptance gates

The reference council found two practical gates: the legacy checkout inventory
does not match the paper's reported family total, and reviewed source trees do
not declare clear repository-wide licenses for code/data. Do not guess missing
models or copy assets merely because they are downloadable. Resolve provenance
and documented reuse terms per asset, or rebuild from permitted genome sources
using independently authored methods. User-supplied packs may support research
use but do not by themselves satisfy the usable-reference requirement.

Before model selection, build a record-level evidence ledger distinguishing
experimentally supported orientation/junction/motif effects, comparative core
boundaries and unvalidated extracted windows. An adaptation experiment is not
automatically a full-leader-boundary label. Published motif-positive sets are
not automatically suitable negative-free classifier benchmarks.

Freeze family/species-group partitions before training or motif discovery;
track overlaps against all reference/build data. Use within-supported-family
holdouts to measure ordinary performance and family-held-out challenges to
measure novelty/abstention. A family excluded from a reference library cannot
be expected to receive a specific reference-family call. Add taxonomic and
time-based tests where records support them.

Measure core detection precision/recall where labels permit, coverage,
abstention, orientation accuracy against independent known-strand labels,
orientation disagreement, endpoint error and interval overlap
separately. Include shuffled/composition-matched flanks, wrong-side controls,
truncated records, reverse complements, circular-origin rotations, shifted
anchors, unrelated repeats and competing families. Treat unannotated natural
and alternative flanks as uncertain controls unless independently adjudicated;
they are not proven biological negatives. Report confidence
intervals and denominators. Exact biological cutoffs and acceptance targets
must be fixed from development evidence before final evaluation, never guessed
in this planning document or tuned on final results.

Compare fixed-window extraction, motif-only screening, sequence-family matching,
profile inference and the legacy method where reproducible/available. Test
the added value of each evidence source through ablations. No superiority claim
is made now; release claims must name their tested population and coverage.

## 6. Execution after confirmation

1. **Repository and contracts.** Recheck ownership/name; create the approved
   repository and local independent checkout at
   `tools/CRISPRleader-v2/`. Add original-code licensing, package skeleton,
   schemas/interfaces, synthetic contract fixtures and worker ownership. Do not
   edit the CRISPRloci tool lock or current runtime at this stage.
2. **Parallel worker development.** Run the three workers specified in
   [WORKERS.md](WORKERS.md), based on the same contract commit. The coordinator
   composes CLI/exports/packaging and reviews contract changes.
3. **Reference/model development.** Resolve asset provenance; build the initial
   supported reference release; freeze partitions, calibration and decision
   policy. Core engineering can proceed with synthetic fixtures, but biological
   completion cannot.
4. **Integrated validation.** Run independent reconstruction/containment tests,
   fresh offline installation, CLI flows and real-data benchmarks. Address
   findings before proposing a scientific release. Publish limitations and
   unsupported families in the repository documentation.
5. **Standalone release readiness.** Push reviewed work to the approved
   repository; preserve source/reference/benchmark identities. Do not publish
   packages or make the repository public unless that visibility/distribution
   was approved. Report the working CLI and measured results to the owner.

Server integration follows a separately reviewed release change: add the pinned
tool/reference identity, Integration adapter, independent backend validation,
safe artifacts and frontend evidence view. That change must run
`scripts/check.sh`, real worker cases and browser tests. No deployment or
production service operations are authorized by this plan.

## 7. Definition of done

The independent repository contains a working documented tool and an actual
default biological reference bundle produced and verified by the documented
setup, with fixed version/digests, supported families and recorded reuse basis.
Acquisition/build stays outside prediction and need not distribute raw genomes.
A recipe or user-supplied pack interface alone is insufficient. Known-family
core predictions on curated biological examples, honest no-call behavior,
verified exports, clean installation outside this monorepo, passing CI and a
frozen benchmark report meeting predeclared scope-specific acceptance gates
are required. If sample sizes or label quality cannot support a claim, narrow
the supported scope or keep the release explicitly prerelease.
Model/reference permissions and the legacy
inventory discrepancy are documented. The five existing tools remain intact.
If references or labels cannot support a biological claim, record that limit;
do not declare the goal complete based only on engineering tests.
