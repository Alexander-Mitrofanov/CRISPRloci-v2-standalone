# Council review: benchmark validity, reproducibility and acceptance

Phase 1 review, 2026-09-27. This is a predeclared evaluation plan for the next
implementation, not a report that it already passes. The same validation council
member should review the implementation and its actual evidence against this
file. This review changes no detector code, historical benchmark outputs, upstream
repository or `old_version` file.

## Evidence inspected and baseline

Inspected `docs/LITERATURE_BENCHMARK.md`, `docs/VALIDATION.md`,
`docs/IMPROVEMENTS.md`, `benchmarks/results/methodology_audit.md`, the canonical
and diagnostic scripts, benchmark-method tests, natural-locus/native-tool tests,
and detection/ranking/export orchestration. The following additional quantities
were computed read-only from retained reports; no new scientific run was made.

The reference comparator is **`benchmarks/final_fixed/complete`**, not the much
older failing `baseline/complete` stage. Its frozen figures are:

| Quantity | Comparator |
| --- | ---: |
| Input loci completed | 79/79 |
| Uniquely mapped references | 78; Ghc2 remains unassessed |
| Localized, correct contig and strand | 68/78 |
| Top 1 / top 3 / top 5 / top 10 | 53 / 61 / 63 / 64 of 78 |
| Mean reciprocal rank; misses contribute zero | 0.7283666082366188 |
| At least 90% reference coverage | 64/78 |
| Raw reported hypotheses | 3,395 |
| Origins | 3,390 anti-repeat; 5 covariance-model |
| Candidates per input: median / p95 / maximum | 12 / 192 / 464 |
| Total input bases; candidates per Mb | 1,364,763; 2,487.61140212623 |
| Candidate lengths: median / p95 | 191 / 202 nt |
| Sum of per-locus elapsed time | 843.9207 seconds |
| Per-locus elapsed time: median / maximum | 6.3041 / 74.1794 seconds |

Here p95 is nearest rank, `sorted(values)[ceil(0.95*n)-1]`. Historical elapsed
values include competing workers and are not batch wall time. They are useful
for estimating run cost, not for claiming speedups. Baseline software evidence
is 78 passing tests, including seven real-tool integration cases.

## Problems the evaluation must prevent

1. **Development data presented as independent validation.** All 79 loci, the ten
   misses and their isolated-RNA alignments have already been inspected. A random
   split of these loci now is not a held-out evaluation. Training-source and
   homologous-sequence overlap remain unresolved. The existing collections are
   development/regression panels, regardless of publication date. Freeze one
   proposed configuration before the full comparison and list every exploratory
   configuration. If results motivate another change, label that as a second
   development iteration with another frozen protocol; never reuse the word
   held-out for it.
2. **More candidates mistaken for better accuracy.** These are positive loci,
   not exhaustively annotated genomes. Unmatched predictions are unverified,
   not established false positives. No precision, specificity, FDR, calibrated
   probability or genome-wide accuracy claim is permitted. Report candidate
   burden, raw hypotheses and distinct intervals separately from localization.
3. **Long windows gaming coverage.** The established rule accepts 20 nt and
   half of the shorter interval, so both short conserved regions and very long
   enclosing windows can pass. Preserve that rule for comparison, but report
   intersection-over-union (IoU), prediction length and both endpoint errors.
   A 90%-coverage gain alone cannot justify a transcript-boundary claim.
4. **Truth-assisted diagnostics leaking into detection.** Supplying the known
   RNA, published repeat, truth coordinates or failure labels is allowed only
   in explicitly diagnostic analyses. Final native runs receive source FASTA
   and ordinary configuration only. The detector must not read literature
   truth tables, catalogue hashes, reference labels or a list of known misses.
5. **Historical evidence overwritten.** `summarize_literature_benchmark.py`
   currently always writes `benchmarks/results/genome_*` and `metric.json`,
   including when `--stage` is supplied. Do not run it for this iteration until
   a destination option exists. `diagnose_literature_failures.py` also writes
   fixed diagnostic/results paths; its existing outputs are read-only here.
6. **A seed that is not stable per specimen.** Existing mononucleotide shuffle
   generation advances one RNG through a sorted, changing pool. Retain those
   historical FASTAs and counts. New controls must derive seeds independently
   per source sequence, model pool, control kind and replicate.
7. **Success or rank asserted from the wrong artifact.** The current scorer
   regards the existence of `report.json` as success and trusts its source
   linkage. New scoring must verify completion, input hashes, expected source
   IDs and manifest membership. Do not substitute cached reports generated
   under another configuration. Preserve or explicitly specify tie ordering;
   a new tie-breaker can change top-1 without adding evidence.

## Preservation and provenance contract

Use the coordinator's append-only iteration directory
`benchmarks/council_round2/` and its pre-change `baseline/` source archive,
79 copied reports, metrics and manifest. Keep new manifests outside historical
`benchmarks/results`. Verify that the preservation coverage includes the
following without rewriting them:

- All files in the 79 `final_fixed/complete` run trees, including raw native
  evidence, logs, exports, input FASTAs and completion/failure artifacts.
- `benchmarks/results/`, `benchmarks/literature/loci.json`, source FASTAs named
  by that manifest, `local.toml`, `upstream-lock.json`, the current source,
  tests and benchmark scripts, and this plan.
- Model/index and native executable hashes, upstream revisions and repository
  cleanliness. There is no Git repository at the `new_version` root in this
  environment, so source checksums/a source archive are required, not an
  invented commit ID.

Retain a source archive for the pre-change detector or use the already frozen
source archive after verifying its recorded hashes. A copied manifest is cheap;
copying every large historical file is unnecessary if its digest and location
are retained and all end-of-run checks confirm no mutation. Do not rely on a
hard link as an independent backup of a mutable file.

Each new run gets a new directory and a manifest with protocol/version, exact
ordered input list and bytes hashes, settings and profile, source-file hashes,
model/native executable hashes, Python/platform/upstream revisions, commands,
worker/thread counts, timestamps, batch wall time, per-input result status and
control-generation parameters. Retain tool outputs and failure artifacts.
Refuse an existing output directory. Hash inputs and protected outputs again
after evaluation. Timing/paths are allowed to differ across repeat runs;
normalized biological outputs must be compared independently of those fields.

## Frozen metrics and coordinate rules

Keep the 78-locus denominator; an execution error contributes zero localization
and reciprocal rank. Keep Ghc2 in the 79-input completion and candidate-burden
panels with null biological outcomes. Keep Ssa and Sth1A as separate natural
loci despite their duplicated reference sequence. Source-row counts, distinct
sequences and natural loci are separate units.

Use the existing canonical localization predicate: identical source contig and
strand, inclusive overlap at least 20 nt and at least half of the shorter
interval. Reference coverage is the best fraction among matching candidates.
Compute top-k and reciprocal rank from the complete raw ranking; also show a
separate grouped-locus ranking if grouping is implemented. Do not quietly
replace the raw ranking denominator with grouped candidates.

Export one paired row per locus, including completion, localized status,
best rank, maximum reference coverage, candidate count, candidates/Mb, origin
counts, elapsed time and changes from baseline. Add:

- Top 1/3/5/10, fixed-denominator MRR, lists of gained/lost detections and ranks;
  rank ties and their deterministic handling.
- Median, nearest-rank p95 and maximum raw candidate count and lengths; total
  candidates/Mb; strongest-evidence subset counts with the rule stated first.
- At the highest-ranked localized candidate: coverage, IoU and signed genomic
  start/end error. Also report best attainable IoU separately, without
  selecting a different candidate for every displayed metric. Convert start/
  end to transcript 5-prime/3-prime interpretation using strand where needed.
- Candidate origins and boundary kinds, accepted-array count, duplicate
  same-strand intervals, new weak-evidence hypotheses, and rescue suppression.
- Counts and ranking metrics per evidence tier if an experimental short-match
  branch introduces tiers. Define tier eligibility before runs, retain every
  raw hypothesis, and show both the combined output and tier-specific views.
  Do not improve a denominator by dropping the weak tier from the count table.
- Run errors/timeouts independently from biological misses and all-control
  counts as descriptive observations, never estimated biological FDR.

Every prediction in all 79 new reports must reconstruct exactly from its
reported source interval and strand. Public intervals are 1-based inclusive;
internal intervals are 0-based half-open. Require `1 <= start <= end <= L`,
`end-start+1 == len(sequence)`, bounded local evidence, correct local/genomic
round trips, and matching CSV/FASTA/GFF/report IDs, strand and sequence. Reverse
complementing a source maps `[s,e]` to `[L-e+1,L-s+1]` and flips strand. A reverse
CM hit within an oriented window must not become forward support by accident.

For candidate consolidation, retain raw rows and their provenance. Tests must
cover same interval/opposite strand, adjacent intervals, chained overlaps,
multiple arrays and competing origins. Preserve each supporting hypothesis;
raw duplicate suppression and a biologically justified locus definition are
different operations.

## Deterministic controls and adversarial checks

Freeze generated control FASTA bytes and a manifest **before inspecting their
results**. Use namespace `crisprtracrrna-council-controls-v1`, master seed
`20260927`, and a per-item integer seed derived from SHA-256 of the UTF-8 string
`namespace|master_seed|model_pool|source_sha256|control_kind|replicate`.
Use explicitly versioned generation code and stable identifiers. Adding,
removing or reordering other inputs must leave an item's controls unchanged.
Do not use Python's process-randomized `hash()`.

A small optional sequence control panel is the 78 distinct Gasiunas references:
one mononucleotide and one exact-dinucleotide shuffle per sequence, screened
alongside their originals with the unchanged Type II model. A dinucleotide
shuffle must preserve length and all ordered adjacent-pair counts; a trivial
or unchanged realization is recorded rather than forcibly altered. Preserve
IUPAC symbols/ambiguous positions explicitly and record whether they partition
canonical runs. Record collisions with the source or another control. Require
unit tests for these invariants and catalogue-order independence. If exact
dinucleotide generation cannot be validated, ship only the correctly named
mononucleotide control and mark the stronger control incomplete; do not call a
window shuffle dinucleotide-preserving. This panel is optional for the bounded
round-two change if model resources and model thresholds remain unchanged;
the genomic controls below are the relevant mandatory stress experiment for
short anti-repeat proposals. Do not rescan the full 42,320-sequence inventory
merely to retest unchanged model recognition.

This small sequence panel tests model agreement and composition dependence,
not the whole pipeline. Native genomic controls are also needed for any new
short-match or weak-model rescue branch. Freeze these 19 development loci:
`Esp2, Bbo, Cca2, Cca1, Cga, Cme1, Ffr, Ghy2, Orh, Psp, Ain, Kki, Lsp1, Lan,
Rsp, Cme3, Ghe, Mse, Sth1A`. They span all known misses, prior native failures,
partial CM regions and a high-burden case. This is a deliberately selected
stress panel, not a representative random sample.

For each, generate one length/composition control retaining accepted array
intervals and shuffling non-array sequence segments with the same versioned
scheme; preserve ambiguous positions. Run the unchanged baseline and proposed
profile on exactly those same control bytes. Record retained array sequence,
newly inferred arrays/Cas context, candidate counts per Mb, strongest evidence,
execution and runtime. Destruction of genuine Cas context is a limitation, so
zero output does not establish specificity. If costs prevent full controls,
run the new native search/CM branch using frozen accepted arrays as a clearly
labeled **conditional component test**; that cannot substitute for a claimed
end-to-end specificity result. Do not tune cutoffs on these controls after
calling them an evaluation panel.

Adversarial software/native checks also include: wrong strand with strong
score; correct strand with a 19-nt hit; exact threshold equality; no arrays;
no CM hits; array-overlapping weak hits; repeated contig IDs across files;
long paths containing spaces; clipped terminal windows; malformed native
output; and candidate order under permuted input records. Repeat a small
native subset twice and compare normalized results to expose nondeterminism.
A complete-pipeline reverse-complement stress run should be reported even if
upstream array inference is asymmetric; parser/coordinate transformations
remain strict invariants, while new upstream biological differences must be
explained rather than hidden.

## Affordable run order and acceptance gates

The summed retained time (about 14 minutes for 79 loci) makes a complete real
comparison practical. Historical runs used two loci in parallel, one tool
thread per locus. Preserve that budget, log actual batch wall time, and avoid
simultaneous unrelated workloads when comparing speed. This estimate is not a
promise that a sensitivity profile or control input has the same cost.

Run in this order after one configuration is frozen:

1. Existing 78 tests plus focused tests for the scientific changes; run all
   seven existing native cases with the real configuration, not skipped.
   Run lint/format/type checks required by the changed files.
2. Native checks on `Ain, Kki, Lsp1, Mse, Cca1`, including a repeated run and
   coordinate/export audits. This catches the known path/parser regressions
   and expensive candidate growth before launching the full panel.
3. A fresh full 79-input default run, even when most changes concern an opt-in
   profile. Rescore it into the new iteration directory and audit every
   prediction/export. Historical files remain byte-identical.
4. One fresh full 79-input run of any proposed sensitivity profile, using the
   same manifests and resources, plus the 19 paired genomic controls above
   (38 native control runs: unchanged baseline and proposed profile). The
   primary budget is 79 default + 79 experimental + 38 control runs. Report branch
   ablation if a new rescue rule claims credit: disable that rule on the same
   inputs and configuration or retain enough branch provenance to establish
   the difference. Retained-raw-hit rescoring is diagnostic, not a native run.
5. Build/install smoke testing only when packaging or distributed resources
   change; otherwise the existing native/source checks are the relevant gate.

**Hard correctness gates for any shipped change:** no unexpected test skips or
failures; no coordinate/export/provenance invariant failure; no protected-file
mutation; all 79 natural inputs complete; correct Ghc2 handling; deterministic
control generation; no truth/catalogue-dependent code path. Any violation is a
failure requiring correction, independent of apparent recall improvement.

**Default compatibility gate:** retain all 68 baseline-detected loci, at least
53 top-1, 63 top-5 and 64 at 90% coverage. If the implementation claims unchanged
default biological behavior, normalized default predictions and ranking should
match the comparator, not merely its aggregates. An intentional default change
requires the council to review each gained/lost rank and all burden differences;
changed defaults cannot be justified by a positive-only recall gain alone.

For exact-repeat IntaRNA grouping, require equality of every candidate's parsed
interaction result (including null, coordinates, energy and structure) and of
downstream evidence/ranking in the 79 fresh default reports. Group on the full
oriented query sequence and any other interaction parameters, not just an
array ID, repeat hash alone or un-oriented genomic sequence. Check identical
sequence reuse across different arrays, repeated target sequence with distinct
candidate IDs, multiple output rows, an empty interaction result, out-of-order
rows, duplicate/unknown/missing IDs according to the native contract, malformed
rows, and any batch/single fallback. Record actual command counts and time;
the reported planning result of 3,390 calls versus 152 exact-repeat groups is
a proposed optimization opportunity, not measured end-to-end performance.
Default behavior equivalence is a hard gate for this engineering change, even
if an experimental profile improves detection. Compare parsed numeric values
and normalized biological dictionaries; textual JSON whitespace/key order and
run-specific paths/timing are excluded. No numeric tolerance is needed for
the claimed exact batching optimization. If a mismatch appears, retain and
explain it; do not introduce a tolerance after seeing the results.

**Opt-in sensitivity-profile gate:** retain all 68 baseline detections and the
53/63/64 aggregate floors, and demonstrate at least one additional localized
reference (69/78 or more) in actual native output. Pure ranking or consolidation
work instead must meet its predeclared aim (at least one top-1 improvement or
at least 20% reduction in displayed grouped entries while preserving raw
hypotheses); do not claim detection gain for either. All new results remain
regression/exploration claims. Leave wider-coverage/weak-CM settings explicitly
opt-in unless independent negative/held-out evidence later supports promotion.

**Operational burden budget for a sensitivity profile:** total raw hypotheses
at most 1.5 times baseline (at most 5,092 across 79 inputs), p95 at most 288,
maximum at most 928, and measured comparable batch wall time at most twice a
fresh baseline run. These are practical guardrails, not biological calibration.
If historical timing alone is available, report performance as descriptive and
do not claim the timing gate passed. If the profile exceeds a guardrail, keep
it experimental or narrow the proposal; do not quietly raise the bound after
observing results. Raw count budgets still apply when a grouped view is added.

**Control review gate:** compare paired baseline/new burden and strongest-
evidence counts for every fixed control. A systematic rise (more than half of
controls gain candidates), a >1.5-fold aggregate burden increase, or any new
highest-support control prediction requires explicit council review of native
evidence and blocks promotion to default. A zero baseline denominator must be
shown as zero versus the observed absolute count, not a finite ratio. Passing
controls permits a bounded experimental feature, not a precision claim.

Council acceptance may record “implementation correct; scientific improvement
not established” or “useful experimental profile; default unchanged.” A failed
scientific target must remain visible even if the code is maintainable. The
final evaluation should link this frozen plan, implementation source hashes,
actual commands/logs and paired machine-readable outcomes; list each gate as
pass, fail, or not measured with its evidence.
