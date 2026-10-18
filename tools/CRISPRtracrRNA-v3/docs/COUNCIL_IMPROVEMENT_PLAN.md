# Council improvement plan: round 2

Frozen on 2026-09-27 before implementation or examination of new-profile results.
Coordinator: root. Implementer: `council_executor`. Independent reviewers, reused
after execution: `council_detection`, `council_engineering`, `council_validation`.

## Objective and evidence

Improve runtime and scientific-output contracts without changing default biological
predictions. Separately evaluate a disabled-by-default Type II short anti-repeat
rescue experiment. The 79 natural inputs are an inspected development/regression
collection, not an independent validation set. Success does not require recovering
every known reference or promoting experimental settings to defaults.

The council's detailed rationale is in [detection](council/01_detection_review.md),
[engineering](council/02_engineering_review.md), and
[validation](council/03_validation_review.md). This document resolves their scope
differences and is the execution contract. The machine-readable counterpart is
`benchmarks/council_round2/protocol.json`.

Baseline: 79 successful inputs, 78 uniquely mapped references, 68 localized;
top 1/3/5/10 = 53/61/63/64; 64 with at least 90% reference coverage; 3,395 raw
hypotheses (3,390 anti-repeat plus five independent genomic CM rescues). Keep
Ghc2 unassessed biologically. Keep Ssa and Sth1A as separate natural loci.

## 1. Preserve and freeze before changing code

- Use `benchmarks/council_round2/baseline/`: pre-change source archive and hashes,
  all 79 copied reports, existing metrics and baseline manifest.
- Coordinator records hashes of the full historical final run trees, results,
  literature inputs, bundled models, local configuration, upstream lock, native
  executables and protected repository states. Verify them again after execution.
- Never write new results to `benchmarks/results`, overwrite a previous run, edit
  `old_version`, change either upstream repository, or modify model resources.
- The executor freezes control generation and its checks before native control
  runs. Each independent shuffle seed derives from SHA256 of
  `namespace|master_seed|model_pool|source_sha256|control_kind|replicate`, namespace
  `crisprtracrrna-council-controls-v1`, master seed `20260927`. Adding/reordering
  inputs must not alter existing controls. Include array identity where relevant.
- Freeze this plan and protocol hashes. Any unavoidable deviation must be recorded
  before inspecting its affected result; do not tune thresholds on this panel.

## 2. Implement exact-repeat IntaRNA batching

Owner: executor. Review owner: engineering council.

- Add `intarna_batch_size`, positive non-boolean integer, default 128, TOML and
  CLI `--intarna-batch-size`. Size 1 is an explicit serial comparison path.
- Implement `interactions(candidates, work, settings, runner)` with every input
  candidate ID represented, including valid no-hits as null. Empty input launches
  no process. Reject duplicate input IDs before writing scientific inputs.
- Group by exact oriented repeat and every interaction constraint. Each batch
  has one query FASTA record and at most 128 target records. Sort deterministically,
  use opaque IDs, request `id1,id2,start1,end1,start2,end2,E,hybridDP` and
  `--outNumber 1`, and retain the existing thread/process timeout limits.
- Parse and validate every row before choosing an interaction: required unique
  columns, expected IDs, finite energy, complete numeric fields, positive ordered
  in-bounds intervals, exactly one ampersand, target characters `.(`, repeat
  characters `.)`, structure lengths matching intervals, equal positive paired
  counts. Unknown IDs and malformed nonwinning rows are hard failures. A header
  with no rows is a valid no-hit result. Extra named columns are harmless.
- Select multiple valid rows by a documented total deterministic key beginning
  with energy; reject exact duplicate/inconsistent records as specified by the
  parser contract. Expose one-based paired positions on full oriented inputs.
  Preserve existing interaction field values and null semantics.
- Keep compatibility entry points where appropriate. Never silently retry failed
  batches or turn corrupt output into missing evidence. Pipeline failures publish
  failure artifacts and do not publish a complete report.
- Preserve standard candidate IDs, ranking, biological evidence and exports.
  Document raw evidence changing from per-candidate to per-batch files and the
  explicit IntaRNA 3.2.0 tested contract; no arbitrary version rejection.

## 3. Implement a bounded experimental detection branch

Owner: executor. Review owner: detection council.

Add `short_anti_repeat_rescue=false` and its explicit CLI switch. It applies to
Type II complete-mode detection. Reject incompatible bundled Type V configuration;
record user responsibility for family validity with custom models. Keep defaults
at identity 0.70, coverage 0.60, array masking enabled, CM-only rescue E <= 1e-5,
existing extension lengths and scoring weights. Models must not change.

The additional native genomic alignment proposals must have query span >=12 nt,
identity >=0.80, and coverage below the standard configured coverage gate. Exclude
accepted array overlap. Test both repeat orientations, admitting complementary
transcript hypotheses only. Reuse real genome-wide search output, never a known
reference RNA, truth coordinate, label, catalog entry or published repeat.

Evaluate an anchored interaction using installed IntaRNA's documented and tested
`seedTRange` contract. Group only identical constraints; individual short-branch
calls are acceptable given the small proposal volume. Admit a short candidate only
when all of the following hold:

1. Energy < 0, at least 10 paired target residues within `anti_local`, and at least
   half of that interval's target residues paired. A bounding span is insufficient.
2. Pair mapping uses ascending target and descending query positions; T is treated
   as U, and admitted pairs are AU, UA, GC, CG, GU or UG. This complementarity gate
   applies to experimental admission, preserving default native energy semantics.
3. At least one independent regional support: forward candidate-window CM with
   E <=0.01, starting at/after duplex start and overlapping duplex by <10 nt; OR
   forward ERPIN terminator E <=0.001 starting downstream of duplex end, fully
   within the oriented candidate window.

These are predeclared engineering hypotheses, not calibrated probabilities or
literature-established universal thresholds. Never relax them after seeing results.
Esp2/Bbo currently lack truth-localizing genomic seeds and may remain missed.

Record every proposal and acceptance/rejection reason, native supporting hits,
paired/anchored positions and support search context. Keep every standard candidate
available with its original ID. Add `admission_branch=standard|short_rescue|
genomic_cm_rescue` and descriptive evidence classes (`sequence_only`,
`sequence_duplex`, combinations with `cm`/`terminator`, and `cm_only`). Missing
evidence is not evidence of absence. Keep calibrated_probability false and search
windows distinct from transcript endpoints.

Do not implement locus consolidation in this round. It is independent of rescue
and batching, adds a new grouping/ranking contract, and is not needed to assess
either change. Also defer model training, Type V expansion, circular replicons,
new search backends and transcript-end inference. This resolves the detection
review's optional P2 grouping suggestion without weakening any raw-count gate.

## 4. Build isolated benchmark and audit tools

Owner: executor. Review owner: validation council.

- Add an explicit destination to the existing scorer or create a round-specific
  harness. Historical output paths must remain untouched. Refuse existing run dirs.
- Each run manifest records ordered input hashes, frozen protocol/source hashes,
  tool/model hashes, settings, environment, upstream revisions, workers, timestamps,
  actual batch wall time and per-input status. Use two loci at once and one native
  tool thread per locus for natural/control panels.
- Score with existing correct-contig/strand overlap >=20 nt AND >=50% shorter
  interval rule, fixed 78 denominator, errors as misses, Ghc2 null. Report all
  top-k ranks, MRR, gains/losses, >=90% coverage, highest-ranked match coverage,
  IoU and signed start/end errors; distinguish best-IoU from highest-ranked match.
- Report raw candidate counts, distinct coordinate/strand intervals, origins and
  admission classes, median/nearest-rank p95/max, per-Mb burden, runtimes and all
  rejection reasons. Do not claim a precision/FDR estimate from these positive loci.
- Audit source reconstruction on both strands, internal/public coordinate round
  trips, local and genomic evidence bounds, export ID/sequence/strand agreement,
  completion and evidence manifests, and expected source hashes.
- Compare every default candidate and rank with the frozen baseline, stripping only
  known additive fields and run metadata. Compare existing numeric fields exactly.
  Experimental results must preserve the original candidate projection, except
  rank shifts from inserting admitted extra hypotheses.

## 5. Execute the frozen validation protocol

Run sequential phases; retain actual logs, failures and source hashes.

1. Unit/contract tests, focused strand/admission/parser tests, Ruff and mypy;
   all seven existing integration cases must run with real tools, not skip.
2. Small native panel Ain, Kki, Lsp1, Mse, Cca1; repeat a small subset and compare
   normalized output. Run a complete reverse-complement stress case and explain
   any upstream asymmetry while enforcing strict coordinate reconstruction.
3. Controlled IntaRNA size-1 versus size-128 comparison on identical frozen Mse
   candidates, plus a multi-repeat case and Ssu repeated sequences. Use identical
   executable, thread count and host without concurrent unrelated benchmarks.
   Compare whole interaction maps, downstream evidence, calls and measured wall time.
4. Fresh complete runs on all 79 natural inputs: default, then experimental profile.
   All affected stages actually execute; do not substitute retained-hit rescoring.
5. Paired genomic stress controls on exactly these 19 inputs:
   Esp2, Bbo, Cca2, Cca1, Cga, Cme1, Ffr, Ghy2, Orh, Psp, Ain, Kki, Lsp1,
   Lan, Rsp, Cme3, Ghe, Mse, Sth1A. Preserve accepted array sequence intervals and
   ambiguous positions, deterministically shuffle non-array segments preserving
   length/composition; freeze FASTA hashes before running. Run both profiles on
   the identical bytes (38 native runs). Destroyed Cas context limits interpretation.
6. Conditional repeat-background controls on all 79 inputs with frozen arrays:
   two exact dinucleotide-preserving repeat shuffles and one cross-locus repeat
   derangement. Use per-source/array stable seeds (replicates 1 and 2); preserve
   length, ambiguous positions and dinucleotide counts; test invariants and flag
   unchanged/degenerate shuffles. Derangement excludes identical and reverse-
   complement-identical donor repeats and records donors. Execute new native search,
   interaction and supporting-evidence stages with source DNA fixed. These are
   explicitly conditional component controls, not end-to-end biological negatives.
   Report short-rescue admissions per input/array/Mb and all failed runs.

The primary end-to-end budget is 79 default +79 experimental +38 controls; small
and conditional runs are additional. Model-only sequence shuffles are omitted
because models and model thresholds are unchanged. Exact dinucleotide shuffling is
required only for the conditional repeat controls. Full genome controls use the
explicit composition-preserving scheme above, avoiding an unvalidated dinucleotide
claim. No timing or biological conclusions are inferred from failed controls.

## 6. Acceptance and independent council evaluation

Hard gates: all tests pass without unexpected skips; all natural runs complete;
all coordinate/export/provenance audits pass; protected files unchanged; no truth
leakage; deterministic controls; exact normalized default output equivalence.

Engineering acceptance: full interaction-map equality, exact pair routing, strict
malformed-row failures and reduced invocation count. Fresh Mse grouped interaction
time must be below serial; 20% reduction is a target, not a scientific correctness
substitute. Report actual whole-pipeline times separately.

Useful experimental sensitivity requires >=69/78 localized, all original 68 retained,
top1 >=53, top5 >=63, >=90%-coverage >=64; <=5,092 raw hypotheses, p95 <=288,
max <=928, no locus adding >100 short candidates, comparable batch time <=2x.
Each repeat-control panel requires median added short candidates <=1 and max <=10.
Genomic control burden rising on >half the controls, >1.5x aggregate, or a new
strongest-support prediction requires explicit council review and blocks default
promotion. Zero denominators remain explicit. No default promotion is planned.

A scientifically unsuccessful but correct, bounded, explicitly experimental feature
may remain disabled if the council accepts it and the negative result is prominent.
Never silently revise cutoffs to make this regression panel pass.

The same three council agents independently inspect implementation and evidence and
write `04_detection_evaluation.md`, `05_engineering_evaluation.md` and
`06_validation_evaluation.md`. Each reports pass/fail/not-measured per gate, concrete
findings, native evidence and limitations. Send blockers back to the executor; the
affected council member rechecks corrections. Coordinator verifies final tests,
source/protection hashes and records a final results report. Complete the goal only
after this full cycle, including honest unresolved scientific outcomes.
