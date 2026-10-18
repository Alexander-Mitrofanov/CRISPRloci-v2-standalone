# Council validation evaluation: final round 2

**Verdict: accept the amended dependency-controlled implementation and the useful,
disabled experimental profile. No unresolved validation blocker remains. The
original stock-native historical exact-equivalence gate remains FAIL.**

Reviewer: the same validation/reproducibility council member who wrote
`03_validation_review.md`. Review date: 2026-09-27. This reviewer ran read-only
artifact checks and nine focused software tests, but launched no native scientific
workloads or implementation changes. Earlier findings, corrections and protocol
decisions remain in `06_validation_evaluation_history.md`; none is erased by this
final verdict.

## Contract and comparison units

The original contract is `docs/COUNCIL_IMPROVEMENT_PLAN.md` and
`benchmarks/council_round2/protocol.json`, with pre-result clarifications. The
append-only amendment is `docs/COUNCIL_REPRODUCIBILITY_AMENDMENT.md` and
`protocol_amendment_v2.json`. Original scientific floors, the historical set of
68 detected references, absolute candidate caps and controls were retained.

The natural panel contains 79 attempted inputs and **78 uniquely mapped reference
loci**. Ghc2 remains biologically unassessed/null. Ssa and Sth1A remain separate
natural loci despite their duplicated reference sequence. Correct contig and
strand, overlap >=20 nt and >=50% of the shorter interval define localization.
Execution errors would count as non-detections; all completed natural panels have
79 successes. These are previously inspected development/regression loci with
training-source overlap, not an independent held-out accuracy study.

The four relevant comparisons are deliberately distinct:

| Comparison | Result |
| --- | --- |
| Historical stock-native outputs versus first new stock-native default | **FAIL** exact equality: 3,395 versus 3,397 candidates; Lsp2 14 versus 16 |
| Historical outputs versus archived detector using pinned native dependency | Dependency changes retained explicitly: Lpn 12→10, Lsp2 14→16; all reference localization/rank totals unchanged |
| Archived versus current default, identical pinned dependency | **PASS** exact original-field prediction/evidence/ranking equality across all 79 inputs and 3,395 candidates |
| Current default versus experimental, identical pinned dependency | **PASS** all standard candidates/evidence retained; 11 admitted hypotheses add seven localized references |

The initial historical discrepancy begins in FASTA36 before IntaRNA: masked DNA,
repeat FASTA, native argv and binary hashes match, but historical output has seven
alignments versus eight. Six retained unseeded repetitions have different
statistical output hashes and the same eight alignment keys; those repeats do
not independently reproduce the seven-hit sample. Archived-detector replay also
produces the new 16-candidate Lsp2 and 132-candidate reverse-Ain outcomes. Aggregate
stability therefore never substitutes for the failed original equality gate.

## Independent evidence and provenance

The reviewer verified all of the following directly:

- The amended native execution manifest fixes master seed 20260927, one native
  thread, two locus workers, and separately named FASTA36 binary digest
  `e53ea2de4b1584388da0841222b6254fab8f1ebd9e5b0863cd42f341132f348e`.
  The original binary, local configuration, upstream repositories and models were
  preserved. Only executable selection and thread count differ semantically in
  `local.reproducible.toml`.
- All 40 extracted original archived files match the baseline checksum map.
  Each archived run's 12 Python-source digests matches the baseline archive.
  Current default, experimental and every control run's 13 source/harness digests
  match the retained source archive
  `current_default_pinned_native_source.zip`, digest
  `126cfe22a36523e87d2ef217febd38fbb35626630f3a8cdf1a11ebf88ce422fe`.
- All three natural panels have the exact expected 79 input identities/hashes,
  completed results, expected settings, pinned executable and model resources.
  Per-record source bytes, reported sequences/coordinates and local/genomic
  evidence round trips were checked independently. CSV/FASTA output IDs, sequence,
  strand and positions agree; coordinator/auditor checks additionally cover GFF.
- Independently rehashed **23,829** evidence files across archived, current-default
  and experimental natural runs: 13,913 + 4,655 + 5,261, all matching.
  Independently checked another **6,614** evidence-file hashes across completed
  genomic and conditional controls, all matching.
- Rehashed **14,902 protected files** and **390 dependency-freeze files** without
  a mismatch. `old_version`, `external/CRISPRidentify-v2` and `external/CasAndra`
  each return empty `git status --short` with exit code zero.
- Native seed repeatability artifacts contain six identical Lsp2 output hashes
  and two identical outputs each for Ain, Kki and Mse, under path/launch-order
  changes with the same ordered input bytes. This validates the declared
  single-thread profile, not arbitrary thread counts or cross-platform behavior.

Engineering independently owns the exact C patch/build review. Its and the
coordinator's source audit restrict modifications to the declared typed RNG
header and two seed/declaration sites; RNG recurrence and scientific scoring
remain unchanged. Preserved build attempts, source archive, patch, binary and
licenses form part of the frozen dependency evidence.

## Natural results and frozen acceptance gates

An independent standard-library recount, without importing the harness matching
function, produced these results from the actual reports:

| Metric | Pinned archived/current default | Experimental |
| --- | ---: | ---: |
| Completed inputs | 79/79 | 79/79 |
| Localized references | 68/78 | 75/78 |
| Top 1 / 3 / 5 / 10 | 53 / 61 / 63 / 64 | 60 / 68 / 70 / 71 |
| Mean reciprocal rank; misses=0 | 0.7283666082 | 0.8181101980 |
| At least 90% reference coverage | 64/78 | 68/78 |
| Raw candidates | 3,395 | 3,406 |
| Candidates: median / nearest-rank p95 / maximum | 10 / 192 / 464 | 10 / 193 / 464 |
| Candidate length: median / p95 / maximum | 191 / 202 / 204 nt | 191 / 202 / 204 nt |
| Current whole-panel wall time | 344.020 s | 347.170 s |

The seven experimental gains are **Cca1, Cga, Cme1, Ffr, Ghy2, Orh and Psp**.
There are no lost detections, and every original historical detection is retained
in archived-pinned, current-default-pinned and experimental-pinned stages. Esp2,
Bbo and Cca2 remain undetected; Ghc2 remains unassessed. The pinned dependency itself
adds no localized reference and changes no reference rank. Its two candidate-count
changes offset in the total; they reduce the median from the historical 12 to 10.

The experimental branch admits 11 of 66 proposals: two with qualifying CM support,
four with qualifying terminator support and five with both. Its original-field
standard-candidate projection is exact; rank-derived export IDs may change on
insertion, as declared before results. Unchanged candidate-length distributions
show that this measured gain did not come from broadening all windows. Coverage
and localization still do not establish transcript endpoints.

| Predeclared gate | Verdict | Evidence |
| --- | --- | --- |
| Original stock-native historical exact equality | **Fail, retained** | Lsp2 discrepancy and initial run preserved |
| Amended default exact equality | Pass | 79 paired reports, 3,395 predictions; full original-field equality |
| Full natural execution/source/coordinate/export audits | Pass | All 237 amended natural runs complete; source/evidence checks above |
| Historical 68 retained; useful experimental localization >=69 | Pass | All 68 retained; experimental 75 |
| Experimental top1>=53, top5>=63, coverage90>=64 | Pass | 60, 70 and 68 |
| Raw candidates<=5,092; p95<=288; maximum<=928 | Pass | 3,406; 193; 464 |
| Added short candidates per locus<=100 | Pass | Maximum one |
| Comparable experimental/default wall ratio<=2 | Pass | 1.0091544643; same pinned settings/resources |
| Exact frozen-candidate interaction replay | Pass | 3,390 exact old interaction maps across 79 loci, 152 calls |
| Mse grouped native interaction faster | Pass | 44.614→29.904 s; 464→6 calls; exact maps |
| Genomic control completion/burden | Pass | 19+19 successes; 60 candidates each; no new admission/burden increase |
| Conditional control completion/median/max limits | Pass | All 3×79 inputs; per-panel median 0, maxima 1/0/1 |
| Control generation/frozen bytes/source provenance | Pass | Exact regeneration and manifest/source/archive checks |
| Protected files/upstream state | Pass | 14,902 + 390 digests unchanged; three clean repositories |
| Final software validation | Pass | 119 tests including all seven native cases; Ruff and strict mypy pass |
| Independent biological precision/FDR/generalization | **Not measured** | Current data/control design cannot support these estimates |

The archive versus current default wall times are 390.527 and 344.020 seconds,
one paired observation, not a replicated runtime estimate. The controlled Mse,
Ssu and Sth1A interaction comparisons improve by 33.0%, 29.4% and 25.5%, with exact
complete maps. Do not substitute these stage timings for genome-wide speed claims.
The full 3,390-interaction replay uses frozen historical candidate inputs and is
explicitly a component equivalence check, not genomic discovery.

## Controls and interpretation

Controls retain manifest digest
`67efd2ede9053e2638de91c2646c08c157c347c13c833a8aa4a66b51725ebb55`.
The original generator version was independently recovered from its archive and
regenerated: two panels of 84 exact-dinucleotide repeat shuffles, 84 cross-locus
repeat assignments and 19 array-preserving genomic compositions all match their
frozen bytes/seeds. Adding or reordering unrelated shuffle inputs does not alter
per-specimen seeds; the derangement donor pool itself is frozen.

Both genomic profiles use identical frozen input bytes. All 38 complete-mode
runs succeed. Default and experimental each return 60 hypotheses; the additional
30 short proposals yield zero admissions, and no individual locus has increased
candidate burden. Biological outcome fields are null, appropriately: these are
composition controls, not validated biological negatives.

| Conditional repeat panel | Complete inputs / frozen arrays | Short proposals | Admitted | Median / max admitted per input |
| --- | ---: | ---: | ---: | ---: |
| Dinucleotide replicate 1 | 79 / 84 | 126 | 2 | 0 / 1 |
| Dinucleotide replicate 2 | 79 / 84 | 79 | 0 | 0 / 0 |
| Cross-locus derangement | 79 / 84 | 100 | 1 | 0 / 1 |

All conditional panels actually run the new native search/interaction/support
stages with the pinned dependency; per-array repeats and original DNA match the
frozen manifests. Counts are recomputed from proposal admission flags, with no
failed loci silently removed. The three admitted controls remain visible:
Tpu (CM) and Eit (terminator) in replicate 1, and Cga (CM) in the derangement.
Their presence is consistent with a descriptive stress experiment; calling them
biological false positives or estimating FDR from these counts would be invalid.
The original median<=1/max<=10 operational limits pass.

The admitted deranged Cga control exposes a limitation of broad-window
localization: its minus-strand window 833–1001 overlaps the reference 810–911
by 79 nt and passes the declared localization rule, but its native seed 983–996
lies entirely outside that reference and its supporting CM hit 910–937 overlaps
by only 2 nt. Window localization alone therefore does not establish the correct
anti-repeat or biological mechanism. In contrast, all seven natural-reference
gains have native seeds overlapping their corresponding references, independently
checked from the retained reports. This is useful supporting geometry, still not
proof of transcript endpoints or biological specificity.

These panels inspect composition dependence and candidate burden. Genomic
shuffling disrupts biological context; conditional controls retain inferred arrays
and native loci. Neither establishes precision or justifies promoting the
experimental branch to default. The experimental setting should remain disabled
unless explicitly selected, with training overlap and boundary uncertainty stated.

## Review fixes and remaining limitations

Five concrete audit/harness/reporting problems were found and corrected during
this same-council cycle: default projections could hide extra short candidates;
provenance links were incomplete; conditional controls ignored explicit config;
amended profiles could accept stock-native metadata; and historical retention
was initially tested against the new default only. Negative tests now reject the
original counterexamples, including a replacement detection that preserves total
count while losing a historical reference. Nine focused tests passed in the
reviewer's independent rerun. Their fixes changed no scientific threshold or
control bytes, and affected conditional runs had not started.

Final retained software logs report 119 tests passed in 18.84 seconds, including
all seven native integration cases, with no skips; Ruff and strict mypy pass.
`executor_final_validation.json` binds these logs, final audit/result files and
the final implementation archive. `final_stage_source_bindings.json` preserves
the exact archived/current source versions for every scientific panel. The
coordinator's final binding check confirms all 390 dependencies, six stock native
executables, three amendment files, 32 final archive files and the archive digest,
plus result/audit/log hashes. Its final preservation check at
`protection_verification_20260927T000259.json` again finds all 14,902 protected
files and upstream repository states unchanged.
Source versions used by scientific runs remain archived even where later audit
or reporting code was corrected. The three missed references, original stock
nondeterminism and reverse-complement upstream array asymmetry remain documented
limitations, not hidden test exclusions.

The three control admissions and seven natural-reference gains do not turn
this into a held-out study. The useful conclusion is bounded: exact batching
under paired native inputs, a reproducible selected single-thread dependency
profile, and improved regression localization with small measured candidate
burden from an explicitly experimental branch. No precision, calibrated confidence,
FDR, transcript-end accuracy or independent-generalization claim is supported.

## Independent verification artifacts

- `06_validation_natural_check.py/.json`: independent natural matching, paired
  projection, source-archive, reconstruction/export and evidence-hash checks.
- `06_validation_control_check.py/.json`: complete control input/repeat/provenance,
  source-archive, evidence-hash, admission recount and burden checks.
- `06_validation_protection_check.json`: final protected/dependency checksums.
- `06_validation_final_checks.json`: frozen interaction replay, seed-repeatability,
  repository status and final software-log evidence.
- `06_validation_localization_limitation.json`: native-seed/reference overlaps
  for all seven natural gains and the Cga derangement counterexample.
- Earlier `06_validation_*audit*.json`, `*_correction_check.json`,
  `*_native_discrepancy.json` and `*_amendment_wiring.json` retain independent
  findings and counterexamples before correction.

These final results can be reproduced from retained manifests and raw outputs
without rewriting historical `benchmarks/results` or the original baseline.
