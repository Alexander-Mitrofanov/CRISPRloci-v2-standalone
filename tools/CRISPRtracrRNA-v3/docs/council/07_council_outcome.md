# Council round two: implementation and evaluation

Date: 2026-09-27. Coordinator: root.

The council designed a bounded improvement round, a separate executor implemented
it, and the same three council members independently evaluated the result. The
[detailed plan](../COUNCIL_IMPROVEMENT_PLAN.md) was frozen before implementation.
The [reproducibility amendment](../COUNCIL_REPRODUCIBILITY_AMENDMENT.md) was frozen
before its affected runs. The original historical comparison failure is retained.

**Final decision: accept the engineering changes and the bounded experimental
detection feature under the documented reproducibility amendment.** The same
detection, engineering and validation reviewers each report acceptance with no
unresolved blocker within this scope. Implementation, native panels, controls,
tests and coordinator preservation checks are complete. The experimental feature
remains disabled by default; the original historical-equivalence gate remains
failed, separately from the successful amended paired comparison.

## Implemented changes

- **Faster interaction evaluation.** IntaRNA batches candidates sharing the exact
  oriented repeat, with a default limit of 128 targets. Every expected candidate
  remains represented, including no-hits. Native rows receive strict ID, numeric,
  coordinate and structure validation before deterministic result selection.
  Size 1 remains an explicit serial comparison path. The compatibility wrapper
  retains evidence in isolated directories across repeated calls.
- **Experimental short anti-repeat detection.** A default-off Type II branch
  tests short native genomic alignments with anchored RNA interaction and separate
  regional model or terminator support. It records all proposals, rejection
  reasons, paired residues and support context. It does not read reference
  sequences, truth coordinates or literature labels during detection. Existing
  thresholds, models and default biological predictions remain unchanged under
  the same dependency profile.
- **More reliable evidence and benchmarking.** Added explicit admission/evidence
  categories, corrected Python-module dependency provenance, and bound benchmark
  outputs to input, code, control-generator, configuration and executable hashes.
  The harness checks candidate preservation, exports, geometry and completion.
- **A separate reproducible native profile.** A minimal FASTA36 source repair
  supplies explicit random seeds. The separately named executable and
  `local.reproducible.toml` preserve the stock executable and `local.toml`.
  Reproducibility was tested with one native thread on this host; arbitrary
  multithread or cross-platform equivalence is not established.

## Natural-locus results

All 79 inputs completed in each paired panel; 78 references map uniquely. Ghc2
remains unmapped and is excluded from the localization denominator. These are
inspected development/regression loci with training-source overlap.

| Measurement | Current default | Experimental short rescue |
| --- | ---: | ---: |
| Localized references | 68/78 | 75/78 |
| Localized at rank 1 | 53/78 | 60/78 |
| Localized within top 3 | 61/78 | 68/78 |
| Localized within top 5 | 63/78 | 70/78 |
| Localized within top 10 | 64/78 | 71/78 |
| At least 90% reference coverage | 64/78 | 68/78 |
| Raw candidate hypotheses | 3,395 | 3,406 |
| Median / 95th percentile / maximum per input | 10 / 192 / 464 | 10 / 193 / 464 |
| Observed whole-panel wall time | 344.02 s | 347.17 s |

Seven previously missed references are localized: **Cca1, Cga, Cme1, Ffr, Ghy2,
Orh and Psp**. Each new match ranks first and has a native short seed overlapping
the mapped reference. Every historical detection is retained. Eleven of 66
proposals are admitted; three overlap references already localized, and one is
an unannotated Eit hypothesis. No input gains more than one candidate.

The remaining mapped misses are **Esp2 and Bbo**, which lack a truth-localizing
native short proposal, and **Cca2**, whose anchored duplex lacks qualifying
regional support. Its terminator E-value of 0.0236 exceeds the frozen 0.001 gate.
The thresholds were not relaxed to fit these references.

Localization requires the correct contig and strand, at least 20 overlapping
bases, and overlap covering at least half of the shorter interval. A matching
search window does not establish a full transcript, biological function or
calibrated precision. Detailed overlap, endpoint and coverage results are retained
in [experimental_geometry.json](../../benchmarks/council_round2/experimental_geometry.json).

## Engineering and reproducibility results

All **3,390 frozen interaction maps** preserve their original fields exactly,
using 152 batched calls. In the controlled Mse stress case, calls decrease from
464 to 6 and interaction wall time from 44.61 to 29.90 seconds, a **33.0% reduction**.
This is an interaction-stage measurement on that case. The archived/current
whole default panels took 390.53/344.02 seconds, an observed 11.9% reduction in
one paired measurement, not a universal speedup claim.

The first new run using stock FASTA36 **failed exact historical equivalence**:
Lsp2 had 16 candidates rather than 14, with the difference beginning in native
FASTA36 output before batching. Archived Python replay reproduced the new count.
Native source inspection found inconsistent seed initialization. Repeated stock
runs showed changing shuffled statistics, but did not reproduce the historical
seven-alignment sample. This original failure remains a failure.

The amended comparison runs the archived and current Python detectors using the
same separately built, seeded executable. Their default biological predictions,
evidence and ranks match exactly across all 79 inputs and 3,395 candidates. Against
historical output, dependency-only changes remain at Lpn (12 to 10 candidates)
and Lsp2 (14 to 16); reference localization and rank totals are unchanged. Equal
aggregate counts do not erase those differences. Six Lsp2 repeats and two each
for Ain, Kki and Mse produced identical native output bytes under the tested profile.

## Controls and limits

| Control panel | Completed inputs | Short proposals | Short admissions |
| --- | ---: | ---: | ---: |
| Genomic shuffle, default | 19/19 | Not enabled | 0 |
| Same genomic shuffle, experimental | 19/19 | 30 | 0 |
| Conditional dinucleotide shuffle 1 | 79/79 | 126 | 2 |
| Conditional dinucleotide shuffle 2 | 79/79 | 79 | 0 |
| Conditional repeat derangement | 79/79 | 100 | 1 |

The paired genomic controls have 60 hypotheses in both profiles. Conditional
panel medians are zero and maxima are 1, 0 and 1, within the frozen burden limits.
All shuffle invariants and native-stage provenance checks pass. These synthetic
controls do **not** provide a biological false-positive rate or specificity estimate.

The shuffle-1 admissions are Eit and Tpu; the derangement admission is Cga.
Shuffled Eit overlaps the same unannotated terminator region as the natural extra
Eit candidate. In deranged Cga, a broad admitted window overlaps the reference,
while its short seed lies outside it and its supporting model overlaps by only
two bases. This demonstrates a limitation of window-based localization, even
when the admission rules are satisfied. The seven natural gains were separately
checked for truth-overlapping native seeds. Independent biological validation and
more complete negative annotation are needed before default promotion.

## Validation and preservation

The final suite passes **119 tests**, including all seven real-tool cases with
no skips. Ruff and mypy pass. The seven final stage audits pass. Independent
reviews also reconstruct candidate sequences and evidence coordinates, check
exports and paired-base identities, verify all original detections, regenerate
controls and verify native evidence hashes.

Coordinator checks confirm that all **14,902 protected files** and protected
repository states are unchanged. All 390 frozen native/profile files, six stock
executables, three amendment files and 32 final archived implementation files
match their recorded hashes. Test logs, final metrics, audit results and the
implementation archive are bound to their recorded hashes as well. The final
implementation archive captures the tested source; later documentation-only
clarifications to README and IMPROVEMENTS are recorded separately. Final handoff
verification confirms all 32 archived members, the 30 unchanged live files and
the two explicit prose deltas; no code, harness or test changed after the suite.

## Council records and reproducibility

- [Detection design review](01_detection_review.md) and [execution evaluation](04_detection_evaluation.md).
- [Engineering design review](02_engineering_review.md) and [execution evaluation](05_engineering_evaluation.md).
- [Validation design review](03_validation_review.md) and [execution evaluation](06_validation_evaluation.md).
- [Final metrics and gate outcomes](../../benchmarks/council_round2/final_metrics.json).
- [Executor test and evidence record](../../benchmarks/council_round2/executor_final_validation.json).
- [Coordinator final hash checks](../../benchmarks/council_round2/coordinator_final_bindings.json).
- [Final handoff and documentation-delta checks](../../benchmarks/council_round2/coordinator_handoff_bindings.json).
- [Protected historical-file and repository verification](../../benchmarks/council_round2/protection_verification_20260927T000259.json).
- [Remaining-miss diagnosis](../../benchmarks/council_round2/remaining_miss_diagnosis.json).
- [Stage-to-source archive bindings](../../benchmarks/council_round2/final_stage_source_bindings.json).

From this checkout, use a fresh output path to try the experimental profile:

```bash
python CRISPRtracrRNA.py run --config local.reproducible.toml \
  --short-anti-repeat-rescue --input genome.fasta --output runs/short-rescue
```

Omit `--short-anti-repeat-rescue` for the default biological detection profile.
The experimental option remains disabled by default. The next scientific work
should test independently curated native references and better negative
annotations, then investigate seed generation for Esp2/Bbo and support/model
coverage for Cca2 without tuning on this same collection.
