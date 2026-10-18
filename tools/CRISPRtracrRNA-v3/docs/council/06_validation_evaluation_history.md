# Council validation evaluation: round 2

Status: **original historical exact-equivalence gate failed; native
reproducibility diagnosis and any prospectively amended comparison pending**.
Same validation/reproducibility reviewer as `03_validation_review.md`.
Review date: 2026-09-27. No native workloads were launched during this review;
implementation and retained run files were read-only. Findings were sent to the
coordinator and `council_executor` for correction without interrupting runs.

## Frozen contract and evidence

The governing contract is `docs/COUNCIL_IMPROVEMENT_PLAN.md`, its corresponding
`benchmarks/council_round2/protocol.json`, and pre-result
`protocol_clarifications.md`. The plan digest still equals
`4928592891d03d70ab06faae313c6b39f3248afa34774aa97a43c4b3d2e86a0d` and matches
the protocol. The controls refer to that same protocol. Clarifications preserve
thresholds and define review of **every** newly corroborated genomic-control
short admission; the evidence tier cannot be selected after viewing results.

Inspected the complete benchmark/audit scripts, focused control tests, frozen
control manifest, archived pre-native generator, small-panel reports/metrics,
coordinator audits, test logs, controlled performance results and reverse stress
comparison. Independent read-only checks are retained in
`docs/council/06_validation_preliminary_audit.json`.

The source generator extracted in memory from `implementation_pre_native.zip`
exactly matches the control manifest's generator digest:
`4171485205794e68db6ac0cf8a4945a616c0cd801acc6ccf04cdaa2eccee489f`.
The controls manifest digest is
`67efd2ede9053e2638de91c2646c08c157c347c13c833a8aa4a66b51725ebb55`.
The generator is preserved separately because later harmless harness/provenance
edits change the current script digest. This correctly preserves reproducibility.

## Findings and correction recheck

**V1 — Default-equivalence audit could hide unexpected extra candidates; fixed.**
At first review `council_audit.projection()` unconditionally excluded rows whose
`admission_branch` is `short_rescue`, even for `experimental=False`. A minimal
in-memory check comparing one standard candidate against the same candidate plus
an extra short-rescue candidate returns equality. This violates the exact-default
acceptance gate. Exclusion is appropriate only for the experimental projection.
The default audit should also explicitly reject any short proposals/admissions,
compare complete counts/order, and have a regression test for this example.

This is a defect in an audit, not evidence that default detection actually
changed. All five observed default cases have identical full candidate counts
and no short candidates; the independent coordinator comparison also examines
all original fields without this blind spot. Nevertheless the reusable audit
was corrected before final acceptance evidence. The revised function filters
only experimental projections and `report_membership()` rejects default
proposals/admissions explicitly. The reviewer re-ran the in-memory counterexample:
it now rejects the added default candidate. A regression test was added.

**V2 — Provenance cross-link enforcement; substantially fixed, final links pending.**
At first review panel provenance recorded source/native/model/protocol/plan hashes,
settings and input hashes, while the audit checks artifact hashes and individual
source reconstruction. It does not yet enforce all cross-links: immutable
plan/protocol agreement; exact expected full-panel membership; stage profile
against report configuration; report source IDs/files against the corresponding
run input; genomic input hashes against the frozen control manifest; conditional
source bytes/control repeats against that same manifest; and completeness before
interpreting conditional statistics. A set of internally consistent artifacts
from the wrong configuration or specimen must not pass the provenance gate.

The existing records contain enough evidence to add strict final checks without
rerunning native stages whose bindings prove correct. Preserve original run
manifests and source archives, including recorded versions preceding audit-only
fixes. Do not compare an old source hash against a newly edited script and
mistake that expected change for an invalid old run. Compare to the retained
source for the run and explicitly check scientific code/configuration consistency.

The executor added immutable protocol/plan/control bindings, exact panel
membership/profile rules, source/model/native checks, and report input/configuration
matching. The reviewer confirmed the correct small manifest passes and in-memory
mutations of profile, membership, protocol hash, input hash, threshold and added
default short rows all fail. Evidence is in
`06_validation_audit_correction_check.json`.

The executor subsequently added exact selected-panel repeat-map/array-count
matching, proposal/admission recounting, and the 79-completion requirement. These
changes were inspected; their completed conditional artifacts do not yet exist
to test. Retained per-run source versions and protected native provenance still
need their final archive linkage. Neither finding establishes corrupt scientific
inputs/results; current native runs can continue.

## Independent preliminary checks

| Gate or observation | Status | Evidence and scope |
| --- | --- | --- |
| Frozen plan/protocol/control protocol linkage | Pass for inspected files | Plan hash matches protocol; control manifest's protocol hash matches |
| Generator retained at recorded version | Pass | Archive member digest matches control manifest |
| Per-specimen shuffle reproducibility | Pass | Regenerated both 84-array dinucleotide panels from archived generator and recorded seeds; exact repeat bytes match |
| Dinucleotide and mononucleotide invariants | Pass | Independent exact length/base/adjacent-pair comparisons; coordinator additionally checks ambiguity positions |
| Cross-locus donor assignment | Pass | Recomputed all 84 assignments, including reversed donor-pool ordering; matches manifest; coordinator verifies donor exclusions and bijection |
| Array-preserving genomic controls | Pass | Reconstructed all 19 sequences from original sources and per-segment seeds; full bytes, source hashes, composition and accepted-array preservation match |
| Small default metrics | Pass for 5 loci | Independently recomputed overlap/strand/rank/coverage for Ain, Lsp1, Mse, Kki, Cca1 without importing the benchmark metric |
| Small default biological equivalence | Pass for 676 predictions | Full counts are 8, 4, 464, 8, 192; all old-field projections match and no short rows are present |
| Reusable exact-default audit | Pass after V1 correction | Counterexample now rejects; explicit default proposal/admission rejection |
| Strict final provenance audit | Partial pass; final V2 links pending | Positive small manifest and six negative mutation checks pass; conditional/source-version checks remain |
| Full default input completion and reference metrics | Pass | 79/79; 68/78 localized; top 1/3/5/10 = 53/61/63/64; 64 at >=90% coverage |
| Historical exact default prediction equivalence | Fail | Initial new panel has 3,397 candidates; historical 3,395. Lsp2 has 16 versus 14 |
| Experimental native localization/ranking/burden gates | Not measured | Full experimental panel pending |
| Paired 19 genomic controls / 3 conditional 79-input panels | Not measured | Generated inputs verified; native outcomes pending |
| Complete all-prediction coordinate/export audit | Not measured for full panels | Coordinator's small audit passes all 676 predictions and 324 evidence files |
| Final protected-file and upstream-state check | Not measured at end of execution | Coordinator's 14,902-file preservation manifest exists; final comparison remains due |

The known localization predicate remains unchanged, with contig and strand
required and inclusive overlap >=20 nt and >=50% of the shorter interval.
The implemented summary correctly assigns null detection/coverage to Ghc2 and
zero contributions for failed mapped loci. For full panels, the final audit
must explicitly establish all 79 expected labels and the 78 mapped denominator;
small exploratory panels correctly use their own declared subset. Ssa and Sth1A
remain different loci. MRR uses zero for missed references. Highest-ranked match
coverage/IoU/end errors and best-IoU are distinct. Nearest-rank p95 is correct.

A minor reporting detail: the current median implementation chooses the upper
middle for even-sized lists. Required completed panels have odd sizes (79 and
19), so their medians are unaffected. Do not interpret an even successful-only
subset of a failed conditional panel as its acceptance statistic. All expected
conditional inputs must complete before applying its median/max gates.

Genomic controls are not biological truth. Their per-row `detected`/`coverage`
fields initially initialized to false/zero. The executor corrected these to null
before any genomic-control runs; biological aggregate metrics remain suppressed.
Their native outcomes measure burden under a composition perturbation.

## Controlled performance and reverse stress

The retained native interaction comparison reports exact complete interaction-map
and old-field equality in all three cases:

| Locus | Serial → grouped calls | Serial → grouped wall time | Reduction |
| --- | ---: | ---: | ---: |
| Mse | 464 → 6 | 44.614 → 29.904 s | 33.0% |
| Ssu | 104 → 2 | 11.543 → 8.148 s | 29.4% |
| Sth1A | 128 → 2 | 15.045 → 11.204 s | 25.5% |

This meets the predeclared Mse speed target and invocation criterion as recorded;
engineering independently reviews native routing/maps. These are controlled
interaction-stage timings, not measured genome-wide speedups or replicated
statistical runtime estimates. The experimental <=2x whole-panel time gate still
requires comparable completed full-panel wall times.

Reverse-complement Ain is explicitly asymmetric at the upstream inference/search
level: the comparison contains 4 unique forward intervals versus 66 reverse-mapped
intervals, with 2 shared. The raw prediction counts are 8 and 132. The coordinator
traced an accepted-array extent change from 2,340 nt to 234 nt; array masking then
exposes different sequence. This is an observed limitation of the existing full
pipeline, not a coordinate-invariance pass. Reconstruction remains mandatory and
passes in the coordinator's audit. Do not state that complete biological output
is reverse-complement invariant or infer a new default regression from a property
the baseline upstream detector never guaranteed.

Initial retained tests show two test-fixture failures, subsequently fixed without
scientific threshold changes; the later full log reports 110 passing tests.
Subsequent code review fixes increased the unit total (105 pass, seven native
cases deselected in `unit-final.log`). Final source-version tests, all seven real
integration cases and lint/type results remain required after the audit fixes;
a unit-only rerun is not a substitute for final full-suite evidence.

## Final review still required

Finish V2 checks against conditional artifacts and retained source versions; V1 and six in-memory provenance
rejection cases have been rechecked. Then independently compare complete 79+79 results with the frozen baseline,
including all original detections, exact default output, experimental original
candidate preservation, gained/lost rank lists, IoU/length/burden summaries and
per-origin/tier counts. Check all 19 paired genomic controls use identical frozen
bytes and explicitly review every new corroborated short admission. Verify each
conditional panel covers 79 inputs/84 frozen arrays, contains native stage evidence,
and meets or fails its median<=1/max<=10 admitted-candidate limits without dropping
failures or unchanged shuffle controls. Finish with source/executable/model hashes,
protected-file checks, final tests and gate-level pass/fail/not-measured reporting.

A correct disabled experimental branch may fail scientific usefulness. That must
be reported prominently; neither this positive-only training-overlapping benchmark
nor any shuffled-control outcome supports a precision/FDR/generalization claim.

## Initial full-panel discrepancy and defensible acceptance path

The first full `natural_default` run completed all 79 inputs in 346.422 seconds.
Localization, all top-k counts, MRR and 90%-coverage counts equal the historical
baseline, but exact prediction equivalence **fails**: 3,397 versus 3,395 raw
candidates. Only Lsp2's candidate count differs, 16 versus 14. Aggregate agreement
does not override this failed predeclared hard gate. Preserve this run, its audit
failure and all historical files; do not replace it with a later matching rerun.

The reviewer independently inspected both Lsp2 native search directories. Source
FASTA, masked target FASTA and repeat-query FASTA bytes are identical. FASTA36
argv (`-T 1 -m 8 repeats.fasta anti-search.fasta`) and native executable digest
are also identical. Historical output has seven alignments, initial new output
eight. The extra row is the reverse-repeat query at genomic 19880–19905, query
26–1, identity 88.46%, span 26 nt. Other E-values/bitscores differ. This establishes
that the first relevant difference precedes IntaRNA batching. It does not alone
prove RNG as the cause; the executor/engineering native diagnosis must establish
the mechanism. Evidence is retained in `06_validation_native_discrepancy.json`.

Subsequent retained diagnostics strengthen the attribution to pre-existing native
behavior. `fasta_reproducibility/results.json` records six fixed-input native
repeats: all emit the same eight alignment keys, but all six full output hashes
differ because their statistical values differ. These repeats establish output
variability; they did **not** recreate the seven-hit historical sample, so they
must not be described as six observations of candidate-count fluctuation.
`archived_baseline_replay/comparison.json` records fresh execution of the archived
detector producing 16 Lsp2 candidates and 132 reverse-Ain candidates, with exact
projected biological equality to the new detector in both cases. This directly
shows both observed behaviors also occur in the archived implementation. These
limited diagnostic comparisons do not retroactively pass the original 79-locus
historical equality gate.

The apparent numerous candidate-ID/score differences must not be interpreted as
many independent batching errors: inserting a native hit changes sequential
candidate IDs and candidate-window search context. Nevertheless this causal
explanation is not a substitute for a paired exact comparison.

Recommended next steps, before any amended whole-panel results are inspected:

1. Keep two explicit ledgers: the **original historical gate remains failed**;
   a separately named amended native-reproducibility comparison may later pass.
   The final report must show both. Existing tests and scientific thresholds
   retain their original acceptance rules.
2. Establish native variability with bounded, fully retained diagnostic repeats
   of the identical Lsp2 native inputs under the archived command, and confirm
   the archived detector has the same exposure. Record executable, environment,
   command, outputs and source documentation/help explaining the RNG behavior.
   Diagnostics are exploratory explanation, not substitute accuracy runs. Do not
   rerun until an output happens to equal history.
3. Before running an amended full comparison, freeze an append-only amendment
   with exact seed mechanism, one seed, native flags/environment, source versions,
   expected run list, comparison fields and unchanged thresholds/burden limits.
   Prefer the already predeclared control master seed `20260927` if the verified
   native seed interface accepts it. Its value is not to be selected for matching
   historical outputs or improving recall. Do not edit the original plan/protocol.
4. Apply the identical deterministic-native protocol to both the archived
   pre-change detector and the new default detector. Extract the archived code
   into a new isolated location; do not edit protected originals. If an external
   wrapper/injection is necessary for the old version, retain and hash it, record
   the underlying executable, and show that it only supplies the declared seed.
   Report this as an instrumented archived-detector comparison, not an untouched
   historical rerun. Confirm repeat native outputs under that seed agree before
   treating the setting as deterministic; pin thread counts and dependencies.
5. Use exact paired original-field equality over all 79 amended default reports.
   Also use the same frozen candidate bank for size-1 versus grouped IntaRNA,
   including both observed Lsp2 candidate sets. Native RNG differences must not
   be removed by tolerances or by dropping unmatched candidates. Candidate-bank
   replay isolates the batching contract; it is not a new genome-sensitivity run.
6. Run the experimental detector and paired controls under the same amended
   deterministic-native protocol. Do not mix its comparator with the initial
   unseeded new run. Preserve the original absolute gate floors/caps (including
   68 original detections, 69 useful experimental localizations, top-k floors and
   5,092-candidate cap); report any amended baseline metric change separately.
   All controlled timings must identify the seed/protocol they actually used.

If no validated deterministic mechanism exists, fixed-input exact component
equivalence can still support the batching implementation, and the initial native
run can establish execution and observed reference metrics. It cannot satisfy
the historical exact full-pipeline gate. Keep that gate failed and the default
invariance claim qualified; scientific usefulness of the experimental branch
requires a separately predeclared repeated/paired design or remains not
established. It is acceptable to finish with an explicitly disabled experiment
and an honestly limited result. It is not acceptable to silently relax the gate.

## Amendment wiring review while archived panel runs

Reviewed `COUNCIL_REPRODUCIBILITY_AMENDMENT.md`,
`protocol_amendment_v2.json`, `native_execution_manifest.json`, dependency/build
freeze manifests and current archived-run/benchmark/audit/results scripts.
Independent evidence is in `06_validation_amendment_wiring.json`. This review
launched no native processes and did not inspect a completed amended accuracy
result. The original historical exact gate remains failed.

The actual frozen artifacts and active archived run pass these checks:

- All 390 paths in `dependency_freeze_sha256.json` retain their digests. The
  amendment/prebuild-clarification digests match the append-only amendment freeze.
- Native execution manifest links the amendment digest, binary digest
  `e53ea2de4b1584388da0841222b6254fab8f1ebd9e5b0863cd42f341132f348e`,
  configuration, master seed 20260927, native threads=1 and locus workers=2.
- All 40 extracted archived files match the original source checksum map.
  The active archived run records the same pinned binary/configuration and 79
  inputs; its imported Python module digests match the archived originals.
- Original control manifest remains byte-identical. The new reproducibility
  configuration differs semantically from protected `local.toml` only in its
  FASTA36 executable and thread count. The installed executable and original
  configuration remain protected.
- Coordinator's independent source audit records only the two declared original
  source-file changes plus the added typed seed header, with RNG recurrence source
  unchanged. Engineering owns the detailed C/build review. The pinned native test
  log reports 115 tests passed; final source-version checks remain due.

**V3 — Conditional benchmark ignores the selected configuration; blocker before
conditional runs.** `conditional(panel_name, config=None)` receives the CLI config
but loads `local.toml` unconditionally. Thus it would use stock FASTA36 even when
invoked with the pinned configuration. It must load the passed configuration and
its run must be rejected unless it declares the amended pinned dependency. This
was found before conditional controls ran; the existing archived run is unaffected.

**V4 — Amended provenance checks can be bypassed by missing metadata/profile
mislabeling; blocker before amended acceptance.** Pinned checks are conditional on
`native_dependency_profile` being present, while expected experimental status is
derived from the manifest itself rather than the declared stage. In an in-memory
negative test, the actual initial **stock-native default** manifest passes
`provenance_audit()` when supplied under either
`current_default_pinned_native` or `current_experimental_pinned_native`. Require
an explicit stage-to-profile map and mandatory amended/native/configuration
bindings for both current panels, genomic panels and every conditional panel.
Validate the archived run's different manifest schema against the same native
execution/dependency freeze, and validate the paired comparator's provenance as
well as its report contents. The fact that the active archived manifest is correct
does not make this permissive gate safe.

**V5 — Historical retention/dependency-delta reporting is incomplete; blocker for
final scientific attribution.** `council_results.py` labels pinned-default to
experimental retention as `all_original_detections_retained`. That is not the
required frozen historical set of 68 references. Build and compare the historical
68-label set explicitly against archived-pinned, current-default-pinned and
current-experimental-pinned results. Retain the numerical floors and absolute
caps unchanged. Add historical-to-archived-pinned gained/lost labels, ranks and
burden as the dependency delta; retain archived-to-current paired equality and
current-default-to-experimental gains as separate comparisons. Also correct the
outdated `audit_natural_experimental.json` filename to the actual amended stage's
audit. No baseline substitution may conceal loss of an original reference.

All three findings were sent to coordinator and executor before their affected
runs/acceptance claims. They are harness/audit/reporting defects; the inspected
active archived native execution remains correctly configured. Required recheck:
conditional CLI-config wiring test, rejection of stock/default manifests under
both amended stage names, archived/current same-dependency binding, and a synthetic
historical-loss example where the pinned default and experiment both miss an
originally detected locus. Final empirical verdict awaits completed panels.

### Amendment correction recheck: proceed

**V3, V4 and V5 are corrected. Proceed with the amended native panels.** This is
approval of the corrected execution/audit wiring, not a final empirical pass.

The reviewer independently ran `tests/test_council_controls.py` with bytecode and
pytest cache writes disabled: **nine tests passed in 0.06 seconds**. The direct
conditional-call test confirms the requested configuration reaches the shared
pinned-settings loader; actual selected settings resolve the frozen binary with
one thread. Stock configuration is rejected for a pinned run.

Independent in-memory negative checks now reject the original stock-native
manifest under all seven amended natural, genomic and conditional stage names.
Valid prospective pinned default and experimental manifests pass. Removing
pinned metadata or substituting a wrong binary, thread count or profile rejects
in both natural stages. The dedicated archived auditor now checks the archived
source/resource digests, pinned execution/configuration, exact ordered 79 inputs
and completed report membership; its final execution depends on completion of
the still-running archived panel.

Historical retention is now separate from pinned-default retention. A synthetic
example with one original detection lost and one replacement gained—equal total
counts, with identical pinned default and experimental results—correctly reports
the original loss in archived, current and experimental stages. Dependency deltas
and experimental gains are written separately, and the experimental audit filename
matches its amended stage. Original floors/caps remain unchanged.

All 390 dependency-freeze files remain unchanged. Evidence and auditor digest are
in `06_validation_amendment_correction_check.json`. Original historical equality
remains **FAIL**; final amended equality, control burden and experimental usefulness
remain **not measured** until their complete native evidence is audited.
