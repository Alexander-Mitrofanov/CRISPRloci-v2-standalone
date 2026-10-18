# Engineering council evaluation — round 2

Date: 2026-09-27. Reviewer: the same `council_engineering` agent that wrote the phase-1 review. **Final verdict: ACCEPT the engineering implementation and the separately declared pinned-native, one-worker dependency profile. No unresolved engineering blocker remains in the reviewed scope.** The experimental branch remains default-off; this verdict does not promote it or establish biological specificity.

The original stock-FASTA full historical-equivalence gate is permanently **FAIL** (Lsp2 14→16). The amended comparison is separately **PASS**: archived and current Python detectors using the same pinned binary produce exactly equal existing prediction fields/order on all 79 loci. The failed run, native diagnosis, source amendment, and all later evidence remain retained.

## Reviewed contract and code

I reread `docs/COUNCIL_IMPROVEMENT_PLAN.md`, `benchmarks/council_round2/protocol.json`, and `protocol_clarifications.md`. I compared the current search, pipeline, configuration, CLI and process code with the frozen `baseline/source.zip`, and inspected the new focused tests and benchmark performance function. I made no implementation changes and launched no native scientific work during the executor's controlled timing.

The final reviewed core snapshot has these SHA-256 values. These match the retained current-stage manifests and live files; the wrapper correction is included:

| File under `src/crisprtracrrna/` | SHA-256 |
| --- | --- |
| `search.py` | `b79a58d279f008c0a561e2c7de37057d1c22e0ae444785b39a7b304020fc7db0` |
| `pipeline.py` | `6d3608a7ca94216b2fbd6375566ea324cf5118ceee70683a36c8c46e0537f583` |
| `config.py` | `af9eceb2af57a02b83e4a2df3b8922e22d650f076197178835d4841509e4c48d` |
| `cli.py` | `35befc0806315d7745955e47093450ff982fe9230cb262e4ec7fc957fbb8ef15` |
| `process.py` | `c3eb0f9e1cad0dadfe081d97b262dd46bed84b31728380dbcbb75892a8561ec6` |

## Gate status

| Engineering gate | Status | Independent evidence or remaining requirement |
| --- | --- | --- |
| Exact intended target/repeat pairs; configured chunk limit | Pass | Fake runner and independent 676-candidate / 14-batch routing audit; final all-79 native replay preserves all 3,390 interactions in 152 calls |
| Every expected candidate retained, no-hit as null, empty input launches nothing | Pass | Focused tests and offline historical parser check, including 457 no-hits |
| Strict all-row validation, including malformed nonwinning rows | Pass | Invalid IDs, columns, coordinates, nonfinite energy, structure characters/length/count, duplicates and incomplete rows are rejected; valid extra named columns accepted |
| Deterministic tie selection and stable routing | Pass for inspected scope | Total key includes energy, both coordinate intervals and structure; group/repeat sorting and candidate sorting are explicit; reversed native rows tested |
| Paired-position mapping and singleton anchored constraints | Pass for engineering contract | Positive equal paired counts; one-based positions relative to full oriented inputs; anchored calls pass each candidate's own public anti-local range and use exactly one target |
| Standard scientific inputs and evidence preserved | Pass for fixed-input batching and amended paired pipeline | Independent all-79 map comparison: 3,390 exact interactions; pinned archived/current default: 3,395 exact ordered predictions |
| Configuration/CLI backward compatibility | Pass for new settings and existing local configuration | Default 128 / false; explicit size 1 and boolean enable/disable work; invalid sizes rejected; experimental incompatible modes rejected |
| Failed scientific stage never publishes complete report | Pass | Injected batch-contract failure test; existing runner/tool failure and timeout tests pass; no silent retry path |
| Repeated legacy `interaction()` calls retain every raw result | Pass after correction | Unique private per-call directories prevent collisions, including repeated/unsafe-looking IDs; focused regression independently passed |
| Configured CasAndra module provenance correction | Pass for current configured absolute PYTHONPATH | Actual source file and expected CasAndra revision recorded separately from interpreter checkout; unresolved module stays explicit |
| Controlled native size-1 / size-128 evidence equality and timing | Pass for three required fixtures | Independently loaded all maps and command logs: exact size1/128/historical interaction equality, same native hash and one thread; Mse stage time reduced 33.0% |
| Original stock-native historical equivalence | **FAIL, retained** | First run: 79 complete reports, 3,397 predictions; Lsp2 historical14/fresh16 differs upstream of batching; the amendment does not rewrite this gate |
| Seeded native byte repeatability | Pass for declared one-worker scope | Independently rehashed all 12 probes: Lsp2 six; Ain/Kki/Mse two each, all exact within locus |
| Amended pinned-native archived/current default equivalence | Pass | Independent direct ordered prediction projection on all 79 loci, 3,395 predictions; same seeded native dependency |
| Experimental preservation and control completion | Pass | Independent comparison retains all existing predictions/order/evidence except expected export renumbering; 11 additions; genomic 19+19 complete, 60→60/zero short; conditional 79 each, admissions 2/0/1 |
| Final tests, lint, types, source/artifact bindings | Pass | 119/119 tests including seven native cases, Ruff and strict mypy (12 files); retained logs rehashed; seven current stages share identical 13-file core/harness source |
| Protected trees/dependencies | Pass by independent validation/coordinator audit | Their final rehash verifies 14,902 protected and 390 dependency files unchanged; this reviewer did not duplicate the full filesystem traversal |
| Arbitrary native thread counts, platforms, import modes; biological generalization | Not measured | Outside the accepted bounded profile; limitations below remain explicit |

## Independent executed checks

### Focused tests

From `new_version`, I ran:

```text
pytest -q tests/test_council_round2.py tests/test_scientific_parsers.py \
  tests/test_process_and_reports.py tests/test_sequences.py -m 'not integration'
```

Result: **72 passed in 0.73 seconds**, with no failures or skips. This includes the process timeout/descendant cleanup regression; no native scientific pipeline was started by this reviewer. The executor's retained `tests-pass1.log` additionally reports **110 passed in 14.47 seconds**; the executor says these include all seven native integration cases. The 72-test result is my independent execution, while the 110-test result is retained executor evidence. The final post-correction log has since been inspected and rehashed: 119 passed, including the seven native cases, in 18.84 seconds. These full-suite results are executor-run evidence, distinct from this reviewer’s focused executions.

### All historical native IntaRNA output files

I traversed the 79 historical `final_fixed/complete/*/work/*/logs/commands.json` files, selected every IntaRNA command, and recovered its exact target length, query length, output CSV path and candidate ID. I parsed the original native CSV using the **new** strict `parse_intarna()` and compared its result with that candidate's frozen `council_round2/baseline/reports/*.json` interaction.

Only the two additive fields `target_paired_positions` and `repeat_paired_positions` were removed from each new non-null result. The comparison covered **3,390 interactions: 2,933 nonempty and 457 null**. **Zero existing-field differences and zero parser errors** occurred. This shows that strengthened structure/row checks accept every retained real interaction while preserving its old coordinates, energy and structure. It is an offline compatibility check, not a substitute for new native batch equivalence.

### Completed small native panel artifact audit

I independently read each completed `small_default/{label}/report.json` and its frozen baseline report. I compared the entire ordered prediction lists after removing only `admission_branch`, `evidence_class`, `target_paired_positions`, and `repeat_paired_positions` from the new predictions. I did not remove coordinates, IDs, export IDs, ranking features/contributions/scores, old evidence, or candidate ordering.

For every `.ids.json` batch map I also read the corresponding target and query FASTAs. I required one query, at most 128 targets, exact agreement between target IDs and map keys, exact candidate target sequence, exact oriented repeat sequence, and a one-to-one multiset match between routed and expected candidate IDs.

| Locus | Standard anti-repeat candidates | Batch invocations | Existing prediction fields/order equal | Exact pair routing |
| --- | ---: | ---: | --- | --- |
| Ain | 8 | 2 | Yes | Yes |
| Kki | 8 | 2 | Yes | Yes |
| Lsp1 | 4 | 2 | Yes | Yes |
| Mse | 464 | 6 | Yes | Yes |
| Cca1 | 192 | 2 | Yes | Yes |
| **Total** | **676** | **14** | **All five** | **Every candidate exactly once** |

This independently establishes reduced launch count and correct routing for completed native small-panel runs. It makes no controlled timing claim.

### Configuration and provenance

I loaded the existing `local.toml` through the unchanged load path with parser-generated overrides. No new CLI flags gave `(intarna_batch_size=128, short_anti_repeat_rescue=False)`; explicit size 1 plus enable gave `(1, True)`; explicit disable gave `(128, False)`. The constructor rejects boolean/non-integer/nonpositive batch sizes, and the scope check rejects experimental model mode and Type V.

I called `command_provenance()` with the actual configured CasAndra command and configured child environment. It records:

- interpreter checkout: `784cb6dc5b38520cc84defcaa07e6461787fe6fc`;
- module source: `new_version/external/CasAndra/src/casandra/__init__.py`;
- module source SHA-256: `03f4a7769bfcaa957d9a8ac365ccfbec434e20b4aa89634d72342fda226ff847`;
- CasAndra checkout: `ba046f90242ac91dc075183ea51dcbedfac28484`;
- explicit resolution description: `static explicit effective PYTHONPATH; no import executed`.

The focused test verifies that package code is not imported and that an unresolved module does not retain the interpreter's revision as its own. This fixes the demonstrated historical attribution error for the configured execution. Static PYTHONPATH inspection is not a general Python import resolver: isolated/environment-ignoring interpreter flags, relative child paths, namespace packages, site packages and import hooks are outside this resolution method. Preserve the explicit limitation; do not claim universal proof of the imported module.

## Concrete finding and correction

**E1 — compatibility wrapper overwrote prior raw evidence; now resolved.** `interaction(candidate, work, settings, runner)` now delegates to a fresh one-candidate `interactions()` call. Every call begins at group 0/chunk 0 and writes `intarna-000000-000000.*`. Two calls with candidates `first` and `second` using the same working directory leave only `{"t000000": "second"}` in the retained ID map. I reproduced this with a header-only fake runner, with no native work.

This is not a defect in the main pipeline's current usage: standard interactions execute once per record directory, and experimental interactions execute once in a distinct `short-rescue` directory. It is nevertheless a regression in the preserved legacy entry point's ability to retain evidence from the former per-candidate calling pattern. The executor corrected the wrapper using `tempfile.mkdtemp(prefix="intarna-single-", dir=work)`. Candidate IDs are retained as mapping values and are not interpolated into filesystem paths. I inspected the correction and independently ran its four-call regression (distinct IDs, repeated ID, and `../unsafe`): **1 passed in 0.03 seconds**. A focused batch/parser selection additionally gave **27 passed, 5 deselected in 0.07 seconds**. The main `interactions()` function and native timing code did not change. The corrected `search.py` SHA-256 is `b79a58d279f008c0a561e2c7de37057d1c22e0ae444785b39a7b304020fc7db0`. This resolves E1 without invalidating the main-pipeline native results already inspected.

No other main-pipeline engineering blocker was found in the code review or final evidence audit. Existing missing-executable pre-launch provenance gaps identified in phase 1 remain deferred; the implementation did not alter `Runner` or claim to resolve them.

## Review chronology

The sections below preserve the native diagnosis, prebuild recommendation and build approval as they occurred. Their then-pending gates are resolved by the final evidence section. They do not supersede the final gate table above.

## Stock FASTA36 reproducibility failure: independent diagnosis

The first `natural_default` run is preserved. I independently compared all 79 ordered prediction lists with their frozen historical reports using the same four-field additive projection as above. **78 match exactly; Lsp2 differs, 14 historical versus 16 fresh candidates.** The fresh total is 3,397. This is an actual failure of the frozen exact-default gate, even though aggregate localization metrics remain unchanged.

For Lsp2, the historical and fresh FASTA36 commands are both `fasta36 -T 1 -m 8 repeats.fasta anti-search.fasta`, with executable SHA-256 `9a6e368d4d955a169d3749b5f752a8e87816d26fcfb4c7ec8e577adc36744d55`. I checked byte hashes independently:

| Native input | Identical SHA-256 in both runs |
| --- | --- |
| `repeats.fasta` | `aa2dc15532b143b2986e0cbd278491df0e534c1e73bd1163fff29db61d8fe692` |
| `anti-search.fasta` | `338cce665ad81f25cfa9ce6c9e530979fe8763f1e8ae1427674037793f8ae03f` |
| `input.fasta` | `f1244ed8c1b7e7f8befc07ad9378d666dd4cc4264327938c1d0b705138fbb343` |

The historical native table contains seven alignments; the fresh table contains eight. The extra row is `repeat000000r`, identity 88.46%, query coordinates 26→1, source coordinates 19880–19905, E-value 0.00059. Other native rows also have changed E-values/bit scores despite identical identity/coordinates. This precedes candidate generation and IntaRNA, so batching cannot cause the added native row. The executor's six controlled stock-binary repetitions retain eight alignments but have **six distinct full output hashes**; all repetitions, including this result, remain archived. The archived pre-change pipeline replay reports Lsp2=16 and reverse Ain=132, exactly matching the corresponding new-pipeline projections (`archived_baseline_replay/comparison.json`). That establishes that the behavior is not introduced by the new detector code.

### Verified upstream source and binary evidence

I downloaded the exact [upstream fasta-v36.3.8g source archive](https://github.com/wrpearson/fasta36/archive/refs/tags/fasta-v36.3.8g.tar.gz) into `/tmp/crispr-council-engineering/fasta-v36.3.8g.tar.gz`. Its SHA-256, `fa5318b6f8d6a3cfdef0d29de530eb005bfd3ca05835faa6ad63663f8dce7b2e`, exactly matches the installed Bioconda package recipe (`fasta3-36.3.8-hec16e2b_7/info/recipe/meta.yaml`). Selected unmodified source files were extracted only to the temporary audit directory; no dependency code was changed or built by this reviewer.

The installed manual documents shuffled significance estimation for pairwise/small-library comparisons. Exact source confirms small-sample shuffling in `comp_lib9.c:1066–1085`, the default DNA secondary-hit E cutoff of 0.001 in `initfa.c:1517`, and conversion of that E cutoff through estimated statistics to `repeat_thresh` in `comp_lib9.c:1525`. Therefore different random statistics can change which secondary alignments are emitted, even though the Python parser does not gate on those E-values.

There is also an upstream C initialization defect. `mrandom.c:42` defines `my_srand(int set)`; it uses a Park–Miller RNG with a positive supplied seed, otherwise a wall-clock-derived seed (`gettimeofday`). Production files declare the function without argument types and call `my_srand()` **without an argument** at `comp_lib9.c:570` and `work_thr2.c:148`. That is an argument/signature mismatch with undefined behavior; it cannot be treated as a supported reproducible seed contract. I inspected the installed executable's `my_srand` and its two call sites with `objdump`: the binary reads the argument register and the calls do not set an explicit seed. This corroborates the source-level diagnosis. It does not justify assigning one exact historical missing row to a known historical RNG seed, which was never recorded.

No documented/configured stock CLI random-seed option was found in installed help, manual or exact production source. `-z -1`, changing `-E`, or raising shuffle counts changes native statistics/secondary-hit behavior and would be a scientific search change, not a harmless reproducibility toggle. Repeatedly rerunning until historical counts happen to match would not repair the failed gate.

### Bounded repair recommendation before the frozen amendment

The coordinator is considering a **separately named source-built native dependency**, restricted to FASTA `-T1`, using fixed untuned master seed **20260927**. This is a new pinned dependency profile, not an alleged stock FASTA36 flag. Preserve the original executable, all old runs, and the original failed gate.

A minimal patch should provide an explicit shared `my_srand(int)` declaration, remove both production untyped declarations, and supply positive seeds at both manager and worker initialization. The manager uses master seed 20260927. Initialize the worker only after determining its index and use a deterministic distinct stream seed such as `20260927 + worker_index + 1`, within 1..2147483646. Preserve the existing RNG recurrence, scoring, E thresholds, statistics, and shuffle counts. `randtest.c` already supplies an explicit argument and is not part of production execution; audit every call site nonetheless.

**A fixed per-worker seed does not prove arbitrary multithread reproducibility.** Workers receive shuffle buffers dynamically. The bounded profile must explicitly pin one FASTA worker; later parallel guarantees require separate scheduling-independent RNG design or evidence. Reproducibility tests should vary process start times, output paths and launch order while preserving the exact ordered input bytes and declared thread count. Do not silently promise invariance to changed FASTA record order or thread count.

Build feasibility is supported by available GCC/make and the same upstream `Makefile.linux64_sse2` used by the Bioconda recipe. Build in a separate staging tree, only the required FASTA target, with documented compiler/flags and no fast-math retuning. Retain the archive hash, full patch, build log, compiler version, binary hash, upstream LICENSE/COPYRIGHT and modification notices. The archive includes Apache-2.0 licensing; retain its complete notices in the vendored dependency rather than stripping them.

Before full new-profile panels: exact repeated-output byte equality on at least six Lsp2 trials plus Ain/Kki/Mse with fixed inputs, one worker, and path/start-time perturbation; malformed-output checks still enforced. Then compare archived and current Python pipelines on all 79 loci using the **same new binary/settings**. Freeze the new protocol, paths, source and binary hashes before examining results. Exact paired-current/archived success would establish the new dependency profile's implementation comparison; it must not replace the failed stock historical gate.

## Completed controlled IntaRNA timing review

I independently loaded `performance/*/interactions.json` and every corresponding command ledger. On Mse (464 candidates), Ssu (104), and Sth1A (128), size 1 and size 128 maps are identical; after removing only paired-position additions, they also equal frozen historical interactions. Every call uses one native thread and executable SHA-256 `4605114ecd777c34d9f2b35f7b1e834f0c135afcb7fa5de25102b832147820b6`.

| Fixture | Serial/grouped calls | Serial/grouped stage wall seconds | Grouped fraction |
| --- | --- | --- | --- |
| Mse | 464 / 6 | 44.614 / 29.904 | 0.6703 |
| Ssu | 104 / 2 | 11.543 / 8.148 | 0.7059 |
| Sth1A | 128 / 2 | 15.045 / 11.204 | 0.7447 |

Command-time sums are consistent with stage wall times: Mse 43.068/29.889 s, Ssu 11.365/8.142 s, Sth1A 14.817/11.197 s; the remaining stage time includes input/mapping/parser/log work. The executor scheduled this controlled phase before the full natural panel; no additional native work was launched by this reviewer. These measurements pass the Mse speed target for the tested host/settings and are separate from the unrelated FASTA36 historical mismatch. The all-79 frozen-candidate replay subsequently completed; its final direct comparison is recorded below.

The coordinator subsequently froze `docs/COUNCIL_REPRODUCIBILITY_AMENDMENT.md`,
`protocol_amendment_v2.json`, and `amendment_prebuild_clarifications.md`, adopting
the separate dependency, explicit manager/worker seeds and one-worker scope.
Build/patch inspection was still pending at that checkpoint; its completed review
follows. The historical stock exact-equivalence gate remains **FAIL**.

## Seeded native build review before paired panels

**Decision: build/source contract approved for the declared Linux pthread, FASTA
`-T1` profile. Proceed to the paired 79-locus stages only after the amendment's
mandatory repeatability probes pass. No source/build engineering blocker was found.
The original stock historical-equivalence gate remains FAIL.**

I reviewed the frozen amendment, machine-readable amendment protocol and prebuild
clarification, then independently audited `native/fasta36-seeded-20260927` without
editing, rebuilding or launching scientific panels.

- Compared all **346 original regular files** in the verified archive with the
  retained source tree: only `src/comp_lib9.c` and `src/work_thr2.c` differ; no
  original files are missing. `src/council_seed.h` is the declared new header.
  All **350** recorded source-file/alias hashes in `build_manifest.json` verify.
- The manager now passes master seed **20260927** through a typed declaration.
  The worker index is assigned before RNG initialization; a range check keeps
  `20260927 + worker_index + 1` within 1..2147483646. Both callers include the
  shared typed header. Unsupported `PCOMPLIB`/MPI builds fail explicitly.
- The RNG implementation `mrandom.c`, recurrence, shuffle routines/counts,
  alignment/scoring code, statistical estimators and thresholds are byte-identical
  to the upstream archive. The installed stock binary remains unchanged.
- Independently inspected the new executable's disassembly: the manager supplies
  `0x135283f` (20260927); the worker passes index plus `0x1352840` (20260928) after
  the bounds check. This verifies the seeds actually reached the compiled calls.
- Verified archive, patch, new binary, builder script, both build logs and explicit
  link-library hashes against the manifest. The accepted binary SHA-256 is
  **`e53ea2de4b1584388da0841222b6254fab8f1ebd9e5b0863cd42f341132f348e`**;
  patch SHA-256 is **`c273f64943434b414b03041e58c39887a5599ff257587f1a4edae6d41dc9f50a`**.
- GCC 13.3.0 / GNU Make 4.3 and the exact compile/link commands are retained.
  The selected upstream `Makefile.linux_sse2` resolves to bytes identical to
  `Makefile.linux64_sse2`; compile options include `-g -O -msse2` and the upstream
  statistics/thread defines, with no fast-math option. This is honestly a new
  source-built dependency, not a binary-identical recreation of the installed
  Bioconda build.
- Attempt 1 compiled the objects and failed only while resolving `-lz`; its exit
  code 2, complete log and initial manifest remain retained. Attempt 2 successfully
  linked the same objects after an explicit versioned system-zlib library argument.
  The second log contains the final link command and success. No source, seed,
  statistics or optimization adjustment occurred between attempts. `ldd` shows
  libc/libm/the loader as runtime dependencies; the manifest correctly describes
  the zlib pathname as a **link** library, without claiming that it remains an ELF
  runtime dependency.
- Upstream LICENSE/COPYRIGHT copies exactly match the retained source tree, and the
  local modification notice explains the seed changes, separate binary and
  one-thread limitation. Source archive, patch, source files and build evidence
  are retained for review.
- `local.reproducible.toml` resolves to this separate executable with `threads=1`.
  I compared `local.toml` with the frozen archive and verified identical bytes;
  stock FASTA36 still hashes to `9a6e368d4d955a169d3749b5f752a8e87816d26fcfb4c7ec8e577adc36744d55`.

The patch does not impose a universal worker-scheduling guarantee, and the notice,
manifest and selected configuration correctly limit the claim to one native worker.
Runtime repeatability is an empirical gate still to be established by the executor's
predeclared probes; source inspection alone does not mark it passed. Paired archived
and current runs must use this exact binary/hash and identical inputs/settings.


## Final independent evidence and verdict

**ACCEPT.** The final state closes the bounded engineering gates. No threshold was
retuned after results, the seed is the predeclared 20260927, and the stock-native
historical failure remains separate from the successful amended comparison. I
made no implementation/build changes or native benchmark reruns during this final
review. Direct artifact-check results are retained in
`benchmarks/council_round2/engineering_final_independent_checks.json`.

### Direct final artifact checks

I loaded the actual report files for all 79 paired loci, flattened each report's
ordered predictions, and compared archived/current default after removing only
the four previously declared additive evidence fields from the current output.
**All 79 ordered lists and all 3,395 existing predictions match exactly.** This
includes IDs, export IDs, coordinates, sequence, native interaction values,
ranking and old evidence. Separately, experimental output retains every default
prediction in the same relative order with identical evidence; only `export_id`
is excluded in that comparison because accepted additions shift displayed ranks.
There are exactly **11 added short candidates**, with no standard loss.

I independently read every `fixed_candidate_batch_replay/*/interactions.json`,
verified the cited frozen report hash, and compared the complete interaction map
against its frozen anti-repeat candidates, removing only the two paired-position
additions. **79 loci / 3,390 interactions match, including no-hits**. Actual
command-ledger lengths sum to **152 calls**, consistent with the declared replay.
This new native execution complements the earlier all-3,390 historical CSV parser
check and the independent 676-candidate exact FASTA/ID-map routing audit; it does
not rely on rediscovering the same candidates with variable stock FASTA36.

I independently rehashed all native repeatability inputs and output tables and
checked each command's successful return, executable hash and `-T 1` argument.
All **six Lsp2** outputs are byte-identical; **two each of Ain, Kki and Mse** are
also identical within locus. Paths and interleaved launch positions vary as
predeclared. The accepted binary remains
`e53ea2de4b1584388da0841222b6254fab8f1ebd9e5b0863cd42f341132f348e`.
This passes the empirical one-worker repeatability prerequisite that was pending
at build review; it does not imply invariance to FASTA record reordering,
multithread scheduling, hardware or compiler changes.

I independently counted complete control artifacts: genomic default and
experimental each have **19 complete reports / 60 predictions / zero short
admissions**. Each conditional stage has **79 complete result files**. Their
proposal/admission totals are **126/2**, **79/0** and **100/1**. Final metrics retain
all 68 historical localized references, record 75/78 experimental localization,
zero natural losses and all frozen burden/length/runtime caps passing. The
experimental/default natural stage wall ratio is **1.0092** (347.170/344.020 s).
These are engineering consistency checks on the reports; the detection and
validation council provide the independent scientific/control interpretation.
The controlled batching result remains the separate, reviewed Mse **33.0% stage
wall-time reduction**, with exact maps also verified on Ssu and Sth1A.

### Source and verification bindings

I independently verified the four retained source archive hashes, every source
hash listed in the eight stage manifests, and their manifest hashes. The archived
stage's 12 module files bind to `baseline/source.zip`; the archived harness also
binds to its pre-run archive. All **seven current natural/control stages use the
same 13 core/harness files**, each identical to current live source and both
retained current source archives. The final reviewed wrapper and provenance fixes
are therefore included in the evaluated code, with no unreviewed core delta.

The immutable `implementation_final.zip` SHA-256 is
`fb0b61ded6486aadc7dc80e4917e7f1e29622c1435ce8be49254a89b3cb0148f`.
I verified all 32 members against live files. Thirty still match exactly; only
`README.md` and `docs/IMPROVEMENTS.md` changed during the final prose cleanup.
`documentation_final_sha256.json` correctly records both archived and final live
hashes, which I independently checked. The changes update obsolete future-batching
wording, final links and the concrete conditional-control localization caveat.
They do not change tested core, harness, tests, configuration or native behavior.
The immutable tested archive was retained, so no repeated scientific/test run is
needed for these prose-only deltas.

I rehashed all **19** final result/audit/test-log files cited in
`executor_final_validation.json`. The actual retained logs show **119 passed in
18.84 seconds**, with three native literature cases plus four native tool cases,
no skip reported; Ruff passes; strict mypy reports no issues in **12 source
files**. This is verification of executor-run evidence, not a claim that I reran
those native tests. My earlier 72-test focused execution, wrapper regression and
27-test batch/parser selection remain independent executed evidence.

Validation/coordinator preservation audits report **14,902 protected files** and
**390 frozen dependency/profile files** unchanged, with the original native
executables and protected repositories preserved. I rely on their independent
full traversal here; my final rehash scope is the source/build/result evidence
above. No protection blocker was reported by either reviewer.

### Remaining limits, without unresolved blockers

- The original stock-native historical exact gate stays **FAIL**. Stock FASTA36
  retains its pre-existing seed/signature problem; reproducible results require
  explicitly selecting the separate pinned profile. One-worker empirical success
  is not a universal platform or parallel guarantee.
- Type II short rescue remains **default-off**. Development/regression localization
  and the prescribed controls do not estimate held-out sensitivity, precision or
  a biological false-positive rate. Predicted windows are not transcript ends.
  In particular, the deranged-repeat Cga window overlaps 79 reference nucleotides
  while its seed is outside the reference and CM overlap is only two nucleotides;
  the Eit shuffle admission revisits an unvalidated natural-candidate region.
  Passing bounded control caps does not erase these scientific caveats.
- The configured CasAndra provenance correction uses static explicit absolute
  PYTHONPATH resolution. Other import modes remain outside that proof. The older
  missing-executable pre-launch ledger gap is deferred, unchanged and not claimed
  as fixed. Neither issue blocks the tested configured execution.
- The batching speed observation applies to the controlled host, native binary,
  fixtures and one-thread settings. No whole-genome memory/scaling or broad
  machine-independent speed claim was measured. Duplicate consolidation and
  wider architecture changes remain outside this round.

No additional engineering change is required before accepting this bounded
implementation and documented dependency profile. Scientific promotion decisions
remain with the same council and require the limitations above to stay visible.
