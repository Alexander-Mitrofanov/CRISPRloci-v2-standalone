# Round 2 amendment: an explicit native dependency comparison

Frozen 2026-09-27 before seeded native results or experimental-profile results.
This adds a reproducibility repair to the original council plan. It does not
replace the original protocol or erase its failed historical-equivalence gate.

## Observed failure and attribution

The first fresh default run completed all 79 inputs with unchanged reference
localization and top-k totals, but produced 3,397 rather than 3,395 candidates.
Only Lsp2 differed: 16 rather than 14 candidates. The first divergence is in
FASTA36 output, before IntaRNA. Masked input, repeat FASTA, arguments and executable
hash are identical to the historical run. All reconstruction/export invariants pass.

Six prospectively bounded native repetitions produced six different statistical
output hashes, but all six had the same eight alignment keys. They demonstrate
statistical variability; they did not reproduce the historical seven-row output.
The archived detector independently produced 16 Lsp2 candidates and 132 candidates
on reversed Ain, exactly matching the new detector's biological projections on
those inputs. The latter exposes upstream array-inference asymmetry, not a coordinate
conversion error. Retain these results and the original exact-comparison FAIL.

The installed FASTA36 36.3.8g source has `my_srand(int set)` but two production
call sites/declarations omit the argument. Its normal unseeded initialization uses
the clock. Source and installed-code inspection found no supported CLI seed option.
Do not invent one, disable statistics, adjust E-value thresholds, or select a seed
that reproduces a desired reference result.

## Bounded dependency repair

The same detection, engineering and validation council members support an explicit
new dependency-controlled comparison. The implementation agent builds a separately
named executable from the exact verified upstream archive, under the writable
project. Keep the installed executable and protected upstream repositories unchanged.

The minimal source patch is fixed before testing:

1. Add a header with a correctly typed `void *my_srand(int set)` declaration and
   master seed **20260927**, chosen from the already frozen control seed.
2. Replace the old no-argument declarations in `comp_lib9.c` and `work_thr2.c`.
   Pass the master seed explicitly at the manager call.
3. Initialize each worker after obtaining its worker index, using master seed plus
   worker index plus one; all seeds must remain in 1..2147483646. Handle the declared
   build path correctly; do not promise unsupported MPI or platform behavior.
4. Preserve the RNG recurrence, shuffling algorithm/counts, alignment/scoring code,
   statistical estimation and all search/reporting thresholds. No fast-math flags,
   upstream version upgrade, model changes or biological parameter changes.

Build only FASTA36 with the verified upstream Linux SSE2 makefile, using available
gcc/make. Preserve the original archive and checksum, exact patch, modified source,
compiler/version, complete build command/flags/log, binary checksum, copyright/license
and modification notice. Freeze these build artifacts and an execution manifest
before any amended full-panel result is inspected. Build failures may be corrected
as engineering failures with retained attempts; no data-driven seed selection.

This is a separately named, explicitly selected native dependency profile. Pin
FASTA36 to **one thread** for all amended comparisons and the recommended local
reproducibility configuration. Do not claim arbitrary multithread reproducibility:
native workers receive dynamically scheduled shuffle buffers. Keep existing local.toml
byte-identical; add a clearly documented separate configuration for the pinned binary.
Stock-FASTA users retain an explicitly documented reproducibility limitation.

## Amended execution and acceptance

All scientific gates in the original protocol remain fixed: 78 mapped references,
historical 68-reference retention, top1 >=53, top5 >=63, coverage90 >=64, useful
experimental localization >=69, and every original absolute burden/runtime/control
limit. No default promotion of the short-rescue experiment is planned.

1. Prove repeatability using at least six exact-byte Lsp2 repetitions. Test repeated
   source/repeat inputs for Ain, Kki and Mse with path/order perturbations, respecting
   the declared single-thread scope. Fail explicitly if exact repeatability fails.
2. Run the archived pre-change detector and current default detector on all 79 inputs
   with the identical new binary, seed, inputs, one thread per locus and two loci in
   parallel. Use new output directories. The archived Python detector may select
   this executable through its existing `fasta36` setting; no seed-related Python
   source patch is needed. If extraction/resource relocation is needed, record it
   and verify original source hashes and resource bytes.
3. Require exact candidate, interaction, evidence and ranking projections between
   those two dependency-controlled stages, allowing only documented additive fields
   and run metadata. Record historical-to-fixed archived gains/losses, burden and
   rank changes separately. Never use a new baseline to excuse losing any of the
   historically localized 68 references or weakening original floors/caps.
4. Separately finish the fixed-input native IntaRNA replay over all 79 historical
   candidate sets (3,390 interactions). This is a component proof, not genomic
   discovery. Preserve the controlled size1/128 timing results already obtained.
5. Only after the paired default gate passes, run the full experimental 79-input
   panel with the same pinned dependency, followed by the original 38 paired genomic
   controls and all three conditional repeat panels. Retain the frozen control bytes,
   donor pool and seeds. Audit every run and preserve all initial failures.
6. Attribute a gain to short rescue only when it is absent from the fixed-dependency
   default and present in the otherwise identical experimental run. Dependency-driven
   changes and experimental-branch changes must be reported separately.
7. Repeat final tests with the selected reproducibility configuration. The same
   council evaluates source/build contracts, exact paired results, historical deltas,
   controls, preservation and scientific limits. Each review shows original historical
   exact equivalence **FAIL** separately from any amended comparison **PASS**.

The additional full baseline run is authorized work toward repairing the reproducibility
defect. This remains a development/regression study, with no independent precision,
FDR, transcript-end accuracy or generalization claim.
