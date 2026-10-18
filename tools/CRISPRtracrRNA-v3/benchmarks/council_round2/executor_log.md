# Executor log

## Pre-result implementation notes

- Implemented exact-repeat grouping, configured default 128, strict ID/structure validation for every native interaction row, deterministic energy/coordinate ties, and retained target-ID maps.
- Short proposals reuse the native genomic alignment file. Standard candidate-window model/terminator calls remain separate from experimental proposal calls, preserving standard CM search context. Genomic CM rescue executes before experimental insertion, retaining standard candidates and IDs.
- IntaRNA 3.2.0 documents `seedTRange` as single-target only. Anchored experimental calls are therefore individually executed; this is allowed by the frozen plan.
- Coordinator identified an existing provenance defect before new panel results: `python -m casandra` inherited the interpreter checkout revision. Corrected metadata records interpreter-checkout revision separately and statically resolves module source/revision from effective explicit PYTHONPATH without importing module code. Unresolved module paths remain explicit. This changes provenance only.
- Initial test run: 104 passed (including all seven integration cases), two new fake-runner tests failed because their fixture treated FASTA iterators as lists. Fixed fixture; no scientific code or threshold adjustment resulted.

## Independent preliminary review corrections

- Detection council found a conditional-harness annotation mismatch: evidence class used any raw forward model/terminator instead of only qualifying short-admission support. Corrected before conditional runs. Admission rules and counts were unaffected.
- Engineering council found that repeated calls to the legacy single-candidate compatibility wrapper reused batch filenames. The production batched pipeline does not use this wrapper. Added a unique private subdirectory per legacy call and a regression test for repeated and unsafe-looking candidate IDs; no changes to batched scientific inputs or controlled timing code.

- Validation council corrected two audit defects before full-panel audit: the default comparator must never remove unexpected short candidates, and provenance needs explicit cross-links to frozen inputs/settings/protocol rather than only internally consistent artifact hashes. Both are enforced. Conditional audits also bind every array's control repeat to its frozen panel, count accepted proposals directly and require all 79 completions. Genomic-control metrics retain null biological outcomes because no truth annotation exists.

## Historical equivalence failure and bounded diagnosis

The first fresh default panel completed 79/79 in 346.42 seconds and retained 68/78
localization, top 1/3/5/10 = 53/61/63/64 and 64 at 90% coverage. It has 3,397 candidates,
versus 3,395 historical. Exact comparison FAILED at Lsp2 (16 versus14), and the
execution sequence stopped before experimental runs. The initial panel is retained.
All 79 reconstruction/export/provenance invariants pass separately.

Lsp2's masked genomic input, repeat-query FASTA and native argv/executable are
identical. First divergence occurs in native FASTA36, which reports an extra
reverse-repeat counterpart at 19880–19905; changing candidate composition then
changes ERPIN background scores. All other 78 prediction projections are exact.
Six predeclared unchanged native repetitions each reported 8 alignments but all
six statistical output hashes differ. The archived baseline detector independently
produced 16 Lsp2 candidates and 132 reverse-Ain candidates, exactly matching the new
detector on both inputs. Thus both the historical stochastic-search difference
and reverse-input asymmetry predate batching.

Engineering independently found no supported seed CLI in exact FASTA36 36.3.8g
source. No undocumented flag, statistics-disabling option, tuned seed or
rerun-until-match was used. Root and the same council are deciding the append-only
protocol response. A native IntaRNA-only replay on all 79 frozen candidate sets is
being retained as a distinct component-equivalence test, not genomic discovery.

## Frozen reproducibility amendment

Root froze `protocol_amendment_v2.json` before seeded native results. The separate
native source build changes only two RNG call/declaration sites plus a typed
header, using master seed 20260927 and worker seed master + index + 1. Both root
and engineering independently reviewed the exact patch and hashes. The installed
binary, protected repositories and `local.toml` remain unchanged. An initial link
failure due to the missing unversioned zlib development symlink is retained;
linking the installed versioned system library succeeded without scientific changes.

Six Lsp2 outputs and paired Ain/Kki/Mse outputs are byte-identical across different
paths and process launch orders. One native thread is required; arbitrary
multithread determinism is not claimed. The execution manifest was frozen before
full amended panels, and the archived Python detector's 40 original files were
verified byte-for-byte after extraction. All 115 tests, including seven real-tool
cases, pass with the selected reproducibility configuration (18.42 seconds).

The archived 79-input panel, current 79-input default panel, strict paired audit,
and then experimental 79-input panel are queued in that order. A failed gate stops
the sequence; no existing panel is overwritten.

Validation reviewed the amended harness before affected controls and found three
wiring issues: conditional execution ignored the requested config; pinned-stage
audits allowed missing profile metadata; and the final summary needed independent
retention of the historical 68 references rather than only the new baseline.
All are corrected with focused tests, and historical-to-dependency differences
are reported separately from experimental gains. No affected conditional run had
started, and no scientific cutoff or frozen native artifact changed.

## Amended natural panels complete

Archived and current defaults completed all 79 inputs with 3,395 hypotheses and
68/78 localized references. The paired biological projection is exactly equal,
including order and existing fields; source/export/evidence audits pass. Archived
and current wall times were 390.53 and 344.02 seconds respectively (one paired
observation). The original stock-native historical comparison remains failed.

The experimental panel completed all 79 inputs: 75/78 localized, top 1/3/5/10 =
60/68/70/71, 68 references with >=90% coverage. It accepted 11 of 66 short
proposals, producing 3,406 total hypotheses (p95 193, max 464); wall time was
347.17 seconds. All standard biological evidence is preserved; inserting new
hypotheses changes 345 rank-derived export IDs as explicitly permitted. Gains are
Cca1, Cga, Cme1, Ffr, Ghy2, Orh and Psp. Esp2, Bbo and Cca2 remain missed; Ghc2
remains unmapped. No historical detection was lost. Paired genomic and conditional
controls are the next required evidence, and the experiment remains disabled.

## Final controls and executor verification

Both 19-input genomic profiles completed and passed strict audits. Each produced
60 hypotheses; 30 short proposals produced zero admissions. Conditional panels
all completed 79 inputs with strict source/array/repeat/admission/native-evidence
audits: dinucleotide 1 admitted 2/126 proposals, dinucleotide 2 admitted 0/79, and
derangement admitted 1/100. Per-input medians are zero; maxima are 1, 0 and 1.
No repeat control is unchanged or identical across its other frozen panels.

`final_metrics.json` reports every amended quantitative gate as passing while
retaining the original historical exact-projection failure. Companion files
include all natural proposals, paired natural/genomic results, conditional
admission evidence, historical dependency deltas and oriented endpoint geometry.
`final_stage_source_bindings.json` verifies exact retained source archives for
every amended stage; all current stages share the same 13 core/harness files.

Final commands: `PYTHONPATH=src CRISPRTRACRRNA_TEST_CONFIG=$PWD/local.reproducible.toml
python -m pytest` passed all 119 tests, including all seven native cases, in 18.84
seconds. `python -m ruff check src tests scripts/council_*.py` and
`python -m mypy src` passed. No code changed after these final checks. The
experiment remains disabled pending independent generalization evidence; the
council's final read-only evaluations and coordinator preservation audit are
separate from this executor's claims.
