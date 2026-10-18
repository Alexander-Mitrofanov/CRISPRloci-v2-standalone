# Engineering council review — phase 1

Date: 2026-09-27. Role: software engineering. Scope: `src/crisprtracrrna`, tests, configuration/CLI, `docs/IMPROVEMENTS.md`, and retained benchmark evidence. No core implementation or protected repository changes were made. The prior 78-test baseline and 79-input benchmark are the starting point, not newly claimed test results.

## Recommendation

Implement **bounded IntaRNA batches grouped by the exact oriented repeat sequence**, with ID-aware parsing and validation of every returned row. This is a complementary engineering change that preserves the scientific inputs, all candidate hypotheses, and output evidence. It can be combined with a separately reviewed detection change without conflating an algorithmic sensitivity improvement with a performance improvement.

Do not implement general scheduling, persistent caches, hypothesis consolidation, or a process-runner rewrite in this cycle. The existing process-group timeout/cancellation behavior should remain in force for every batch.

## Demonstrated findings

### 1. IntaRNA is the largest measured tool-time component

`pipeline.py:129–130` invokes `search.interaction()` once per candidate; `search.py:347–368` starts a new process for each target/repeat pair. Both the target and query are single sequences, so the configured IntaRNA thread count cannot parallelize multiple target/query combinations within these calls.

I inspected all 79 retained `benchmarks/final_fixed/complete/*/work/*/logs/commands.json` files:

| Component | Processes | Sum of logged seconds |
| --- | ---: | ---: |
| IntaRNA | 3,390 | 408.951 |
| CasAndra | 79 | 194.918 |
| FASTA36 | 79 | 110.949 |
| CRISPRidentify | 79 | 92.698 |
| Infernal | 151 | 21.188 |
| ERPIN | 72 | 1.846 |

IntaRNA accounts for about 49% of the 830.55 seconds summed across these commands. These are retained per-command timings, not a fresh end-to-end wall-time trial.

Only 42 calls reuse an identical `(candidate_sequence, oriented_repeat)` pair within a record. A pure exact-pair cache saves little here. Grouping candidates by their exact oriented repeat gives 152 batches overall; a maximum of 64 targets per process gives 166 batches, and a maximum of 128 gives 152. No biological hypotheses need to be removed.

The installed IntaRNA 3.2.0 `--fullhelp` confirms that `--target` and `--query` accept FASTA files; `--outCsvCols` supports `id1` and `id2`; `--threads` parallelizes target/query combinations; and the default `--outNumber` is one interaction per pair. A multi-target file plus **one** query sequence therefore evaluates exactly the intended pairs. A multi-target/multi-query Cartesian product would not.

### 2. A native batching probe preserved every Mse interaction

A temporary, non-pipeline probe reconstructed the 464 Mse target/repeat pairs from retained commands, wrote six target FASTA files and six single-record query FASTA files, and ran the installed binary with otherwise matching scientific settings and two threads. The retained Mse calls used one thread. Candidate IDs were FASTA IDs; the requested columns added `id1,id2` to the existing columns. Each returned row was routed by ID and normalized through the existing single-pair parser.

All 464 interaction dictionaries matched the retained serial result exactly, including candidates with no interaction. The six groups contained 52, 52, 90, 90, 90, and 90 targets; they returned 33, 39, 89, 90, 90, and 90 interaction rows respectively. Thus missing rows are an actual supported no-interaction result and must map to `None` for the expected candidate.

The fresh grouped run took 15.65 seconds; retained serial commands total 46.64 seconds. This establishes feasibility and a promising timing direction, **not** a controlled threefold-speedup claim: the old and new timings were collected at different times and used different thread counts. The probe and outputs were confined to `/tmp/crispr-council-engineering`.

### 3. Current IntaRNA validation can hide malformed scientific output

`search.py:329–340` reads all rows, selects the minimum energy, and only then validates the selected row's coordinates. A valid row at energy `-10` followed by a row at energy `-9` with target coordinates `999..1000` is accepted for a 20-base target. This was reproduced during review. Batched parsing must not perpetuate this behavior: validate every row before selecting evidence.

The selection key currently contains only energy. Equal-energy valid rows select whichever arrives first. Introduce a documented total tie key, such as `(energy, target_start, target_end, repeat_start, repeat_end, hybrid_dot_bracket)`, in the common row-selection routine. Then reordered CSV rows cannot change the selected coordinates and downstream consistency features. This is particularly relevant when batch execution may reorder output.

### 4. Runner provenance has a real but separate pre-launch gap

`process.py:41–47` resolves and hashes the executable before adding a record or entering the persistence `finally`. A nonexistent binary raises the appropriate `ToolError`, but leaves no attempted-command entry in `commands.json`; this was reproduced. The enclosing pipeline still writes `failure.json`, so this does not turn a failed tool into a scientific negative. Treat comprehensive pre-launch provenance and atomic command-ledger writes as a later bounded task; they are not needed to achieve batching safely.

### 5. Existing scope limits remain justified

The current CLI/config and immutable coordinate types provide a small, clear interface. Preserve the short FASTA36 filenames, array masking, independent genomic model rescue, opaque IDs, source reconstruction, and explicit failure behavior. Changing strand hypotheses or deduplicating overlapping predictions as a side effect of this engineering work would be an unreviewed scientific change. `docs/IMPROVEMENTS.md` is correct to defer hypothesis consolidation, circular genomes, persistent caching, and resume until their contracts are defined.

## Proposed implementation contract

### Files and interfaces

1. **`src/crisprtracrrna/search.py`**: add `interactions(candidates, work, settings, runner) -> dict[str, dict[str, Any] | None]`. Return an entry for every supplied candidate, with `None` for a valid absent hit. An empty candidate list launches no process. Reject duplicate candidate IDs before writing inputs. Keep the existing `interaction()` and `parse_intarna()` interfaces as thin compatibility entry points or preserve their behavior with a shared validated-row routine.
2. Group by exact `candidate.anti.repeat`, including orientation; keep every candidate even if sequences are identical. Sort groups and targets deterministically. Write opaque group/query IDs and candidate IDs into FASTA headers. Each process takes one single-record query file and at most the configured number of target records. Add `id1,id2` to CSV output and explicitly request `--outNumber 1`.
3. Add a parser receiving expected target IDs/lengths and the expected query ID/length. Require all named columns; validate every row's IDs, complete fields, finite energy, ordered positive coordinates within both sequences, and the strict hybrid structure contract below. Unknown IDs and malformed rows are `ContractError`, including malformed nonwinning rows. Header-only output is a valid set of absent hits. Preserve support for harmless extra columns. If the common parser accepts multiple valid interactions, select deterministically using the total tie key rather than file order; document that normal pipeline calls request one interaction.
4. **`src/crisprtracrrna/pipeline.py`**: call `interactions()` once before the evidence loop, then look up each candidate by ID. Keep scoring and result ordering otherwise unchanged. Do not catch `ToolError`/`ContractError` and convert them to an empty evidence map.
5. **`src/crisprtracrrna/config.py`, `cli.py`, `config.example.toml`**: add a positive, non-boolean integer `intarna_batch_size`, default 128, exposed as `--intarna-batch-size`. Size 1 provides a serial comparison/debugging path using the same ID-aware contract. Continue using the existing `threads` and per-process `timeout`; avoid adding worker pools or multiplying the CPU budget. Document that the timeout applies per batch and memory grows with active IntaRNA threads.
6. **`tests/test_scientific_parsers.py` plus a focused batch test module**, and opt-in **`tests/test_real_tools.py`**: cover the contracts below. Update **`docs/IMPROVEMENTS.md`, `docs/MIGRATION.md`, and `README.md`** with actual behavior and measured acceptance results. Do not advertise a speedup from the exploratory comparison alone.

### Compatibility and migration

- Preserve candidate IDs, all hypotheses, `interaction` dictionary field names/semantics, `None` for no interaction, coordinates, scores, predictions, and all export columns. No report-schema major bump is needed for this additive configuration/performance change.
- Existing configurations continue to work with the default batch size; size 1 is the comparison path. Invalid sizes fail as user-readable configuration errors.
- Raw evidence files will become per-group/per-chunk CSV and FASTA rather than one CSV per candidate. Document this internal artifact-layout change; raw rows carry stable candidate IDs and remain covered by the manifest. Do not reconstruct misleading per-candidate process records.
- Additive paired-position fields can support the opt-in detection branch below; preserve every pre-existing interaction field. The new deterministic equal-energy tie rule can select a different coordinate interval only where the previous parser depended on row order. Document this limited intentional correction and test it directly.
- Unsupported native CLI behavior or malformed output must fail explicitly. Avoid implicit serial retry after scientific-tool failure, which could conceal a broken batch contract.

## Integration with the council's opt-in short-anti-repeat branch

The detection council proposes an experimental, default-off Type II branch requiring actual paired bases inside `candidate.anti_local`, plus independently specified structural corroboration. Engineering should expose validated pairing positions without choosing its biological thresholds.

The official [IntaRNA CSV documentation](https://github.com/BackofenLab/IntaRNA#customizable-csv-rna-rna-interaction-output) specifies that `hybridDP` contains an interaction-site string for the target, then `&`, then one for the query. Target paired positions use `(`, query paired positions use `)`, and dots are unpaired bases. Coordinates identify the interaction within each full sequence. This is a source-format contract, not evidence that any particular minimum paired-base count establishes a true tracrRNA.

Implement a shared structure validator/mapping helper with these checks:

1. Exactly one `&`; target side contains only `.` and `(`; repeat side only `.` and `)`.
2. Target-side length equals `end1 - start1 + 1`; repeat-side length equals `end2 - start2 + 1`.
3. A positive number of paired positions, with equal counts on both sides. Reject zero-pair and unequal-pair encodings.
4. The already validated interval bounds place every derived position inside its full target or query. Derive `target_paired_positions = [start1 + i for each '(' at index i]` and `repeat_paired_positions = [start2 + j for each ')' at index j]`.

These additive arrays use **one-based positions in the full oriented candidate sequence and full oriented tested repeat**, respectively, not positions relative to the interaction site or input forward genome. Existing interval, structure, and energy fields retain their values. The detection predicate counts `anti_local.start <= p - 1 < anti_local.end` for each target paired position `p`; it must not substitute interval overlap, which includes unpaired bulges. For a reverse-strand candidate, these remain transcript-oriented coordinates; use `Candidate.to_genome()` if genomic presentation is needed. No sequence-complementarity enforcement is required here; the native energy model already admits its supported pair types.

I inspected all **2,933 nonempty interaction rows** retained in the 79-locus baseline: every row has these structure characters, matching side lengths, and equal pair counts. This therefore tightens corrupt-output detection without rejecting these real outputs. Some synthetic existing parser-test strings have inconsistent lengths and should be replaced with meaningful valid structures when enabling validation. Add bulged examples in which the interaction span overlaps the anti-repeat but its paired positions do not; add boundary tests for the first and last anti-repeat bases and reverse-oriented candidates.

The supported/tested CLI contract is IntaRNA **3.2.0**, the installed version. Request `id1,id2,start1,end1,start2,end2,E,hybridDP` and `--outNumber 1`; retain existing energy/accessibility/seed defaults. Size 1 is an explicit serial comparison path using the same parser, not an automatic fallback. Do not invent broader version support without tests or change parameters to accommodate malformed output. Stable future compatible output can be accepted by the contract; a hard-coded exact-version rejection is unnecessary.

The parent has frozen baseline material under `benchmarks/council_round2/baseline/`; implementation validation should consume that snapshot and write separate round-two output directories. Original measurements above came from the corresponding retained `final_fixed` outputs before the freeze.

## Meaningful acceptance tests

1. A fake runner receives candidates sharing one repeat and candidates with a different or reverse-complement repeat. Assert the exact target/query pairs and maximum chunk size, not just the number of calls. At `batch_size=2`, five targets in a repeat group produce three chunks and no cross-group query comparisons.
2. Parse output deliberately reordered across candidates; equal-energy valid rows reordered within a candidate yield identical evidence. Include absent-hit candidates, header-only output, and multiple valid rows.
3. Fail on unknown target/query IDs, duplicate candidate IDs in input, missing/duplicate required columns, truncated required fields, non-finite energy, out-of-range coordinates, invalid coordinate ordering, empty or inconsistent structure, and a malformed nonbest row. Extra named columns remain accepted.
4. Inject a batch tool failure or malformed result and verify the pipeline emits `failure.json` and never publishes a complete `report.json`. Existing process-group timeout tests remain applicable because every batch still uses `Runner`.
5. A real-tool comparison should run the exact same frozen candidate inputs in batch sizes 1 and 128, compare the complete interaction maps including no-hits, and compare downstream rankings/exports after normalizing execution metadata only. Mse is the high-count fixture; add a multiple-repeat fixture and Ssu for repeated target sequences. The Mse exploratory probe is not a substitute for this committed test/validation record.
6. All existing unit/contract tests and relevant real-tool tests pass. If accompanied by a scientific change, compare batching on identical generated candidates separately so biological differences cannot obscure an interaction-contract regression. Preserve complete execution on all 79 loci; the existing 68/78 localization is the minimum regression floor unless the council explicitly accepts a documented scientific tradeoff.

## Performance acceptance

Use the same installed executable, thread count, input candidates, and host for fresh size-1 and size-128 runs. Report subprocess count, summed process wall time, total interaction-stage wall time, and peak memory if available. Run serial and grouped trials in alternating order if more than one trial is needed.

The deterministic scaling target is `sum(ceil(group_size / batch_size))` invocations, with each group defined per record and exact oriented repeat. The retained baseline gives **152 instead of 3,390** at size 128. Fresh Mse grouped wall time must be lower than serial on a controlled comparison with identical evidence; a suggested target is at least 20% lower interaction-stage time, subject to host variance. Do not promise an overall pipeline speedup when upstream runtimes dominate a particular input. No extra processes should run concurrently beyond the existing sequential pipeline, and IntaRNA keeps the configured thread ceiling.

## Evaluation commitment

The same engineering council member should independently inspect the implemented diff, run the malformed-output and pair-routing tests, and verify the controlled equivalence/performance record. Final acceptance requires all three: unchanged intended pairs and evidence, explicit failure on corrupted outputs, and a demonstrated reduction in launch count. Timing is secondary to correctness.
