# Bug-fix log

Changes delivered in **v3.0.0**, grouped by their effect on users. Work was performed
on 2026-09-26–27. These phases are development milestones, not separate releases.

## 1. Correct sequence processing and evidence

The v2 rewrite replaced shared temporary files and loosely coupled parsing with
typed intervals, isolated records and explicit tool adapters.

| Problem | Corrected behavior |
| --- | --- |
| Multiple contigs could be concatenated and IDs truncated at a dot | Process each record independently; retain original IDs, headers and source hashes |
| Inclusive endpoints, negative starts and reverse projection were inconsistent | Use one half-open interval model internally, clip windows and explicitly convert public coordinates |
| Array exclusion checked only the hit start | Exclude every anti-repeat interval overlapping an accepted array |
| Gapped alignment length inflated coverage; exact cutoffs were rejected | Use repeat-coordinate span and inclusive identity/coverage limits |
| Model evidence could attach to the wrong candidate | Key hits by candidate/model and select supporting evidence explicitly |
| ERPIN read only the first hit and confused score with E-value | Parse all hits, retain both values and rank stronger significance correctly |
| Missing or reverse-oriented evidence could receive support bonuses | Score only available evidence supporting the reported transcript hypothesis |
| String `False` could enable a boolean option | Parse explicit CLI/TOML booleans |
| Failed computations, no-hit runs and stale output could be confused | Fail explicitly, refuse existing output paths and emit valid empty exports |

The rewrite also removed executable text parsing, corrected Type V motif
coordinates, and added shell-free execution, process-group timeouts, quoted CSV
and versioned JSON/GFF3/FASTA output. These repairs can change v2 candidate counts
and ranks; v3 does not promise equivalence to erroneous v2 output.
See the [migration guide](docs/MIGRATION.md) for detailed behavior and CLI mapping.

## 2. Literature regression fixes

- **Long FASTA36 paths:** use short local filenames inside isolated working
  directories to avoid corruption by the legacy native filename parser.
- **Repeat copies hiding anti-repeats:** mask accepted arrays in a separate search
  sequence while retaining original bases and coordinates.
- **Malformed redundant searches:** remove the whole-genome-as-query FASTA36 pass;
  retain both repeat orientations in the regular direction.
- **Model-supported loci missed by sequence search:** add independent genomic
  covariance-model rescue with explicit origin and null unavailable evidence.

Completion rose from 54/79 to 79/79 inputs and localization from 45/78 to 68/78
mapped references in the inspected panel. This is development regression evidence,
not held-out accuracy. See the [benchmark log](BENCHMARK_LOG.md).

## 3. Interaction batching and reproducibility

- Batch IntaRNA targets sharing an exact oriented repeat, with bounded batch size
  and deterministic selection of equal-energy alternatives.
- Validate every returned row, including nonwinning results; reject malformed
  IDs, structures and coordinates while preserving valid missing interactions.
- Distinguish interpreter provenance from the package invoked through `python -m`
  without importing that package for inspection.
- Preserve the observed stock-FASTA36 reproducibility failure and provide a
  separate explicit-seed build for controlled one-thread comparisons.

Short anti-repeat rescue is a **new experimental feature, disabled by default**.
It is not a claim that weaker alignments are universally valid. Gates and controls
are in the [scientific plan](docs/COUNCIL_IMPROVEMENT_PLAN.md) and benchmark log.

## 4. Code quality and failure recovery

Four expert councils agreed on plans, separate agents implemented them, and a new
three-member council independently reviewed the integrated result.

| Area | Fix |
| --- | --- |
| Configuration | Reject numeric booleans, malformed weights and invalid path shapes before incidental Python errors |
| Coordinates and strands | Require genuine integer coordinates and valid strands; unknown strands no longer silently reverse-complement |
| Executable lookup | Resolve commands and relative PATH entries in the child directory; hash, record and execute the same symlink-preserving path |
| Failure diagnostics | Read a bounded stderr tail, retain full logs and preserve legitimate empty subprocess arguments |
| Upstream contracts | Strict JSON/JSONL decoding rejects nonfinite values and overflow; validate consumed IDs/counts and provide record context |
| Experimental policy | Define existing rescue thresholds once, retaining all values, comparison boundaries and public fields |
| Publication | Assign IDs before saving records, count only persisted records and preserve original errors when failure-marker writing fails |
| Atomic JSON | Use unique temporary siblings; clean abandoned writes without damaging an existing destination |
| Dependency updates | Resolve annotated tags to commits; stage locks; restore branch/HEAD on ordinary failures and describe failed restoration honestly |
| Dependency sync | Publish only completed checkouts, allowing retry after clone/checkout failure |
| Packaging | Verify installed wheels outside the source tree, including entry points, hashes, model indexes and executable permissions |

The final council found an additional defect: cleanup of a staged dependency lock
could hide the original error and recovery instructions. The executor fixed it,
added five regression cases, and the reviewer reran the original reproduction.
Cleanup errors now attach a diagnostic note to the primary failure.

**Final checks:** 360 unit tests, seven native-tool tests, strict mypy, Ruff and
formatting passed. The source distribution passed the same 360 unit tests, and
the installed wheel verified 27 resources. All 104 native prediction records in
five comparison reports remained unchanged.

See [the final council outcome](docs/code_quality/08_final_outcome.md) for plans,
implementation records, the initial review finding and independent sign-offs.

## Compatibility and limits

The latest quality round preserves valid scientific behavior and final report
schemas; earlier rejection of malformed values is intentional. Per-record JSON
now contains the same export IDs as the final report. Public genomic coordinates
remain one-based and inclusive.

Dependency recovery handles ordinary caught failures, not crash-atomic
Git/filesystem transactions or concurrent updates. Scientific ranking, transcript
boundaries, model coverage and native reproducibility retain the limitations
explained in the README and benchmark log.

## 5. Scientific report boundary audit (2026-10-11)

- Reject duplicate JSON keys at every nesting level in upstream reports and
  JSONL records. Conflicting status, coordinates or scientific evidence must
  fail rather than silently select the final field.
- Require numeric CRISPRidentify certainty scores. Booleans and numeric strings
  no longer become array support through Python's float conversion; valid integer
  and floating-point scores retain their values.
- Reject repeated CasAndra cassette IDs even when the supplied counts and
  checksums are consistent. A cassette reference must identify one evidence row.
- Reject Infernal hits whose strand contradicts the direction of the reported
  sequence coordinates. Sorting endpoints must not conceal inconsistent strand
  evidence.

Eighteen additional regression cases preserve valid numeric confidence and cover
these malformed-report boundaries. All 378 portable tests and seven native-tool
integrations pass, together with lint, formatting and strict typing. Report
schemas, valid scientific outputs, CLI flags and defaults are unchanged; no models
or upstream dependency revisions were changed.
