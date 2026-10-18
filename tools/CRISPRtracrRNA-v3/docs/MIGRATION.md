# V2 to V3 migration and bug audit

The baseline is BackofenLab/CRISPRtracrRNA tag 2.0. All v2 tracked files remain unchanged
in the sibling `old_version` checkout. V3 replaces the old 1,735-line modules directory with
an installable `src` package, immutable domain values, validated configuration, explicit
orchestration, typed tool adapters, and separate evidence/report layers.

## Corrected behavior

| V2 issue | V3 behavior |
| --- | --- |
| Input FASTAs concatenated, including multiple contigs | Parse records independently and preserve full identifiers and source hashes |
| Accessions truncated at the first dot | Run-local opaque tool IDs map to full original contig IDs and file names |
| Anti-repeat slices omitted the inclusive final base | One internal half-open interval representation; explicit public conversion |
| Negative window starts wrapped to the end of DNA | Clip each window to its own contig boundaries |
| Reverse model-hit coordinates caused an extra strand flip | Normalize endpoints while retaining Infernal's explicit strand |
| Multiple-model hits flattened in model order and paired with wrong candidates | Key every hit by candidate identity and model; select the best matching hit |
| Always-true string conditional disabled CMSCAN parsing | Parse the supported Infernal table, including all records and models |
| Overlap exclusion checked only hit start | Reject every anti-repeat interval overlapping any accepted array |
| Coverage used gapped alignment length and could exceed one | Divide aligned repeat coordinate span by repeat length |
| Strict `>` threshold unexpectedly rejected boundary values | Similarity and coverage thresholds are inclusive |
| Shared CWD temporary filenames and unchecked shell commands | Isolated per-record directories, argv commands, return-code checks, timeouts and logs |
| Exceptions swallowed during failed runs | Explicit failure artifact; publish completion only after all stages succeed |
| CSV assembled by string concatenation and parsed by splitting commas | Standard CSV readers/writers and structured JSON |
| Summaries could include stale or unrelated output files | Export only the current in-memory run; refuse existing output directories |
| No-hit cases could omit a summary header or crash | Stable header-only CSV and explicit empty result lists |
| `type=bool` parsed the string `False` as true | Explicit boolean parsing and TOML booleans |
| ERPIN parser consumed only the first hit per header | Parse all forward and reverse hits; keep all raw evidence |
| ERPIN E-value mislabeled as score and rewarded larger values | Preserve raw score and E-value separately; rank stronger significance higher |
| Missing evidence received consistency bonuses; containment overlap failed | Missing evidence has zero bonus; compute interval intersection consistently |
| IntaRNA and model coordinates had different upstream-window offsets | All evidence uses the same oriented candidate sequence with genomic projection |
| Type V motif flag did not call the search | Execute optional motif annotation in both modes |
| IUPAC regex classes included literal `|`; motif endpoints/distances were wrong | Valid nucleotide classes, overlapping motif search, inclusive endpoints, true intervening bases |
| `eval` used to decode interval text | No executable text parsing; typed intervals and structured files |

## Deliberate output and scientific changes

`report.json` is the authoritative v3 output. It has schema `CRISPRtracrRNA-report/3.0.0`.
The flattened `summary.csv` is an export, not a drop-in replacement for v2 columns. Complete-mode
predictions report the actual search window, its strand, all evidence, and the ranking feature
contributions. A search window does not establish mature transcript endpoints.

The same v2 covariance models and the default FASTA36 bidirectional searches are retained.
The following changes can alter candidate counts and ranks: CRISPRidentify v2's accepted
array set, full overlap exclusion, corrected endpoint/strand handling, inclusive thresholds,
span-based coverage, all ERPIN hits, and corrected consistency/evidence scoring. This is not
a promise of numerical identity to buggy v2 output.

Ranking keeps the ten configurable weight roles but defines the feature values explicitly:

- array score, identity, coverage, and identity × coverage;
- interaction: `max(0, -IntaRNA_energy / 10)`;
- covariance model: `max(0, Infernal_score / 10)`;
- terminator: `max(0, -log10(max(ERPIN_E_value, 1e-300)))`;
- orientation, anti-repeat/tail order and overlap, and tail/terminator order and overlap:
  Boolean consistency evidence, with no bonus when required evidence is absent.

The weighted sum is **uncalibrated**. Reverse hits inside an oriented candidate are retained
but cannot support that candidate's forward transcript hypothesis. Cas cassette proximity is
context only, filtered to Type II or V as selected; it is not a new evidence weight.

## CLI mapping

| Old | New |
| --- | --- |
| `python CRISPRtracrRNA.py` | `crisprtracrrna run` (source launcher remains available) |
| `--input_folder` | `--input` or `--input-folder` (old alias works) |
| `--output_folder` | `--output` (old alias works) |
| `--run_type complete_run/model_run` | `--mode complete/model` (old values work) |
| `--model_type II/V` | `--model-type II/V` (case normalized) |
| `--perform_type_v_anti_repeat_analysis True/False` | `--type-v-motifs [true/false]` |
| `--temp_folder_path` | Isolated intermediates always live under the run's `work/` |
| `--output_summary_file` | Fixed `summary.csv` inside the output directory |
| Ten `--weight_*` flags | Explicit `[weights]` table in TOML |

Removed options are rejected by argparse; they are never silently ignored. All examples use
fresh output directories. Dependency errors return a nonzero status and retain diagnostic logs.

## Council round-two additions

Existing configuration files remain valid. `intarna_batch_size` defaults to 128
and accepts positive integers; `1` selects serial execution with the same parser.
The tested native contract is IntaRNA 3.2.0. Compatible versions are accepted by
output contract rather than an arbitrary version rejection. Process timeout and
thread limits apply per batch; memory use grows with active IntaRNA threads.

Raw interaction evidence now uses `intarna-GROUP-CHUNK.csv`, target/query FASTAs
and `.ids.json` maps instead of one process/CSV per candidate. These artifacts are
retained and hashed. A valid missing native row means no reported interaction;
unknown IDs, invalid structures, duplicated/inconsistent records and malformed
nonwinning rows fail explicitly. Equal-energy alternatives are selected by target
start/end, repeat start/end and structure, making row order irrelevant.

Interaction records add one-based paired positions in the full oriented target
and repeat. Existing interval, energy and structure fields retain their meaning.
Predictions and CSV exports add `admission_branch` and `evidence_class`.
Candidate IDs remain stable. Export IDs remain rank-derived and may shift when an
experimental hypothesis is inserted into the ranking.

`short_anti_repeat_rescue = false` is the default. Enabling it requires Type II
complete mode; bundled Type V use is rejected. Custom models remain the user's
responsibility for Type II validity. It admits only complementary short proposals
with query span >=12 nt and identity >=0.80 below the ordinary coverage gate.
Native seed-constrained interactions must have energy <0, >=10 paired target
bases and >=50% paired target residues within the anti-repeat, canonical/wobble
antiparallel pairing, and forward downstream CM E<=0.01 or terminator E<=0.001
support. The exact geometry and frozen validation are in the council plan.
These are uncalibrated engineering hypotheses. Every proposal and rejection reason
is retained in `short_rescue_proposals`.

For `python -m` commands, provenance now distinguishes the interpreter checkout
from a package statically resolved through effective explicit `PYTHONPATH`.
Unresolved module paths are labeled; no module code is imported to identify it.

## Optional pinned native reproducibility profile

`local.reproducible.toml` selects the separately named
`native/fasta36-seeded-20260927/fasta36-council-seed20260927` binary and sets one
thread. No new biological threshold or model is introduced. This dependency build
corrects two missing RNG seed arguments in the exact upstream FASTA36 36.3.8g
source, using the council's already frozen master seed. Stock-FASTA configurations
remain usable but retain the observed shuffled-statistics variability.

The source archive, two-file patch plus typed header, compiler flags, build
attempts, license and checksums are retained alongside that native build. The
configuration and paths are local to this checkout; other machines must rebuild
and record their own dependency hashes. The wheel does not silently install or
replace this native executable. See the reproducibility amendment for the explicit
one-thread scope and the historical comparison that failed before this repair.
