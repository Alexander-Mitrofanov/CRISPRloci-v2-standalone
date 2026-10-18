# Changelog

All notable CRISPRidentify v2 changes are recorded here.

## Unreleased

### Added

- Add `--reconstruction_policy conservative` for balanced runs.
  Separate Bona-fide core confidence from reconstructed geometry: freshly
  score both, validate each candidate's source and applicable negative checks,
  preserve original core evidence, and abstain on ambiguity or accepted-array
  collisions. Selected geometry retains its own model score.
- Separate reconstruction selection, candidate validation and atomic category
  updates into focused modules, with replayable report evidence and a
  transactional decision sidecar.

### Changed

- Enable the benchmarked conservative reconstruction policy by default for
  `balanced`, including runs with no preset flag. Keep an explicit `legacy`
  opt-out and preserve the legacy policy for `full`, `region-screen`, and
  architecture-rescue workflows. Record either resolved policy in effective
  configuration and run provenance. The reconstruction sidecar now reports
  `experimental: false`; biological-accuracy establishment remains false.
- Document the complete 4,000-assembly comparison, including the substantial
  completeness gain, largely unchanged locus recovery, 37.3% higher observed
  native wall time than historical OLD Identify, and four memory failures.
- In conservative mode, preserve ordinary competition groups and restrict
  repair-only components that touch several ordinary groups. Retain conflicts
  without allowing repair bridges to merge ordinary groups. Preserve initial
  geometry-tie evidence in reports, retained candidates and reconstruction
  proofs without changing score bands or treating observations as selection
  authority.
- Recover uniquely supported omitted internal repeat copies with bounded,
  complete source-graph analysis and a checked single-core replacement.
  Preserve ambiguous boundaries, unknown restrictions, ordinary model bands
  and neighboring accepted events; expose replayable internal-completion
  evidence and actual search-work counters.

- Repair source-supported weak repeat boundaries before structural pruning,
  with deterministic source-only alignment and unchanged biological filters.
- Complete narrowly attested exact-repeat paths over accepted short cores
  when fresh scores meet the existing category thresholds. Preserve the old
  core as evidence and report the distinct
  `source_backed_terminal_repair_completion` classification basis.
- Bound and reuse immutable source-alignment computations while keeping
  candidate state, source validation and eligibility decisions independent.

## 2.3.2 - 2026-08-27

Likely-tandem soft-quarantine and report-integrity patch release.

### Added

- Add a source-replayed `likely_tandem` quarantine ledger for final accepted
  arrays with 2–8 canonical spacers when every spacer pair has global edit
  identity at least 0.95. The self-versioned report extension preserves the
  original score, category, source geometry, sequences, hashes, and
  complete-link evidence.

### Changed

- Move qualifying `Bona-fide` or `Possible` rows to the public `Low score`
  category before circular-origin analysis and expose the provenance-bound
  `tandem` tag in canonical JSON, HTML, low-score text, legacy JSON, and pickle
  output. This is an intentional scientific-output change; the optional field
  keeps the public report schema compatible, and the distinct scientific
  disposition remains `likely_tandem`.
- Preserve the likely-tandem ledger under accepted-only projection and bind
  every full-report ledger record one-to-one to its tagged Low-score array.
  The writer exhaustively reruns the transaction over every pre-quarantine
  accepted row, replays positive records against the original FASTA and exact
  core architecture, and publishes no legacy or canonical artifact on a
  mismatch. HTML and benchmark consumers independently rescreen accepted rows.
- Invalidate stale `report.json` completion markers before the first pipeline
  stage and publish the replacement strictly last.
- Render a visible tandem-demotion audit table in full and accepted-only HTML
  reports.
- Add `tags` to CRISPRidentify-300 call tables without changing the frozen run
  envelope version; union and delta-audit readers accept historical tagless
  tables as well as new tables.

### Fixed

- Reject boolean-for-numeric substitutions and incomplete transition
  provenance in canonical, HTML, benchmark, and batch/resume consumers.
- Preserve historical two-field category rows when a qualifying row is padded
  with tandem metadata and projected into a full report.
- Validate every batch completion against source-reconstructed core arrays and
  the complete likely-tandem ledger before accepting it as resumable output.

## 2.3.1 - 2026-08-26

Archive-first documentation patch release.

Runtime code, models, dependency locks, report schemas, benchmark evidence, and
scientific behavior are unchanged from v2.3.0.

### Changed

- Replace clone-based setup with exact archive verification, extraction,
  installation, and smoke-test instructions.
- Record the release version and archive SHA-256 as the supported source
  identity for citations and support requests.
- Document how to select the optional IS-element environment after installation.

## 2.3.0 - 2026-08-24

Scale-oriented screening, restartable batch execution, and an audited release
contract.

### Added

- Add named `balanced`, `region-screen`, and `full` presets plus effective
  configuration printing and machine-readable run provenance.
- Add a bounded, restartable batch runner with deterministic sharding, retries,
  timeouts, atomic completion records, and failure ledgers for multi-VM runs.
- Publish the report JSON Schema, output/category contract, cluster guide,
  measured performance evidence, and checksum-pinned 100/1,000-genome cohort
  manifests.
- Add citation, contribution, security, support, conduct, and GitHub workflow
  documentation.

### Changed

- Run the complete standalone test suite in CI, add first-party Ruff checks,
  and pin every third-party workflow action to an immutable commit.
- Make high-volume console telemetry opt-in while retaining concise progress
  and explicit effective-configuration output.

### Fixed

- Make competitive-component splitting default to enabled so the documented
  high-recall fast profile is internally consistent with its default bounded
  rescue. The previous contradictory defaults caused every ordinary CLI run
  to stop during argument validation.

### Performance

- Add native path-feature ABI v4. Successive exact interfaces aggregate
  first-encounter-ordered counts, construct normalized packed sparse rows,
  perform sparse dot products, and build/hash ASCII path neighborhoods in C,
  while retaining bounded native and Python fallbacks.
- In the machine profile, skip source-phase geometry scores when every
  competing geometry has already failed the mandatory boundary-evidence gate;
  this preserves the candidate group and accepted output by construction.
- On the paired speed-100 benchmark these changes reduced wall time from
  381.78 s to 270.94 s (1.409x throughput) and CPU time by 33.77%, with all
  164 `report.json` files byte-for-byte identical.
- Bulk-copy packed native rows, use mutation-sensitive owned fuzzy-match
  snapshots, cache bit-exact Jeffreys source-phase scores by canonical base
  counts, and scan FASTA symbols as exact sets. These preserve accepted output
  while removing high-frequency Python allocation and scoring work.
- Add explicit `--fuzzy_boundary_paths {True,False}` with a default of `True`.
  On speed-100, `False` completes in 124.95 s versus 270.94 s for the accepted
  exact run (2.168x), retaining overlapping calls for 555/560 exact-fast loci.
  Because it loses five Possible arrays, downgrades one Bona-fide locus, and
  changes other boundaries, it is an opt-in screening tradeoff rather than an
  exact default.
- On the frozen release-10 inputs, the current boundary-on fast stack takes
  25.14 s (19.38x versus the 487.16 s historical full reference); boundary-off
  screening takes 19.82 s (24.58x). Both ratios are same-input cross-version
  agreement benchmarks, not biological sensitivity estimates.

## 2.2.0 - 2026-08-21

High-recall fast-default and machine-output scale release.

### Changed

- Make the balanced high-recall fast path the CLI default. Historical full
  search remains available explicitly with `--fast_run False`.
- Add `--output_profile machine`, which transactionally publishes only the
  canonical `report.json` required by database builders and removes stale
  presentation artifacts.
- Increase order-preserving native fuzzy batches from 32 to 128 patterns and
  omit Vmatch index tables unused by the `-s leftseq` query.

### Performance

- Replace general tuple construction/hash dispatch in compact path
  vectorization with exact Python-3.7-compatible three- and four-integer
  low-16-bit recurrences.
- Deduplicate and cache identical RNAfold and BLAST scoring requests in bounded,
  process-local, failure-atomic LRU caches.
- On the final frozen ten-genome comparison, default fast mode completed in
  48.25 s versus 487.16 s for historical normal mode (**10.10x**), at
  169,108 KiB peak RSS. It was byte-identical to the prior explicit
  high-recall fast profile; preceding final-candidate runs were 49.90--52.71 s.

### Correctness and limits

- Stream Vmatch result records and request one sentinel above the one-million
  accepted-record ceiling. Repeat-dense inputs now fail explicitly if the
  candidate domain exceeds that ceiling instead of silently accepting
  `-best` truncation.
- The ten-genome normal-mode comparison is an agreement test, not an
  independent truth set: fast mode matched 39/46 accepted normal loci, including
  6/6 Bona-fide and 33/40 Possible loci; 19/46 borders were exact. These numbers
  are release evidence, not a claim of equal recall.

## 2.1.2 - 2026-08-19

CasAndra cassette-hierarchy contract hotfix.

### Fixed

- Normalize the independent final cassette class/type/subtype winners after
  both direct-profile aggregation and architecture-model classification: type
  remains the winning type, class is derived from that type, and an
  incompatible subtype is masked rather than replaced by an unreported
  lower-ranked label.
- Preserve each provider's unnormalized final cassette claim in
  `raw_classification` evidence, mark repaired cassettes with
  `hierarchy_conflict_masked`, and report a validated per-run cassette repair
  count for cohort-scale quality control.
- Apply the genome evidence gate to the normalized cassette claim and
  independently replay the raw-to-normalized transformation, route evidence,
  gate result, and repair count in the CRISPRidentify adapter and report
  validator.
- Cover the observed direct-profile Class 1/Type V and Class 1/Type VI
  conflicts, the architecture route, and adversarial provenance claims with
  regression tests.

### Compatibility

- The bundled CasAndra package is now
  `0.1.0.dev0+crispridentify.3`; model artifacts and schema version 4 are
  unchanged.

## 2.1.1 - 2026-08-19

CasAndra protein-hierarchy contract hotfix.

### Fixed

- Normalize independent protein class/type/subtype head winners into a
  hierarchy-consistent classification: type remains the winning type, class is
  derived from that type, and a conflicting subtype is masked instead of being
  replaced by an unreported lower-ranked label.
- Preserve the unnormalized head winners in canonical per-gene
  `raw_classification` evidence, mark every repaired record with
  `hierarchy_conflict_masked`, and report a validated per-run repair count for
  cohort-scale quality control.
- Recompute the normalization and conflict flag inside the CRISPRidentify
  adapter so malformed output and false provenance claims remain fail-closed.
- Cover the observed type-I/V-M and type-I/V-F conflicts with regression tests.

### Compatibility

- The bundled CasAndra package is now
  `0.1.0.dev0+crispridentify.2`; model artifacts and schema version 4 are
  unchanged.

## 2.1.0 - 2026-08-18

CasAndra-default provider-contract release.

### Added

- Vendored, offline CasAndra runtime and frozen deployment model.
- Versioned Cas provider metadata, source-coordinate cassette/gene evidence, and
  CasAndra profile-family evidence in machine-readable reports.
- Public `--cas_gene_mode {auto,single,meta}` control.
- Fail-closed CasAndra bundle and nested-artifact integrity verification.
- Reproducible source/model provenance, full Apache-2.0 license material,
  downstream modification notice, and model/reference-data attribution.
- Exact CasAndra dependency pins in both supported Conda profiles and the Linux
  core lock.
- Pinned Conda-forge GCC 14.4 toolchain, standalone install/launch wrappers, and
  clean-clone GitHub Actions coverage for repeat setup and CasAndra smoke tests.

### Changed

- `--cas True` now uses CasAndra by default instead of
  CRISPRcasIdentifier.
- Setup now reconciles existing environments, installs the vendored package
  offline and non-editably, runs `pip check`, builds the native matcher, and
  runs an expanded environment/model doctor.
- Setup now rejects unsupported platforms before solving and builds the native
  matcher only with the compiler supplied by the active locked environment.
- The launcher now resolves relative input and output paths from the caller's
  working directory, including when invoked by absolute path.
- CRISPRidentify v2 is a tracked product in this repository instead of an
  ignored local export.
- Type VII candidates with no supported subtype can pass a type-aware genome
  evidence gate without weakening gates for other types. These type-only calls
  are marked experimental because independent whole-genome Type VII cassette
  validation is not yet available.

### Compatibility

- Cas annotation remains opt-in through `--cas True`; only the provider used
  by that option changed.
- The supported runtime is Python 3.12 on Linux x86-64.
- CasAndra uses a downstream local version suffix to distinguish the integrity
  and output-contract changes from the upstream `0.1.0.dev0` snapshot.

## 2.0.0

- Introduced the standalone v2 runtime with the compact ten-feature array
  classifier, portable strand prediction, source-coordinate JSON/HTML reports,
  and reproducible local setup.
