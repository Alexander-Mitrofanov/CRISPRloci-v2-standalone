# 2.1.1

- Reject duplicate native association/support records and alignments containing
  double-gap columns, while preserving distinct target intervals and occurrence
  projection. Reject estimated calibration without two target/control scores.
- Bind input/group hashes to the exact bytes parsed, and close catalogue handles
  deterministically on successful and failed association transactions.
- Parse multiline genomes in linear time, require an indexed catalogue accession
  key, cap directory inventories during enumeration, and reject fractional or
  boolean native resource parameters.
- Replace repeated linear FDR score lookups with binary search. Exact native
  output/calibration parity is checked against the prepatch implementation for
  calibrated, repeated-score and unavailable-control populations.
- Preserve reference and association schema v2, with explicit 2.1.0/2.1.1
  reference compatibility. Existing operator-provisioned 2.1.0 indexes need no
  rebuild; the executable and report identities advance to 2.1.1.
- Include upstream examples/assets and privacy scripts in source distributions,
  and pin the test toolchain and transitive packaging dependency for repeatable
  checkout and unpacked-source verification.

# 2.1.0

- Add broad immutable reference contract v2 with compact public provenance/counts,
  SHA-256/size-receipted indexes and a SQLite accession catalogue. Stream database
  checksums and load only matching source genomes in an eight-entry cache.
- Add command-length-safe genome-list/tar import, operator source metadata and
  complete reference verification. Preserve ACGTN target sequence content. Fail
  instead of truncating at explicit record, sequence and index-size ceilings.
- Bound operator native build steps separately from submitted search deadlines.
- Correct the vendored native `createdb` 16-bit source-index overflow above 65,535
  genomes by using 32-bit source indexes. Search/scoring/calibration algorithms
  are unchanged; regression tests cover all 65,537 native source mappings.
- Version association reports as 2.0.0 to replace the full target list with compact
  reference counts, source provenance, scope and manifest digest.

# 2.0.2

- Enforce aggregate 50 Mb reference coverage and exact complete database inventory at load, including native target/control base indexes. Reject a symlinked index directory and cap index storage at 1 GiB. Scientific algorithms and outputs are unchanged.

# 2.0.1

- Include every required upstream tracked regular source and matrix file in Git archives. The2.0.0 publication lacked files hidden by inherited ignore rules; its tag remains unchanged.
- Add upstream regular-file inventory and archive completeness regression.

# 2.0.0

- Modern C++17 first-party utilities, shared checked parsing and scoring modules.
- Typed offline Python package, explicit group maps, sequence-bound versioned results.
- Immutable target/control panel preparation with checksummed source and index inventory.
- Explicit unavailable calibration for undefined small/null-control edge cases.
- Tested calibrated native numerical parity and original biological support parity.
- Isolated bounded native execution and descendant cleanup; packaging and CI.
