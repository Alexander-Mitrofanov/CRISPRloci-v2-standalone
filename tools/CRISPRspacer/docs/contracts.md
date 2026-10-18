# Version 0.3.1 contracts

Every JSON artifact identifies its schema. Runtime validators are implemented
in the relevant input, artifact, database and viral modules; incompatible schemas fail explicitly.
Native reports are parsed into source-bound records without importing a detector.

The 0.3.1 patch retains all schema versions, scoring policies and defaults.
Replay distinguishes boolean values from numeric counts and scores. Corrupt
deflate payloads produce structured input errors. Indexed searches validate
normalized target identities and require linear traversal. SQLite metadata is
size/type bounded before Python materialization, and viral catalog SQL honors
the operation's deadline/cancellation checkpoints. Existing valid 0.3.0 results
remain replayable; indexed results from 0.2.0 are also accepted when replay agrees.

## Viral references and reverse search (added in 0.3.0)

Operator acquisition uses `crisprspacer.ncbi-viral-plan/1.0.0` and
`crisprspacer.ncbi-viral-download/1.0.0`. A metadata plan freezes the complete
genomic member inventory and official release catalog. Downloads verify official
MD5, add SHA256, preserve acquisition metadata and never replace existing output.

`crisprspacer.viral-reference/1.0.0` binds source receipts, normalized FASTA,
SQLite sequence/metadata inventory and every native BLAST database file. Its
identity is SHA256 of `reference.json`. Native volumes/aliases stay inside the
reference directory. Operator validation roundtrips every native subject against
the normalized inventory. Operators make the validated directory immutable;
searches verify its complete file inventory and hashes before use.

`crisprspacer.viral-search/1.0.0` is a distinct heuristic report with
`execution_completed: true` and `retrieval_completeness: "heuristic"`. Its manifest
binds `result.json`, `inputs.json`, `matches.tsv`, `blast-hsps.tsv` and
`blast.stderr.txt`. Query occurrences retain repeated original identifiers and
sequence duplicates. Eligible queries are unambiguous ACGT, 18–80 bases.

The recorded `blastn-short` profile searches both strands with word size 7,
masking disabled, E-value 1000, ungapped scoring reward 1/penalty -3, and a subject
limit equal to the reference record count. No native per-subject HSP cap or
best-hit/culling filter is used. Acceptance requires each individual HSP to span
the full query with 0–2 substitutions (default 2); source coordinates and aligned
bases reconstruct exactly from the pinned subject catalog. Native retrieval may
miss matches or clip terminal mismatches. Output/deadline/HSP caps fail the entire
publication; they never convert a partial run into a no-match result.

Counts distinguish input occurrences, unique sequences, accepted HSPs and
per-accession locations. Candidate summaries are sequence associations only.
Taxonomy, infection, host range, PAM and biological targeting remain unassigned.
Validation replays retained HSP filtering, source reconstruction and summaries;
it does not rerun native retrieval or prove omitted HSPs never existed.

These commands remain independent of the CRISPRloci v2 web adapter, which uses
isolated workers and its own independently validated, sanitized result contract.
Standalone executable provenance may contain operator paths.
See [methods and reproduction](viral-reference-methods.md).

## Indexed spacer collections (added in 0.2.0)

`database-index` publishes `spacers.sqlite3` and `reference.json`. All spacer
occurrences remain in the SQLite table; eligible occurrences have one forward
and one reverse seed, or one `both` seed for palindromes. Construction and
`validate-database` prove the entire seed inventory corresponds exactly to the
spacer table. Operators validate, pin the database SHA-256, and make the snapshot
immutable before exposing it to workers. Per-job use checks those pinned bytes
and the internally bound reference metadata. Only operators can provision data.

`crisprspacer.database-search` schema `1.0.0` is a separate exact-only report:
tool identity, `status`, `search_complete`, `outcome`, `reference`, fixed `policy`,
target identities, counts and occurrence-expanded `matches`. Coordinates are
1-based inclusive on linear forward targets. The `inputs.json` bundle contains
normalized targets and the database identity, never the whole database.
`manifest.json` schema `1.0.0`, kind `crisprspacer.database-search`, binds exactly
`analysis.json`, `inputs.json` and `matches.tsv`. Validation replays against an
explicitly supplied matching database. TSV strings are spreadsheet-escaped.

Database matching nominates candidates by a 16-base exact prefix, then checks
the entire oriented spacer. Every eligible exact match necessarily has that
prefix. Targets are scanned in bounded batches of 400 starts; rows are streamed
from SQLite. Both Python work and SQLite execution have deadlines. Resource
errors do not publish completed output. Prefix ambiguity counts are diagnostics,
not per-spacer coverage. No PAM or array-overlap inference is performed.

See [source, eligibility, limits and deployment instructions](ncbi-spacerome.md).

## Inventories

`crisprspacer.inventory` version `1.0.0` contains `sources`, `arrays`, `queries`,
`annotated_records` and input `provenance`. Each source has unique `id`,
uppercase IUPAC `sequence`, `length`, normalized `sha256`, nullable `assembly_id`
and traversal `topology`. Arrays have inventory-scoped ID, original native ID,
source ID, inclusive interval, source category and nullable orientation.

Queries are occurrences, with unique ID, observed/deletion kind and source
metadata. Source-bound observed queries reconstruct exactly from their reported
forward-source slice. A deletion has `sequence`, `sha256`, `start` and `end`
all null, and a 0-based interbase `boundary`. Spacer FASTA queries have no
invented source record, array, ordinal, interval or category.

The CRISPRidentify adapter supports schema `1.2.0`, forward-source coordinates
and `uppercase_sequence` SHA-256 scope. It validates the complete reported
repeat/spacer architecture, including native counts and ordinal gaps. It retains
all reported array intervals independently of query category selection. It does
not assert that supplied reports include every real array or that reported
boundaries/classifications are biological truth. Native report digests preserve
identity; proprietary filesystem paths and arbitrary native extensions are not
copied into output. Inputs are consumed from the same byte snapshot that is hashed.

## Search semantics

`crisprspacer.analysis` version `1.0.0` records tool version, complete effective
policy, input inventory and target digests, execution/outcome, query statuses,
per-pair coverage, every alignment and aggregate counts. The retrieval profile
is exhaustive full-length Hamming comparison over eligible ACGT queries and
unambiguous target windows. Default `substitutions=0`; supported values are 0–2.

For `k < query_length`, partition the query into `k+1` nonempty blocks. A window
with at most `k` substitutions must contain one block unchanged. All exact seed
occurrences nominate starts; candidates are deduplicated and verified across
the full query. For `k >= query_length`, enumerate all starts. Perform this
independently for the query and its reverse complement, merging palindromic
orientation equivalence only. Caps fail the whole publication, so successful
retrieval is not silently truncated.

Origin-crossing targets use ordered nonwrapping `target_segments`. Start/end
refer to traversal endpoints and can have start greater than end. Their strand
is relative to the stored spacer sequence, independent of source-array
orientation. Left/right flanks are in the forward target frame, and circular
flanks never exceed the nonaligned remainder of a lap.

`same_source_record` can be assessed only by the source-target workflow; a name
or hash match in external inputs is insufficient. `same_declared_assembly`
compares explicit user declarations and is null if either is missing. These
fields do not infer organism or sample identity.

Array classification is `source_array`, `other_array`, `array_padding`,
`outside_reported_arrays`, or `unknown`. The first three are excluded from the
candidate table and retained in exclusions. All reported array categories are
used for overlap; `--array-padding` is a declared user policy, not a universal
scientific cutoff. Missing annotations leave overlap unknown.

PAM status is `not_evaluated`; score/probability fields are null. There is no
aggregate host assignment, guide design, functional-targeting decision or
implicit CRISPR subtype-specific mismatch rule.

## Reference identity and verification

`crisprspacer.reference` version `1.0.0` is an immutable local normalized
sequence inventory with source-file and normalized-inventory hashes, packaging
time, a declared source/release label, record/base counts and normalization
policy. It is not a seed index or database downloader. The exact normalized
targets and reference identity travel with results.

Manifest schema `1.0.0` fixes bundle kind and exact filenames, sizes and hashes.
Unexpected files, symlinks, unsafe paths and altered bytes are rejected. Bundles
have an aggregate read ceiling of 512 MiB. Verified bytes are the bytes parsed;
a changed file is not reopened as if it were the verified snapshot.

Result validation first checks artifact integrity, inventory source binding and
direct alignment reconstruction. It checks retained reference identity against
the target inventory, then replays deterministic discovery/classification and
compares the full report and regenerated TSVs. This catches omitted/altered
matches even if a manifest was regenerated. It cannot authenticate a deliberately
replaced, internally consistent input collection without an independently
trusted original digest. Replay currently requires this report-producing tool
version and obeys the retained resource policy.

## Deferred work

- Broader performance evaluation beyond the declared collection and query budgets.
- A separately versioned adapter for additional array-caller schemas.
- Optional calibrated partial/gapped acceptance profiles.
- Scientifically validated PAM/context annotation and aggregate associations.
- Reference acquisition for additional collections and separately pinned taxonomy metadata.

These are not silently approximated by version 0.3.1.
