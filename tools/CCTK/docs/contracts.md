# Structured observed arrays

Input schema is exactly `{"name":"cctk.arrays","version":"1.0.0"}` and an
`arrays` list. Each array has `id`, `source_id`, `array_id`, `strand` and `spacers`.
The execution ID is unique; `(source_id,array_id)` is also unique. Strand is
`null`, `+` or `-`, and never changes the supplied sequence or order. Each spacer
has `id`, contiguous 1-based `ordinal`, `sequence` and `status`. `observed` requires
1–1000 A/C/G/T bases; `deletion` requires null sequence. Occurrence IDs are unique
within an array. No extra fields are accepted. Source names are annotations,
not verified host identity.

Limits: 200 arrays, 500 slots per array, 10,000 slots total, 32 MiB input/output
files and 20 million plot pixels. Complete output must also fit the 32 MiB file budget; oversized reports fail explicitly.
Missing observations must be repaired at source; they cannot be quietly skipped.

Result schema `cctk.comparison` version 1.0.0 carries tool identity, exact input
SHA-256, options, echoed arrays/slots with derived DNA-equivalence symbols,
retained edges with supporting symbols, connected components and isolated nodes.
A symbol is `s_` followed by SHA-256 of the exact DNA (or lexicographically smaller
of DNA/reverse-complement under the explicit alternative policy). Deletions have
null symbol and participate in no set metric. Duplicate sequences may share a
symbol but retain distinct occurrence identities and positions.

`min_shared` is an integer 1–500; `identity` is `exact` or `reverse-complement`;
`order` is `supplied` or `id`. Network edges exist only when the set intersection
meets the positive threshold. A zero-overlap pair is never an edge. Components
include every array once; transitivity does not imply pairwise sharing.

Each array reports `observed_spacer_count` and `evidence_status` (`observed` or
`no_observed_spacers`). An all-deletion array remains a node with no DNA evidence
and cannot contribute an edge. Plot labels longer than 33 characters are
abbreviated; full identities remain in JSON. Above 10,000 rendered connections
both plot interfaces omit connector lines with an explicit note; the JSON/TSV
network remains complete. The observed viewport reserves all slots and labels.
