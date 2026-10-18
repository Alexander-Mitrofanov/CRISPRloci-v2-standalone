# Input, reference and report contract

The runtime performs no network requests and launches no external scientific
processes. Python APIs use the same input validation as the CLI.

## FASTA

Inputs are bounded regular files; symlink traversal is rejected. ASCII FASTA
supports blank lines and multiline sequences. The first header token is a
unique ID, preserved up to 500 printable nonwhitespace characters. Descriptions
are retained up to 2,000 characters. Sequence case is normalized, U becomes T,
and allowed IUPAC bases are checked. Internal sequence whitespace is invalid.

The parser caps input at 64 MiB, 100,000 records, 10 million bases and 200 bases
per record. Search additionally caps queries at 10,000 occurrences. Ambiguous
queries are `unsupported`; ambiguous reference entries remain in the manifest
but are excluded from matching. Their exclusion is reflected in eligible counts.

## Reference schema 1.0.0

`schema.name` is `crisprmap.reference`. `provenance` requires `namespace`,
`source_url`, and `source_sha256`; values are bounded strings or null. `entries`
bind unique IDs, descriptions, normalized sequences, sequence checksums and
nullable `family`, `motif`, `accession`. Unexpected reference fields are rejected.

`content_sha256` covers the complete reference except that checksum field,
serialized as sorted compact ASCII JSON with a terminal newline. It covers
annotations as well as sequences. Load rebuilds the reference and compares
exactly; duplicate JSON keys and nonfinite numbers are errors. These checks
establish content consistency. Trust in biological annotations comes from
external provenance, not the mere existence of a checksum.

## Search schema 1.0.0

`schema.name` is `crisprmap.search`. A successful result contains
`execution_completed: true`, method, reference identity, policy, query/work
counts and ordered `results`. CLI reports additionally include software version,
input byte SHA-256 and reference file SHA-256/provenance. No paths identifying the
local workstation are stored in result metadata.

Each result contains the normalized query and sequence SHA-256, status/reason,
minimum edit distance **within the configured retrieval limit**, tied best
occurrence count, displayed count and truncation flag. `no_hit` distance is null:
the engine does not claim a globally nearest distance outside the search limit.

Annotation fields each contain `status`, all known best-tie `values`, nullable
assigned `value`, and `missing_reference_count`. Statuses distinguish
`historical_association`, `candidate` (family only), `reference_evidence` (near
motif only), `conflicting`, `incomplete`, `unannotated`, and `not_assessed`.
Only unanimous complete exact associations populate `value`. Labels on near-hit
motifs are never query motif assignments.

Hits identify reference occurrence, family/motif/accession evidence, distance,
normalized edit similarity, and one deterministic complete alignment.
`orientation` is forward/reverse/both relative to the stored reference.
For both-orientation ties, `aligned_orientation` declares which optimal
orientation is displayed. Traceback prefers diagonal, then deletion, then
insertion. Alignment columns are strings with `-` gaps, not genome coordinates.

Results are deterministic apart from user-chosen policy values. All closest ties
within the limit are pooled before a deterministic ID-sorted `top_k` display cap.
Caching repeated sequences does not collapse query occurrences. Comparison
reservation is conservative; the actual executed comparison count is also given.
The complete serialized search report is capped at 32 MiB during assembly;
repeated queries cannot bypass the cap through caching. Exported artifacts
together are capped at 128 MiB.
Deadline/work/output limits fail the entire operation rather than silently
returning a partial success.

## Artifact publication

`report.json` is the completion artifact, published atomically after HTML and TSV
are written. Destination directories must be new; existing results are refused.
Normal failures clean up newly created artifacts. A killed process may leave an
incomplete directory without completion JSON; choose a new path when retrying.
Reference construction also publishes exclusively and atomically.

HTML escapes caller-supplied content, has no scripts/CDN assets and uses a
restrictive content policy. TSV neutralizes leading spreadsheet formula
characters. JSON preserves original IDs without spreadsheet escaping.
