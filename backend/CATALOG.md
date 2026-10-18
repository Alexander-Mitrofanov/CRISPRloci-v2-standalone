# Public reference catalogue

The catalogue serves published CRISPRcasDb annotations independently of private
analysis jobs. Configure an operator-owned snapshot with
`CRISPR_WEB_CATALOG_DATABASE=/absolute/path/to/catalog.sqlite3`. Without it,
`GET /api/v1/catalog` reports `available: false`; job submission remains available.
Anonymous requests cannot download, import or replace a snapshot.

The overview reports imported entity counts separately from campaign coverage.
Failed and unavailable assemblies are analysis outcomes, not evidence that an
assembly lacks CRISPR arrays or Cas genes. Numeric IDs are snapshot-local; cursors
are bound to the snapshot, entity and filters. Retain the snapshot ID with exported
record IDs. Scores retain their source meaning and are not calibrated probabilities.

## API

- `GET /api/v1/catalog`: snapshot identity, coverage, provenance and limitations.
- `GET /api/v1/catalog/assemblies`, `/arrays`, `/repeats`, `/spacers`,
  `/cas-systems`, `/cas-genes`: `{items, next_cursor, limit}` pages.
- `GET /api/v1/catalog/assemblies/{accession}` and `/arrays/{array_id}`: one record.

Pages accept `limit` (1–100, default 25), `cursor`, `assembly_accession` and `q`
(at most 1,000 characters). Repeat/spacer pages also accept `array_id`. Follow
`next_cursor` until it is null, including when a page has no items. Child records
are paginated separately rather than silently truncated inside array detail.

Assembly search is a literal accession prefix; `assembly_accession` selects an
exact accession. Scientific queries match exact repeat consensus, repeat/spacer
sequence, array call ID or Cas system subtype. DNA letter case is ignored;
IUPAC ambiguity symbols must match literally and do not act as wildcards.
Returned source sequence text is preserved. Cas queries ignore case and match
the complete source annotation/profile or an explicit gene-name token in it:
`C25_Cas5_cmr3_1` supplies `cas5` and `cmr3`, while `Cas5f1` does not supply
`cas5`. Searches perform no family/subtype inference, alignment, reverse
complementation or similarity inference. The API exposes source coordinates as 1-based inclusive positions on
the forward source sequence, preserving unknown strand and all missing fields.
Repeats and spacers are occurrences; array consensus is a distinct field. A spacer
deletion remains `kind: deletion` with its reported boundary and null sequence.
Sharing a contig does not imply that a Cas system belongs to a particular array.

## Global exact-search sidecar

Global scientific queries use `exact-search.sqlite3` beside the immutable
registry, or the operator override `CRISPR_WEB_CATALOG_SEARCH_INDEX`. The offline
builder indexes distinct `(kind, first 16 SHA-256 bytes of normalized term,
shard_id)` triples in a `WITHOUT ROWID` table. It indexes observed repeat/spacer
sequences and source gene labels, not every Cas call ID. Source-label aliases and
DNA normalization are shared by the builder and reader in `catalog_terms.py`.

The sidecar metadata must mark a completed build and match the snapshot ID,
registry SHA-256, published shard-digest identity and term schema version. The
registry digest is cached only while its file identity remains unchanged. A
global query reads at most five indexed candidate IDs, visits at most four
candidate shards per page, and verifies every match against the original scalar
values or explicit source-name tokens. Hash collisions cannot yield false
scientific matches. Empty pages can still have a continuation cursor.

Missing, incomplete or mismatched indexes produce a generic 503 for global
scientific queries. Genome/array-scoped queries remain usable without the
sidecar. The public summary replaces the old genome-scope limitation only after
the sidecar validates; the published registry and manifest remain unchanged.

## Storage and safety

Small imports may use one normalized SQLite database. Comprehensive imports use
an accession registry plus checksummed gzip-compressed SQLite shards. The root
holds `export_metadata` (the `summary` JSON record), `assemblies`,
`assembly_shards(assembly_pk, shard_id)` and
`shards(id, relative_path, sha256, uncompressed_bytes, counts_json)`.
Shard row IDs use `shard_id * 1_000_000_000 + source_id`.

Only the offline importer writes data. Reads use immutable/query-only SQLite
connections, allowlisted columns and parameterized SQL. No source JSON, private
paths, commands or error text are exposed. Symlinks and traversal paths are
rejected. Compressed shards have a 64 MiB input limit, a verified SHA-256 and a
256 MiB expansion limit. Streaming decompression fills an exact-size buffer;
the complete compressed file is never held in memory. Up to two decompressed
shard images are cached within a combined 128 MiB limit. At most two shard
readers run concurrently, with weighted admission reserving twice the decoded
image size plus 8 MiB per reader within a 544 MiB active budget. This permits
256 MiB shards while queueing another large read for at most two seconds, and
keeps image allocations plus cache within 672 MiB under the API's 1 GiB service
ceiling. SQLite pager caches are limited to 1 MiB per connection. Each page
visits at most four relevant shards. SQL has both an instruction budget and a
deadline, and returned pages have a byte budget. A malformed, missing or
oversized configured snapshot produces a generic 503 without affecting the job
database.

Import and transfer instructions are maintained with the offline snapshot tools.
Tests construct synthetic data in this repository and never require predecessor
checkouts, network access or a scientific runtime.
