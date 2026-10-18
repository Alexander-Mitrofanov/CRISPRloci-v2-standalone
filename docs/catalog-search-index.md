# Building the catalogue's exact-search index

The scientific catalogue is immutable. A separately published SQLite sidecar
maps exact search terms to candidate shards, allowing the API to find repeats,
observed spacers and Cas labels across the complete snapshot without scanning
the corpus during a request. The reader verifies the actual scientific values
inside each candidate shard before returning a record.

From this repository, after publishing a complete catalogue:

```bash
backend/.venv/bin/python scripts/index_catalog.py \
  .catalog/snapshots/crisprcasdb-20260927 --workers 4
```

This administrator command requires no network access or predecessor checkout.
It never changes `catalog.sqlite3`, `manifest.json` or scientific shard files.
It creates `exact-search.building.sqlite3`, then atomically renames that file to
`exact-search.sqlite3` only after all checks succeed. The API automatically
detects that sibling filename; `CRISPR_WEB_CATALOG_SEARCH_INDEX` can select an
explicit location instead. An already published sidecar is not overwritten.
On a host with sufficient available RAM, `--cache-mib 3072` allows a 3 GiB
offline SQLite cache to reduce indexing I/O. The default is 512 MiB and the
explicit maximum is 4 GiB; this option does not alter the API's memory limits.
For a large catalogue, `--commit-batch 1000` reduces repeated index-page writes
by sorting up to 1,000 shards' candidate rows per transaction. The default is
50; this affects transient operator memory and repeated work after interruption,
but preserves the same atomic receipts and search semantics.

## Search meaning

The builder and API share
[`catalog_terms.py`](../backend/src/crispr_web/catalog_terms.py). Changing these
semantics requires a new term version and rebuilt sidecar.

| Kind | Indexed source values |
| --- | --- |
| 1: arrays | Repeat consensus and source array call identifier. |
| 2: repeats | Observed repeat sequence. |
| 3: spacers | Sequence of an observed spacer; deletion boundaries supply no sequence. |
| 4: Cas systems | Source subtype. |
| 5: Cas genes | Complete annotation and profile labels, plus explicit gene-name tokens in those labels. |

DNA strings containing IUPAC nucleotide symbols are normalized to uppercase.
This does not interpret ambiguity symbols as wildcards, perform similarity
search, trim a sequence, or match its reverse complement. Cas labels are
case-insensitive. For example, `C25_Cas5_cmr3_1` supplies its complete label and
the tokens `cas5` and `cmr3`. A variant such as `Cas5f1` does not supply `cas5`.
These are source-label aliases, not inferred biological classifications.
Gene call identifiers are not indexed globally. Genome accessions use the
existing exact accession registry.

## Searching from array results

The Arrays results tab exposes observed repeat and spacer sequences under each
array. Each reported DNA sequence has a **Search in DB** button; the same action
is available in database array details. Missing sequences and deletion records
cannot be searched, and consensus sequences are never substituted for missing
repeat occurrences. Newly generated summaries include these sequences; older
saved summaries without them must be rerun to populate the controls.

The repeat and spacer endpoints accept `match=similar` with a DNA `q` of 1–1,000
bases. This opt-in search finds equal-length sequences with **at most one literal
nucleotide substitution**, considering both the query and its reverse complement.
IUPAC symbols are literal symbols, not wildcard matches. Insertions, deletions,
local alignments and more distant homologues are outside this search's scope.
The default `match=exact` behavior remains unchanged.

The reader probes the existing index for the query's bounded one-substitution
neighbourhood, keeping only the next five distinct candidate shard IDs. Every
probe reads at most five index entries; probing has a two-second and two-million
SQLite-instruction budget. At most four scientific shards are verified per page
under the existing memory, time and response-size limits. Search never changes
or builds the index. A missing or mismatched sidecar returns an unavailable error,
not a negative biological result. Cursors bind the snapshot, query, filters and
match mode; even empty pages may have a continuation cursor.

Results retain every matching occurrence, in database order. The interface shows
the source sequence, accession, contig, array and source coordinates, together
with full-length literal identity, substitution count and relative strand. These
are sequence-comparison values, not confidence or inferred biological function.

## Storage and bounded execution

The `hits` table stores unique `(kind, key_hash, shard_id)` rows as its primary
key without an additional row identifier. `key_hash` is the first 16 bytes of
SHA-256 over the normalized UTF-8 term. Repeated occurrences of a sequence in
the same shard share one candidate row, while every original occurrence remains
in the scientific shard. Hashes nominate candidates only: a collision cannot
establish an exact scientific match.

Metadata pins the schema `crisprloci.exact-search.v1`, shared term version,
snapshot identifier, full registry SHA-256 and ordered shard-inventory SHA-256.
Only a `status=complete` sidecar with matching identity is served. The reader
checks at most four candidate shards per page and uses its bounded shard cache.

Up to four worker processes verify and parse independent shards. Only one
process writes the sidecar. It commits search keys and per-shard completion
receipts in bounded atomic batches, so interruption may repeat an uncommitted
batch but cannot retain its keys without its receipts. Running the same command
resumes committed work only when the snapshot and search semantics still match.
An exclusive lock prevents concurrent writers. SQLite uses full synchronous
transactions and bounded WAL checkpoints. The builder checks for a 2 GiB free
disk reserve plus conservatively estimated room for the current index's WAL
and new keys before each commit, and stops safely when that reserve is reached.

## Publication checks

Each source shard is checked against the registry's compressed SHA-256, declared
uncompressed size and every entity row count. Decompression is bounded, source
paths cannot traverse symlinks, and oversized scalar values are rejected. At
publication the builder verifies all source shards are represented, aggregate
source counts agree with the catalogue summary, candidate counts match committed
receipts, SQLite passes its integrity check, and the registry checksum has not
changed. It checkpoints and closes SQLite, synchronizes the sidecar, atomically
renames it, and synchronizes its directory. Source records and their identifiers
remain unchanged.

## Verified complete snapshot

On 27 September 2026 the builder published the index for
`crisprcasdb-20260927`: **43,665 source shards**, **55,547,358 candidate rows**,
and **1,632,821,248 bytes** (approximately 1.52 GiB). Its source totals matched
3,331,847 assemblies, 9,583,665 arrays, 55,023,435 repeat occurrences,
45,439,770 spacer records, 9,237,875 Cas systems and 220,117,236 predicted Cas
genes. Candidate rows are not occurrence counts: repeated source values within
one shard share a candidate entry, and one Cas label may supply several aliases.

The completed sidecar SHA-256 is
`bf17ea6f3289840c2a18c8514f5af3196ab0bae32bc27d5d04f697ad27f31aec`.
The pinned registry SHA-256 is
`f9aa6b9601f1b7bde7a4902984e08cbc5a9130a6c9816cb7b2513794daef08dc`.
All source checks, aggregate counts and SQLite integrity checks passed. Actual
global reader checks also returned exact repeat/spacer matches, source Cas
labels, `Cas9`/`cas1` aliases and a real 131-base spacer. The local verification
receipt is `.catalog/search-index-verification.json`; it includes metadata,
source totals and measured query outcomes. More than 4.7 GiB of free disk
remained after publication and verification.
