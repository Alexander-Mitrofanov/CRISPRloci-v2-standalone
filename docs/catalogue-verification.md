# Full catalogue verification — 27 September 2026

Snapshot `crisprcasdb-20260927` was published locally at
`2026-09-27T16:29:41.917903+00:00` after the complete source transfer and verification.
No production service or release was deployed.

## Data

All **3,331,847 unique versioned accessions** reconcile with the completed
ClusterCRISPRCasDb campaign: **3,331,684 complete, 151 failed and 12 unavailable**.
Unsuccessful outcomes retain null feature counts. The imported records comprise
9,583,665 arrays, 55,023,435 repeat occurrences, 45,439,770 spacer occurrences,
9,237,875 Cas systems and 220,117,236 predicted Cas genes. Their coordinate
context comprises 98,367,031 feature-bearing source contigs. Assembly contig
counts retain the full upstream totals.

Every one of the **43,665 compressed shard files** passed SHA-256 verification
before atomic publication. All receipts reconcile with registry identifiers,
counts, sizes and digests. The 43,661 bundle receipts share the expected source
profile file digest; the four retained/raw groups carry the same pinned analysis
profile identity. Registry integrity and foreign-key checks pass. Array, system
and gene totals match the independent source summaries exactly.

| Identity | SHA-256 |
| --- | --- |
| Scientific profile | `5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540` |
| Ordered shard identity | `788048d4b340660f37b8324f221ff4318dc1f53ec89ea129593f9f2b10be4285` |
| Published registry file | `f9aa6b9601f1b7bde7a4902984e08cbc5a9130a6c9816cb7b2513794daef08dc` |

The registry is 412,823,552 bytes and the compressed shards total 23,416,865,282
bytes. The largest decoded shard is 13,692,928 bytes. Keep the registry, shards,
manifest and receipts together at
`.catalog/snapshots/crisprcasdb-20260927/`.

The snapshot excludes protein sequences, bulky native reports and feature-free
contig records. It preserves the requested feature annotations and all observed
repeat/spacer nucleotide sequences. Scientific limitations and coordinate
semantics appear in the manifest, API and Database view.

## Application checks

`bash scripts/check.sh` passed with **474 Python tests and 145 frontend tests**,
along with linting, formatting, source boundary checks, production build and
public export scanning. Tests use repository-local fixtures and require no
predecessor project or source cluster.

The independent full-snapshot API audit passed **128 GET checks** plus a rejected
POST, with a
1.672 ms median and 19.035 ms maximum in the local test environment. It checked
all entity routes and totals, pagination without duplicate records, accession
prefix paging, exact scientific searches, array/repeat/spacer relationships,
concurrent reads, and complete/failed/unavailable genomes across the snapshot.
All 163 unsuccessful outcomes have null scientific counters; a successfully
analysed genome with no features correctly retains zeroes. Invalid bounds and
cursors are rejected and write methods return 405. The published registry was
unchanged by these reads. These timings describe the local audit, not a
production load benchmark.

Browser acceptance against the complete published snapshot passed for all six
displayed totals, genome pagination, array details, repeat/spacer occurrences,
CSV and consensus FASTA downloads, Cas gene evidence, Back navigation and a
390-pixel mobile viewport. No JavaScript errors occurred. The temporary loopback
API/frontend processes were stopped after verification; no worker was started
and no analysis job was submitted.

Bounded-memory testing reached 667.805 MiB peak RSS with queued 254 MiB shard
reads and a 120 MiB cache, within the API service's 1 GiB limit. A separate index
comparison returned identical records/cursors for 335 scoped queries while
reducing sample compressed shard sizes by 15.87%.

Local evidence is retained in:

- `.catalog/final-check-all.log`
- `.catalog/audit/completed-snapshot-verification.json`
- `.catalog/audit/final-api-smoke-20260927T163257Z.json`
- `.catalog/audit/catalogue-memory-audit.json`
- `.catalog/audit/catalogue-index-audit.json`
- `output/playwright/catalog-full/acceptance.json` and adjacent screenshots

## Search-only interface revision

The follow-up interface uses one search-type selector, one input and an explicit
Search button. It displays no records, counts or empty-result state before a
valid submission. Editing the input, changing its type or choosing New search
clears previous output and cancels pending requests. A complete accession opens
genome-scoped context; repeat, spacer and Cas queries use the global exact-search
index. Provenance is collapsed beneath results.

The updated repository gate passed **489 Python tests and 149 frontend tests**.
Additional focused checks cover atomic index-batch interruption/resumption and
case-insensitive source DNA matching with literal ambiguity symbols. Initial
desktop/mobile browser acceptance verified the absence of record requests while
typing or changing search type, along with exact accession results, array
details, downloads and reset behavior. A real 131-nucleotide spacer is included
in the final global-search browser checks; input is explicitly validated at
1,000 characters without truncation.

Index construction and its snapshot identity checks are documented in
[catalogue search index](catalog-search-index.md). The original registry,
manifest and scientific shards remain unchanged.

The complete `exact-search.sqlite3` was published at
`2026-09-27T19:39:28.539273+00:00`. It covers all 43,665 shards, occupies
1,632,821,248 bytes, and reconciles every scientific source count. Its candidate
rows comprise 12,863,975 array terms, 6,721,695 repeat terms, 25,738,990 observed
spacer terms, 745,507 system terms and 9,477,191 Cas annotation/profile terms.
These are term-to-shard mappings, not biological occurrence counts.

The published-index API audit passed 14 HTTP checks covering exact global
repeat/spacer searches, Cas aliases and full profiles, a real 131-nucleotide
spacer, paging without duplicates, no-match and substring-negative controls,
oversized input and mismatched cursors. Registry and index file identities
remained unchanged. Evidence:
`.catalog/audit/global-search-api-smoke-20260927T193951Z.json`.

Final browser acceptance passed against the published index for all four search
types, original Cas profiles, the `Cas9` alias, and a real 131-base spacer with
pagination. The initial screen and input edits made no record requests. Both
desktop and 390-pixel mobile results passed overflow and JavaScript-error checks,
including long-sequence and Cas-gene results. Temporary verification servers
were stopped. Evidence: `output/playwright/catalog-search-final/acceptance.json`
and adjacent screenshots. The index verification receipt is
`.catalog/search-index-verification.json`; its file SHA-256 is
`bf17ea6f3289840c2a18c8514f5af3196ab0bae32bc27d5d04f697ad27f31aec`.

## Server configuration

Use a backend release containing this feature and set:

```bash
export CRISPR_WEB_CATALOG_DATABASE=/home/alex/CodexProjects/CRISPRloci-v2/.catalog/snapshots/crisprcasdb-20260927/catalog.sqlite3
```

The normal release launcher also accepts this path through `--catalog`.
Production installations should copy the whole snapshot to an operator-managed
read-only directory accessible to the API service, following
[the deployment instructions](../deploy/README.md). The reference data must
never be copied into the public frontend export. See
[catalogue setup](catalogue.md) and [API routes and limits](../backend/CATALOG.md).
