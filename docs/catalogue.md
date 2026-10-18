# CRISPRcasDb reference catalogue

The Database tab (`?view=database`) opens with a search form. Choose an accession
number, repeat, spacer or Cas gene, enter a value, and submit it. No record
requests, tables, totals or empty-result messages appear before submission.
Editing the query or changing its type clears results and cancels pending reads.

Accession searches require the complete versioned genome accession. Its results
link arrays, repeat/spacer occurrences, Cas systems and predicted Cas genes.
Repeat and spacer searches match exact nucleotide sequences across the imported
catalogue. Cas searches accept an original annotation/profile or a source-derived
gene name such as `Cas9`; they do not infer protein families from similarity.
Genomic context and array details open from matching records. Dataset provenance
is collapsed beneath results. CSV and nucleotide FASTA exports cover the
displayed page only. Searches do not perform alignment or sequence similarity
analysis.

## Imported snapshot

The complete `crisprcasdb-20260927` snapshot was fetched and verified locally on
27 September 2026. Its registry is
`.catalog/snapshots/crisprcasdb-20260927/catalog.sqlite3`; the adjacent manifest,
receipts and 43,665 compressed shards belong to the same immutable snapshot.

| Records | Imported count |
| --- | ---: |
| Genome outcomes | 3,331,847 |
| Arrays | 9,583,665 |
| Repeat occurrences | 55,023,435 |
| Spacer occurrences | 45,439,770 |
| Cas systems | 9,237,875 |
| Predicted Cas genes | 220,117,236 |

The compressed feature shards occupy 23.42 GB and the registry 0.413 GB
(decimal units), excluding receipts. Every shard checksum was verified before
publication. Registry integrity, foreign keys, unique accessions, terminal
outcomes, source profile identities and independent source feature totals
reconcile. The linked Database interface and read-only API have been exercised
against this full snapshot. See [verification evidence](catalogue-verification.md)
and the [council decisions](catalog-council.md).

## Source and scientific scope

The source is the completed **ClusterCRISPRCasDb / CRISPR Ensemble 0.6.2** campaign
over the frozen 2026-09-02 GenBank bacterial and archaeal catalogue
`c26a75991f66d0c8`. Its final 2026-09-22 reconciliation has **3,331,847** accessions:
3,331,684 successful, 151 failed and 12 unavailable. Failed and unavailable
genomes have unknown scientific counts, not measured zeroes.

Independent source summaries report 9,583,665 arrays, 9,237,875 Cas systems and
220,117,236 predicted Cas genes. These are feature occurrences across all source
assemblies, including related/redundant assemblies, not unique gene families or
unique biological loci. A full import must reconcile these counts exactly.

The pinned scientific profile is
`5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540`.
It uses CRISPRidentify 2.3.2 region-screen and independent whole-assembly
CasAndra. Accepted Bona-fide and Possible arrays are retained. Missing strands
remain unknown; source scores are not calibrated probabilities. Cas annotations
are predictions, including unassigned genes; array/Cas proximity does not imply
a functional association. The original system membership is preserved.

Coordinates are 1-based inclusive on the forward source sequence. Repeat
consensus is distinct from repeat occurrences. Spacer deletion boundaries are
distinct from observed nucleotide spacers. The compact catalogue preserves all
requested feature rows and their contig context, but omits feature-free contigs,
protein sequences, raw native JSON, commands and internal source paths from the
public API. Assembly contig counts retain the upstream totals. This is not a
complete genome or protein FASTA repository.

## Immutable storage and bounded serving

A snapshot directory contains `catalog.sqlite3`, `shards/*.sqlite3.gz`,
`receipts/*.json` and `manifest.json`. Keep these files together. The registry
indexes all accessions and their shard locations. Stable numeric row identifiers
are scoped to a snapshot and preserve every foreign-key relation. The server
verifies compressed checksums, rejects unsafe paths, limits decompression and
SQLite work, and caches a bounded amount of shard data. Anonymous requests cannot
fetch, import or update a catalogue. The job queue remains separate.

The API contract and limits are in [backend/CATALOG.md](../backend/CATALOG.md).
Snapshots under `.catalog/` are ignored by source control. Only `frontend/` is
exported publicly; never copy reference receipts, SSH configuration or the
private registry into the static frontend deployment.

Global scientific searches use the adjacent immutable `exact-search.sqlite3`.
It maps normalized exact terms to candidate shards and is pinned to the registry
and source shard digests. The API verifies actual source values before returning
matches and visits at most four candidate shards per page. Hashes never stand in
for biological evidence. A missing or mismatched index produces an explicit
unavailable response for global scientific queries; exact genome context remains
available. Search values are limited to 1,000 characters, and the interface
rejects oversized input without truncating it.

Build the sidecar offline after publishing a snapshot:

```bash
python3 scripts/index_catalog.py /absolute/snapshot --workers 4
```

The API discovers the adjacent index automatically. An operator can instead set
`CRISPR_WEB_CATALOG_SEARCH_INDEX` to an exact absolute path. See
[index construction and verification](catalog-search-index.md). Keep this index
with the snapshot when configuring the web server.

## Import

The scripts use Python's standard library and are independent of predecessor
project directories. They accept an operator's SSH configuration as an explicit
argument; no SSH keys or host settings are embedded in application code.

1. Run `scripts/catalog_source_inventory.py --node N --gzip` on each source host
   through its configured SSH connection, saving the gzip stdout locally. The
   script only reads the original and recovered published result trees. It
   verifies collection hashes and the scientific profile, and identifies raw
   retained results and unavailable receipts.
2. Reconcile outcomes and create the transfer plan:

   ```bash
   python3 scripts/catalog_source_inventory.py --merge /path/to/vm*.json.gz \
     --expect-total 3331847 --gzip > /path/to/selected.json.gz
   python3 scripts/plan_catalog_transfer.py /path/to/selected.json.gz \
     /path/to/transfer-plan.json.gz
   ```

3. Fetch and normalize into a new private directory:

   ```bash
   python3 scripts/import_catalog.py \
     --output .catalog/snapshots/crisprcasdb-YYYYMMDD \
     --snapshot-id crisprcasdb-YYYYMMDD \
     --plan /path/to/transfer-plan.json.gz \
     --ssh-config /path/to/operator-ssh-config \
     --expected-assemblies 3331847
   ```

The configured host aliases are recorded only in `.ops-private/`. Only bounded
read/transform processes run there; no production services or scientific results
are modified. A full raw copy is unnecessary: verified published SQLite exports
are normalized in memory and compressed before transfer. The local importer
checks each transferred checksum, commits one shard at a time, reserves 5 GiB
free disk space, and refuses to overwrite a published snapshot. Rerun with the
same plan/output to resume an interrupted import. The finished registry is
published atomically only after its expected coverage and SQLite integrity pass.

For an explicitly partial local snapshot, provide `--bundle /path/to/bundle.tar.gz`
and `--sha256 DIGEST` instead of the plan/SSH arguments. Its metadata and UI
identify its imported counts separately from the full campaign.

## Use with the server

Set the exact absolute registry path in the API environment:

```bash
export CRISPR_WEB_CATALOG_DATABASE=/absolute/snapshot/catalog.sqlite3
```

The release launcher also accepts:

```bash
python3 scripts/serve.py --catalog /absolute/snapshot/catalog.sqlite3
```

The chosen backend release must include the catalogue feature; follow the normal
release build/gate process before activating it. Installing/importing reference
data does not activate releases or start production services. Without a configured
catalogue, the API and UI show an explicit unavailable state and existing job
workflows continue to operate.

Verify `/api/v1/catalog`, genome search, all feature types and a genome-scoped
array detail before selecting a new snapshot. Keep a previous immutable snapshot
for rollback. Restart the API with the new explicit path after verifying it.
