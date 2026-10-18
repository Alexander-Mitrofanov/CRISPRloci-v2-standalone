# CrisprCasDb catalogue: council decisions

Recorded 27 September 2026. The data/provenance, API/storage, and frontend
agents reviewed the source campaign and this web server independently, then
agreed on the representation below. These are decisions for the catalogue
integration; the import manifest records what has actually been fetched.

## Source and coverage

Use the completed **ClusterCRISPRCasDb exhaustive campaign**, not its earlier
three-caller benchmark or the incomplete 11 September reporting snapshot.
The frozen 2 September NCBI GenBank bacterial and archaeal catalogue contains
**3,331,847 assemblies**. The final 22 September receipt reports **3,331,684
successful, 151 failed, and 12 unavailable**. Retries do not increase this
denominator. A failed or unavailable assembly is never a negative biological
observation.

The pinned scientific profile is
`5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540`:

- CRISPRidentify **2.3.2**, accepted-only `region-screen`, with array strand
  inference and embedded Cas analysis disabled.
- Standalone CasAndra **0.1.0.dev0+crispridentify.3**, run independently on each
  complete assembly multi-FASTA. Its ORF selection is assembly-wide.
- Workflow **0.6.2**. This is not a three-caller consensus collection.

Source catalogue ID: `c26a75991f66d0c8`; frozen catalogue manifest SHA-256:
`f434c57b2279f7868e282eb657f5df1a15c7104c98e292086d08fe9822438d99`.

Authoritative campaign records consulted in the CrisprCasDb project:

- `recovery/2026-09-22/retry-tail/status/20260922T095116Z/summary.json`
  and `accession-status.json`: final outcomes.
- `recovery/2026-09-20/result-locations.json`: original, recovery,
  redistributed, rebalanced, tail, and retained-source locations.
- `src/crispr_ensemble/exporter.py`: publication SQLite schema.
- `README.md`: coordinates, scientific profiles, publication and compaction
  guarantees. Some older README sections still describe 2.3.1; use each
  publication's pinned analysis profile for the imported release.

The independently completed protein export's inventory was also inspected.
It records the same **3,331,684 unique successful assemblies**, **220,117,236
predicted Cas genes**, and eight byte-identical repeated successful records
removed during reconciliation. These counts describe predictions, not a
curated set of experimentally established Cas proteins.

## Scientific representation

| Entity | Representation and interpretation |
| --- | --- |
| Assembly | Versioned NCBI accession, original terminal outcome, counts, and source publication identity. |
| Sequence | Contig/replicon ID scoped to its assembly; source length and SHA-256. |
| Array | One caller prediction with source call ID, coordinates, category, original certainty score, repeat consensus, and counts. |
| Repeat | An ordered occurrence within an array, preserving its observed sequence and coordinates. A repeat consensus is a separate array attribute. |
| Spacer | An ordered occurrence within an array. Preserve `kind=observed` versus `kind=deletion`; a deletion's boundary is not an observed sequence. |
| Cas system | One independently predicted cassette/system with caller type and subtype, including unresolved classification. |
| Cas gene | One protein-coding prediction, with original profile, score and score kind, strand, translation table, and optional system membership. |

Source coordinates are **1-based, inclusive, on the input forward strand**.
Orientation is separate. Preserve missing coordinates and unknown strand.
Zero-width deletion boundaries use the source's separate **0-based boundary**
convention; do not display them as ordinary intervals or count them as observed
spacers. Array and repeat/spacer sequences retain the source orientation; do
not silently reverse-complement or collapse reverse-complement sequences.

IDs must be stable and collision-free across shards: source publication
identity plus source table primary key is suitable. A caller's `call_id` alone
is not globally unique. Sequence grouping, when offered, must distinguish
unique sequence counts from occurrence counts and retain links to every source
occurrence.

Keep missing values as `null`. In the inspected real publication, array strand
was null as expected; Cas system normalized confidence, confidence kind, and
experimental fields were also null. Native CasAndra classification includes
vote fractions and sometimes further final-classification evidence. These
must not be promoted into a generic confidence percentage or substituted for
missing normalized values. Original profile scores are not calibrated
probabilities.

Cas genes can occur outside accepted Cas systems. Preserve nullable system
membership. No source array-to-system foreign key establishes a biological
association: co-location or proximity can be displayed as context with an
explicit distance, not asserted as proven linkage.

## Storage, fetching, and serving

Use a **read-only, versioned catalogue independent of the job queue**. Only an
explicit administrator import fetches data. Requests cannot update scientific
tools, initiate imports, or contact the source cluster. The catalogue and
source connection details stay outside the public frontend export.

The selected local representation is compressed, immutable per-publication
SQLite shards plus a small accession registry. A bounded decompression cache
supports indexed access to one shard at a time. An inspected 100-assembly
source SQLite was 12.4 MB; retaining normalized scalar fields while excluding
bulky JSON yielded approximately 2.0 MB, or 0.49 MB gzip. This demonstrates
feasibility, not a final whole-catalogue size guarantee. The importer must
measure actual size, reserve free disk space, and report complete coverage
only after verifying its final manifest.

Retain normalized sequences, classifications and scores; preserve source
checksums and profile identity to trace every record. Store contig metadata
only when the contig has at least one array, Cas system, or Cas gene. This
retains the complete coordinate context of all requested features while
avoiding millions of feature-free contigs. Assembly `sequence_count` remains
the original source contig count, not the number of retained contig records.
An assembly with no reported features remains in the accession registry.
Full native JSON,
raw diagnostic paths, and protein sequences are separate large artifacts.
If they are omitted from the lean store, state that explicitly rather than
providing reconstructed or invented protein/whole-genome downloads.

The standalone [source inventory script](../scripts/catalog_source_inventory.py)
finds committed publications, raw terminal records from unfinished batches,
the retained VM3 subset, and unavailable-input job receipts. It verifies
collection-index hashes, the scientific profile, safe source paths, and
declared bundle sizes. The fetcher must additionally verify the full bundle
hash, export artifact hashes, SQLite integrity, and relationships before
publishing a locally derived shard.

Merge inventory outcomes by versioned accession. A validated successful
outcome takes precedence over failed attempts; conflicting successful
checksums fail the import. For remaining terminal outcomes, recovery stage
and recorded time determine the latest result. The full inventory must
reconcile to the final campaign totals; an imported subset must be labeled
as a subset, with its actual genome/entity counts shown separately.

The eight-machine read-only inventory performed for this integration reconciled
exactly to **3,331,847 unique outcomes** with the recorded final status counts.
The raw-outcome adapter was also exercised against the real retained VM3
results and unavailable receipts. Inventory completeness is distinct from
download completeness; only the final local manifest certifies the latter.

Queries are parameterized and bounded. The main browse unit is the assembly:
accession lookup, paginated records, then sequence/array/system detail.
Record tables should show caller evidence without implying validation.
Global filters must be backed by the registry or another explicit index;
an anonymous search must never decompress and scan the whole corpus.

## Later false-positive investigation

The remote `casandra-fp-full-20260923` directory was inspected read-only.
The inspected VM2 worker contains a one-protein Gins15 validation example,
context evidence, and worker provenance. Its `suspected_false_positive`
label still lists unresolved original-profile, provenance, and context gates.
Its manifest explicitly calls this discovery evidence, not a false-positive
rate estimate, and says this stage assigns no negative label.

This directory is **not evidence of a corrected complete catalogue**. Preserve
the original campaign predictions and their uncertainty. A later verified
adjudication release could be added as a separately versioned annotation
layer with per-call provenance; do not silently remove original predictions
or extrapolate the small validation example to the full collection.


## Search-only revision

The subsequent interface revision opens with only a search type selector, input
and submit button. It fetches no records until a complete accession, exact
repeat/spacer sequence, or Cas gene/profile is submitted. Editing a query clears
its output and cancels pending reads. Genome context tabs appear only after a
specific accession is selected; dataset totals and provenance are collapsed
beneath submitted results.

An offline, snapshot-pinned exact-search index maps normalized terms to candidate
shards. The API verifies actual source values in at most four candidate shards
per page. Repeat/spacer strings retain literal IUPAC semantics. Cas queries also
match complete gene-name tokens explicitly present in source annotations or
profiles, with case ignored; this does not establish family membership by
similarity. Gene call identifiers are not a global search key. The source
registry and scientific shard files remain immutable.
