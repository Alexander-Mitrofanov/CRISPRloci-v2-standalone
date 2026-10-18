# Spacer group associations with SpacePHARER

`spacer_association` consumes DNA spacers as FASTA or a complete detection
`spacers.json` sequence-set. It combines translated and nucleotide search evidence
against an immutable operator-provided target/control reference database. The standalone
[SpacePHARER 2.1.1 repository](https://github.com/Alexander-Mitrofanov/SpacePHARER/tree/4f264252846acef7e0d291209e3e5e3f6faf306c)
was refactored and validated before integration; `tools/tools.lock.json`
pins that revision. It is a fork of upstream commit
`d1e94c04025be6e04fd8204f4f8bb1ea782e993c`, not a new upstream release.

## Explicit grouping and biological interpretation

The API requires `association_grouping`; the interface starts with no selection.

- `single_set`: the user explicitly asserts all submitted spacers belong together.
- `per_record`: each spacer occurrence forms its own group (at most 100).
- `per_source`: a complete sequence-set must supply origins for every occurrence.
  The pair `(source_id, source_sha256)` defines separate groups.

Source records and contigs are not automatically hosts. Supplied origin metadata
is retained as an annotation, never promoted into independent source verification.
The server accepts at most 1000 occurrences, 200 bases each and 100 groups, also
subject to configured request/sequence budgets. Duplicate occurrences retain IDs
and origins; identical sequences within a group are searched once and projected
onto all their occurrences. Group membership and support remain separate, while empirical target/control calibration pools score populations across the entire run. Adding unrelated groups can therefore affect estimated FDR; batch invariance and independent per-group calibration are not promised.
Ambiguous IUPAC sequences and sequences shorter than 27 bases are unsupported,
not negative evidence. A group–target similarity association does not establish
active infection, cleavage, a PAM or a complete host identity.

```json
{
  "mode": "spacer_association",
  "sequence": ">spacer1\nACGTACGTACGTACGTACGTACGTACGTACGT\n",
  "association_grouping": "single_set",
  "tool_options": {"spacepharer": {"--fdr": 0.05}}
}
```

The native combined score is not a probability. FDR is an empirical estimate
from native target/control score populations. A versioned native calibration
receipt supplies `target_score_count`, `control_score_count` and its status.
Genome and control-score counts alone do not guarantee an FDR estimate; the native
target/control score curve must support its cutoff estimate.
`native_fdr` stays null when calibration is unavailable; these rows have status
`uncalibrated_candidate` and are never accepted predictions. Numeric FDR above
the requested cutoff remains `above_cutoff_candidate`; only estimated FDR at or
below the cutoff is `association`. All-unsupported runs are `not_run` with a null
calibration receipt. No missing scientific value is substituted with zero.

## Operator preparation and reference scope

The version-2 reference supports up to 100,000 single-record genomes, 5 Mb per
genome and 10 Gb total, with explicit index and catalogue bounds. These are
engineering ceilings, not a target database size. There is no 100-genome cap.
The [council decision](spacepharer-council.md) explains the database and scaling
choices. Target and reversed-translated-ORF control indexes are separate from
CRISPRspacer BLAST indexes and cannot be interchanged.

Build a dated NCBI snapshot through the operator-only acquisition and selection
script. It freezes the complete result of its documented query, downloads
versioned records with original metadata, verifies every requested record, and
accounts for exclusions and exact-sequence duplicate aliases. GenBank and RefSeq
records are included; RefSeq representatives take precedence for identical
sequences. The collection is all eligible records under the query and policy,
not every phage in nature. Missing annotations, RNA and segmented genomes, partial
records and unsupported sequences limit coverage. Isolation-host metadata is
not treated as independently established infection host range.

The discovery query is:

```text
txid10239[Organism:exp] AND "complete genome"[Title]
AND (phage[Title] OR bacteriophage[Title]
     OR "vhost bacteria"[Filter] OR "vhost archaea"[Filter])
```

The downloaded record must independently have a taxonomy beginning with
`Viruses`. Search relationships can return nonviral records, which are accounted
for as `not_viral_taxonomy`. “Complete genome” is a title annotation; it is not an
independent completeness or assembly-quality measurement. Eligible records have
one unsegmented genomic-DNA source feature, 1,000–5,000,000 bases, only `ACGTN`, and
at least one called base. RNA, segments, multiple source features,
partial/satellite/vector/synthetic title annotations, unsupported ambiguity codes
and all-`N` sequences are excluded. Each excluded accession receives one reason
from the ordered policy; exclusion counts are not independent overlapping totals.

The acquisition requires an exact enumeration of all candidate UIDs and fails if
the query exceeds 100,000 candidates. It does not download a capped subset. This
limit is checked before exclusions and deduplication; the separate reference
capacity limits apply to selected sequences. Query results are frozen before
batch retrieval, and raw-response hashes, accession versions, and acquisition
start/completion times are retained. The result is a reproducible dated
acquisition, not an instantaneous snapshot of all NCBI records.

```sh
python scripts/prepare_phage_reference.py acquire --snapshot /data/ncbi-phage-DATE \
  --download-workers 3
python scripts/prepare_phage_reference.py select --snapshot /data/ncbi-phage-DATE
/path/to/spacepharer prepare-reference \
  --genome-list /data/ncbi-phage-DATE/selected/genomes.txt \
  --source-metadata /data/ncbi-phage-DATE/selected/source.json \
  --scope broad_phage_reference --name ncbi-phage-DATE \
  --threads 8 --build-timeout 14400 --output /data/references/ncbi-phage-DATE
/path/to/spacepharer verify-reference --reference /data/references/ncbi-phage-DATE
sha256sum /data/references/ncbi-phage-DATE/reference.json
```

Downloads use one worker by default or up to three workers sharing a request
budget below two requests per second. Resume verifies completed batch hashes and
UID coverage receipts without repeating the frozen queries. New requests use
`gbwithparts` XML to obtain expanded sequences; existing verified `gb` batches
remain valid. Every accepted batch must contain exactly its requested UIDs.
Truncated XML, sequence-length mismatches and missing viral sequences
fail acquisition. An explicit sequence-less nonviral `CONTIG` can instead be
recorded as excluded, preserving declared length and a missing sequence hash.
No missing sequence is replaced by an empty or fabricated genome.

Keep the selection SQLite metadata, accession/alias list, selection receipt and
source download hashes alongside the sealed reference. The public selection
policy includes the selection receipt digest. `prepare-reference` also accepts
explicit `--genome` arguments for small operator-supplied panels; their scope
remains `targeted_panel`. Database preparation and downloads never run in jobs.

Exact-sequence deduplication preserves all accession aliases and prefers RefSeq,
then accession ordering, for the representative. It does not cluster species,
similar strains, reverse-complement equivalents or circularly rotated genomes.
Neither the broad scope label nor a calibrated association establishes complete
known-phage coverage, infection, or independently verified host range. A missing
match applies only to the installed reference and the submitted spacer evidence.

Set `CRISPR_WEB_ASSOCIATION_REFERENCE` to the sealed reference directory and
`CRISPR_WEB_ASSOCIATION_REFERENCE_SHA256` to its manifest digest. SpacePHARER 2.1.0
must be pinned in the release. A compact manifest binds all index files and an
SQLite genome catalogue. The tool, adapter and independent backend gate stream
complete index/catalogue checksums and validate catalogue identities. Sequences
are loaded only for emitted matches, with a bounded cache; operator preactivation
verification checks every source sequence. Symlinked references or members are
rejected. `N` bases in reference genomes remain unknown, not inferred bases.

Readiness rechecks the compact manifest digest. Each queued job pins that digest;
a reference change requires draining the old queue, building and testing a new
reference/release, and atomically changing the operator configuration. Retain the
old reference and runtime for rollback. Never mutate an active reference.

`provision_tools.py` exports the pinned Git source, installs the locked build
requirements in a separate environment, builds the first-party C++17 engine and
vendored MMseqs with SSE4.1 (AVX2 disabled), runs native CTest, and installs the
native binary beside the Python CLI. Native version/digest, the environment and
all source files including required vendored zstd build recipes are receipted.
Earlier tool receipts retain their inventory semantics. Scientific manifest 8.0.0
requires SpacePHARER and the core tools; unrelated optional tools are not required.
Reference acquisition never activates a production runtime.

## Reports and validation

The worker invokes `associate` with bounded threads, a 2048 MiB native split-memory
budget, job deadlines and descendant containment. The Python wrapper terminates
its native process group on interruption, timeout and leader exit. The durable
worker additionally tracks and terminates descendants, including separate groups.

The adapter and independent backend gate check record order/sequences/hashes,
explicit group partitions and supported/unique occurrence counts, target identity,
calibration state/counts, FDR interpretation, finite score bounds, hit counts and
complete duplicate projection. Both reconstruct query and target alignment slices
from original sequences. Coordinates are 1-based inclusive and can descend for
reverse-complement orientation; gaps are preserved. No flank is labeled a PAM.
The publication gate does not rerun the native scoring algorithm.

Version-2 reports publish compact reference counts, scope and provenance instead
of the entire genome inventory. Historical version-1 result exports remain readable.
The browser previews at most 100 associations and 50 alignments per association,
with truncation stated. Complete `spacer-associations.json`, association-level
`spacer-associations.tsv`, `input-sequences.json`, compact `result.json` and the ZIP
use existing capability-protected downloads. Compact JSON joins all submitted
occurrences and origins to groups/support. Native working reports, indexes, logs
and input copies under `association/` and `association-input/` are private.
See [compact example](examples/module-results/spacer_association.json) and
[module schema](schemas/module-result-v1.json).

## Scientific preservation and verification

The standalone refactor reorganizes first-party native statistics, checked parsers,
RAII readers, workflow execution and build code. It preserves vendored MMseqs
scientific code and upstream GPL/authorship. Oracle tests compare well-defined
scores/support and calibrated filtering with an independently built original
upstream engine. A frozen 300-target/300-control score fixture tests the calibrated
accepted/rejected paths. Upstream biological fixtures supply actual alignments and
reverse-strand coordinates; the small three-phage panel legitimately produces two
uncalibrated candidates with six supporting alignments. It is a test fixture, never
advertised as a production collection.

Documented correctness fixes include short-control undefined behavior and replacing
false zero FDR with an unavailable sentinel and receipt. The upstream order-sensitive
combined-score fallback is preserved and explicitly tested. Tests establish
computational preservation, not new biological performance claims.

Run `SPACEPHARER_EXECUTABLE=/absolute/path/to/spacepharer
REPEAT_TYPER_EXECUTABLE=/absolute/path/to/repeattyper scripts/check.sh` with both
variables on the same shell command line. Real local HTTP tests cover queued
native execution, per-source isolation, duplicate projection, detection handoffs,
unsupported and negative controls, malicious outputs and authenticated downloads.
`frontend/scripts/browser_spacer_association.py` separately tests the UI with
intercepted API responses generated by the actual native HTTP test. Operator
`scripts/gate_release.py` adds the real three-genome positive/unsupported fixtures
to a complete staged runtime gate, without public database downloads.

A full staged toolset build, native release gate and resource-constrained HTTP
acceptance against the actual reference are operator prerequisites before deployment.
Local correctness tests alone do not establish broad-database operational capacity.

Please cite [Zhang et al., Bioinformatics (2021)](https://doi.org/10.1093/bioinformatics/btab222).

## Historical targeted-panel verification recorded 10 October 2026

The following record describes the earlier 2.0.2 targeted-panel integration; it
does not verify the capacity or biological coverage of a broad reference.

The final backend gate passed **775 Python tests, one existing skip**, including
actual native HTTP jobs and protected JSON/TSV/ZIP exports. Ruff passed. The final
frontend gate passed **254 tests**, source boundaries, lint/format checks, the
production build and public-export scan. An earlier full run identified a missing
documentation schema enum and the frontend line-length gate; these were fixed and
the affected contract/frontend checks passed. The generic all-modes form test now
explicitly selects grouping after confirming it starts disabled.

The published standalone v2.0.2 passed 25 Python/native-oracle tests and active
Release CTest assertions; its clean GitHub CI build passed. Parent acceptance
built backend 2.9.0 and Integration 0.11.0 wheels, verified all 116 packaged
Python/JSON source files against the checkout, and executed the actual 2.0.2 native
engine from extracted wheels on 71 occurrences across three groups. This retained
two candidates, six alignments, one unsupported occurrence and duplicate IDs;
independent validation rejected a changed alignment. Browser acceptance used the
actual HTTP job's public evidence/input fixture, checked desktop and 390px widths,
explicit grouping, source metadata, calibration uncertainty, duplicate support,
reload recovery, protected TSV download and absence of JS/CSP errors.
