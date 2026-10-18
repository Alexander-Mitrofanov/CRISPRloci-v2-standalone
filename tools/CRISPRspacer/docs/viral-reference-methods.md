# Viral RefSeq search: methods and reproducible provisioning

CRISPRspacer's `viral-search` mode searches supplied spacer sequences against an
operator-provisioned viral nucleotide reference. Its results are candidate
**viral reference sequence matches**. A match alone does not establish infection,
host range, active immunity, or a functional CRISPR target. This workflow is
separate from `database-search`, which searches supplied target DNA against the
2017 NCBI-hosted spacerome collection.

## Reference scope

The first reference scope is every `viral.N.N.genomic.fna.gz` member in a pinned
[NCBI viral RefSeq release](https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/).
This includes phage reference records; it is not restricted to viruses of human
or animal hosts. It also is not all viral sequences in GenBank. NCBI describes
this location as the complete viral RefSeq release. RefSeq can contain incomplete
records, and segmented genomes have separate sequence records, so a record count
is not a virus or species count.
[NCBI Virus documentation](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/)
explains the collection and distinguishes it from the larger all-nucleotide
download.

Keep the original accession **including its version**, the original FASTA title,
and the exact source bytes. Do not infer taxonomy, isolation host, genome
completeness, segment grouping, or circular topology from a title. Those values
remain unknown unless provided by separately pinned, validated metadata. NCBI
recommends citing accession.version and the RefSeq release number for a dataset
obtained from the release FTP site.
[RefSeq citation guidance](https://www.ncbi.nlm.nih.gov/books/NBK50679/)

During the initial operator inspection on 2026-09-28, the official release number
was `237`, and the genomic FASTA inventory contained
`viral.1.1.genomic.fna.gz` with official compressed-file MD5
`104124827b73be6bfd87af42596bf4be`. These are observations of that snapshot, not
constants used to choose future input files. The plan enumerates **all** members
and fails if the release catalog and directory listing disagree.

The completed VM build contains **19,625 sequence records and 583,708,944 bases**.
Its [snapshot receipt](viral-refseq-237-snapshot.json) records source/software
hashes, the normalized content identity, complete reference identity and successful
native and search checks. The full reference bundle occupies 1,343,792,830 bytes;
compressed sources and software are stored separately. No reference sequence data
are included in the receipt or this repository. Rebuilt native metadata may have
different dates or executable paths, so compare retained source and normalized
content identities as well as recording the new bundle identity.

## Acquire sources on the VM

Run acquisition only as an explicit operator action on the VM designated to hold
the complete reference. Ordinary searches have no network or update step. Local
development and tests use small synthetic fixtures; do not copy the full viral
collection to the development checkout.

Use an installed, version-pinned CRISPRspacer environment. Choose a new working
directory on the VM, with at least 20 GiB free before acquisition:

```bash
mkdir -p /operator-chosen/reference-work
cd /operator-chosen/reference-work
crisprspacer viral-plan --output viral-plan.json
crisprspacer viral-download --plan viral-plan.json --output viral-sources
```

The path above is a placeholder to replace with the VM's configured reference
storage directory. The two commands do the following:

1. `viral-plan` fetches only bounded metadata: `RELEASE_NUMBER`, the viral release
   listing, and `release-catalog/release{N}.files.installed`. Each catalog/listing
   response is capped at 8 MiB. It requires a nonempty, identical genomic FASTA
   inventory in the listing and catalog, then fetches all three metadata resources
   again to detect a publication in progress. It records their hashes and every
   selected file's fixed official URL and MD5. The catalog contains MD5 and
   filename, not exact byte size; planned source bytes are therefore `null`.
2. `viral-download` validates the plan's fixed NCBI URLs and filenames, retrieves
   the pinned release catalog, and verifies both the catalog hash and complete
   file inventory. It streams each compressed source, checks its official MD5,
   computes SHA256 and actual byte count, and fully reads the gzip stream to check
   integrity. It rechecks the catalog before publication. Redirects and arbitrary
   supplied download URLs are rejected.
3. Publication creates a new directory atomically. Existing plans and source
   directories are never replaced. Failure removes the incomplete directory.
   `download.json` records source hashes, sizes, source release, official URLs,
   catalog identity, and the SHA256 of the preserved original `plan.json`.
   `release.files.installed` preserves the catalog bytes too.

The download bounds are 8 GiB per compressed member, 16 GiB compressed total,
20 GiB decompressed total, a 24-hour acquisition deadline, and 30-second network
socket timeouts. The 20 GiB initial free-space check is a provisioning minimum,
not a prediction of the final BLAST index size. Resource limits fail explicitly;
they do not publish a truncated reference.

NCBI documents release checksums in the `release###.files.installed` catalog.
MD5 is used to check agreement with the publisher; independently computed SHA256
binds the retained source files and downstream reference.
[NCBI release checksum documentation](https://ncbiinsights.ncbi.nlm.nih.gov/2022/11/15/refseq-release-215/amp/)

## Build and verify an immutable BLAST reference

Install a pinned BLAST+ distribution in an operator-owned environment and retain
the distribution checksum and version output. The initial VM provisioning uses
BLAST+ `2.17.0+`. Local tests may use other supported recorded versions; identical
results across different BLAST versions are not assumed.

For a Linux x86-64 VM, the original installation can be reproduced in a new
operator-owned software directory using the pinned official distribution:

```bash
curl --fail --show-error --proto '=https' --tlsv1.2 \
  -o ncbi-blast-2.17.0+-x64-linux.tar.gz \
  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz
printf '%s\n' 'bdec166721de3b55f90a3badc83538e8  ncbi-blast-2.17.0+-x64-linux.tar.gz' \
  > distribution.md5
md5sum -c distribution.md5
sha256sum ncbi-blast-2.17.0+-x64-linux.tar.gz > distribution.sha256
tar -xzf ncbi-blast-2.17.0+-x64-linux.tar.gz
./ncbi-blast-2.17.0+/bin/blastn -version
```

The MD5 above was read from NCBI's adjacent `.md5` file on 2026-09-28. Preserve
both the archive and computed SHA256. A mismatched archive must not be extracted.
Use the matching official distribution for another operating system or CPU;
this archive is specifically Linux x86-64. NCBI's
[BLAST installation documentation](https://www.ncbi.nlm.nih.gov/books/NBK52640/)
describes the supported installation alternatives.

From the same VM working directory, read the release value in `viral-plan.json`
and supply that exact value. The following command illustrates the inspected
release `237`; replace it only when building from a different explicit plan:

```bash
crisprspacer viral-index \
  --sources viral-sources/viral.*.*.genomic.fna.gz \
  --source-manifest viral-sources/download.json \
  --name ncbi-viral-refseq-237 \
  --source-release 237 \
  --source-url https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ \
  --makeblastdb /operator-chosen/blast/bin/makeblastdb \
  --blastdbcmd /operator-chosen/blast/bin/blastdbcmd \
  --output viral-refseq-237

/operator-chosen/blast/bin/blastdbcmd \
  -db viral-refseq-237/db -dbtype nucl -info
```

The additional NCBI `blastdbcheck` test needs its separate taxonomy lookup files.
In BLAST 2.17.0, its metadata check queries taxonomy even when taxid testing is
disabled. Without `taxdb.bti`/`taxdb.btd`, it can fail with `Taxid 9606 not found`
while individual sequence checks pass. This is a checker dependency; it does not
assign taxonomy or hosts to the viral records. The
[NCBI checker source](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blastdbcheck_8cpp_source.html)
and [taxonomy loader](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/objtools/blast/seqdb_reader/seqdbtax.cpp)
define this behavior.

Provision these lookup files on the VM in a separate new directory, preserving
the archive, official MD5 and local SHA256 receipt:

```bash
mkdir -p /operator-chosen/taxonomy-check
cd /operator-chosen/taxonomy-check
curl --fail --show-error --proto '=https' -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz.md5
curl --fail --show-error --proto '=https' -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
md5sum -c taxdb.tar.gz.md5
sha256sum taxdb.tar.gz > taxdb.tar.gz.sha256
tar -xzf taxdb.tar.gz
BLASTDB=/operator-chosen/taxonomy-check \
  NCBI_DONT_USE_NCBIRC=1 NCBI_CONFIG__BLAST__DATA_LOADERS=none BLAST_USAGE_REPORT=false \
  /operator-chosen/blast/bin/blastdbcheck \
  -db /operator-chosen/reference-work/viral-refseq-237/db -dbtype nucl -full
```

The taxonomy URL is mutable too; retain its exact archive for historical reruns.
Keep these files outside the immutable viral reference. Ordinary CRISPRspacer
searches do not use this lookup or infer taxonomy from it.

The builder streams the large reference separately from the bounded interactive
query parser. It verifies each supplied source against `download.json`, preserves
sequence identity and metadata in the reference bundle, and constructs a local
nucleotide BLAST database. It checks every stored sequence by round-tripping it
through `blastdbcmd`; the additional `blastdbcheck -full` command exercises NCBI's
own database integrity checks. See the
[NCBI integrity-check guide](https://www.ncbi.nlm.nih.gov/books/NBK603843/?report=printable)
and its documented
[`-full` implementation](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blastdbcheck_8cpp_source.html).

The published BLAST prefix is `viral-refseq-237/db`; its reference identity is the
SHA256 of the exact `viral-refseq-237/reference.json` bytes. The reference includes
the normalized subject FASTA, subject metadata index, native BLAST files and a
copy of the supplied source receipt. Retain all files together, the compressed
originals, and the acquisition metadata. Pin the reference identity in operator
configuration before searches use it. Build a new directory for every update;
never replace an index beneath a running search.

BLAST database version 5 uses LMDB and can reserve hundreds of gigabytes of virtual
address space during construction; that is not the same as physical RAM use.
Provision the operator environment accordingly without weakening worker process
containment.
[NCBI makeblastdb options](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.makeblastdb_application_opt/?report=objectonly)

## Search method and interpretation

```bash
crisprspacer viral-search \
  --spacers supplied-spacers.fasta \
  --reference viral-refseq-237 \
  --blastn /operator-chosen/blast/bin/blastn \
  --output viral-match-results
```

The search uses the short-nucleotide BLAST task. NCBI describes `blastn-short` as
intended for short nucleotide queries; its word size is 7. The implementation
records the complete search settings and executable identity, with explicit strand,
masking and scoring settings. The permissive retrieval E-value is not a biological
confidence threshold. Short-query search settings differ from those used for long
genomes.
[NCBI task descriptions](https://www.ncbi.nlm.nih.gov/books/NBK569839/?report=classic),
[blastn options](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.blastn_application_options/),
[NCBI short-query guide](https://blast.ncbi.nlm.nih.gov/blast/BLAST_guide.pdf)

The default acceptance rule retains individual, ungapped BLAST alignments that
cover the full supplied spacer and contain at most two substitutions. An exact
profile uses zero substitutions. Coverage is tested on each alignment, not by
combining disconnected local alignments. Coordinates and relative strand are
reported explicitly; original BLAST evidence is retained.
Rejection reasons are nonexclusive: one rejected HSP can contribute to several
reason counts, while the rejected-HSP total counts that HSP once.

Retrieval remains **heuristic**. BLAST can fail to nominate a candidate, or return
a local alignment that clips a terminal mismatch and consequently fails the
full-span filter. Thus the two-substitution profile does not promise to enumerate
every reference window within Hamming distance two. No-match results mean that no
accepted match was found under the recorded method, reference and limits; they do
not establish absence of a viral relationship.

Retain multiple loci and ties rather than selecting one virus from a short match.
Record limits on subjects, alignments, output size and execution time. Hitting a
limit cannot be interpreted as an exhaustive search. NCBI's default subject limit
is 500, with ties resolved by database order, and setting `max_hsps` to one retains
only one alignment per query/subject pair. These defaults must not silently stand
in for an all-matches promise.
[NCBI BLAST output and limit definitions](https://www.ncbi.nlm.nih.gov/books/NBK279684/)

PAM compatibility and biological function are not evaluated. Array overlap is not
established by searching an external viral reference. References are searched as
linear sequence records; matches spanning the stored end/start boundary of a
circular molecule are not represented unless a future explicit topology-aware
method handles them. Ambiguous sequence characters and incomplete reference
coverage introduce additional uncertainty. A full viral RefSeq snapshot remains
a reference sample of known diversity, not a comprehensive virome.

## Reproduce or update

For a rerun, reuse the recorded reference, source files, exact CRISPRspacer and
BLAST versions, spacer input bytes and search parameters. For a new release,
repeat plan/download/build into new destinations, inspect the new manifest and
integrity results, then update the operator's pinned reference explicitly.

The current FTP location is not treated as an immutable archive. A historical
catalog may remain available after the source file at the live URL has changed.
In that situation the pinned old plan fails its checksum check. Reproduction
requires the retained original bytes; the downloader will not silently substitute
the newer file. Accession.version alone also does not freeze descriptive metadata,
so preserve the source receipt and files even when accession lists appear equal.

Neither reference acquisition nor updates occur inside anonymous requests or
queued scientific jobs. Full reference storage stays on the designated VM;
development fixtures remain synthetic and small.
