# SpacePHARER broad-reference integration decision

The user requested a council review followed by implementation and deployment.
The scientific/reference, architecture, and operations reviews agreed that the
100-genome panel was an integration limit, not a scientific requirement. It
remains useful for regression tests but is inadequate for general discovery.

## Reference scope

Use a dated NCBI Nucleotide snapshot of all records matching the documented
complete-genome and phage/prokaryotic-virus query, subject to explicit sequence
eligibility rules. Include both GenBank and RefSeq; collapse exact sequence
duplicates, prefer the RefSeq representative, and retain accession aliases.
Do not sample an arbitrary number of genomes. Every acquired record must be
accounted for as selected, excluded with a reason, or a duplicate.

The query is a discovery step, not sufficient proof of record eligibility. The
downloaded GenBank taxonomy must start with `Viruses`; nonviral records returned
through search relationships are explicitly excluded. The complete-genome
requirement is a title annotation, not an independent assembly-quality assessment.
An acquisition exceeding 100,000 candidate records fails before download rather
than retaining only the first 100,000. This discovery ceiling applies before
sequence exclusions and deduplication.

The NCBI Datasets isolation-host field cannot alone establish a virus's host:
phages isolated from human samples can have human isolation-host annotations.
Likewise, a title containing “phage” is annotation evidence, not an independently
verified host assignment. The selection combines explicit phage titles with
NCBI's bacterial/archaeal `vhost` annotation and states that coverage depends on
those annotations. The snapshot does not claim to contain every known phage.

DNA sequence ambiguity `N` is retained without guessing bases. Eligible sequences
contain only `ACGTN`, contain at least one called base, and are 1,000–5,000,000 bases
long. RNA, segmented records, multiple source features, partial or
satellite/vector/synthetic title annotations, other ambiguity codes and all-`N`
records are excluded and counted. Each excluded record receives the first
applicable reason in the ordered policy. Input spacer ambiguity remains unsupported.
Target and reversed-ORF
control indexes must be built from the exact same selected collection. The
selection query, accession versions, sequence hashes, metadata, exclusion counts,
download receipts and source database remain available to the operator.

Full-sequence retrieval uses GenBank XML with `gbwithparts`, checked against the
requested nucleotide UIDs and declared lengths. Previously verified `gb` batches
remain valid and retain their original hashes. An explicitly sequence-less
nonviral `CONTIG` record can be accounted for as excluded, with its declared
length retained and its sequence hash missing. A viral record without
its full sequence stops acquisition; missing DNA is never fabricated or silently
treated as a negative search result. Nonempty length mismatches and truncated
responses also fail.

The frozen query results and retained responses make the acquisition reproducible;
they do not constitute an instantaneous snapshot of all NCBI state. Acquisition
start and completion times are recorded. Exact duplicates are collapsed by the
hash of the supplied uppercase sequence, with RefSeq preference and accession
ordering as a deterministic tie-breaker. Reverse-complement equivalents, circular
rotations and near-identical strains are not clustered.

## Prepared reference snapshot

The 2026-10-10 acquisition returned 33,898 candidate records. The selected
reference, `ncbi-complete-dna-phage-20261010`, contains 26,952 distinct genomes
and 1,855,809,866 bases; 5,677 exact duplicate accession records remain mapped
to their representatives. The 1,269 exclusions comprise 598 unsupported sequence
alphabets, 349 satellite/vector/partial/synthetic title annotations, 287 records
failing the single genomic-DNA qualification, 33 nonviral taxonomies, one missing
complete-genome title and one segmented/multiple-source record. These are ordered,
nonoverlapping exclusion reasons.

The built reference has 99 index files totaling approximately 16.8 GiB. Full
operator verification checked all index hashes and every selected genome. Its
manifest SHA-256 is
`645d20f79334bfad1ce6f6f52618bcd40f362fbfcae1fb423e670e0edb9016c5`;
the selection receipt SHA-256 is
`638d60daa9bb0878c18116f9c210626cfabccdcfbcb1a58078cb03325e5d5d22`.
Runtime activation remains a separate operation governed by the criteria below.

## Validation of the prepared snapshot

The candidate release passed isolated API and worker acceptance against this
exact reference under copies of the production resource limits. Reference-derived
positives, an upstream biological support subset, unsupported inputs, a seeded
random-input control, cancellation with descendant cleanup, subsequent job
execution, authenticated exports and API restart recovery all passed.

The six-spacer biological support subset, kept in its two original groups,
produced 34 associations and 56 support alignments with estimated native FDR in
approximately 264 seconds. The reference-derived control produced ten candidates
with unavailable FDR, preserved as candidates rather than accepted associations.
The eight-spacer random control produced no candidates or accepted associations.
The three native-search cases each took approximately 264–271 seconds; the worker
peaked at 3.65 GiB and did not restart. These small-input checks establish
computational operation and preservation of native results, not independent
biological accuracy, FDR validity or performance at the maximum input size.

The software suites, complete runtime release gate and frontend fixture/browser
checks also passed. The corresponding SpacePHARER and frontend revisions were
subsequently published with passing GitHub CI. Production activation and a real
public-site browser submission passed, including reference provenance, native
FDR missingness, desktop/mobile presentation, reload recovery and authenticated
ZIP/TSV downloads. Detailed operational receipts remain in private deployment
evidence.

## Software contract

SpacePHARER 2.1.0 uses a compact version-2 reference manifest and indexed SQLite
genome catalogue. Shared indexes remain immutable and are never copied into
jobs. Only genomes needed to validate reported alignments are read into memory,
with a bounded cache. Index and catalogue hashes are streamed and independently
checked in the tool, Integration and backend. Operator verification additionally
checks every source genome before activation. The native score calculation and
calibration semantics are preserved; a 16-bit source-index capacity defect is
fixed and tested beyond 65,535 sources.

The public result contains the reference digest, count, total bases, collection,
snapshot, scope and selection policy. It does not repeat the complete genome
inventory. Historical version-1 results remain displayable. The 1,000-spacer,
100-group and bounded result/preview limits are separate from reference capacity
and remain in force. Resource exhaustion fails explicitly rather than silently
truncating a scientific result. Missing FDR stays missing and candidates without
calibration are not accepted associations.

## Activation criteria

Build the reference and a new pinned tool/application release separately from
the live runtime. Check exact source/index identities, native preservation tests,
large-catalogue and corruption regressions, the complete release gate, and browser
presentation. Benchmark the actual collection under copied production service
limits, including real native execution, source-bound positive controls, negative
and unsupported controls, grouping, duplicate projection, protected downloads
and recovery. Synthetic controls demonstrate computational operation, not
biological prediction accuracy.

Preserve the previous reference and runtime, drain queued jobs before changing
the configured digest, take a matched queue/data backup, and verify the new
service before reopening submissions. Reference refresh is an explicit operator
action; a submitted job can never download or replace its reference. Production
identifiers and deployment evidence belong only in ignored `.ops-private/`.

## Primary references

- [SpacePHARER database preparation and control databases](https://github.com/soedinglab/spacepharer#creating-databases)
- [NCBI Virus data and metadata documentation](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/)
- [NCBI E-utilities documentation](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
