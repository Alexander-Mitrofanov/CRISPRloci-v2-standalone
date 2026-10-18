# Independent methodology audit

Audit date: 2026-09-26. Scope: the literature catalog and locus builders,
sequence/genome benchmark runner, both canonical summarizers, original source
table headers, and retained benchmark results. No scientific analyses were
rerun during this audit. Small read-only consistency checks used retained
sequences, hits, source files, and reports.

## Verified results

The canonical final genome result is unchanged: 68 of 78 uniquely mapped
reference loci localized on the correct contig and strand, 53 with a matching
top-ranked candidate, and 64 with at least 90% reference coverage. All 79
attempted final reports match their current source FASTA file hashes. Every
mapped reference reconstructs exactly from the stated interval and strand.
The unmatched Ghc2 reference is excluded from the fixed localization denominator.

The canonical sequence summaries correctly exclude Type III and VI annotations
from the Type II and V pools. Deduplication is within model family using the
normalized sequence hash. The audited canonical counts were 14,925 recognized
among 41,787 unique Type II/unknown sequences and 112 among 185 unique Type V
sequences. The latter pool includes exploratory screens outside the bundled
V-K model scope. These results measure model recognition, not universal
biological sensitivity or genome-wide precision.

## Findings and corrections

1. The genome runner still used an exploratory overlap rule while the canonical
   scorer required at least 20 nt overlap and at least 50% of the shorter
   interval. The coordinating agent changed the runner to use the canonical
   function and added an optional `--stage` to the scorer. Historical run
   counters remain preserved; the canonical comparison uniformly rescores the
   retained reports.
2. The raw sequence runner's `wrong_strand_only` flag was also set when a
   sequence had positive-strand model hits shorter than 20 nt and no eligible
   positive-strand hit. This affected 684 unique queries in the original Type II
   scan pool: 683 in the corrected Type II/unknown pool and one Type III-D query.
   The 683 canonical sequences correspond to 2,242 catalog records: Faure
   predictions 163, Alexander 2023 3, Dooley predictions 376, and ATLAS
   predictions 1,700. This is a raw diagnostic-label error only; canonical
   detection counts are unchanged. The coordinating agent corrected the flag
   to require absence of positive-strand hits. Historical raw files remain
   unchanged.
3. Dooley's `organism` field contained the source identifier from column C,
   rather than a taxonomic organism name. The coordinating agent cleared that
   field while retaining the identifier as accession. Ziemann's sequence
   locator now names the actual sheet, `Dataset_7(Data_of_CAST_systems)`.
4. The coordinating agent added explicit `--batch` selection and retained input
   hashes to sequence summarization. The control documentation now explains
   that the seeded random stream depends on pool membership and ordering;
   reproducing the exact historical shuffled-control counts requires the
   retained `cm-input.fasta` inputs identified by their hashes. No control
   rescan was needed for this reporting correction.
5. The coordinating agent changed detailed exports for unmatched Ghc2 from
   `detected=false` to null/not assessed. The denominator and aggregate
   localization results are unchanged.
6. `tests/data/literature/README.md` incorrectly cited Gasiunas Supplementary
   Data 3. It now cites Supplementary Data 2, confirmed by the original article
   XML's label for `41467_2020_19344_MOESM5_ESM.xlsx`.

## Interpretation limits

Source-record totals are not unique-sequence totals: Gasiunas contains 79
records and 78 unique sequences; Ziemann contains 112 records and 111 unique
sequences. The duplicate Gasiunas sequence occurs in distinct Ssa and Sth1A
loci, so locus-level scoring retains both. Compilation, predicted, engineered
component, and curated-region evidence remain distinguished in the catalog.
Candidate localization and reference coverage do not establish exact native
transcript boundaries. Positive-only locus truth cannot establish genome-wide
precision. The historical model-training overlap remains a limitation of the
regression benchmark.

Corrections above were implemented by the coordinating agent except for the
fixture citation correction. This audit handoff changed only this document and
that citation; it did not rerun tests or scientific benchmarks.

## Final inventory extension

After this audit, the official Chyou supplement added 521 predicted records and
the original Briner Figure S2B added 41 comparative references. New real-model
scans were joined by sequence hash, preserving first-batch duplicate controls.
The final canonical inventory is 304,447 records / 42,712 unique sequences;
the Type II/unknown pool is 42,320 unique sequences, 15,124 recognized and
532 shuffled controls with any hit. Type V remains 185 / 112 / 0. The
68/78 natural-locus result and its denominator are unchanged. Original audit
counts above describe the earlier inventory, not the expanded final one.
Briner extraction has separate source-specific font/coordinate QC in
`literature/raw/alternate_retrieval/briner2014_extraction_qc.json`.
