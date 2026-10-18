# CRISPRmap-v2: council decision

Initial decision: 30 September 2026. Modern extension decision: 1 October 2026.

The user requested an agent council, an implementation, a comparison with the
original CRISPRmap, and a repository under their GitHub account. The coordinator
and three agent reviewers considered scientific relevance, historical reference
and benchmark recovery, and software architecture. This is an internal design
review, not independent experimental validation.

The initial decision produced the historical-reference mapper. The later
2.1.0a1 extension implements a modern occurrence catalogue, indexed search,
operational sequence groups, optional folding and a separate benchmark pipeline.
These scopes are distinguished below; implementation of a pipeline does not
assert that every full-corpus analysis has finished.

## Initial decision: historical-reference mapping

Build a standalone, offline tool for mapping prokaryotic CRISPR repeat sequences
to an explicit reference release. Preserve exact historical annotations, and
provide transparent sequence-family candidates for nonexact sequences. Deliver a
batch command line interface, machine-readable results, bounded input handling,
and a reproducible comparison with recoverable original evidence.

The useful contribution is an inspectable, maintainable reference analysis that
can run locally in contemporary workflows. This release does not claim to
recreate every original clustering analysis or update the underlying biological
classification to cover all diversity known in 2026.

## Historical identity

CRISPRmap already had a version 2 in 2014. The official documentation distinguishes
the 2013 release with 3,527 consensus repeats, 40 sequence families and 33
structural motifs from the 2014 release with 4,719 repeats, 24 families and 18
motifs. `CRISPRmap-v2` is the name of this 2026 software project; it is not a claim
to be the Freiburg 2014 release. Reference names and historical labels must
identify their actual source release. The name also does not refer to the
unrelated optical-screening method named CRISPRmap.

Source: [CRISPRmap official help](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=CRISPRmap).

## Agreed historical-mapper scope

1. Normalize submitted DNA/RNA repeats explicitly, recording the original input
   file's hash and the transformation policy. The report contains normalized
   sequences; keep the source FASTA to reproduce its exact formatting and case.
   Do not invent bases for ambiguous symbols.
2. Search both sequence orientations when appropriate. Use full-length edit
   distance, with an explicit scoring definition and deterministic alignment.
   Preserve the original reference records when identical search sequences are
   deduplicated.
3. Treat an exact match as reference membership evidence. Preserve missing and
   conflicting historical annotations. Compare all equally supported references
   before applying any display cap.
4. For a nonexact match, report similarity evidence and possible historical
   sequence-family associations. A configured maximum distance is an operational
   retrieval filter, not calibrated confidence or a probability of family
   membership. The historical defaults are fixed policies, not thresholds
   optimized on held-out data.
5. Keep motif annotations of reference records distinct from motifs assigned to
   the query. The historical `classify` command does not infer a structural motif
   for a nonexact query, run folding, or reproduce LocARNA structural clustering.
   The modern extension reports folding through a separate command and evidence
   schema; it does not convert those folds into historical motif assignments.
6. Include reference checksums, provenance and method parameters with results.
   Reference acquisition/building is explicit; a normal analysis does not update
   its reference or contact an external service.
7. Provide installation instructions, automated correctness checks and a
   reproducible benchmark report that states exactly which comparisons ran.

## Council discussion and resolution

| Review | Main concern | Resolution |
| --- | --- | --- |
| Scientific relevance | Family, structural motif, subtype and orientation could be conflated. | Report them as separate evidence; no subtype, host, transcription-strand or nonexact motif inference. |
| Historical evidence and benchmarking | Historical releases differ, and executable availability may limit comparison. | Identify every reference/output release. Use authentic original evidence and disclose unavailable comparisons. |
| Engineering | A new predictive model would outrun its training provenance and calibration. | Start with deterministic local reference mapping, explicit ambiguity and bounded machine-readable results. |
| Coordinator | A wrapper alone would add little; unsupported biological claims would undermine the result. | Implement a standalone alignment and evidence workflow, then measure it against the evidence actually recovered. |

The reviewers agreed that exact retrieval alone would be too narrow. Full-length
near-reference search, conflict handling, alignment evidence and a reproducible
benchmark add useful analysis while keeping the interpretation supportable.
They rejected transferring structural motifs from sequence similarity alone and
calling a reference-relative reverse complement a transcription-strand
prediction.

## Benchmark agreement

- Authentic original outputs are suitable for output-concordance tests. They do
  not measure the original executable's runtime or establish a speedup.
- An original executable comparison is valid only when the executable,
  parameters, database release and workload are identified. A newly written
  baseline is not the original tool.
- Historical labels can measure reference recovery or agreement with a historical
  partition. They cannot establish biological accuracy on current diversity.
- Sequence-disjoint held-out evaluations must keep duplicate sequences and
  reverse complements together. An accession split can retain sequence overlap;
  report that coverage separately from an exact/RC-purged evaluation. Exact
  sequence separation is weaker than homology, species or family separation;
  the report must state which separation it uses.
- The fixed candidate policy is evaluated as specified. No claim of calibration
  or development-selected thresholds follows from that evaluation.
- Missing historical labels are not negative labels. Synthetic controls test
  defined behavior, not biological detection specificity.

## Gates for stronger biological claims

The modern importer implements accession-level provenance, retained occurrence
records, source hashes and a reproducible build. Curated biological annotations
and data redistribution rights still require their own evidence; the software
license does not supply either.

Predictive family models need transparent training labels, grouped development
and test partitions, performance at accepted coverage, and absent-family
evaluation. Structural motif assignment needs a structure-aware classification
method and independent validation beyond an individual thermodynamic fold.
Strand and Cas-subtype methods each require separate labels, comparators and
overlap audits. Implemented sequence components and folding hypotheses do not
by themselves satisfy these predictive validation requirements. See
[scientific interpretation](SCIENCE.md) for the evidence and limits.

## Modern catalogue extension — 1 October 2026

The user authorized using the completed ClusterCRISPRCasDB campaign and required
that its existing VMs perform data processing and benchmarking. The council
converged on a frozen provenance-preserving catalogue, separate quality tiers,
disk-backed exact deduplication and lossless indexed search. The implemented
grouping method constructs operational connected components from eligible
Bona-fide consensus nodes at fixed 85/90/95% edit similarity.
Possible-only and observed-variant-only nodes can map to every qualifying core
group but cannot merge groups. Chaining is measured rather than concealed.

The optional ViennaRNA workflow supplies separate hypotheses for both
orientations; it does not create motif IDs or select transcription strand.
Benchmarks distinguish coverage after holding out accession stems, full
retained-cohort exact/RC purging including observed variants, and sample-scoped
homology-purged abstention stress. The latter cannot measure family generalization at the same purge/retrieval threshold. Independent modern
family accuracy and original-executable speed remain unavailable.

The implementation and reproducible operator commands are documented in
[modern reference policies](MODERN_REFERENCE.md) and [build instructions](BUILD_MODERN.md).
Executed release totals and scientific results must be backed by the relevant
complete manifests and VM reports. A completed catalogue is separate from
completed indexing, clustering, folding and benchmark runs; this decision
document does not establish completion of those stages.

## Executed release

The full VM build, grouping, folding and benchmarks completed. The
[release report](../benchmarks/modern-20261001/RELEASE_REPORT.md) records the
verified counts, threshold sensitivity, measured coverage and limitations.
The [validation manifest](../benchmarks/modern-20261001/release-validation.json)
binds these results to their data and implementation revisions.
