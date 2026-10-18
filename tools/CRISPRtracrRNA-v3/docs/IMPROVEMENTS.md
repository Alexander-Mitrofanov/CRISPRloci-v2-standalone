# Scientific and engineering improvements

Implemented improvements include isolated tool runs, full multi-record FASTA support, portable
source IDs, validated machine-readable contracts, dependency pinning/update tooling, checked
coordinate projection, correct terminator significance, Type V motif annotation, explicit
candidate hypotheses, readable feature contributions, GFF3/FASTA exports, input/model/executable
hashes, and failure/completion artifacts. No model retraining was performed.

## Literature-driven changes

The [literature inventory and benchmark](LITERATURE_BENCHMARK.md) adds 304,447 source records,
42,712 distinct normalized sequences, original table locations, evidence labels, source hashes
and explicit training-overlap annotations. Most records are computational predictions; the
catalogue preserves that distinction. It includes the original 41-sequence Briner alignment
and all 521 Chyou predictions, with the paper's 275-record conservative subset labeled separately.
Published native-locus references exposed four concrete
problems that have now been addressed:

1. FASTA36 receives short local input filenames inside each isolated working directory, avoiding
   corruption of long absolute library paths by its legacy parser.
2. Accepted CRISPR arrays are masked in a separate search copy before anti-repeat alignment,
   preventing perfect repeat copies from suppressing weaker external matches. Original bases
   and coordinates are preserved.
3. The redundant whole-genome-as-query FASTA36 pass was removed after it emitted malformed
   native alignments at Kki. Both repeat orientations are still queried in the regular direction,
   and malformed output continues to fail explicitly.
4. Complete mode independently scans genomic sequence for strong covariance-model hits. Hits
   with E-value ≤`1e-5` can be retained even when no anti-repeat candidate was generated. They
   have explicit `candidate_origin="covariance_model"` and
   `boundary_kind="model_supported_region"`, with absent anti-repeat/array/interaction/terminator
   evidence represented as null. An existing same-strand candidate or accepted array prevents
   duplicate/overlapping rescue in the defined cases.

With an unchanged 70% identity and 60% repeat-coverage threshold, these changes improved
localization from 45/78 to 68/78 exactly mapped references and execution from 54/79 to 79/79
inputs. Top-1 localization is 53/78. These references overlap the training-source collection;
the result establishes a regression improvement, not independent generalization or precision.
Ten mapped references remain missed, and model-supported regions do not establish full
transcript boundaries. The four bundled Type V models have demonstrated Type V-K coverage;
none recognized the 14 assessed Cas12f components.

## Highest-value next scientific work

1. **Independent biological validation.** Extend the sourced inventory with experimentally
   supported native transcripts and exhaustively annotated negative regions. Partition by
   species/clade and model-training provenance, and predeclare coordinate matching rules before
   assessing a held-out collection. The current positive-only regression set cannot estimate
   genome-wide precision. Later publication dates alone do not exclude homologous training data.
2. **Calibrate ranking.** The current score is an interpretable heuristic. Estimate feature
   scaling and operating thresholds on held-out data. Test whether correlated identity,
   coverage, and identity × coverage features are overweighted. Stratify Type II and Type V.
3. **Resolve transcript boundaries.** Current complete-mode intervals are anti-repeat search
   windows or model-supported regions. RNA-seq/TSS evidence, independent terminator validation
   and structural support can distinguish full transcripts from mature processed tracrRNAs.
   Neither the 150-base extension nor the covariance-model span is a validated endpoint.
4. **Circular replicons.** V3 clips at contig boundaries. Circular wrapping needs explicit
   topology metadata and split-interval output, with source reconstruction tests, before use.
5. **Candidate consolidation.** Keep raw hypotheses and add a separately defined locus layer
   that groups overlapping predictions across repeat orientations and models. Avoid discarding
   informative alternatives before validating a clustering rule.
6. **Search sensitivity comparisons.** FASTA36 remains the default; BLAST is optional. Compare
   their short anti-repeat sensitivity and hard-negative behavior with predeclared thresholds.
   Explore local Cas/array windows as an optional speed setting only after measuring recall loss.
7. **Broader RNA-family models.** Recover missing primary sequence tables and curate additional
   subtype-specific native references before training new covariance models. Use separate
   training, threshold-selection and test partitions. Extending Type V-K coverage to Cas12f or
   other Cas12 architectures is a model-development task, not a change of subtype label.

## Highest-value scaling work

Exact-repeat IntaRNA batching is implemented and validated in council round two below.
Next, cache upstream outputs using full input, model and configuration hashes; implement
checkpoint validation before resume.
For large cohorts, schedule independent records with an explicit global CPU/memory budget.
Current execution is sequential across records and bounded by the configured tool thread count.
Stream completed records into final exports for very large result sets instead of retaining all
prediction objects until completion. These are future improvements, not implemented guarantees.

## Scope and scientific sources

The preserved baseline is [CRISPRtracrRNA 2.0](https://github.com/BackofenLab/CRISPRtracrRNA/releases/tag/2.0).
The original [CRISPRtracrRNA publication](https://academic.oup.com/bioinformatics/article/38/Supplement_2/ii42/6702007)
describes the combined evidence approach. The supplied model families and their training
provenance should be checked against that work and the
[Type V-K V2 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12416283/) before any new performance claim.
Dependency classifications likewise retain the uncertainty and limitations documented by each
upstream tool. No change here turns classifier scores into biological probabilities.

## Council round two: bounded implementation

The [frozen council plan](COUNCIL_IMPROVEMENT_PLAN.md) separates a default-preserving
engineering change from a disabled experimental detection branch. Exact-repeat
IntaRNA batching retains every target/repeat pair and adds strict all-row structure
and ID validation. The short anti-repeat branch uses native genomic proposals,
actual anchored paired residues and independent regional support; it never reads
reference RNA, coordinates or literature labels at runtime. Candidate-window
model scans for the new branch are isolated from standard candidates, and genomic
CM rescue occurs before experimental insertion, preserving default candidates and
search-context E-values. No model or default biological threshold changed.

The isolated harness is `scripts/council_benchmark.py`; all new artifacts are under
`benchmarks/council_round2`. It records source/input/model/tool hashes, actual wall
time and independent per-locus status. The controls use stable per-item SHA256
seeds. Exact dinucleotide shuffles use randomized Euler trails and preserve
ambiguous positions and run endpoints; they are not claimed to be uniform samples.
The frozen cross-locus donor pool supports a deterministic bijective derangement.
Conditional repeat controls retain inferred arrays and therefore do not estimate
end-to-end biological specificity. Full genomic controls preserve arrays and base
composition while disrupting non-array sequence and Cas context. Neither set is a
validated biological-negative collection.

The initial full default comparison exposed a pre-existing native reproducibility
problem: exact FASTA36 input bytes produced different shuffled-statistics outputs,
and Lsp2 changed from 14 historical candidates to 16. The archived Python detector
also produced 16, confirming that this difference precedes batching. A minimal,
separately built native dependency now supplies explicit RNG seeds; the
[append-only amendment](COUNCIL_REPRODUCIBILITY_AMENDMENT.md) preserves the original
failure and requires a fresh paired old/new detector comparison using identical
pinned dependency bytes. Six Lsp2 and paired Ain/Kki/Mse probes are byte-identical
under the declared one-thread profile.

The amended paired 79-input run preserves every existing default biological
prediction under the same pinned dependency: 3,395 hypotheses and 68/78 mapped
references localized. On the same inspected collection, the experimental branch
localizes 75/78, with top-1/top-5 counts of 60/70 and 68 references covered by at
least 90%. Eleven of 66 short proposals are admitted, giving 3,406 hypotheses.
The additional detections are Cca1, Cga, Cme1, Ffr, Ghy2, Orh and Psp; every
historical detection remains available. These results measure this development
collection and do not establish generalization or calibrated specificity.

Esp2 and Bbo still have no truth-localizing native short proposal. Cca2 has a
truth-localizing, 16-base anchored duplex but no qualifying CM support; its
forward terminator E-value of 0.0236 exceeds the predeclared 0.001 limit. The
threshold is unchanged. Ghc2 remains unmapped and outside the 78-reference
localization denominator. The full remaining-miss evidence is retained in
`benchmarks/council_round2/remaining_miss_diagnosis.json`.

All 3,390 frozen candidate interaction maps match the old fields exactly under
batching. In the controlled Mse comparison, calls fall from 464 to six and
interaction wall time falls from 44.61 to 29.90 seconds (33.0%). The complete
paired archived/current default panels took 390.53/344.02 seconds; this is one
measured panel comparison. Experimental wall time was 347.17 seconds. Inserting
experimental hypotheses shifts 345 rank-derived export IDs while preserving
standard candidate IDs and biological fields.

The paired 19-input genomic controls complete under both profiles with 60
hypotheses each and zero short admissions from 30 experimental proposals. All
three 79-input conditional panels complete: the two dinucleotide controls admit
two and zero candidates, and repeat derangement admits one. Their per-input
medians are zero and maxima are one, zero and one. All control shuffles satisfy
their frozen invariants, with no unchanged repeat or cross-panel identical
repeat control. The prescribed coordinate, export, admission and provenance
audits pass. These controls do not support a biological false-positive rate.

The first repeat shuffle admits an Eit window overlapping the same unannotated
region as the natural experiment's extra Eit candidate. That region remains
unvalidated. The deranged-repeat Cga window overlaps 79 nt of the reference even
though its seed lies outside the reference and its supporting CM overlaps by only
2 nt. This exposes a limitation of broad-window localization. All seven natural
gains were separately checked to have reference-overlapping native seeds; see the
[detection council's control review](council/04_detection_evaluation.md#completed-conditional-repeat-controls).

Final validation passes 119 tests, including seven real-tool cases, plus Ruff and
mypy. The quantitative amended development-panel gates pass; the original
stock-native historical equivalence gate remains failed and explicitly reported.
No default promotion is made. Detailed measurements, gate outcomes and all
retained proposals are in `benchmarks/council_round2/final_metrics.json` and its
paired-results/conditional-control companion files.
The [final council outcome](council/07_council_outcome.md) records the independent
reviews, acceptance scope and remaining limitations.
