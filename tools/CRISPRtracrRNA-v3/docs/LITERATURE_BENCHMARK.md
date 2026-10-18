# Literature inventory and detection benchmark

The literature search (cutoff 2026-09-26) and local benchmark produced **304,447 source records
and 42,712 distinct normalized sequences**. Complete-mode localization on 78 exactly mapped
published Type II reference loci improved from **45/78 (57.7%) to 68/78 (87.2%)**, and completed
runs increased from 54/79 to 79/79. The retained data and executable benchmark distinguish
sequence recognition, genomic localization, ranking and transcript coverage.

This is an inventory of sequences extractable from the named retrieved sources. It is not a
claim to enumerate every tracrRNA ever published or every biological tracrRNA. In particular,
**303,646 records are computational predictions**. The other 801 records include mixed
reference compilations, active systems, engineered guide components and curated family
regions; they are not 801 independently validated native transcripts. Identical sequences can
have multiple source records, and different sequence lengths can represent processing states
or constructs from the same system.

## Sequence inventory

Start with [reference_sequences.tsv](../benchmarks/literature/reference_sequences.tsv) or
[reference_sequences.fasta](../benchmarks/literature/reference_sequences.fasta) for the 801
records outside the purely computational-prediction category. Each TSV row gives the sequence,
repeat where available, source, evidence category, subtype, organism, original table location,
training-overlap annotation and sequence SHA-256. This subset still requires its evidence labels
to be read.

The complete inventory is available as [catalog.tsv](../benchmarks/literature/catalog.tsv),
[catalog.jsonl](../benchmarks/literature/catalog.jsonl) and
[unique_sequences.fasta](../benchmarks/literature/unique_sequences.fasta).
[catalog_summary.json](../benchmarks/literature/catalog_summary.json) contains the machine-readable
counts. Source files were preserved under `benchmarks/literature/raw/`, with
[input checksums](../benchmarks/literature/input_checksums.json).
[sources.json](../benchmarks/literature/sources.json) records publication metadata, retrieval
status, evidence notes and source licenses;
[download_manifest.json](../benchmarks/literature/download_manifest.json) records the frozen
public file URLs and hashes. There are 18 documented source entries and 123 frozen download
files in these manifests. Accession identifiers remain distinct from organism names; where the
source table does not provide an organism, the catalogue leaves that field empty.

| Source and retrieved material | Records | Evidence interpretation |
| --- | ---: | --- |
| [Chylinski et al. 2014](https://doi.org/10.1093/nar/gku241), Supplementary Table S4 | 63 | Published collection mixing prediction and prior validation |
| [Briner et al. 2014](https://doi.org/10.1016/j.molcel.2014.09.019), original Supplemental Figure S2B | 41 | Comparative reference collection of predicted native RNA regions; boundaries not experimentally established for every row |
| [Gasiunas et al. 2020](https://doi.org/10.1038/s41467-020-19344-1), `Cas9_database` | 79 | Biochemically supported systems; native transcript boundaries are not independently established for every row |
| [Faure et al. 2018](https://doi.org/10.1080/15476286.2018.1493331), Supplementary Data 1, `DB` | 215 | Reference compilation with mixed original evidence |
| Same publication, Supplementary Data 1, `Barrangou` | 144 | Reference compilation with mixed original evidence |
| Same publication, Supplementary Data 3 | 3,353 | Computationally predicted tracrRNA-like sequences |
| [Chyou and Brown 2019](https://doi.org/10.1080/15476286.2018.1498281), official Figshare S4 CSV | 521 | Predictions: 275 in the conservative single-prediction subset and 246 alternatives from genomes with multiple predictions |
| [Dooley et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8336202/), public `AllCas9_TracrRNA_Systems.xlsx` | 41,741 | Predicted sequences, not a new set of experimentally validated transcripts |
| [Alexander et al. 2023](https://doi.org/10.1089/crispr.2022.0090), Table 2 | 96 | Active systems tested with engineered sgRNAs |
| [Ziemann et al. 2025](https://doi.org/10.1093/femsml/uqaf020), Dataset 3 | 112 | Predicted Type V-K RNAs, with selected expression validation in the study |
| [Ruffolo et al. 2025](https://doi.org/10.1038/s41586-025-09298-z), public CRISPR-Cas Atlas v1.0 | 258,031 | Predicted annotations retained from 1,246,088 operons |
| [Cas12f comparative study 2026](https://doi.org/10.1038/s41594-026-01788-6), Supplementary Table 1 | 14 | Explicit tracr components of engineered sgRNAs; not assertions of native transcript boundaries |
| [Rfam RF02348](https://rfam.org/family/RF02348), seed alignment | 37 | Curated family regions, including one duplicate normalized sequence |

After source-specific extraction, catalogue normalization only removes whitespace, uppercases
letters and converts U to T. IUPAC ambiguity codes are preserved. No inferred nucleotides,
spacer removal by guesswork or sequence repair
was used. Forty-one rows with missing or non-sequence values were excluded and retained in
[excluded_rows.json](../benchmarks/literature/excluded_rows.json). This includes 16 literal `NA`
cells in Ziemann Dataset 3; that workbook has 128 rows but only 112 extractable tracr sequences.
Faure's representative distal-tail FASTA was not treated as full tracrRNAs.

The [original Briner supplement](https://ars.els-cdn.com/content/image/1-s2.0-S1097276514007515-mmc1.pdf)
contains the 41 RNAs as an alignment in Figure S2B, comprising 38 distinct sequences. Its PDF
text layer misreads printed gap marks as `N`. The
[source-specific extractor](../scripts/extract_briner_figure.py) checks the frozen PDF hash,
extracts both positioned sequence blocks and restores the printed gaps. An independent pass
using the font's stored character encoding confirms every aligned sequence and gap; visual
inspection confirms the figure. Both raw and corrected alignments are retained in the
[extracted rows](../benchmarks/literature/raw/alternate_retrieval/briner2014_figure_sequences.json),
with [verification details](../benchmarks/literature/raw/alternate_retrieval/briner2014_extraction_qc.json).
The published RNA regions are 85–140 nucleotides long and contain no ambiguous residues after
the verified alignment gaps are removed. This extraction does not establish native transcript
boundaries beyond the authors' predictions.

The [official Chyou Figshare deposit](https://doi.org/10.6084/m9.figshare.6807203.v3), version 3,
provides 521 S4 rows with 498 distinct sequences in the explicitly named `TRACR_RNA_SEQ` field.
`N_PRED=1` selects exactly the paper's 275-record conservative collection. The 246 other rows
are multiple-prediction alternatives. The `CHY` flag identifies 55 rows that the source associates
with earlier references or experimental evidence; their sequences are still this predictor's
versions, so the flag is not promoted to validation of their boundaries. Source probabilities,
terminator flags and subset membership are retained. The deposit's license is CC BY 4.0;
[metadata](../benchmarks/literature/raw/alternate_retrieval/chyou_figshare_article.json) and an
[extraction plan](../benchmarks/literature/raw/alternate_retrieval/extraction_plan.json) document
the files and fields. Its S3 file contains paired sequence comparisons, not ordinary wrapped
FASTA, and its tail models were not used to evaluate the same training collection.

The public Atlas was streamed in full; the retained extraction contains the original operon
ID, metadata, subtype annotation, tracr annotation and CRISPR repeats. Its download size,
generation, SHA-256 and extraction are recorded in
[atlas_manifest.json](../benchmarks/literature/raw/atlas_manifest.json). The Atlas data license
is **CC BY-NC 4.0**, separate from this software's license. Other data retain their original
source terms. An operon subtype attached to a predicted tracr-like annotation does not itself
establish a functional tracrRNA for that subtype.

## Literature coverage and gaps

The search followed the original [CRISPRtracrRNA paper](https://doi.org/10.1093/bioinformatics/btac466)
and its training references, early Type II comparative studies, later biochemical system
collections, public prediction resources, Rfam and subtype-specific Cas12 literature. Primary
article XML and supplementary data were retrieved from Europe PMC/PMC where available,
with the authors' public repositories used for the Dooley and Atlas resources. The Dooley
repository was fixed at commit `07c5d28d1c1afdc03e81938beaafd1069ce6fe15`.
Follow-up searches recovered Briner's original supplement from the official Elsevier archive
and Chyou's supplements from the publisher's Figshare deposit, closing both initial retrieval
gaps while retaining their original failed-access history.
A broad Europe PMC metadata query for `tracrRNA` or `trans-activating crRNA` in titles/abstracts,
from 2011 through the search date, retrieved 179 records; the query and metadata are retained in
[broad_search_2026-09-26.json](../benchmarks/literature/raw/broad_search_2026-09-26.json).
Selected data-rich sources and their cited references were inspected in depth. This is not a
claim that all 179 full texts were reviewed.

The following gaps remain explicit:

- [Cas12l](https://doi.org/10.1038/s41422-022-00771-2) supplementary sequences could not be
  retrieved through the available public endpoints.
- The inspected [Cas12b](https://doi.org/10.1038/s41421-018-0069-3),
  [Cas12e](https://doi.org/10.1038/s41467-024-54491-9),
  [Cas12 reprogramming](https://doi.org/10.1038/s41467-024-50243-x) and
  [2021 Cas12f](https://doi.org/10.1038/s41467-021-26469-4) supplements contained primers,
  constructs or fused guide sequences rather than a clean, separately identified native
  tracrRNA table suitable for automatic extraction. These were not silently reclassified as
  isolated native tracrRNAs.

Publication search coverage and sequence accessibility are separate limitations. More papers,
additional processing states and nonstandard guide architectures can extend this inventory.

## Model recognition of isolated sequences

The sequence benchmark runs the actual bundled Infernal covariance models against normalized
published sequences. It collapses exact duplicates within a model pool, searches the oriented
sequence and records a detection only when a reported hit is on the expected strand and spans
at least 20 nucleotides. A hit need not cover the full RNA. This answers whether a model recognizes
a sequence region; it does not test recovery from a genome or reconstruction of transcript ends.

Selected reference collections give the following source-row counts:

| Collection | Recognized / tested | Interpretation |
| --- | ---: | --- |
| Chylinski 2014 | 37 / 63 | Early mixed reference collection; possible training overlap |
| Briner 2014 | 41 / 41 | 38 unique sequences; original Type II model-training source |
| Gasiunas 2020 | 50 / 79 | Source collection used in model construction |
| Faure reference `DB` | 160 / 215 | Mixed reference compilation |
| Faure `Barrangou` | 123 / 144 | Mixed reference compilation |
| Chyou 2019, all predictions | 169 / 521 | 498 unique sequences; includes the conservative subset and multiple-prediction alternatives |
| Alexander 2023 | 29 / 96 | Postdates the 2022 models, but homologous training overlap has not been excluded |
| Ziemann 2025 Type V-K | 112 / 112 | Closely tied to the bundled Type V-K model source; not independent validation |
| Rfam RF02348 seed | 37 / 37 | 36 unique regions; likely overlap with early Type II references |
| Cas12f 2026 tracr components | 0 / 14 | Type V-F is outside the demonstrated coverage of the bundled V-K models |

Within Chyou's collection, 116/275 conservative single predictions and 53/246 alternatives are
recognized. The 55 rows carrying the source's prior-reference/experimental `CHY` flag have
30 recognized predicted sequence versions. These subsets overlap: the flagged rows belong
to the conservative subset and must not be added again to the total.

The larger prediction collections showed agreement with the bundled models for 2,353/3,352
screened Faure rows, 13,056/41,741 Dooley rows and 119,197/256,797 Atlas rows. These are agreement
counts against predictions from other resources, not independently measured biological recall.
The corrected pooled results are:

| Model pool | Distinct sequences screened | Recognized distinct sequences | Shuffled controls with any hit |
| --- | ---: | ---: | ---: |
| Type II plus subtype-unknown exploratory records | 42,320 | 15,124 | 532 |
| Type V annotations, including subtypes outside V-K | 185 | 112 | 0 |

Across the catalogue, 303,212 source records were assigned a model screen and 1,235 records
with other subtype annotations were marked not assessed. Counts of distinct sequences are
deduplicated within each pool and must not be summed as if the pools were disjoint biological
families. All 91 V-K_V1 and 21 V-K_V2 source records were recognized. Screening other Type V
annotations is a scope check, not a claim that the bundled V-K models cover those RNA families.

Canonical source-row outcomes are in
[sequence_detection.tsv](../benchmarks/results/sequence_detection.tsv), with the corrected
aggregation in [sequence_summary.json](../benchmarks/results/sequence_summary.json).
Retained raw model hits are in
[final/sequences](../benchmarks/final/sequences/),
[final/additional_sequences](../benchmarks/final/additional_sequences/),
[final/chyou_sequences](../benchmarks/final/chyou_sequences/) and
[final/briner_sequences](../benchmarks/final/briner_sequences/). The later batches add Rfam,
Cas12f, Chyou and Briner records retrieved after the larger scan began. Sequence hashes join
all four batches and prevent duplicate counting within each model family. The original model
resources and thresholds were not retrained or optimized against these sequences.

An audit caught an initial pooling error: the raw benchmark's prefix test for `II` also selected
`III`, and its prefix test for `V` selected `VI`. The canonical summary rejoins retained hits by
sequence hash using exact major-type grouping and excludes those annotations from the assessed
denominators. The benchmark script was corrected for future runs; the original raw outputs
remain available for audit.

One deterministic mononucleotide shuffle per unique sequence was generated with seed
`20260926`. Shuffle hits are a composition-control observation, not validated biological false
positives. The random stream depends on pool membership and order; reproducing the exact
historical control counts requires the retained `cm-input.fasta` files whose hashes are in
`sequence_summary.json`. A new scan with the completed catalogue creates different shuffled
controls because the pool membership has changed. The positive collections, duplicated source
rows and predictions cannot
establish genome-wide precision or a calibrated false-discovery rate. Unknown and unsupported
subtypes must remain separate from confirmed Type II or Type V-K references.

## Natural-locus benchmark and matching rule

The 79 Gasiunas RNA sequences were matched against the corresponding natural-locus FASTAs in
the authors' Dooley reference collection. Both strands were searched for an exact sequence
match. Seventy-eight reference RNAs mapped uniquely; **Ghc2 did not map**, with a longest exact
matching block of only 10 bases in its supplied locus. This source discrepancy was retained in
[loci.json](../benchmarks/literature/loci.json) and excluded from the fixed 78-locus detection
denominator. The pipeline still attempted all 79 input loci.
Detailed results mark Ghc2's detection as null/not assessed rather than treating it as a
biological failure.

A localized prediction must have the same contig and strand as the reference, overlap at least
20 nucleotides and overlap at least 50% of the shorter interval. Execution failures count as
non-detections within the fixed mapped denominator. Top-1 and top-5 refer to the position of a
matching prediction in the complete output ranking. Reference coverage measures the fraction
of reference bases overlapped by the best matching candidate; 90% coverage is reported separately.
Neither criterion implies accurate transcript endpoints or the absence of false candidate loci.

All retained stages were rescored using the same
[metric](../benchmarks/results/metric.json) and
[scoring script](../scripts/summarize_literature_benchmark.py).
The original exploratory counters used a looser overlap rule and are preserved in the raw
run directories. Use [genome_summary.json](../benchmarks/results/genome_summary.json) and
[genome_results.tsv](../benchmarks/results/genome_results.tsv) for the final comparisons below.

| Pipeline stage | Completed inputs / 79 | Localized / 78 | Top 1 / 78 | Top 5 / 78 | ≥90% reference coverage / 78 | Reported candidates |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Frozen pre-fix complete mode | 54 | 45 | 35 | 42 | 45 | 2,448 |
| Genomic models alone | 79 | 45 | 44 | 45 | 0 | 208 |
| Short FASTA36 input paths | 78 | 61 | 48 | 57 | 61 | 3,404 |
| Add independent genomic model rescue | 78 | 65 | 51 | 61 | 61 | 3,411 |
| Add accepted-array masking | 78 | 67 | 52 | 62 | 63 | 3,439 |
| Final: remove redundant whole-genome query pass | 79 | **68** | **53** | **63** | **64** | 3,395 |

The final figures correspond to **87.2% localization**, **67.9% top-1**, **80.8% top-5** and
**82.1% at least 90% reference coverage**. Candidate counts include all completed input loci and
all retained hypotheses, including alternative orientations and windows. They are not counts
of validated tracrRNA genes, and the unannotated candidates cannot simply be called false
positives. The increase from the failing baseline reflects both restored execution and greater
candidate sensitivity. It is not evidence of an improvement in precision.

This collection overlaps the sources used to train the Type II model, as described in the
[2022 CRISPRtracrRNA methods](https://doi.org/10.1093/bioinformatics/btac466).
The natural-locus experiment is a useful regression benchmark, **not an independent held-out
accuracy estimate**. Even the later Alexander collection needs clade/homology-based exclusion
before it can support a generalization claim.

## Failure diagnosis and implemented fixes

**Long input paths broke the native FASTA36 library-name parser.** A natural reference run
produced a corrupted library filename when an absolute path reached the legacy parser's buffer
limit. The pipeline now invokes FASTA36 inside each isolated record directory using short input
basenames. This restored 24 of the 25 initially failed input runs without editing FASTA36 or any
upstream repository.

**Perfect repeats inside arrays hid weaker anti-repeat matches.** Post-search exclusion could
not recover hits that FASTA36 had never emitted. Accepted array intervals are now replaced by
`N` in a separate search copy before either FASTA36 or BLAST is run. Original sequence,
coordinates and the later overlap exclusion remain intact. This recovered the Ain and Esp1
reference loci in the measured complete-mode sequence of changes.

**Anti-repeat detection had gated otherwise strong model evidence.** Complete mode now scans
the genomic sequence independently and rescues model-supported regions with E-value ≤`1e-5`.
Hits overlapping accepted arrays or fully contained in an existing same-strand candidate are
excluded. Only overlapping or directly adjacent eligible model hits are merged. Rescue candidates carry an explicit
origin and boundary kind; no array, anti-repeat, interaction or terminator evidence is invented.
Model-supported regions may be shorter than a full RNA, so localization and reference coverage
remain separate metrics.

**A redundant FASTA36 direction caused native alignment corruption.** The former extra search
used the entire genome as a query against short repeats and emitted invalid coordinates and a
negative identity value at Kki. Both repeat orientations were already queried in the regular
direction. Removing the extra whole-genome-query pass avoided this native failure without
relaxing parser checks. A diagnostic Smith–Waterman option did not consistently prevent the
failure and was not adopted. The final run completed all 79 inputs.

The fixes leave default anti-repeat identity at 70% and coverage at 60%. No model retraining or
benchmark-specific sequence insertion was performed. Short or divergent pairing regions and
RNA families absent from the covariance models remain biological/model-scope limitations.
The ten exactly mapped references still missed in the final run are **Esp2, Bbo, Cca2, Cca1, Cga,
Cme1, Ffr, Ghy2, Orh and Psp**. Recovering these indiscriminately by lowering thresholds would
need a separate evaluation with biological negatives and clade-separated validation.

The [per-locus failure analysis](../benchmarks/results/failure_analysis.md) and its
[tabular evidence](../benchmarks/results/failure_diagnosis.tsv) distinguish these residual
misses from software execution errors. All ten runs succeeded, all had an accepted array and
all reference intervals lie outside the masks. None had a returned prediction meeting the
localization rule, so these are not merely ranking failures. Both-orientation genomic BLAST
diagnostics also recovered none at the unchanged identity/coverage gates.

| Residual evidence | Loci | Consequence |
| --- | --- | --- |
| No localizing genomic covariance-model hit reported | Esp2, Bbo, Cca2, Cga, Ffr, Ghy2 | Independent model rescue has no qualifying region to add |
| Localizing model hit weaker than the `1e-5` rescue gate | Cca1 (`0.0045`), Cme1 (`0.012`), Orh (`0.0038`), Psp (`0.00077`) | The current rescue gate intentionally excludes these hits |

Identity-passing repeat alignments in truth-assisted isolated-RNA diagnostics cover less than
60% of the repeat at these missed loci. Those diagnostics help explain the coverage failure;
they do not count as genome detections. Four inferred repeats (Bbo, Cca1, Ffr and Psp) also differ
from the published repeat even after considering both orientations. The retained JSON includes
every inspected alignment and distinguishes a missing native hit from a threshold rejection.

Four successful model rescues illustrate why localization and boundaries must be separated:
Lan covers 58/99 reference bases, Rsp 60/97, Cme3 58/93 and Ghe 61/96. These conserved regions
pass localization but fail the 90% coverage criterion. Their explicit `model_supported_region`
boundary kind prevents them from being mistaken for recovered full transcripts.

## Reproduction and retained evidence

The frozen source used before these changes is retained in
[baseline/source.zip](../benchmarks/baseline/source.zip), with
[source checksums](../benchmarks/baseline/source_sha256.json). Each genome run retains the
effective configuration, raw upstream outputs, native command logs, timing, reports or failure
marker, and export checksums. Baseline failures have not been overwritten.

After retrieving the source files and installing the documented scientific runtimes, install
the pinned [benchmark preparation dependencies](../benchmarks/requirements.txt). They include
`openpyxl`, `python-docx` and `xlrd` for tables, `ijson` for the Atlas stream, and `pdfplumber`
with `pypdf` for Briner's figure extraction. They are separate from runtime package dependencies.

```bash
python -m pip install -r benchmarks/requirements.txt
python scripts/fetch_literature_sources.py --verify-only
python scripts/extract_briner_figure.py
python scripts/build_literature_catalog.py
PYTHONPATH=src python scripts/build_literature_loci.py
PYTHONPATH=src python scripts/benchmark_literature.py sequence \
  --config local.toml --output benchmarks/reproduction/sequences
PYTHONPATH=src python scripts/benchmark_literature.py genome \
  --config local.toml --mode complete --workers 2 \
  --output benchmarks/reproduction/complete
python scripts/summarize_literature_benchmark.py --stage reproduction/complete
python scripts/summarize_sequence_benchmark.py --batch reproduction/sequences
```

To obtain missing public source files, run `fetch_literature_sources.py` without `--verify-only`.
The separate `fetch_tracr_atlas.py` streams the 5.27 GB public Atlas and retains only the
tracr-associated extraction; it verifies the complete source hash before replacing that
extraction. These preparation scripts download data explicitly. Normal pipeline analysis does
not download dependencies or literature data.

Use a fresh output directory for every run. The sequence scan used two Infernal threads;
genome runs used two parallel loci with one tool thread each. With no arguments, the genome
scorer regenerates the retained comparison stages. To score another run, pass a path relative
to `benchmarks`, for example `--stage reproduction/complete`; this rewrites the aggregate
`benchmarks/results/genome_*` files for the selected stages. The sequence summarizer accepts
one or more scan paths through `--batch`; without arguments it reads the four retained batches.
The current benchmark runner uses the canonical localization rule. Original exploratory
counters preserved in earlier run directories used a looser rule and should not be substituted
for the uniformly rescored results. See [VALIDATION.md](VALIDATION.md) for software and native-tool
regression checks, the [independent methodology audit](../benchmarks/results/methodology_audit.md)
for reporting corrections identified before the final source additions, and
[IMPROVEMENTS.md](IMPROVEMENTS.md) for the scientific work still needed.
