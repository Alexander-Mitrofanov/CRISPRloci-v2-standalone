# Benchmark log

This log separates **biological development benchmarks**, **runtime measurements**
and **software regression checks**. Results were collected on 2026-09-26–27.
Repository preparation did not rerun the full literature panel.

## Scope and interpretation

The inventory contains **304,447 source records and 42,712 distinct normalized
sequences**. Of these records, 303,646 are computational predictions. The other
801 mix reference collections, tested systems, engineered guide components and
curated regions; they are not 801 independently validated native transcripts.
See the [inventory](docs/LITERATURE_BENCHMARK.md) and
[machine-readable counts](benchmarks/literature/catalog_summary.json).

The Type II native panel contains **79 inspected inputs and 78 uniquely mapped
references**. Ghc2 is unmapped and excluded from localization denominators. These
references overlap model-training material and were inspected during development.
Results measure regression behavior, not held-out generalization, genome-wide
precision, specificity or false-discovery rates.

Localization requires the correct contig/strand, at least 20 overlapping bases,
and overlap covering at least half the shorter interval. A matching search window
does not establish a complete transcript boundary.

## 1. Literature-driven default fixes

| Development stage | Completed inputs | Localized references | Top-1 | Top-5 |
| --- | ---: | ---: | ---: | ---: |
| Initial literature-benchmark run | 54/79 | 45/78 (57.7%) | 35/78 | 42/78 |
| After long-path repair | 78/79 | 61/78 (78.2%) | 48/78 | 57/78 |
| After genomic model rescue | 78/79 | 65/78 (83.3%) | 51/78 | 61/78 |
| Completed default fixes | 79/79 | 68/78 (87.2%) | 53/78 | 63/78 |

The initial row is an earlier v3 development run, **not a fresh benchmark of the
official v2 release**. Final fixes include short FASTA36 filenames, array masking,
removal of the malformed redundant search direction and genomic model rescue.
Default anti-repeat limits remain 70% identity and 60% repeat coverage.
[Stage results](benchmarks/results/genome_summary.json) retain intermediate runs,
including the incomplete 67/78 stage before the last native-search repair.

## 2. Controlled default versus experimental comparison

Stock FASTA36 varied its shuffled-background statistics for identical inputs. The
first exact historical comparison **failed**: Lsp2 had 16 rather than 14 candidates.
Archived Python code reproduced the new result, locating the difference before
interaction batching. This failure remains recorded.

The amended experiment compared archived/current Python detectors using the
**same separately seeded FASTA36 binary with one native thread**. Their default
biological predictions and ranks matched exactly across 79 inputs and 3,395
hypotheses. Dependency-only differences from the original historical run remain
at Lpn (12 to 10 candidates) and Lsp2 (14 to 16); equal totals do not establish
historical identity. See the [reproducibility amendment](docs/COUNCIL_REPRODUCIBILITY_AMENDMENT.md).

| Measurement under the seeded profile | Default | Experimental short rescue |
| --- | ---: | ---: |
| Completed inputs | 79/79 | 79/79 |
| Localized references | 68/78 (87.2%) | 75/78 (96.2%) |
| Localized at rank 1 | 53/78 | 60/78 |
| Localized within top 5 | 63/78 | 70/78 |
| At least 90% reference coverage | 64/78 | 68/78 |
| Candidate hypotheses | 3,395 | 3,406 |
| Median / 95th percentile / maximum per input | 10 / 192 / 464 | 10 / 193 / 464 |
| Short proposals admitted | 0 | 11 of 66 |
| Observed panel wall time | 344.02 s | 347.17 s |

All original detections are retained. New localized references are Cca1, Cga, Cme1,
Ffr, Ghy2, Orh and Psp. Esp2/Bbo lack a truth-localizing short proposal. Cca2 lacks
qualifying independent support; its terminator E-value of 0.0236 exceeds the
unchanged 0.001 experimental limit.

The option remains **disabled by default**. No model was retrained and thresholds
were not relaxed to fit misses. [Final metrics](benchmarks/council_round2/final_metrics.json)
retain the complete measurements, rejection counts and failed original gate.

## 3. Interaction batching performance

All **3,390 frozen interaction maps** preserved their original fields exactly,
using 152 batched calls across the panel. Controlled stage measurements:

| Locus | Candidates | Serial → batched calls | Serial → batched interaction time |
| --- | ---: | ---: | ---: |
| Mse | 464 | 464 → 6 | 44.61 → 29.90 s |
| Ssu | 104 | 104 → 2 | 11.54 → 8.15 s |
| Sth1A | 128 | 128 → 2 | 15.04 → 11.20 s |

Mse interaction time fell by 33.0% in this measurement. Archived/current whole
default panels took 390.53/344.02 seconds, an observed 11.9% reduction in one paired
run. These are local measurements, not universal speedup guarantees.
[Final metrics](benchmarks/council_round2/final_metrics.json) include timings and
interaction-map equality checks.

## 4. Controls and remaining limitations

| Control | Completed inputs | Short proposals | Short admissions |
| --- | ---: | ---: | ---: |
| Genomic shuffle, default | 19/19 | Disabled | 0 |
| Genomic shuffle, experimental | 19/19 | 30 | 0 |
| Conditional dinucleotide shuffle 1 | 79/79 | 126 | 2 |
| Conditional dinucleotide shuffle 2 | 79/79 | 79 | 0 |
| Conditional repeat derangement | 79/79 | 100 | 1 |

Genomic controls have 60 hypotheses in both profiles. Conditional controls retain
inferred arrays and are not end-to-end negative genomes. One deranged-repeat
window overlaps a reference while its seed lies outside it, illustrating the
limits of window localization. These controls cannot give a biological
false-positive rate.

The bundled Type V models cover Type V-K families. None recognized the 14 assessed
Cas12f components; the `V` label does not imply coverage of all Cas12 systems.
See the [scientific council outcome](docs/council/07_council_outcome.md) and
[remaining-miss diagnosis](benchmarks/council_round2/remaining_miss_diagnosis.json).

## 5. Software quality regression round

The later engineering round did **not** rerun the full literature panel or claim
new detection accuracy. It checked preserved behavior and failure handling.

| Check | Result |
| --- | --- |
| Unit suite | 112 before the round; **360 passed** after repairs/new regressions |
| Native scientific integrations | **7 passed** |
| Native before/after comparison | **104 prediction records in five reports identical** |
| Exports | FASTA/GFF3 identical; CSV identical after normalizing temporary input-directory prefixes |
| Experimental decision probe | 1,156 admission decisions plus standard/adapter/search cases match baseline bytes |
| Independent scientific probe | 697 additional comparison cases match baseline bytes |
| Source distribution | Builds; **360 unit tests passed** from the extracted archive |
| Wheel | Isolated imports/entry points pass; **27 resources** and all 14 modules verify |
| Static checks | Ruff, formatting and strict mypy pass |

All three independent reviewers accepted the changes after a dependency-cleanup
error was repaired and rechecked. See the
[engineering outcome](docs/code_quality/08_final_outcome.md),
[validation summary](docs/code_quality/validation/summary.json) and
[native comparison](docs/code_quality/validation/native-comparison.json).

## Evidence and reproduction

[Artifact scope](benchmarks/README.md) distinguishes compact evidence in Git from
larger workspace-only inputs/runs. Historical metadata can contain the original
machine's absolute paths; it is evidence, not portable configuration.

```bash
python -m pip install -e '.[dev]'
pytest -q -m 'not integration'
CRISPRTRACRRNA_TEST_CONFIG="$PWD/config.toml" pytest -q -m integration
```

The final command requires installed scientific runtimes and your configuration.
Full panel scripts under `scripts/` additionally require downloaded literature
data, generated loci, frozen contexts and archived runs; cloning this compact
repository alone does not reproduce those panels.
