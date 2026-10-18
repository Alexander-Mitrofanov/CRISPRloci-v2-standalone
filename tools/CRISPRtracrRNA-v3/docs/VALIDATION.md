# Validation evidence

Validation was performed on this Linux x86-64 workstation on 2026-09-26–27. The exact upstream
commits are recorded in `upstream-lock.json`. The locally configured runtimes were:

- Python 3.12.7 for the v3 package;
- Infernal 1.1.4, FASTA36 36.3.8g, IntaRNA 3.2.0, and the preserved ERPIN executable;
- CRISPRidentify v2 checkout `9706d9b2959f74df967f64b41645b11363ace22e`, using its existing
  Python 3.12.13 core environment and its own freshly built native matcher;
- standalone CasAndra checkout `ba046f90242ac91dc075183ea51dcbedfac28484`, runtime 0.3.0.dev2,
  with its packaged cap15 model and the installed compatible scientific dependencies.

## Automated checks

**78 tests passed** (71 unit/contract tests and seven real-tool integration tests) in 20.85 seconds
on the literature-benchmark changes. Ruff lint/format checks for source, tests and scripts,
and strict mypy checking of all 12 package modules passed. The sdist and Linux x86-64 wheel
built successfully.
All 13,913 retained evidence/export hashes across the 79 final natural-locus runs were verified;
the count and result are recorded in
[`evidence_verification.json`](../benchmarks/results/evidence_verification.json).
Machine-readable verification and logs are retained in
[`benchmarks/results/`](../benchmarks/results/), including
[`validation.json`](../benchmarks/results/validation.json) and
[`pytest.log`](../benchmarks/results/pytest.log).

Run the complete test suite using the machine-specific configuration:

```bash
CRISPRTRACRRNA_TEST_CONFIG="$PWD/local.toml" pytest -q tests
ruff check src tests scripts CRISPRtracrRNA.py
ruff format --check src tests scripts CRISPRtracrRNA.py
mypy
python -m build --no-isolation
```

The suite covers FASTA/gzip and identifier handling, IUPAC reverse complements, forward/reverse
coordinate round trips and boundary clipping, malformed configuration, virtual-environment
symlinks, similarity/coverage thresholds, full array overlap exclusion, multi-record Infernal
parsing, ERPIN multi-hit parsing and E-value ordering, named IntaRNA columns, missing-evidence
ranking, Type V motifs, source hashes, additive upstream schema changes, rejection of breaking
schemas, CasAndra checksum validation, shell metacharacters, process-group timeouts, failure
markers, safe CSV/GFF output, and dependency updates against local Git repositories.
The new regression checks also cover short FASTA36 input names under long paths containing
spaces, masking accepted arrays without changing the original sequence, one regular search
with both repeat orientations, model rescue thresholds, forward/reverse rescue coordinates,
exclusion of array overlaps and already covered hits, explicit candidate origins and boundary
kinds, and nullable evidence for candidates supported only by a genomic covariance model.
Nine benchmark-method checks protect exact major-type grouping and the contig, strand,
minimum-overlap and interval-fraction requirements of the localization metric.

Seven opt-in real-tool tests exercise:

1. Type II screening followed by a multi-record forward/reverse-complement reconstruction check.
2. The full CRISPRidentify → CasAndra → FASTA36 → Infernal → ERPIN → IntaRNA pipeline.
3. All four Type V models and the enabled motif analysis.
4. An explicit model-only zero-hit case with a valid header-only CSV.
5. Three natural reference fixtures (Ain, Kki and Lsp1) that reproduce the literature-discovered
   FASTA36 failures. Each requires a correctly stranded anti-repeat window covering at least
   90% of the published reference sequence, with long working paths containing spaces.

The three natural fixtures have source provenance and SHA-256 checksums under
[`tests/data/literature/`](../tests/data/literature/). Their test isolates the native anti-repeat
search and checks source reconstruction; it complements the full pipeline benchmark rather
than repeating every pipeline run in CI.

The real tests are skipped unless `CRISPRTRACRRNA_TEST_CONFIG` is set. Standard CI intentionally
runs unit/contract tests and packaging without access to the private upstream repositories.

## Retained smoke runs

The workspace's ignored `runs/` directory retains the earlier refactor smoke runs and their raw
outputs and logs. These predate the literature-driven masking and rescue changes, so their
candidate counts are historical execution checks, not expected counts for the current defaults:

| Run | Result |
| --- | --- |
| `smoke-model-II` | 442,855 nt upstream Type II example; four merged candidates |
| `smoke-model-V` | 61,319 nt upstream Type V example; one merged candidate |
| `smoke-complete-built` | One accepted array, five Cas calls, one Cas cassette, 100 transcript hypotheses |
| `smoke-blastn` | Optional BLAST backend; 96 transcript hypotheses |
| `wheel-complete` | Installed wheel invoked from `/tmp`; complete pipeline, 100 hypotheses |

For the complete example, 89 hypotheses had IntaRNA interactions, four had forward-oriented
model hits, and ten had forward-oriented terminator hits. These counts refer to transcript
hypotheses, not distinct validated tracrRNAs. The complete smoke run took about 16 seconds on
this machine with two tool threads; that is an observed example, not a general benchmark.

The first fresh-checkout attempt failed because CRISPRidentify's native library had not been
built. It retained `failure.json` and no completion report. Running the upstream native build
script resolved that installation issue without changing tracked upstream files.

The wheel is built with a Linux x86-64 platform tag, installed into a clean local venv without
project-source imports, and exercised from outside the source directory. Model resources and
the ERPIN executable are included; upstream repositories, local configuration and run outputs
are excluded from distributions. The sdist includes tests, examples and dependency tooling.

Every retained evidence/export checksum in these five run manifests was verified.
The original v2 checkout and both upstream dependency checkouts have no tracked or
untracked changes according to Git (upstream-ignored native build products are permitted).

## Literature benchmark

The [literature report](LITERATURE_BENCHMARK.md) documents a separate biological-reference
regression experiment. Seventy-eight of 79 published Gasiunas reference RNA sequences map
exactly and uniquely to their supplied natural loci; Ghc2 remains an explicit source discrepancy.
Using the same contig/strand and overlap rule on every retained stage:

| Result | Frozen complete-mode baseline | Final complete mode |
| --- | ---: | ---: |
| Successful input runs | 54/79 | 79/79 |
| Localized reference loci | 45/78 (57.7%) | 68/78 (87.2%) |
| Reference localized at rank 1 | 35/78 | 53/78 |
| Reference localized within top 5 | 42/78 | 63/78 |
| At least 90% reference coverage | 45/78 | 64/78 |

The fixed mapped denominator includes execution failures as non-detections. The benchmark
retains raw native outputs, baseline failures, final predictions, stage comparisons and a
sequence-level inventory; it is not a test based only on synthetic fixtures. A separate
covariance-model scan covers the retrieved sequence collections and explicitly records
unsupported subtype annotations. See
[`genome_summary.json`](../benchmarks/results/genome_summary.json) and
[`sequence_summary.json`](../benchmarks/results/sequence_summary.json).

The final catalogue contains 304,447 source records and 42,712 distinct normalized sequences.
Of these, 303,646 records are computational predictions; the other 801 retain mixed reference,
system-support, engineered-component or curated-region evidence labels. The model scan covers
303,212 source records; 1,235 records with unsupported subtype annotations remain not assessed.
The Type II/unknown pool recognizes 15,124 of 42,320 distinct sequences, and the exploratory
Type V pool recognizes 112 of 185. These pooled figures measure model agreement across mixed
evidence collections, not independent biological recall.

The recovered Briner figure was checked visually and by independent PDF character decoding;
its 41 rows contain 38 distinct sequences, all recognized by the bundled Type II model.
Chyou's official supplement contributes 521 predictions, of which 169 are recognized; 116 of
the paper's 275-record conservative subset are recognized. These source collections overlap
earlier reference/model-development work. Their inclusion closes two retrieval gaps without
turning predicted intervals into experimentally established boundaries.

## Limits of this evidence

Software tests validate implementation and supported integration contracts. The natural-locus
benchmark demonstrates improved localization on this specific collection, which overlaps the
Type II model's training sources. It does not establish independent biological sensitivity,
specificity, precision, score calibration or experimental transcript boundaries. Ten mapped
reference loci remain undetected; the bundled Type V-K models also miss the 14 assessed Cas12f
tracr components. The older smoke examples originate in v2 and are not independent held-out
data. V3 deliberately changes results where parsing, coordinates, evidence, search or ranking
behavior needed correction; equality to v2's full ranked output is not a validation criterion.
