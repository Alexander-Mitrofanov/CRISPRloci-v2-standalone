# CRISPRrepeat

Reproducible CRISPR repeat folding and array-context evidence. This standalone
Python CLI uses ViennaRNA 2.7.0 numeric APIs for isolated MFE/partition-function
folding and RNAplfold-equivalent local pairing/accessibility. It accepts repeat
FASTA or original genomic FASTA with canonical CRISPRidentify reports.

```sh
python -m venv .venv
.venv/bin/pip install .
.venv/bin/crisprrepeat --input examples/repeats.fasta --output repeat-results
.venv/bin/crisprrepeat --input genome.fasta --reports identify-results --output array-results
.venv/bin/crisprrepeat --input repeats-rna.fasta --molecule RNA --output rna-results
```

DNA inputs retain both orientation hypotheses by default. RNA is supplied in
transcribed 5′→3′ orientation. `--orientation forward|reverse` is an explicit
user assumption, never chosen using the more stable structure. All genomic
coordinates in this tool are zero-based half-open on the source forward strand.
Detector one-based inclusive intervals are independently reconstructed.

The default model uses Turner 2004 RNA parameters, 37°C, dangles=2, salt=1.021 M,
window=150 nt and maximum span=80 nt (including both endpoint bases). Parameters
and their checksum travel with results. Short contexts explicitly reduce the
effective window/span. `--temperature`, `--window`, `--span`, and
`--no-lonely-pairs` support documented sensitivity analyses.

For every observed repeat/orientation, results contain isolated MFE and ensemble
energy, exact numeric pair values, an isolated MFE pair-set hypothesis, and its
mean pair support in isolation and context. Their difference is descriptive;
different ensembles are compared. Empty motifs and span-ineligible comparisons
remain missing with reasons. Local accessibility uses ViennaRNA's native
single-position calculation. Local pair sums are not pairing probabilities.
RNAplfold does not produce a unique full-array dot-bracket structure.

`repeat-analysis.json` is the complete report. TSV, RNA FASTA and isolated
dot-bracket text are also written. Pair coordinates are zero-based; absent
eligible pairs represent numeric zero because calculations have no reporting
cutoff. Computation or output limits cause an explicit error and no completed
JSON report, never silent context shortening.

Optional `--reference frozen-reference.json` performs full-length unit-edit
alignments to a provenance-bearing local reference, in both orientations,
preserving occurrence records and ties (report cap 50 with total counts).
See [the contract](docs/CONTRACT.md). No reference database or calibrated
strand/subtype/family classifier is shipped: unavailable stages are explicit.
The package neither contacts external servers nor modifies scientific checkouts.

These outputs prioritize structural hypotheses for follow-up. They do not
establish processing, immunity, host identity or array authenticity. Protein
binding, pseudoknots, actual transcription ends and kinetics are not modeled.
Weakly folding/non-hairpin repeats remain valid outputs. Separate instances'
MFE pair sets need not describe the same motif.

## Code navigation

`src/crisprrepeat/` follows the evidence flow: `inputs.py` reconstructs bounded,
source-bound sequences; `folding.py` owns the pinned ViennaRNA model and numeric
validation; `references.py` retrieves frozen-reference alignments; `analysis.py`
composes orientation-specific evidence; `exports.py` publishes bounded outputs;
`cli.py` handles execution. Keep scientific interpretation and missing-value
rules in [the contract](docs/CONTRACT.md) and [scientific scope](docs/SCIENCE.md).
`tests/test_science.py` contains native numerical and coordinate controls;
`test_native_boundaries.py` tests invalid native results and process-global
model contamination; `test_parameter_ownership.py` checks owned Turner2004
tables across global replacement, mixed settings and native callbacks;
`test_reference_bounds.py` compares retrieval to full edit alignments, including
duplicate origins, ties and the report cap. Published substrates live in
`tests/fixtures/` and run through `benchmarks/reimann_2017.py`.

## Verification

```sh
.venv/bin/pip install '.[test]'
.venv/bin/pytest
```

Tests exercise actual ViennaRNA numerics, source tampering, negative-strand
mapping, full-window/global equivalence, context competition and missing-value
semantics. Synthetic examples establish engineering/numeric behavior, not
biological accuracy. See [scientific scope](docs/SCIENCE.md).

## Published experimental substrates

Two attributed Reimann et al. transcripts exercise a context-sensitive R6/R7
contrast and a challenge case that was cleaved despite low structural support:

```sh
.venv/bin/python benchmarks/reimann_2017.py --output reimann-results.json
```

The offline script compares default and whole-fragment settings and records
sequence/source checksums, experimental labels and current model scores in JSON.
The [comparison report](docs/REIMANN_REPRODUCTION.md) documents the substantial
window effect and includes a recorded run. These selected substrates do not
establish biological accuracy, and current inferred-MFE scores do not exactly
reproduce the historical scores. The published-fixture tests run with `pytest`.

Scientific background: [ViennaRNA](https://doi.org/10.1186/1748-7188-6-26),
[CRISPRloci](https://doi.org/10.1093/nar/gkab456),
[context-sensitive Cas6 processing](https://doi.org/10.1093/nar/gkw786), and
[non-stem-loop Cas6 processing](https://doi.org/10.1016/j.str.2016.02.009).
