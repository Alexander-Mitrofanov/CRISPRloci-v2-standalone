# Published-substrate comparison

This check evaluates two exact experimental RNA substrates from
[Reimann et al. (2017)](https://doi.org/10.1093/nar/gkw786). It demonstrates a
context-dependent contrast, parameter sensitivity, and an experimental
counterexample. It does not establish biological prediction accuracy or exactly
reproduce the paper's numerical scores.

Run from the CRISPRrepeat checkout after installing the test dependencies:

```sh
.venv/bin/pip install '.[test]'
.venv/bin/pytest -q tests/test_published_reimann.py
.venv/bin/python benchmarks/reimann_2017.py --output reimann-results.json
```

The script uses the installed, pinned ViennaRNA 2.7.0 API without downloading
anything. Omitting `--output` prints JSON. Its default fixture path is relative
to the script, so the command also works from another working directory. The
[recorded run](benchmark-results/reimann-viennarna-2.7.0.json) includes source,
fixture, script, sequence and energy-parameter checksums, effective settings,
experimental labels, all motif pairs, current scores and elapsed time. Runtime
and timestamp fields naturally vary between runs; numeric assertions allow a
small floating-point tolerance. The short elapsed time is not a scalability
benchmark.

## Source binding

The [fixture](../tests/fixtures/reimann_2017.json) contains RNA copied from
Supplemental Table S2 and cleavage labels copied from Supplemental Table S3 of
`gkw786_supp.docx`, retrieved from the
[Europe PMC supplementary archive](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5741207/supplementaryFiles).
The downloaded document's SHA256 is
`707a800dbfd82c357dce3208ab9ff61ef54e0896daae108cf896a63bfc15ca7e`.
The original document is not bundled. Whitespace and cleavage slashes were
removed from sequences; the synthetic initial `GG` was retained. These are
complete experimental transcripts, not automatically inferred genomic array
boundaries. No genomic coordinates are invented.

| Transcript | Repeat | Transcript interval, zero-based half-open | Cas6-1 cleavage | Cas6-2a cleavage | Published model score |
| --- | --- | --- | --- | --- | --- |
| CRISPR2_VIII_190nt | R6 | [38,75) | No | No | 0.004 |
| CRISPR2_VIII_190nt | R7 | [114,151) | Yes | Yes | 0.9773 |
| CRISPR2_I_121nt | R3 | [42,79) | Yes | Yes | 0.0001 |

All three repeat instances have the same 37-nt sequence:

```text
GUUCAACACCCUCUUUUCCCCGUCAGGGGACUGAAAC
```

Sequence lengths, hashes, repeat intervals and retained prefixes are checked
before folding. Experimental cleavage labels remain separate from both the
paper's model scores and the current calculations.

## Current model and observed results

The current implementation infers the following isolated MFE hypothesis:

```text
................(((((....))))).......
```

Its zero-based pair set is `(16,29), (17,28), (18,27), (19,26), (20,25)`.
This is not a recovered curated motif from the paper. Scores below are the mean
support for these five pairs in the respective context model. They do not give
the probability that the complete motif forms or that cleavage occurs.

Three configurations keep the comparison explicit:

- `default_local`: W150/L80, lonely pairs allowed, the tool's default.
- `local_no_lonely_pairs`: W150/L80, `noLP=true`.
- `whole_fragment_no_lonely_pairs`: W equals the exact fragment length,
  L80, `noLP=true`, matching the paper's stated window/span/lonely-pair policy.

All use the tool's declared Turner 2004 parameters, 37°C, dangles=2 and salt
setting. The historical ViennaRNA version and motif definition differ, so the
last preset is an experimental-fragment comparison rather than an exact rerun
of the publication. The 121-nt case uses an effective window of 121 nt under
every preset; its requested default window of 150 remains recorded.

| Repeat | Default local support | Local support with noLP | Whole-fragment support with noLP |
| --- | ---: | ---: | ---: |
| R6 | 0.560092 | 0.560064 | 0.004978 |
| R7 | 0.947893 | 0.951908 | 0.944794 |
| R3 | 0.000213 | 0.000203 | 0.000203 |

The isolated support is 0.959673 with lonely pairs allowed and 0.961181 with
`noLP=true`. R6 has lower context support than R7 under all three configurations.
Holding the lonely-pair policy fixed exposes a large effect of the window on
R6. This makes parameter sensitivity visible instead of treating one score as
an intrinsic property of a repeat.

R3 was experimentally cleaved despite very low support in both the published
and current calculations. The benchmark retains this challenge case; it does
not turn the successful R6/R7 contrast into a universal processing rule.

## What the tests establish

Tests verify source binding, native-engine execution, the R6/R7 contrast,
window sensitivity with other settings held fixed, preservation of the R3
counterexample and finite machine-readable output. Numeric regression anchors
apply only to these exact sequences and the pinned model; they are not
biological thresholds. No classifier, calibrated probability, accuracy
percentage or processing decision is derived.

These related substrates come from one study and were selected because their
behavior was already known. Generalization still requires independent
experimental datasets, locus/family/study separation, comparisons with isolated
folding alone, and evaluation of unsupported processing mechanisms.
