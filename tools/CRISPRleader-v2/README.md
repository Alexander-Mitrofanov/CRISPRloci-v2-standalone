# CRisprleader v2

Standalone CRISPR-array context and leader-analysis development project.
**Leader prediction is not available yet.** The biological reference bundle,
inference policy and independent biological evaluation are incomplete. The
working components below are engineering infrastructure; they do not establish
leader presence, functional activity or biological accuracy.

The owner approved the council plan on 27 September 2026. Development uses an
independent repository and three workers for genomic inputs, references/numerics,
and runtime/evaluation. See [current status](docs/STATUS.md) and
[interfaces](docs/CONTRACT.md) and [research basis](docs/RESEARCH.md).

## Install and run

Linux and Python 3.11 or newer. There are no runtime package dependencies.
Atomic output publication requires Linux `renameat2`; subprocess containment
requires POSIX process groups and `waitid` with `WNOWAIT`.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install .
.venv/bin/crisprleader --help
.venv/bin/crisprleader validate --fasta examples/synthetic.fa --arrays examples/synthetic.arrays.json
.venv/bin/crisprleader extract --fasta examples/synthetic.fa --arrays examples/synthetic.arrays.json --flank-length 12 --output context-output
.venv/bin/crisprleader validate --fasta examples/synthetic.fa --result context-output/result.json
```

The example is synthetic software test data. `extract` retains both adjacent
windows without selecting a leader side. Window size is a user policy, not a
predicted biological boundary. The result declares `analysis_mode=context_extraction`
and has an empty predictions list. Existing output directories are never replaced.

## Commands and formats

- `extract`: source-verified flank contexts and JSON, TSV, FASTA, GFF3 and BED.
- `adapt identify`: canonical CRISPRidentify v1 reports to native annotations.
  Requires the original FASTA and verifies the source hashes and coordinates.
  Bona-fide and Possible arrays are included; other categories are recorded as
  excluded. Joined-origin report extensions require native annotations.
- `validate`: source reconstruction and structural validation of annotations or
  results. It does not independently establish biological topology or correctness
  of upstream scientific claims.
- `references verify`: strict local reference manifest/file verification. No
  default pack ships with this snapshot; omitting `--bundle` reports an error.
- `benchmark`: descriptive evaluation of supplied transcription-orientation labels
  and tool calls. Retains missing runs, abstentions and errors; it does not run
  predictors or establish scientific superiority.
- `predict`: reserved interface; exits with a clear execution error until the
  required scientific implementation and evaluation are complete.

Native coordinates and BED are zero-based half-open. GFF3 is one-based inclusive.
Circular arcs use ordered source-forward segments, with a separate strand. GFF3
and BED split wrapped features into linked pieces; FASTA contains each complete
oriented arc. `manifest.json` records output sizes and SHA-256 checksums.
`leaders.tsv` and `motifs.tsv` contain headers only for extraction mode.

Inputs are plain regular-file multi-FASTA and [native annotation JSON](docs/CONTRACT.md).
All FASTA sources are bound by ID, length and the SHA-256 of uppercase sequence.
Circular topology requires an explicit provenance assertion. Reads, array counts,
window lengths and runtime output are bounded. Missing values remain missing;
errors are distinct from completed analyses. The CLI returns 0 for success,
2 for invalid input/unavailable functionality, and 130 for interruption.

The Python API also provides a [generic PWM scanner and reference container](src/crisprleader/data/REFERENCE_SCHEMA.md).
Its synthetic numerical tests establish software behavior only. A statistical
motif match is not a probability of biological function or a leader decision.

## Development and evidence

```sh
.venv/bin/python -m pip install -r requirements-dev.lock
.venv/bin/python -m pip install --no-build-isolation --no-deps -e .
PYTHON=.venv/bin/python scripts/check.sh
```

The [recorded local checks](docs/ENGINEERING-CHECKS.md) include 320 passing tests.
The check runs tests, lint, type checks, a wheel build, and an offline installation
smoke test outside the source tree. CI runs the same checks on Python 3.11–3.13.
The [benchmark protocol](benchmarks/GATES.md) was recorded before biological
prediction outputs. Its gates remain unevaluated, and the candidate partition is
not a curated benchmark. See [benchmark status](benchmarks/README.md).

Original project code is MIT licensed. No legacy CRISPRleader models or third-party
biological tables are redistributed. External assets would retain their own
terms. This repository is not activated in the CRISPRloci server; server integration
and deployment require a separate release decision.

## Scientific comparison work

The [literature-backed comparison protocol](benchmarks/SUPERIORITY_PROTOCOL.md)
and [competitor task audit](benchmarks/COMPARATORS.md) define how a meaningful
comparison must be performed. Conserved-core localization is the proposed primary
scientific task. The implemented reporter addresses transcription orientation
only, as a secondary aid. Margins, independent biological labels, powered sample
size and a final study freeze are still missing.

```sh
.venv/bin/crisprleader benchmark --dataset examples/benchmark.synthetic.dataset.json --predictions examples/benchmark.synthetic.predictions.json --output benchmark-report.json
```

These examples contain invented outcomes for fictional tools, not measured
CRISPRleader or competitor performance. The report always declares
`scientific_superiority_established=false`. See [input/report semantics](docs/EVALUATION_FORMAT.md).
