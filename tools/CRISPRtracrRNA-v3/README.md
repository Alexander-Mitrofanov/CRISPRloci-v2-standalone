# CRISPRtracrRNA v3

**Reproducible tracrRNA candidate prediction with explicit scientific evidence.**

CRISPRtracrRNA combines CRISPR arrays, anti-repeat similarity, RNA–RNA interactions,
terminators, covariance models and Cas-system context. Version **3.0.0** replaces
the v2 implementation with an installable, typed Python package, validated tool
adapters and traceable results.

Predictions are research hypotheses. Scores are not calibrated probabilities,
and reported intervals do not establish full native transcript boundaries.

**Read next:** [Bug-fix log](BUGFIX_LOG.md) · [Benchmark log](BENCHMARK_LOG.md) ·
[V2 migration](docs/MIGRATION.md) · [Dependency management](docs/DEPENDENCIES.md)

## Features

- **Two analysis modes:** covariance-model screening, or complete analysis with
  independently installed CRISPRidentify v2 and standalone CasAndra.
- **Reliable sequence handling:** FASTA, multi-record FASTA, gzip FASTA and input
  directories; full contig identifiers and file provenance are preserved.
- **Explicit evidence:** candidate origins, strand-aware coordinates, raw hits,
  ranking contributions, tool commands and resource checksums.
- **Safe result publication:** isolated working directories, subprocess timeouts,
  refusal to overwrite existing runs, and a completion report written last.
- **Efficient interactions:** bounded IntaRNA batches share an exact oriented
  repeat, with every returned row validated.

The bundled Type II and four **Type V-K** covariance models are retained from
the original project. The Type V option does not cover every Cas12 subtype.

## Installation

Supported scientific runtime: **Linux x86-64, Python 3.11+**. The orchestration
package uses Python's standard library; scientific programs run as external
commands. The bundled historical ERPIN executable is Linux x86-64.

```bash
git clone https://github.com/Alexander-Mitrofanov/CRISPRtracrRNA-v3.git
cd CRISPRtracrRNA-v3
conda env create -f environment.yml
conda activate crisprtracrrna-v3
```

The environment installs this package in editable mode and provides Infernal,
FASTA36, IntaRNA and BLAST. In an existing compatible environment, use
`python -m pip install -e .` and supply those programs yourself.

### Model-only quick start

Model mode uses Infernal and the bundled models, without CRISPRidentify or CasAndra.

```bash
crisprtracrrna doctor --mode model
crisprtracrrna run --mode model --input examples/model-II.fasta --output runs/model-II
crisprtracrrna run --mode model --model-type V --type-v-motifs --input examples/model-V.fasta --output runs/model-V
```

### Complete analysis

Install the upstream applications in their own environments. Exact tested
revisions are in [upstream-lock.json](upstream-lock.json). GitHub access to both
repositories is required; they were private when validation was performed.

```bash
python scripts/dependencies.py sync
cd external/CRISPRidentify-v2
./install.sh
cd ../..
python -m venv .envs/casandra
.envs/casandra/bin/python -m pip install ./external/CasAndra
cp config.example.toml config.toml
```

Set these commands in `config.toml`, adapting them to your installations:

```toml
identify_command = ["external/CRISPRidentify-v2/crispridentify"]
casandra_command = [".envs/casandra/bin/casandra"]
```

```bash
crisprtracrrna doctor --config config.toml
crisprtracrrna run --config config.toml --input examples/complete.fasta --output runs/complete
```

`doctor` checks paths/resources; actual execution validates upstream runtime
contracts. Analysis never downloads dependencies. Commands are argument arrays,
not shell snippets. See [dependency setup and updates](docs/DEPENDENCIES.md) for
environment overrides, native builds and update recovery.

## Inputs and configuration

`--input` accepts a FASTA file, multi-record/gzip FASTA, or directory of FASTAs.
Each record is analyzed independently. Duplicate IDs within one FASTA are rejected;
identical IDs in different files are supported. Use a fresh output path per run.

Defaults retain 70% anti-repeat identity, 60% repeat coverage, array masking and
independent genomic model rescue at E-value ≤ `1e-5`.

| Option | Purpose |
| --- | --- |
| `--mode model` | Covariance-model screening |
| `--model-type II` / `V` | Select Type II or Type V-K families |
| `--search-backend blastn` | Explicit alternative to FASTA36; results may differ |
| `--intarna-batch-size 128` | Bound interaction batches; `1` selects serial comparison |
| `--no-model-rescue` / `--no-mask-arrays` | Disable those stages for controlled comparisons |
| `--short-anti-repeat-rescue` | Experimental Type II complete-mode branch; **off by default** |

The experimental branch requires a short native alignment, an anchored duplex and
independent model or terminator support. Its fixed gates and limitations are in
the [benchmark log](BENCHMARK_LOG.md). It is not promoted to a default setting.

```python
from pathlib import Path
from crisprtracrrna.config import load_settings
from crisprtracrrna.pipeline import run

settings = load_settings(Path("config.toml"), {"threads": 4})
report = run(Path("genome.fasta"), Path("runs/python-api"), settings)
```

Run `crisprtracrrna run --help` for all options. The source launcher
`python CRISPRtracrRNA.py` supports documented legacy aliases; v3 deliberately uses
a new output schema. See [migration details](docs/MIGRATION.md).

## Results and interpretation

| File | Contents |
| --- | --- |
| `report.json` | Canonical versioned result, published only after success |
| `summary.csv` | Quoted table, including a header for zero candidates |
| `predictions.fasta` | Candidate sequences in reported transcript orientation |
| `predictions.gff3` | Features in forward-reference genomic coordinates |
| `configuration.json` | Effective settings; environment values are redacted |
| `manifest.json` | SHA-256 hashes of retained evidence and exports |
| `work/` | Per-record results, upstream reports, raw outputs and command logs |
| `failure.json` | Failure diagnostic when storage permits; no successful report |

Public genomic coordinates are **1-based and end-inclusive on the original
forward strand**, with strand separate. Local evidence coordinates refer to the
oriented candidate. Export IDs map to original file/contig/header/digest metadata.

Complete-mode `anti_repeat` candidates are search windows, including up to five
upstream and 150 downstream bases by default. `covariance_model` rescue candidates
are model-supported regions; unavailable array, interaction and terminator
evidence is null. Neither origin establishes a complete transcript boundary.
Missing Cas context does not suppress search. Unverified predictions are not
automatically false positives.

## Validation and benchmarks

The latest engineering review passed **360 unit tests and seven native-tool
integrations**, strict typing, lint, packaging and isolated wheel checks. It
preserved all **104 prediction records** in five native comparison reports. Three
independent reviewers accepted the changes after an additional failure-path fix.
See [the validation report](docs/code_quality/08_final_outcome.md).

Earlier development-locus experiments localized **68/78 references with defaults**
and **75/78 with the experimental branch**. These inspected references overlap
model-training material. The [benchmark log](BENCHMARK_LOG.md) explains the
counts, controls, timings and original reproducibility failure; these figures are
not independent accuracy or specificity estimates.

Stock FASTA36 showed variable shuffled-background statistics. The controlled
comparison used a separately seeded build with **one native thread**. See
[native reproducibility artifacts](native/README.md). Ordinary configurations
remain supported; arbitrary multithread reproducibility is not claimed.

## Development

```bash
python -m pip install -e '.[dev]'
pytest -q -m 'not integration'
ruff check src tests scripts CRISPRtracrRNA.py
ruff format --check src tests scripts CRISPRtracrRNA.py
mypy
python -m build
python scripts/check_wheel.py dist/*.whl
CRISPRTRACRRNA_TEST_CONFIG="$PWD/config.toml" pytest -q -m integration
```

The last command requires configured scientific runtimes. The wheel check uses a
fresh environment outside the checkout and verifies both entry points, all 27
resources, model indexes, the typing marker and ERPIN executable permissions.
CI runs portable checks on Python 3.11–3.13; native integrations are opt-in.

This repository includes source, tests, examples, documentation, compact benchmark
evidence and optional native reproducibility artifacts. Local configurations,
environments, external checkouts and raw benchmark runs are excluded. See
[benchmark artifact scope](benchmarks/README.md).

## Origin and licensing

V3 derives from BackofenLab/CRISPRtracrRNA **2.0**, commit
`d4cbbbb2172d4114c28fdb8a93c8abae871140f3`. V3 does not require the separate historical
checkout. The original [MIT license](LICENSE) is retained.
[THIRD_PARTY.md](THIRD_PARTY.md) records bundled resources and dependency licensing.
No models were retrained for these changes.
