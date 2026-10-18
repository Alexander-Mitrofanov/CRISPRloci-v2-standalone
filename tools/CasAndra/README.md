# CasAndra

CasAndra is an offline, CPU-only command-line tool for detecting Cas proteins,
grouping nearby calls into Cas cassettes, and assigning CRISPR-Cas class, type,
and subtype. The frozen deployment model is included in the Python package, so
normal inference does not depend on the benchmark code or research datasets.

This standalone repository packages runtime **0.3.0.dev2** with the selected
**cap15** model, `casandra-full361-v3-cap15-2026-09-20-v1`: 361 core HMMs,
9 inherited competitor HMMs, 26 additional competitor HMMs, and the cassette
architecture classifier. The runtime Python modules and selected model files
are preserved byte-for-byte from the reviewed September 20, 2026 handoff.

## Scope

CasAndra accepts amino-acid FASTA for protein-level prediction and nucleotide
FASTA, GenBank, or FASTA plus GFF3 for genome-level prediction. It reports Cas
proteins and Cas cassettes. CRISPR arrays, spacers, repeats, and leaders are not
predicted.

## Install

Clone the repository and install with Python 3.11 or newer:

```bash
git clone https://github.com/Alexander-Mitrofanov/CasAndra.git
cd CasAndra
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install .
casandra --version
casandra inspect-model
```

The default model occupies about 76 MB. Installation downloads Python
dependencies; inference is local and offline after installation. To use the
original runtime's tested dependency versions, install with
`python -m pip install -c requirements-tested.txt .`. Scikit-learn 1.8.0 is
pinned because the bundled Joblib estimator requires that serializer version.

## Run

Predict Cas genes and cassettes from a genome:

```bash
casandra predict-genome \
  --genome genome.fna \
  --output casandra_output \
  --gene-mode auto \
  --threads 6
```

For short or metagenomic contigs, use `--gene-mode meta`. To reuse existing
annotations, use `--genbank annotations.gb`, or provide `--gff3
annotations.gff3` together with `--genome genome.fna`.

Classify proteins directly:

```bash
casandra predict-proteins \
  --input proteins.faa \
  --output cas_proteins.jsonl \
  --threads 6
```

For a provenance-bound batch annotation suitable for services and retained
results, use the atomic directory command:

```bash
casandra annotate-proteins \
  --input proteins.faa \
  --output protein_annotation \
  --threads 6
```

Every record carries a literal `result` field: its curated Cas family (for
example, `Cas3`, `Cas9`, or `Csx32`) or the exact negative label `no cas`.
System class, type, and subtype remain supplementary. Cap15 preserves baseline
evidence and marks calls withheld by its specificity rules with
`abstained_from_baseline_Cas_call`. A withheld call is an abstention, not a
claim that the protein is biologically non-Cas. The directory contains
`protein_predictions.jsonl`, `run.json`, and `manifest.json`; run metadata binds
the input digest, model and bundle identities, counts, timing, and result
semantics, while the manifest checksums the prediction and run files.

The bundled frozen model is used by default. An alternative compatible bundle
can be selected with `--model PATH`. Only load trusted bundles: the cassette
classifier uses Joblib/Pickle. `python -m casandra` provides the same CLI.

Classify one cassette from an ordered set of supplied proteins:

```bash
casandra classify-cassette \
  --input putative_cas_genes.faa \
  --output cassette_output \
  --threads 6
```

FASTA record order is treated as cassette gene order. The command reports a
classification (or an unresolved result), all per-protein calls, run
provenance, and a SHA-256 manifest. It does not invent or report nucleotide
coordinates for protein-only input.

Genome output includes TSV summaries, JSONL records, GFF3 annotations, an HTML
report, rejected cassette candidates, run metadata, and SHA-256 manifests.
CasAndra refuses to overwrite an existing output directory.

## Validate the checkout

```bash
python -m pip install -e '.[dev]'
python -m pytest -q
python -m pip install build
python -m build
python -m pip install --force-reinstall --no-deps dist/*.whl
python scripts/verify_installation.py
```

Tests and fixtures run without the original research workspace. The installation
verifier checks the actual wheel's runtime/model identities and exercises
protein, cassette, and genome analysis, including metagenomic CDS calling.

## Python API

```python
from pathlib import Path
from casandra import AnalysisOptions, Analyzer, GenomeInput

analyzer = Analyzer.from_bundle()
result = analyzer.analyze_genome(
    GenomeInput(genome_fasta=Path("genome.fna")),
    Path("genome_results"),
    AnalysisOptions(gene_mode="auto", threads=4),
)
print(result.summary)
```

Reuse one `Analyzer` per process to avoid repeatedly loading the model.

## Provenance and limits

See [provenance](provenance/README.md) for exact source identities and packaging
changes. This is pre-release research software. The included deployment refit
reuses source splits and cannot establish independent accuracy on those
records. Architecture confidence is not a calibrated biological probability.
No new biological performance claim is made by this repository extraction.

The model's historical `local_successor_candidate` / `local_candidate_not_live`
metadata is retained. Publishing this repository does not change the running
web deployment. The service interface lives in the separate
[CasAndra web repository](https://github.com/Alexander-Mitrofanov/casandra-web).

License: Apache-2.0.
