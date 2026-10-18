# CasAndra

> This directory is the downstream CRISPRidentify snapshot, version
> `0.1.0.dev0+crispridentify.3`. It includes integrity, evidence, Type VII
> gate, protein/cassette hierarchy normalization, packaging, and dependency
> constraints documented in `VENDORING.json`. It is not an official upstream
> release.

CasAndra is an offline, CPU-only command-line tool for detecting Cas proteins,
grouping nearby calls into Cas cassettes, and assigning CRISPR-Cas class, type,
and subtype. The frozen deployment model is included in the Python package, so
normal inference does not depend on the benchmark code or research datasets.

## Scope

CasAndra accepts amino-acid FASTA for protein-level prediction and nucleotide
FASTA, GenBank, or FASTA plus GFF3 for genome-level prediction. It reports Cas
proteins and Cas cassettes. CRISPR arrays, spacers, repeats, and leaders are not
predicted.

## Install

Install CasAndra as part of CRISPRidentify from the CRISPRidentify v2 root:

```bash
./scripts/setup_local.sh
```

Setup creates or reconciles the supported Python 3.12 environment and installs
this snapshot non-editably, with no package-index access and exact validated
dependencies. The model adds about 62 MB. Inference is local and does not
require network access.

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

The bundled frozen model is used by default. An alternative compatible bundle
can be selected with `--model PATH`.

Genome output includes TSV summaries, JSONL records, GFF3 annotations, an HTML
report, rejected cassette candidates, run metadata, and SHA-256 manifests.
CasAndra refuses to overwrite an existing output directory.

## Integrity and provenance

CRISPRidentify setup runs its environment doctor after installation. The doctor
validates the installed version, exact dependencies, model manifests, artifact
hashes, all 341 HMM profiles, and the cassette architecture model.

Exact source provenance and checksums are recorded in `VENDORING.json` and
`SOURCE_MANIFEST.sha256`; downstream changes are summarized in `NOTICE`. The
bundled deployment refit is not an independent evaluation model and must not
be used to re-estimate accuracy on its training records.

License: Apache-2.0.
