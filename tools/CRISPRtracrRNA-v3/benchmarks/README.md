# Benchmark evidence in this repository

Start with the [benchmark log](../BENCHMARK_LOG.md) for results and limitations.
The larger development workspace contains approximately 2.5 GB of downloaded
literature, source corpora, generated controls and scientific-tool runs. Git keeps
compact records of those experiments, not the entire workspace.

## Included

- `literature/`: inventory counts, source metadata, retrieval checksums, mapped
  locus metadata, excluded-row records and the 801-record reference subset.
- `results/`: stage summaries, per-locus results, methodology/failure analyses,
  validation records and the small historical source archive.
- `council_round2/`: final metrics, paired comparisons, control/audit summaries,
  proposal/geometry records, source archives and hash manifests; the frozen
  control manifest is also retained.
- [`../tests/data/literature/`](../tests/data/literature/): three small, sourced
  natural regression fixtures used by the opt-in native tests.
- [`../docs/code_quality/`](../docs/code_quality/): engineering plans, independent
  reviews, validation logs, probes and the pre-change source snapshot.

These are original retained records. Absolute paths, stage names and historical
test counts identify the environment and point in development when each record
was produced; they are not portable runtime settings or current status claims.

## Not included

The full `catalog.jsonl`/`catalog.tsv`, complete unique-sequence FASTA, original
downloaded supplements, large sequence-detection table, per-genome native output
trees, control FASTAs, installed Python dependencies and environments remain local.
Two historical source archives that embed workstation configuration also remain
local (`implementation_final.zip` and `pinned_native_python_source.zip`); their
recorded hashes and comparison summaries are retained.
Historical reports may link to those workspace-only paths. Their absence from Git
does not turn aggregate metrics into a complete independently runnable dataset.

The extraction, locus-building, benchmark and audit programs are in `../scripts/`.
Full scientific reproduction requires their original input data and frozen
contexts in addition to the pinned scientific runtimes. Download records retain
source URLs and hashes; upstream access and source-specific licenses still apply.
The standard unit suite is self-contained and does not require those large data.

## Interpreting the evidence

The original exact stock-native history comparison failed. The separately seeded
dependency comparison succeeded under its stated one-thread profile. Both outcomes
remain visible in `council_round2/final_metrics.json`; neither is relabeled by this
repository packaging step. Positive reference panels and synthetic controls do
not establish biological precision, specificity or native transcript endpoints.
