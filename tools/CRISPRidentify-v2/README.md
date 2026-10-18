# CRISPRidentify v2

CRISPRidentify detects and classifies CRISPR arrays in bacterial and archaeal
DNA. This release combines the CRISPR array workflow with an optional,
bundled CasAndra runtime for offline Cas-gene and CRISPR-Cas cassette
annotation. The current GitHub source enables core-supported reconstruction
by default. The published 2.3.2 archives below predate this source change.

## Quick start from the release archive

CRISPRidentify v2 supports Linux x86-64 and requires Conda, Mamba, or
Micromamba. Installation needs internet access for the locked Conda packages;
normal analysis is local. Distribute these two files together:

- `CRISPRidentify-v2-2.3.2.tar.gz`
- `CRISPRidentify-v2-2.3.2.tar.gz.sha256`

From the directory containing both files, verify and extract the release:

```bash
sha256sum --check CRISPRidentify-v2-2.3.2.tar.gz.sha256
tar -xzf CRISPRidentify-v2-2.3.2.tar.gz
cd CRISPRidentify-v2-2.3.2
./install.sh

./crispridentify \
  --file examples/NC_017040.1.fasta \
  --result_folder runs/example \
  --html_report True
```

The run creates one result directory per input sequence. The canonical
`report.json` is always written; the example also requests a portable
`report.html`.

Check the installed environment or inspect all command-line options with:

```bash
.envs/crispridentify-core/bin/python scripts/check_environment.py --profile core
./crispridentify --version
./crispridentify --help
```

## Core-supported reconstruction

Normal runs use `balanced` with `--reconstruction_policy conservative` by
default. No extra flag is needed to enable the benchmarked reconstruction
behavior:

```bash
./crispridentify --file genomes.fasta --result_folder runs/reconstruction \
  --retain_candidate_evidence True
```

Conservative mode preserves ordinary candidate competition and keeps repair
proposals from joining separate ordinary groups. Conflicting repair proposals
remain available as evidence. Initial geometry ties are recorded.

The reconstruction stage tests directly compatible Alternative candidates
around a Bona-fide core, freshly scores the core and proposed geometry, and
applies source, structural, and applicable tandem checks to each. A unique
qualifying reconstruction may extend the ends or restore omitted internal repeats. Ambiguous best geometries
and collisions with other accepted arrays cause abstention. The displaced
core remains Alternative evidence.

When a donor's repeat boundaries differ slightly from its core, added repeats
can use a phase reference established unanimously by the donor copies mapped
to that original core. This keeps the identity threshold unchanged and avoids
using added repeats to justify themselves. The source-bound proof records the
reference; all full-array checks still apply.

The new geometry keeps its **own** classifier score, which may be below 0.75;
its Bona-fide classification records the separately supported core. The score
is not a calibrated probability. `reconstruction_evidence.json` records all
decisions; selected arrays carry source-bound evidence in `report.json`.
See [the output contract](docs/OUTPUTS.md#core-supported-reconstruction)
for thresholds, negative-check scope, and replay limits.

After the Bona-fide reconstruction pass, an accepted Possible core can also
authorize a checked longer retained path. This second pass keeps the result
Possible and requires the longer path's own fresh score to reach 0.60, along
with the same source, structure, negative, occurrence and collision checks.
It preserves the original core and unselected alternatives as Possible
Discarded evidence. Tied distinct geometries remain unresolved.
`possible_reconstruction_evidence.json` records this separate pass; the full
path's score never replaces or raises the original core's confidence tier.

A separate contextual check can retain a generated whole path containing one
internal block of identical spacers. The full path and both maximal flanks
must each reach an own score of 0.50 and pass their applicable checks. The
result remains Possible even when its score exceeds 0.75. Tied top paths with
shared repeat starts retain their complete endpoint alternatives; no flank
score breaks the tie. See [localized spacer context](docs/OUTPUTS.md#localized-identical-spacer-context)
for the exception's scope and evidence.

Reconstruction can improve the completeness of an existing array without
detecting a new locus. Its classification basis and retained alternatives
should be considered alongside the model score; default status does not
establish biological precision or recall.

Use `--reconstruction_policy legacy` to disable the reconstruction stage and
restore its previous competition policy. The `full` and `region-screen`
presets, and explicit `--architecture_rescue True` runs, retain `legacy` unless
an incompatible reconstruction request is supplied, in which case the CLI
rejects it. Both resolved policies are recorded in the effective configuration.
Candidate-evidence retention remains independently opt-in.

## Choose a preset

CRISPRidentify exposes three named, provenance-recorded presets:

| Preset | Use it for | Important consequence |
| --- | --- | --- |
| `balanced` | General high-recall analysis with standard reports. This is the default. | Uses optimized fast search and conservative core-supported reconstruction, retaining orientation prediction, degenerated-array search, repair paths, and explicit flag customization. |
| `region-screen` | High-throughput biological-region screening and database ingestion. | Emits machine-profile JSON for accepted arrays only and disables optional work. Candidate representation, category, and exact boundaries may differ from `balanced` or `full`. |
| `full` | Historical full-search comparison. | Runs the substantially more expensive non-fast search with standard output. |

`balanced` is selected when `--preset` is omitted. The `region-screen` and
`full` presets are atomic: an explicitly supplied preset-controlled option that
conflicts with their contract is rejected instead of being silently
overridden. Cas and insertion-sequence annotations remain opt-in for
`balanced` and `full`.

For an individual high-throughput screening check:

```bash
./crispridentify \
  --file examples/NC_017040.1.fasta \
  --result_folder runs/example-screen \
  --preset region-screen
```

For a cohort, use the restartable batch runner rather than a shell loop or the
detector's ad hoc folder mode:

```bash
./crispridentify-batch manifest \
  --input-root genomes \
  --output manifests/genomes.tsv

./crispridentify-batch run \
  --manifest manifests/genomes.tsv \
  --run-dir runs/region-screen \
  --preset region-screen \
  --workers 12 \
  --cpus-per-worker 2
```

Re-running the exact `run` command safely validates and skips completed
samples, retries bounded failures, and resumes the same immutable
configuration. Use a new run directory when the configuration changes.
`--workers * --cpus-per-worker` must fit the CPU affinity granted by the
scheduler.

The parallel native matcher starts fresh Python worker processes to avoid
copying retained arrays from the parent. Direct Python callers of `Pipeline`
or `Detection` must use an importable script and put execution inside
`if __name__ == "__main__":`. For notebooks or other interactive callers, use
the CLI or disable parallel matching with `--parallel False`.

For four hosts, use the same frozen manifest, `--shard-count 4`, one unique
`--shard-index` per host, and a distinct run directory for each shard. See the
[cluster execution guide](docs/CLUSTER.md) for the complete restart, resource,
failure, and CasAndra contracts. See [Resources and
performance](docs/PERFORMANCE.md) before selecting a production preset or
concurrency level.

## Add Cas annotation

Cas annotation is optional and disabled by default. Enable the bundled
CasAndra provider with:

```bash
./crispridentify \
  --file examples/NC_017040.1.fasta \
  --result_folder runs/example-with-cas \
  --cas True \
  --html_report True
```

No separate Cas tool installation is required. CasAndra uses pinned
Pyrodigal, PyHMMER, and scikit-learn dependencies plus the bundled protein
profiles and cassette model. After installation, prediction runs offline.

`--cas_gene_mode auto` is the default. Use `single` for one complete genome or
`meta` for short or metagenomic contigs. Each CasAndra genome annotation has a
one-hour limit by default. Override it with a positive number of seconds, for
example `--cas_timeout_seconds 7200`. A timeout terminates the complete
CasAndra process group.

CasAndra Type VII cassette calls without subtype support are explicitly marked
experimental. The bundled model is a deployment refit, not an independent
evaluation model, and its training records must not be reused to estimate its
accuracy. `--architecture_rescue True` controls a separate, source-backed
CRISPR array analysis; it is not required for CasAndra.

## Understand the results

The report model uses five canonical category names; a particular report only
contains records emitted by its preset and output profile:

| Category | How to use it |
| --- | --- |
| **Bona-fide** | Accepted, highest-confidence CRISPRidentify calls. This is not experimental proof of an active CRISPR system. |
| **Possible** | Accepted calls retained when biological-region recall matters. |
| **Alternative** | An alternate representation associated with a candidate region; do not automatically count it as a distinct biological locus. |
| **Possible Discarded** | A non-accepted candidate. |
| **Low score** | A non-accepted review candidate. |

For candidate-selection research, add `--retain_candidate_evidence True`.
Each sequence then also receives `candidate_evidence.jsonl.gz`, retaining
all five final internal categories, including Alternative and Possible
Discarded, with scores, source coordinates, alignments and structural evidence.
These are competing research representations, **not additional accepted
arrays**. Normal reports and accepted counts are unchanged. The sidecar also
works with `--accepted_only True`; see the [research evidence contract](docs/OUTPUTS.md#candidate-research-evidence).

Low-score arrays may carry optional provenance-bound tags. The `tandem` tag
identifies a final accepted-row demotion under the complete-link
low-spacer-diversity rule; its scientific disposition is `likely_tandem`.
Current reports preserve a source-replayed audit ledger for these decisions,
and the HTML renders it as a visible **Tandem demotion audit** table even when
accepted-only projection omits the Low-score row. The writer exhaustively
reruns the transaction over all pre-quarantine accepted rows, and independent
HTML/benchmark consumers rescreen accepted architectures before trusting an
empty ledger.

The canonical accepted set is Bona-fide plus Possible. When comparing tools,
collapse overlapping representations into biological regions under a
predeclared matching rule rather than treating every representation as a
separate locus. `certainty_score` is a model score, **not a calibrated
probability**. Score bands are only the initial partition; structural
filtering, within-region competition, and documented source-backed
transactions can affect the final category.

Source-backed boundary repair can remove weak terminal alignment columns
before structural rejection. A uniquely supported, longer path of identical
repeats can then replace its accepted short core when its freshly computed
score meets the incumbent category's ordinary threshold. This decision uses
the `source_backed_terminal_repair_completion` classification basis and
retains the original core as lower-category evidence. See the
[output contract](docs/OUTPUTS.md#source-backed-repeat-boundary-repair).

Internal repeat completion searches anomalous spacers for a uniquely supported
source path. A checked replacement preserves the existing score thresholds,
all structural filters and neighboring accepted events. Ambiguous boundaries
remain unresolved. See the [internal completion contract](docs/OUTPUTS.md#source-backed-internal-repeat-completion).

Public source intervals are one-based and end-inclusive on the input forward
strand. Predicted orientation is a separate field. A report validation status
of `valid` means the array was reconstructed from the supplied source sequence;
it does not establish biological truth, classification correctness, or score
calibration.

Start with these fields when consuming `report.json`:

- `schema`: report name and schema version;
- `source`: sequence identity, length, and uppercase-sequence SHA-256;
- `run.configuration`: release, named preset, preset provenance, and effective
  scientific configuration;
- `run.parameters`: the complete detector parameter mapping;
- `summary`: array and category counts;
- `arrays`: coordinates, repeats, spacers, category, optional tags, score, and provenance;
- `cas`: `not_requested`, failure, or computed CasAndra results;
- `extensions`: declarative interpretation and optional topology/locus views.

An array with `boundary_family.status = boundary_uncertain` reports shared
repeat positions with several admissible end assignments. Its coordinates,
sequences and score show one named actual member. Read the complete coupled
alternatives in `classification_basis.boundary_family` before treating that
view as a unique segmentation; they remain available with candidate retention
disabled.

Conservative reports also disclose exact terminal competitors in
`extensions.competing_terminal_extents`. These are existing unaccepted
Possible Discarded observations, with their own scores and explicit missing
checks or blockers. They do not change accepted arrays or repeat counts.
Read the [competing-extent contract](docs/OUTPUTS.md#unaccepted-competing-terminal-extents)
before using them; the full evidence is present even with candidate retention
disabled.

See [Outputs and schema](docs/OUTPUTS.md) and the machine-readable
[`CRISPRidentify-report` 1.2.0 JSON Schema](schemas/report-1.2.0.schema.json)
before building a database parser. A detector run publishes `report.json`
transactionally after its other artifacts. It invalidates any stale copy
before the first pipeline stage and writes the replacement strictly last. For
`crispridentify-batch`, use
the validated per-sample `completion.json` described in the
[cluster execution guide](docs/CLUSTER.md) as the completion marker.

With the standard profile, accepted-array text/CSV reports are also produced.
`--html_report True` adds an offline report. `--fasta_report True` adds FASTA
exports. With `--cas True`, the result also contains `Cas_systems.csv` and
`Cas_genes.csv`, and the JSON `cas` section records the provider, model bundle,
coordinates, classifications, and evidence.

## Input and optional analyses

Exactly one input selector is required:

- `--file PATH` for a FASTA file, including a multi-record file;
- `--input_folder PATH` for a directory of single-record FASTA files;
- `--input_folder_multifasta PATH` for a directory of multi-record FASTA files.

Relative input, output, and manifest paths are resolved from the directory in
which the launcher is called. The launcher itself can therefore be invoked by
absolute path from anywhere.

Array orientation prediction is enabled by default. Optional insertion-sequence
annotation requires a larger environment:

```bash
./install.sh --profile is-elements
export CRISPRIDENTIFY_ENV="$PWD/.envs/crispridentify-is-elements"
```

Keep `CRISPRIDENTIFY_ENV` set when launching either `crispridentify` command
from that installation.

Use a strict topology manifest only when its entries are backed by the frozen
INSDC source records expected by `--sequence_topology_manifest`. Sources not
covered by the manifest retain unknown topology.

## Installation and reproducibility

The extracted release is approximately 90 MiB. The locked core environment is
approximately 2.3 GiB; allow at least 5 GiB for the environment and Conda
package cache. The core lock is Linux x86-64 only and includes the validated
GCC toolchain, so no host compiler is required.

`./install.sh` creates or reconciles `.envs/crispridentify-core` from
`locks/conda-linux-64-core.lock`, verifies the vendored CasAndra source
manifest, installs CasAndra without contacting a Python package index, builds
the native matchers, runs `pip check`, and validates the required executables,
databases, models, and Cas profiles. Use `./install.sh --no-lock` only when you
intentionally accept a fresh solve from `environment.yml`.

For reproducible work, record all of the following with the result bundle:

- `./crispridentify --version` and the release-archive SHA-256;
- the SHA-256 of `locks/conda-linux-64-core.lock`;
- the input accession or file and its sequence checksum;
- the complete command line and hardware allocation;
- the `report.json` schema version and `run.configuration` object.

The compact array classifier uses ten inputs: repeat length and count, repeat
similarity and A/T fraction, mean spacer length, spacer similarity,
repeat-consensus mismatch count, spacer-length evenness, consensus-repeat RNA
folding energy, and similarity to the bundled verified-repeat BLAST database.
It does not use the former ORF, tandem-HMM, or second-BLAST classifier features.

To uninstall a release-local environment, remove
`.envs/crispridentify-core` and any result directories you created. A custom
environment is selected with `CRISPRIDENTIFY_ENV=/absolute/prefix`.

## Provenance, citation, and licensing

CRISPRidentify v2 is a downstream continuation of the original
[BackofenLab/CRISPRidentify](https://github.com/BackofenLab/CRISPRidentify)
software. Cite the original method paper and report the exact v2 release version
and archive SHA-256 used:

> Mitrofanov A, Alkhnbashi OS, Shmakov SA, Makarova KS, Koonin EV, Backofen R.
> CRISPRidentify: identification of CRISPR arrays using machine learning
> approach. *Nucleic Acids Research*. 2021;49(4):e20.
> <https://doi.org/10.1093/nar/gkaa1158>

Machine-readable citation metadata is available in [`CITATION.cff`](CITATION.cff).

The bundled CasAndra directory is a modified downstream snapshot, not an
official upstream release. The supplied source identified itself as
`0.1.0.dev0` but contained no VCS metadata, immutable revision, or published
artifact locator; this repository therefore does **not** assert an upstream
commit or release. The integrated version is
`0.1.0.dev0+crispridentify.3`. Known origin limitations, downstream changes,
model provenance, and checksums are recorded in
[`vendor/casandra/VENDORING.json`](vendor/casandra/VENDORING.json),
[`vendor/casandra/SOURCE_MANIFEST.sha256`](vendor/casandra/SOURCE_MANIFEST.sha256),
and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

CRISPRidentify is distributed under the [MIT License](LICENSE). CasAndra is
licensed under Apache-2.0; its complete license and modification notices are
retained below `vendor/casandra/`. Review dependency licenses before
redistributing a complete Conda environment.

## Help and contributing

- Installation and use questions: [`SUPPORT.md`](SUPPORT.md)
- Reproducible bug reports and changes: [`CONTRIBUTING.md`](CONTRIBUTING.md)
- Security disclosures: [`SECURITY.md`](SECURITY.md)
- Project history: [`CHANGELOG.md`](CHANGELOG.md)
- Community expectations: [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md)
