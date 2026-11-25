# Working in CRISPRidentify

This directory is the Git repository. Run commands here. Sibling benchmark,
profiling, output, and research-handoff directories are outside the repository;
do not move, clean, or edit them as part of tool maintenance.

## Find the code

Read [docs/CODE_MAP.md](docs/CODE_MAP.md) for task-to-file and test mappings.
Start in the relevant package, rather than reading the entire pipeline:

| Task | Start here |
| --- | --- |
| Flags, defaults, presets | `components/cli/` |
| Pipeline stage ordering | `components/pipeline.py` |
| Repeat discovery and fuzzy matching | `components/detection/` |
| Candidate objects, filters, repeat repair | `components/refinement/` |
| Features, classifiers, score caching | `components/scoring/` |
| Lineage, competition, occurrence protection | `components/evidence/` |
| Core-supported reconstruction and extent alternatives | `components/reconstruction/` |
| Source-backed rescue, tandem checks, circular joins | `components/rescue/` |
| Cas, strand, leader, insertion-sequence annotation | `components/annotation/` |
| JSON/HTML/CSV output and publication validation | `components/reporting/` |
| FASTA, coordinates, source validation, topology | `components/sequence/` |
| CPU/process limits and folder scheduling | `components/runtime/` |
| Restartable cohort runner | `scripts/run_crispridentify_batch.py` |

Use scoped searches, for example:

```bash
rg -n 'def |class ' components/reconstruction
rg -n 'reconstruct_possible_arrays' components/reconstruction tests
```

## Keep the structure useful

- Import canonical package paths. `components/_compat.py` only serves historical
  flat imports and old pickles; do not add new implementation there.
- Keep package `__init__.py` files lightweight. Native CPU limits must be set
  before importing NumPy or starting the pipeline; spawned workers must not
  parse arguments, change directories, or launch the CLI on import.
- Use `components/paths.py` for checkout assets. Never infer the repository root
  from the depth of a feature module or from the caller's working directory.
- Keep scientific decisions in their owning package and stage ordering in
  `pipeline.py`. Preserve thresholds, coordinates, evidence and transaction
  checks during structural refactors.
- Update the code map when moving or adding a subsystem. Tests use canonical
  imports; `tests/test_package_layout.py` covers legacy import compatibility.
- Keep generated outputs in a temporary directory or ignored `runs/` directory.
  Vendored CasAndra and bundled model/data files have separate integrity
  contracts; see [CONTRIBUTING.md](CONTRIBUTING.md).

## Verify

Use the existing locked interpreter; installation instructions are in the
README if it is missing. Run focused tests while editing, then:

```bash
.envs/crispridentify-core/bin/python -m unittest discover -v
ruff check .
git diff --check
```

For runtime, CLI, asset-path, or reporting changes, also run the bundled-genome
smoke command in CONTRIBUTING.md with a fresh output directory. Do not run or
rewrite external benchmarks for a structural cleanup.
