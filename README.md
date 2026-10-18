# CRISPRloci v2

The complete publication source repository is
[Alexander-Mitrofanov/CRISPRloci-v2](https://github.com/Alexander-Mitrofanov/CRISPRloci-v2).
It includes the frontend, API, isolated worker, Integration adapters and all
17 scientific modules at the exact revisions in `tools/tools.lock.json`.
See [source availability and reproduction](docs/source-availability.md) for
the module inventory, licenses, installation and archival instructions.

Every mode provides a [compact `result.json`](docs/module-results.md).
Every mode also has an [independently reproduced precomputed example](docs/precomputed-examples.md),
including separate RNA-repeat and Type V tracrRNA examples.
See [composable analysis outputs](docs/module-exchange.md) for typed sequence
ports and detection → repeat analysis/mapping or viral search handoffs.
The [component connectivity audit](docs/component-connectivity-audit.md) records
module contracts, supported connections, repairs and verification limits.

A standalone CRISPR analysis server: Vue frontend, FastAPI API, durable SQLite
queue, isolated scientific worker and updateable GitHub tool releases.
It preserves the preceding server's detection, reconstruction, orientation,
recovery links, cancellation, retention and result downloads, and adds Cas and
tracrRNA analyses, source-verified leader contexts and repeat structure evidence.

Public web server: [CRISPRloci](https://alexander-mitrofanov.github.io/crispr-evolution-web/).
The frontend runs on GitHub Pages and the backend on de.NBI; see the
[production installation record](deploy/PRODUCTION.md).

| Analysis | Scientific tools | Results |
| --- | --- | --- |
| Observed array comparisons | CCTK | Ordered spacer occurrences, explicit gaps, sharing metrics, components and isolates from detection arrays.json |
| Repeat subtype evidence | RepeatTyper | Frozen repeat-model predictions with native scores, uncertain and unsupported sequences |
| Spacer-group associations | SpacePHARER | Explicit groups, versioned phage reference databases, native target/control calibration and supporting alignments |
| Repeat reference mapping | CRISPRmap v2 | Exact annotation recovery and nearby repeat candidates with explicit uncertainty |
| Target DNA → reference spacers | CRISPRspacer | Exact full-length observed-spacer matches on both strands |
| Spacers → viral references | CRISPRspacer + BLAST+ | Full-length ungapped matches with a bounded mismatch policy |
| Complete locus annotation | CRISPRidentify v2 + CasAndra + CRISPRtracrRNA v3 + CRISPRleader v2 | Arrays, Cas cassettes, tracrRNA evidence and both adjacent context windows |
| Leader context | CRISPRidentify v2 + CRISPRleader v2 | Source-verified flank windows; leader prediction unavailable |
| Repeat sequences | CRISPRrepeat + ViennaRNA | Isolated DNA/RNA folds; explicit strand alternatives |
| Repeats in array context | CRISPRidentify v2 + CRISPRrepeat + ViennaRNA | Per-instance isolated/local pair support and accessibility |
| Cas systems | CasAndra | Cas protein calls and classified cassettes |
| tracrRNA screening | CRISPRtracrRNA v3 | Type II or Type V-K covariance-model candidates |
| Array detection | CRISPRidentify v2 | Categorized arrays and native reports |
| Evolution and reconstruction | CRISPRidentify v2 + SpacerPlacer | Related-array groups, ancestral spacer states and events |
| Orientation-aware evolution | CRISPRidentify v2 + evOr | Input/reverse comparison and selected reconstruction |
| Locus neighborhood comparison | clinker-v2 | Source-bound translated genes, ordered protein alignments and synteny links |
| Prophage context | Phigaro-v2 | Pinned-model candidate intervals and explicit unassessed short records |
| Insertion-sequence context | ISEScan-v2 | Source-bound insertion intervals, native family evidence and terminal repeats |
| Defense-system context | DefenseFinder-v2 | Pinned defense-model systems/components, supporting HMM hits and separate Cas context |
| RNA family analysis | RNAalifold-v2 | Source-bound alignment, consensus folding and compensatory evidence |

No scientific packages are imported into the API. SpacerPlacer and evOr have
separate environments. A complete locus run reuses the array and Cas evidence
produced by CRISPRtracrRNA's complete workflow. Backend validators independently
bind coordinates, sequence checksums, counts and retained evidence to the input.
Scientific uncertainty and missing values remain explicit.

Every analysis has a collapsed **Analysis options** menu with applicable tool
flags, defaults, effects and CLI documentation. Selected overrides survive queue
recovery and appear in result provenance. See the [tool flags and modes reference](docs/tool-options.md)
for all seventeen pinned CLIs and the distinction between scientific settings and
operator-managed arguments.

The **Database** page searches a separate CRISPRcasDb reference catalogue by
complete accession, repeat, spacer or Cas gene. Records appear only after a
query is submitted and link genomes, arrays, occurrences and Cas annotations.
The read-only API uses checksum-verified compressed shards and an indexed genome
directory; scientific jobs never download or modify reference data. See
[catalogue import and operation](docs/catalogue.md) and the
[agent council's scientific/design decisions](docs/catalog-council.md).

## Run locally

Requirements: Linux x86_64, Python 3.12, Node 24, Git, Conda/Mamba, a C++ build
toolchain. The complete publication repository bundles all scientific sources;
the separate development layout requires GitHub access to their origins. Build in the final
project location: installed environments contain absolute paths.

```bash
# Only for the separate development layout without bundled sources:
# gh auth login
# gh auth setup-git

bash scripts/setup.sh
python3.12 scripts/provision_tools.py
python3.12 scripts/manage_release.py build
# The build prints its absolute release directory; use it below.
python3.12 scripts/gate_release.py /absolute/path/to/release
python3.12 scripts/manage_release.py activate /absolute/path/to/release

# Terminal 1: API and scientific worker, using .runtime-data/.
python3.12 scripts/serve.py

# Terminal 2: frontend.
cd frontend
npm run dev -- --host 127.0.0.1
```

Open `http://127.0.0.1:5173`. The backend listens on `127.0.0.1:8000`.
The seventeen-tool release gate defines 37 real scientific jobs through HTTP
across all 20 modes, plus 5 input rejections before dispatch. It covers
positive/negative controls, both tracrRNA model families, repeat DNA/RNA inputs,
empty results, authenticated downloads and archive integrity.
It stores local evidence under `.runtime/gates/` and approves only that exact
release. `serve.py` verifies source, installed packages, native executables and
models before startup. It refuses an ungated release.

Sources are pinned in [tools/tools.lock.json](tools/tools.lock.json).
See the [October 13 refinement](docs/tool-refinement-2026-10-13.md) for the current review and verification scope, and the [October 11 review](docs/tool-review-2026-10-11.md) for the earlier release. Built
scientific environments live in `.runtime/toolsets/`; backend and Integration
wheels live in `.runtime/releases/`. Jobs use these snapshots, never a moving
checkout. Source checkouts are restored by `python3 scripts/sync_tools.py`.

## Update a tool

```bash
python3.12 scripts/update_tools.py --tool 'evOr@main'
# Use the candidate lock path printed by that command:
python3.12 scripts/provision_tools.py --lock /absolute/candidate.json
python3.12 scripts/manage_release.py build --lock /absolute/candidate.json
python3.12 scripts/gate_release.py /absolute/new-release
# Let jobs finish, then stop the API and worker.
python3.12 scripts/manage_release.py activate /absolute/new-release
python3 scripts/sync_tools.py
python3.12 scripts/serve.py
```

Updates require operator GitHub access. Candidate resolution checks the upstream
CRISPRtracrRNA dependency lock; incompatible CRISPRidentify/CasAndra revisions
must be updated as a compatible set. Activation refuses while managed services
are running, retains `.runtime/previous`, and updates the reference lock.
`manage_release.py rollback` selects the prior gated v2 release after services
stop. Never run pip/conda installs inside an active release.

## Verification and deployment

```bash
bash scripts/check.sh                  # Python + frontend regression gates
python3 scripts/sync_tools.py --check  # Exact origins, commits, clean checkouts
python3.12 scripts/manage_release.py verify .runtime/current
```

- [Replacement audit](docs/replacement-audit.md): predecessor compatibility, added capabilities and verification limits.
- [Private operations](docs/private-operations.md): local connection records and commit/push safeguards.
- [Operations runbook](docs/operations.md): release, recovery and routine checks.
- [Deployment and migration](deploy/README.md): service settings, backups,
  existing jobs, TLS/Pages and the predecessor-folder retirement checklist.
- [Backend API and limits](backend/README.md).
- [Scientific contract](docs/scientific-contract.md).
- [Architecture](docs/architecture.md).
- [Verification evidence](docs/verification.md): regression, real-science, browser and independence checks.
- [Tool repositories and source provenance](tools/README.md).

Neither predecessor source folder is required at runtime. Their provenance is
recorded in [docs/migration-source.json](docs/migration-source.json). This project
does not delete those folders or replace existing live services automatically.

Repeat folding and array-context workflow: [implementation and operator guide](docs/repeat-analysis.md).

See [Protospacer search and VM spacer reference](docs/protospacer-integration.md) for the new queued workflow,
immutable reference identity and three-record release gate.

See [Bidirectional spacer searches](docs/spacer-searches.md) for the web interface, viral reference search, result interpretation and operator configuration.

See [tool integration audit](docs/integration-audit-2026-10-03.md) and [repeat reference mapping](docs/repeat-mapping.md) for current coverage, remaining workflows and CRISPRmap configuration.

RepeatTyper is integrated as the queued `repeat_type` mode, consuming repeat DNA
FASTA or the existing repeat occurrence handoff. See [repeat typing](docs/repeat-typing.md) for immutable model identity,
uncertainty semantics, the isolated runtime and validated exports.

SpacePHARER is integrated as `spacer_association`: explicitly grouped spacer occurrences, an operator-owned targeted phage panel, native empirical calibration and protected evidence exports. See [spacer association](docs/spacer-association.md).

CCTK is integrated as queued `array_compare`: detection arrays.json feeds observed
spacer order and sharing comparisons. See [array comparison](docs/array-comparison.md).
