# Guided standalone setup

The complete public repository is
[CRISPRloci-v2-standalone](https://github.com/Alexander-Mitrofanov/CRISPRloci-v2-standalone).
It includes all scientific source modules as ordinary files, so cloning or
unpacking the source does not require access to any private repository.

## Start the graphical installer

```bash
git clone https://github.com/Alexander-Mitrofanov/CRISPRloci-v2-standalone.git
cd CRISPRloci-v2-standalone
bash install.sh
```

Requirements before opening the wizard: Linux x86-64, Python 3.10 or newer,
a browser, network access and at least 35 GiB free disk space. Git is needed only
when cloning. Unsupported platforms are identified before installing anything.
An x86-64 Linux environment in WSL2 can run the scientific stack; native Windows,
macOS and ARM environments are not supported by these pinned scientific packages.

The GUI uses Python's standard library and local HTML/CSS/JavaScript. It does not
require Node, Conda, a web framework or an external website to open. Its address
is printed in the terminal if automatic browser opening is unavailable:

```bash
bash install.sh --no-browser
bash install.sh --no-browser --port 8876
```

The installer listens only on `127.0.0.1`. Its session URL contains a random
fragment required to authorize setup actions. Keep the installer terminal open.
The scientific job API cannot install software or change reference databases.

## Choose references and resources

Every installation includes the web interface, durable local queue, isolated
worker and all 17 scientific modules. Choose CPU threads per scientific job and
whether to include developer dependencies and regression checks. Required
scientific models, including pVOGs and DefenseFinder models, are part of the
software installation and are checksum-verified.

| Reference | Bundled example | Public download | Local construction/import |
| --- | --- | --- | --- |
| Observed spacers | Three original sample records | Shmakov et al. 2017 spacerome, pinned SHA-256, 720,391 source occurrences, about 13 MiB compressed | Build a SQLite spacer index from your FASTA or gzipped FASTA |
| Viral DNA | Three synthetic targets | Freeze current NCBI viral RefSeq metadata, download its complete supported genomic FASTA inventory, verify official checksums and build a BLAST index | Build a BLAST reference from your FASTA or gzipped FASTA |
| Spacer-group associations | Three bundled phage genomes | Acquire complete DNA phage/prokaryotic-virus records with the documented NCBI selection policy and build the SpacePHARER target/control reference | Build from a directory of `.fa`, `.fna` or `.fasta` files, optionally gzipped |
| Repeat mapping | Two synthetic reference entries | Original 2013 CRISPRmap assets, checksum-verified before importing | Use a compatible CRISPRmap `reference.json` |
| Searchable catalogue | A clearly labelled synthetic example, generated and indexed locally | No automatic mirror of the hosted production catalogue | Use a compatible complete snapshot directory, or import a CRISPRcasDb export bundle with its SHA-256 and build the search index |

Example collections are functional demonstrations, not comprehensive biological
references or scientific benchmarks. Their limited scope remains explicit in
reference provenance. There is no requirement to download large optional datasets.
For local inputs, enter an absolute path. Catalogue snapshot directories must
retain the registry, manifest, compressed shards and receipts together. A source
export bundle is the format accepted by `scripts/import_catalog.py`, not an
arbitrary SQLite or tar file. Local reference paths must not traverse symlinks.
The installer validates each completed reference using its scientific tool.

Large NCBI viral/phage downloads may need many additional GiB and hours of work.
Acquisition code retains disk reserves and refuses changed snapshots rather than
silently switching versions. The review screen describes selected data sources
before any download begins. Original methods, source datasets, models and native
dependencies retain their own [licenses and attribution](../THIRD_PARTY.md).

## What the pipeline does

1. Check platform, free disk, exact scientific source pins and stopped local
   services. No production service is contacted or stopped.
2. Download the checksum-pinned Micromamba 2.3.3 binary and create a private
   Python 3.12, Node 24, Git and C++ build environment under `.installer/`.
3. Install application/interface dependencies and build every scientific tool
   in its own environment. SpacerPlacer and evOr remain separate because both
   own `sp_model`; IntaRNA also keeps its separate native environment.
4. Build and verify an immutable application release, then construct and validate
   the selected reference databases.
5. Run the real-science HTTP acceptance gate for the release. The wizard cannot
   activate a failed or ungated release.
6. Start the candidate on an isolated local port and confirm all 20 analysis
   modes and the selected catalogue are available through the API.
7. Activate the validated local release and save its reference paths/checksums,
   CPU setting and data directory in private installer state.

A download bar reports actual received bytes when an HTTP total is available.
Without a total it displays bytes received with an active bar. Spacer indexing
reports source record counts, phage acquisition reports records acquired,
catalogue indexing reports committed shards, and scientific provisioning reports
completed tool environments. Native operations such as BLAST construction do
not provide a reliable total and use an active bar. Overall progress counts
completed stages, not an estimated elapsed-time percentage. Details remain
available in the bounded installation log.

## Run the local application

```bash
bash start.sh
# Do not open a browser automatically:
bash start.sh --no-browser
# If the default local ports are already used:
bash start.sh --api-port 8001 --frontend-port 5174
```

The local Vite interface uses `127.0.0.1:5173`, and the API uses
`127.0.0.1:8000`. Both bind to loopback only. The launcher starts the validated API
and isolated worker, waits for readiness and opens the interface. It discards
inherited deployment/reference configuration and uses the choices saved by the
installer. Press Ctrl+C to stop all its child processes. The graphical installer
also provides explicit Start/Stop controls for an application it started.

Scientific jobs retain capability-protected downloads, cancellation, bounded
resources and retention. Their private files are in `.runtime-data/standalone/`.
This launcher is for a local workstation, not a public deployment. Existing
production deployment/export commands and HTTPS build rules remain unchanged.

## Cancel, retry and diagnose

Use **Cancel installation** to terminate the current command and its process
group. Verified completed downloads, tool environments and databases are retained.
Run `bash install.sh` again, select the same choices and retry. Completed assets
are verified before reuse; an incomplete application build can resume for the
same exact source identity. A retry after the managed bootstrap has completed
requires at least a 5 GiB free-space reserve; unfinished downloads and remaining
environments may require more. A different source lock starts a fresh plan.
Failed or cancelled runs never mark an installation
ready. Changed cached downloads, tool receipts or local inputs are rejected.

Stop a running local application before changing its installation. Do not move
the checkout after installation: native environments and receipts contain absolute
paths. If software installation fails, open the installation log to identify the
failed package, unavailable source or disk-space limit. `.installer/last-install.json`
retains the latest bounded log/state, and `.installer/application.log` records
startup messages when the GUI starts the app. These files are ignored and must
not be published. Do not delete an active `.runtime/` directory to retry.

Check prerequisites without installing:

```bash
bash install.sh --check
```

## Headless use with the same pipeline

Prepare a JSON choices file. Omitting a reference retains its bundled example.

```json
{
  "cpu": 4,
  "developer": false,
  "references": {
    "spacer": {"mode": "download"},
    "viral": {"mode": "sample"},
    "association": {"mode": "sample"},
    "repeat": {"mode": "historical"},
    "catalog": {"mode": "sample"}
  }
}
```

```bash
bash install.sh --headless /absolute/path/choices.json
bash start.sh --no-browser
```

Allowed reference modes are those in the table above: `sample`, `download` or
`local` for spacer/viral/association; `sample`, `historical` or `local` for repeat;
`sample`, `local` or `import` for catalogue. A `local`/`import` entry requires
`path`; catalogue `import` also requires `sha256`. Unknown options and arbitrary
commands are rejected.

## Developer verification

```bash
bash scripts/check.sh
backend/.venv/bin/python -m playwright install chromium
backend/.venv/bin/python scripts/browser_installer.py
```

Browser acceptance tests exercise the real local wizard and explicitly stub
installation actions. They do not claim to perform a fresh 17-tool install.
Native reference builders, downloader integrity, HTTP isolation, cancellation,
release gating and failure handling have separate regression tests. A complete
fresh scientific installation requires adequate disk/network resources and is
validated by its own real-science release gate before activation.

The repository also includes a manually dispatched **Complete standalone
installation** GitHub Actions workflow. It installs the private bootstrap,
all 17 scientific environments and sample references on a disposable Linux
runner, runs the real scientific release gate, checks all analysis modes through
the local API and activates the result. It does not use an existing runtime or
contact a production service. The ordinary CI run tests the graphical wizard
separately in Chromium.
