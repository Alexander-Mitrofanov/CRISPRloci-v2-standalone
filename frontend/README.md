# CRISPRloci v2 frontend

Static Vue 3 frontend for CRISPRidentify v2, CasAndra, CRISPRtracrRNA v3,
SpacerPlacer, evOr and connected context workflows. It is built with Vite and contains no server-side runtime or secrets.

## Genome exploration

Completed genomic jobs with `result.json` and `input-sequences.json` offer a
**Genome & locus** tab. A single, lazy-loaded igv.js 3.8.10 integration provides
sequence zoom, coordinate navigation, strand-aware tracks and SVG figures.
Feature selection opens the existing scientific evidence alongside the tracks.
The original tables, array comparisons and evolutionary views remain the owners
of detailed interpretation; no second genome-browser library is installed.

Supported modes are detection, complete loci, Cas, tracrRNA, leader context,
protospacer search, prophage context, insertion context and defense context.
After the viewer starts successfully, it replaces the static annotation overview
strip. Precomputed examples, older results and jobs outside the viewer's limits
retain their existing views.

`src/features/genome/` verifies complete artifacts, binds source hashes to actual
target DNA, and converts inclusive annotations and half-open flanks once. Inputs
are bounded to 10 million bases, 50,000 placed features and 1,000 overlapping
features per track. The feature chooser lists at most 100 search matches; it
does not truncate the scientific tracks. Empty or unlocated windows are counted
explicitly. Unknown topology never implies circular adjacency. See the existing
evidence tables for module scope, raw scores and unassessed records.

Only the API client receives the job capability. IGV receives local `File` data
and in-memory features, with external genome defaults, remote searches and URL
session loading disabled. Context menus that offer external BLAT searches are
disabled. No new endpoints, scientific environments or production changes are
required. The MIT license is distributed in `public/licenses/igv.txt`.

The production CSP permits two exact static stylesheet hashes from the pinned
IGV distribution; arbitrary inline styles/scripts, workers and frames remain
blocked. `scripts/igv-vendor.mjs` verifies the upstream source checksum and an
isolated output chunk before the public scan exempts its built-in public genome
identifiers. Application chunks and example data retain the strict identity scan.
Dependency updates require reviewing both this checksum and the CSP hashes.

Run `npm test` and the production build. In the private monorepo, run
`scripts/check.sh frontend` and `scripts/browser_genome.py` against a locally
served build configured with `VITE_API_BASE_URL=https://analysis.example.org`.
The browser test uses synthetic intercepted artifacts and creates no live jobs.

## Version 2.12.0 context workflows

Five additional methods appear when advertised by the analysis service:

- **Compare gene neighborhoods:** clinker compares source-verified proteins from 2–8 related genomic windows, up to 50,000 nt each.
- **Inspect prophage and spacer context:** Phigaro candidate regions are joined to exact source-spacer matches.
- **Inspect insertion-sequence context:** ISEScan predictions retain partial elements, boundary evidence and Cas-cassette proximity.
- **Survey the defense repertoire:** DefenseFinder preserves system calls, component assignments and profile hits with historical model release 1.2.3.
- **Compare an RNA family:** 2–50 homologous, consistently oriented DNA or RNA records, each 30–500 nt, are aligned with MAFFT and folded with RNAalifold. At least two distinct sequences are required.

RNAalifold and the ViennaRNA folding engine are credited to the ViennaRNA authors and the Institute for Theoretical Chemistry of the University of Vienna. The Methods section links the original RNAalifold and ViennaRNA publications.

Results retain missing evidence, scope limitations and source coordinates. These methods have fixed scientific policies and no public tool options in this release. RNA-family homology and transcript boundaries remain researcher supplied. Comparisons, model predictions and proximity do not establish function or activity.

With a build configured for `VITE_API_BASE_URL=https://analysis.example.org` served locally, run `python scripts/browser-context-modes.py --origin http://127.0.0.1:8491 --output /tmp/context-browser` for Chromium acceptance checks. The script uses intercepted synthetic API fixtures and submits no live jobs.

## Structure

The implementation intentionally uses small files so independent changes and
reviews stay local:

- `src/App.vue` composes the page and owns no scientific rendering details.
- `src/components/submission/` contains input, mode, readiness, and option controls.
- `src/components/jobs/` contains capability-link recovery and progress flows.
- `src/components/catalog/` presents the read-only reference database, with state in
  `src/composables/useCatalog.js` and scientific formatting in `src/features/catalog/`.
- `src/components/results/` contains vertical overview, orientation,
  reconstruction, history, and publication slices behind `ResultsView.vue`.
- `src/components/shell/` and `src/components/common/` contain layout and shared UI.
- `src/composables/` owns service configuration and in-memory job session state.
- `src/utils/` contains framework-neutral formatting, Newick, result, and download helpers.

Bearer capabilities stay in memory and in the explicit `#job=` recovery fragment. URL fragments
are client-side and are not included in HTTP requests, so API calls still carry the capability only
in the Authorization header. The complete recovery link is a bearer secret and must be kept private.
Capabilities are never stored in cookies, local storage, or session storage.

## Spacer searches

Choose **Spacer searches** on the main page, then select a direction:

- **Find spacers in my sequence:** submit target DNA to find full-length exact
  reference-spacer matches on both strands, using the NCBI-hosted 2017 research
  spacer collection.
- **Find viruses for my spacers:** submit one DNA spacer per FASTA record and
  choose 0, 1 or 2 substitutions. The server searches viral RefSeq with BLAST+
  and retains full-length ungapped alignments. Eligible spacers are 18–80 nt;
  shorter or ambiguous records are reported as skipped. The displayed service
  record and request limits also apply.

Both modes appear only when their server reference is configured. The two forms
keep separate drafts while switching direction. Results show coordinates, strands,
reference identity, and complete JSON/TSV downloads. Viral results can be filtered
by spacer, accession or title and summarize support from distinct spacer sequences.
BLAST retrieval is heuristic; matches do not establish infection, host range or
functional targeting. PAM compatibility and CRISPR-array overlap are not evaluated.
Reference databases stay on the server and are not bundled with this frontend.

`python scripts/browser_spacer_searches.py --origin http://127.0.0.1:4187 --output
/tmp/spacer-browser` verifies both UI flows against an already running build using
intercepted synthetic API responses. Build with `VITE_API_BASE_URL=https://analysis.example.org`
and serve the Vite preview on port 4187 first. This browser check makes no live jobs.

## Reference database

The main page's **Search the database** option and the **Database** navigation link
open `?view=database` while preserving any active job recovery
fragment. Public catalog requests use the same configured API origin and never include
job credentials. The dataset is served by the API; database files are not bundled in this
static public frontend.

The initial page shows only a search form. Select accession number, repeat, spacer or
Cas gene, enter a specific value, and submit it. No records are fetched before submission;
editing the input or changing the search type clears results and cancels pending reads.
Accessions must be complete and versioned. Repeat/spacer queries match exact nucleotide
sequences globally. Cas queries match original annotations/profiles or explicit gene-name
tokens from those labels, ignoring case. The API's snapshot-pinned search index supports
these global searches. Every result uses bounded pagination; an empty intermediate page
can still have a next page. Coverage, provenance and scientific limits are collapsed
beneath submitted results.

Arrays link to ordered repeat and spacer occurrence pages. Unknown values stay unknown;
deletions use explicit 0-based interbase boundaries rather than invented intervals or
sequences. Cas evidence retains caller scores and flags without treating predictions as
experimentally established. Downloads contain only the currently displayed page: CSV
includes all public record fields, and FASTA includes observed sequences (or explicitly
labelled array repeat consensus). Repeats and spacers are occurrences, not unique families.

From the private monorepo, `backend/.venv/bin/python scripts/browser_catalog.py --origin
http://127.0.0.1:5173` checks the running frontend against a configured, imported catalog
with its exact-search index. It checks the empty initial screen, explicit searches,
genome context, details, downloads, browser history and mobile layout, saving screenshots under
`output/playwright/catalog/`. It starts no services and submits no analysis jobs.

## Bundled examples

The 15 previously released methods have exact inputs and real precomputed outputs in
`public/examples/<method>/`. RNA repeats and Type V tracrRNA have separate
examples selected by their molecule/model option. Load example fills the input
and restores the recorded options. View precomputed result displays the normal
result components without submitting sequences, including when the API is offline.
The five context methods do not offer example buttons until validated native examples are published.

All examples are generated with the pinned tools and reproduced by a second
independent worker job with the same input, options and references. Only execution
durations may vary. Input bytes, normalized sequences, counts and all submitted
options are checked before a stored result is used. Edited inputs or options use
the normal live submission path, which still requires the advertised capability.
The result includes input/output downloads and exact tool provenance.

Reference-dependent examples use small, explicitly labeled demonstration panels;
a fresh run using larger or different references can produce different matches.
Scientific uncertainty and missing values remain explicit. The older masked
orientation teaching assets remain validated for compatibility.

## Independent users

Anonymous jobs use capability authorization instead of shared browser sessions. Every submission receives a random 128-bit job locator and a separate 256-bit bearer token; only the token digest is stored. There is no job-list endpoint. Status, cancellation, artifacts, and bundles all require the token bound to that job, and a wrong job/token pair returns the same 404 as an unknown job.

The frontend writes the active credential to a strict `#job=<locator>.<token>` fragment so reloading
or copying the current page can recover the job. The fragment is never sent to GitHub Pages or the
API, but it is visible in browser history and to scripts on the same origin; anyone holding the full
link can access that job until expiry. API requests omit cookies and use no-store and no-referrer
policies plus an Authorization header. No capability enters an API URL, local storage, or session
storage. Two jobs can coexist in the durable queue and complete independently, including when users
share one NAT address under the configured per-client allowance.

## Development

```bash
npm ci --ignore-scripts
npm test
npm run dev
```

Build the deployable Pages artifact with an exact HTTPS API origin and project
base path:

```bash
VITE_API_BASE_URL=https://analysis.example.org \
VITE_BASE_PATH=/YOUR-FRONTEND-REPOSITORY/ \
npm run build
```

`VITE_API_BASE_URL` must be an origin without a path, query, credentials, or
fragment; production builds fail closed when it is missing or invalid.
`VITE_BASE_PATH` is normalized to one leading and trailing `/`, so both
`YOUR-FRONTEND-REPOSITORY` and `/YOUR-FRONTEND-REPOSITORY/` produce the same base.

The public repository contains the contents of this directory at its root.
Consequently, `.github/workflows/pages.yml` here is intentionally shaped for a
frontend-only checkout; it is not the private monorepo workflow. Export through
`../deploy/export-pages-frontend.sh` so dotfiles are included and private
backend/scientific directories cannot cross the publication boundary.

The build also runs scripts/scan-public-example.mjs. The scanner validates every example and its exact input binding in public and dist, permits only the declared sequence assets, and rejects source identities, capabilities and private paths. Synthetic viral identifiers and null reference accessions are allowed only within validated example contracts.

Production origin is chosen at deployment; no production site is changed by this migration.

API origin: configured by the operator at build time.
