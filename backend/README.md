# CRISPRloci v2 backend

Backend 2.11.0 supports compact, independently validated SpacePHARER 2.1.0 broad-reference results and bounded lazy source validation. See [spacer associations](../docs/spacer-association.md).

Backend 2.6.0 adds a primary [compact `result.json`](../docs/module-results.md)
for all 12 modes. Full validated evidence supplies the export independently
of preview limits. Detailed artifacts and sequence handoffs remain available.

Backend 2.5.0 adds [typed JSON handoffs](../docs/module-exchange.md): compact
array evidence, reusable repeat/spacer sequence sets, and JSON inputs through
the existing job submission endpoint. Scientific tool contracts are unchanged.

Backend 2.2.0 supports the optional queued `protospacer` mode. It searches
submitted DNA against an operator-installed observed-spacer collection with
exact full-length comparisons on both strands. Set both
`CRISPR_WEB_SPACER_DATABASE=/absolute/reference-directory` and
`CRISPR_WEB_SPACER_DATABASE_SHA256=<SHA-256 of spacers.sqlite3>` in the operator
environment. The scientific release must pin CRISPRspacer (manifest schema
5.0.0); its isolated command belongs in the operator runner configuration.

Without this configuration the mode is absent from the advertised capabilities
and explicit submissions return 503. The public request schema accepts neither
database paths nor executable commands. Submission pins the configured digest
in the durable queue; a later worker configuration change fails that job rather
than silently changing its reference. The standalone worker verifies the actual
database bytes. Reference refreshes are offline operator actions using a new
directory and digest; pending jobs retain their original identities.

This mode accepts at most 5,000,000 input bases, subject to the server's smaller
general limits, and inherits queue admission, per-job capabilities, cancellation,
process containment, retention and storage limits. Native search is capped at
300 seconds and 10,000 matches. Resource-limited searches fail explicitly.
Independent backend validation reconstructs each match from immutable target DNA
before publishing `protospacer-analysis.json` and `protospacer-matches.tsv`.
The preview is limited to 1,000 rows and 700,000 encoded bytes; its truncation
does not limit the downloads. Native retained inputs and internal manifests are
excluded from publication. Sequence evidence does not establish biological
targeting; PAMs remain unevaluated and array overlap remains unknown.

Production-oriented FastAPI and durable worker boundary for the repository's
CRISPRidentify, CasAndra, CRISPRtracrRNA, SpacerPlacer and evOr integration. This directory is
self-contained: it installs none of the scientific tools into its own
environment and does not modify their source trees.

## Architecture

Run the API and worker as separate processes under the unprivileged
`crisprsvc` account:

```text
GitHub Pages -> HTTPS reverse proxy -> FastAPI -> SQLite WAL queue
                                                |
                                     one polling worker
                                                |
                         shell-free crispr-tools Integration CLI
                                                |
                         private job root + reproducibility ZIP
```

SQLite is the source of truth. Queue claims use an immediate transaction and a
renewed lease. A stopped or unexpectedly dead worker can retry only up to
`CRISPR_WEB_MAX_JOB_ATTEMPTS`; exhaustion is terminal rather than an infinite
crash loop. Every attempt gets a new empty output directory, which satisfies
the Integration contract without reusing partial scientific output.

The API, worker, and cleanup command refuse to run as root. The data root,
SQLite database, normalized input, logs, and results use owner-only
permissions. Configure the systemd services with `User=crisprsvc`,
`Group=crisprsvc`, and `UMask=0077`.

## Immutable installation

Use the project-level [setup and release workflow](../README.md) and
[deployment guide](../deploy/README.md). `scripts/serve.py` selects a gated
release and supplies its exact runner and wheel identities. Every release
contains backend 2.0.0 and Integration 0.3.0 wheels plus the five pinned
scientific source and environment receipts. Public provenance exposes repository
URLs and Git commits without exposing installation paths.

The API runs on loopback with Uvicorn proxy-header interpretation disabled.
Configure exact frontend origins and trusted proxy networks in the service
settings. New deployments use `.runtime/releases/` under the new project root.
Existing queue data can be preserved with the transactional schema-3 migration.

## Public API v1

Service endpoints:

- `GET /api/v1/health` — database, worker, storage, and exact scientific
  manifest readiness; returns 503 unless all are healthy. It includes all three
  release digests and six sanitized scientific identities.
- `GET /api/v1/config` — public limits, modes, and fixed policies.
- `GET /api/v1/version` — backend/Integration contract versions, backend release
  and wheel digests, scientific manifest digest, and path-free tool identities.
- `GET /api/v1/openapi.json` — typed OpenAPI contract. Interactive Swagger
  and ReDoc pages are disabled so the API can retain its strict script-free CSP.

Job endpoints:

- `POST /api/v1/jobs` — submit a bounded JSON job.
- `GET /api/v1/jobs/{job_id}` — status and structured summary.
- `DELETE /api/v1/jobs/{job_id}` — cancel a queued/running job.
- `GET /api/v1/jobs/{job_id}/result` — complete results ZIP.
- `GET /api/v1/jobs/{job_id}/artifacts/{artifact_id}` — one registered file.

Submission example:

An optional `tool_options` object accepts exact shipped CLI flag names, grouped
by tool. For example, a `repeat_map` submission can include
`"tool_options": {"crisprmap": {"--max-distance": 0, "--top-k": 5}}`.
Only the [documented bounded options](../docs/tool-options.md) that apply to the
selected mode are accepted; arbitrary paths, commands and resource overrides
are rejected. Omitted flags retain the workflow defaults. The transactional
schema-7 migration adds an empty-options default to existing jobs. Selected
options are returned under `job.options.tool_options` and retained in provenance.

The following example uses default scientific flags:

```json
{
  "sequence": ">isolate_A\nACGTRYSWKMBDHVN...\n>isolate_B\nACGT...\n",
  "filename": "related-isolates.fasta",
  "mode": "orientation",
  "category_policy": "bona_fide_possible",
  "spacer_distance": 1,
  "bias_corrections": true
}
```

These modes are accepted:

- `loci` maps to `--tools loci`: complete array/Cas/tracrRNA evidence.
- `cas` maps to `--tools casandra`: CasAndra genome analysis.
- `tracrrna` maps to `--tools tracrrna`: covariance-model screening.

For `loci` and `tracrrna`, `tracr_model_type` is `II` (default) or `V` (V-K).
Array/Cas annotations and all candidates are exported in `annotations.json`, TSV
and candidate FASTA. Complete loci also retain native array evidence in
`array-evidence.json`. UI previews are bounded; these exports contain all
validated records. Genomic coordinates are 1-based and end-inclusive.

The original evolutionary modes remain available:

- `detection` maps to Integration `--tools identify`.
- `reconstruction` maps to `--tools spacerplacer`.
- `orientation` maps to `--tools crispr-evor`. CRISPR-evOr includes the
  reported SpacerPlacer reconstruction; its `crispr_evor/0_results.csv` is the
  definitive supported history for a decisive call and the input-order reporting
  default for an unresolved call.

When native event exports are available, the bounded public summary adds
`orientation.reconstructions` in this shape:

```json
{
  "group": "repeat_group_001",
  "hypothesis": "input",
  "newick": "(...);",
  "spacer_order": [1, 2, 3],
  "nodes": [{
    "name": "Inner1",
    "spacers": [1, 3],
    "gains": [1],
    "loss_blocks": [[2, 3]],
    "contradictions": [],
    "duplications": [],
    "rearrangements": [],
    "reacquisitions": [],
    "independent_gains": [],
    "other_duplication_events": []
  }],
  "acquisition_count": 1,
  "deletion_count": 2
}
```

`hypothesis` is `input` or `reverse`; both are retained when available.
Spacer IDs are positive canonical numbers that remain stable between the two
hypotheses for one group. Original spacer labels and DNA are withheld from the
public summary.
Node arrays, gains, loss blocks, and special events are model reconstructions,
not directly observed history. A decisive comparison selects its supported
hypothesis. ND uses the input-order history only as the reporting default and
does not support that direction. `reconstructions_truncated` reports a summary
group cap. The browser may further cap visual detail; downloadable artifacts
retain the detailed tool output.

Category policy is either `bona_fide_possible` or `bona_fide_only`. Raw
category counts remain available, while the primary `detection.arrays` list and
selected counts contain only categories chosen for analysis. Spacer distance is
an integer from 0 through 2. Bias correction is a boolean.
No arbitrary arguments, categories, filesystem paths, executables, CPU counts,
or timeouts are accepted. CPU, CRISPRidentify folder parallelism, and timeouts
are operator-owned settings; folder jobs are positive and cannot exceed the
global worker CPU budget. If the operator forces bias corrections off, the
requested and effective values plus a warning are persisted and returned.

The submission response is:

```json
{
  "access_token": "<shown once; 256 random bits>",
  "token_type": "Bearer",
  "job": {"job_id": "...", "status": "queued"}
}
```

Only a SHA-256 digest of this high-entropy capability is stored. Send the token
as `Authorization: Bearer <token>` for every status, cancellation, result, and
artifact request. Losing it makes the anonymous job inaccessible.

## FASTA and scientific result contract

Input accepts multi-record IUPAC DNA FASTA or one raw DNA sequence. This is a
bounded genome/contig service, not an unlimited whole-genome upload endpoint;
the configured request, total-base, per-record, record-count, and header limits
are authoritative. Headers are normalized to safe source identifiers.
Duplicate original identifiers and identifiers that collide after sanitization are rejected
instead of silently renamed. The normalized combined FASTA, per-record FASTA
files, normalized source IDs, and normalized SHA-256 are retained.

Canonical CRISPRidentify reports must use the exact `CRISPRidentify-report`
schema name and a strict semantic version with supported major `1`. Additive
minor/patch releases (for example `1.2.0`) are accepted so the planned
CRISPRidentify update does not require a synchronized backend release; malformed
versions and unknown majors fail closed. Every report must bind its source ID,
sequence length, lowercase SHA-256 of the whitespace-free uppercase sequence,
and `sha256_scope: uppercase_sequence` to the exact normalized per-record input.

Successful exit status is not sufficient. The backend independently validates
the Integration run schema, adapter group/array counts, result group sets,
finite and domain-safe numeric fields, the LRT model enum, insertion/deletion
counts, and CRISPR-evOr likelihood/orientation consistency. Reverse predictions
use the selected reverse likelihood/tree while preserving the documented
forward-minus-reverse delta; Forward and ND select the forward tree. For
SpacerPlacer and every CRISPR-evOr tree policy, the unique leaf set must exactly
equal the file stems of that adapter group's `arrays[*].ccdb_file`. Every
non-root edge needs a finite nonnegative branch length, and the explicit branch
sum must match the CSV tree length within five-decimal serialization tolerance.
Unrelated, missing, extra, duplicate, or length-inconsistent trees are rejected.

`completed_no_eligible_groups` is a successful scientific result: detection
worked, but no repeat/overlap group had enough consistently oriented related
arrays for evolutionary reconstruction. Adapter skip counts, including
`orientation_not_determined`, are surfaced.

CRISPRidentify category is the primary classification. The API exposes the raw
model/certainty score with `model_score_is_probability: false`; it must never
be displayed as a percentage or calibrated probability.

The ZIP contains normalized input, Integration and tool logs, safe detailed
CSV/JSON/text/Newick/FASTA outputs, and versioned `provenance.json`. Provenance
records input hash/source IDs, exact backend release and wheel digests, the
exact scientific manifest digest, and all six sanitized scientific identities
without installation paths. It also records requested and effective options,
schemas and available versions, sanitized logical stage status/durations,
citations, and SHA-256 hashes for bundled scientific files. Raw stage commands
and `integration_result.json` are never published. The adapter manifest is
copied only after recursive path-field sanitization; any other path-bearing
public JSON is omitted. Pickle, HTML, symlinks, hidden files, and `work_folder`
content are excluded. Individual downloads use a narrower database-registered
allowlist, are rehashed before serving, and are attachments resolved inside the
private job root.

If safe outputs or ZIP overhead exceed the exact archive cap, the ZIP is omitted
and provenance carries an explicit warning; the writer cannot cross the cap.
Public summaries are field-allowlisted, symlink-safe, bounded while reading, and
bounded again after JSON serialization.

## Security and operations

- Exact-origin CORS only; no wildcard and no credential cookies.
- Per-IP submissions are rate limited in SQLite. Forwarded addresses are
  honored only from explicitly trusted proxy networks.
- Queue depth and restart attempts are bounded transactionally; even queued
  jobs expire when workers remain unavailable.
- The job wall-time budget begins on the first claim and remains fixed across
  worker restarts and retries; an attempt never receives a fresh timeout.
- Each job has aggregate byte and entry ceilings, bounded stdout/stderr capture,
  a minimum free-space reserve, and private per-job temporary/cache directories.
  Low storage blocks both submission and worker claims and degrades readiness.
- Jobs use independent 128-bit locators and 256-bit bearer capabilities.
- Worker commands are operator-owned argv arrays executed with no shell.
- Cancellation and timeouts terminate the Integration process group and
  descendant scientific process groups, including stages that start sessions
  or double-fork.
  If any attributed process survives SIGKILL, the job fails with a containment
  error, the worker refuses another claim, and exits so systemd can kill the
  complete service cgroup before restarting it.
- Submissions hold a private filesystem lock until their SQLite row commits.
  Cleanup quarantines and removes only old 32-hex job directories with no row
  and no live submission lock, so power-loss orphans cannot leak disk and fresh
  submissions cannot be raced.
- Strict CSP, HSTS, frame denial, MIME sniffing denial, no-referrer,
  Permissions-Policy, and no-store headers are emitted.
- Terminal jobs expire after the configured retention period.
- Submitted input is temporary. The frontend should warn users not to submit
  sensitive or personally identifying sequences.

Back up neither the transient job tree nor capability tokens. Monitor
`/api/v1/health`, filesystem capacity, worker service state, and the systemd
journal. A stale worker makes health return 503.

## Tests

```bash
cd backend
python -m pip install -e '.[test]'
pytest -q
ruff check src tests
```

The suite uses a fake Integration executable and covers FASTA provenance,
token-digest-only storage, durable queue bounds/recovery, fixed option mapping,
all three scientific outcomes, category-aware/path-sanitized summaries,
mode-specific output validation, safe bounded bundles, hash/traversal/symlink
checks, exact CORS and trusted-proxy handling, retry/retention ceilings,
authorization, strict current/future adapter compatibility, report source/hash
binding, UI-consumed numeric fields, selected-orientation likelihood mapping,
exact tree leaves and branch-length totals, persistent retry deadlines,
low-disk/log/disk-flood controls, three-way release provenance, manifest-tamper
readiness, crash-orphan reconciliation, bounded errors and summaries, and
fail-closed detached/double-fork process-tree cleanup.

RepeatTyper is integrated as the queued `repeat_type` mode, consuming repeat DNA
FASTA or the existing repeat occurrence handoff. See [repeat typing](../docs/repeat-typing.md) for immutable model identity,
uncertainty semantics, the isolated runtime and validated exports.

SpacePHARER is integrated as `spacer_association`: explicitly grouped spacer occurrences, an operator-owned targeted phage panel, native empirical calibration and protected evidence exports. See [spacer association](../docs/spacer-association.md).

## Context analysis release 2.12.0

The five context modes are advertised and admitted only when their required scientific artifacts are present in the release manifest. `locus_compare` requires clinker and CasAndra; `prophage_context` requires Phigaro, CRISPRidentify and CRISPRspacer; `insertion_context` and `defense_context` each require their own engine plus CasAndra; `rna_family` requires RNAalifold. Scientific release manifest version 10.0.0 records these independent source identities.

Locus comparison accepts 2–8 sequences of at most 50 kb each. Genomic context modes accept at most 20 records and 2 Mb total. RNA-family analysis accepts 2–50 records of 30–500 unambiguous bases, requires at least two distinct sequences, and accepts either DNA or RNA in the supplied biological orientation. These are additional bounds below the operator-wide limits. RNA-family JSON input uses the explicit `rna` sequence port; a repeat port does not assert RNA-family homology.

Complete protected exports preserve source identities, native method/model scope and missing scientific values. The compact module-result schema includes every new mode. Coordinate-only Cas-context joins retain up to 32 equally near cassettes per interval and expose total tied counts and truncation; the separate full Cas annotations remain available. Source and native-result directories, command lines and operator model paths remain excluded from public artifacts.
