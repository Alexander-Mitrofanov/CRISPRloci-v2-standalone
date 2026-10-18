# Deployment, updates and retirement

This project runs independently of both predecessor source trees. No running
service, DNS record, TLS edge or public Pages repository is changed by the setup
scripts. Build and test a new installation before redirecting traffic.

## Install at the final path

The example systemd units use `/srv/crisprloci-v2`, the existing unprivileged
`crisprsvc` service account and `/srv/crispr/jobs/web` for persistent data.
Provision at that final path; do not move Conda/venv environments after building.
For another path, edit all three unit WorkingDirectory/ExecStart entries.

Keep SSH/VM details only in `.ops-private/`; see
[private operations safeguards](../docs/private-operations.md). Exclude
`.ops-private/`, `output/` and `spacer analisis/` from every source transfer.

Copy the project source, tool lock and independent Git checkouts (or clone them
with `scripts/sync_tools.py`). Exclude `.runtime/`, `.runtime-data/`, local `.env`
files, `backend/.venv`, `node_modules`, caches and build output. Follow the root
README's setup, provision, build and real release gate. Select the passing release
with `manage_release.py activate` before starting services.

The production service account must be able to read the new project and runtime.
The project, tools, runtime snapshots, launch scripts and release pointers should
be owned by an operator account and unwritable by `crisprsvc`. The service account
needs write permission only to the persistent data root. Production units enforce
read-only scientific/application code with `ProtectSystem=strict` and retain
process, memory, CPU, disk, temporary-directory and network containment settings.
Local development has the normal permissions of the current user.

Copy `backend/.env.example` to `/etc/crisprloci.env`, owner root, mode 0600.
Set these deployment-specific fields explicitly:

- `CRISPR_WEB_DATA_ROOT=/srv/crispr/jobs/web`
- `CRISPR_WEB_DATABASE=/srv/crispr/jobs/web/queue.sqlite3`
- `CRISPR_WEB_CORS_ORIGINS=https://YOUR-EXACT-FRONTEND-ORIGIN`
- Optionally, `CRISPR_WEB_CATALOG_DATABASE=/srv/crispr/catalogs/SNAPSHOT/catalog.sqlite3`
  for the read-only Database explorer. Copy the complete verified snapshot,
  including its `shards/` directory, to operator-owned storage readable by
  `crisprsvc`. Keep it unwritable by the service account and outside home
  directories (`ProtectHome=yes`). See [catalogue operation](../docs/catalogue.md).
- `CRISPR_WEB_TRUSTED_PROXY_CIDRS=127.0.0.1/32` only for the supplied local proxy
- CPU, input, queue, retention and storage limits appropriate for this machine.

Remove old runner commands, release hashes and scientific-manifest paths from
imported environment files: `serve.py` derives them from the selected release.
Do not copy a predecessor `.env` wholesale. The reference defaults retain the
original admission limits: 2 million bases, 20 records, 8 active jobs, 2 active
jobs per client, 24-hour retention and a 2 GB disk budget per job.

Install `deploy/systemd/crisprloci-{api,worker,cleanup}.service` and the cleanup
timer in `/etc/systemd/system/`. Review resource limits against available RAM,
CPUs and mount points. Run `systemctl daemon-reload`, then start the new API and
worker after the migration steps below. Enable the cleanup timer. Workers have
no network access under systemd; updating and building tools is an operator task.

The supplied Nginx template retains the previous local TLS-edge arrangement:
PROXY protocol on `127.0.0.1:8080`, a configured `__PUBLIC_HOST__`, forwarded
headers overwritten, body/rate/connection limits, and API on port 8000.
It requires a trusted local TLS edge sending PROXY protocol. For a different
TLS termination arrangement, adapt that listener and client-address handling
before use. Never expose the API's loopback endpoint directly.

## Preserve existing jobs

1. Disable new submissions at the edge and allow current jobs to finish, or
   cancel them through their authorized API. Stop the old API, worker and cleanup
   timer. Never run an old and new worker against the same queue concurrently.
2. Preserve the entire persistent data root, including `jobs/` and the SQLite
   database. Also preserve any existing frontend recovery links and exact API
   origin if old links must keep working. Tool checkouts and source folders do
   not contain the durable queue unless explicitly configured that way.
3. Back up the stopped database with the included SQLite backup helper. Its
   destination directory must be root:root, mode 0700, with enough free space:

   ```bash
   sudo python3 deploy/backup-database.py backup \
     --source /srv/crispr/jobs/web/queue.sqlite3 \
     --destination /absolute/private-backup/queue-before-v2.sqlite3 \
     --release before-crisprloci-v2
   ```

   Back up the matching job directories as well; a database-only backup cannot
   restore downloads. Keep backup manifests alongside their databases.
4. Start the v2 API against the existing data root, wait for schema initialization
   and fresh health, then start the worker. Startup migrates queue
   schema 2 through the current schema 7 in one transaction, preserving job IDs, capability hashes,
   artifacts, options, timestamps and rate-limit salt. The migration has a
   regression test with a completed job and registered downloads.
5. Check `/api/v1/health/fresh` locally and `/api/v1/version`; verify a saved job
   recovery link and download, then submit one new job and restore edge traffic.
   The active backend and scientific release IDs must match the tested release.

Schema 7 retains per-tool CLI overrides. Older backends cannot honor these overrides; rollback to a release without schema-7 support requires its matching data snapshot. Rolling back to the predecessor
requires stopping v2 and restoring the pre-migration database and matching job
snapshot; jobs created after that snapshot are not present in it. Rollback among
compatible CRISPRloci v2 releases uses `manage_release.py rollback`, with services
stopped. Only use that rollback with a schema-compatible release; otherwise
restore the matching data/configuration snapshot.

## Frontend

Build `frontend/` with Node 24 and an exact HTTPS `VITE_API_BASE_URL`.
Use `VITE_BASE_PATH=/` for a root site or `/<repository-name>/` for Pages.
Only the frontend directory is public-exportable. The nested Pages workflow
uses the repository variable `VITE_API_BASE_URL` and derives its own base path.

`deploy/export-pages-frontend.sh /absolute/clean-frontend-checkout
https://github.com/OWNER/FRONTEND.git` exports only frontend files into an
explicit, separate, clean repository. It does not commit or push. Review and
publish that diff using the intended frontend repository. No production URL is
hard-coded in the application metadata.

## Tool updates

Use the root README's candidate → provision → build → gate → stop → activate
workflow. Running jobs keep their immutable toolset. The API cannot fetch GitHub
code, choose executables or update a release. GitHub credentials belong to the
operator and never go in `VITE_*` or job data. Installed-file inventories include
Python libraries, native executables and bundled scientific models; verification
refuses edits to a staged runtime. Keep previous releases for rollback.

`provision_tools.py --resume` resumes an interrupted, not-yet-sealed toolset.
It never repairs a sealed active runtime. Create a new release when code or
scientific dependencies change. Exact Conda package URLs and pip inventories
are retained in each toolset as rebuild/diagnostic evidence; platform packages
must still be available when provisioning a new machine.

## Before deleting a predecessor folder

- The selected release verifies and its real seventeen-case gate passes entirely
  from this project. No service ExecStart/WorkingDirectory, environment file,
  symlink, tool launcher or frontend export command points to a predecessor.
- Any persistent data stored inside that old folder has been copied and tested
  through v2 recovery/downloads, or intentionally retired after retention.
- Production units and the published frontend have been switched and verified
  if that folder still serves production. Keep independent source/data backups.

The development replacement can be used immediately without either predecessor
folder. Deleting source folders does not itself switch an existing live service.

The production units set `ProtectHome=yes`: install their base Python and Conda
outside home directories (for example system Python 3.12 and `/opt/conda`).
The local development runtime may use the current user's Python installation.

See [Protospacer reference provisioning and activation](../docs/protospacer-integration.md) for the new queued workflow,
immutable reference identity and three-record release gate.

### Bidirectional spacer searches

CRISPRloci backend 2.3.0 / Integration 0.7.0 expose two independently advertised
queued modes in **Spacer searches**. `protospacer` searches submitted linear DNA
against the pinned NCBI-hosted 2017 spacer collection with full-length exact
matching on both strands. `viral_search` searches submitted spacer FASTA against
a pinned viral reference with CRISPRspacer 0.3.0 and BLAST+. The latter accepts
`viral_max_mismatches` 0, 1 or 2 (default 2), retains full-query ungapped alignments,
and reports ineligible short/ambiguous spacers separately. Its 80 nt per-record
limit and the deployment's ordinary record/request limits apply. Viral retrieval
is heuristic; sequence matches do not establish infection or host range.

Operators configure `CRISPR_WEB_SPACER_DATABASE` and
`CRISPR_WEB_SPACER_DATABASE_SHA256` (SQLite bytes), plus
`CRISPR_WEB_VIRAL_REFERENCE` and `CRISPR_WEB_VIRAL_REFERENCE_SHA256`
(`reference.json` bytes). Both references are read-only, outside public exports,
and prepared before starting services. Queue schema 6 retains reference digests
and mismatch policy per job. A worker refuses a changed reference. Native BLAST
is part of the inventoried tracrRNA environment and runs only in the isolated
worker. No API request downloads data or changes a reference.

Use the standalone CRISPRspacer repository's `docs/viral-reference-methods.md`
and `docs/viral-refseq-237-snapshot.json` to reproduce the full viral reference
on the server. Local and release tests build only the bundled three-record
samples. `scripts/gate_release.py` exercises both directions, result binding,
private-file exclusion and complete JSON downloads before activation. For a
schema upgrade, drain jobs and snapshot the full queue/job directory before
activation; rollback requires that matching snapshot as well as the old release.
