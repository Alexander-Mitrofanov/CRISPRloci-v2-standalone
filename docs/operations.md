# Operations runbook

Connection details and dated execution records belong exclusively in the local
`.ops-private/` directory. Follow [private operations](private-operations.md).
These procedures require an explicit deployment request before changing a live
service. The development checks never establish current remote service health.

## Release and verify

1. Run `bash scripts/check.sh`, including the privacy and compatibility checks.
   Use the browser scripts for changed workflows with disposable local data.
2. Prepare pinned scientific sources and independent environments using the
   root README's candidate, provision, build and real-science gate sequence.
   Keep SpacerPlacer and evOr isolated. Use the gate for the exact candidate;
   earlier reports do not approve changed sources or environments.
3. Before a production switch, close submissions, drain work, stop API, worker
   and cleanup, and back up the complete queue/job directory and private
   configuration together. Authenticate backups with `deploy/backup-database.py`.
   Preserve the selected and previous releases and the corresponding schema.
4. Follow [deployment and migration](../deploy/README.md). Start API first and
   wait for database initialization and fresh health before starting the worker.
   Verify the installed package and scientific identities, worker readiness,
   SQLite integrity, retained job recovery/download and a new scientific job.
5. Verify effective systemd/cgroup limits, loopback-only API, exact CORS origin,
   unauthorized negatives, cancellation, restart recovery and retention. Keep
   ingress closed if rollback or health verification fails.
6. Export only the frontend with `deploy/export-pages-frontend.sh`, review the
   diff, run public scanning and verify the browser against the approved origin.
   Record versions, checksums, outcome and evidence paths privately.

## Routine checks and incident recovery

- Check API/worker/cleanup unit status, cleanup timer, fresh health, queue length,
  stale leases, failures, disk bytes/inodes, SQLite/WAL sizes and oldest job age.
- Keep OS updates, dependency review, encrypted backups and restore drills on an
  operator schedule. Record each completed check without input sequences or
  bearer tokens. This document does not create a monitor or scheduled task.
- On low storage, stop new submissions, terminate or drain work, run supported
  retention cleanup and inspect orphaned jobs. Never recursively delete the
  queue or job root. Retain emergency reserves and enforce filesystem quotas.
- To roll back an incompatible schema, stop all writers and restore the matching
  database, job files and configuration together. Preserve later jobs separately;
  an older snapshot cannot contain them. Never start old and new workers together.
- Existing predecessor edge transaction/redaction guards may still protect the
  deployed host. Keep their installed files and state until a separately approved
  replacement has passed failure/recovery tests. Do not copy or bypass an old
  installer against v2's schema or delete a guard to reopen traffic.

The v2 release manager owns immutable local releases; it is not a replacement for
every predecessor host bootstrap/edge transaction command. Host installation and
those retained controls need their own deployment verification. The residual
shared-UID, quota and edge risks remain in [the security model](security.md).
