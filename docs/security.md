# Security and privacy model

## Public request controls

- JSON request bodies require exactly one bounded `Content-Length`; ambiguous
  `Content-Length`/`Transfer-Encoding` framing, byte-count mismatches, and
  disconnected partial submissions fail before JSON parsing.
- FASTA accepts only IUPAC DNA, bounded records, headers, per-record bases, and
  aggregate bases. Empty, duplicate, and sanitization-colliding identifiers are
  rejected.
- Pydantic forbids unknown fields. Public options are enums, booleans, and a
  small integer range; requests cannot provide paths, commands, trees, archives,
  pickle flags, or output-reuse flags.
- Queue capacity, salted-HMAC per-source limits, and an independent global
  submission budget are persisted and enforced atomically in SQLite. Rotating
  source addresses, restarting Nginx, or reaching the loopback API does not
  reset the global admission budget. Forwarded addresses are trusted only from
  explicitly configured loopback proxies.
- Nginx adds per-IP and global request/submission budgets, dedicated job-status
  budgets, and authenticated-download concurrency/rate budgets. Unknown-host
  traffic traverses the same pre-access limits and short idle timeouts.
- CORS uses the exact GitHub Pages origin. Credentials are disabled.

## Job authorization and data

The browser receives an independent random 256-bit job capability. The server
stores only its SHA-256 digest and compares digests in constant time. A failed
job locator/token pair returns the same not-found response to avoid existence
leaks. The capability is sent to the API only as an authorization header. For
reload and cross-tab recovery, the frontend also writes it to a strict `#job=`
fragment. Fragments are not part of HTTP requests, but they persist in browser
history and are readable by same-origin scripts. The complete recovery link is
therefore a bearer credential and must be protected until expiry.

Uploads and results are research data. They are private to the bearer capability
but are not end-to-end encrypted from the service operator. The production
policy deletes terminal jobs after 72 hours. Operators should not submit
personal, clinical, controlled-access, export-controlled, or unpublished
sensitive sequences unless institutional policy explicitly permits this
service.

## Storage and overload controls

Admission requires both byte and inode headroom: the configured emergency
reserve plus the complete worst-case input or worker allocation. The worker
checks free bytes, free inodes, visible entry count, apparent size, and allocated
filesystem blocks every second. Queued cancellation atomically leaves a small
database tombstone while its private input tree is quarantined and removed
immediately.

The SQLite main database has a fixed 1 GB page ceiling and bounded WAL retention.
Release rollback refuses a database above that ceiling and refuses to copy it
onto the root filesystem unless the destination can retain the copy plus a
2 GB and 1,024-inode reserve. Release source authentication and a 25 GB /
250,000-inode persistent-volume preflight happen before candidate bytes are
copied. Nginx request-body spooling uses a private 192 MB, 4,096-inode tmpfs.

## Process isolation

Both API and worker currently run as `crisprsvc`, never root. systemd
applies a strict read-only filesystem, private devices/tmp, no new privileges,
empty Linux capability sets, namespace and kernel protections,
task/file-descriptor limits, explicit restart/log budgets, zero swap and zswap,
and memory/CPU/I/O cgroups. The worker has no network address families and
`IPAddressDeny=any`; the API is reachable only on loopback.

Tool environments and source are operator-owned and read-only. Job output is on
a dedicated persistent volume. Commands are fixed JSON arrays loaded from
root-owned configuration. Every subprocess belongs to a killable process group.

## Response and browser controls

API responses set no-store, nosniff, frame denial, no-referrer, restrictive
permissions policy, CSP, and HSTS headers. Generated content is returned only as
an attachment; HTML is excluded. Artifact paths are opaque registrations, not
request paths. The static site uses no third-party analytics or cookies. Job
capabilities are prohibited from paths and queries; the explicit recovery
fragment is the only URL-bearing exception. Synchronous in-memory latches reject
submit, example-load, cancel, and download button reentry before Vue's reactive
render cycle completes.

## Residual risks

Application checks and cgroups are not filesystem quotas. An open-but-unlinked
file is invisible to tree traversal, native code can create entries between
one-second checks, and the shared `crisprsvc` identity can write outside the
current job tree. Before treating this service as safe for hostile public input,
the de.NBI deployment **must** enforce hard byte and inode quotas on job storage
and should run each scientific job under a separate UID/mount namespace or
container with only its own directory writable. The dedicated `/srv/crispr`
mount protects the VM root filesystem, but without those quotas an exploited
scientific process can still make the service volume unavailable.

Per-IP and global limits bound work at this deployment edge but cannot guarantee
fair access against a volumetric or widely distributed denial-of-service attack.
One anonymous principal can also consume the scientific worker continuously
within the accepted queue and timeout policy. Use an upstream managed WAF with
total request/response duration controls plus authenticated institutional
CPU-hour and retained-storage quotas for sustained operation. Never expose the
loopback API to untrusted local processes.

Scientific tools process complex adversarial input and do not provide a formal
sandbox. systemd/cgroup confinement limits host impact but is not equivalent to
a dedicated VM per job. Because the API, worker, and scientific child processes
share one Unix identity, a compromised scientific process has a cross-job
filesystem blast radius within the service data root; accept only non-sensitive
research data until per-job identities and stronger filesystem isolation are
deployed.

GitHub project Pages deployments under one account share an origin. Recovery
fragments are readable by same-origin script running in the open workbench page,
so a compromised dependency or sibling-project navigation is a capability risk.
The site therefore ships no third-party script and uses a restrictive CSP, but a
dedicated custom hostname is still the preferred browser isolation boundary. Tailscale Funnel is a beta edge
transport with non-configurable bandwidth limits. For sustained institutional
operation, prefer a university-managed reverse proxy/WAF, routine OS and dependency
updates, external monitoring, and periodic restore and incident-response tests.

Operator connection records and deployment evidence are private local data.
See [publication safeguards](private-operations.md) for ignore rules, staged-blob
and outgoing-history checks. See [operations](operations.md) for recovery and
[the replacement audit](replacement-audit.md) for the scope of verification.
Historical predecessor security reviews are retained locally under
`.ops-private/predecessor-security/`; they are historical evidence, not a fresh
security certification of this release.
