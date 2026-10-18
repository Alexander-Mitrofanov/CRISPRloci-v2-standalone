# Verified replacement — 27 September 2026

Historical evidence for the release below. For the current source comparison,
regressions and publication safeguards, see the
[28 September replacement audit](replacement-audit.md).

The installed CRISPRloci v2 backend and Integration executable sources match the
workspace. The selected immutable release is
`e54e11421cf705f0c2c733274a0dcbcfa38443c71df3a3db8e111cb96a8e0b08`.
All five source repositories remain clean at the commits in `tools.lock.json`.

- 357 Python tests passed: inherited backend/Integration behavior, new annotation
  validation and API options, schema-2 migration with a retained completed job,
  capabilities/download records, full exports beyond preview limits, update
  activation locks and runtime integrity.
- 88 frontend tests passed; formatting, lint, source-boundary checks, production
  build and public-export scanning passed.
- Seven real HTTP → queue → worker → scientific tool → authenticated archive
  cases passed on the final release, covering all six modes and both RNA model
  families. Cas/locus examples produced 5 Cas calls; the complete locus example
  produced 97 RNA candidates. Type II and V-K model examples produced 4 and 1.
- Chromium passed actual locus submission, annotation/provenance rendering,
  authenticated ZIP download, recovery after reload and 390-pixel mobile layout.
  The compiled production frontend also passed the existing evolutionary
  example and a real locus job with no CSP violations or JavaScript exceptions.
  Its test HTTPS API origin was forwarded locally by Playwright; no public
  deployment was contacted or changed.
- The standalone audit checked 160 runtime code/configuration files and 7,514
  runtime/tool symlinks without finding references to either predecessor folder.
  Historical source-provenance documents intentionally retain their lineage.
- An actual GitHub fetch resolved evOr's `main` into a candidate lock without
  changing the active runtime. The previous gated application release is kept
  for rollback; activation rejects a running managed server.

Machine-readable evidence is in [verification.json](verification.json).
Detailed local gate jobs, logs and browser screenshots remain under `.runtime/`.
The scripts `gate_release.py`, `browser_acceptance.py`, `browser_production.py`
and `audit_standalone.py` reproduce the corresponding checks. Browser scripts
use the bundled tool examples and require the local server; the production
browser check additionally requires a preview on port 4173 built with the test
origin `https://analysis.example.org` and backend CORS allowing that preview.

The local application runs entirely from this project. Production services,
public Pages, TLS configuration and existing production job data were not
changed. Follow [deployment and migration](../deploy/README.md) before retiring
any predecessor folder still referenced by a separate live deployment.
