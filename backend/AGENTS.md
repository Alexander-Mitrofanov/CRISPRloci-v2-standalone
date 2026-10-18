# Backend agent guide

The backend is an installable `src/`-layout package. Public compatibility
entry points are `crispr_web.api:create_app`, `crispr_web.worker:main`, and
`crispr_web.cleanup:main`.

## Ownership

- `http/middleware.py`: ASGI security, framing, and bounded-error middleware.
- `http/routes/service.py`: health, config, and version routes.
- `http/routes/jobs.py`: submission, capability-protected lifecycle, downloads.
- `api.py`: composition, readiness cache, dependencies, compatibility exports.
- `db.py`: atomic persistence and state transitions.
- `config.py`: stable settings façade; environment parsing and release-manifest
  identity checks remain together so startup validation has one fail-closed path.
- `service.py`: job-level policy over persistence and private storage.
- `worker.py`: orchestration, process containment, and compatibility façade;
  focused helpers are mapped in `src/crispr_web/execution/README.md`.
- `summary.py`: stable façade over the focused modules documented in
  `src/crispr_web/summaries/README.md`.
- `artifacts.py`: stable façade over the independent publication modules
  documented in `src/crispr_web/publication/README.md`.

## Rules

- Do not import Integration validators into backend validation; duplication is
  intentional defense in depth.
- Do not split a transaction across repositories or helper modules.
- Keep paths, commands, tokens, and private manifests out of public responses.
- Preserve bounded reads, symlink rejection, process-group containment,
  cancellation checkpoints, deadlines, and storage accounting.
- `db.py` is intentionally transactionally cohesive. `Worker._execute` keeps
  checkpoint, launch, containment, validation, and finalization order visible;
  extract only leaf helpers without distributing ownership of that sequence.
- Keep current façade imports working while internals move.
- Runtime dependencies must stay compatible with
  `requirements-production.lock`. Test dependencies live in the `test` extra.

## Checks

- HTTP/API: `cd backend && .venv/bin/pytest -q tests/test_api.py tests/test_api_hardening.py tests/test_admission_readiness.py`
- Worker: `cd backend && .venv/bin/pytest -q tests/test_worker.py tests/test_worker_hardening.py tests/test_integration_contract.py`
- Publication: `cd backend && .venv/bin/pytest -q tests/test_output_hardening.py tests/test_resource_exhaustion_regressions.py`
- Backend gate: `scripts/check.sh backend`

FastAPI TestClient uses an asyncio thread bridge. In restricted sandboxes that
block cross-thread loop wakeups, run the same test command in the approved local
execution environment rather than weakening the tests.
