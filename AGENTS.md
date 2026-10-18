# CRISPRloci v2

This repository owns the API, durable worker, integration adapters, frontend and
runtime release management. Scientific sources live in independent GitHub
checkouts under tools/; their revisions are pinned in tools/tools.lock.json.

Read the nearest nested AGENTS.md when working in backend, Integration or frontend.

- The API queues work; only isolated workers invoke scientific tools.
- Keep backend output validation independent from Integration's validation.
- Preserve atomic queue operations, per-job capabilities, bounded reads,
  cancellation, retention and process containment.
- Missing scientific values remain missing; uncertainty remains explicit.
- SpacerPlacer and evOr must use separate environments because both own sp_model.
- Never update a tool checkout or runtime from an anonymous request or running job.
- Public deployment exports only frontend/. Never publish backend/config/job data.
- Tests and runtime setup must not depend on predecessor project directories.
- Run scripts/check.sh for relevant changes, including browser tests for UI flows.
- Do not operate production services without an explicit deployment request.
- Keep SSH/VM identifiers and connection journals only in ignored `.ops-private/`;
  record steps there without credentials. Run `scripts/check_private_data.py`
  before publication and install `.githooks` in every real checkout. Never copy
  `output/`, `.ops-private/` or local deployment evidence into a GitHub repository.
