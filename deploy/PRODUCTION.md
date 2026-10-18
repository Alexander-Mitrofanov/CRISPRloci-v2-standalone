# Production operations records

Deployment-specific connection details, VM names, SSH configuration, exact
runtime locations and dated rollout evidence are recorded locally in
`.ops-private/production-history.md` and `.ops-private/connections.md`.
Those records must never be committed, including to a private GitHub repository.
See [private operations and publication safeguards](../docs/private-operations.md).

The public frontend remains linked from the root README. Configure the API
origin through deployment settings; source examples use `https://analysis.example.org`.

Use [deployment and migration](README.md) and the [operations runbook](../docs/operations.md)
for generic procedures. Existing edge-maintenance guards must remain installed
until their replacement has been explicitly deployed and verified. A source
comparison does not establish the current state of a remote service.
