# CRISPRloci v2 replacement

1. Preserve the maintained server's durable queue, capability authorization,
   cancellation, admission limits, independent scientific validation, artifact
   policy, recovery flow, and evolution visualization.
2. Add independently configured CasAndra and CRISPRtracrRNA stages, including a
   complete locus workflow. Validate and present real machine-readable results.
3. Build immutable runtimes from the five GitHub commits in tools/tools.lock.json.
   Stage updates separately; verify provenance and compatibility before activation;
   retain the previous release for rollback. Never pull tools inside a job.
4. Rebrand the frontend and expose each implemented workflow and its result views.
5. Provide standalone setup, local launch, deployment, migration, backup and update
   documentation; eliminate dependencies on either predecessor source folder.
6. Verify regression suites, real scientific jobs, browser submission/results,
   fresh runtime setup, and deletion readiness without removing the user's folders.

Production activation is a separate operator action. This implementation prepares
and tests its release commands without changing the running de.NBI service.
