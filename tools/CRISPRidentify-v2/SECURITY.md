# Security policy

## Supported versions

Security fixes are applied to the latest tagged release and, when appropriate,
to the current development branch. Older tags are not guaranteed to receive
backports.

| Version | Security support |
| --- | --- |
| Latest tagged release | Supported |
| Current development branch | Best effort before release |
| Older releases | No guaranteed backports |

## Report a vulnerability privately

Do not disclose vulnerability details in an ordinary issue, pull request,
discussion, benchmark log, or FASTA attachment.

1. If GitHub shows **Security → Advisories → Report a vulnerability** for this
   repository, use that private form.
2. If private vulnerability reporting is unavailable, contact a repository
   administrator through an already established private channel.
3. If no private channel exists, open an issue that contains only a request for
   private security contact. Do not include technical details in that issue.

Include the affected tag or commit, threat model, reproduction steps, impact,
and any proposed mitigation. Remove credentials and confidential sequence data.
Use a minimal synthetic input whenever possible.

Maintainers will acknowledge the report through the same private channel,
validate the scope, coordinate a fix and release, and agree on disclosure
timing with the reporter. Response time depends on severity and maintainer
availability; this policy does not promise a fixed service-level agreement.

## Security scope

Examples that belong in a private report include:

- command or argument injection;
- unsafe archive, path, symlink, or temporary-file handling;
- unbounded resource use that permits a practical denial of service;
- bypass of model, manifest, checksum, or report-provenance validation;
- dependency or GitHub Actions supply-chain compromise;
- unintended disclosure through logs, reports, or temporary files;
- processes surviving a documented timeout or cancellation boundary.

Unexpected biological predictions, score disagreements, and benchmark
regressions normally belong in the scientific-result or performance issue
forms unless they also create a security impact.

## Coordinated disclosure

Please allow time to reproduce the issue, prepare tests, update affected
artifacts, and notify authorized users before publishing details. The project
will credit reporters who request credit and will not name reporters who prefer
anonymity.

No security policy can make untrusted sequence analysis risk-free. Run large or
untrusted inputs with operating-system resource limits and least-privilege file
access, and keep production credentials out of the analysis environment.
