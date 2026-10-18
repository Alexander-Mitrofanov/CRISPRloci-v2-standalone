# Component connectivity audit

Local source audit, 11 October 2026. No production deployment or runtime update
was performed. The [sequence contract](module-exchange.md) defines biological
roles; [compact results](module-results.md) describe findings and are not general
sequence inputs.

## Coverage

| Component | Checks |
| --- | --- |
| Frontend and API contract | All 15 modes, advertised references, request fields, upload previews, input bounds, response normalization, source boundaries, production build and public export scan |
| API and queue | Durable admission, capacity, capabilities, CORS/security policy, cancellation, retention, isolation, concurrent state transitions and recovery |
| Worker and Integration | Input hashes, stage commands/options, process containment, deadlines, independent scientific validation, bounded publication and compact results |
| Results and subsequent jobs | Actual publisher output through consumer admission, exact sequence/identifier preservation, source grouping, incompatible-role rejection and an HTTP download → submission → worker → download round trip |
| Browser | JSON upload/download and direct continuation, old recovery preserved until submission, incompatible/incomplete input rejection, desktop/mobile layout and JavaScript errors |
| Runtime/release support | Existing release management, reference preparation, deployment configuration, private-data and systemd resource-policy regression tests |

## All module inputs and reusable outputs

Every successful module supplies its original input through `input-sequences.json`.
This permits rerunning compatible analyses with different options without treating
findings as observed sequences. Each destination still enforces its own limits.

| Module | Input | Additional reusable output |
| --- | --- | --- |
| `detection` | Target DNA | Arrays, repeat occurrences, spacer occurrences |
| `orientation` | Target DNA | Arrays, repeat occurrences, spacer occurrences |
| `reconstruction` | Target DNA | Arrays, repeat occurrences, spacer occurrences |
| `cas` | Target DNA | Native/compact Cas findings; no current Cas-result consumer |
| `tracrrna` | Target DNA | Candidate findings and native FASTA; candidates are not repeat inputs |
| `loci` | Target DNA | Arrays, repeat occurrences, spacer occurrences; combined locus findings |
| `leader` | Target DNA | Arrays, repeat occurrences, spacer occurrences; adjacent-window FASTA and coordinates |
| `repeat_context` | Target DNA | Arrays, repeat occurrences, spacer occurrences; contextual folds |
| `protospacer` | Target DNA | Exact reference-match evidence; original target DNA remains reusable |
| `repeats` | Repeat DNA or RNA | Original repeat port; predicted folds remain findings |
| `repeat_map` | Repeat DNA | Original repeat port; reference annotations remain findings |
| `repeat_type` | Repeat DNA | Original repeat port; subtype predictions remain findings |
| `viral_search` | Spacer DNA | Original spacer port; match evidence remains findings |
| `spacer_association` | Spacer DNA and explicit grouping | Original spacer port with any supplied origins |
| `array_compare` | Structured arrays | Original arrays and all observed spacers, including supplied source origins |

Repeat ports connect to structure, reference mapping and subtype evidence.
Spacer ports connect to viral search and explicit group association. Array ports
connect to observed comparison. Target ports connect to the nine modes accepting
genomic DNA. RNA ports connect only to RNA repeat structure analysis. Missing,
incomplete, ambiguous or oversized inputs are rejected according to the relevant
consumer; the handoff never silently trims records or substitutes consensus DNA.

Evolutionary analysis still requires its existing comparable-array eligibility
checks. A spacer list or a comparison network does not constitute evolutionary
input. Transcript candidates, inferred histories, model labels and folding
scores are not interchangeable with observed repeat/spacer occurrences.

## Repairs

- Made FASTA identifier normalization idempotent. Canonical identifiers ending
  in `_` or `-`, including truncated headers, retain their IDs across worker
  verification and subsequent submissions.
- Restored comparison spacer origins and added `spacers.json`, so comparison
  results can feed source-grouped association without losing occurrences,
  array ordinals or source identity. Deletion slots remain explicit in arrays.
- Aligned the sequence-origin source-ID limit with the array contract.
- Added a bounded, authenticated browser handoff into a fresh destination form.
  The new job is submitted only through its normal form and receives a separate
  capability. Failed, empty and incomplete outputs cannot start a handoff.
- Corrected source grouping for BOM-prefixed JSON and the frontend limits for
  target searches and repeat-context cohorts.
- Corrected stage display: CCTK follows validation only in observed comparison;
  unrelated analyses no longer display that stage.
- Updated the UI guidance to include repeat subtype, spacer association and
  observed-array comparison consumers.

## Verification and limits

`scripts/check.sh` passed with **1,015 backend/Integration/repository tests** and
**298 frontend tests**, plus lint, source boundaries, production build and public
export scanning. The added backend matrix exercises 96 input-port connections,
36 detection-output connections, all mode/role acceptance pairs, RNA boundaries,
source preservation and an authenticated HTTP/worker round trip.

The default run skipped 67 executable-dependent tests. A separate run using the
existing BLAST tools and the current CRISPRspacer source passed all 39 viral and
protospacer tests, including the opt-in native lifecycle check. CCTK, RepeatTyper
and SpacePHARER executables are unavailable in this workspace; their fresh native
runs are not certified by this audit. Their admission, fixtures, output contracts
and frontend coverage remain exercised. Scientific accuracy is not established
by software contract tests.

The headless Chromium handoff check passed against a local production build with
an intercepted fixture API. A second browser run against a real temporary local
API and durable worker passed submission, authenticated download, reload recovery,
return-to-results navigation, service-outage recovery, mobile layout and production
CSP checks. These runs use fixture scientific execution and are not a production
deployment health check. Browser evidence is local-only under the system temporary
directory.

Git metadata is inaccessible in this source workspace. The guard installer
reported source-archive mode; install `.githooks` in the real checkout before
publication. The private-data scanner passed. No public publication was made.
