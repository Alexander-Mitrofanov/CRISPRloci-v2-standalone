# Initial standalone verification

Recorded 28 September 2026 for CRISPRspacer 0.1.0.

## Automated tests

`bash scripts/check.sh` completed with **79 tests passed**, plus a passing Ruff
check. Python 3.12.7 was used locally. CI also specifies Python 3.10 and 3.13;
their status is available in the repository's workflow runs rather than assumed
from the local result.

The tests cover:

- 2,100 seeded cases comparing exact/one/two-substitution discovery with an
  independent exhaustive window enumerator, across linear and circular targets.
- Forward/reverse and palindromic matches; duplicate spacer occurrences.
- Source-array and other-category array exclusions; unknown annotation coverage.
- Explicit assembly membership and unknown relationships across records.
- CRISPRidentify source hashes, repeat/spacer geometry and deletion ordinal gaps.
- Ambiguous sequence handling, no-match/no-query/no-window outcomes.
- Malformed FASTA/JSON, duplicate IDs, nonfinite/oversized numeric tokens and symlinks.
- Candidate, hit, pair, work, evidence and deadline limits, including zero-window cases.
- Atomic output failure/overwrite refusal, altered manifests and rehashed result tampering.
- Reference identity binding, source snapshot consistency and complete deterministic replay.

## Package independence

Both a wheel and source archive were built successfully using the standard
Python build frontend. The wheel was installed with `--no-index --no-deps` in a
new virtual environment. Its import location was verified inside that
environment's `site-packages`.

Nine CLI invocations then ran from an unrelated temporary directory with no
`PYTHONPATH`: version reporting; direct search and validation; same-source
analysis and validation; reference creation and validation; reference-backed
substitution search and validation. Only copied synthetic fixtures were used.
No CRISPRloci, Integration, scientific-tool checkout or source-directory import
was required at runtime.

## Synthetic benchmark

[Machine-readable measurements](benchmark-synthetic-v1.json) come from
`python scripts/benchmark.py`: seed 20260928, one 100,000-base random ACGT
target and 100 observed 32-base query substrings. Each run recovered 100
alignment occurrences. Measured elapsed times were:

| Target traversal | Substitutions | Seconds |
| --- | ---: | ---: |
| Linear | 0 | 0.0422 |
| Linear | 1 | 0.0652 |
| Linear | 2 | 0.0888 |
| Circular | 0 | 0.0424 |
| Circular | 1 | 0.0657 |
| Circular | 2 | 0.0885 |

These timings describe a small synthetic workload on this machine. They are
not estimates of catalogue-scale capacity, biological sensitivity, specificity
or host-prediction accuracy. Larger/repetitive collections may reach declared
budgets; those fail explicitly rather than yielding a silently partial result.

## Council review

The scientific reviewer found the search, coordinate and uncertainty semantics
consistent with the declared scope. The integration reviewer identified resource
accounting, reference binding, single-snapshot input, JSON numeric parsing and
circular-copying issues. Each was fixed with regression coverage. A final
focused review found no material blockers among those findings and independently
reran the 79-test suite successfully.

No biological input analysis, external reference download, v2 integration or
production deployment was performed as part of the initial 0.1.0 release.

## Version 0.2.0 database checks

The suite now includes indexed-search versus exhaustive-oracle comparisons on
160 generated targets and 92 spacer occurrences, original three-record sample
CLI round trips, both strands, palindromes, duplicate sequences, ambiguous tails,
separate records, eligibility counts, failure atomicity, pinned hash rejection
and five rehashed-but-semantically-invalid database mutations. Independent
review led to complete seed-inventory validation during import and explicit
SQLite deadline interruption. Full downloaded data and VM measurements are
kept on the VM; development retains only the three-record source sample.

## Version 0.3.1 review regressions

The standalone gate passes **225 tests** and Ruff with local BLAST 2.12.0.
The patch adds corrupted-deflate CLI/reader checks, rehashed boolean-versus-number
evidence mutations, exact-target identity/topology checks, bounded metadata
materialization, rejection of recursive metadata views, cancellation/deadline
interruption during catalog scans and evidence SQL, and prior-version replay.
All result schemas and scientific scoring defaults remain unchanged.
The viral/reference/CLI and review regression suites also pass **109 tests**
using the separate pinned toolset's BLAST 2.17.0 executables.
An offline-built 0.3.1 wheel installed in a fresh virtual environment also passes
the indexed three-record search/validation and generated viral reference/search
roundtrips at both zero and two substitutions, from an unrelated temporary
directory with no `PYTHONPATH`, using BLAST 2.17.0.

## Version 0.3.0 viral-reference checks

The full standalone gate passes **198 tests** and Ruff locally. New tests use
generated DNA and mocked network replies; no viral corpus is downloaded by CI
or local tests. Actual BLAST 2.12.0 builds and searches run against these small
fixtures. The reference builder also passed a full roundtrip with BLAST 2.17.0.

Coverage includes release-plan consistency, checksum and gzip failures, fixed
source URLs, source receipt binding, complete native-file inventories, full
sequence roundtrips, both search strands, repeated locations, duplicate spacer
occurrences, 0/1/2/3 substitution boundaries, clipped terminal mismatches,
ambiguous queries, empty outcomes, process-group termination, resource caps,
atomic failure and rehashed evidence/provenance tampering. Metadata acquisition
was also tested through import and native validation as one workflow.

The built 0.3.0 wheel was installed with `--no-index --no-deps` in a new virtual
environment. From an unrelated temporary directory with no `PYTHONPATH`, its
public CLI built and validated a synthetic reference, searched three generated
spacers, then validated the retained HSP evidence. It required no CRISPRloci
installation. These tests establish implementation consistency, not biological
sensitivity, specificity or host-prediction accuracy.

The separate VM verification used the complete viral RefSeq release 237:
19,625 records and 583,708,944 bases. Source checksums, complete sequence
roundtrip and `blastdbcheck -full` passed. Exact and two-substitution CLI modes
both recovered seven known interval matches (forward/reverse plus a duplicate
occurrence), with two ineligible queries retained as skipped. A further run as
the service account with only a loopback network interface passed against the
read-only reference. Its [metadata receipt](viral-refseq-237-snapshot.json)
contains hashes and counts; source-derived test sequences remain on the VM.
The optional native checker needed separate NCBI taxonomy lookup files, now
documented in the reproduction method. GitHub CI passed on Python 3.10, 3.12
and 3.13 for the implementation commit. No production services were changed.
