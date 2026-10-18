# Independent upstream dependencies

| Component | Repository | Tested contract |
| --- | --- | --- |
| Baseline v2 | https://github.com/BackofenLab/CRISPRtracrRNA | tag 2.0, preserved separately |
| CRISPRidentify v2 | https://github.com/Alexander-Mitrofanov/CRISPRidentify-v2 | CLI `--file --result_folder --preset --cas False --cpu`; report schema 1.x |
| CasAndra standalone | https://github.com/Alexander-Mitrofanov/CasAndra | CLI `predict-genome`; genome output schema 5 |

The exact revisions are in `upstream-lock.json`. Source checkout locations and runtime commands
are independent. Configure an installed executable, an upstream launcher, or a Python module
command as an argument array. V3 neither imports upstream internals nor changes upstream source.
The dependencies do not need to share Python versions or environments with v3.

CRISPRidentify's canonical `report.json` is validated for schema major, coordinate semantics,
source identity/digest, validation status, array bounds and scores, and available source
sequence reconstruction. New optional fields and schema 1.x minor additions are preserved.
Only `Bona-fide` and `Possible` arrays seed predictions. Never parse its display CSV or pickle
files for the integration.

CasAndra's `run.json`, `manifest.json`, `proteins.jsonl`, and `cassettes.jsonl` are consumed.
Input FASTA digest, output checksums, counts, contig identities and bounds are checked. Empty
Cas results remain distinct from failed computations. Unknown schema versions fail explicitly;
CasAndra currently declares an integer schema rather than a major/minor compatibility promise.

## Fetch or update

```bash
python scripts/dependencies.py sync
python scripts/dependencies.py update CRISPRidentify-v2 --ref main
# Or choose a release tag / exact commit instead of main.
python scripts/dependencies.py update CasAndra --ref main
```

An explicit update checks the configured remote and refuses local changes or a lock mismatch.
It accepts a tag, branch, or commit, including annotated tags, and records the resolved commit
and previous commit. Fetch refspecs that write local branches are rejected. The new lock is
staged before switching the checkout and published only after checkout succeeds. For caught
failures before publication, the original lock remains intact and the original branch or
detached HEAD is restored. If restoration fails, the error includes both failures and the
original revision for manual recovery. This is ordinary failure recovery, not crash atomicity
across Git and the lock file; do not run concurrent dependency updates.

For a missing checkout, `sync` stages its clone and pinned checkout before publishing the
destination. A failed clone or checkout is cleaned up so a corrected retry can proceed.
Existing directories are preserved; untracked build products in an already pinned checkout
do not prevent a no-op sync.

The maintenance commands do not edit upstream files or run network operations
during analysis. A revision change may require that upstream's installer, a native rebuild, or
reinstalling CasAndra into its own environment. Use its instructions, then check integration:

```bash
crisprtracrrna doctor --config config.toml
CRISPRTRACRRNA_TEST_CONFIG="$PWD/config.toml" pytest -q -m integration
```

A clean test run verifies real sequence processing and parsing for that configuration. To undo a
candidate update, check out `previous_commit` and restore that value as `commit` in the lock.
Do not promote a dependency revision based only on `doctor`, which checks paths/resources.

No adapter can guarantee compatibility with arbitrary future breaking CLI or schema changes.
Stable output-contract updates require no edits to CRISPRidentify or CasAndra. A breaking change
belongs in this package's adapter and tests, never in a patch to the external tool.

The local validation reused the already installed CRISPRidentify core environment and invoked
CasAndra from the freshly downloaded standalone source through its public `python -m casandra`
entry point. The fresh CRISPRidentify native library was built with its own
`scripts/build_native_matcher.sh`; all upstream tracked files remained clean.
