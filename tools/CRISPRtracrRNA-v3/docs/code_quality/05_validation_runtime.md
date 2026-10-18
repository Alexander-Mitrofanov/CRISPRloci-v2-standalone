# Independent validation: runtime and publication

Reviewer: `validation_runtime`, a new validation-council member who did not plan or
implement the changes. Review completed against the frozen source in
`/tmp/tracrrna-quality-baseline` on 2026-09-27.

## Decision

**Accept the runtime and reporting changes. No blocking or actionable regression
was found.** The implementation follows the two relevant consensus plans and
addresses their reproduced defects without broadening scientific behavior.
Allowing empty non-executable argv elements is the correct compatibility
adjustment: an empty argument is valid subprocess input.

## Independent review

Compared the baseline and current `config.py`, `models.py`, `sequence.py`,
`process.py`, and `pipeline.py`; inspected the extracted `reporting.py`, CLI
behavior, and focused regression tests.

- Numeric, coordinate, path-shape, command-shape, and strand validation changes
  intentionally reject malformed values at their public boundaries. Valid default
  values, override behavior, config-relative data paths, executable symlink
  identity, and environment-value redaction are retained.
- Executable lookup interprets both explicit relative paths and relative/empty
  `PATH` components in the effective child directory. The absolute path passed to
  `Popen` is also the one hashed and recorded. Timeout cancellation policy is
  unchanged; stderr diagnostics now read a bounded UTF-8 tail.
- JSON serialization keeps its prior bytes and finite-number requirement.
  Exclusive unique sibling creation prevents publication attempts from sharing a
  temporary file. Normal umask permissions, same-directory replacement, and
  cleanup preserve the intended single-file atomicity contract.
- Export ID construction is centralized and applied before retained-record
  serialization. The exported IDs, prediction order, columns, quoting, and
  escaping are preserved. Accepted rescue proposal aliases receive matching IDs.
- `run` now counts records after their publication succeeds; setup failures enter
  the same failure boundary. Secondary failure-marker errors cannot replace the
  primary exception. Existing output refusal remains before writes, and final
  report publication remains after the exports and manifest.
- The explicit pipeline re-exports and calls through its module globals preserve
  existing imports and monkeypatch seams. The extraction introduces no dependency
  or report-schema change.

## Independent execution

The additional probe is retained at
`/tmp/tracrrna-quality-runtime-validation.py`; its machine-readable result is
`/tmp/tracrrna-quality-runtime-validation.json`.

Command from the repository root:

```bash
PYTHONPATH=new_version/src python3 /tmp/tracrrna-quality-runtime-validation.py
```

All five checks passed:

1. A mixed absolute/relative/empty `PATH`, including a nonexecutable shadow entry,
   was used through a symlinked child directory containing `..`. Direct native
   `subprocess.run` and `Runner` selected the same executable; the recorded
   absolute symlink-preserving path and SHA-256 matched the launched file.
2. A real child deliberately ignored SIGTERM. The configured timeout still
   stopped it within the expected bounded interval, raised `ToolError`, and
   persisted a timeout record equal to the in-memory record.
3. Eight concurrent writers performed 80 JSON replacements. Every observed
   destination parsed completely and contained a self-consistent writer payload;
   attempts cleaned their temporary files and preserved an unrelated pre-existing
   fixed-name `.tmp` sibling.
4. A second retained-record write was failed deliberately, both with working and
   failing failure-marker storage. Exactly one record remained published; the
   successful failure marker counted one; neither run published `report.json`;
   both raised the original exception object. The secondary storage failure was
   logged as expected.
5. A TOML boolean timeout was rejected with a setting-specific `TracrError`, while
   a valid CLI-style numeric override still took precedence.

Also independently ran:

```bash
cd new_version
python3 -m pytest -q tests/test_runtime_contracts.py -k 'cli_'
```

Result: **2 passed, 129 deselected**. These cover the legacy direct-flag entry,
aliases, explicit booleans, and the doctor path. Broader tests were not needlessly
repeated: the lead's integrated unit/native, lint, formatting, typing, packaging,
and frozen prediction/export comparisons are recorded separately under
`/tmp/tracrrna-quality-final`.

## Accepted limits and council discussion

The documented exclusions remain appropriate: this is atomic single-file
publication, without fsync or an export-wide transaction; a storage failure may
prevent diagnostic publication or temporary cleanup; the original command-log
publication/cancellation policy remains outside this bounded change. No new
promise about those behaviors is introduced.

Findings and explicit acceptance were sent to the scientific validation member
and lead. The scientific reviewer independently reported no cross-component
blocker, unchanged native command construction, and preserved policy thresholds,
rejection ordering, and ranking. The final council synthesis belongs to the lead
once the independent tooling review is also complete.

### Tooling blocker adjudication

After the runtime/publication review, the separate tooling member reproduced a
cross-component failure at `scripts/dependencies.py:96-98`: `update`'s final
staged-file `unlink` can mask either the primary lock-publication exception or
the combined publication/restoration error, including its required checkout
recovery instructions. I inspected both the code and the two-case probe/output
in `/tmp/tracrrna-tooling-cleanup-probe.py` and
`/tmp/tracrrna-tooling-cleanup-probe.txt` and agree this is a **blocking P2
correctness issue for overall council acceptance**.

The agreed bounded correction is to preserve the exception already propagating,
including its cause and recovery text, while attaching or logging secondary
cleanup details and the staged-file path. Add focused tests for publication plus
cleanup failure and publication plus restoration plus cleanup failure. Assert
that the disk lock and caller mapping remain unchanged, that the checkout state
matches whether restoration succeeded, and that required recovery guidance
remains visible. A separate executor will implement the correction; the tooling
reviewer must revalidate it.

Runtime/reporting acceptance remains unchanged. Overall council acceptance is
withheld until this tooling blocker is corrected and independently verified.

### Final council disposition

The tooling member independently re-ran the original failure probe after the
separate executor's correction and accepted the result. Publication failures now
preserve their original `OSError`; failed restoration preserves the complete
`RuntimeError` and its cause/recovery instructions. Cleanup failures add notes
with the residual staged-file path. The expected checkout, disk-lock, and caller
mapping states remain correct. Its 30 focused tests and scoped Ruff/format checks
passed, with branch and detached-head cases covered.

**Final concurrence: ACCEPT the complete change set.** The runtime, scientific,
and tooling members have exchanged their conclusions, and the sole blocking
finding is corrected and independently verified. No blocker remains. This final
disposition supersedes the temporary conditional acceptance above. Unchanged
runtime/reporting tests were not repeated; their earlier independent evidence and
the lead's native baseline comparisons remain applicable.
