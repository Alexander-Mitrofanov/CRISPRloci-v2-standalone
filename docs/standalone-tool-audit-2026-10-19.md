# Standalone tool bug audit — 19 October 2026

Targeted engineering review of all 17 bundled tool entrypoints and their bounded
scientific workflows. Confirmed defects received one separate fix commit per
tool, with regression tests and matching bundled-source/runtime pins. Four tools
had no actionable finding in this pass. This is not a claim that these projects
are free of bugs or a new measurement of biological accuracy.

The tool suites passed **3,547 tests**, plus **2,284 CRISPRidentify subtests**.
Seven tracrRNA native integration tests were deselected; four optional
SpacePHARER comparisons requiring external original/calibration oracle binaries
were skipped. SpacePHARER's installed native engine and frozen native-output
controls ran. CRISPRmap's optional ViennaRNA controls also ran.

| Tool | Fix or reviewed boundary | Passed | Standalone fix commit |
| --- | --- | ---: | --- |
| CRISPRidentify-v2 | Terminate descendants after normal exits and SIGTERM; reject non-ASCII native feature inputs before FFI dispatch. | 481 | `b65db9f72eb5` |
| CasAndra | Replace six racy directory renames with atomic publication that preserves a concurrent output. | 73 | `cc8d0d3c8406` |
| CRISPRtracrRNA-v3 | Stop native descendants on success, nonzero exit and SIGTERM before validating outputs. | 442 | `4d6b81c77f6c` |
| SpacerPlacer | Reap MAFFT process groups on all exits and cancellation; retain a bounded diagnostic tail. | 167 | `5bb1ab697c78` |
| evOr | Reap MAFFT process groups on all exits and cancellation; retain a bounded diagnostic tail. | 155 | `161181582365` |
| CRISPRleader-v2 | Reviewed array geometry, circular flank extraction, projection and source binding; no actionable finding. | 320 | Unchanged |
| CRISPRrepeat | Reviewed input bounds, independent folding parameters and source-bound exports; no actionable finding. | 86 | Unchanged |
| CRISPRspacer | Prevent a concurrent empty directory from being replaced at the atomic publication step. | 254 | `48b3bb0ac5e0` |
| CRISPRmap-v2 | Own Turner2004 parameter tables and explicit folding settings; ambient defaults no longer alter results. | 238 | `75a48684d753` |
| RepeatTyper | Reviewed input bounds, immutable model/feature/label binding and exclusive publication; no actionable finding. | 68 | Unchanged |
| SpacePHARER | Preserve concurrent outputs and canonical LF source bytes so fresh clones and Git archives match their receipts. | 101 | `eed481c2e941` |
| CCTK | Reviewed observed/deletion identity, order, sharing metrics and bounded/exclusive I/O; no actionable finding. | 160 | Unchanged |
| clinker-v2 | Refuse replacement of a result directory created while comparison exports are being written. | 717 | `a3261491d758` |
| ISEScan-v2 | Contain and reap the native group; timeout raises a structured failure instead of signaling the wrapper itself. | 86 | `bf330606a600` |
| DefenseFinder-v2 | Reject impossible coverage fractions, negative component E-values and negative system completeness. | 55 | `bf9f9afb6f1b` |
| RNAalifold-v2 | Contain both MAFFT probes and alignments on success, failure, timeout and SIGTERM. | 60 | `4f0bba3a7d98` |
| Phigaro-v2 | Bind the input digest to the parsed byte snapshot and contain native descendants on every exit. | 84 | `8fd06872275a` |

The process regressions cover successful exits, nonzero exits, timeouts, SIGTERM
and calls from worker threads. Publication controls preserve the inode and
contents of competing directories, files and symlinks. Folding controls compare
exact outputs before and after changing process-global native defaults and
thermodynamic tables. DefenseFinder controls reject malformed native evidence;
Phigaro controls mutate the input after parsing and verify the original digest,
source identity and normalized sequence remain bound together.

`bash scripts/check.sh backend` passed: **1,523 passed, 70 skipped**, followed by
a passing platform Ruff check. The skips require unavailable optional
executables or managed runtime resources in this standalone source checkout.
The check ran in the
local execution environment required by the documented FastAPI test thread
bridge restriction. Changed Python files passed focused Ruff checks; tracrRNA's
changed runner passed strict mypy. `scripts/sync_tools.py --check` verified all
17 bundled trees. Fresh Git-archive and local-clone checks also verified all
17 source identities, including SpacePHARER's vendored Windows build files.
An additional 84 focused provenance, runtime and CLI-contract repository tests
passed after the final pin updates.
The private-data scanner reported zero findings.

The original upstream commits and `upstream_tree` values remain provenance.
Each patched tool's `bundled_tree` in `tools.lock.json` matches its exported
`tree` in `bundled-sources.json`. These pins change the runtime identity and
reject mismatched or missing bundles, preventing reuse of a pre-fix toolset.
Package versions, scientific thresholds, model assets and report schemas remain
unchanged. Existing installations must be rebuilt and release-gated explicitly
to use these patches.

Verification reused isolated local environments and added disposable test
dependencies outside the checkout. ISEScan selection controls used Python 3.12
with NumPy 1.26.4, SciPy 1.13.1, Biopython 1.85 and fastcluster 1.2.6; its deployment
recipe remains Python 3.11. DefenseFinder checks covered its shipped offline
context adapter. Retained legacy renderers, training workflows and unrelated
upstream commands were outside this focused audit. A fresh complete installation,
full native biological release gate and browser acceptance were not run. No UI
flows changed, no production service was operated, and no commits were pushed.
