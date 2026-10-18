# CRisprleader v2 — approved council plan and implementation status

Date: 27 September 2026. The owner confirmed the plan with “aprove”.
Standalone repository creation and implementation are authorized.

The private repository is [Alexander-Mitrofanov/CRISPRleader-v2](https://github.com/Alexander-Mitrofanov/CRISPRleader-v2).
Its checkout is `tools/CRISPRleader-v2`. Three independent workers delivered
source-bound genomic input/mapping, reference-container and generic numerical
infrastructure, and bounded runtime execution. The coordinator added the CLI,
independent output validation, exports, packaging and engineering checks.

**The scientific tool remains incomplete.** Automated safety restrictions stopped
biological research in two worker turns. No default biological reference bundle,
leader inference policy or independent biological evaluation was completed.
The working context-extraction command makes no leader prediction; `predict`
returns an explicit execution error. See the [current implementation status](../../tools/CRISPRleader-v2/docs/STATUS.md).

Review:

1. [Approved scientific and delivery plan](PLAN.md).
2. [Independent worker ownership and gates](WORKERS.md).
3. [Literature, reference audit and benchmark policy](EVIDENCE.md).
4. [Scientific-value and comparative-performance follow-up](SCIENTIFIC-VALUE.md).

The original completion criteria remain: working prediction, an actual usable
biological reference bundle with reviewed provenance, and independent biological
results meeting the predeclared gates. Passing synthetic tests or a motif scanner
does not satisfy those criteria. No biological accuracy claim is made.

No existing scientific tool or server runtime was modified. Identify, CasAndra,
CRISPRtracrRNA, evOr and SpacerPlacer retain their current responsibilities.
Integration into CRISPRloci and deployment are separate decisions; no tool lock,
active runtime or production service was changed.
