# Tool integration audit — 3 October 2026

Nine scientific tools are pinned in `tools/tools.lock.json`. Eight already had
queued web workflows. This change adds the missing CRISPRmap reference-mapping
workflow, so all nine pinned tools now have API, isolated-worker, validation,
result-view and download integration in source. This is a local implementation,
not a production deployment or verification of the currently deployed release.

## Evidence reviewed

Project chats reviewed include **Implement CRISPRmap v2**, **Assess CRISPR loci
portal methods**, **Plan protospacer search tool**, **Prepare tools repository
folder**, **Plan CRISPRloci improvements**, **Test web server for bugs**, **Add
database search for arrays**, **Simplify interface and results tabs**, and the
repeat-analysis and CRISPRcasDb catalogue implementation chats. Their proposals
were reconciled against the current lock, API modes, runner commands, independent
validators, release builder and frontend, rather than treating chat completion
messages as deployed-state evidence.

## Pinned tools

| Tool | Web integration after this change | Scientific boundary |
| --- | --- | --- |
| CRISPRidentify v2 | Detection, locus annotation, comparative evolution, repeat/leader context | Array categories and strand evidence retain uncertainty. |
| CasAndra | Cas analysis and complete locus annotation | Predicted proteins/cassettes, not experimentally established activity. |
| CRISPRtracrRNA v3 | Type II/V-K screens and complete locus workflow | Candidates and retained evidence, not confirmed transcripts. |
| SpacerPlacer | Reconstruction and evOr-selected histories | Requires eligible comparable arrays; separate environment. |
| evOr | Orientation-aware evolution | Acquisition-order evidence differs from transcript strand; separate environment. |
| CRISPRleader v2 | Leader context and locus workflow | Both adjacent windows only; no validated leader prediction. |
| CRISPRrepeat / ViennaRNA | Isolated repeat and detected-array context folding | Descriptive structural evidence; no inferred cleavage positions. |
| CRISPRspacer / BLAST+ | Reference-spacer and viral-reference searches | Reference-dependent modes require operator configuration. Self/supplied-target CLI workflows are still not web modes. |
| CRISPRmap v2 | **New `repeat_map` mode:** DNA repeats to a frozen annotated JSON reference | Historical exact annotation recovery and nearby sequence candidates; no subtype, transcription-strand or structural-motif prediction. |

CRISPRcasDb is already exposed through the separate catalogue API and Database
view, including per-repeat/spacer searches. CRISPRstrand, MAFFT, IntaRNA, Infernal
RNAfold/RNAplfold (through ViennaRNA bindings), and the other native tools used by the pinned workflows are dependencies, not
additional standalone web submission modes.

## Other named methods and remaining gaps

The planning chats mention **RepeatTyper/CRISPRCasTyper**, **Spacer2PAM**,
**CRISPRclassify**, **CRISPRcasIdentifier**, **CCTK** (including CRISPRdiff/Constrain), **clinker**, **geNomad**, **VIBRANT**, **PHASTER**,
**iPHoP**, **CRISPRscope**, **CRISPRTarget**, and historical **CRISPRmap/LocARNA**
methods. These are proposals, alternatives, literature precedents or historical
baselines; they are not all installed/pinned executable integrations. None is
advertised as running in this change. Their scientific and licensing evaluation
remains as recorded in the [portal assessment](../improvenments/05-portal-council-2026-09-29.md)
and [connected pipelines plan](../improvenments/06-connected-scientific-pipelines.md).

The newer standalone CRISPRmap **2.1.0a1 modern catalogue/index/MAP2 group search**
is also distinct from the **2.0.0a1 reference mapper currently pinned by this
server**. Its large catalogue lives outside the local runtime and uses another
output/storage contract. This change does not upgrade the scientific pin, copy
that database, or claim modern catalogue coverage for the historical web mode.

The next concrete gaps are CRISPRspacer self/supplied-target web workflows,
modern CRISPRmap catalogue integration, detected-array-to-analysis handoffs,
feature-overlap evidence, comparison matrices and validated leader prediction.

## New workflow and release contract

See [repeat reference mapping](repeat-mapping.md) for operator setup and limits.
The API stores the operator-pinned reference digest in the same atomic queue
transaction as the job. A changed reference fails the queued job. Commands stay
shell-free and scientific execution stays within the existing isolated worker,
with cancellation/deadline/storage checkpoints and capability-protected exports.

Integration checks completion, source identities, policy, hit bounds and
annotation assignments. The backend independently reads the pinned reference,
recomputes nearest global-edit ties on both orientations with a different
bounded dynamic-programming algorithm, validates all label conflicts (including
undisplayed ties), and checks alignments and counts. It imports neither the
scientific mapper nor Integration validators.

The UI distinguishes exact, near, no-hit and unsupported evidence. Missing labels
remain null. Full query-level JSON/TSV downloads are generated from validated
native output, never from a capped summary. Native inputs, report paths and
commands remain private. The top-20 native alignment limit is disclosed; full
annotation aggregation still includes every tied reference occurrence.

Application versions are backend/frontend **2.4.0**, Integration **0.8.0**.
Scientific manifests use **6.0.0** when CRISPRmap is present, preserving older
manifest support. The lock's new runtime profile creates a new toolset identity
without changing scientific commits or mutating an existing runtime. The release
gate now includes real repeat-mapping controls and public-archive checks.

## Verification

Final `scripts/check.sh` result: **601 Python tests passed, 1 opt-in installed-
tool lifecycle test skipped; 197 frontend tests passed**. Lint, source boundaries,
production frontend build and private/public data scans passed. Both candidate
wheels built successfully and imported with the new mode in isolated Python.

- Real CRISPRmap subprocess tests cover exact/reverse/near/no-hit/ambiguous input,
  omitted ties, conflicting/missing labels, changed reference/input, invalid
  alignments, booleans, duplicate JSON keys, symlinks and stage/count tampering.
- HTTP-to-real-worker test covers admission limits, durable reference identity,
  completion, authenticated JSON/TSV/ZIP/provenance, access isolation and failure
  after a queued reference changes.
- The exact pinned 2.0.0a1 source was exported without changing the checkout and
  tested against the real 3,527-entry historical reference: 2 exact, 1 near,
  1 no-hit and 1 unsupported query passed the independent validator.
- Browser acceptance covers capability advertisement, limits, submission,
  result filtering, reload recovery, authenticated download, missing evidence,
  expanded alignments at 390px width, and no JavaScript/CSP errors.
- Local evidence is under ignored `output/integration-audit-20261003/`; built
  candidate wheels are under ignored `.runtime/candidates/repeat-map-20261003/`.

No production service was operated. Building/activating a full nine-tool release
and configuring its chosen biological reference remain deployment steps. The
source archive has no accessible parent Git metadata, so parent publication
hooks cannot be installed here; install `.githooks` in the real checkout before
publication. The privacy scan also found and removed a pre-existing host alias
from `docs/catalogue.md`, replacing it with the private operations record pointer.
