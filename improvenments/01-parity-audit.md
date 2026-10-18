# Feature parity and current-code audit

Read-only audit dated 27 September 2026. “Existing” means supported by inspected
source/contracts, not independently rerun in this task. “Latent” means an
upstream capability or retained evidence lacks a complete web workflow.
“Gap” means no supported public path was found in the scoped audit; it does not
claim exhaustive absence from every vendored file.

Historical source keys: **A** = [2021 paper](https://doi.org/10.1093/nar/gkab456);
**S** = [supplement, including visually inspected Tables S1/S2](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8265192/supplementaryFiles);
**H** = [official server help](https://rna.informatik.uni-freiburg.de/Help.jsp?toolName=CRISPRloci).
The current help corroborates the interface but is not by itself a dated 2021
snapshot. Literature IDs below resolve in [the evidence register](02-literature.md).

| ID | Historical capability | v2 evidence/status | Proposed package |
| --- | --- | --- | --- |
| F01 | Genomic array detection; categories and orientation (A) | Existing; narrow previews and category policies [C1–C3]. | Preserve; P1/P7 |
| F02 | Leaders and downstream context (A,H) | Latent 200-nt flank extraction; no structured conserved-leader prediction found [C4]. | P3 |
| F03 | Repeat/spacer sequences and features (H) | Native canonical records and complete-locus `array-evidence.json`; limited UI [C2–C3]. | P1 |
| F04 | Repeat MFE and array-context structure comparison (A) | Native detector MFE feature; no equivalent per-repeat structure workflow [C5]. | P4 |
| F05 | Repeat-only orientation, subtype and reference matches (A,H) | Gap: no repeat input mode [C1]. | P4 |
| F06 | CRISPRmap family/motif exploration (H) | Gap: no supported public handoff or local equivalent found. | P4 |
| F07 | Cas annotation, cassette boundaries and subtypes (A) | Existing; genomic Cas genes and cassette intervals/classifications [C6]. | Preserve; P1 |
| F08 | Protein-only and ordered cassette analysis (A,H) | Latent CasAndra entry points; web is DNA-only [C1,C7]. | P2 |
| F09 | GenBank input, supplied CDS, completeness setting (S,H) | Latent CasAndra GenBank/GFF support; no web parser or explicit completeness input [C1,C7]. | P0/P2 |
| F10 | Cas modules, incomplete/unseen cassettes, missing-protein hypotheses (A,S) | Gap: public counts/classifications do not establish roles/completeness or these inference capabilities [C6]. | P7; E14–E16 |
| F11 | Self-target matches outside arrays (A) | Gap: no matching workflow found. | P5 |
| F12 | Prophage/Cas/genomic target context (S,H) | Cas intervals exist; prophage/target overlap layer absent. | P5/P6 |
| F13 | Viral input, protospacer mapping and host evidence (A,H) | Gap: no viral mode or host-spacer reference contract [C1]. | P5/P6 |
| F14 | IS interruption and irregular/partial arrays (S,H) | Detector has repair/degenerated-array paths; complete web exposure and exceptional-case parity remain unverified [C10]. | P7 |
| F15 | Expert scientific controls (H) | Current bounded API options cover category, evolutionary matching, bias and RNA family; not legacy detector/Cas controls [C1]. | P0/P7 |
| F16 | Interactive genomic overview and detail views (A,H) | Existing array/Cas/RNA tracks and history explorer; incomplete general locus navigation [C8]. | P1 |
| F17 | GFF and sequence/report downloads (A,H) | Existing safe ZIP/JSON/TSV/RNA FASTA; native Cas GFF is supported. Unified feature export and protein FASTA remain gaps [C9]. | P1 |
| F18 | Standalone availability (A,S) | Existing local setup/Integration; verify equivalent CLI paths for every addition [C11]. | All packages |

The supplementary comparison explicitly marks multi-FASTA DNA unsupported in
the original server. Preserve v2's multi-record support as an existing
improvement, while adding explicit organism/assembly grouping where needed.
Original complete/draft terminology does not justify assuming multiple input
records belonged to one genome. The paper's future-work discussion included
tracrRNA; this is already part of v2 and is not a new proposal.

## Code evidence

**C1 — Modes and accepted inputs.**
[API models](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/models.py:11)
define `detection`, `orientation`, `reconstruction`, `cas`, `tracrrna` and `loci`.
The [submission contract](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/models.py:34)
accepts DNA with a small bounded options set. These six modes already route to
the five tools recorded in
[the tool lock](/home/alex/CodexProjects/CRISPRloci-v2/tools/tools.lock.json).

**C2 — Native array detail.**
[Detector canonical records](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/components/report_model.py:1049)
contain consensus, repeats, spacers, architecture and features.
[Complete-locus publication](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/publication/annotations.py:44)
retains array evidence. Reuse these validated records instead of recomputing
sequences in the browser.

**C3 — Preview/category limitations.**
[Detection projection](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/summaries/detection.py:59)
is intentionally narrow; [row selection](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/summaries/detection.py:138)
excludes categories outside analysis policy.
[The runner](/home/alex/CodexProjects/CRISPRloci-v2/Integration/crispr_integration/runner.py:1092)
can request accepted-only output. A full-category explorer therefore requires
an output-policy review. Alternative candidate representations are not
necessarily distinct loci; see
[detector output semantics](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/docs/OUTPUTS.md:18).

**C4 — Leader naming trap.**
[LeaderSeqSearch](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/components/components_non_array_computations.py:512)
extracts an orientation-directed window. The legacy text label
[“Leader region”](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/components/components_output_maker.py:412)
does not turn it into a biological leader-boundary predictor. It can support
accurately labeled flank context, pending P3.

**C5 — Energy is not a structure report.**
[MFE extraction](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/components/compact_feature_extraction.py:329)
retains finite energy values as features. That alone does not provide
dot-bracket structures, ensemble uncertainty or array-context comparisons.

**C6 — Cas already has real public results.**
[Gene projection](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/summaries/annotations.py:152)
and [cassette projection](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/summaries/annotations.py:167)
expose intervals and classification evidence.
[Native cassette construction](/home/alex/CodexProjects/CRISPRloci-v2/tools/CasAndra/src/casandra/cassette.py:60)
uses gene evidence and a gap policy. Neither gene counts nor confidence alone
prove completeness, activity or association with a nearby array.

**C7 — Input opportunities already upstream.**
[CasAndra documentation](/home/alex/CodexProjects/CRISPRloci-v2/tools/CasAndra/README.md:54)
describes GenBank, supplied annotation, protein annotation and ordered-cassette
paths. [Ordered cassette output](/home/alex/CodexProjects/CRISPRloci-v2/tools/CasAndra/src/casandra/cassette_input.py:81)
explicitly lacks genomic coordinates. The current
[web adapter](/home/alex/CodexProjects/CRISPRloci-v2/Integration/crispr_integration/annotations.py:163)
uses genome prediction with automatic gene mode.

**C8 — Extend the existing interface.**
[Annotation tracks](/home/alex/CodexProjects/CRISPRloci-v2/frontend/src/components/results/annotations/AnnotationResults.vue:64)
already show arrays, Cas genes and RNA candidates. Their source metadata
currently comes through RNA evidence. Generalize source-centric navigation
without replacing the existing evolutionary explorer or claiming a new viewer
where a usable one already exists.

**C9 — Preserve the publication boundary.**
[Annotation exports](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/publication/annotations.py:33)
and [bundle rules](/home/alex/CodexProjects/CRISPRloci-v2/backend/src/crispr_web/publication/contracts.py:28)
already support safe downloadable data.
[CasAndra writes GFF3](/home/alex/CodexProjects/CRISPRloci-v2/tools/CasAndra/src/casandra/reporting.py:84),
which is bundle-eligible; complete-locus internal `work/` paths are excluded.
Raw protein JSONL and `.faa` are not automatically public. Generate dedicated
validated CDS/protein exports rather than weakening the allowlist wholesale.

**C10 — Detector settings need a focused audit.**
[Current profiles](/home/alex/CodexProjects/CRISPRloci-v2/tools/CRISPRidentify-v2/README.md:49)
include orientation, degenerated-array search and repair paths; the checkout
also contains an optional IS-element environment. Installed support, effective
runtime settings and every irregular-array case were not executed in this
planning review. “Not evaluated” must remain visible until release tests exist.

**C11 — Existing reproducibility.**
[Project setup](/home/alex/CodexProjects/CRISPRloci-v2/README.md),
[architecture](/home/alex/CodexProjects/CRISPRloci-v2/docs/architecture.md) and
[scientific contracts](/home/alex/CodexProjects/CRISPRloci-v2/docs/scientific-contract.md)
already define standalone operation, snapshots, source binding and uncertainty.
The [public exporter](/home/alex/CodexProjects/CRISPRloci-v2/deploy/export-pages-frontend.sh)
copies only `frontend/`; this internal planning folder is outside that source.

## What counts as parity completion

Every F-item needs a demonstrated user workflow, validated evidence contract,
export, negative/uncertain fixture and documented local execution path where
applicable. Mark the outcome as retained, restored, validated alternative, or
explicitly changed/deferred. A changed/deferred original feature remains an
open parity item; it cannot disappear under a broad “complete annotation” label.
In particular, observed Cas component coverage alone does not replace an
evaluated missing-component hypothesis feature.
