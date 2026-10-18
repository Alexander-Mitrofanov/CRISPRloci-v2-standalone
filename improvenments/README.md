# CRISPRloci v2 — internal improvement plan

Original plan: 27 September 2026. Latest research update: 12 October 2026.
Status: historical research and planning, with subsequent implementation tracked separately.

The October 12 implementation request is tracked in [context integrations](../docs/context-integrations.md), including license decisions, limits and verification. The research reports below preserve the evidence and recommendations at their original dates.

Latest follow-up: [literature, biological workflows and community benefit — 12 October 2026](10-literature-workflow-council-2026-10-12.md)
tests the five additions against primary application studies and defines the
biological question, expected benefit and finished result for each. It revises
the ease-first priority and repeat-first RNAalifold scope in the
[engineering council](09-tool-expansion-council-2026-10-12.md), which reviews the
current twelve-tool implementation. Both follow the
[10 October council](07-tool-expansion-council-2026-10-10.md), whose RepeatTyper,
SpacePHARER and CCTK recommendations are now integrated.

Restore the original server's useful scientific capabilities while keeping v2's
isolated execution, explicit uncertainty and reproducible releases. This folder
uses the spelling requested by the project owner. It sits outside `frontend/`
and therefore outside the existing frontend-only public export.

The original planning council used three parallel agent reviewers (publication
parity, current-code audit, scientific evidence) and a coordinating reviewer.
Individual follow-up reports record their own reviews and discussions.
Recommendations are proposals for future work, not independent experimental
validation or permission to implement or deploy.

CRISPRidentify v2, CasAndra, CRISPRtracrRNA v3, evOr and SpacerPlacer are already
integrated. Preserve them. Reuse their supported outputs and entry points before
considering another scientific dependency.

Read the plan in this order:

1. [Feature parity and current evidence](01-parity-audit.md): original scope,
   existing functionality, latent capabilities and actual gaps.
2. [Literature and interpretation](02-literature.md): primary sources, what
   they support, and limits on the claims v2 should make.
3. [Delivery roadmap](03-roadmap.md): proposed work packages, dependencies,
   acceptance criteria and release gates.
4. [Council decisions](04-council-decisions.md): recommendations, disagreements
   resolved during review and remaining research questions.

Detailed P4 follow-up: [repeat-analysis implementation decision](repeat-analysis/DECISION.md)
and [evidence/council record](repeat-analysis/EVIDENCE.md). The repeat council
selects per-instance context evidence first, then reference retrieval, with
strand/subtype/family inference subject to separate validation gates.

Recommended order:

- Establish typed input/evidence contracts and expose existing array/Cas
  detail, protein workflows, interoperable exports and better locus navigation.
- Restore leader analysis, repeat structure in array context, and repeat-only
  analysis with reference provenance.
- Add a shared spacer-target analysis layer for self-targeting and viral input,
  then complete contextual reporting with locally run prophage annotation.
- Evaluate PAM-context analysis only after matching, orientation and subtype
  validation. Defer automatic anti-CRISPR discovery and unrelated expansion.

GenBank input, IS-element evidence, exceptional array structures, Cas module
interpretation and the legacy CRISPRmap workflow are explicitly tracked; they
must not disappear between milestones. Full parity is a final audit outcome,
not a claim attached to the first delivery.

The historical baseline is Alkhnbashi et al., *CRISPRloci: comprehensive and
accurate annotation of CRISPR–Cas systems*, NAR 49(W1), W125–W130 (2021),
[doi:10.1093/nar/gkab456](https://doi.org/10.1093/nar/gkab456), checked against
its supplementary material and the authors' server help. Published tool
benchmarks motivate candidate approaches; they do not establish performance of
this v2 integration. No scientific benchmarks or application tests were run for
this documentation-only task.
