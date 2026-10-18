# Council decisions and review record

Convened 27 September 2026 in response to the project owner's request.
Three independent agent tasks researched bounded topics alongside a coordinating
reviewer. They then reviewed the concrete roadmap and returned corrections.
This is an internal agent council, not an external expert panel or experimental
replication. No numerical vote or confidence score is inferred from agreement.

| Member/task | Evidence reviewed | Recommendation |
| --- | --- | --- |
| `publication_parity` | Original article; visually inspected supplement; authors' interface; original tool papers | Restore original capabilities; retain exceptional arrays, prophage context and local reproducibility. |
| `v2_capabilities` | Repository guides, public contracts, adapters, publication path, frontend and pinned tool code | Expose existing capabilities first; preserve independent validation and release boundaries. |
| `scientific_evidence` | Primary method papers, experimental counterexamples, host/matching benchmarks and updated classification | Support the staged plan with explicit abstention and independent evaluation before stronger inference. |
| Coordinator | Cross-checks of source/code evidence, export boundary, check scripts and all reviewer findings | Adopt P0–P7 for a parity program; gate P8 separately; deliver planning documents only. |

## Decisions

1. **Keep the five existing integrations.** Add adapters and evidence views
   around them; introduce a new scientific dependency only for an actual gap.
2. **Sequence by reuse and prerequisites.** P0 establishes contracts; P1/P2
   improve existing evidence and inputs; P3/P4 restore locus/repeat analyses;
   P5/P6 restore target/host context. Reference release work begins early, and
   P7 can proceed alongside these once relevant contracts exist.
3. **Full parity includes conditional work.** Every F-item remains on the
   ledger. A gated or deferred feature is not considered implemented. Legacy
   accuracy claims and software versions need not be reproduced to preserve
   a useful capability, but any replacement must be validated and named.
4. **Separate biological concepts.** Preserve native scores, uncertainties,
   source coordinates, input types, reference identities and unsupported cases.
   Avoid turning a flank into a leader, a fold into activity, a match into an
   interaction, or a cassette count into completeness.
5. **Use updated literature without overstating coverage.** The 2026 matching
   benchmark informs candidate selection and the 2025 classification informs
   vocabulary management. Neither automatically improves an installed model.
6. **Defer unrelated expansion.** PAM inference remains conditional;
   anti-CRISPR discovery, novel effectors, guide design and broad defense mining
   need separate objectives and evidence.

## Corrections incorporated after council review

- Self-target defaults exclude all annotated array intervals, with array
  overlaps retained separately; excluding only the originating spacer is
  insufficient.
- Repeat reports must retain per-instance context effects and inspectable
  pairing probabilities. One consensus MFE or array-level score is insufficient.
- Exceptional-array fixtures explicitly include complete spacer deletions at
  different positions and partial repeats.
- Prophage validation measures biological region/boundary accuracy, not merely
  correctness of overlap arithmetic.
- Protein exports require translation or supplied-protein provenance. Genomic
  slicing alone cannot validate an amino-acid sequence.
- GFF3 missing-value notation follows column/attribute rules. Existing native
  Cas GFF is acknowledged; the gap is consistent, unified export.
- Retaining native Cas evidence does not permit publishing excluded JSONL or
  raw work directories.
- Missing-component hypotheses remain distinct from observed coverage; removing
  that inference is an explicit unresolved parity decision.
- PAM motif enrichment must not be described as an established interference
  requirement.

## Questions to resolve in future work

| Question | How it will be resolved | Blocks |
| --- | --- | --- |
| Which leader/reference implementation is maintained and reproducible? | Evaluate candidate tool/data licenses, installation and independent leaders. | P3 selection |
| How should multiple records map to organisms, assemblies and topology? | Design explicit input metadata and ambiguous-input behavior. | P2/P5 semantics |
| Which repeat classifier/family resource has adequate independent coverage? | Compare native assets, legacy families and subtype baselines with family holdouts. | P4 selection |
| Which host-spacer snapshot and search engine fit the server budget? | Curate sources; measure recall, ambiguity, storage, RAM and runtime. | P5/P6 selection |
| Which local prophage predictor covers intended host lineages? | Independent region/boundary benchmark and database/license audit. | Prophage parity |
| Can missing-component hypotheses be validated beyond observed coverage? | Mask known components; test incomplete/mixed/rare systems. | F10 closure |
| Which optional detector paths are actually present in the release? | Check pinned runtime dependencies/settings and irregular-array fixtures. | F14/F15 closure |
| What numerical acceptance targets and resource limits are realistic? | Freeze them using development data before a separate final evaluation. | Scientific release gates |

No questions block delivery of this plan. They are concrete research and design
gates for subsequent authorized work, not requests to begin implementation.

## Scope and verification of this planning task

Repository writes are confined to Markdown files in `improvenments/`. No code,
tool checkout, dependency lock, runtime, database, job or service was changed.
The supplement was inspected from a temporary research download; its stable
retrieval link is recorded in the literature register rather than vendoring
the paper into this project.

Document links, source anchors, cross-references and coverage were checked.
Application tests and runtime benchmarks were not run because no application
implementation was requested. The roadmap records the required future checks,
including explicit browser tests beyond `scripts/check.sh`. The workspace root
does not expose a Git repository, so no root Git diff/commit was produced.
