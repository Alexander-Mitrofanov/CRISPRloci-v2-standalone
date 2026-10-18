# CRISPRloci module integration council report

Review date: 14 October 2026. **This is the completed pre-implementation report.** All 15 requested councils finished before any implementation began. Later changes and their verification belong in the separate implementation record; they do not retroactively alter these findings.

The supported scientific workflows are substantially integrated, but several independent validation and information-preservation defects need repair. Two native environment overrides also undermine declared execution identity. No evidence supports indiscriminately replacing every reference with a larger dataset. The most important capability limitation is that CRISPRleader currently extracts flanks rather than performing the original leader-prediction workflow. The most important fresh scientific uncertainty is ISEScan's forward/reverse sensitivity under the tested compatibility environment.

## Scope, method and evidence

The repository pins 17 independent scientific sources. CasAndra and CRISPRidentify-v2 were excluded as requested; their existing interfaces were read only where a reviewed module consumes them. Each of the other 15 tools was reviewed **sequentially** by a chair and two specialist subagents. One specialist evaluated author intent and scientific execution, another evaluated data/output/handoffs, and the chair reviewed platform boundaries. Every council held an explicit cross-challenge and rebuttal round, resolved competing claims, and finished both specialists before the next council began. The detailed reports retain those discussions and evidence.

At the start, all 17 scientific checkouts were clean, matched `tools/tools.lock.json`, and had privacy hooks configured. The parent workspace's Git metadata was inaccessible, so original content and hashes were saved for review. Before this report was frozen, all **1,331 pre-existing source, test and documentation files** still matched their initial hashes; no implementation, tool update, reference replacement, runtime activation or production operation had occurred.

The initial `scripts/check.sh` gate passed **1,414 Python tests and 396 frontend tests**, plus lint, source-boundary, build and privacy checks. **70 optional native tests skipped** without explicitly selected executables. The local Chromium JSON handoff baseline passed upload, malformed-input recovery, incompatible/incomplete rejection, submission, compact-result and authenticated downloads, direct continuation, recovery and mobile layout. It used an intercepted fixture API, not production or a complete native pipeline. A sandbox FastAPI thread-bridge hang was resolved by running the same tests in the approved local execution environment, without weakening them.

Each report distinguishes inspected source, malformed-output probes, fresh native execution, historical evidence and biological evaluation. Native tests demonstrate execution and computational consistency; they do not establish new sensitivity, specificity or biological calibration. Test counts overlap across councils and must not be added into a purported independent total. A qualified pass means the supported method works within its declared scope while listed repairs or evidence limits remain.

## Four-axis verdicts

Each module name links to its substantive council discussion, source citations, primary literature and acceptance criteria.

| Module | Author-intended integration | Required or beneficial data | Compact scientific output | Downstream readiness |
| --- | --- | --- | --- | --- |
| [CRISPRtracrRNA-v3](module-council-2026-10-14/01-crisprtracrrna.md) | Qualified: supported model/complete paths are correct; independent gates need strengthening. | 27/27 bundled resources verified; broader subtype models need evaluation. | Repair lost IntaRNA energy and model-only interpretation labels. | Original DNA replay works; derived RNA-family transfer requires curated homology/boundary provenance. |
| [SpacerPlacer](module-council-2026-10-14/02-spacerplacer.md) | Qualified: grouping, biological order and native reconstruction are sound; correction/metric gates need repair. | No missing biological database; informative related cohorts and independent trees are user inputs. | Keep valid extremes/missing estimates; correct standalone orientation/tree wording. | Observed-array ports and detailed histories exist; inferred ancestors are not observed sequence inputs. |
| [evOr](module-council-2026-10-14/03-evor.md) | Qualified: both hypotheses, trees and ND threshold behave correctly; selected metric/status binding needs repair. | No missing biological database; same cohort/evaluation requirements as reconstruction. | Keep both hypotheses and unresolved selection; correct rho/LRT explanations. | Selected histories are inspectable evidence; retain input/reverse identities and uncertainty. |
| [CRISPRleader-v2](module-council-2026-10-14/04-crisprleader.md) | Pass for flank extraction; original leader prediction is not implemented. | Extraction needs no predictor database; motifs alone cannot create a validated predictor. | Explicit prediction_performed:false is correct; repair TSV strand tokens. | Source-bound flanks are reusable; neither side is automatically a validated leader. |
| [CRISPRrepeat](module-council-2026-10-14/05-crisprrepeat.md) | Pass for supported descriptive thermodynamics; strengthen impossible-pair gates. | Pinned ViennaRNA/Turner parameters present; no missing biological database. | Shared folds and sparse compact evidence are appropriate; full probabilities remain available. | Repeat/context evidence is reusable; folding does not establish homologous RNA families or processing. |
| [CRISPRspacer](module-council-2026-10-14/06-crisprspacer.md) | Qualified: exact and viral search semantics are sound; exact reported-hit/reference binding needs repair. | Two operator-pinned references required; production coverage unknown. Larger scoped collections merit evaluation. | Preserve exact strand TSV; simplify duplicate viral overview while retaining detailed evidence. | Source-derived spacer matches already feed Phigaro; ordinary ports replay input, not matched subjects. |
| [CRISPRmap-v2](module-council-2026-10-14/07-crisprmap.md) | Pass for historical offline mapping; distinguish independent fork from historical service. | Historical assets complete; modern catalogue integration is separate scientific work. | Repair effective-policy loss and misleading display truncation above 20 ties. | Repeat inputs and full tied mapping evidence are reusable; labels do not establish family homology. |
| [RepeatTyper](module-council-2026-10-14/08-repeattyper.md) | Qualified: frozen 37-label score fidelity passes; add guaranteed ACGT RC-consistency checks. | Model assets complete; training overlap and calibration provenance need independent evaluation. | Closed-set scores and unknown threshold remain explicit; correct stale TSV/provenance documentation. | Preserve repeat occurrences/origins; subtype alone is insufficient for RNA-family or evolutionary grouping. |
| [SpacePHARER](module-council-2026-10-14/09-spacepharer.md) | Pass for native combined evidence and grouping, subject to run-wide calibration interpretation. | Pinned matched target/control panels required; expansion must re-evaluate calibration. | Remove unconditional calibrated wording; clearly state run-wide pooled populations. | Groups and occurrence provenance retained; association is not infection or verified host assignment. |
| [CCTK](module-council-2026-10-14/10-cctk.md) | Pass for observed sharing/equality/components; no evolutionary inference claimed. | No reference database for compare; standalone MinCED/BLAST workflows are separate. | Restore coordinate and uncalibrated-score semantics; uploaded annotations remain client-supplied. | Observed arrays retain occurrences/deletion slots and replay; fix stale offline/live browser test. |
| [clinker-v2](module-council-2026-10-14/11-clinker.md) | Qualified pass for bounded all-CDS global protein alignment; upstream clustering/figures are outside scope. | No missing reference database; users must select related bounded loci. | Current full alignment evidence matches v1; repair leading-gap TSV escaping. | Original DNA and checked protein evidence available; ambiguity-normalized caller compatibility needs evaluation. |
| [ISEScan-v2](module-council-2026-10-14/12-isescan.md) | Faithful pipeline; native orientation sensitivity needs exact-release evaluation. Repair double-gap gate drift. | 266 HMMs+355 singleton queries complete; new ISfinder acquisition has permission/compatibility constraints. | Full TIR/ORF/partial-call evidence and typed TSV are appropriate; preview limits explicit. | Target replay and source-bound IS evidence work; proximity does not establish disruption. |
| [DefenseFinder-v2](module-council-2026-10-14/13-defensefinder.md) | Supported historical repertoire method; isolate native config and strengthen native rank/length binding. | Pinned 1.2.3 assets complete; later releases need rights, compatibility and scope evaluation. | Full evidence preserved; fix numeric strand normalization and component-preview truncation notice. | Original targets replay; system/profile/proximity evidence is not function or a derived protein port. |
| [RNAalifold-v2](module-council-2026-10-14/14-rnaalifold.md) | Supported pinned comparative folding; isolate MAFFT helper selection. | No missing biological database; curated homologues and independent structures matter. | Retain promised columns/probabilities; repair pair-type normalization and clarify metadata/version documentation. | Curated RNA/DNA import and replay work; consensus/candidate sets are not automatic homologous-family inputs. |
| [Phigaro-v2](module-council-2026-10-14/15-phigaro.md) | Substantially faithful; disclose joint Prodigal training and preserve native/source GC distinction. | Full 9,518-profile pVOG specification pinned; local full-model rerun unavailable. PHROGs needs recalibration/evaluation. | Retain genes/regions/spacer evidence; repair TSV and add effective search interpretation policy. | Enable existing validated array/repeat/spacer ports; current export gate permits only target replay. |

## Agreed implementation checklist

The following is the deduplicated engineering scope for the separate implementation agent. P1 prioritizes trustworthy execution/evidence binding, P2 scientific field fidelity and interpretation, and P3 documentation or presentation maintenance; these are implementation priorities, not claims of a public exploit or proven wrong biological inference. The detailed reports provide precise reproductions and acceptance criteria.

**A01 · P2 · tracrRNA.** Preserve native energy and evidence labels; independently verify oriented candidate sequence, optional fractions and finite energy. Default labels only for absent keys in validated model-only mode; preserve explicit nulls and complete-mode labels.

**A02 · P1 · SpacerPlacer, evOr.** Validate published numeric/event fields, effective correction policy and selected-hypothesis correction status independently. Keep legitimate missing estimators, huge finite alpha and negative log-likelihoods; do not change formulas or scientific defaults.

**A03 · P2 · SpacerPlacer, evOr.** Correct standalone orientation/tree/null-count presentation and rho/LRT explanations, including native help. Standalone means input orientation not compared; rho correction does not recompute optimizer likelihood.

**A04 · P2 · CRISPRrepeat.** Reject impossible positive pair identities and minimum-loop distances at both independent gates. Preserve valid GU and crossing ensemble marginals and model-specific local probabilities; no refolding.

**A05 · P1 · CRISPRspacer exact search.** Bind reference metadata and reported spacer identity/sequence/eligibility to the immutable pinned catalogue using bounded integrity checks and indexed lookups. No exhaustive search replay; acceptance of reported hits does not prove no omitted matches.

**A06 · P1 · RepeatTyper.** Enforce prediction consistency for verified unambiguous ACGT reverse-complement-equivalent repeats independently. Keep occurrences, literal sequence counts, input order and origins; do not extend to ambiguous/unclassified rows without a guaranteed invariant.

**A07 · P2 · CRISPRleader, CRISPRspacer, clinker, Phigaro.** Preserve typed strand enums in public TSV and validated leading-gap protein alignments in clinker public/native TSV. Arbitrary identifiers/titles retain formula protection; no blanket escaping removal; native clinker correction needs its own source identity.

**A08 · P2 · CRISPRmap.** Preserve effective distance/top-k policy and accurately show counts/display truncation, including more than 20 native ties. Distinguish native reporting limits from browser preview limits; do not invent historical labels.

**A09 · P2 · CCTK.** Restore declared coordinate_system and uncalibrated score_semantics in public/compact comparison JSON; explain uploaded annotations. No new confidence or verification claim; preserve explicit unknowns and zero values.

**A10 · P2 · SpacePHARER.** Remove unconditional calibrated wording and explain run-wide pooled target/control calibration while groups remain distinct. No batch-invariance, independent-group-FDR or biological accuracy claim; do not change scores.

**A11 · P3 · leader, map, RepeatTyper, SpacePHARER.** Correct stale present-tense pins/versions, inaccurate TSV provenance descriptions and historical/fork attribution. Retain clearly dated historical verification records; links to source lock remain authoritative.

**A12 · P2 · CCTK browser acceptance.** Update stale visibility assertion and test that unadvertised offline examples cannot enqueue arbitrary/live jobs; retain advertised desktop/mobile/download flows. Do not weaken feature availability enforcement or bypass browser checks.

**A13 · P2 · ISEScan.** Reject native-forbidden double-gap TIR alignment columns in both platform gates. Boundary hardening, not a native scientific error claim; retain native partial-call semantics.

**A14 · P3 · CRISPRspacer viral overview.** Use a concise overview with complete detailed evidence still available in the dedicated view and downloads. Retain uncertainty, heuristic completeness, truncation and missingness; avoid hiding information required for interpretation.

**A15 · P1 · DefenseFinder.** Reject consumed undeclared MacSyFinder configuration with a DefenseFinder-local native preflight. An inherited setting changed score 2 to 6 without declared policy change. Match pinned config resolution; preserve clean native results and other tool environments. Empty --cfg-file is insufficient. Propagate native source identity.

**A16 · P1 · DefenseFinder.** Independently bind native global protein FASTA rank and exact protein length in backend checks. Global rank 64/65 differs from per-source position 22/23; preserve both, shared assignments and ties. No scientific rerun or wholesale Integration duplication.

**A17 · P2 · DefenseFinder.** Normalize native numeric strand locally and show component-preview truncation. Preserve unknowns; do not change unrelated shared string-strand policies or complete downloads.

**A18 · P1 · RNAalifold; shared managed MAFFT launchers.** Remove inherited MAFFT_BINARIES from both RNA native subprocess calls and future generated managed launchers. Preserve caller environment and temporary-directory settings; transparent override marker must not run; retain pinned 7.526 outputs. Shared launcher also serves isolated SpacerPlacer/evOr; never edit active runtimes.

**A19 · P2 · RNAalifold.** Preserve bounded valid pair-type count objects in frontend normalization; clarify compact metadata exception and historical versus managed MAFFT versions. No extra display required, no lost promised probability/column evidence, no invented counts or missing values.

**A20 · P2 · Phigaro context.** Publish existing typed arrays/repeats/spacers from already validated detection evidence and declare their available output ports. Use existing builders and source/native origins; preserve accepted-category filtering, destination eligibility, incomplete/cap states and browser validation. No detector rerun or invented grouping.

**A21 · P2 · Phigaro context.** Publish effective source-spacer scope, query eligibility, exactness, all-array exclusion, completion and PAM/targeting limitations in full/compact evidence. Propagate independently through both validation contracts and fixtures. Test no hits, excluded queries, short sources and palindromes; do not promise excluded raw work rows in downloads.

**A22 · P2 · Phigaro context.** Explain joint single-mode Prodigal training, code 11 and appropriate genome/assembly input before submission and in methods. Multiple unrelated genomes can alter predictions for unchanged records. Keep native gc and source_gc distinct; 20 kb is eligibility, not validated accuracy. Do not change calls, thresholds or caller mode.

The exact-reference repair must independently bind the catalogue's internal and external metadata and reported rows to the pinned identity. Strictly byte/time-bounded streamed integrity hashing is distinct from an unbounded SQL scan or rerunning the entire search. Preserve cancellation, immutable read-only access, indexed lookups and row-size limits. Checking reported matches **does not prove that no matches were omitted**; do not relabel that residual trust boundary as exhaustive independent verification.

TSV repairs must be column-aware. Validated `+`/`-` enums and restricted protein/gap alphabets need exact scientific values, while arbitrary IDs, names and titles still need spreadsheet-formula protection. The ISEScan and DefenseFinder TSVs did not show the same strand defect; do not apply indiscriminate escaping changes.

No accepted action changes alignment/scoring formulas, probability models, biological thresholds, model bytes, training labels or reference scope. Correcting declared configuration ownership preserves the intended baseline. Native source changes require their own tested local commits, exact lock/source references and package identities; no push or activation is implied.

## Reference-data decisions

Required runtime assets, optional reference breadth, evaluation truth sets and user-supplied comparative inputs are different requirements. A larger collection is useful only with compatible biological scope, provenance, rights, identity, deduplication, calibrated interpretation and resource bounds. Jobs must continue to use immutable operator-provisioned references and must never acquire or update them themselves.

- **tracrRNA:** all 27 bundled resources matched their hashes. The large Atlas inventory and the relevant Rfam seed were already represented in development documentation. Useful next work is subtype-stratified, overlap-aware evaluation on native transcripts and negatives. Type V remains scoped to V-K; additional Cas12 subtypes require evaluated models.
- **SpacerPlacer and evOr:** there is no missing biological lookup database. Appropriate related/diverse array cohorts and justified independent trees are scientific inputs. Larger unrelated cohorts do not repair model assumptions or determine acquisition orientation.
- **Leader extraction:** no predictor database is missing from the implemented extraction operation. Restoring original leader prediction requires an inference method, appropriate curated references and held-out validation. New scoped leader motifs are not a universal drop-in predictor.
- **Repeat folding and RNA families:** pinned thermodynamic/native resources support the calculation. Experimental processing or structural truth sets would support new claims; they are not required lookup databases for the current descriptive method. RNA-family expansion requires curated homologues, comparable boundaries, verified orientation and useful diversity rather than duplicate abundance.
- **Exact spacer and viral search:** two different reference bundles are required. The review shell did not establish full production reference coverage. Documented historical corpora and tiny release fixtures are not proof of current deployed coverage. A quality/ecosystem-defined JGI SpacerDB subset or a scoped viral GenBank extension merits a pilot with accession/release identity, filtering/deduplication, added-match evidence, resource measurements and recalibrated interpretation where applicable. Current import limits rule out treating the largest public collections as drop-in replacements.
- **Historical repeat mapping:** the required historical catalogue is present. A broader modern catalogue is a separately labelled and evaluated capability; its presence in standalone source does not mean the web workflow already consumes it.
- **RepeatTyper:** frozen 37-label model assets and computational score fidelity are present. Training overlap, calibration and broader label coverage remain evaluation/retraining questions. Model scores must not become calibrated probabilities.
- **SpacePHARER:** target and control panels must remain matched and explicitly identified. Calibration pools submitted run-wide score populations; per-group output does not imply independently calibrated groups or batch invariance. More sequences require renewed calibration and biological evaluation.
- **CCTK comparison and clinker:** these modes need appropriate observed arrays or related loci, not a universal downloaded reference database. Standalone optional discovery workflows, clustering, synteny ordering and figures are separate capabilities.
- **ISEScan:** 266 HMMs and 355 singleton queries are complete for the pinned method. A new ISfinder acquisition has permission and model-compatibility constraints. More models cannot be assumed to repair the observed translated-peptide/HMM sensitivity.
- **DefenseFinder:** model 1.2.3 assets are complete. Later model metadata introduces unresolved rights declarations as well as changed coverage. Resolve rights, engine/rule compatibility and evaluation before a versioned expansion. The model does not assess CasFinder or all later systems.
- **Phigaro:** provisioning defines the complete 9,518-profile pVOG bundle, but the prior full local asset was absent, so a fresh full HMM search was unavailable. This is not evidence of a deployed-resource defect. PHROGs is a different scientific reference, requiring profile-semantic compatibility, recalibration and held-out host/taxon, GC and contig-length evaluation. Preserve native GC used by the inherited calibration separately from exact source GC. Coherent genome/assembly input matters because Prodigal trains jointly across eligible records.

The detailed reports link the primary dataset/method sources and distinguish freshly verified facts from retained historical receipts. No council established that automatic reference replacement would improve biological accuracy.

## Output and pipeline interpretation

A compact `result.json` is validated scientific evidence. A typed sequence packet is submission input. `input-sequences.json` ordinarily replays the original normalized input; it is not a derived finding merely because it accompanies a result. Existing explicit derived connections, including source-bound spacer evidence feeding Phigaro context, are evaluated separately.

A valid downstream connection needs preserved source identities, sequence/coordinate frame, strand and origin semantics, eligibility for the destination, and any curator decision about homology, grouping or orientation. Subtype labels, shared family names, coordinate proximity, candidate windows, a folded consensus or a similarity link do not supply those decisions. A new protein, annotated-locus, candidate-RNA or inferred-ancestor port requires a versioned contract and a real consumer, not a cast into an unrelated existing port.

Missing values remain missing; zero observations remain zero; unassessed coverage, unavailable models and unresolved hypotheses remain explicit. Native detail in protected full artifacts is not lost merely because the browser preview is smaller. Conversely, a full artifact does not excuse a mislabeled or silently truncated preview. Preserve concise qualifications needed to interpret results rather than treating all explanatory text as clutter.

The councils rejected unsupported reductions of clinker alignments and RNA-family columns/probabilities because those are promised v1 scientific contents. They also rejected a blanket finite/nonnegative rule across evolutionary CSV fields: deliberately unavailable estimators, negative likelihoods and finite extreme alpha values require field-specific semantics. No new confidence, homology, function, activity or organism-wide absence claim is justified.

## Scientific follow-ups that are not automatic code changes

1. Reproduce the ISEScan forward/reverse and fragmented-input controls in the exact release environment and compare the original author implementation. The measured sensitivity arose before coordinate serialization, at translated ORF/HMM evidence, under a compatibility environment differing in Python/BLAST versions. Define sensitivity, false-positive, boundary and resource criteria before altering caller settings, combining orientations or updating models.
2. Evaluate clinker compatibility with ambiguity-normalized caller proteins using real interface cases and an explicit original-versus-analyzed sequence policy. The synthetic `MK` versus `MX` discrepancy is a conditional compatibility finding; the retained real ambiguity fixture passed. Do not alter excluded CasAndra or silently treat X as a universal wildcard.
3. Design bounded independent exact-search completeness verification separately from reported-hit membership checks. Establish shared work/candidate/time budgets and measure costs before adding another exhaustive search.
4. Evaluate new leader prediction, subtype/model coverage, calibrated biological claims, alternative alignment/folding strategies and larger reference panels against independent, overlap-aware truth sets. Do not use predictions generated by the same method as independent truth.
5. Design optional derived ports, standalone history exploration, sanitized GenBank export and upstream synteny/cluster displays around explicit consumer requirements and scientific interpretation. None is a missing prerequisite of an existing demonstrated pipeline edge.

## Implementation and release acceptance

Only after this report and every detailed council report are complete should a new dedicated subagent execute the engineering checklist. It must record every action's disposition in `module-council-implementation-2026-10-14.md`, preserve these pre-change findings, and distinguish completed repairs from evaluated scientific follow-ups and concrete blockers.

Keep Integration and backend validators independent, scientific imports out of the API/backend, and execution in isolated workers. Preserve atomic queue transitions, cancellation/deadlines, bounded reads, symlink rejection, per-job capabilities, retention and process containment. SpacerPlacer and evOr must retain separate environments. CasAndra and CRISPRidentify sources and pins remain unchanged.

Propagate changed contracts through independent validation, public and compact projections, frontend models/views, schemas/examples/fixtures, documentation and coherent package versions/wheels/source identities. Run focused regressions, relevant native controls, `scripts/check.sh`, affected browser flows and the privacy scanner. Fresh malformed-output rejection must coexist with valid native baselines, null/zero behavior and uncertainty. Native package tests and fixture browser tests must remain labelled as such.

Local disk had approximately 7 GB free during review, and the existing immutable toolsets did not match the current 17-tool source lock. Assess complete candidate staging before allocating large environments. Do not mutate/relabel an older toolset, fabricate readiness, remove unrelated user data or activate a partial candidate. Changed wheels and source identities can be verified independently; any uncompleted full native release gate must remain an explicit deployment prerequisite. No production deployment or publication was requested.

## Detailed reports and frozen evidence

- [CRISPRtracrRNA-v3 council](module-council-2026-10-14/01-crisprtracrrna.md)
- [SpacerPlacer council](module-council-2026-10-14/02-spacerplacer.md)
- [evOr council](module-council-2026-10-14/03-evor.md)
- [CRISPRleader-v2 council](module-council-2026-10-14/04-crisprleader.md)
- [CRISPRrepeat council](module-council-2026-10-14/05-crisprrepeat.md)
- [CRISPRspacer council](module-council-2026-10-14/06-crisprspacer.md)
- [CRISPRmap-v2 council](module-council-2026-10-14/07-crisprmap.md)
- [RepeatTyper council](module-council-2026-10-14/08-repeattyper.md)
- [SpacePHARER council](module-council-2026-10-14/09-spacepharer.md)
- [CCTK council](module-council-2026-10-14/10-cctk.md)
- [clinker-v2 council](module-council-2026-10-14/11-clinker.md)
- [ISEScan-v2 council](module-council-2026-10-14/12-isescan.md)
- [DefenseFinder-v2 council](module-council-2026-10-14/13-defensefinder.md)
- [RNAalifold-v2 council](module-council-2026-10-14/14-rnaalifold.md)
- [Phigaro-v2 council](module-council-2026-10-14/15-phigaro.md)

The adjacent report manifest records report checksums, original tool pins, baseline test counts and the source-hash inventory identity. It marks the report-before-implementation boundary; later implementation evidence is separate.
