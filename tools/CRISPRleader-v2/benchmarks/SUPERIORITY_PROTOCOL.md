# Scientific value and comparison protocol

Version: 2026-09-27, designed before any v2 biological prediction results.
**Status: study design, not completed preregistration or demonstrated superiority.**
The existing small `partition.json` is a source-audit/prerelease candidate panel;
it is not a powered superiority benchmark. Its frozen history is retained.

## What would establish scientific value?

The useful research question is whether v2 can provide more correct, independently
supported annotations at an acceptable error rate, especially outside the most
represented reference families. Better packaging and coordinate handling improve
reliability but cannot establish improved biological prediction.

Priority hypotheses for future implementation/evaluation are:

1. Repeat-family-aware conserved-core evidence with explicit uncertainty improves
   localization or useful coverage over the original leader tool. This follows
   the task distinction in [CRISPRleader](https://doi.org/10.1093/bioinformatics/btw454).
2. A separately reported motif channel adds annotation evidence without turning
   every match into a leader call; applicability is family/context dependent,
   consistent with the [2025 motif study](https://doi.org/10.1089/crispr.2024.0093).
3. Preserving both orientations until evidence resolves them, and abstaining on
   conflict or unsupported families, reduces misleading annotations at a measured
   coverage cost. This is a testable proposal, not a literature-proven improvement.
4. Independent orientation evidence can expose disagreements worth reviewing.
   Acquisition orientation from [evOr](https://doi.org/10.1371/journal.pcbi.1013706)
   must not be silently relabeled as transcription direction or as ground truth.

These hypotheses are priorities, not algorithms or calibrated thresholds already
implemented. The biological reference, predictor and validation remain incomplete.

## Separate questions and labels

Conserved-core localization is the proposed primary scientific task because it
is the purpose of this standalone leader tool. Orientation is a secondary task;
the first implemented report covers only that secondary task and cannot establish
a win on leader boundaries. The primary overlap endpoint requires comparable
exact core intervals; uncertain endpoint labels belong in a separately declared
interval-aware analysis until their evaluation rule is frozen.

| Track | Eligible independent truth | Primary measurement | Prohibited substitution |
| --- | --- | --- | --- |
| Transcription orientation | Reviewed, single-direction transcription labels in a declared source frame | Equal-group mean of correct calls / all eligible loci, jointly with error and coverage | Agreement with another predictor, acquisition direction alone, or pseudo-labels |
| Conserved-core localization | Independently reviewed core boundaries with the same semantic definition | Group-averaged strand-aware interval overlap over all eligible loci; abstentions/errors have zero recovered overlap | Fixed extraction windows, functional-assay intervals, or labels produced by the method being evaluated |
| Distal core endpoint | Exact or interval-valued endpoint labels with recorded resolution | Error on called loci plus coverage; interval distance is zero only inside the admissible interval | Treating the known array junction as proof of upstream-endpoint accuracy |
| Motif/annotation evidence | Independent motif-class/annotation labels, with unverified cases separate | Task-specific precision/coverage where both positives and negatives are justified | Numerical PWM significance treated as experimental activity |
| Acquisition orientation | Independently reviewed acquisition labels or explicitly labeled simulations | A separate group-level analysis | Scoring these calls against transcription labels |

V2 leader-side evidence also is not automatically a transcription-direction
prediction. Every admitted tool output needs a documented semantic mapping to
the track's target, independent of whether the tool is v2 or a comparator.
Never pool simulated and natural biological cases in the same accuracy denominator
or uncertainty estimate; run and report them separately.

Bidirectional/conflicting transcription observations do not become arbitrary
binary truth. Preserve them for a separate analysis. Keep experimental and
reviewed-comparative evidence strata visible; computational labels are excluded
from accuracy denominators. More predicted loci are not automatically more true
positives. Unlabeled alternatives are not established biological negatives.

## Sampling, leakage and independence

Follow the data/optimization/model/evaluation reporting structure of
[DOME](https://doi.org/10.1038/s41592-021-01205-4). Define the intended population,
eligibility, source availability and evidence review before predictions. Include
bacterial/archaeal and common/rare-family strata where sufficient truth exists;
report unsupported or unsampled strata instead of extrapolating into them.

Split and audit before any feature discovery, model/reference construction,
calibration or threshold selection. Group related genomes and homologous
repeat/core loci so that close relatives and repeated copies cannot straddle
training and final evaluation. [GraphPart](https://doi.org/10.1093/nargab/lqad088)
demonstrates why biological sequence similarity must inform partitioning.
The exact similarity procedure, thresholds, taxonomy snapshot, group membership
and exclusions need to be frozen in the final study manifest. They are not yet
selected or executed here.

Audit comparator training/reference overlap too. If it cannot be reconstructed,
report the uncertainty and distinguish deployment performance from clean
unseen-family generalization. Do not claim that accession-only splitting removes
homology leakage. Keep all transformed/fragmented derivatives with their parent;
they test robustness and add no independent biological sample size.

Curators should be blind to model predictions; discrepancies need a recorded
second review/adjudication. Independent-group assignment is a scientific claim
requiring audit, not a fact established by supplying `group_id` in a JSON file.

## Operating point and complete accounting

At a threshold frozen on development/calibration data, publish all counts:
correct calls, wrong calls, abstentions, unsupported cases, execution failures,
missing results and unavailable tools. Use the same prespecified eligible locus
set for every tool; never compare each tool on only its own successful cases.
Report micro totals, equal-group summaries and label/family strata.

For N eligible loci, correct C and wrong W:

- Correct-call yield = C/N.
- Call coverage = (C+W)/N.
- Error among calls = W/(C+W), undefined with no calls.
- Wrong-call yield = W/N.

Risk and coverage must be considered together; a tool can appear highly accurate
by declining difficult cases. [Selective classification](https://papers.nips.cc/paper_files/paper/2017/hash/4a8423d5e91fda00bb7e46540e2b0cf1-Abstract.html)
provides the relevant conceptual framework. Use development data to select the
operating point. Final risk–coverage curves are descriptive, not a reason to
retune the reported threshold. An all-locus measure includes failures operationally;
also show completed-run scientific accuracy separately so failures are not
misdescribed as biological false calls.

## Confirmatory superiority rule — not yet activated

Before looking at final predictions, freeze:

- Named task/population and included comparator implementations.
- Primary endpoint and practical improvement margin delta.
- Maximum acceptable error rho and minimum acceptable coverage kappa.
- Whether experimental-only or a specified evidence mixture is primary.
- Independent grouping, resampling/statistical method, power target and sample size.
- Multiplicity correction, stopping rule and planned sensitivity analyses.

Choose delta/rho/kappa using intended scientific use and pilot variability,
without inspecting the final labels/predictions. Do not invent universal values
from the papers or choose values that make an observed result pass. The values
remain null in `study-design.json`, so the study cannot currently make a
confirmatory superiority claim.

Compare on the same units with paired group-level uncertainty. Related loci are
not independent Bernoulli trials; ordinary per-array Wilson intervals cannot
establish group-level generalization. With too few independent groups, report
inconclusive. A bootstrap does not compensate for poor labels or too few groups.

For the proposed orientation secondary analysis, define group-weighted coverage
as the mean call fraction across all eligible groups, including zero-call groups.
Define group-weighted conditional error as mean wrong-call fraction divided by
mean call fraction, with null when that denominator is zero. Do not average
conditional errors only over groups that happened to receive calls. The final
study must explicitly bind rho and kappa to its declared weighted quantities;
core-localization error semantics remain unresolved until labels and the claim
are fixed.

A claim against multiple competitors requires a predeclared family of comparisons
and familywise error control, for example an appropriate Holm-adjusted planned
test family or simultaneous intervals. [Demšar](https://www.jmlr.org/papers/v7/demsar06a.html)
motivates care in multiple-method comparisons; its cross-dataset procedures are
not automatically a test for dependent loci. The statistical unit and assumptions
must match this study. The implemented report's unadjusted descriptive intervals
are not the confirmatory test.

A scope-specific superiority claim requires the prespecified lower confidence
bound for improvement to exceed delta, with error and coverage meeting their
prespecified bounds. Name the comparator, endpoint, dataset/scope and uncertainty.
If any required eligible comparator is untested, no "better than all prespecified eligible
alternatives" claim is available. No finite inventory establishes superiority
to every tool that exists or will exist.

## Fair execution and attribution

Run the registry tracks in COMPARATORS.md. Hold hardware, threads, resource
limits and input information fixed within each track; randomize execution order
and prescribe repeats. Record initialization/model loading separately from
analysis time, cold/warm runs, peak memory, output validity and completion.
Unavailable tools are reported, not assigned zero accuracy or infinite time.

Ablate reference matching, motifs and optional context on the same splits. This
checks whether each component contributes scientifically; it does not replace
external comparisons. Preserve negative and inconclusive results. The novelty
claim must explain the demonstrated contribution, not the number of components.

## Present implementation boundary

`crisprleader benchmark` summarizes supplied orientation labels and tool calls;
it neither runs competitors nor produces labels/predictions. Its deterministic
metrics and paired equal-group bootstrap are engineering infrastructure. It
retains synthetic provenance and always reports that scientific superiority is
not established. Core/motif/acquisition evaluation, blinded biological curation,
competitor execution and the confirmatory statistical study are not completed.
The earlier restriction on biological reference construction has not been
bypassed by this methodological literature review.
