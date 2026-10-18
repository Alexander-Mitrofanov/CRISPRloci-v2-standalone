# Scientific acceptance gates — frozen before final predictions

Initial declaration: 27 September 2026. W3 independent evaluation worker.
No final predictions inspected. Record any amendment and its timing here.

Status: **not evaluated**. Biological work was halted by the coordinator after
a reported safety restriction. All requirements below remain requirements,
not achievements. The pre-output source-audit partition correction is recorded
in `partition.json`; no acceptance gate was relaxed to accommodate it.

## Integrity and honest claims (all required)

1. Every source has accession.version, original record URL, source identity,
   independently computed sequence SHA-256, retrieval date and rights note.
2. Every array geometry, observed repeat and predicted feature reconstructs
   exactly from that source; transformed fixtures reconstruct to the same
   oriented biological sequence. Required success: 100% of tested records.
3. Final genomes and close K-12 relatives are absent from fitting and
   calibration; all shared sequence identities are reported. Exact duplicated
   reference cores are reported as homology/regression cases, not evidence of
   generalization. The full set of fitting and calibration loci is auditable.
4. Experimental labels retain their true scope. Strand/junction experiments
   cannot become boundary truth. No full-leader/promoter/activity claim may
   be emitted from extracted windows, reference cores or motif matches.
5. Unsupported families, insufficient sequence and orientation conflict retain
   scientific no-call evidence. A broken reference or process execution is an
   error, never a biological negative.

## Initial prerelease functional gates

These are engineering/curated-panel acceptance gates, not a claim of general
biological sensitivity or a population accuracy target.

- Known-family final loci with independently reviewed orientation must not
  produce a confident call to the experimentally contradicted side. Report
  every abstention and its reason; do not count abstention as correct detection.
- At least one nontraining known-family natural locus must produce supported
  homologous-core evidence on the reviewed side. Otherwise the default bundle
  is not demonstrated usable beyond its own training examples.
- Every final novel-family locus must retain unknown/unsupported family
  coverage; a confidently misassigned specific family blocks the prerelease.
- All shuffled controls and uncertain natural alternative sides must be
  reported with their own denominators. A shuffled call requires adjudication
  and blocks a clean specificity claim; it cannot be removed post hoc.
- Reverse-complement, circular-origin rotation, source fragmentation and
  supplied-strand conflict transformations must preserve scientific evidence
  where the same bases remain available; losses of necessary sequence must
  become explicit partial/no-call results.

## Measurements and quantitative release claims

Report coverage, supported-side accuracy and false-call rate with numerator,
denominator and Wilson 95% intervals where a binomial interpretation is
justified. Also report abstention by family, natural-genome count and unique
core count. Small samples cannot establish broad accuracy.

Core endpoint error and interval overlap are evaluated only against independently
reviewed comparative-core labels of the same meaning. Adaptation-essential
regions and arbitrary upstream windows are excluded from that metric. If no
adequate boundary labels exist, endpoint accuracy is **not estimated**; the
release must say so, rather than reporting zero error against a chosen window.

Compare fixed-window extraction, motif-only evidence, sequence-family matching
and DNA profile inference. A route unavailable because of licensing, missing
models or missing executable is explicitly unavailable and cannot count as a
win. Thresholds must be fixed using development/calibration only. Quantitative
superiority or universal accuracy claims require a larger independent panel
and a separately preregistered study; this prerelease makes neither claim.

## Additive comparison-study design — 27 September 2026

The original prerelease gates above are unchanged. They do not establish
scientific superiority. The owner's request for stronger comparative performance
is addressed by `COMPARATORS.md`, `SUPERIORITY_PROTOCOL.md`, `comparators.json`
and `study-design.json`. These separate task meanings and require independently
reviewed labels, leakage audits, complete outcomes, practical margins and an
appropriate powered comparison before a scope-specific performance claim.

The candidate primary scientific endpoint concerns conserved-core localization.
The first implemented report evaluates supplied transcription-orientation calls
only, as a secondary descriptive aid; it is not a core predictor or a completed
biological evaluation. Practical margins, final sample size and confirmatory
statistics are not frozen. The design is therefore not completed preregistration.
