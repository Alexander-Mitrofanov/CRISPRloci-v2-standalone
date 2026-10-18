# Proposed delivery roadmap

All items below are future work. No application changes, scientific-tool updates,
database builds, service operations or deployments are part of this plan.
Package IDs map to the parity audit. Priorities express dependency order, not
calendar estimates; runtime and database costs still need measurement.

## P0 — Define evidence, inputs and benchmark ownership

Agree on a versioned contract before adding modes. Distinguish genomic DNA,
protein sets, ordered proteins with explicit grouping, consensus repeats and
viral nucleotide sequences. File order alone must not establish genome order.
Preserve original IDs, input digests, molecule type, source record, coordinate
convention, strand and transformations such as reverse complementation or U/T
normalization. Add genome/assembly grouping: separate FASTA records are not
automatically separate organisms, nor automatically contigs of one organism.

Model each analysis stage with explicit `not_requested`, `completed`, `no_hit`,
`unsupported`, `failed` and `truncated` semantics, mapped deliberately to the
existing job lifecycle. These are proposed evidence states, not a request to
replace the durable queue state machine. Each result needs tool/model/database
identity, options, source binding and a reason for unavailable values.

Keep array strand, predicted transcription direction, leader interval and evOr
acquisition-order evidence separate. Record raw method scores in their native
units. Store native Cas labels with their classification vocabulary version;
any newer taxonomy mapping is a separately versioned interpretation.

**Acceptance:** reviewed schemas and fixtures cover missing, conflicting,
unsupported and no-hit evidence; input alphabets cannot be confused; negative
strand and circular intervals have a documented reversible representation.
Create a feature-to-fixture ledger for every parity item before implementation.
P0 must not silently change current modes or current evolutionary grouping.

## P1 — Make existing evidence useful and portable

Extend the existing result views and publication path, reusing existing tools.
Expose validated repeat/spacer sequences, consensus, category, raw detector
features and optional native evidence. Let users inspect excluded/low-scoring
array categories without admitting them to evolutionary analysis. This needs
an explicit detector output-policy audit as well as UI filtering; alternative
array representations must not be counted as separate biological loci. Display
native Cas protein sequences and explicit gene-to-cassette membership where
source output supports them.

Provide a common locus browser with contig selection, zoom, feature selection
and linked detail tables. Linear views are the default for fragments; use a
circular view only for records whose topology is known. Preserve the existing
evolutionary history explorer and its comparison semantics. Add stable
downloadable JSON/TSV/FASTA and a unified GFF3 export for supported genomic
features. First audit native GFF availability; a file extension allowlist is
not proof that a job generates the file.

Use one documented ID and coordinate mapping across array, repeat, spacer,
Cas, RNA and future target features. Genomic features and model predictions
must remain distinguishable. Use `.` only in GFF3 columns that permit it; omit
unavailable optional attributes or follow a documented attribute schema.
Protein-only results must not acquire
fabricated nucleotide coordinates. Include citations and evidence provenance
in every applicable bundle, beyond browser preview caps.

**Acceptance:** genomic sequences reconstruct from validated input intervals;
proteins bind to validated CDS/translation provenance or the supplied protein
FASTA; GFF3 round trips through an independent parser; detail IDs join to full exports;
preview truncation is visible; unknown and zero remain distinct; keyboard,
mobile, download, recovery and no-hit browser flows pass. Cross-linking is
proximity evidence unless a supported association method says more.

## P2 — Restore protein and annotated-genome input

Expose the pinned CasAndra protein annotation and ordered-cassette entry points
through new adapters. Keep unordered protein identification distinct from
cassette classification. Require explicit protein order/group provenance before
claiming cassette boundaries. Check which upstream inference assumptions need
to be surfaced, and do not call arbitrary FASTA ordering genomic context.

Add a bounded GenBank parser and normalization contract. Preserve accession,
topology, feature IDs, joins/complements and translated CDS provenance. Define
behavior for partial CDS, alternative translation tables, pseudogenes and
translation mismatches. Separate supplied annotation from re-predicted genes;
report conflicts instead of silently choosing a translation. Inspect current
tool support before choosing whether adapters or an upstream release own each
new path.

**Acceptance:** unordered proteins never generate positional cassette claims;
ordered positive fixtures and shuffled controls behave as specified; GenBank
and equivalent FASTA agree on source sequence identity; compound coordinates,
partial records and invalid translations produce defined results. Preserve
multi-record DNA support already present in v2. Depends on P0; P1 supplies the
presentation and export surface.

## P3 — Leaders and orientation evidence

Audit detector flank fields first. A fixed upstream window may be shown as a
flank but is not a conserved leader prediction. Evaluate CRISPRleader or a
maintained equivalent against curated leaders before selecting an implementation
and reference set. Return predicted conserved-core boundaries, supporting
family/model, native score, clipping at contig edges and unsupported cases.

Show detector transcription evidence, leader evidence and evOr's existing
group-level acquisition-order comparison side by side. Report disagreement.
Do not combine their scores into a probability or silently reverse the genome,
rewrite source coordinates, override evOr ND, or infer a leader from evOr alone.

**Acceptance:** independent known-strand/leader fixtures, reverse complements,
divergent leaders, edge-truncated sequences and conflicting evidence; report
leader detection coverage and boundary error separately. Use clade/family
holdouts where data permits. No hit is 'leader not detected', not biological
absence. Depends on P0/P1. Scientific basis: E3–E4 in the literature register.

## P4 — Repeat structure and repeat-only analysis

The [repeat-analysis council decision](repeat-analysis/DECISION.md) specifies
the implementation order, metric semantics, current reusable assets and
scientific release gates for this package. P4a and the repeat-only folding entry point are implemented in the separate CRISPRrepeat tool and web integration; see [implementation status](../docs/repeat-analysis.md). Calibrated strand/subtype/family inference and a curated public reference bundle remain pending.

Treat three outputs separately: isolated-repeat folding; repeat folding in full
array context; sequence/structure-family and subtype assignments. Evaluate
pinned ViennaRNA tools for MFE/dot-bracket output and local pairing probabilities,
retaining temperature, energy parameters, window/span, input orientation and
sequence. The historical isolated-versus-context comparison is a parity target;
a detector MFE feature or a tracrRNA pairing output does not complete it.

Include per-repeat-instance context effects and inspectable pairing-probability
plots/matrices, not only one array-level aggregate. Specify the
context-comparison metric and its handling of unpaired repeats
before coding. Use a descriptive label such as 'predicted pairing consistency';
do not label a model comparison experimental structure accuracy. For unresolved
orientation, retain clearly labeled alternatives when meaningful. A weak fold
must not veto an array.

Add a repeat-only workflow with an explicit DNA/RNA normalization policy,
orientation and subtype evidence, reference-sequence retrieval, alignment
details, genomic origins and database version. Evaluate reuse of existing
repeat/strand assets before selecting another classifier. Benchmark CRISPRmap
family assignment and repeat-subtype methods as distinct capabilities. Preserve
a user-controlled CRISPRmap handoff or provide a scientifically validated local
equivalent; do not silently send sequences to another server.

**Acceptance:** repeated sequences yield deterministic outputs with pinned
settings; structures are syntactically valid and sequence-bound; strand and
indel mappings survive export; non-hairpin and no-reference-hit cases stay
valid; family-disjoint tests measure classifier abstention as well as errors.
Depends on P0/P1, with P3 supplying optional additional strand evidence.
Scientific basis: E5–E8.

## P5 — One target-matching layer, two user workflows

Design one versioned spacer–target contract reused by self-targeting and viral
input. Separate matching from biological interpretation. Record query spacer
identity and origin, target source, both coordinate frames, alignment strand,
aligned length, query coverage, mismatches, indels, scoring policy, reported
hit count and any reporting cap. Define Hamming distance and edit distance
explicitly; do not reuse the current evolutionary `spacer_distance` option.

Benchmark suitable short-sequence search tools against an exhaustive oracle on
small fixtures. Evaluate the 2026 matching study's candidates and configurations,
including Bowtie1, indelfree and an appropriately configured alignment baseline,
before selecting an engine. Pin all-hit/reporting behavior. Handle ambiguity,
low complexity, reverse complements, duplicates, multi-mapping and circular
targets deliberately. Include complete-query and partial-alignment policies as
separate concepts rather than one permissive similarity threshold.

For **self-targeting**, search the defined genome/assembly, exclude overlaps
with all annotated CRISPR-array intervals from the default self-target set,
retain those overlaps in a separate inspectable category, and retain
Cas/prophage/other-region overlaps as potentially simultaneous annotations.
Multiple uploaded genomes need an explicit within-genome versus cross-genome
policy. Keep uncertain array categories inspectable but separate from the
default target-evidence cohort. Output 'candidate genomic match', not confirmed
autoimmunity or anti-CRISPR activity.

For **viral input**, search a frozen host-spacer reference release. Report
candidate host associations with unique spacer counts, source arrays,
accessions, alignment evidence, taxonomy level and shared/ambiguous matches.
Deduplicate reference copies without discarding occurrence provenance. No
matching host means no supported association in that release. RNA input must
retain its molecule identity; U/T normalization does not establish a DNA-targeting
mechanism. Viral sequences carrying CRISPR arrays require explicit masking or
overlap annotation to avoid treating spacer-to-spacer identity as ordinary
protospacer evidence.

**Acceptance:** exact expected hits on seeded synthetic fixtures, controls for
array self-matches and duplicate inflation, sensitivity/precision on held-out
known associations, abstention at unsupported taxonomic ranks, and resource
limits on repetitive input. Database-related false positives are evaluated
with length/composition-matched controls. Thresholds must be chosen on a
development partition and frozen before final evaluation; no universal
biological accuracy percentage is claimed. Depends on P0/P1 and P6 reference
releases for database-backed viral analysis. Scientific basis: E9–E12.

## P6 — Reference releases and prophage context

Create operator-managed immutable releases for leader, repeat and host-spacer
references and any prophage models/databases. Record accession versions,
source URLs, retrieval/build dates, licenses, taxonomic snapshot, detector
release, array quality policy, deduplication, checksums and index parameters.
Set measured storage/RAM limits and update cadence before deployment. Tool and
database identities must both travel with results. Preserve old releases for
reproducibility and rollback.

Evaluate a locally runnable prophage method; PHASTEST is a literature-supported
candidate, not a selected dependency. Compare its scope, database costs,
redistribution terms and benchmark performance with suitable alternatives.
Preserve boundary uncertainty and prediction quality when overlaying targets.
Distinguish no overlap from annotation not run, unsupported host lineage and
annotation failure. Bacterial prophage benchmarks cannot establish archaeal
virus coverage. Do not make completed locus annotation depend on a live public
prophage service or automatically upload submitted genomes externally.

**Acceptance:** rebuild/verify a frozen miniature database offline; a running
job keeps its original release; provenance resolves every reported hit;
prophage overlaps pass independent interval checks; measure prophage region
detection and boundary accuracy on independent curated genomes, stratified by
supported lineage and assembly fragmentation; stage failure does not
masquerade as absence. Research reference ownership in P0 and implement release
support before P3/P4/P5 require it. Scientific basis: E13.

## P7 — Close exceptional-array and Cas interpretation parity

Audit IS-element detection, degenerate/terminal and partial repeats, repeat
indels, complete spacer deletion at the beginning/middle/end, and interrupted
arrays against current pinned detector options and outputs. Expose
supported results with explicit 'disabled/not evaluated' states; unsupported
paths require upstream work and a new scientific release. Do not reimplement
the detector in the web layer. Expert controls should be bounded scientific
parameters or validated presets, with effective settings retained in output.

Add a versioned, curated mapping from supported Cas calls to module roles where
justified. Keep unassigned/ambiguous roles. Observed component coverage and
contig truncation can support a completeness description only relative to a
stated subtype model. A model-imputed missing component is a hypothesis, never
a detected gene with invented sequence/coordinates. Preserve native CasAndra
labels and validated scientific fields through approved artifacts alongside
any newer classification mapping; this does not authorize publishing native
JSONL or internal work trees. Evaluate a separately labeled missing-component
hypothesis report against masked-protein benchmarks. If only observed coverage
is delivered, record the historical inference feature as changed/deferred,
not equivalent or completed.

Explicitly evaluate incomplete and previously unseen cassette configurations
and retain interpretable component-association evidence where a validated model
provides it. Unknown-family candidates may be flagged with supporting evidence;
they must not be announced as novel effectors without independent follow-up.

**Acceptance:** irregular-array fixtures cover the original supplementary
examples; inactive optional dependencies are visible; observed and inferred
Cas components cannot be confused; rare/derived systems can remain unassigned;
mixed systems and fragment boundaries are represented. Full original parity
requires P1–P7 and an item-by-item audit, including verified alternatives for
historical functionality. Scientific basis: E2, E14–E16.

## P8 — Conditional later additions

Consider subtype-aware PAM/flank extraction and motif summaries after P3/P5.
Require correct protospacer orientation, adequate independent hits, background
controls and an explicit distinction between observed flanks, inferred motifs
and experimentally supported PAMs. Do not apply one DNA-PAM rule to all CRISPR
systems. Numeric sample-size/score thresholds require a benchmark, not a guess.
Flank enrichment may reflect acquisition and survival biases; it is not an
experimentally established interference PAM. Retain motif uncertainty and
abstain when support is sparse or orientation unresolved.

Defer automatic anti-CRISPR discovery, new effector discovery, broad defense
island mining, guide design and unrestricted metagenomic-scale search. They need
separate scientific objectives and validation and do not close the immediate
parity gaps. Preserve existing tracrRNA II/V-K coverage and evolutionary analysis
without claiming that they cover all system types. Scientific basis: E15–E17.

## Architecture and release gates applying to every package

- The API validates bounded input and queues jobs. Only isolated workers invoke
  operator-owned, shell-free commands. Preserve deadlines, containment,
  cancellation, retention, quotas and atomic queue transactions.
- Add focused scientific adapters; do not distribute stage-order ownership or
  place scientific models in the API/frontend. Reuse validated outputs within
  a job rather than rerunning the five existing tools unnecessarily.
- Keep backend validation independent of Integration validation. Coordinate,
  sequence, count, source and checksum checks must hold on both sides.
- Propagate contract changes through Integration, backend publication and
  summaries, frontend models/views, fixtures, documentation, wheels and release
  identities. Keep SpacerPlacer and evOr in separate environments.
- Safe artifacts must remain allowlisted, bounded, input-bound and free of
  internal paths or credentials. Generated HTML/pickle cannot bypass the
  publication boundary. Keep `improvenments/` out of public exports.
- No automatic tool/database update during a request or running job. Any change
  to independent scientific checkouts is an explicit upstream/release task.
- Preserve standalone reproducibility: every new scientific workflow needs a
  documented local CLI route with equivalent effective settings, normalized
  inputs, references and result schema, without a live web service or any
  predecessor project directory. Do not claim CLI parity until verified.

For future implementation, run `bash scripts/check.sh` for relevant changes.
It runs Python/Integration tests, lint, frontend tests and a production build;
it does **not** run the browser scripts. Add/update and run
`scripts/browser_acceptance.py` and `scripts/browser_production.py` for changed
UI flows in a controlled local setup. Extend `scripts/gate_release.py` with real
science cases for each new input/mode and negative/uncertain outcomes. Do not
reuse old release-gate results as evidence for a changed scientific release.

Scientific evaluation must distinguish software correctness from biological
validity. Assemble curated positives, hard negatives, species/family-held-out
examples, reverse-complement and fragmented variants, and database-leakage
checks. Freeze datasets, accession versions, labels, parameters and acceptance
targets before the final run. Report per-feature coverage, errors, abstentions,
resource usage and uncertainty; avoid one misleading aggregate accuracy score.
Publication benchmarks and agreement with the legacy server are baselines,
not ground truth. Production deployment remains a separate explicit request.

Proposed ownership: scientific maintainers own method/reference selection and
biological benchmarks; Integration owns adapters and stage contracts; backend
owns independent validation and safe publication; frontend owns evidence
presentation; release maintainers own pinned environments and gate evidence.
P0 should assign actual maintainers and measure storage/runtime budgets. No
delivery dates or numerical accuracy promises are justified before that work.
