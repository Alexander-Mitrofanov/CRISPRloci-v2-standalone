# Output, categories, and database contract

Every completed sequence result has one canonical `report.json`. The
`machine` output profile keeps this file and suppresses presentation
artifacts. Standard output can also contain text, CSV, FASTA, HTML, and
optional CasAndra tables.

## Coordinates

The report declares its coordinate system. Public non-empty intervals are
1-based and end-inclusive on the input forward strand. Orientation is a
separate field. Empty spacers use a 0-based interbase boundary and never
pretend to be a one-base interval.

Always read `coordinate_system` and `source.sha256`; do not infer
coordinates from file names.

## Array categories

| Category | Database interpretation |
|---|---|
| `Bona-fide` | Highest CRISPRidentify confidence category |
| `Possible` | Accepted candidate retained for high-recall screening |
| `Alternative` | Alternate representation associated with a candidate region, not automatically a distinct biological locus |
| `Possible Discarded` | Non-accepted candidate |
| `Low score` | Non-accepted candidate |

The canonical accepted set is `Bona-fide + Possible`, which is also reported
as `summary.accepted_array_count`. Use only `Bona-fide` when a smaller,
higher-confidence set is required. Do not compare exact representation counts
as if they were biological-region counts: overlapping candidates can be
different boundaries or repeat phases for one region.

Extract accepted regions as TSV:

```bash
jq -r ".source.id as \$source | .arrays[] | select(.category == \"Bona-fide\" or .category == \"Possible\") | [\$source, .source_interval.start, .source_interval.end, .category, .certainty_score, .id] | @tsv" report.json
```

## Unaccepted competing terminal extents

Conservative reconstruction reports can expose an existing longer **Possible
Discarded** observation beside an accepted **Possible** array. This is a
disclosure of a competing path, with its own score, rather than another accepted
array or an accepted extension. It does not change `arrays`, accepted counts,
classification, scores or reconstruction proofs.

The optional `extensions.competing_terminal_extents` uses schema
`CRISPRidentify-competing-terminal-extents` version `1.0.0`. Its coordinates are
**zero-based, half-open**. Rule `exact_terminal_competing_extent_v1` requires
the same final search key and an exact contiguous copy of the entire accepted
repeat-interval path within the longer observation. Additional repeats can be
on either flank. A search key is routing context, not a biological locus or a
limit on the observed source extent. Shifted boundaries, internal insertions,
cross-key candidates and structurally failed or ineligible observations are
outside this narrow view.

The observation registry preserves each linked final row's own score, features,
alignment, source geometry and available provenance once. Links identify the
final accepted array and complete competing path. Distinct aliases and repeated
inventory locations remain inspectable; they are not independent support votes.
Scores are `existing_final_row`: no classifier rerun or feature-refresh
certification is implied. Missing negative assessments remain `not_recorded`.

The diagnostic `routing_exclusions` describes the Bona-fide/Alternative
reconstruction lane only. It does not explain why a path remains after the
separate Possible reconstruction pass. Actual fresh Possible checks and
abstentions are recorded in `possible_reconstruction_evidence.json`; the
disclosure does not infer them from a final category or stored score.

An `unaccepted_observation` can still have unassessed biological specificity.
A `blocked_observation` records a known negative, original-occurrence or scope
conflict, or overlap with another accepted array. Collision checks include the
full union extent of a boundary family. Blocked observations remain evidence,
not recommended extensions. Malformed or stale source, score, provenance,
membership or accepted-array bindings prevent publication.

The full disclosure remains in `report.json` when optional candidate retention
is off, including accepted-only output. Direct APIs without a complete supplied
inventory cannot claim that a complete search found no competing extents.
Publication checks correspondence to the supplied final inventory, not every
path generated earlier or biological truth.

HTML and text label these paths **Unaccepted competing terminal extent** and
separate accepted from observed repeat counts and own scores.
`Competing_terminal_extents.csv` is a companion table; its rows are linked
observations, not accepted arrays. The primary `Summary.csv`, accepted FASTA
headers and sequences, and legacy JSON/pickle scientific contents are unchanged.
The companion is published transactionally with `report.json` last and removed
when absent or when the machine output profile suppresses presentation files.

## Candidate research evidence

`--retain_candidate_evidence True` writes a deterministic gzip JSON-lines
sidecar, `candidate_evidence.jsonl.gz`, for each completed input sequence.
It is opt-in and does not change candidate generation, scoring, category
assignment, canonical `report.json`, accepted counts, text/CSV/HTML array
lists, or annotations. It is an output observation rather than a scientific
parameter. Its configuration is recorded in the sidecar and in the optional
`research_outputs` part of `--print-effective-config`.

The inventory is captured after final negative/circular decisions and before
output projection. It therefore retains Bona-fide, Alternative, Possible,
Possible Discarded and Low score even with `--accepted_only True`. This is
the final inventory, not every intermediate geometry ever generated.

The `CRISPRidentify-candidate-evidence` schema is version `1.0.0`:

- The first `header` binds the inventory to source ID, sequence SHA-256,
  length, effective scientific configuration and capture stage.
- `comparison_scope` records group rows by their original search context
  and link their accepted candidate IDs. Such a scope is **not a proven
  biological locus**, and scope membership does not establish that candidates
  can safely be merged.
- `candidate` records contain final category, genuine existing score and
  feature values, eligibility, repeat intervals, gapped repeats, fresh
  structural-filter evidence, compact source-path lineage, explicit available
  negative evidence, and links to canonical accepted IDs where applicable.
  Coordinates are **zero-based, end-exclusive**, as declared by the header.
  Repeat and spacer sequences can be recovered from the bound source sequence.
- `source_geometry_id` identifies the source hash and ordered repeat
  intervals. `alignment_id` identifies the gapped alignment. `candidate_id`
  identifies the scored representation/category/search context; identical
  observations can share it. `observation_id` distinguishes inventory rows.
  Multiple alignments or seed histories must not be counted as independent
  arrays or biological support.
- The final `footer` records completion, category/observation counts and a
  SHA-256 over every preceding uncompressed JSON line including its newline.

`source_valid` certifies sequence reconstruction, not biological confidence.
`structural_filter.decision` gives the full conjunction result; its
`structural_valid` field alone only checks basic candidate structure.
Negative evidence marked `not_recorded` is unknown, not clearance. Scores and
features are retained unchanged from final rows; the writer does not rerun
the classifier or independently certify score freshness. The optional
feature-repeat-count consistency field is a diagnostic, not a score repair.
Unaccepted candidates do not inherit orientation, Cas or IS annotations from
an accepted neighbor. Full recursive family provenance is summarized to avoid
duplicating large histories; the header declares that limitation.

The sidecar is staged in the same rollback-capable transaction as the report;
`report.json` is published last. A run with retention disabled removes a stale
sidecar from a reused output directory. A failed run must not be interpreted
as a completed empty inventory. Research consumers should require successful
run status, the current canonical report, a valid gzip/footer and matching
source binding. Never append these rows to the ordinary accepted-array list.

## Core-supported reconstruction

`--reconstruction_policy conservative` is the default policy for `balanced`,
including when the preset and policy flags are omitted. It is mutually
exclusive with architecture rescue. Conservative mode
preserves ordinary candidate competition while restricting repair proposals
that would join distinct ordinary groups. Its later reconstruction stage runs
after ordinary enhancement/repair and before final tandem quarantine. The
explicit `legacy` opt-out adds no reconstruction stage and keeps its existing
competition behavior. `full`, `region-screen`, and runs explicitly selecting
architecture rescue resolve to `legacy` by default; an explicit incompatible
`conservative` request is rejected. Both policies are recorded under
`effective_config.reconstruction_policy` and
`run.parameters.param_reconstruction_policy`, including legacy runs.

The reconstruction decision sidecar reports `experimental: false` for this
default implementation. Its separate `biological_accuracy_established` field
remains false: default status does not establish biological precision or recall.

Initial categorization keeps the configured ordinary candidate components:
strict source-envelope overlap components when component splitting is enabled,
or the existing unsplit comparison scope otherwise. It forms separate
components of terminal-repair-origin candidates, using zero-based, half-open
source-envelope overlap. A repair component touching one ordinary component
may compete within it; a repair component touching none competes independently.
A repair component touching several ordinary components is retained as
noncanonical conflict evidence and cannot merge them, even if a repair has a
higher model score. This includes indirect repair-to-repair chains. It can
defer a valid single-group repair when another connected repair touches a
second group; this is an explicit conservative tradeoff.

Compact, nonrecursive origin witnesses bind repeat coordinates, repeat/spacer
sequences and alignment. A terminal repair receives a repair birth bound to
its source hash and generated-child attestation; derived candidates inherit
origins, and deduplication preserves independent origins. A candidate also
observed through an ordinary origin remains ordinary. Missing or stale origin
bindings reject conservative competition. Ordinary groups and score bands
remain unchanged. These groups do not establish biological loci, and the rule
does not claim to detect all fused arrays in the ordinary inventory.

The optional `initial_competition` observation records `ranked_scope` or
`repair_scope_conflict`, the candidate geometry and score, the scope span,
ordinary and repair component descriptors, and the repair groups' ordinary
owner IDs. A conflict retains its own valid score and structural eligibility;
noncanonical retention is a competition decision, not a biological-invalidity
label. Full candidate sidecars also retain `candidate_origins` summaries.

Exact score ties keep the legacy stable input order. Distinct geometries
within 1e-12 of the winner's score appear in `tied_candidates`. Exact-geometry
aliases do not create a phase ambiguity. `unique_phase: false` records
multiple competing geometries; `null` does not certify a unique biological
phase. Tied geometries can be disjoint within a transitive ordinary component.
Observation validation checks internal consistency and source bounds, not the
complete earlier competitor inventory or independent biological correctness.

An initial tie does not automatically block reconstruction: one larger donor
may explain several tied short cores. A selected reconstruction preserves
historical observations as `anchor_initial_competition` and
`donor_initial_competition` in its proof, bound to the original scores and
geometries. The original core's observation is checked against retained core
evidence. These observations do not authorize selection; donor uniqueness is
conditional on the chosen core and does not establish biological resolution
of every initial competing parse.

Each Bona-fide anchor must retain a fresh model score of at least 0.75.
Alternative donors must have their own fresh score of at least 0.60, pass the
ordinary structural conjunction, reconstruct exactly from the source, retain
classification eligibility, and pass each applicable negative check. This
does not lower the ordinary classifier thresholds or independently admit an
Alternative as a new locus. Accepted array count is preserved by this stage.

Compatibility requires a unique ordered one-to-one mapping of every anchor
repeat into the donor: each boundary differs by at most three bases, with at
least 80% reciprocal repeat overlap. Added repeats require at least 80%
identity to one core-derived reference and adjacent spacers of at most 140 bases.
Terminal additions and internal insertions are recorded separately. There is
no synthesis of new geometry or transitive merging through overlapping spans.
Donors overlapping another Bona-fide or Possible array are rejected; competing
proposals that collide both abstain. Selection ranks qualifying geometries by
repeat count, then their own fresh model score. Distinct best geometries tied
within 1e-12 remain unresolved. When the entire tied set shares ordered repeat
starts, count and original-core mapping, the tool can report the common extent
as a boundary-uncertain family, as described below. Other ties abstain;
equivalent alignment aliases are recorded separately.

An ordinary reconstruction moves only its selected donor row into the accepted
output. Nonselected alignment or provenance observations of that geometry stay
as complete Alternative rows, with their original scores and evidence. They do
not add accepted loci or votes. The displaced original core also remains an
Alternative.

The original core consensus is tested first. If an added repeat fails that
identity test, a narrow phase fallback can use the donor's copies mapped to
the original core. All original core repeat strings must equal their original
consensus, every mapped donor copy must have the same start/end offset pair,
and all mapped strings must be identical canonical DNA. The existing 3 bp and
80% correspondence limits still apply. That shared mapped string supplies one
reference for every added repeat. Added copies and the donor's overall
consensus never determine this reference; each copy cannot choose whichever
reference gives the better score. An inconsistent core, mixed offsets or a
weak added copy leaves the donor ineligible under this fallback. The final
whole path still requires every other score, source, structure, occurrence,
negative and collision check.

Only relationships that need this fallback include `added_repeat_reference`,
with the versioned rule `mapped-identical-original-core-v1`. Its evidence binds
the original consensus, mapped original and donor intervals, offsets, shared
reference and the old/new edit distances and denominators for each added
repeat. Publication reconstructs this evidence from the source. A relationship
without this discriminator replays the original consensus rule, including
older reports; unknown or malformed rules are rejected. An original-reference
passing relationship retains its previous evidence shape. Each member of a
boundary-uncertain family is checked independently against its own original
core mapping.

The existing small-array tandem check applies to 2–8 spacers. The existing
coding-tandem check applies to equal in-frame period runs and obtains source
CDS predictions when applicable. An inapplicable check is explicitly distinct
from an applicable check that could not run; the latter blocks selection.
Known negative annotations also block selection. These scoped checks do not
establish general absence of tandem artifacts or biological authenticity.
Cas and IS support are not selection requirements.

The selected array remains Bona-fide through
`classification_basis.method = source_backed_candidate_reconstruction`.
Its `certainty_score` is its own freshly computed score, potentially below
0.75. The basis separately records the original and fresh core scores,
original donor and final features, source geometry, candidate checks, and
selection decision. The original core moves to Alternative with its original
score. Consumers must read the classification basis rather than infer the
category from score alone.

Publication replays source reconstruction, feature-to-model scores,
candidate-specific checks, the chosen direct relationship and collisions with
other accepted arrays. It binds original/fresh geometry and retained core
evidence. It does **not** independently reproduce the original full donor
ranking universe, rerun external CDS prediction, or repeat full feature
extraction. The producer performs batched fresh feature extraction before
selection; the report replay has the narrower integrity scope just described.

`reconstruction_evidence.json` (schema version `1.0.0`) records source SHA-256,
policy, counters and every selected, retained or abstained anchor decision.
It uses zero-based, half-open coordinates, is published transactionally with
the canonical report, and is removed on a subsequent legacy run in the same
output directory. Full candidate retention remains a separate opt-in sidecar.
The current 500-assembly comparison reuses the inspected development cohort;
it cannot establish independent biological precision or recall.

### Checked extent of a Possible core

The conservative policy runs a separate Possible reconstruction pass after
the unchanged Bona-fide pass. It considers existing Possible Discarded paths
around accepted Possible cores. It does not recover a disjoint Low score
candidate or change the initial 0.50/0.75 score bands. The original Possible
core must remain eligible under fresh checks; a proposed longer path needs
its own freshly computed score of at least 0.60 and every applicable source,
structural, negative and ordinary-occurrence check. Its original core must
map uniquely under the same three-base/80% limits. The existing added-repeat
identity and spacer limits apply.

`classification_basis.method = source_backed_possible_reconstruction` keeps
the replacement in Possible. The proof separates original/fresh core scores
from the full path's own score and records the actual retained donor and
checked relationship. Only the chosen observation moves to Possible; its
original core and all nonchosen donor observations remain Possible Discarded.
Alignment aliases do not provide extra votes. Equally ranked distinct
geometries abstain in this pass, including tied repeat-end alternatives.

Previously accepted Bona-fide extents, including the complete union span of
boundary-uncertain families, are collision barriers. A Possible proposal
cannot veto or merge a completed Bona-fide reconstruction. It also cannot
overlap another accepted Possible core; colliding new Possible proposals
both abstain. Accepted locus count is preserved by this pass.

`possible_reconstruction_evidence.json` records the source binding, policy,
decisions and counters separately from `reconstruction_evidence.json`.
It is published transactionally with the report, remains available when
candidate retention is disabled, and is removed if a later run does not
produce that audit. The classification basis and compact per-array summary
remain in the canonical report. As with Bona-fide reconstruction, publication
replays the declared source/model/proof checks; it does not independently
rerun full feature extraction, external CDS prediction or the original
entire donor search. Passing these scoped checks does not establish
biological authenticity or calibrate either score as a probability.

### Localized identical-spacer context

The conservative policy separately assesses source-valid generated paths with
one internal maximal block exceeding the ordinary consecutive identical-spacer
limit. It captures these paths before ordinary filtering and enhancement.
The duplicated spacer must occur only in that block, with no unrelated
identical-spacer burden. Both maximal flanks retain their boundary repeat and
must contain at least three repeats (or the configured larger minimum).

The complete path and both flanks are freshly scored on their own features and
must each reach 0.50. The full path must pass every structural predicate except
the specifically reviewed identical-spacer burden. Both flanks must pass the
ordinary structural filter. All three receive the applicable existing coding
and tandem checks. An inapplicable screen is recorded as such; it is not a
general negative result. Original eligibility, lineage, known negative
evidence, occurrence protection and accepted-array collision barriers remain
binding. Previously accepted arrays are preserved.

`classification_basis.method = source_backed_localized_spacer_context` caps
the accepted interpretation at **Possible**, including full-path scores above
0.75. The flanks are correlated diagnostic context, not independent votes or
substitute confidence scores. Selection ranks qualified paths by repeat count
and the whole path's own score; flank scores cannot resolve a tie. Every
equally ranked pending peer is then checked, including nonqualifiers. A
higher-ranked path that failed its checks remains in the audit but does not
veto a qualified lower-ranked path. Equally ranked safe
geometries with identical ordered starts can form one complete boundary
family. Mixed starts or an unsafe same-rank peer cause abstention. End choices
from different members cannot be recombined.

The compact `contextual_spacer_context` summary names the displayed actual
member, all member IDs, own-score range, union span and variable repeat
indexes. This family has no accepted short core. Every captured original
observation remains Low score evidence, including nonselected paths and
alignment aliases. The selected Possible view is additional and retains its
own score. The union span describes the alternatives' envelope, not a
synthetic array to export as sequence.

With FASTA reporting enabled, `Contextual_Arrays.fasta`,
`Contextual_Repeats.fasta` and `Contextual_Spacers.fasta` export the selected
actual member on the input forward strand. Headers identify the Possible
confidence ceiling, member, own score and boundary uncertainty. The existing
`Arrays.fasta`, `Repeats.fasta` and `Spacers.fasta` retain their Bona-fide-only
scope. Contextual FASTA files are absent in machine mode or when no contextual
array was accepted, and stale files are removed on a later run.

`contextual_spacer_evidence.json` records capture, assessment and selection
separately. It is validated and published with the canonical report even when
candidate retention is disabled, and removed on a later run with no audit.
Publication checks source, feature-to-model scores and transaction proofs;
it does not repeat feature extraction or external CDS prediction. This route
does not establish biological authenticity or population recall.

### Shared extent with uncertain repeat ends

`classification_basis.method = source_backed_boundary_uncertain_reconstruction`
identifies one accepted occurrence family with several admissible complete
geometries. Every member passes the existing source, structure, fresh own-score,
negative and original-core checks. An equally ranked direct peer that fails
the occurrence-preservation, collision or candidate checks prevents promotion;
filtering cannot silently turn a failed tie into a resolved result. Different
ordered starts or core mappings also preserve abstention. Lower-ranked failures
do not veto a qualified winner.

The array's `boundary_family` descriptor reports `status: boundary_uncertain`,
the shared repeat count, original core count, actual member IDs, display-member
ID, member score range and conservative union span. Its coordinates and variable
repeat indexes are explicitly **zero-based, half-open**, unlike the report's
public one-based inclusive `source_interval` fields. The full source-bound
proof is under `classification_basis.boundary_family` in mandatory `report.json`,
including when optional candidate retention is disabled.

The ordinary array coordinates, repeat/spacer sequences, consensus, features,
orientation and `certainty_score` describe one explicitly named **actual full
member**. This deterministic display choice does not settle the repeat ends.
The original core and member observations remain retained. Each member has its
own score and checks; the original core's confidence is separate. All members
must directly preserve the protected original ordinary occurrences carried by
the core and every tied member or alignment alias. Preserving only a member's
own witness is insufficient, and an intermediate geometry cannot transfer
authority. Collision checks use the family's union span. A family counts once
as an accepted array and once when counting added repeat occurrences.

Per-repeat observed end sets summarize the complete members. They never permit
mixing ends into a new segmentation: only the listed full geometries are
admissible. The existing boundary assessment marks repeat-phase uncertainty
with the reason `shared_ordered_starts_uncertain_repeat_ends`; its biological
occurrence-path assessment remains unassessed. Shared coordinates and exact
source reconstruction do not establish biological completeness or truth.

Publication checks the recorded same-rank members against retained original
observations and the emitted transaction decision. It does not independently
rescore every earlier donor or reconstruct the complete upstream search.
Consistently rewriting the entire provenance record is outside this integrity
check; the producer's complete-set selection and regression tests establish
the supported pipeline behavior.

Text and HTML identify uncertainty before showing the member's sequences.
FASTA headers label the family and actual member. Optional legacy JSON includes
the compact descriptor and a reference to the matching source's canonical
report; pickle retains the complete row metadata. `Summary.csv` appends eight
columns: `Geometry status`, `Boundary family ID`, `Display geometry ID`,
`Admissible geometry count`, `Supported repeat count`, `Original core repeat
count`, `Extent start` and `Extent end`. The extent columns are one-based
inclusive union coordinates; the existing Start/End columns still describe
the actual displayed member. These columns are present for ordinary and empty
components too, so folder aggregation keeps a consistent schema. An ordinary
row's `reported_geometry` status makes no biological resolution claim.

## Source-backed repeat-boundary repair

Before structural pruning, the detector may trim weak terminal columns from
a source-canonical repeat alignment when every other structural predicate
passes. Repeat and spacer sequences are rebuilt from the input sequence;
classifier and structural thresholds remain unchanged. Restricted lineage
parents cannot become eligible through boundary trimming.

A separate checked transaction can complete an accepted short core through
a unique boundary-repaired path with identical repeat copies, at least two
additional copies, ordinary individual spacer lengths, and no second
overlapping accepted event. The longer candidate is scored independently
and must meet the incumbent category's threshold: 0.75 for `Bona-fide`, 0.50
for `Possible`. Its actual score determines its final category. Relative
score ranking alone need not select this more complete representation.

The classification basis `source_backed_terminal_repair_completion` records
the original core, parent and replacement geometries, genuine old/new scores,
source identity, repair evidence and transaction guards. Run-level repair
attempts, abstentions and bounded-search counters are available under
`run.parameters.param_repeat_boundary_repair`. The displaced core remains
lower-category evidence internally; presentation profiles do not necessarily
serialize all five categories. These records justify a source reconstruction,
not a claim of biological validity.

## Source-backed internal repeat completion

An anomalous spacer can hide omitted repeat copies. The bounded source search
holds the parent's existing anchors fixed and tests every allowed interval
against its dominant source motif using the unchanged fuzzy edit limits.
Exact analysis of the complete interval graph must establish one inserted-copy
count and one minimum-cost geometry in every anomalous gap. Equal-cost boundary
alternatives, unresolved gaps and exhausted search budgets cause abstention.
A bounded diagnostic witness list is explicitly separate from graph completeness;
exact path counts are serialized as decimal strings.

Derived candidates remain pending until a checked transaction can replace one
accepted core on the same original repeat path. Every original restriction is
replayed; only an attested broken-edge restriction can be discharged for the
new geometry. Fresh scores must meet the incumbent's existing category band,
all ordinary structural filters still apply, and another overlapping accepted
event anywhere in the source blocks replacement. Independently attested stage-1
terminal trimming may compose with this internal operation.

When a core needs a separate boundary witness, the complete available family
of overlapping equal-count stage-1 children, current rows and reattested pending
paths participates in the existing source-phase comparison. More than 32 distinct
geometries, contrary restrictions, incomplete evidence or an unresolved phase
blocks this route. A successful witness authorizes only the exact corrected
local core map to retained parent anchors; it transfers no whole-path geometry,
eligibility or score. Report replay checks the family against the independently
retained run audits. The completed geometry may already exist as an internal
Alternative; replacement does not imply discovery of a new biological array.
The internal run audit includes `current_authority_inventory`, captured before
replacement from every current category row overlapping a pending target.
It preserves duplicate derivations and absent, malformed or explicit authority
attributes. Both completed-array and retained-core replay require the exact
relevant row census; missing or changed members invalidate the proof.

The `source_backed_internal_repeat_completion` classification basis records
source, parent, retained core, inserted copies, complete search evidence and
actual scores. A displaced core keeps its original score and is retained as
lower-category evidence; its basis identifies that retained role when needed.
`run.parameters.param_internal_repeat_completion` records proposal abstentions,
actual search work, cache hits and transactions. Complete source reconstruction
is not biological validation, and an ambiguous occurrence count does not
certify exact repeat boundaries or improved accepted coverage.

## Likely tandem quarantine

Before output projection, every final source-reconstructed `Bona-fide` and
`Possible` row is screened for the narrow low-spacer-diversity phenotype. A
row receives the disposition `likely_tandem` only when it has 2–8 non-empty
canonical A/C/G/T spacers and **every** spacer pair has full-length global
Levenshtein identity of at least `0.95`, where identity is
`1 - distance / max(lengths)`. This is complete-link similarity: one similar
pair or a chain of near neighbours is insufficient. Individual spacers are
not reverse-complemented for comparison.

The decision uses no external-tool call, model-score override, or benchmark
label. A qualifying row moves to the existing public `Low score` category and
receives the sorted array tag `tandem`; `likely_tandem` remains its scientific
disposition in classification provenance. It is therefore absent from
`Possible` and from `summary.accepted_array_count`, but remains visible in a
full canonical report and `Low_Score_Candidates.txt`. Optional HTML,
legacy-JSON, and pickle writers preserve the same tag when enabled. Its
separate source-hash-bound record at
`extensions.likely_tandem` binds one-to-one to that Low-score array and
preserves the original category and score, source interval, spacer sequences,
pairwise evidence, policy, and transaction counts. With
`--accepted_only True`, the row is projected out and the ledger explicitly
records that omission; current canonical output requires this ledger even
when no tandem rows exist. For a non-empty ledger, the writer replays each
record against the original FASTA, verifies the whole-source hash, and
requires an exact repeat/spacer/core match to its report row. The writer also
reruns the transaction across every pre-quarantine accepted row and requires
the complete positive ledger and exact transaction counts. HTML rendering and
benchmark ingestion independently rescreen accepted core architectures and
replay positive evidence, so a self-consistent but incomplete ledger is
rejected. `report.html`
shows the same evidence in a **Tandem demotion audit** table, including for
accepted-only output where the quarantined core row is intentionally absent.

`likely_tandem` is an auditable soft quarantine, not a claim of proven
biological falsehood. Authentic CRISPR arrays can duplicate spacers; therefore
an array is never filtered merely because any one spacer pair is identical.
The eight-spacer maximum is a bounded small-array policy matching the screened
scope, not a biological upper bound on tandem-repeat evidence.

The pipeline removes a stale `report.json` before its first fallible stage and
publishes the replacement strictly last, after all requested artifacts have
validated. Its presence therefore denotes a completed per-sequence result,
not merely a prior successful invocation in the same directory.

For a simple within-tool region inventory, sort by source and start, then
collapse overlapping accepted intervals on the same source. For comparative
recall, freeze a matching rule before looking at results. Specify target
categories, how multiple representations form a region, coordinate conventions,
overlap thresholds and failure denominators. Cross-tool agreement does not
establish biological truth, and locus overlap does not establish repeat or
spacer completeness.

## Stable database keys

The human-readable array `id` includes display order and can change when a
different preset adds, removes, or reorders representations. Do not use it as
a cross-run primary key.

A database row should retain at least:

- tool release and report schema version;
- effective preset and configuration SHA-256;
- `source.id`, `source.sha256`, and source length;
- interval start/end, category, orientation, and consensus repeat;
- the original human-readable array ID for traceability.

Read the named preset from `run.configuration.preset` and its digest from
`run.configuration.preset_provenance.effective_config_sha256`. The preset
digest covers the complete effective scientific configuration, including
explicit options outside the preset-controlled speed path. Adding custom
detector arguments therefore creates a distinct, recorded configuration; it
must not inherit the release benchmark claim without separate validation.

For a run-local deterministic key, hash a canonical serialization of source
SHA-256, interval, consensus repeat, category, and effective configuration
SHA-256. Keep representations separate from a region table so several arrays
can point to one collapsed biological region.

## CasAndra

`--cas True` annotates Cas genes and systems with the bundled, integrity
checked CasAndra provider. The JSON contract is under `cas`; standard output
also writes `Cas_systems.csv` and `Cas_genes.csv`. Check `cas.status`
before consuming systems:

- `computed`: annotation completed and passed adapter validation;
- `not_requested`: Cas annotation was disabled;
- an explicit failure status or process error must not be interpreted as an
  empty biological result.

Cas classification contains provider, model-bundle, route, evidence, and
hierarchy-normalization provenance. Type VII type-only calls are experimental
and are labeled as such. CasAndra is opt-in and should normally be a second
stage after region screening at multi-million-genome scale.

## JSON schema and version policy

The current report contract is `CRISPRidentify-report` version `1.2.0`.
Validate it with
[`schemas/report-1.2.0.schema.json`](../schemas/report-1.2.0.schema.json).
The runtime also performs stronger source-slice and provenance validation that
a structural JSON Schema cannot express.

Within report schema major version 1:

- patch changes clarify validation without changing valid structure;
- minor changes can add optional fields;
- removing a field, changing coordinates, or changing field meaning requires a
  new major schema version.

Consumers should reject an unknown schema name or major version and preserve
unknown optional fields.
