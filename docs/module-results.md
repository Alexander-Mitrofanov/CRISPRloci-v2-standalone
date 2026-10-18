# Compact results for every module

Backend/frontend 2.12.0 publishes `result.json` for all 20 analysis modes and
places it first among individual downloads. The shared, versioned envelope is:

```json
{
  "schema": {"name": "crisprloci.module-result", "version": "1.0.0"},
  "module": "cas",
  "status": "completed",
  "complete": true,
  "data": {
    "coordinate_system": "1-based-inclusive-source-forward",
    "score_semantics": "uncalibrated-model-evidence",
    "evaluated_gene_count": 10,
    "genes": [],
    "cassettes": []
  }
}
```

This example means ten genes were evaluated and no Cas genes or cassettes were
reported. Missing values stay null, genuine zeroes stay zero, and empty findings
stay empty. Unrequested modules, UI state, execution paths, command lines, logs,
runtime measurements and redundant row counts do not belong in this file.

The [JSON Schema](schemas/module-result-v1.json) describes the common envelope
and each module's data collections. It deliberately leaves individual scientific
records open: the field meanings below and independently validated native
contracts remain authoritative. It is not a replacement scientific validator.
Consumers dispatch on the schema version and `module`, and check `complete`
before reading `data`.

## Module data

The linked examples are generated from repository test fixtures. Evolutionary
examples use the synthetic worker, and the locus example composes fixture
array/flank evidence with empty Cas/tracrRNA results; these are contract examples,
not new biological analyses.

| Module and example | Compact scientific contents |
| --- | --- |
| [detection](examples/module-results/detection.json) | Accepted arrays with categories, raw scores, source identities, repeat/spacer occurrences and boundary uncertainty |
| [cas](examples/module-results/cas.json) | Evaluated gene count, Cas gene intervals/classifications and cassette classifications; no empty tracrRNA/detection sections |
| [tracrrna](examples/module-results/tracrrna.json) | Sources, candidate sequences/intervals, evidence origin, boundary kind, model/ranking scores and anti-repeat/interaction evidence |
| [loci](examples/module-results/loci.json) | Detection, Cas, tracrRNA and flank components with one shared source table |
| [leader](examples/module-results/leader.json) | Arrays and both adjacent flank sequences/segments; `prediction_performed: false` explicitly preserves the absence of a leader prediction |
| [repeats](examples/module-results/repeats.json) | Repeat occurrences referring to shared isolated folds, RNA sequences, structures, energies, support measures and unsupported-classification reasons |
| [repeat_context](examples/module-results/repeat_context.json) | The same fold data plus array/instance/context links, consensus folds and local support changes; analysis windows are not biological boundaries |
| [protospacer](examples/module-results/protospacer.json) | Target identities, a shared observed-spacer table, exact full-length matches, relative strands, reference identity and search completeness |
| [viral_search](examples/module-results/viral_search.json) | Query occurrences, a shared subject table, accepted matches, substitutions, scores, relative strands and heuristic retrieval/saturation limits |
| [spacer_association](examples/module-results/spacer_association.json) | Explicit groups, occurrence origins, target/control panel identity, calibration state and score populations, native combined scores and source-bound alignment support |
| [array_compare](examples/module-results/array_compare.json) | Observed ordered spacer occurrences, deletion slots, scoped source identities, explicit equality policy, sharing edges, components and isolates |
| [repeat_type](examples/module-results/repeat_type.json) | Occurrence-linked repeat subtype evidence, explicit unknowns, native model scores, threshold and immutable model vocabulary/provenance |
| [repeat_map](examples/module-results/repeat_map.json) | Query results, a shared reference-label table, tied alignments/distances, total best-hit counts, hit caps and conflicting/incomplete labels |
| [reconstruction](examples/module-results/reconstruction.json) | Group members, an input-order tree per group, estimated rates/likelihoods/event totals, correction statuses and excluded arrays |
| [orientation](examples/module-results/orientation.json) | Group members, input/reverse likelihoods and trees, decision threshold, selected reconstruction metrics, correction statuses and exclusions |

| [locus_compare](examples/module-results/locus_compare.json) | Verified CDS/protein neighborhoods, native protein alignments, source hashes and separate Cas-coordinate context |
| [prophage_context](examples/module-results/prophage_context.json) | Source-derived spacer matches, array exclusion evidence, candidate prophage regions and explicit assessed-source scope |
| [insertion_context](examples/module-results/insertion_context.json) | IS intervals, transposase ORFs, native terminal-repeat evidence and coordinate-only Cas proximity |
| [defense_context](examples/module-results/defense_context.json) | Pinned model coverage, defense systems, components, isolated profile hits, translated genes and separate Cas context |
| [rna_family](examples/module-results/rna_family.json) | Original-to-representative crosswalk, alignment, column counts, consensus structure, probabilities and descriptive covariation |


Array IDs are scoped to their source. Use `(source_id, array_id)` when joining
array, flank, gene or other locus evidence. Fold IDs refer to entries in `folds`;
reference IDs and subject accessions refer to their corresponding shared tables.
Viral query IDs identify occurrences, so repeated queries are not merged merely
because they have the same sequence. Matching queries include their sequence;
all original queries are available in `input-sequences.json`.

Coordinates are declared per component. Detection, annotations and search
intervals use 1-based inclusive coordinates. Flanks and repeat structure use
zero-based half-open coordinates. Sequence frames and relative strand remain
explicit; consumers must convert deliberately when joining different components.
Flank sequences follow the reported strand/segment semantics. Fold RNA is the
instance's analyzed orientation, not necessarily the original DNA strand.

For orientation, `input` and `reverse` refer to acquisition order, not genomic
transcription direction. An `ND` decision keeps the input-order reconstruction
with `selection_basis: "unresolved_default"`; this is not a positive orientation
call. Both independently estimated trees remain available. Reconstruction event
counts and ancestral histories are model estimates, not observations.

## Compactness and detail

Shared metadata is recorded once. Repeated sequences in fold/subject/reference
records are replaced by identity references, while distinct scientific
occurrences remain separate. Dense pair-probability/accessibility arrays, native
feature dumps, ancestral branch-event matrices and raw search HSP bookkeeping
stay in the existing detailed artifacts and ZIP. RNA-family is an explicit exception: its compact projection retains every reported pair probability and alignment-column evidence, including bounded pair-type count objects. Effective interpretive policies
remain in compact data when needed to understand results; software versions,
execution options and citations remain in provenance.

`complete: true` means the defined compact projection is fully exported. It does
not mean every native detail is included, that model evidence is certain, or
that heuristic searches are exhaustive. In particular, `hits_truncated` and
`best_hit_count` remain explicit when a search's reporting policy caps ties;
viral `retrieval_completeness` remains `heuristic`. Unknown overlap and unevaluated
PAM/host inference are retained in search policies. No rows are silently dropped
to fit a browser preview or export limit.

## Limits and module connections

The isolated worker publishes compact results after independent scientific
validation, using full evidence rather than capped UI summaries. Source reads
reject symlinks, duplicate JSON keys and non-finite JSON values, with a 20 MB
per-file and 40 MB combined read budget. Publication is atomic and bounded by
the smaller of 20 MB and the configured artifact limit, with cancellation
checkpoints, private-path screening and the existing per-job download capability.

If evidence cannot be safely projected or the compact result exceeds those
limits, `result.json` contains `complete: false` and
`reason: "compact_export_unavailable"`, omits `data`, and the job receives a
warning. Even this small status file is omitted if the artifact limit or storage
prevents writing it. Failed/cancelled jobs do not generate new compact results.
Successful evolutionary jobs with no eligible groups keep their explicit status,
empty groups and any excluded-array reasons.

Use the [typed sequence ports](module-exchange.md) for supported job-to-job
handoffs. A module-result file describes findings; it is not automatically a
valid sequence input or an executable workflow. `input-sequences.json` accompanies
successful modes, and detection-containing modes also supply `arrays.json`,
`repeats.json` and `spacers.json`. Existing detailed filenames and scientific tool
contracts remain available for consumers that need more evidence.

`array_compare` exports [complete observed comparisons](examples/module-results/array_compare.json):
ordered occurrences and deletion slots, scoped original IDs, source annotations,
explicit equality policy, thresholded sharing edges, components and isolates.
It declares one-based inclusive source-forward coordinates and uncalibrated native model-score semantics. Uploaded source/category/score/boundary annotations and claimed producers remain client-supplied and unverified. Explicit nulls and zero scores remain distinct. It does not infer an alignment or evolutionary tree.

The context examples are projections of native scientific fixtures, not deployment records. RNA-family and locus-comparison fixtures are synthetic contract examples; their structures or alignments are not biological claims. The insertion example uses ISEScan’s upstream NC_012624 fixture, the defense example uses the documented 40 kb upstream DefenseFinder fixture, and the prophage example uses a verified bacteriophage lambda source. Cas-context tables can be empty when the accompanying Cas fixture supplies no cassettes; an empty proximity list does not establish biological absence.

tracrRNA interaction energy is the native `energy_kcal_mol` value, exported as `interaction_energy` in kcal/mol. Optional missing values remain null; zero remains zero. Model-only absent evidence labels default to covariance-model/CM-only/merged-model-region provenance, while explicit nulls and complete-mode labels are retained.
