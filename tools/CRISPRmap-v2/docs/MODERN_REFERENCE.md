# Modern repeat reference

The modern catalogue pipeline imports frozen ClusterCRISPRCasDB campaign
observations. It is independent of the historical CRISPRmap reference and its
family/motif vocabularies. The pipeline is implemented for release 2.1.0a1;
executed release statistics and validation belong in the benchmark release
report, not in this specification.

## Frozen source and evidence

The selected campaign starts from the NCBI GenBank bacterial/archaeal catalogue
frozen on 2 September 2026. Its authoritative outcome inventory contains
3,331,684 complete assemblies, 151 failed assemblies and 12 unavailable ones.
The complete campaign reports 9,583,665 arrays and 55,023,435 observed repeats.
These are source totals, not an assertion that every row passed this importer.
The build must reconcile selected outcomes and scientific record totals before
its reference release is accepted.

The scientific profile SHA-256 is
`5f24613d42b75c8477ec8f4988355e428b848c8a3f8f75078dca9d6d2a74e540`.
Arrays were predicted with CRISPRidentify 2.3.2, accepted-only region-screen,
without strand inference or embedded Cas analysis. Independent CasAndra calls
have no authoritative array-to-Cas linkage and are not imported as repeat labels.
The catalogue represents computational observations, not curated ground truth.

Each source archive is checksum pinned. The importer independently validates
its database checksum, schema, scientific profile, collection index and selected
assembly identities. Private source paths are input configuration only. Public
provenance contains hashes and accession identities; it does not contain host
addresses, SSH routing, connection records or raw metadata with contact details.
Only allowlisted taxonomy fields are retained from original assembly metadata.
The software MIT license does not relicense upstream data; preserve source
attribution and review upstream terms before redistributing a data release.

## Schema and normalization

`assemblies`, `contigs`, `arrays` and `repeat_occurrences` preserve occurrence
origins. Accession versions and stems are both retained. Only contigs supporting
arrays are included, with their original sequence checksums when present.
Failed and unavailable assemblies remain explicit outcomes with no invented
scientific values. Arrays require a complete, source-validated owning assembly;
an invalid source transaction is rolled back. The benchmark independently
checks assembly eligibility rather than trusting global reference counters.

`nodes` deduplicates exact canonical sequences after uppercase and U-to-T
normalization. A sequence and its reverse complement share a node. Occurrences
retain the orientation relative to that canonical representation, with `0` for
palindromes. This orientation is not a transcription-strand prediction. Source
strand remains null when absent. Non-normalized raw strings are stored when
needed to reconstruct source values. Missing sequences remain null; ambiguous
or out-of-policy sequences remain represented but cannot enter ordinary
unambiguous sequence search.

Four separate counts distinguish Bona-fide consensus, Possible consensus,
Bona-fide observed repeats and Possible observed repeats. An observed variant
is not silently promoted into a new consensus. Neither frequency nor a detector
score constitutes calibrated confidence in a family or functional repeat.

## Operational sequence groups

Core nodes have an eligible ACGT sequence and at least one Bona-fide consensus
occurrence. Only these nodes can form connected components. For threshold `s`,
an edge exists when the smaller full-length unit edit distance over the two
orientations obeys `distance <= floor((1-s) * max(lengths))`. All comparisons
use exact rational thresholds. Default sensitivity thresholds are 0.85, 0.90
and 0.95; these are fixed operational policies, not fitted biological boundaries.

The partition index is lossless within its declared range. Splitting a reference
into `k+1` nonempty segments ensures that at most `k` edits leave at least one
segment intact. Searching every feasible shifted position and both orientations
retrieves a superset of true neighbors, which are then exactly verified. If a
budget is exceeded, the operation fails instead of publishing partial success.

Candidate-index schema 2.0.0 stores the core-membership flag in the posting key,
so graph queries seek core postings without scanning variant-only postings.
Deep validation checks every flag against the bound catalogue. Indexes built
with schema 1.0.0 must be rebuilt with `catalog-index`; they are rejected rather
than silently interpreted using the new layout. The catalogue schema is unchanged.

Parallel workers own deterministic query lanes. Each retains a verified
spanning forest for each threshold. Combining every lane produces the same
components as the complete qualifying graph. Stored forests are not full edge
lists: they do not establish graph degree or full-graph diameter. Membership
hashes and policy identify the new `MAP2-S…` groups. IDs are independent of
worker count and node insertion order; they can change when membership changes.
They are not historical `F` family IDs or validated biological families.

Connected components permit chaining. Reports count members outside the
threshold from a deterministic lexicographic representative and give component
sizes and representative-distance histograms. A component does not promise
all-pairs similarity above its threshold.

Possible-only and variant-only nodes remain searchable. With cluster evidence
requested, the mapper performs a separate complete core search and retains all
qualifying operational cluster candidates. Such nodes cannot bridge components.
Nearest reference matches, direct core membership and threshold-based cluster
candidates are different report fields.

## Structural evidence

Optional ViennaRNA 2.7.0 folding computes both orientation hypotheses at 37 °C
with dangles=2 and the engine's default Turner 2004 parameters. Evidence includes
dot-bracket MFE structure, MFE and ensemble free energies, MFE probability,
ensemble mean base-pair distance and normalized ensemble defect to the MFE
structure. A lane manifest states exact coverage; completion of one lane is not
completion of the full core set.

No stability-based strand choice is made. These predictions do not establish
shared structural motifs, covariation, processing or cellular function. Stronger
motif claims require independent structure-aware alignment and validation.

Implementation reference: [ViennaRNA Python API](https://viennarna.readthedocs.io/en/latest/api_python.html).

## Evaluation boundaries

The modern benchmark separates genome-held-out coverage with measured overlap,
genome-held-out coverage after complete exact/RC sequence purging (including
all observed variants), and sample-scoped homology-purged novelty stress.
Contig checksum overlap and unknown checksums are explicitly handled. Taxonomy
availability is reported; an accession split is not a species-disjoint split.
Detector training overlap remains unknown unless training manifests are supplied.

At a 90% retrieval cutoff, removing every training sequence at 90% similarity
to a test query forces zero coverage for that query. The homology purge is an
abstention stress check, not evidence of family generalization. New clusters
computed on a complete catalogue do not provide independent held-out labels.
Modern family/motif biological accuracy remains unavailable without independent
labels. Historical output concordance remains a separate benchmark and cannot
establish original-executable runtime or speedup.
