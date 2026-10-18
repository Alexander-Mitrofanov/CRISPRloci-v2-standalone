# Scientific interpretation

CRISPRmap-v2 provides two related workflows for prokaryotic CRISPR repeat
analysis. The historical `classify` command maps queries to a frozen reference
with its original family and motif annotations. The modern catalogue extension
imports provenance-bearing detector observations, searches consensus and
observed repeat sequences, constructs operational sequence groups, and can
report separate RNA folding hypotheses. Neither workflow detects arrays from
genomes or assigns modern Cas-system subtypes.

The modern pipeline is implemented in 2.1.0a1. Its methods are described here;
completion and measured results of each full-corpus stage require that stage's
manifest and executed release report.

## Different forms of evidence

| Question | Historical reference | Modern catalogue |
| --- | --- | --- |
| Is this sequence in the reference? | Exact membership with the corresponding records and provenance. | Exact canonical membership with retained array and repeat occurrences. |
| Which sequence families or groups resemble it? | Candidate historical family associations, including disagreement and missing labels. | Operational sequence-component candidates in a separate namespace; these are not validated biological families. |
| Which structural motif does it form? | Nonexact motif identity is not inferred; original reference annotations remain attributed to their records. | Optional thermodynamic folding hypotheses for both orientations, without motif IDs or experimental structure claims. |
| Which Cas subtype, organism or transcription strand does it establish? | None. | Source accessions and available taxonomy describe observation origins; they do not assign the query a host, subtype or transcription strand. |

In the historical workflow, an exact match can carry the same annotation as a
reference record, but does not validate that annotation experimentally. Several records may share
a sequence and disagree in labels or origins. Null annotations mean that the
source did not provide the annotation; they do not establish biological absence.

## Sequence comparisons and uncertainty

Full-length edit distance counts substitutions, insertions and deletions needed
to transform one complete sequence into another. It differs from a local
alignment that can ignore unmatched ends. Each edit has cost one. The shared
similarity definition is `1 - distance / max(query_length, reference_length)`;
it is not a confidence score.

The historical mapper filters and ranks by raw edit count, returning all
closest ties before a display cap. The modern catalogue search also separates
nearest-reference hits from operational cluster evidence. Its optional cluster
search retains every qualifying core group under the declared similarity
threshold, including groups beyond the nearest reference ties. The component
builder and modern coverage benchmark use explicit relative-similarity
thresholds. These outputs have different selection rules and must not be
interpreted as one common nearest-hit list.

Input normalization uppercases sequences and converts U to T. Mixed T/U is
accepted under this policy; the tool does not infer the molecule type. Results
contain normalized sequences and the original FASTA file's SHA-256 hash, not a
copy of its original formatting and case. Preserve the input file for exact
reproduction.

Search in both orientations addresses how DNA sequence strings were recorded.
A reverse-complement hit identifies an orientation relative to a reference
string. It does not determine which genomic strand produces mature crRNA.
Neither the search score nor a tie-breaking rule supplies strand evidence.
Palindromic repeats and equally supported orientations can remain unresolved.

The configured distance or similarity cutoff is a fixed operational filter. It
is not a posterior probability, a calibrated error rate or a boundary between functional
and nonfunctional repeats. A result outside the filter means that the configured
search found no qualifying reference, not that the repeat is biologically novel.
Likewise, resemblance does not prove common function or establish host taxonomy.

Ambiguous bases remain explicit. Neither workflow chooses an arbitrary base,
treats an ambiguous character as an exact wildcard match, or silently shortens
the query. The modern catalogue preserves ambiguous and out-of-policy source
sequences as observations while excluding them from ordinary ACGT search.
Report unsupported sequence evidence separately from a completed search with
no qualifying match.

## Why historical labels are separate

The original 2013 study used independent sequence-family and structural-motif
analyses. Its examples show that sequence conservation and conserved RNA
structure need not coincide. A nearest sequence match therefore cannot
substitute for structure-aware motif classification. The historical mapping
command does not infer motif identity for a nonexact query. The modern catalogue
extension can report separate thermodynamic folding hypotheses; these are not validated motif assignments.

Source: [Lange et al., CRISPRmap (2013)](https://doi.org/10.1093/nar/gkt606).

The historical 2014 update used CRISPRstrand orientation predictions and changed
the reference partition. It reported revised orientations for 536 of the 3,527
previous repeats. Family numbers from the two releases consequently cannot be
treated as interchangeable identifiers, and reference direction alone is not
experimental strand evidence.

Source: [Alkhnbashi et al., CRISPRstrand (2014)](https://doi.org/10.1093/bioinformatics/btu459).

Cas-system classification has since expanded. The 2025 evolutionary update
describes two classes, seven types and 46 subtypes. This taxonomy is based on
system evidence and does not convert historical CRISPRmap repeat families into
current Cas-subtype calls.

Source: [Makarova et al., updated classification (2025)](https://doi.org/10.1038/s41564-025-02180-8).

Repeat-only subtyping also has intrinsic ambiguity: CRISPRclassify discusses
closely similar repeats shared between some type I and type III systems through
shared adaptation machinery. A repeat reference tool should retain such
conflicts rather than force one subtype. RepeatTyper and CRISPRclassify are
possible future comparators for a separately scoped subtype method.

Sources: [Nethery et al., CRISPRclassify (2021)](https://doi.org/10.1089/crispr.2021.0021),
[Russel et al., CRISPRCasTyper (2020)](https://doi.org/10.1089/crispr.2020.0059).

## Reading the historical benchmark

See the [executed historical benchmark report](../benchmarks/results/REPORT.md)
for its workloads, counts and limitations. It evaluates the historical mapper
and recovered original evidence; it is not a report of the modern catalogue
pipeline. The rules below interpret those results without adding unmeasured
performance claims.

Reference self-recovery verifies indexing, orientation handling and annotation
preservation. It is not a test on unseen sequences. Concordance with authentic
original output measures agreement with that output and its release; it is not
independent ground truth. If the original executable was not run, its runtime
and a speedup comparison remain unavailable.

A held-out reference analysis must remove every occurrence of the same sequence
and its reverse complement from the search reference. Report its unit of
analysis, split rule, denominator, missing labels and conflicts. Such a split
tests recovery of historical associations for different sequence strings; close
homologs can still occur on both sides. Do not call it homology-held-out or
independent biological validation unless those stronger separations were
actually enforced.

Report accepted coverage alongside agreement among accepted cases. Coverage
alone can reward overassignment; agreement alone can conceal widespread
abstention. Per-family results are needed because abundant families can dominate
pooled counts. Distinguish lack of a historical annotation from disagreement
with an available one. Fixed candidate thresholds must be described as
uncalibrated policies, including when they happen to perform well.

Synthetic mutations, shuffled sequences, and length/composition-matched controls
are useful for specific correctness and sensitivity checks. Their inherited or
constructed labels do not turn them into experimentally established biological
examples. An unavailable comparison should be reported as unavailable, not
estimated from a superficially similar task.

## Reading the modern benchmark

See the [executed modern release report](../benchmarks/modern-20261001/RELEASE_REPORT.md)
for verified full-catalogue counts, grouping, folding and benchmark results.

The implemented modern evaluation samples unique eligible Bona-fide consensus
sequences from a cohort held out by accession stem. Its search reference contains
Bona-fide training consensuses, rather than every searchable Possible or
observed-only node in the catalogue. Assembly eligibility requires complete
outcomes and valid upstream output validation; that validation does not
experimentally confirm an array or its biological function.

Coverage means finding a training reference at similarity at least 90%, using
full-length unit edit distance normalized by the longer sequence and checking
both orientations. The accession split keeps accession versions together and
excludes cross-fold identical contig hashes and unknown contig hashes. It does
not guarantee species, strain, or near-homology separation. Taxonomy availability
and overlap are audited but taxonomy is not used to define these folds.

Ordinary accession-held-out coverage may include identical sequences shared
across different accessions. The separate exact/RC-purged evaluation removes
training nodes matching any retained held-out consensus or observed repeat,
including unsampled arrays and both detector tiers. The sample-scoped homology
stress additionally removes training nodes at the same 90% threshold from
sampled queries and their eligible associated observed variants. Excluded
ineligible variants must be reported. Zero coverage is expected by construction
in that stress test; it is not a measure of family generalization or a claim of
full-cohort homology separation.

The optional comparison with the frozen historical reference uses the same
coverage cutoff. Its genome overlap is unknown, and coverage differences do not
establish biological accuracy or predictive superiority. The independent
dynamic-programming oracle checks a bounded reference subset; it is separate
from full index-inventory validation. Original-executable speed remains
unavailable. Detector training overlap and independent modern family labels
are also unavailable.

A complete catalogue, a pilot run, or completion of one processing lane does
not establish completion of the full index, groups, folds or benchmark. Read
those results from their executed release reports and complete manifests.

## Gates for stronger biological claims

1. **Scientific curation beyond catalogue coverage.** The modern importer
   implements accession versions, array locations, occurrence-preserving
   deduplication, source hashes and available taxonomy. Recorded sequence
   orientation is not transcription-strand evidence. Curated labels, independent
   assessment of detector errors and missing metadata need further evidence.
   Data license terms remain separate from the software license.
2. **Predictive family assignment.** Define the family vocabulary and label
   evidence. Partition development and held-out data by sequence homology and
   genome/strain, with sensitivity to grouping criteria. Choose thresholds only
   on development data, then measure accepted error versus coverage, per-label
   behavior and uncertainty on the untouched test data.
3. **Unknown-family behavior.** Remove entire families from the reference and
   evaluate abstention. A family absent from the vocabulary cannot be correctly
   assigned to that unavailable label. This is a different question from
   generalization within represented families.
4. **Structural motifs.** Recover authentic motif alignments and provenance,
   select a structure-aware method, and validate motif mapping independently.
   The implemented thermodynamic folds remain model hypotheses and do not
   establish cellular processing, protein-bound structure or immune activity.
5. **Strand and subtype.** Treat each as a separate prediction task with its own
   experimentally or contextually supported labels, baselines, ambiguous cases,
   and calibration. Audit training/reference overlap before claiming an
   independent external test of a pretrained model.

The implemented catalogue, groups and folding methods address reference
analysis and evidence reporting. They do not establish that these stronger
predictive and biological validation requirements have been met. The
[council decision](COUNCIL.md) records both implementation scopes.

## Modern catalogue extension

See [modern reference policies](MODERN_REFERENCE.md) for quality tiers, canonical
sequence provenance, operational sequence components, independent folding
evidence and leakage-aware evaluation. Modern groups have their own namespace
and cannot replace historical family labels. Accession and source taxonomy
record where an observation arose; sequence similarity does not assign a host.
