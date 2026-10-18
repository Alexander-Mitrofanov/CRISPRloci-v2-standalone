# Scientific interpretation contract

This contract defines what the public interface may claim. It is intentionally
narrower than what a user might infer from words such as “confidence” or
“orientation.”

## Detection

CRISPRidentify v2 detects and categorizes candidate arrays in bacterial or
archaeal genomic DNA. The interface treats category as the primary result.
The emitted model/certainty score is displayed as a raw model score, never
converted to a percentage and never described as a calibrated probability.

The normal policy admits `Bona-fide` and `Possible` calls to evolutionary
preflight. Expert strict mode admits only `Bona-fide`. A stricter policy can
remove the only comparable partner, so a valid detection job can yield no
eligible evolutionary group.

The supplied v2 schema is versioned. The summarizer tolerates additive fields
and known aliases, while the Integration adapter fails closed when required
sequence, coordinate, category, or orientation semantics are absent. Updating
CRISPRidentify requires rerunning canonical, future-alias, and production-worker
contract tests.

## Comparable groups

SpacerPlacer and CRISPR-evOr compare homologous arrays, not arbitrary arrays in
one upload. Integration groups arrays conservatively by reverse-complement
canonicalized direct-repeat identity and transitive spacer overlap. Spacer edit
distance is an explicit bounded matching policy (0–2 in the public API).

A group needs at least two usable arrays. Three or more related isolates are
recommended because a two-array reconstruction contains little phylogenetic
information. Source-frame arrays whose strand/orientation was not determined
are excluded: CRISPR-evOr reverses an entire group and cannot repair one
independently misoriented member.

`completed_no_eligible_groups` is an informative success. Detection remains
valid, but no reconstruction or orientation claim is made. The response and
bundle report all exclusion reasons.

## Reconstruction

SpacerPlacer estimates ancestral spacer gains and losses under its evolutionary
models. Results depend on array detection, homology/grouping, spacer matching,
tree estimation or a valid independent tree, model fit, and numerical
optimization. They are estimates rather than observed historical events.

Zero reconstructed deletions is valid output, but deletion-rate comparisons and
some orientation evidence may then be poorly identified. Warnings remain visible
and travel with the result bundle.

## Acquisition-order orientation

CRISPR-evOr evaluates two group-wide hypotheses: the supplied spacer order and
the completely reversed spacer order. For each group:

```text
ΔlnL = BDM lnL(input order) − BDM lnL(reverse input order)
ΔlnL >  5  → input order supported
ΔlnL < −5  → reverse input order supported
otherwise  → unresolved (ND)
```

This is an evidence threshold, not a p-value, posterior probability, or error
rate. ND is a valid result and must never be coerced to a direction.

“Input” and “reverse input” refer only to array order. They do not by themselves
identify transcription direction, genomic strand, leader location, PAM
orientation, phage/protospacer targets, CRISPR activity, or genome-editing
guides.

The distance representation depends on spacer order. Therefore, when no truly
independent rooted phylogeny is supplied, forward and reverse trees are
estimated separately and both are retained. A tree estimated by a standalone
input-order SpacerPlacer run is not an independent tree and is never reused for
both hypotheses. CRISPR-evOr’s selected `0_results.csv` is the definitive
reconstruction in orientation mode; no provisional standalone reconstruction
is published.

The bounded public summary retains an `orientation.reconstructions` row for
each available `input` and `reverse` hypothesis. Each row contains the group,
Newick tree, canonical numeric `spacer_order`, ordered node states, branch
`gains` and `loss_blocks`, special-event lists, and acquisition/deletion
totals. Canonical spacer IDs have meaning within a group and remain consistent
across its hypotheses; original spacer identity strings and DNA are withheld
from the public summary. Node arrays and every branch event are model estimates,
not observations.

For a decisive comparison, the selected history is the supported hypothesis.
For ND, the input-order history is reported only as the default reference; this
must not be described as input-order support. Both hypothesis reconstructions
remain available when the tool emitted both. Public-summary and browser display
caps may limit interactive detail, while downloadable artifacts retain the
detailed tool output; `reconstructions_truncated` signals a summary group cap.

## Reproducibility

Every terminal job preserves:

- normalized input and its SHA-256 digest;
- source record identifiers and effective bounded options;
- detector categories, preflight selection/exclusion reasons, and warnings;
- tool/report/contract versions and stage durations;
- the reported reconstruction (supported hypothesis when decisive; input-order
  default when unresolved) and, for orientation mode, both hypothesis trees;
- an artifact manifest with sizes and SHA-256 digests.

Raw operator paths, secrets, bearer tokens, and unsafe serialized Python objects
are never public artifacts.

## CRISPRloci v2 annotation contracts

Integration 0.3.0 adds `casandra`, `tracrrna` and `loci` workflow selections while
retaining the existing integration-run schema and evolutionary contracts. The
backend 2.0.0 public modes are `cas`, `tracrrna` and `loci`. Type II and Type V-K
are distinct model families; `tracr_model_type` is `II` or `V`. The experimental
short anti-repeat rescue remains disabled.

Standalone CasAndra consumes each normalized FASTA record with gene mode `auto`.
The native schema-5 run manifest and JSONL checksums must match. Gene coordinates,
strand, nucleotide sequence and nucleotide hash bind every call to the input;
Cas counts and cassette counts must agree with the retained tables. Confidence
and score fields are model evidence, not proof of biological activity.

CRISPRtracrRNA schema 3.0.0 binds source IDs, lengths, sequence hashes and FASTA
file hashes. Candidate sequences must reconstruct exactly from their reported
1-based inclusive intervals and strand. Export IDs are unique, counts cover all
input records, and each retained record matches the report and its manifest.
Complete mode also validates the nested CRISPRidentify report and CasAndra
evidence. Its Cas context must equal the retained native cassette/run records.
Backend validation does not import the Integration or tool validators.

The complete workflow reports evidence-supported candidates with ranking scores.
Model screening reports covariance-model bit scores and E-values. Neither is a
calibrated probability. No-hit output remains a successful analysis and does not
establish absence of a Cas system or tracrRNA. Type V models cover V-K; the UI
does not imply screening of all Type V systems.

Public summaries contain narrow fields only: Cas genes/cassettes, RNA candidate
intervals, evidence types, scores, sequences and source metadata. Browser rows
are bounded and truncation is explicit. `annotations.json`, TSVs and candidate
FASTA are generated independently of the preview cap; complete-locus jobs also
export native array details in `array-evidence.json`. Raw internal work folders
and paths are excluded from publication.

Scientific manifest schema 2.0.0 records five Git repository/commit identities
and the Integration wheel. Operator runtime receipts additionally hash installed
packages, native binaries and bundled models. Backend artifacts retain these
public identities with effective options and normalized-input hashes.
