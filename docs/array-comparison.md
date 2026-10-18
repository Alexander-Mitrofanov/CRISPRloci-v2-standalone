# Observed array comparisons with CCTK

`array_compare` connects detection `arrays.json` directly to ordered spacer
comparison and sharing networks. The standalone [CCTK 2.0.1 repository](https://github.com/Alexander-Mitrofanov/CCTK/tree/6fb81fec60f4951bf23259199ea3d9bc8f8a3d22)
was modernized, tested and published before integration. It preserves the
scientific algorithms of upstream 1.0.3 while adding bounded inputs, safe optional
executables and reproducible oracle tests. The original publication is
[Collins et al., The CRISPR Journal (2023)](https://doi.org/10.1089/crispr.2022.0080).

## Input and interpretation

Submit the JSON text in the existing `sequence` field:

```json
{
  "mode": "array_compare",
  "sequence": "<the contents of detection arrays.json>",
  "tool_options": {"cctk": {"--identity": "exact", "--min-shared": 1}}
}
```

The [input schema](schemas/arrays-v1.json) is `crisprloci.arrays` 1.0.0, with no `complete` flag. A typed
contract validates unique sources, source descriptors, scoped `(source_id,id)`
array identities, intervals, occurrence identifiers, contiguous ordinals,
sequence lengths, source bounds and nonoverlapping observed units. Coordinates
remain 1-based inclusive source-forward. Repeated DNA keeps each occurrence and
position. Array identifiers and occurrence IDs are never DNA-equivalence symbols.

Unknown strand and `boundary.exact=false` or unassessed boundary dimensions are
accepted. Comparison does not use the separate evolutionary eligibility gate.
Observed spacers need source intervals and unambiguous A/C/G/T DNA. Missing,
not-reported, unverified or ambiguous spacer evidence is rejected explicitly.
An explicit deletion has null sequence and a source boundary; it is retained as
×, contributes no identity, and is never replaced with dummy DNA. Deletion-only
arrays can accompany observed arrays and remain isolates with
`evidence_status=no_observed_spacers`. The submission requires at least one actual
spacer overall. Incomplete repeat evidence is retained as an annotation when its
fields are consistent; repeats do not participate in sharing.

Uploaded source hashes and coordinates are annotations. Without the original
genomes, this job cannot independently verify those genomic origins. Immediate
spacer DNA and the structured array map are independently pinned to the job.

Exact DNA equality is default. `--identity=reverse-complement` explicitly merges
a sequence with its reverse complement for matching only; it does not infer
strand or acquisition order. Array and spacer order always follows the supplied
JSON. `nshared` is the size of the set intersection; Jaccard is intersection size
divided by union size. Duplicate occurrences never increase either set size.
`--min-shared` filters both edges and connected components. All arrays, including
singletons, remain represented. A component need not share a spacer pairwise.

## Limits and interface

Structural limits are 200 arrays, 500 slots per array, 10,000 slots total and
1–1000 bases per observed spacer. The web service additionally enforces its
advertised request/base/record limits: **every observed spacer occurrence is one
record**. The default service record limit is 20; the environment setting accepts
up to 100. The standalone 10,000-slot ceiling is not the web admission capacity.
A complete protected JSON export must fit 20 MB; oversized outputs fail explicitly.

The browser preview shows up to 50 arrays, 50 slots per array and 100 edges, with
explicit truncation. Colors denote shared DNA across shown spacer slots in different arrays; repeated
occurrences and deletion marks remain visible. Original IDs, unknown strand,
counts and a pairwise sharing table accompany a graph. The graph is a sharing
view, not a phylogeny or evidence of horizontal transfer.

The worker runs native CRISPRdiff plotting with bounded dimensions/connector
work; its PNG stays private. The browser renders independently validated data,
never raw tool SVG. Native directories, command logs and private input maps are
excluded from protected archives.

## Outputs and execution

Complete capability-protected `array-comparison.json`, `array-sharing.tsv`,
unchanged reusable `arrays.json`, compact `result.json` and `results.zip` are
published after independent backend checks. JSON preserves original identifiers;
TSV escapes spreadsheet formula prefixes. The backend recomputes every DNA
symbol, overlap, Jaccard, component, isolate, occurrence order and input digest.
A malformed or tampered result cannot become a successful job.

The API only queues. Isolated CCTK workers run offline, with cancellation,
process containment, private job capabilities, atomic queue operations and
retention unchanged. Runtime source and environment receipts bind the Git pin
and separate CCTK environment. The immutable release candidate gate includes a
real detection singleton and a multi-source duplicate/gap fixture under both
identity policies. Provisioning is operator-owned; anonymous jobs never download
or update tools. SpacerPlacer and evOr retain their independent environments.

Tests cover the installed native engine, actual HTTP submission through worker
and protected exports, hostile report/input changes, scoped ID collisions,
unknown/reverse orientation, malformed/partial arrays, duplicate DNA, deletion
slots, deterministic output and browser desktop/mobile behavior.

The consolidated 12-tool immutable runtime was not rebuilt or activated on the
current workstation because available disk space is insufficient. The isolated
CCTK environment and installed server wheels are tested separately; provisioning
and candidate-gate recipes are ready for an operator with adequate storage. No
production service or deployment was changed.
