# CRISPRspacer

Standalone analysis of CRISPR spacer sequences and **candidate protospacer
matches**. Python 3.10+; no Python runtime dependencies. Searching is offline.
The optional reverse viral search requires local NCBI BLAST+ executables.

Version 0.3.1 keeps the existing result schemas and search policies while
hardening typed evidence replay, corrupt gzip handling, exact-target identity
checks, bounded SQLite metadata reads, and viral catalog deadlines/cancellation.
Previously valid 0.3.0 results remain replayable.

The core engine searches both DNA orientations using exhaustive full-length
exact or bounded substitution matching. A separate BLAST workflow searches a
viral reference with explicitly heuristic retrieval. Both preserve spacer occurrences, validate
source coordinates, and export
reproducible JSON and TSV bundles. It runs independently of CRISPRloci and does
not import its API, worker, Integration package or other scientific tools.

This is an initial computational release. Matches are sequence evidence, not
proof of functional targeting, infection, immunity or host range. PAM inference,
automated host prediction and array detection are not implemented. The BLAST
workflow retains raw local alignments but accepts only full-query ungapped HSPs.
CRISPRidentify reports are accepted as explicitly supplied input.

## Install

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
crisprspacer --version
```

From the private GitHub repository, clone using your authenticated GitHub access
first. Installation itself needs only Python packaging tools. Core `search`,
`self` and `database-search` need no external aligner. `viral-index` requires
`makeblastdb` and `blastdbcmd`; `viral-search` requires `blastn`. BLAST versions
and executable checksums are recorded. Searches never download data; only the
explicit operator commands `viral-plan` and `viral-download` access NCBI.

## Search your spacers against viral references

Version 0.3.0 adds the reverse direction: **spacer FASTA → candidate viral
reference matches** using local `blastn-short`. See the reproducible
[NCBI viral database methods](docs/viral-reference-methods.md) for the full
metadata-plan → checksum-verified download → `makeblastdb` → validation workflow.
Run full downloads on the data VM. No viral genome collection is included here
or downloaded automatically during installation/tests.

```bash
crisprspacer viral-search --spacers my-spacers.fasta \
  --reference /path/to/immutable-viral-reference \
  --expected-reference-sha256 '<SHA-256 of reference.json>' \
  --output results/viral-matches
crisprspacer validate-viral-result results/viral-matches \
  --reference /path/to/immutable-viral-reference
```

The default accepts full-query, ungapped BLAST HSPs with at most two
substitutions; use `--max-mismatches 0` for exact reported HSPs. BLAST can miss
alignments, including terminal-mismatch cases it clips. This is not an exhaustive
Hamming-distance search. Results retain accession versions, source titles,
coordinates, strand, substitutions and BLAST scores. Duplicate spacer occurrences
do not count as independent sequence support. Taxonomy and host assignments are
not inferred. Matching a segment does not establish a complete viral genome.

| Command | Direction / purpose |
| --- | --- |
| `search` | Supplied spacers against supplied target DNA; exhaustive bounded substitutions. |
| `self` | Source-bound spacers against their source, auditing reported-array overlap. |
| `database-search` | Uploaded target DNA against an indexed spacer collection; exact only. |
| `viral-search` | Uploaded spacers against a viral BLAST reference; heuristic retrieval. |

A [council design record](docs/viral-council.md) explains the scope and decisions.

## Search supplied spacers against target sequences

```bash
crisprspacer search \
  --spacers examples/spacers.fasta \
  --targets examples/targets.fasta \
  --output results/search

crisprspacer validate results/search
```

The default is exact matching. `--substitutions 1` or `--substitutions 2` enables
full-length Hamming comparison: substitutions only, no insertions, deletions or
partial alignments. These are computational policies, not biological tolerance
thresholds. IDs identify **occurrences**; different IDs with identical sequence
remain separate in output while sharing search computation.

The same command supports the reverse legacy use case: supply known/collated
spacers as queries and an uploaded target sequence as the target. Results link
back to the original spacer IDs. A supplied spacer collection is not evidence of
complete genome or database coverage, and associations are not host assignments.

## Check targets against the NCBI-hosted spacer collection

Version 0.2.0 adds `database-index` and `database-search`: an immutable SQLite
spacer index and an offline, exact-only target check. See the
[NCBI download and usage instructions](docs/ncbi-spacerome.md). Full data belong
on the data VM; this repository includes only three original NCBI sample records.
CRISPRloci's queued web mode invokes the same scientific command and pins the
reference digest for each job. The viral-search web adapter is maintained
independently in CRISPRloci v2.

## Analyze same-source matches

```bash
crisprspacer self \
  --source examples/source.fasta \
  --identify examples/identify.json \
  --output results/self
```

The included examples are synthetic fixtures, not biological benchmarks.
`self` accepts canonical `CRISPRidentify-report` schema **1.2.0** and checks
source ID, sequence length/hash, coordinate frame, repeat/spacer geometry,
counts and sequence reconstruction. Every supplied report must bind to one
FASTA record. Multiple reports are supplied after `--identify`.

By default, each spacer searches its originating record only. Matches inside
its originating array, another reported array, or an explicitly padded array
region are exported as exclusions. Every reported array participates in the
overlap audit even when its category is excluded from query selection. Omitted
or undetected arrays remain unknown. A source match outside reported arrays is
a **same-record candidate**, not a claim of active biological self-targeting.

`--category` can be repeated. Its default is `Bona-fide` and `Possible`.
Category exclusions, observed ambiguous spacers, deletion boundaries and
unsupported lengths are retained in the query inventory with reasons. Original
spacer ordinals are preserved, including gaps created by deletion events.

An occurrence inventory can also be prepared once and searched later:

```bash
crisprspacer prepare --source examples/source.fasta \
  --identify examples/identify.json --output results/prepared
crisprspacer search --inventory results/prepared \
  --targets examples/targets.fasta --output results/from-inventory
```

## Assembly membership and circular targets

Multi-FASTA membership does not imply that records belong to one genome.
Declare assembly membership explicitly in a metadata JSON file:

```json
{
  "schema_version": "1.0.0",
  "records": {
    "contig1": {"assembly_id": "assembly-A", "topology": "linear"},
    "plasmid1": {"assembly_id": "assembly-A", "topology": "circular"}
  }
}
```

For `self`, pass `--metadata metadata.json --scope same-assembly` to search
across contigs explicitly assigned to the same assembly. Missing query-source
assembly declarations cause an error. `--scope all` searches all supplied source
records while retaining unknown assembly relationships. Unannotated target
records have `array_overlap: "unknown"`.

For supplied targets, pass `--target-metadata metadata.json` to `search`.
For a local reference, pass `--metadata` when building it. Assembly IDs are
user declarations in a shared namespace, not inferred or externally verified.
An external target with a matching record name/hash is not automatically treated
as the spacer's source record.

Without metadata, analysis uses **linear traversal** and no assembly membership;
linear traversal is a search convention, not an inferred biological topology.
Circular traversal requires an explicit declaration. Queries may span the origin
but never exceed one target lap. Ordered target segments describe wraparound;
records are never concatenated into artificial multi-contig targets.

## Immutable local reference packages

```bash
crisprspacer index --targets examples/targets.fasta \
  --name synthetic-v1 --source-label synthetic-fixture --release v1 \
  --output results/reference
crisprspacer validate-reference results/reference
crisprspacer search --spacers examples/spacers.fasta \
  --reference results/reference --output results/reference-search
```

`index` creates a versioned normalized sequence inventory and checksummed
manifest. It is a reusable reference bundle, not a large-corpus search accelerator.
Input bytes and normalized records have separate digests. Every use verifies
the package; metadata cannot be overridden at search time. Reference identity
and the exact target sequences are retained in the result for replay.

The `index` command uses local input only. The separate `viral-plan` and
`viral-download` operator commands acquire NCBI viral RefSeq sources explicitly.
Choose and document data sources, releases and actual coverage before constructing
biological reference collections. No full reference collection ships with the tool.

## Evidence and outputs

All search output directories are new, atomic publications. Existing paths are
never overwritten. For the core `search` and `self` modes, a successful directory contains:

| File | Content |
| --- | --- |
| `analysis.json` | Versioned report, policy, query statuses, coverage, all alignments including exclusions, counts and warnings. |
| `inputs.json` | Exact normalized queries/sources/targets, declared metadata and configuration for replay. |
| `queries.tsv` | Every observed/deletion occurrence, its status, source identity, interval and sequence when available. |
| `matches.tsv` | Candidate alignments after reported-array/padding exclusions. |
| `exclusions.tsv` | Auditable array-overlap and padding exclusions. |
| `coverage.tsv` | Possible, assessed and ambiguous windows per unique query sequence/target. |
| `manifest.json` | Exact file inventory, byte sizes and SHA-256 digests. |

`viral-search` instead writes `result.json` with query statuses, accepted matches,
per-accession candidate summaries, method settings and limits; `inputs.json`;
`matches.tsv`; retained `blast-hsps.tsv` and `blast.stderr.txt`; and `manifest.json`.
Its reference bundle stays separate. `validate-viral-result --reference ...`
rechecks retained HSPs against source sequence intervals and recomputes filtering
and summaries. This checks evidence consistency; it does not make BLAST exhaustive.

Coordinates are 1-based inclusive on the forward target. `relative_strand` is
`+`, `-`, or `both` for a palindromic query; a palindromic match is reported once
per location. For nonpalindromic queries, both orientations can legitimately
match one location under a substitution policy and remain separate evidence.
Target-relative strand does not establish transcription direction. Flanks are
forward-target sequence context; no PAM is inferred. Linear end clipping is
explicit through actual retained length and `null` for absent flanks.

Unknown values are `null` in JSON and blank in TSV. The authoritative JSON
preserves identifiers; potentially executable spreadsheet strings are prefixed
with an apostrophe in TSV. Do not treat that display escape as part of an ID.
Normalized input sequences are included in results; protect the directory
according to the sensitivity of your input data.

## Completeness and limits

The following exhaustive guarantees apply to the core `search` and `self` engine.
The separate `viral-search` retrieval is heuristic, as described above.

Every returned alignment spans the whole query. For `k` substitutions, the
engine uses `k+1` disjoint exact seeds followed by full-window verification.
At least one seed must be unchanged for any qualifying Hamming match; queries
shorter than the seed partition are handled by direct start enumeration.
The exact case uses one full-query seed. This gives exhaustive discovery within
the declared scope when execution completes without a resource error.

Only ACGT query sequences are searched. IUPAC ambiguity is accepted in input,
but ambiguous queries are skipped and target windows containing ambiguity are
unassessed and counted. `N` never acts as a wildcard. Circular targets assess at
most one lap; a query longer than a target has no possible window.

The report distinguishes matches, no candidate matches, no eligible queries
and no assessable windows. A successful no-hit report means no qualifying
candidate in the assessed windows under that policy. It does not establish
biological absence. `validate` checks artifacts, source reconstruction,
alignment metrics and deterministic replay, including omitted hits; it is
computational consistency checking, not independent biological validation or
cryptographic proof of who produced the data.

| Budget | Default |
| --- | --- |
| Submitted FASTA bases | 50 million per FASTA / normalized record inventory |
| Input file or JSON artifact | 256 MiB; FASTA decompressed bytes also bounded |
| Records / query occurrences | 100,000 |
| Query lengths searched | 1–1,000 bases; configurable up to 10,000 |
| Occurrence-expanded alignments, including exclusions | 100,000 (`--max-hits`) |
| Seed candidates | 2 million (`--max-candidates`) |
| Query-sequence/target pair iterations | 100,000 (`--max-pairs`) |
| Retained evidence accounting | 64 MiB (`--max-evidence-bytes`) |
| Conservative work units | 500 million (`--max-work`) |
| Search/publication time | 300 seconds (`--timeout`) |

Budget accounting is conservative and includes zero-window pairs, skipped
scopes and array-overlap work. These bounds do not substitute for operating
system memory/CPU limits when embedding the tool in a service. This release is
for bounded local workloads, not a claim of whole-catalogue scalability.

Exit codes: `0` completed (including a valid no-hit outcome), `2` input or
integrity error, `3` resource limit, `130` keyboard cancellation. Errors are
JSON on stderr for runtime failures; argument syntax errors use standard argparse
diagnostics. Successful summaries are JSON on stdout. A cap or cancellation
does not publish a completed result. A process forcibly killed by the operating
system may leave a hidden staging directory and `.partial` reservation; only
remove those after confirming that its writer has stopped.

## Development and verification

```bash
python -m pip install '.[test]'
bash scripts/check.sh
python scripts/benchmark.py
```

Tests use synthetic data only, including an independent brute-force oracle,
source-coordinate tampering, circular boundaries, missing assembly metadata,
array exclusions, budget failures, manifest tampering and replay. The benchmark
is reproducible synthetic computational evidence, not measured biological
sensitivity. See [verification](docs/verification.md) and
[contracts and scope](docs/contracts.md).

CRISPRloci v2 consumes the standalone commands through its isolated workers and
independent Integration/backend validators. This package does not import the
server, queue or frontend.

## Background

The legacy feature directions are described in
[CRISPRloci, Alkhnbashi et al. (2021)](https://doi.org/10.1093/nar/gkab456).
The importance of distinguishing reported-array overlaps from other matches is
discussed in [Nobrega et al. (2020)](https://doi.org/10.1093/nar/gkaa1071).
These references motivate the feature; their biological thresholds and reported
performance have not been transferred to this implementation.
