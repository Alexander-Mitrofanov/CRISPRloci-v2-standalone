# SpacePHARER 2.1.1

A modernized offline SpacePHARER engine with explicit grouping, versioned evidence,
and immutable target/control references. Native code retains the published
protein/nucleotide search and combined-score algorithms; `python/spacepharer`
provides the validated API and atomic CLI. Both are required. This is a fork of
[SpacePHARER](https://github.com/soedinglab/spacepharer), not a new upstream release.
See [third-party attribution](THIRD_PARTY.md) and [original manual](docs/upstream-manual.md).

## Build and install

Linux x86-64 with SSE4.1, a C++17 compiler, Python 3.11+, POSIX shell, awk and coreutils:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.lock
PATH="$PWD/.venv/bin:$PATH" cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DHAVE_SSE4_1=1 -DHAVE_AVX2=0
PATH="$PWD/.venv/bin:$PATH" cmake --build build --parallel 4
PATH="$PWD/.venv/bin:$PATH" ctest --test-dir build --output-on-failure
cp build/src/spacepharer-native .venv/bin/
.venv/bin/pip install .
PATH="$PWD/.venv/bin:$PATH" spacepharer --version
```

The Python wheel is dependency-free and contains the typed interface; install the
native binary beside its entry point or set `SPACEPHARER_NATIVE` to its absolute
path. The CLI verifies native version 2.1.1. A wheel alone does not contain or
silently download a native executable. No GPU, service or network access is used
in jobs. An explicit operator command prepares immutable references:

```sh
spacepharer prepare-reference --genome phageA.fasta --genome phageB.fasta \
  --name my-target-panel --output reference
spacepharer associate --input spacers.fasta --groups groups.json \
  --reference reference --output results --threads 1 --memory-mb 2048 --fdr .05
```

References support both targeted panels and broad, dated phage collections. They
are bounded to 100,000 single-contig genomes, 5 million bases per genome and 10
billion bases total. These are operational ceilings, not a requested sample size.
No eligible records should be silently truncated to a ceiling: preparation fails
if a bound is exceeded. A collection can be complete under a documented selection
policy, but cannot represent every phage in nature or every unsequenced phage.

For large collections the operator supplies a newline-separated genome path list
(relative paths are resolved against the list file), avoiding command-length limits:

```sh
spacepharer prepare-reference --genome-list genome-paths.txt \
  --name ncbi-phages-snapshot --source-metadata source.json \
  --scope broad_phage_reference --threads 8 --build-timeout 14400 \
  --output reference
spacepharer verify-reference --reference reference
```

`source.json` has exactly four nonempty plain-text fields of at most 1000 characters:
`collection`, `release`, `retrieved_at`, and `selection_policy`. Record the source
release or retrieval date, exact accession selection and exclusions; bind a detailed
operator selection receipt by putting its SHA-256 in `selection_policy`. The tool
performs no acquisition or taxonomic selection itself. Do not call a RefSeq-only
collection all of NCBI, or inferred host annotations independently verified hosts.
ACGTN target genomes retain unknown bases unchanged; other ambiguous symbols must
be handled explicitly by the operator. User spacer ambiguity remains unsupported.

Native target/control indexes are built together; control translated ORFs are
reversed. Streaming tar import preserves each genome's source set without passing
tens of thousands of filenames to a subprocess. Every index has a SHA-256 and size
receipt (16 GiB per member, 64 GiB aggregate). A checksummed SQLite catalogue binds
accession, length, sequence digest and local FASTA path. Every job streams the full
index checksums, validates the catalogue, then reads/hashes only genomes needed to
reconstruct reported alignments (at most eight remain cached). `verify-reference`
also hashes every source genome and is required before operator activation.
Version 2.1.1 accepts references prepared with 2.1.0 or 2.1.1; the database format
and scientific scoring/calibration algorithms are unchanged. Existing verified
2.1.0 references do not need rebuilding for this patch.

`spacepharer.reference/2.0.0` and the public association contract carry compact
counts, scope, source provenance and the manifest digest; the complete private
catalogue is not repeated in every result. References and outputs are atomic new
directories, never overwritten. Keep deployed references operator-owned and
service-unwritable. Native builds have an operator-configurable timeout from 60
seconds to four hours per native step; submitted searches retain their 900-second
native step deadline and 2048 MiB default split-memory limit. Benchmark the full
collection under the deployment's actual cgroup limits before activation.

No toy panel is installed as a production reference. Historical native `downloaddb`
remains operator-only legacy functionality; the typed API never invokes it. BLAST
indexes are not interchangeable with these native target/control databases.

`groups.json` partitions every occurrence exactly once:

```json
{"schema_version":"spacepharer.groups/1.0.0","groups":[{"id":"sampleA","members":["spacer1","spacer2"]}]}
```

Grouping is an explicit scientific assertion by the caller. A contig/source ID is
not automatically a whole-host identity. Records from unrelated hosts must not be
pooled. Up to 1000 occurrences (1–200 bases) and 100 groups are accepted. Ambiguous
DNA and sequences shorter than the engine's 27-base ORF minimum are preserved as
unsupported, never counted as negative evidence. Exact duplicate sequences within
a group are searched once; all occurrence IDs are projected back onto support.
Reverse complements are not silently collapsed. Scores are calculated from unique
sequences; the native CLI itself retains its historical supplied-record behavior.

## Scientific result contract

`report.json` uses `spacepharer.association/2.0.0`; `associations.tsv` is a rectangular
summary. The original hierarchical `native.tsv` and `engine.log` are local evidence.
Each result binds its input/groups hashes, immutable panel identity, groups and
occurrences. Hashes bind the same input byte snapshots that are parsed. Support
includes native best-hit p-values and the complete nucleotide
alignments with 1-based inclusive, potentially descending coordinates. Both aligned
sequences are checked against the original inputs and target genomes.
Repeated association/support records and columns with a gap in both sequences
are rejected; distinct target intervals remain separate support.

Native combined scores are not probabilities. Native FDR estimates are reported
only when the empirical calibration procedure is estimable. A versioned native
calibration receipt supplies status and actual target/control score counts; missing
receipts fail closed. Estimated calibration requires at least two target and two
control scores. Uncalibrated rows are `uncalibrated_candidate` with null
`native_fdr`, never accepted associations. `association` requires the native calibrated filter and a reported FDR no greater
than the requested cutoff. The historical threshold can retain a boundary row with
a larger separately reported FDR; these remain `above_cutoff_candidate`. it is an association hypothesis, not proof of active
infection, cleavage, or an experimentally validated PAM. No-match output is scoped
to the supplied panel and is not proof of no biological interaction.

## Modernization and compatibility

All first-party native utilities now use shared checked numeric/column parsing and
scoped database-reader lifetime. Per-query target indexes preserve native HSP row order
and final nucleotide-match selection while avoiding repeated whole-result scans.
Native searches and version probes share process-group cleanup on success, failure,
deadlines and interruption. Combined-score, geometric e-value and empirical-tail
logic are extracted into testable C++ modules; workflows share typed execution
helpers and fail if `exec` returns. CMake requires modern tools, uses C++17 for
first-party code, and keeps native assertions enabled in Release-mode tests. Legacy
spacer imports use bounded null-terminated input. Native flank and missing taxonomy
handling avoid invalid memory accesses. Vendored MMseqs scoring/search algorithms are unchanged. In 2.1.0 the native
`createdb` source-index container changes from 16 to 32 bits, preventing silent
genome-set collisions above 65,535 source records. A native regression imports
65,537 distinct FASTAs and verifies every source mapping.

Two deliberate correctness fixes are versioned: small control populations formerly
read beyond a slope vector and could report false zero FDR; these now yield an
unavailable calibration receipt and native sentinel -1. Missing taxonomy support
now stays unclassified. Native flank reads beyond genome boundaries yield unavailable
bases. The typed API disables historical broad PAM-pattern guesses. The published
order-sensitive combined-score fallback is intentionally preserved and tested.

The original engine at commit d1e94c04025be6e04fd8204f4f8bb1ea782e993c was built
separately with only `-include cstdint` for current-compiler compatibility. Exact
native calibrated FDR output parity is tested on frozen seeded distributions of
300 target and 300 control scores. Upstream biological examples preserve two
candidate associations, all six aligned supports, reverse-query coordinates,
p-values and combined scores. Those small examples do not establish empirical FDR
and are correctly labeled uncalibrated. These tests establish computational
preservation, not an independent biological accuracy benchmark.

```sh
.venv/bin/pip install -r requirements-test.lock
SPACEPHARER_NATIVE="$PWD/build/src/spacepharer-native" \
SPACEPHARER_ORACLE=/path/to/original/spacepharer pytest
ruff check python tests
ruff format --check python tests
mypy python
python -m build
python scripts/check_private_data.py workspace
```

Tests cover references larger than 100 genomes, lazy sequence binding, catalogue
corruption, bounded provenance, input/group/reference tampering, missing calibration, negative and
unsupported evidence, duplicate occurrences and native descendant cleanup. Native
and real biological tests require the native executable; original parity additionally
requires the separately built upstream oracle. CI builds both from their exact
sources. Exit codes: 2 invalid input/reference, 3 engine/deadline/interruption, 4 I/O.

Please cite Zhang et al., Bioinformatics (2021),
[doi:10.1093/bioinformatics/btab222](https://doi.org/10.1093/bioinformatics/btab222).
