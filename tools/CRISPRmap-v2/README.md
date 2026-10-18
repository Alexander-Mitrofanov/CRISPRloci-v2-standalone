# CRISPRmap-v2

**Offline prokaryotic CRISPR repeat reference mapping, 2026 implementation.**

Map repeat FASTA to versioned references, recover exact historical annotations,
and inspect nearby sequence evidence with explicit uncertainty. Analysis runs
offline after references are prepared. The historical mapper uses the Python
standard library; large-catalogue search and RNA folding have separately pinned
optional dependencies. Historical results include JSON, TSV and portable HTML;
modern catalogue search streams JSONL and preserves complete ties in SQLite.

This independent implementation is distinct from Freiburg's **CRISPRmap v2
(2014)** and the unrelated optical-screening tool called CRISPRmap. The
historical baseline release, **2.0.0a1**, implements historical reference mapping.
It does not recreate
the original LocARNA structural clustering, superclass inference, or
CRISPRstrand model. Near matches are candidates, not calibrated classifications.

The **2.1.0a1** extension adds a provenance-preserving modern catalogue builder,
disk-backed indexed search, versioned operational sequence groups and separate
RNA folding evidence. The completed ClusterCRISPRCasDB reference contains **2,184,322 canonical
repeat sequences** from **9,583,665 arrays**, with **206,946 eligible Bona-fide
consensus sequences** forming the grouping core. The full import, search index,
three group thresholds, 413,892 folding predictions and final CLI checks passed
validation on the campaign VMs. See the [full release report](benchmarks/modern-20261001/RELEASE_REPORT.md).
See [modern build instructions](docs/BUILD_MODERN.md) and
[scientific policies](docs/MODERN_REFERENCE.md).

The **2.1.0a2** bugfix release rejects duplicate or incomplete historical metadata
columns and rows, and rejects duplicate keys, nonfinite numbers and oversized
numeric tokens in modern source JSON. Explicitly empty annotations remain
missing. Catalogue merging now opens immutable partition URLs correctly across
SQLite builds. Historical `classify` schemas, retrieval policies and defaults are unchanged.

## Install and try

Historical mapping requires Python 3.10 or later. Modern search and folding use
the optional dependencies shown below.

```bash
git clone https://github.com/Alexander-Mitrofanov/CRISPRmap-v2.git
cd CRISPRmap-v2
git config core.hooksPath .githooks
python3 -m venv .venv
.venv/bin/pip install .
.venv/bin/crisprmap classify \
  --input examples/repeats.fasta \
  --reference examples/reference.json \
  --output runs/demo
```

Open `runs/demo/report.html`. The bundled example is entirely synthetic; its
`DEMO` labels are engineering fixtures without biological meaning. Output paths
must be new. A completed `report.json` is published only after successful
analysis and companion exports. Existing artifacts are never replaced.

## Modern catalogue search

Install the pinned scalable-search and folding dependencies, then use a completed
catalogue, its index and optional sequence groups:

```bash
.venv/bin/pip install '.[scale,structure]'
.venv/bin/crisprmap catalog-search \
  --input your-repeats.fasta --catalog data/catalog --index data/index \
  --clusters data/clusters --output runs/modern-search
```

`results.jsonl` contains displayed nearest hits and operational group evidence;
`results.sqlite3` retains every equal-best match. Possible-only and observed-only
sequences remain searchable but cannot connect the Bona-fide consensus groups.
Modern groups have separate `MAP2` identifiers. Folding results are independent
hypotheses for both orientations, without historical motif assignments.

See [build instructions](docs/BUILD_MODERN.md) for the source-plan format,
partitioned import, index construction, grouping, folding and benchmark commands.

## Historical reference

The optional importer retrieves **3,527 original 2013 records** and the original
16-query example input/output. URLs, byte lengths and SHA-256 digests are pinned.
This is an explicit operator step; analysis never downloads or updates data.

```bash
.venv/bin/python scripts/fetch_legacy.py
.venv/bin/crisprmap inspect-reference --reference data/legacy/2013/reference.json
.venv/bin/crisprmap classify \
  --input your-repeats.fasta \
  --reference data/legacy/2013/reference.json \
  --output runs/your-analysis
```

Raw historical assets remain in ignored `data/` because their separate
redistribution terms were not established. All reports name the label namespace;
2013 and 2014 families are different partitions and must not be interchanged.
See [source provenance](docs/LEGACY_PROVENANCE.md).

## Reading historical evidence

| Result | Meaning |
| --- | --- |
| `exact` | Sequence or reverse complement occurs in the reference. A family or motif value is reported only when all tied occurrences have the same nonmissing annotation. |
| `near` | Closest full-length references within the requested edit limit. Family labels are candidates; motif labels describe those references only. Assigned values stay null. |
| `no_hit` | No eligible reference within the limit; this does not establish biological novelty. |
| `unsupported` | Ambiguous IUPAC query bases prevent unambiguous matching; no wildcard or fabricated match is used. |

Every best tie contributes to annotation evidence before `--top-k` limits the
display. Conflicting and incomplete labels remain explicit, including in exact
matches. Unannotated references are not negative examples. Reference-relative
orientation is not transcription strand, and no Cas subtype or host is inferred.

The engine uses exact unit-cost global Levenshtein distance, searches both
orientations, indexes duplicate/reverse-complement groups, and uses safe length
lower bounds. One deterministic alignment is reported for each displayed hit.
Similarity is `1 - edits / max(query length, reference length)`, not a probability
or alignment identity. The default `--max-distance 3` is a fixed, uncalibrated
retrieval filter. Set `--max-distance 0` for exact recovery alone.

FASTA may be multiline; case is normalized and U becomes T. Full identifiers
are preserved and must be unique. Input hashes bind reports to original bytes.
Limits include 200 nt per repeat, 10,000 query occurrences, bounded input/output,
a default five-million-comparison reservation and a 300-second analysis deadline.
Exceeding a limit fails explicitly. See [contract](docs/CONTRACT.md).

## Bring another frozen reference

Supply FASTA plus TSV columns `id`, `family`, `motif`, optionally `accession`.
Empty fields mean missing. The builder binds metadata and sequence checksums;
inventories never merge themselves automatically.

```bash
.venv/bin/crisprmap build-reference \
  --fasta reference.fasta --metadata annotations.tsv \
  --namespace your-reference-release-2026-01 \
  --source-url https://example.org/your-reference-release \
  --output your-reference.json
```

This supports current curated references without treating old family names as
modern Cas subtypes. Input reference quality remains the curator's responsibility.

## Benchmark and verification

The [recorded benchmark](benchmarks/results/REPORT.md) compares against authentic
published outputs from original CRISPRmap, then separately measures reference
recovery, duplicate/RC-group holdout, family-absent challenges and synthetic
controls. The original service is offline and no executable was recovered:
**original runtime and speedup are unavailable**. Historical agreement does not
establish improved biological accuracy on modern genomes.

The [modern benchmark](benchmarks/modern-20261001/RELEASE_REPORT.md) evaluates
1,000 unique heldout core consensus sequences at 90% global edit similarity,
considering both orientations:

| Reference or evaluation | Coverage |
| --- | ---: |
| Modern training core, accession-stem holdout | 76.1% |
| Modern training core after full exact/RC overlap purge | 66.1% |
| Frozen historical 2013 reference | 18.5% |

These are reference-coverage measurements using the new mapper. The accession
split does not establish species or homology independence, historical overlap
is unknown, and the figures do not measure biological family accuracy.
All **168 tests** and lint passed on a campaign VM; CI also passes on Python
3.10, 3.12 and 3.13 and with pinned scientific extras.

```bash
.venv/bin/pip install '.[test]'
PYTHON=.venv/bin/python bash scripts/check.sh
.venv/bin/python benchmarks/run_benchmark.py
```

Optional browser smoke test for the HTML report:

```bash
.venv/bin/pip install playwright
.venv/bin/python -m playwright install chromium
.venv/bin/python scripts/check_report_browser.py runs/demo/report.html
```

The [code map](docs/CODE_MAP.md) links workflow boundaries to their regression tests.
The [agent council consensus](docs/COUNCIL.md) records scope and rationale.
[Scientific interpretation and future validation gates](docs/SCIENCE.md)
describe the work needed before claiming modern biological classification:
current independent labels, homology/taxonomic holdouts, calibration, and separate
structure-based motif validation. The council is an engineering/scientific
review aid, not independent experimental validation.

## Citation and license

Cite the original [Lange et al. (2013)](https://doi.org/10.1093/nar/gkt606) when
using its historical reference, and
[Alkhnbashi et al. (2014)](https://doi.org/10.1093/bioinformatics/btu459) when
discussing historical orientation/class updates. This new code uses the MIT
license; third-party data retain their own terms. See `CITATION.cff`.
