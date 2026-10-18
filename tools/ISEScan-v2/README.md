# ISEScan-v2

A maintained ISEScan 1.7.3 derivative that produces bounded, source-bound JSON
for inspecting insertion-sequence context around CRISPR loci. It retains the
upstream detection thresholds and supplied transposase models. It does not infer
loss of function from an interval overlap.

## Install and run

Create a dedicated runtime from `environment.yml`, then install this checkout:

```sh
conda env create --prefix /your/runtime/isescan -f environment.yml
/your/runtime/isescan/bin/python -m pip install .
/your/runtime/isescan/bin/isescan-v2 --input genome.fasta --output results --threads 1
```

A C compiler is required to build the bundled MIT-licensed SSW library. Installation
is operator-owned; execution is offline and never downloads programs or models.
Keep `FragGeneScan`, `hmmsearch`, `phmmer`, `blastn`, `blastp` and `makeblastdb` on
the runtime PATH. Their files and reference models must remain immutable during
jobs. Pin the resolved environment with a package manifest before deployment.
The package includes Python sources and reference models in its wheel; it does
not depend on the original checkout after installation. Supported build/pilot
platform: Linux x86-64. Other architectures require a separate native test.

The output directory must be new. The CLI accepts at most 200 DNA records and
2,000,000 total bases; native identifiers are generated privately and mapped back
to source IDs. Duplicate IDs and malformed DNA are rejected. The standalone
process has a bounded timeout (default 3,600 seconds) and a separate process
group; when run by the web worker it uses the already isolated group.

## Results

- `isescan.json`: typed predictions, source/sequence hashes, model hashes and
  explicit uncertainty. Coordinates are 1-based inclusive in the source-forward
  coordinate system.
- `native.tsv`: unchanged native table with private generated source identifiers.
- `work/`: private native logs, source, translations and search intermediates.

`native_type` retains upstream `c` (complete) or `p` (partial); contig-edge status
is separate and does not itself mean truncation. `transposase_orf` is null when
not identified; its strand is the transposase strand, not an independently
inferred direction of the whole insertion element. `evalue_this_copy` is separate
from the best E-value among related copies. Missing negative/zero sentinels are
normalized to null only in their documented contexts. Terminal inverted repeat
alignments are reconstructed from source DNA and validated. Target-site
 duplication is always null with status `not_assessed`.

A completed zero-call result is possible. Missing outputs or failed executables
are errors, not evidence of absence. Finding no predicted insertion sequences
with these models does not establish their biological absence. No CRISPR
inactivation or functional conclusion is computed. DNA-level interrupted-homolog
analysis remains outside this method.

## Audited fixes

The native algorithm is maintained in `isescan_v2/native/`. This release fixes:

- shell interpretation and broken argument boundaries in FragGeneScan/HMMER/BLAST
  calls, including filenames containing spaces;
- FragGeneScan's Perl launcher reintroducing shell interpretation and failing to
  propagate the executable's exit code (the same native options are invoked
  directly instead);
- ignored failed streaming BLAST invocations;
- reference/library discovery depending on the process entrypoint or working
  directory.
- truncated cached HMM tables causing a negative seek instead of a new search;
- BLAST query-length inconsistencies being ignored, and missing self-alignment
  or source ORFs silently assigning unrelated copy evidence;
- temporary BLAST subject files surviving failed executions and double-gap
  terminal-repeat alignment columns being accepted by the report validator.

HMM search preparation shares one completion check. ORF lookup and ranked
boundary deduplication use indices; overlap discovery compares only intervals
that intersect. The original average-linkage distances, cutoffs, representative
ranking, exact external-program arguments and model files remain unchanged.

A separate report layer validates bounded input, native evidence, source hashes,
coordinates and terminal-repeat alignments. These are audited defects and tested
boundaries; this does not claim every possible upstream defect is eliminated.

## FragGeneScan compatibility pilot

The upstream README says 1.31 output is incompatible and recommends <=1.30;
Bioconda's 1.7.3 recipe requires >=1.32. On 2026-10-12 a real 1.32 pilot completed
successfully with HMMER 3.4 and BLAST 2.16.0. Using upstream `NC_012624.fna`
(42,245 bases), it reported one IS5 partial prediction at 30355–31343, with a
transposase ORF at 30399–31298 and independently verified terminal-repeat
alignment. Unmodified upstream 1.7.3 with the same dependencies produced a
byte-identical native TSV. A 1,000-base repetitive negative input completed with
zero calls. This validates this exact version combination and fixture, not all
possible inputs or equivalence to FragGeneScan 1.30.

The pilot used Python 3.12, NumPy 1.26.4, SciPy 1.17.1, Biopython 1.86 and
fastcluster 1.2.6. The deployment environment file uses conservative exact pins
and must receive its own release smoke test. Source archive hashes:

- FragGeneScan 1.32: `54db8dab643f791aae7b05c310fe0d88ecb07db306211185c4c8524734f334c4`
- HMMER 3.4: `ca70d94fd0cf271bd7063423aabb116d42de533117343a9b27a65c17ff06fbf3`

Run focused regression tests with `python -m pytest -q tests` after installation.

## Code navigation

- `isescan_v2/__main__.py` owns bounded input, native execution and atomic JSON
  publication; `report.py` independently binds TSV intervals and TIRs to DNA.
- `native/isPredict.py` prepares FragGeneScan and HMM searches; `is_analysis.py`
  invokes those programs and finds/restores TIR alignments.
- `native/pred.py` owns ranked HMM hits, copy discovery, interval clustering,
  complete/partial classification and exports. `_intersected_ids` retains the
  original observation indices; `_orf_index` retains the first ranked ORF.
- `native/tools.py` parses FASTA/BLAST and implements upstream coordinate and
  boundary measures. `constants.py` and `pHMMs/` define the scientific method.
- `tests/test_native_selection.py` checks independent interval selection,
  source binding, strand/absence exports and HMM reuse. Native equivalence must
  also be tested with the exact pinned release environment before deployment.

## Licenses and attribution

See `LICENSE`, `NOTICE` and `UPSTREAM.json`. ISEScan uses Apache-2.0; the SSW and
kseq components retain their MIT notices. The upstream-distributed pHMM files
have no separate restrictive notice and remain associated with the exact
Apache-licensed upstream distribution. A separate ISfinder database is not
included. Cite Xie and Tang (2017), DOI: [10.1093/bioinformatics/btx433](https://doi.org/10.1093/bioinformatics/btx433).

The exact FragGeneScan 1.32 source README licenses redistribution under the GNU
GPL without specifying a version. This conflicts with Bioconda's BSD metadata;
use the source notice, retain its source availability and notices when
redistributing a runtime, and do not describe it as BSD. FragGeneScan is invoked
as a separate program. HMMER 3.4 uses BSD-3-Clause with its own third-party notices;
NCBI BLAST+ retains its public-domain and third-party notices. These separate
programs are not redistributed in this source repository.
