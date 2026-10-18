# clinker-v2 for CRISPRloci

This derivative retains clinker and adds a bounded, offline machine-readable
entry point for comparing annotated genomic neighborhoods. Its protein
alignments use the original global BLOSUM62 method and identity/similarity
definition. Similarity is not orthology, function or an organism phylogeny.

Install into a dedicated Python 3.12 environment:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements-runtime.lock
.venv/bin/pip install --no-deps .
.venv/bin/clinker-v2 compare --input sources.fasta --features features.json --outdir results
```

`features.json` has schema `{name: "clinker.features", version: "1.0.0"}`,
`coordinate_system: "1-based-inclusive-source-forward"` and `loci`. Each locus
contains `id`, `source_id`, `start`, `end` and `genes`. Each gene contains `id`,
`start`, `end`, `strand` (`+` or `-`), `phase` (0–2), `translation_table`,
`partial_5prime`, `partial_3prime` and `translation`. Genes must be ordered by
`(start,end,id)`. Source identities and gene identities preserve the supplied
evidence; the server supplies CasAndra's complete CDS set rather than only Cas hits.

Each translation is checked against oriented source DNA under the stated code.
Phase skips bases at the oriented 5′ end. Incomplete terminal codons are permitted
only with an explicit partial 3′ flag. Internal stops are rejected. Translation
follows CasAndra's ordinary codon translation without forcing initial residues
to methionine. This contract accepts linear, fully contained CDS intervals;
origin-spanning features must be prepared as explicitly linear source windows.

The fixed bounds are 2–8 loci, 50,000 nt per locus, 120 genes overall, 3,000 aa
per protein, 6,000 cross-locus protein pairs and 250 million alignment cells.
Oversized inputs fail before alignment; evidence is never silently truncated.
No genes and no qualifying links are valid negative results. The CLI uses one
process, does not download data and does not start the upstream HTML server.
`--identity` defaults to 0.3; every cross-locus protein pair is assessed.

Outputs are `comparison.json`, `alignments.tsv`, `proteins.faa` and `loci.gbk`.
The JSON preserves original coordinates and source/translation hashes plus
actual gapped alignments and scores. GenBank feature coordinates are local to
each exported window, with original source coordinates in its comment. Only
sanitized, independently validated evidence is exported by the web server.

Audited fixes covered by regression tests include configuration mutation between
alignments, synteny weight shadowing, generator-valued adjacency pairs,
single-linkage bridge merging, stale group assignments, alignment identity loss
on deserialization, and loss/crashes in circular-feature part merging. The
GenBank fallback now honors translation table/codon start and avoids evaluating
an unused fallback when a translation was explicitly supplied. These are fixes
for identified defects, not a claim that every possible defect has been removed.

Run `.venv/bin/python -m pytest -q tests` after installing pytest. The synthetic
fixtures test reverse-strand source binding, alternative codes, partial phase,
tampering, empty calls, limits, exact alignments and the upstream regressions.
Application-level biological accuracy still depends on gene calling and the
quality/comparability of the chosen neighborhoods.

## Code navigation and invariants

- `clinker/bounded.py`: input/schema bounds, source-oriented CDS verification,
  work accounting, ordered protein comparison, and atomic four-file export.
  The 128-entry alignment cache is local to one comparison; duplicate proteins
  retain separate gene evidence and the full declared work counts.
- `clinker/align.py`: upstream alignment, identity/similarity classes, homology
  groups and synteny. Scores and first-optimum alignment ties remain Biopython's;
  similarity classes exclude extended residues unless they are identical.
- `clinker/classes.py`: native GenBank/GFF parsing and serialization. Native
  GenBank ranges filter source coordinates while retaining the full source DNA
  for CDS fallback translation.
- `tests/test_bounded.py` and `tests/test_upstream_regressions.py`: original
  source-binding, alignment and circular-feature controls.
  `tests/test_refinements.py`: independent residue arithmetic, duplicate evidence,
  range extraction, deserialization and functional-group bridge regressions.

Translation retains Biopython 1.85 ambiguity resolution: unambiguous IUPAC codons
can resolve to one residue, D/N to B, E/Q to Z, I/L to J, and unresolved codons to
X. No initiator-M substitution, biological cutoff, scoring model, output schema
or provenance hash policy is changed by caching or accounting optimizations.

## License and redistribution

The upstream code is MIT; `LICENSE` and its attribution are retained. The source
revision is recorded in `UPSTREAM.json`. Exact installed dependency notices and
source URLs are in `licenses/` and `licenses/runtime-manifest.json`. Biopython's
license and BSD/MIT/Apache dependency notices must accompany redistribution.
The upstream PyDNA-attributed helper's BSD-style notice and D3 notice are also
retained. The bounded CLI does not distribute HTML through CRISPRloci.

`argh` is LGPL-3-or-later, a transitive dependency of the unchanged `gffutils`
package. Keep it a separate, replaceable unmodified package and retain its
LGPL/GPL notices and upstream source link. NumPy/SciPy binary wheels contain
additional runtime notices (including GCC runtime exceptions); preserve their
complete notices if redistributing environments. Installing the original
packages separately for server operation does not require relicensing the
web application. Do not describe the entire environment as MIT-only.

Version 2.0.1 preserves validated leading-gap protein alignment strings verbatim
in `alignments.tsv`. Arbitrary identifiers remain spreadsheet-formula escaped;
JSON identities, scoring, alignment strings and thresholds are unchanged.
