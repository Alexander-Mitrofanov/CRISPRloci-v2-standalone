# Contracts and limits

`CRISPRrepeat-report` schema 1.0.0 includes tool/model identity, parameters,
source hashes, input kind and molecule, source-bound arrays, local contexts,
consensus folds and observed repeat instances. Each instance preserves its
source sequence, oriented RNA/hash, source interval (null for repeat-only),
context ID/offset, isolated result, comparison and reference stage.

Pair triples are `[i,j,p]` with zero-based `i < j`. Isolated pairs are global
marginals. Context pairs are pair-specific window averages, which must not be
row-normalized or constrained to sum to one. Unpaired arrays come from the
appropriate engine. All numbers are finite. Genomic intervals are half-open.
Sequence transformations must be reversible; ambiguous substrates are unsupported.

Supported limits: 1000 FASTA records, 200 nt per repeat, 5000 nt per array,
20000 total array bases, 1000 observed repeats, 250000 pairs per context,
20 MB final JSON. Window is 10–300 nt, span is 5–window, temperature 0–100°C.
Limits are operational bounds, not biological cutoffs. Normalized model values
are retained for each execution. A resource error must not produce a completed
report. The outer worker must additionally bound runtime, RAM and total storage.
The direct folding API applies the same repeat/context length bounds before
native allocation. Native structures, energies and local row coverage are
checked before completed evidence is returned. Model defaults are explicitly
pinned even when ViennaRNA process-global defaults have been changed. Each
Folder owns temperature-scaled Turner2004 energy and Boltzmann tables; every
native fold receives copies before computation, retaining its effective local
window/span even if another library later replaces global parameter tables.

Optional reference schema:

```json
{
  "schema": {"name": "CRISPRrepeat-reference", "version": "1.0.0"},
  "release": "operator-owned-immutable-id",
  "license": "reviewed redistribution terms",
  "source_url": "https://reference-owner.example/release",
  "build_date": "YYYY-MM-DD",
  "entries": [{
    "id": "reference-id", "sequence": "ACGT", "sha256": "sequence-sha256",
    "occurrences": [{
      "accession": "accession.version", "source_url": "https://source.example/record",
      "evidence": "origin and annotation evidence"
    }]
  }]
}
```

Maximum 1000 reference entries, each at most 200 nt. Occurrences may add
documented array intervals, taxonomic information and evidence types; they
are reference origins, not inferred query origins. Search uses complete-query
and complete-reference edit alignments, maximum distance 3. A near match is
not a family or subtype assignment. Reference release and content hash are
retained with search results. No database is fetched at job time.
An exact banded distance filter rejects distant entries before full traceback;
shared reference sequences reuse alignments while retaining separate origins.
