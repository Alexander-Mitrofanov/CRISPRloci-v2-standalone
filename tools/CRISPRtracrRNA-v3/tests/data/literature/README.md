These three unmodified natural reference loci reproduce failures found in the
literature benchmark: Ain (repeat copies suppress a weaker anti-repeat hit),
Lsp1 (FASTA36 input path truncation), and Kki (FASTA36 corruption when a whole
genome is used as the query against a short repeat library).

The FASTA files come from the Dooley et al. reference collection at the frozen
[TRACR_RNA commit](https://github.com/skDooley/TRACR_RNA/tree/07c5d28d1c1afdc03e81938beaafd1069ce6fe15/data/assemblies/Corteva).
The corresponding tracrRNA sequences are from the `Cas9_database` sheet in
Supplementary Data 2 of [Gasiunas et al. (2020)](https://doi.org/10.1038/s41467-020-19344-1).
The public repository lists no separate license file at this commit; no new
license grant for the source data is implied here.

`fixtures.json` records the exact download URL, file and normalized sequence
SHA-256, reference sequence, and its unique exact genomic match. Coordinates
are 1-based inclusive on the input forward strand, with strand stored
separately. The accepted CRISPRidentify array coordinates and consensus repeats
are frozen from the benchmark, together with a hash of each original upstream
report. They are predicted arrays; they are not independent biological truth.

The tests isolate the anti-repeat search with those fixed arrays, so they do not
rerun CRISPRidentify or CasAndra. The full natural loci are retained because
removing the array or genomic sequence context can remove the original failures.
Real-tool tests are opt-in through `CRISPRTRACRRNA_TEST_CONFIG`.
