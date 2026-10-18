# Terminal-boundary repair source fixtures

`repeat_boundary_repair.json` contains small excerpts of the public NCBI
accessions NC_017958.1, NZ_CP014859.1 and NC_014962.1 already present in the
project's local benchmark input collection. No detector output was used to
invent source sequence: every excerpt was sliced from the corresponding
`benchmarking/data/genomes/*.fna` input and compared with the source-validated
saved occurrence representation.

Each record identifies the accession, original benchmark call, input sequence
SHA-256, excerpt SHA-256, excerpt offset in the original sequence, and ordered
repeat start/end coordinates relative to that excerpt. Coordinates are
zero-based half-open. The expected consensus length is diagnostic metadata,
never an input to the repair algorithm.

These fixtures are deliberately chosen development failures and are not an
independent biological test set. Historical model-training overlap is unknown.
