# Upstream provenance and third-party code

SpacePHARER native algorithms and workflows originate from Ruoshi Zhang,
Milot Mirdita, Eli Levy Karin, Clovis Norroy, Cyril Galiez and Johannes Söding.
Upstream: https://github.com/soedinglab/spacepharer
Exact imported commit: d1e94c04025be6e04fd8204f4f8bb1ea782e993c (8 May 2024).
Paper: https://doi.org/10.1093/bioinformatics/btab222 (2021).

The full GPLv3 license and upstream notices are retained. The first-party code in
src/ and offline Python package are modernized here; lib/mmseqs is vendored
third-party code from that exact upstream tree and is not claimed as refactored.
Its own licenses and bundled third-party notices apply. No GitHub release date is
inferred from the date of the upstream commit. Frozen example biological inputs
and their accessions are those distributed by upstream. Example databases are
validation fixtures, not a current comprehensive viral database.

Historical Azure deployment pipelines were excluded from both source trees; they
are unrelated to compiling or running the algorithm and contained upstream server
connection commands. Scientific vendored MMseqs source remains unchanged.
