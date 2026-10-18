# Sources and retained resources

The software rewrite derives from BackofenLab/CRISPRtracrRNA release **2.0**, commit
`d4cbbbb2172d4114c28fdb8a93c8abae871140f3`. The upstream MIT license is retained as `LICENSE`.
The development workspace preserves the original source separately in `old_version/`.
That historical checkout is not included in this standalone v3 repository and is not
required at runtime; the upstream revision above identifies the baseline.

`src/crisprtracrrna/data/` contains byte-for-byte copies of:

- `tools/cov_matrix/s2_b14/cov_matrix.cm` and four Infernal indexes, renamed `II.cm*`;
- `tools/cov_matrix/type_5/cluster{1,2,3,4}.top30.cm` and indexes, renamed `V-{1,2,3,4}.cm*`;
- `tools/erpin/rho-indep.epn` and the historical Linux x86-64 `erpin` executable.

Only the copied ERPIN file's executable permission is changed. The resource SHA-256 manifest
is checked by the preflight command. The upstream repository supplied these third-party ERPIN
assets without a separate ERPIN license file; this refactor does not relabel their authorship
or independently resolve their redistribution terms. The wheel is tagged Linux x86-64 because
it includes that executable. `erpin` and `erpin_model` can select a separately installed version.

The example FASTAs are copied from the same v2 release. They are execution/regression fixtures,
not an independent biological accuracy benchmark.

CRISPRidentify v2 and standalone CasAndra are downloaded into `external/`, excluded from the
v3 wheel and sdist, and retain their own licenses and provenance. Exact tested upstream commits
are recorded in `upstream-lock.json`. No vendored Python imports, monkey patches, source edits,
or copied CRISPRidentify/CasAndra runtime modules are used by v3.

`native/fasta36-seeded-20260927/` retains an optional, separately named FASTA36
reproducibility build, its exact source archive, patch, build provenance, Apache 2.0
license and copyright. It is not installed by the Python wheel. See
[native/README.md](native/README.md) for its tested one-thread scope and build limitations.
