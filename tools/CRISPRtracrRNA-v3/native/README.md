# Optional FASTA36 reproducibility profile

The normal pipeline can use an installed FASTA36. The controlled development
comparison instead used the separately named executable in
[`fasta36-seeded-20260927/`](fasta36-seeded-20260927/), preserving the stock installation.

This repository retains the exact upstream source archive, Apache license and
copyright, explicit-seed patch, build scripts/logs/manifests, and tested Linux
x86-64 executable. The expanded upstream source/build tree is excluded from Git;
it can be recreated from the retained archive. The modification notice describes
the original complete development workspace.

The patch fixes missing RNG seed arguments. It does not change alignment scoring,
shuffling counts, models or reporting thresholds. Only the declared Linux SSE2
pthread profile with **one native thread** was tested. MPI and arbitrary-thread
reproducibility are not claimed.

For that profile, adapt your own `config.toml`:

```toml
threads = 1
fasta36 = "native/fasta36-seeded-20260927/fasta36-council-seed20260927"
```

Keep the upstream application commands/environment in the same configuration.
The workstation's `local.toml` and `local.reproducible.toml` are intentionally not
published. The optional binary is not installed by the Python wheel.

To inspect or rebuild the original profile, start from a **fresh extraction** of
the retained source archive and read `build_seeded.py`, `explicit-seed.patch`,
`build_manifest.json` and `MODIFICATION_NOTICE.md`. The recorded build uses a
host-specific zlib path. Adapt compiler/library paths on another host and record
new hashes; do not reapply the script to an already patched extraction or overwrite
the historical artifacts. The expected source archive SHA-256 is
`fa5318b6f8d6a3cfdef0d29de530eb005bfd3ca05835faa6ad63663f8dce7b2e`.

See the [reproducibility amendment](../docs/COUNCIL_REPRODUCIBILITY_AMENDMENT.md)
and [benchmark log](../BENCHMARK_LOG.md) for the failed stock-native comparison and
the separately controlled rerun.
