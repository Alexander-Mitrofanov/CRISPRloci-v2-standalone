# FASTA36 council reproducibility build

This directory retains the exact upstream `fasta-v36.3.8g` source archive,
its license/copyright, an RNG-only patch, modified source, build logs and a
separately named executable. It does not replace the installed FASTA36.

The council amendment fixes missing RNG seed arguments with explicit manager
seed 20260927 and worker seed 20260927 + zero-based worker index + 1. The RNG,
alignment/scoring algorithm, shuffling count and statistical/reporting cutoffs
remain unchanged. Only the Linux SSE2 pthread profile at **one native thread**
is evaluated. MPI builds fail explicitly; arbitrary multithread reproducibility
is not claimed. The separate `local.reproducible.toml` selects this binary and
one thread. Stock-FASTA configurations retain the documented variability.

Build with the preserved upstream Linux SSE2 makefile. The initial compilation
succeeded but linking lacked the unversioned zlib development symlink. Attempt2
links the installed `/usr/lib/x86_64-linux-gnu/libz.so.1` explicitly via `LIB_M`;
no optimization, scoring or statistics flags changed. Both logs and compiler,
source, patch, library and binary digests are in `build_manifest.json`.

`build_seeded.py` describes the patch and successful command for a fresh extraction
of the verified archive. It intentionally asserts original source fragments and
must not be reapplied to an already patched extraction. These local build paths
are host-specific; users on another host must rebuild and record their own hashes.
