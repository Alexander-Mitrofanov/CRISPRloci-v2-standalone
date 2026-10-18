# Code map

The CLI in `src/crisprmap/cli.py` separates two offline workflows. Read
`COUNCIL.md`, `SCIENCE.md` and the relevant contract before changing evidence.

| Path in `src/crisprmap/` | Responsibility | Regression entry points |
| --- | --- | --- |
| `io.py`, `reference.py` | Bounded snapshots, strict JSON, normalized historical inventories and exclusive publication | `test_io.py`, `test_reference.py`, `test_cli.py` |
| `core.py`, `exports.py` | Complete nearest historical ties, deterministic alignments, annotation consensus and escaped reports | `test_core.py`, `test_cli.py` |
| `source_import.py`, `catalog.py` | Pinned source receipts, retained occurrences, canonical nodes, resumable imports and merges | `test_source_import.py`, `test_catalog.py` |
| `catalog_index.py`, `indexed_search.py` | Lossless partition candidates, bounded disk scratch, verified distances and all-best SQLite evidence | `test_catalog_index.py`, `test_indexed_search.py` |
| `clusters.py` | Bona-fide consensus components, resumable spanning forests and independently checked assembly | `test_clusters.py` |
| `structure.py` | Separately pinned folding hypotheses for both orientations | `test_structure.py` |
| `modern_benchmark.py` | Accession holdout, full retained exact/RC purge, scoped homology stress and bounded DP oracle | `test_modern_benchmark.py` |
| `__init__.py`, `__main__.py` | Package identity and CLI entry point | `test_cli.py` |

`benchmarks/run_benchmark.py` replays checksum-pinned published historical
outputs; it does not execute the original tool or measure its runtime.
`scripts/fetch_legacy.py` is an explicit acquisition step, never a search hook.

Keep retrieval complete before display truncation. Labels, missing annotations,
operational components, folding and reference-relative orientation have separate
meanings. Catalogue source keys and native occurrence keys must survive import;
repeat traversal follows each array's ordinal rather than native primary-key
order. A manifest/report completion marker is written only after successful work.

Run `PYTHON=.venv/bin/python bash scripts/check.sh`. Small synthetic index and
benchmark tests explicitly request a zero free-space reserve; production keeps
its two-GiB default, tested independently. Native folding tests require the
pinned ViennaRNA extra. Indexed search also tests the pure-Python distance
fallback when optional RapidFuzz is absent.
