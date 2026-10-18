# SpacePHARER navigation and preservation

- `python/spacepharer/io.py`: bounded inputs and atomic new-directory outputs.
- `reference.py`: immutable checksummed native indexes, SQLite catalogue and eight-genome lazy cache.
- `pipeline.py` / `report.py`: exact occurrence grouping, native calibration receipt, sequence-bound evidence.
- `engine.py`: version-pinned subprocesses and process-group cleanup on every exit.
- `src/core/Statistics.h`: published order-sensitive fallback, geometric e-values and empirical tails.
- `src/core/RecordIndex.h`: borrowed per-query target rows in original order; reader must outlive the index.
- `src/util`: first-party native utilities; `data/*.sh` and `src/workflow` compose them.
- `lib/mmseqs`: pinned vendored search engine. Do not mechanically reformat or revise scientific algorithms.

Preserve native arithmetic, target/control populations, tie order, final nucleotide-match selection,
gapped alignments and unavailable calibration. Missing values remain missing. Do not pool groups,
collapse reverse complements or modify thresholds without an explicit scientific algorithm revision.
References are operator-owned and immutable. Never replace broad references with small fixtures.

After native edits build a candidate executable, run CTest, then point the complete Python suite at it:

```sh
ctest --test-dir build --output-on-failure
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=python \
  SPACEPHARER_NATIVE="$PWD/build/src/spacepharer-native" \
  SPACEPHARER_ORACLE=/path/to/original/spacepharer \
  SPACEPHARER_CALIBRATION_BASELINE=/path/to/prepatch/spacepharer-native \
  python -m pytest -q -p no:cacheprovider --basetemp=/tmp/spacepharer-tests
ruff check --no-cache python tests
ruff format --check --no-cache python tests
mypy python
python -m build
python scripts/check_private_data.py workspace
```

`tests/native/records.cpp` compares indexed rows against an independent full scan.
`tests/fixtures/indexed-records.json` freezes pre-refinement 2.1.1 bytes for interleaved
protein/nucleotide rows and all summary formats. It includes final-match selection,
missing targets, gapped alignments, descending coordinates and optional FDR/PAM columns.
`tests/test_native_records.py` checks these bytes at one and two threads, including
missing HSP keys at the beginning, middle and end of each match list. Baseline executable
checks establish fixture provenance; only a rebuilt candidate validates native changes.
