# CRISPRspacer

Standalone, dependency-free Python scientific tool. Keep it independent of the
CRISPRloci API, Integration layer and other tools. Read README.md and
`docs/contracts.md` before changing search or evidence contracts.

- `inputs.py`: normalized records, source-bound occurrences and the strict
  CRISPRidentify 1.2.0 input adapter. Every occurrence survives deduplication.
- `search.py`: exhaustive full-query Hamming search, both strands, circular
  coordinates, array exclusions and conservative resource accounting.
- `database.py`: immutable SQLite spacer catalog; 16-base prefixes nominate
  exact full-length matches. Imported snapshots require operator validation.
- `viral_reference.py` and `blast_process.py`: immutable source/native receipts,
  bounded catalog scans and contained offline BLAST subprocesses.
- `viral_search.py`: heuristic BLAST retrieval and source-bound full-query,
  ungapped HSP filtering/replay. Never label retrieval exhaustive.
- `io.py` and `artifacts.py`: bounded snapshots, canonical JSON, atomic bundles
  and deterministic replay. Parse/replay the verified bytes without reopening
  mutable input paths. JSON booleans and numeric values have distinct types.
- `ncbi_viral.py`: explicit operator-only pinned acquisition; never call from
  searches. `cli.py` owns commands and structured error/limit/cancel exits.

Run `bash scripts/check.sh` with local BLAST+ on PATH for the complete native
suite; otherwise native controls skip. Use synthetic fixtures only. The
independent brute-force oracle in `tests/test_search.py` checks 0–2 substitutions,
ambiguity, both strands and both topologies. Keep scientific policies, schema
versions, coordinates, deterministic ordering and unknown values explicit.
