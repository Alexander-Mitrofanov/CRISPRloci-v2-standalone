# DefenseFinder context adapter

The production entrypoint is `python -m defense_context` (or `defense-context`).
Read `CRISPRLOCI.md` for the scientific, source-identity and license contracts.

- `defense_context/pipeline.py` owns bounded input, exact model verification,
  Pyrodigal calls, six-family MacSyFinder execution, native parsing and report export.
- Preserve MacSyFinder 2.1.4, Pyrodigal 3.7.1, HMMER 3.4, model release 1.2.3 and
  `models.lock.json`. Model/rule assets are scientific evidence; never rewrite them.
- Protein FASTA order gives native MacSyFinder **global** ranks. Report gene
  `position` is **local** to its source replicon. Keep both conventions distinct.
- Shared proteins may belong to multiple systems. Keep every component assignment;
  deduplicate system membership only. Profile ties keep the first highest score in
  family order and sorted extract-path order.
- Hash CDS bases from oriented submitted DNA, preserving IUPAC ambiguity. Hash
  exact native proteins, including X translations. Keep linear topology, partial
  evidence, missing values and activity uncertainty explicit.
- `tests/test_context*.py` are the derivative regressions. Run with the pinned
  context environment, bytecode/plugin autoload/cache disabled and temporary
  outputs outside the checkout. `scripts/native_context_smoke.py` scans the real
  pinned 40 kb control; also verify empty-gene and multi-source controls.
- The retained `defense_finder*` packages are upstream functionality. The classic
  CLI requires additional dependencies plus separately installed CasFinder and
  DefenseFinder models and may perform online latest-version discovery. It is
  outside the offline context distribution; do not fetch/update it during web jobs.
- Packaging lives in `setup.py` / `requirements-context.lock`; only
  `defense_context` is shipped. Use `scripts/check_private_data.py workspace` and
  `git diff --check` before publication. Keep operations evidence in ignored
  `.ops-private/`, never in scientific outputs or public source.
