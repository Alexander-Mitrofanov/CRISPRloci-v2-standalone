# Modernization and deliberate behavior changes

All original first-party modules have moved into the installable `src/cctk`
package. Algorithms execute directly in this package, with no subprocess wrapper
around unchanged original scripts. Scientific sequence alignment, module parsing,
parsimony event weighting and tree procedures retain their original formulas and
tie choices. Preserving those procedures is intentional; modernization does not
claim a new or independently validated evolutionary model.

Changes across the full standalone surface:

* `file_handling` uses bounded regular-file I/O, strict duplicate/empty input
  rejection, typed observed-array state and checked BLAST column counts.
* `sequence_operations` retains alignment metrics, adds a matrix-size bound,
  eliminates shared mutable defaults, provides typed network interfaces, replaces
  shell-based BLAST pipelines, and fixes empty/tiny multiprocessing batches.
* `network` emits singletons and applies the threshold before components.
* `CRISPRdiff` bounds exhaustive work, defaults to supplied order and uses an
  isolated seeded generator. Its geometry is parity checked against the original.
* `CRISPRtree`, `constrain`, `array_parsimony` preserve inference algorithms while
  introducing typed command gates, resource limits, independent mutable array
  state and targeted exceptions. Constrain rejects unresolved polytomies.
* `tree_operations` retains node identities and scaling, rejects negative branch
  lengths, and preserves all-zero lengths instead of dividing by zero.
* `evolve` honors seed zero, isolates randomness and rejects impossible event
  configurations rather than looping forever. Simulations remain hypotheses.
* `colour_schemes` isolates random state; `plotting` checks dimensions and path
  safety and supports explicit deletion marks for the observed JSON interface.
* `blast`, `spacerblast`, `minced` use argument vectors, checked exit codes,
  bounded output/deadlines and descendant cleanup. A fresh isolated interpreter
  applies file-size limits before exec, avoiding threaded preexec hazards. BLAST batch FASTA and stdin
  replace shell substitution/printf/awk. Filename iteration is sorted. BLAST
  accession filters now use Python regex, not grep-specific BRE/ERE syntax.
* `quickrun` keeps the detector/comparison/tree workflow but reports no eligible
  clusters as a successful empty result rather than terminating the interpreter.
* The installed CLI rejects unknown arguments and returns explicit nonzero
  failure codes; all command Python entry points validate their own options.

New `comparison`, `io`, `process` and `policy` modules own typed observed-array
contracts, bounded I/O/execution and shared command limits. Upstream extraction
commands remain optional and are not part of the web observed-array workflow.

Attribution and GPL terms remain intact. Archive completeness is tested from
Git's tracked tree, then an sdist/wheel is built and the installed wheel exercised.
No tools download scientific references or modify installations during analysis.

Patch 2.0.2 rejects malformed structured options/strands and non-scalar Unicode
identifiers with domain errors. Observed plot labels are literal text; their
viewport no longer renders unused legacy label-measurement figures. Repeated
spacer connectors are drawn once per occurrence pair, and the 10,000-connection
limit counts those pairs. Pairwise alignment rejects nonfinite or overflowing
scores before traceback can stall. Completed comparison directories are
published with Linux `renameat2(RENAME_NOREPLACE)` so a concurrently created
output is never replaced. Retained FASTA workflows now complement every admitted
IUPAC ambiguity code when reversing DNA. Comparison schema, order, occurrence provenance,
identity policies, set metrics and evolutionary score defaults are unchanged.
