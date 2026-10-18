# CCTK 2

The CRISPR comparison toolkit, originally developed by **Alan Collins**, provides
observed-array comparisons, spacer-sharing networks, maximum-parsimony tree
inference, constrained reconstruction, simulation, and optional BLAST/MinCED
workflows. This GPL-3.0 modernization preserves the scientific algorithms from
upstream 1.0.3 while replacing unsafe execution and ambiguous input boundaries.
`UPSTREAM.json` records the source and file hashes. Cite the original publication
in *The CRISPR Journal* (2023): https://doi.org/10.1089/crispr.2022.0080 .

Requires Python 3.11+ on Linux. Install the pinned dependencies and package:

```sh
python -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install . --no-deps
.venv/bin/cctk --help
```

The complete scientific stack is locked. BLAST+ (`blastn`, `blastdbcmd`) and
MinCED remain optional operator-provided executables via PATH or the existing
MinCED path option. They are never downloaded by an analysis. External execution
is shell-free, time limited, output bounded and contained in a process group.

## Observed comparisons

```sh
cctk network -i tests/fixtures/arrays.txt -o network-results --min-shared 1
cctk crisprdiff -a tests/fixtures/arrays.txt -o arrays.png --preordered --seed 2
cctk compare -i examples/arrays.json -o comparison-results --plot
```

Text inputs contain an array ID followed by ordered spacer symbols. These are
literal identities, not a sequence alignment. Duplicate occurrences remain in
the array; shared counts and Jaccard use sets. Thresholds apply to both edges and
components; isolated arrays are explicit.

The JSON input (`docs/contracts.md`) retains DNA occurrences, scoped source/array
IDs, unknown/+/- strand and explicit deletion slots. Only unambiguous A/C/G/T is
accepted. Exact equality is default; reverse-complement equality is an explicit
option. Source-forward order is not inferred acquisition order. `compare` creates
a new output directory atomically with complete JSON and TSV edges. `--plot`
adds a PNG with crosses for deletion slots. Original IDs remain in JSON; TSV text
starting with spreadsheet formula characters receives an apostrophe.

CRISPRdiff produces a plot, not a structured alignment. Shared spacers are colored
and linked; unique spacers are black. Default row order is now supplied order;
`--optimize-order` explicitly enables original bounded optimization. The fixed
seed defaults to 2. A sharing graph is not a phylogeny or proof of transfer.

## Retained methods

`crisprtree`, `constrain`, `evolve`, `blast`, `spacerblast`, `minced`, and `quickrun`
remain usable standalone, with `--help` for every installed command. Original
method documentation is in `docs/upstream/`; changed defaults and limits are in
`docs/modernization.md`. For example:

```sh
cctk crisprtree -a tests/fixtures/arrays.txt -x --seed 3 -r 2 A B C
cctk evolve -n 12 -s 3 -o existing-output-directory
```

Constrain requires a bifurcating tree and explicit array-to-genome mapping;
unsupported polytomies are rejected. Evolutionary procedures retain their
published assumptions and remain separate from observed comparison. No tree is
inferred by `compare`.

## Verification

Use `requirements-test.lock`, then `python -m pytest`, `ruff check src tests`, and
`python -m build`. Frozen upstream outputs cover network metrics, order, 50
pairwise alignments/module decompositions, the actual CRISPRdiff PNG, tree
reconstruction and simulation. Further tests exercise all subcommands, optional
executable boundaries, malformed inputs, bounds and descendant cleanup. CI tests
an installed wheel from a clean checkout. Reproduce goldens with
`docs/oracle.md` and `scripts/regenerate_oracle.py`. The CLI selects headless Agg
and uses a private temporary Matplotlib cache unless MPLCONFIGDIR is supplied. No result establishes CRISPR activity.


## Code navigation

- `comparison.py` validates structured DNA occurrences, preserves deletion slots
  and source identity, and publishes observed JSON/TSV/PNG outputs atomically.
- `sequence_operations.py` owns sequence alignment and sharing metrics. Network
  membership is cached per run; `NetworkEdge.shared_spacers` is a set, while the
  input spacer lists retain occurrence order and multiplicity.
- `CRISPRdiff.py` scores neighbouring rows using cached pair overlaps. Exhaustive
  search supports up to eight arrays; larger groups use the seeded swap heuristic.
  Strict improvements retain the original tie and random-walk behavior. Groups
  with no sharing return supplied order and score zero.
- `array_parsimony.py`, `CRISPRtree.py`, `constrain.py`, `evolve.py` and
  `tree_operations.py` retain upstream evolutionary assumptions and event costs;
  `plotting.py` renders both observed and evolutionary views. Simulation rejects
  deletion from fewer than two occurrences and trailer loss from an empty array
  explicitly; valid event draws and distributions remain unchanged.
- `io.py`, `policy.py` and `process.py` bound file, input and external-tool work;
  `blast.py`, `spacerblast.py`, `minced.py` and `quickrun.py` are optional workflows.
- `tests/test_science.py` compares frozen upstream alignments, modules, network,
  plot pixels, trees and simulation; `tests/test_comparison.py` checks ordered
  occurrence identity, missing evidence and publication contracts. Keep these
  scientific controls when refining code; parity does not validate biology.
