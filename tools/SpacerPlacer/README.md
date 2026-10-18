# SpacerPlacer

Standalone fixed distribution of SpacerPlacer **1.0.1.post3**, derived from the
fixed server extraction of 2026-09-27. It reconstructs ancestral CRISPR spacer
arrays, estimates spacer deletion parameters, and can infer array orientation.

The scientific implementation, bundled runtime data, examples, and regression
tests are preserved from the fixed server source. See [import provenance](PROVENANCE.md)
and the [full inherited documentation](README_UPSTREAM.md).

## Installation

Supported server target: Linux x86-64, Python 3.10–3.12. MAFFT must be on `PATH`;
Graphviz is needed for plots. Use a dedicated environment: evOr also provides a
package named `sp_model` and must be installed separately.

```bash
git clone https://github.com/Alexander-Mitrofanov/SpacerPlacer.git
cd SpacerPlacer
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -c constraints-tested.txt .
spacerplacer --help
```

The constraint file records the dependency versions used by the original server
verification. System dependencies can be provisioned through Conda or the host's
package manager. For headless plotting, set `QT_QPA_PLATFORM=offscreen`.

## Run

```bash
spacerplacer example_datasets/tomida_datasets/CRISPR1/G1.fa runs/G1 \
  --input_type spacer_fasta \
  --no_plot_reconstruction --no_plot_order_graph
```

Use `spacerplacer --help` for input formats, tree inputs, alignment workers, and
orientation settings. Outputs include detailed reconstruction CSVs and tree data.
The inherited documentation describes the output files and scientific options.

For web jobs, invoke this tool in a separate process and a fresh output directory.
Keep uploaded pickle inputs disabled and enforce CPU, memory, and wall-time limits
in the server. Sparse groups may require an explicit
`--rho_bias_failure_policy fallback_uncorrected`; the default remains strict and
fallback is recorded in the output. Do not suppress that provenance.

Input groups retain the full filename stem, including dots; duplicate stems are
rejected before alignment. An `--aligned_groups_path` cache must match every
observed input array and its spacer order. MAFFT output must preserve the same
array identities and observed spacers; unknown symbols cannot become gaps.
The gap markers `-` and `--` are reserved for alignments, not observed spacer IDs.

Supplied Newick files must contain one tree with unique named leaves and finite,
non-negative lengths for every non-root branch. Its leaves must match the arrays
being reconstructed. Missing or empty tree mappings cannot silently select an
inferred tree. The exported tree represents the actual reconstruction after
branch extension and internal node naming, so its lengths agree with the result
CSV and event exports. An absent incoming root edge is represented by length zero.

Version 1.0.1.post3 adds these integrity checks and rejects non-finite numeric
CLI options. Historical extraction verification in `provenance/validation.json`
continues to describe version 1.0.1.post2. The unused historical
`additional_data/230329_death_lh_up_to_68_lambdifyed_old.pickle` was removed;
both active likelihood resources in
`sp_model/bdm_likelihood_computation/sympy_bdm_lh_fct/` remain bundled.

## Verification

```bash
python -m pip install '.[test]'
python -m pytest -q
```

Run this suite separately from evOr to keep their `sp_model` imports isolated.

## Attribution and license

Derived from [fbaumdicker/SpacerPlacer](https://github.com/fbaumdicker/SpacerPlacer),
with the fixed server implementation preserved. Original authorship and the
[GPL-3.0 license](LICENSE) are retained.

Please cite Fehrenbach et al., *SpacerPlacer: ancestral reconstruction of CRISPR
arrays reveals the evolutionary dynamics of spacer deletions*, Nucleic Acids
Research (2024), DOI [10.1093/nar/gkae772](https://doi.org/10.1093/nar/gkae772).
