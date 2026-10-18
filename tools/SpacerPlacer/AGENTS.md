# SpacerPlacer maintenance

Scientific formulas, model defaults, uncertainty, spacer order, coordinates and
bias-correction provenance are compatibility contracts. Preserve exact behavior
before compacting code; use independent equivalence controls for optimizations.
Legacy `sp_model` imports and trusted-pickle/simulation entry points remain public.

## Navigation

- `spacerplacer.py`: CLI options, filesystem validation and input dispatch.
- `sp_model/input_parser.py`, `crisprfinder_to_spacerfasta.py`, `helpers/sequence.py`:
  observed spacer adapters, DNA handling and transitive clustering.
- `sp_model/mafft_integration.py`, `experiments_tools.py`: alignment, caches, tree
  preparation, simulation parameters and pooled estimation.
- `sp_model/run_experiments.py`, `run_reconstruction.py`: workflow orchestration,
  orientation comparison, bias policy and output provenance.
- `sp_model/model_tools.py`, `helpers/misc.py`: likelihood arithmetic, alignment
  repair and spacer order graphs. Preserve first-observed identities and LIFO
  topological ties; both affect reconstruction numbering.
- `sp_model/reconstruction_tree.py`, `evolution.py`, `distance_tree_constructor.py`:
  ancestral inference, transition models and phylogeny estimation.
- `sp_model/run_simulation.py`, `simulation_tree.py`, `coalescent_tree_sim.py`:
  scientific simulation and random sampling.
- `sp_model/data_classes/`: groups, arrays, tree metadata and parent links.
- `sp_model/orientation/`, `summary.py`, `spacer_visualization/`: orientation
  selection, reports and plots. Reports must be reusable without consuming inputs.
- `sp_model/bdm_likelihood_computation/`: symbolic likelihood resources; preserve
  both active packaged pickle resources and cross-Python fallback behavior.
- `tests/`: adapters, numerics, scientific provenance and reconstruction controls.
  `test_order_graph_and_summary.py` contains graph equivalence and reuse controls.

## Verification

Run `PYTHONDONTWRITEBYTECODE=1 QT_QPA_PLATFORM=offscreen python -m pytest -q
-p no:cacheprovider` in an environment containing this package's dependencies and
pytest. Keep temporary outputs outside the source tree. SpacerPlacer and evOr
both install `sp_model`; use separate environments and never combine their
source directories on `PYTHONPATH`. Do not mutate immutable installed runtimes.
